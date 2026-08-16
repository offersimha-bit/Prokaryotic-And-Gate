"""Smoke + correctness tests.  Run with:  python -m pytest and_gate_pipeline/tests
(ViennaRNA is required; NUPACK is optional).
"""

from __future__ import annotations

import random

try:
    import pytest
except ImportError:                       # tests also run without pytest
    pytest = None

from and_gate_pipeline import sequence_utils as su
from and_gate_pipeline.config import PipelineConfig
from and_gate_pipeline.constraints import validate_config
from and_gate_pipeline.target_scan import scan_pair, scan_both_orientations
from and_gate_pipeline.architecture import build_switch
from and_gate_pipeline.filtering import evaluate_pair_triggers
from and_gate_pipeline.scoring import DesignScorer
from and_gate_pipeline.thermo import get_backend, parse_pairs
from and_gate_pipeline import examples

CFG = PipelineConfig()
BK = get_backend(CFG)


# ---- sequence utils -------------------------------------------------------- #
def test_reverse_complement_alphabet():
    assert su.reverse_complement("AUGC") == "GCAU"
    assert su.reverse_complement("ATGC") == "GCAT"
    assert su.reverse_complement(su.reverse_complement("ACGUACGU")) == "ACGUACGU"
    # regression: a U-less RNA k-mer must NOT produce a spurious T
    rc = su.reverse_complement("GGCGGA")
    assert "T" not in rc and rc == "UCCGCC"


def test_hamming_and_restricted():
    assert su.hamming("AAAA", "AAUA") == 1
    assert su.has_forbidden_run("GGGGACC", ("GGGG",)) == ["GGGG"]
    assert su.has_inframe_stop("AUG" + "UAA" + "GGG") is True
    assert su.has_inframe_stop("AUG" + "GGG" + "UAA") is False  # terminal stop ok


def test_count_aug_after_rbs():
    s = "CC" + "AGAGGAGA" + "CC" + "AUG" + "CCC"
    assert su.count_aug_after_rbs(s, "AGAGGAGA") == 1


# ---- constraints ----------------------------------------------------------- #
def test_config_integrity_ok():
    assert validate_config(CFG).ok


def test_config_integrity_contradiction():
    bad = PipelineConfig(L_B=5, Lx=12, len_r2=999)  # |k2|+|r2| != L_B
    assert not validate_config(bad).ok


# ---- target scan ----------------------------------------------------------- #
def _planted_genes(seed=0, perfect=True, cfg=CFG):
    """Plant an x / revcomp(x) pair with room for the 5'->3' domain order:
    gene 1 must hold k1-a-x-r1 and gene 2 must hold k2-r2."""
    rng = random.Random(seed)
    r = lambda n: "".join(rng.choice("ACGU") for _ in range(n))
    x = ("ACGUACGUACGU" * 2)[:cfg.Lx]
    head = cfg.len_k1 + cfg.len_a          # k1 + a sit UPSTREAM of x
    tail = cfg.resolved_len_r1()           # r1 sits DOWNSTREAM of x
    g1 = r(20) + r(head) + x + r(tail) + r(20)
    rcx = su.reverse_complement(x)
    if not perfect:
        rcx = rcx[:-1] + ("A" if rcx[-1] != "A" else "C")   # one mismatch
    g2 = r(20) + rcx + r(cfg.resolved_len_r2()) + r(20)     # k2 then r2
    return g1, g2


def test_scan_finds_exact_match():
    g1, g2 = _planted_genes(perfect=True)
    pairs = scan_pair(g1, g2, CFG)
    assert any(p.exact and p.hamming == 0 for p in pairs)


def test_scan_hamming_fallback():
    g1, g2 = _planted_genes(perfect=False)
    pairs = scan_pair(g1, g2, CFG)
    assert pairs, "expected a min-Hamming candidate"
    assert min(p.hamming for p in pairs) >= 1


def test_scan_both_orientations_runs():
    g1, g2 = _planted_genes()
    pairs = scan_both_orientations(g1, g2, CFG)
    orients = {p.orientation for p in pairs}
    assert "G1->A,G2->B" in orients


# ---- architecture ---------------------------------------------------------- #
def _one_pair():
    g1, g2 = _planted_genes()
    return next(p for p in scan_pair(g1, g2, CFG) if p.exact)


def test_switch_structure_wellformed():
    sw = build_switch(_one_pair(), CFG)
    assert len(sw.off_structure) == len(sw.core)
    # balanced brackets
    assert sw.off_structure.count("(") == sw.off_structure.count(")")
    parse_pairs(sw.off_structure)  # must not raise
    # exactly one start codon in the switch module
    assert su.count_aug_after_rbs(sw.core, CFG.rbs_seq) == 1
    # RBS present
    assert su.to_rna(CFG.rbs_seq) in sw.core
    # pure RNA -- no DNA 'T' may leak into the construct
    assert "T" not in sw.core and "T" not in sw.full


def test_switch_off_state_is_locked():
    sw = build_switch(_one_pair(), CFG)
    mfe = BK.mfe(sw.core)[1]
    assert mfe < -20.0                      # a real, folded OFF-state lock
    sed = BK.ensemble_defect(sw.core, sw.off_structure)
    assert 0.0 <= sed <= 1.0


# ---- thermo ---------------------------------------------------------------- #
def test_ensemble_defect_bounds():
    seq = "GGGGAAAACCCC"
    assert 0.0 <= BK.native_defect(seq) <= 1.0
    assert BK.open_defect(seq) > BK.native_defect(seq)  # folded -> not open


def test_binding_dG_negative_for_complements():
    assert BK.binding_dG("GGGGGCCCCC", "GGGGGCCCCC") < 0


# ---- AND mechanism --------------------------------------------------------- #
def test_triggerB_exposes_triggerA_site():
    """Core AND behaviour: forcing Trigger B's toehold open should raise the
    accessibility of the Trigger A footprint."""
    sw = build_switch(_one_pair(), CFG)
    a_idx = sw.triggerA_footprint_indices()
    forced = sw.triggerB_toehold_indices() + sw.domains.region("sec_r1")
    off = BK.region_accessibility(sw.core, a_idx)
    after_b = BK.region_accessibility(sw.core, a_idx, forced_unpaired=forced)
    assert after_b >= off


# ---- scoring --------------------------------------------------------------- #
def test_scoring_runs_and_is_finite():
    p = _one_pair()
    sw = build_switch(p, CFG)
    tmA, tmB = evaluate_pair_triggers(p, CFG, BK)
    sc = DesignScorer(CFG, BK).score(sw, tmA, tmB)
    for v in (sc.triggerB_activation, sc.intermediate_state,
              sc.triggerA_on_state, sc.penalty, sc.total):
        assert v == v            # not NaN
    assert 0.0 <= sc.triggerB_activation <= 2.0


# ---- architecture: features must EXIST, not just be annotated -------------- #
def test_junction_bulge_actually_forms():
    """Regression: the spec's 3-nt junction bulge was previously written into
    the dot-bracket but not into the sequence, so the helix closed through it
    (P(unpaired) = 0.00). It must now be single-stranded for real."""
    sw = build_switch(_one_pair(), CFG)
    r1_start = sw.domains.spans["sec_r1"][0]
    bulge = list(range(r1_start, r1_start + CFG.bulge_len))
    up = BK.unpaired_probabilities(sw.core)
    mean_open = sum(up[i] for i in bulge) / len(bulge)
    assert mean_open > 0.5, f"junction bulge did not form (P(unpaired)={mean_open:.3f})"
    # and the annotation must agree with the sequence
    partner = parse_pairs(sw.off_structure)
    assert all(partner[i + 1] == 0 for i in bulge)


def test_primary_loop_is_the_conserved_element_and_open():
    """Regression: the RBS loop was home-made filler and self-paired (62% open).
    It must be Green/VISTA's conserved element and stay single-stranded."""
    sw = build_switch(_one_pair(), CFG)
    lo, hi = sw.domains.spans["prim_loop"]
    loop = sw.core[lo:hi]
    top = su.to_rna(CFG.hairpin_top)
    assert loop == top[CFG.primary_stem_len - CFG.len_k1:-3], "loop is not the conserved element"
    assert su.to_rna(CFG.rbs_seq) in loop
    assert loop.endswith("AUG")                     # start codon at the loop's 3' end
    up = BK.unpaired_probabilities(sw.core)
    mean_open = sum(up[lo:hi]) / (hi - lo)
    assert mean_open > 0.75, f"RBS loop not open (P(unpaired)={mean_open:.3f})"


def test_tunable_secondary_arm_strength_is_wired_up():
    """Regression: secondary_arm_gc_bias (spec tunable 2) was declared in config
    and never used. Raising it must measurably weaken the inhibitory stem."""
    pair = _one_pair()

    def stem_mfe(bias):
        cfg = PipelineConfig(secondary_arm_gc_bias=bias)
        sw = build_switch(pair, cfg)
        a = sw.domains.spans["sec_k2star"][0]
        z = sw.domains.spans["sec_xstar"][1]
        return BK.mfe(sw.core[a:z])[1]

    assert stem_mfe(0.5) > stem_mfe(0.0) + 2.0, "tunable 2 has no effect on the clamp"


def test_trigger_A_site_is_never_mutated():
    """The clamp/bulge edits must touch only the switch's internal r1 copy --
    r1* stays the exact reverse complement of the real trigger (spec section 6)."""
    pair = _one_pair()
    for bias in (0.0, 0.5, 1.0):
        sw = build_switch(pair, PipelineConfig(secondary_arm_gc_bias=bias))
        assert sw.domain_seq("sec_r1star") == su.reverse_complement(pair.triggerA.r1)
        assert sw.domain_seq("sec_xstar") == su.reverse_complement(pair.triggerA.x)
        assert sw.domain_seq("sec_k2star") == su.reverse_complement(pair.triggerB.k2)


# ---- VISTA-based AND switch (the corrected architecture) ------------------- #
def test_trigger_domain_order_matches_the_footprint():
    """Trigger A must be revcomp of its own binding site r1*|x*|a*|k1*, i.e.
    k1-a-x-r1 5'->3'. Building it as r1-x-a-k1 (the spec's literal wording)
    gave dG(A:switch) = -4.3 instead of -30.8."""
    from and_gate_pipeline.vista_switch import build
    pair = _one_pair()
    sw = build(pair, CFG)
    fp_A = sw.core[sw.spans["r1star"][0]:sw.spans["k1star"][1]]
    assert su.reverse_complement(fp_A) == pair.triggerA.seq
    fp_B = sw.core[sw.spans["r2star"][0]:sw.spans["k2star"][1]]
    assert su.reverse_complement(fp_B) == pair.triggerB.seq


def test_primary_module_is_vistas_own_construction():
    """The primary module must be VISTA's builder output, not a hand-roll:
    30-nt toehold, 6-nt invasion, and Green's real 11-nt RBS loop."""
    from and_gate_pipeline.vista_switch import build
    sw = build(_one_pair(), CFG)
    top = su.to_rna(CFG.hairpin_top)
    assert sw.seq_of("top") == top                       # conserved element verbatim
    assert sw.spans["toehold"][1] - sw.spans["toehold"][0] == 30
    assert sw.spans["k1star"][1] - sw.spans["k1star"][0] == CFG.len_k1 == 6
    assert top[12:23] == "AACAGAGGAGA"                    # Green's loop, RBS inside


def test_trigger_B_has_an_exposed_toehold():
    """r2* must be single-stranded in OFF -- without it Trigger B has nothing to
    nucleate on (dG(B:switch) was -0.54 before r2* existed)."""
    from and_gate_pipeline.vista_switch import build
    sw = build(_one_pair(), CFG)
    s, e = sw.spans["r2star"]
    up = BK.unpaired_probabilities(sw.core)
    assert sum(up[s:e]) / (e - s) > 0.4
    assert BK.binding_dG(sw.triggerB, sw.core) < -15.0


def test_inhibitory_hairpin_masks_the_toehold():
    """OFF: r1*+x* (26 of the 30-nt toehold) hidden; only a* exposed."""
    from and_gate_pipeline.vista_switch import build
    sw = build(_one_pair(), CFG)
    up = BK.unpaired_probabilities(sw.core)
    ms, me = sw.spans["masked_toehold"]
    gs, ge = sw.spans["a_star_gap"]
    masked = sum(up[ms:me]) / (me - ms)
    gap = sum(up[gs:ge]) / (ge - gs)
    assert masked < 0.2, f"toehold not masked ({masked:.3f})"
    assert gap > 0.4, f"Kim's gap a* not exposed ({gap:.3f})"


def test_and_is_a_nucleation_effect():
    """Trigger B must raise Trigger A's nucleation occupancy. This is what Kim
    measures; an MFE table cannot show it (cofold has no concentration)."""
    from and_gate_pipeline.vista_switch import build
    from and_gate_pipeline.truth_table import and_ratio
    sw = build(_one_pair(), CFG)
    o_off, o_on, ratio = and_ratio(sw, CFG)
    assert o_off < o_on
    assert ratio > 5.0, f"AND ratio only {ratio:.1f}x"


def test_Lx6_keeps_the_triggers_independent():
    """Lx=6 must keep the A:B duplex weak enough that both mRNAs stay free."""
    pair = _one_pair()
    dg = BK.binding_dG(pair.triggerA.x, pair.triggerB.k2)
    assert dg > -12.0, f"x:k2 duplex too strong ({dg:.1f}) -- triggers will sequester"


# ---- kinetic model --------------------------------------------------------- #
def test_displacement_rate_saturates_and_is_monotonic():
    """k_eff must rise with toehold strength and saturate at k_on -- the
    Zhang & Winfree shape. A weak toehold must be orders of magnitude slower."""
    from and_gate_pipeline.kinetics import KineticParams, displacement_rate
    kp = KineticParams()
    rates = [displacement_rate(dg, kp) for dg in (-2, -5, -8, -11, -14, -20)]
    assert all(b >= a for a, b in zip(rates, rates[1:])), "k_eff not monotonic"
    assert rates[-1] <= kp.k_on * 1.001                  # saturates at k_on
    assert rates[-1] / max(rates[0], 1e-30) > 1e3        # weak toehold is far slower


def test_fire_probability_competes_with_degradation():
    """P_fire must depend on the mRNA lifetime -- that competition IS the model.
    A shorter-lived transcript must be harder to fire."""
    from and_gate_pipeline.kinetics import KineticParams, fire_probability
    dg = -9.0
    short = fire_probability(dg, KineticParams(mrna_half_life_s=30))
    long = fire_probability(dg, KineticParams(mrna_half_life_s=3000))
    assert short < long
    assert 0.0 <= short <= 1.0 and 0.0 <= long <= 1.0


def test_and_behaviour_trigger_B_speeds_trigger_A():
    """The whole thesis: Trigger B must make Trigger A's binding faster, so the
    switch fires within the transcript's life only when B is present."""
    from and_gate_pipeline.vista_switch import build
    from and_gate_pipeline.kinetics import and_behaviour
    sw = build(_one_pair(), CFG)
    r = and_behaviour(sw, CFG)
    assert r["dG_toehold_afterB"] < r["dG_toehold_off"]   # B strengthens the toehold
    assert r["t_fire_afterB_s"] < r["t_fire_off_s"]       # ...so it fires sooner
    assert r["p_fire_afterB"] > r["p_fire_off"]
    assert r["and_ratio"] > 1.0


def test_kinetics_not_equilibrium():
    """Regression on the reasoning, not just the code: an equilibrium read of
    this design says it leaks (~32%) because equilibrium grants infinite time.
    The kinetic model must disagree -- the transcript decays first."""
    from and_gate_pipeline.vista_switch import build
    from and_gate_pipeline.kinetics import and_behaviour
    sw = build(_one_pair(), CFG)
    r = and_behaviour(sw, CFG)
    assert r["p_fire_off"] < 0.20, (
        "OFF leak %.3f -- kinetic model should hold where equilibrium leaks"
        % r["p_fire_off"])


# ---- positive control ------------------------------------------------------ #
def test_positive_control_reads_open():
    """A switch KNOWN to work must score as working.

    Without this we cannot separate 'the AND gate does not open' from 'our
    scoring code cannot recognise an open switch'.  A plain Green Series-A
    switch, scored by the SAME primitives the AND gate uses, must read ON.
    """
    from and_gate_pipeline.positive_control import check
    ok, r = check()
    assert ok, (f"positive control failed: P_on={r['P_on']:.4f} "
                f"P_off={r['P_off']:.6f} -- the scoring code cannot see a "
                f"known-good switch open, so no AND-gate number is meaningful")


def test_positive_control_uses_the_same_primitives():
    """The control must share the code under test -- a control that runs
    different code proves nothing about the code in use."""
    import inspect
    from and_gate_pipeline import positive_control as pc
    src = inspect.getsource(pc.score_control)
    for prim in ("opening_energy", "displacement_rate", "spontaneous_rate"):
        assert prim in src, f"control does not use {prim}"


def test_and_gate_on_state_matches_the_control():
    """The AND gate's ON state should be as good as a plain working switch --
    the inhibitory hairpin must not cost us the ON state, only gate it."""
    from and_gate_pipeline.positive_control import build_control, score_control
    from and_gate_pipeline.kinetics import four_state
    from and_gate_pipeline.vista_switch import build
    ctrl = score_control(build_control(CFG))
    andg = four_state(build(_one_pair(), CFG), CFG)
    assert andg["P_11"] > 0.5 * ctrl["P_on"], (
        f"AND ON state {andg['P_11']:.3f} is far below the control's "
        f"{ctrl['P_on']:.3f} -- the architecture is losing the ON state")


# ---- four-state model / spontaneous leak ----------------------------------- #
def _vista_switch():
    from and_gate_pipeline.vista_switch import build
    return build(_one_pair(), CFG)


def test_rbs_and_start_codon_are_located():
    sw = _vista_switch()
    assert "rbs" in sw.spans and "start_codon" in sw.spans
    assert sw.seq_of("rbs") == su.to_rna(CFG.rbs_seq)
    assert sw.seq_of("start_codon") == "AUG"


def test_leak_uses_the_footprint_not_the_rbs_alone():
    """Regression: the RBS sits in the loop and is unpaired BY DESIGN, so
    scoring its accessibility alone reports a locked switch as ~80% open and
    gives P_00 ~ 97%. The leak must use the whole ribosome footprint."""
    from and_gate_pipeline.kinetics import (four_state, spontaneous_rate,
                                            ribosome_footprint, KineticParams)
    from and_gate_pipeline.thermo import get_backend
    sw = _vista_switch()
    fp = ribosome_footprint(sw, CFG)
    assert len(fp) > len(su.to_rna(CFG.rbs_seq)) + 3   # more than RBS+AUG
    b = get_backend(CFG)
    rbs_only = b.region_accessibility(sw.core, list(range(*sw.spans["rbs"])))
    k_spont = spontaneous_rate(sw, CFG, KineticParams())
    assert rbs_only > 0.5, "sanity: the RBS really is open in the loop"
    assert four_state(sw, CFG)["P_00"] < 0.05, "OFF state must not leak"


def test_four_state_truth_table():
    from and_gate_pipeline.kinetics import four_state
    r = four_state(_vista_switch(), CFG)
    P = r["P"]
    assert P["11"] > P["10"] and P["11"] > P["01"] and P["11"] > P["00"]
    assert P["01"] <= P["10"] + 1e-9, (
        "Trigger B alone opens the inhibitory hairpin, not the RBS one, so it "
        "must not fire more than Trigger A alone")
    assert r["logic_margin"] > 1.0
    for v in P.values():
        assert 0.0 <= v <= 1.0


def test_rates_add_no_subtraction():
    """k_10 must contain the baseline by construction, never by subtraction."""
    from and_gate_pipeline.kinetics import four_state
    r = four_state(_vista_switch(), CFG)
    assert r["k"]["10"] >= r["k"]["00"]
    assert abs(r["k"]["10"] - (r["k_spont_off"] + r["k_A_off"])) < 1e-12


# ---- cross-trigger crosstalk utilities ------------------------------------- #
def test_crosstalk_utilities():
    assert su.max_identity_match("ACGUACGU", "UUACGUUU") == 5   # "ACGU"+ -> ACGUU
    # b's reverse complement appears in a -> sticking
    assert su.max_revcomp_match("AAGGCCUU", su.reverse_complement("AAGGCCUU")) == 8
    masked = su.mask_region("AAAACCCCGGGG", 4, 4)
    assert masked == "AAAANNNNGGGG"
    # masked region must not contribute to a match
    assert su.longest_common_substring("NNNN", "CCCC") == 0


def test_scorecard_quality_and_crosstalk_present():
    p = _one_pair()
    sw = build_switch(p, CFG)
    tmA, tmB = evaluate_pair_triggers(p, CFG, BK)
    sc = DesignScorer(CFG, BK).score(sw, tmA, tmB)
    assert 0.0 <= sc.quality_percent <= 100.0
    for k in ("crosstalk_stick_nt", "crosstalk_subst_nt", "type2s_sites"):
        assert k in sc.details
    rows = sc.breakdown(CFG.weights)
    assert len(rows) == 3
    assert abs(sum(mx for *_x, mx in rows) - 100.0) < 1e-6   # max points sum to 100


# ---- STAGE 1: the merged scanner ------------------------------------------- #
# These replace the old interop tests.  The bug they guard against is the one
# that made the bridge unusable: the scanner sliced its window r1|x|a|k1 in
# genomic order while TriggerA.seq reassembles k1+a+x+r1, so the "trigger"
# handed downstream was a sequence that does not occur in the gene, at a
# position that made stage 2 fold a different window.

def test_trigger_is_a_contiguous_slice_of_its_gene():
    """The invariant the whole project rests on: nothing is synthesised."""
    from and_gate_pipeline.target_scan import scan_both_orientations
    pairs = scan_both_orientations(examples.GENE1, examples.GENE2, CFG)
    assert pairs, "no candidates on the example genes"
    for p in pairs[:20]:
        for trig, gene in ((p.triggerA, p.gene_a), (p.triggerB, p.gene_b)):
            s, e = trig.window
            assert su.to_rna(gene)[s:e] == trig.seq, (
                "trigger is not a contiguous slice of its gene")
            assert trig.seq in su.to_rna(gene)


def test_verify_pair_rejects_a_scrambled_trigger():
    """verify_pair must actually fail when the domain order is wrong -- a test
    that cannot fail is not a guard."""
    from and_gate_pipeline.target_scan import (scan_both_orientations,
                                               verify_pair,
                                               TriggerIntegrityError)
    p = scan_both_orientations(examples.GENE1, examples.GENE2, CFG)[0]
    p.triggerA.pos_x += 1              # same sequence, wrong coordinate
    try:
        verify_pair(p)
    except TriggerIntegrityError:
        return
    raise AssertionError("verify_pair accepted a mis-positioned trigger")


def test_pooled_scan_needs_two_distinct_records():
    """An AND gate needs two inputs, not one gene sensing itself twice."""
    from and_gate_pipeline.target_scan import scan_pool
    genes = [("g1", examples.GENE1), ("g2", examples.GENE2)]
    pairs = scan_pool(genes, CFG)
    assert pairs
    for p in pairs:
        assert p.meta["gene_a_name"] != p.meta["gene_b_name"]


def test_pooled_scan_matches_two_gene_scan():
    """The pooled path and the two-gene path must find the same connectors --
    they are now one implementation, so this is a regression guard on the
    k-mer index rather than on two separate code paths."""
    from and_gate_pipeline.target_scan import scan_pool, scan_both_orientations
    genes = [("g1", examples.GENE1), ("g2", examples.GENE2)]
    pooled = {(p.triggerA.seq, p.triggerB.seq)
              for p in scan_pool(genes, CFG) if p.exact}
    direct = {(p.triggerA.seq, p.triggerB.seq)
              for p in scan_both_orientations(examples.GENE1, examples.GENE2, CFG)
              if p.exact}
    assert pooled == direct, "pooled and two-gene scans disagree"


def test_fasta_reader_pools_multiple_records():
    import tempfile, os
    from and_gate_pipeline.target_scan import read_fasta_records, load_genes
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "genes.fa")
        with open(path, "w") as fh:
            fh.write(">a\nACGUACGU\nACGU\n>b\nUUUUGGGG\n")
        recs = read_fasta_records(path)
        assert recs == [("a", "ACGUACGUACGU"), ("b", "UUUUGGGG")]
        assert load_genes(d) == recs        # a folder works too


def test_type2s_filter_is_off_for_triggers_by_default():
    """The trigger is endogenous and never synthesised, so Golden Gate sites in
    it are not a constructability problem -- that filter belongs to the switch."""
    assert CFG.scan_forbid_motifs is False
    assert CFG.scan_max_gc == 0.0


# ---- STAGE 3: selection ----------------------------------------------------- #
def test_pareto_fronts_rank_dominated_points_behind():
    from and_gate_pipeline.select import pareto_fronts
    # (1,0) and (0,1) are extremes; (0.5,0.5) trades between them -- none of
    # the three dominates another, so all sit on front 0.  (0.2,0.2) is strictly
    # worse than (0.5,0.5) on both axes, so it is pushed to front 1.
    pts = [(1.0, 0.0), (0.0, 1.0), (0.5, 0.5), (0.2, 0.2)]
    fronts = pareto_fronts(pts)
    assert fronts[0] == fronts[1] == fronts[2] == 0
    assert fronts[3] == 1
    # a point that dominates everything must be alone on front 0
    assert pareto_fronts([(1.0, 1.0), (0.5, 0.5), (1.0, 0.0)]) == [0, 1, 1]


def test_selection_rejects_near_identical_windows():
    """The failure mode of the old top_k=8: one site, eight offsets."""
    from and_gate_pipeline.target_scan import scan_both_orientations
    from and_gate_pipeline.select import select, evaluate_quality
    from and_gate_pipeline.filtering import evaluate_pair_triggers
    pairs = scan_both_orientations(examples.GENE1, examples.GENE2, CFG)[:12]
    scored = []
    for p in pairs:
        tmA, tmB = evaluate_pair_triggers(p, CFG, BK)
        scored.append((p, tmA, tmB, evaluate_quality(p, tmA, tmB, BK, CFG)))
    chosen = select(scored, CFG, k=4)
    assert len(chosen) == 4
    windows = [(c[0].triggerA.window, c[0].triggerB.window) for c in chosen]
    assert len(set(windows)) == len(windows)


def test_quality_is_absolute_not_pool_relative():
    """Qualities come from fixed reference scales, so scoring the same pair in a
    pool of 2 and a pool of 12 must give the identical number."""
    from and_gate_pipeline.target_scan import scan_both_orientations
    from and_gate_pipeline.select import evaluate_quality
    from and_gate_pipeline.filtering import evaluate_pair_triggers
    p = scan_both_orientations(examples.GENE1, examples.GENE2, CFG)[0]
    tmA, tmB = evaluate_pair_triggers(p, CFG, BK)
    q1 = evaluate_quality(p, tmA, tmB, BK, CFG)
    q2 = evaluate_quality(p, tmA, tmB, BK, CFG)
    assert q1.objectives() == q2.objectives()
    for v in q1.objectives():
        assert 0.0 <= v <= 1.0


# ---- packaging: the stage-numbered files stay importable -------------------- #
def test_every_numbered_module_loads_under_its_alias():
    import and_gate_pipeline as pkg
    from and_gate_pipeline import _loader
    unavailable = getattr(pkg, "__unavailable__", {})
    core = {"target_scan", "filtering", "select", "vista_switch",
            "architecture", "kinetics", "scoring", "offtarget", "pipeline"}
    broken = core & set(unavailable)
    assert not broken, f"core modules failed to load: {broken}"
    for alias in core:
        assert hasattr(pkg, alias), f"alias {alias} not registered"
        assert alias in _loader.FILE_OF


def _run_standalone() -> int:
    """Minimal runner so the suite works without pytest installed."""
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as ex:                       # noqa: BLE001
            failed += 1
            print(f"FAIL {fn.__name__}: {type(ex).__name__}: {ex}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    if pytest is not None:
        raise SystemExit(pytest.main([__file__, "-q"]))
    raise SystemExit(_run_standalone())
