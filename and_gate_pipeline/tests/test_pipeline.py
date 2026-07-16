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


# ---- interop with the standalone scanner ----------------------------------- #
def _scanner_available():
    try:
        from and_gate_pipeline.interop import load_scanner
        load_scanner()
        return True
    except Exception:
        return False


def test_interop_shim_is_not_the_broken_backend():
    """The scanner's own NUPACK unpaired_probs() silently returns 0.5 for every
    base; the pipeline must never inherit that when it drives the scanner."""
    if not _scanner_available():
        return
    from and_gate_pipeline.interop import _ScannerBackendShim
    shim = _ScannerBackendShim(BK, CFG.temperature_c)
    up = shim.unpaired_probs("GGGAAACUUCGGAUCCGAAGUUUCCC")
    assert not all(x == 0.5 for x in up), "shim fell back to the 0.5 stub"
    assert any(x > 0.0 for x in up)
    struct, e_mfe, e_ens = shim.fold("GGGAAACUUCGGAUCCGAAGUUUCCC")
    assert e_ens <= e_mfe + 1e-6            # ensemble energy is never above MFE
    assert 0.0 <= shim.structure_probability("GGGAAACUUCGGAUCCGAAGUUUCCC", 37.0) <= 1.0


def test_interop_config_roundtrip_satisfies_equations():
    if not _scanner_available():
        return
    from and_gate_pipeline.interop import load_scanner, config_from_params
    mod = load_scanner()
    cfg = config_from_params(mod.Params())
    rep = validate_config(cfg)
    assert rep.ok, rep.errors
    p = mod.Params()
    assert cfg.Lx == p.star_len and cfg.len_k1 == p.K1_len
    assert cfg.resolved_len_r1() == p.R1_len
    assert cfg.resolved_len_r2() == p.R2_len


def test_interop_window_conversion_maps_domains():
    if not _scanner_available():
        return
    from and_gate_pipeline.interop import load_scanner, windows_to_pair
    mod = load_scanner()
    gene_a = "AAAACCCCGGGGUUUUACGUACGUACGUACGUACGUACGUAAAA"
    gene_b = "UUUUGGGGCCCCAAAAUGCAUGCAUGCAUGCAUGCAUGCAUUUU"
    w1 = mod.Trigger1Window(gene="A", start=2, seq="", r1="AACCCC",
                            star="GGGGUU", a="UUAC", k1="GUACGU",
                            mean_unpaired=0.5)
    w2 = mod.Trigger2Window(gene="B", start=1, seq="", r2="UUUGGG",
                            k2=su.reverse_complement("GGGGUU"),
                            mean_unpaired=0.5)
    pair = windows_to_pair(w1, w2, {"A": gene_a, "B": gene_b})
    assert pair.triggerA.x == "GGGGUU"          # star -> x
    assert pair.triggerA.pos_x == 2 + len("AACCCC")
    assert pair.triggerB.pos_k2 == 1 + len("UUUGGG")
    assert pair.hamming == 0 and pair.exact
    assert pair.triggerA.gene == su.to_rna(gene_a)   # sequence, not the name
    assert pair.meta["gene_a_name"] == "A"


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
