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
def _planted_genes(seed=0, perfect=True):
    rng = random.Random(seed)
    r = lambda n: "".join(rng.choice("ACGU") for _ in range(n))
    x = "ACGUACGUACGU"
    g1 = r(25) + "GGA" + x + "CAC" + "UUCAGG" + r(25)
    rc = su.reverse_complement(x)
    if not perfect:
        rc = rc[:-1] + ("A" if rc[-1] != "A" else "C")   # one mismatch
    g2 = r(20) + rc + r(20)
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
