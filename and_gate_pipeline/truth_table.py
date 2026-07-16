"""The AND behaviour, measured as nucleation probability.

IMPORTANT -- why this is not an MFE table.  Kim 2019's AND is a *nucleation*
effect: with a short exposed gap ('a'=4) the primary trigger has "reduced
binding probability" and essentially never grabs on; once the secondary trigger
opens the inhibitory hairpin, enough toehold is exposed that it does.  The final
bound state is favourable either way.

MFE cofold has no concentration term -- it places every strand in one complex
and returns the lowest-energy structure -- so it will *always* show Trigger A
binding, regardless of how little toehold is exposed.  An MFE truth table
therefore cannot express this AND and will always look like a leak.

So we measure what Kim actually means: the equilibrium occupancy of the
*nucleation* step at a realistic cellular concentration, before vs after
Trigger B.

    python -m and_gate_pipeline.truth_table
"""

from __future__ import annotations

import math

import RNA

from . import sequence_utils as su
from .config import PipelineConfig
from .thermo import get_backend

_RT = 0.0019872 * 310.15          # kcal/mol at 37 C


def occupancy(dG: float, conc_M: float) -> float:
    """Fraction of switch bound at ``conc_M`` trigger, for a nucleation dG."""
    Kd = math.exp(dG / _RT)
    return conc_M / (conc_M + Kd)


def and_ratio(sw, cfg: PipelineConfig | None = None, conc_M: float = 10e-9):
    """Nucleation occupancy for Trigger A before vs after Trigger B.

    OFF  : only ``a*`` (Kim's exposed gap) is available.
    +B   : Trigger B has invaded ``k2*``, additionally freeing ``x*``.
    Returns (occ_off, occ_afterB, ratio).
    """
    cfg = cfg or sw.cfg
    b = get_backend(cfg)
    ta = sw.pair.triggerA
    off_th, off_tr = sw.seq_of("a_star_gap"), ta.a
    onB_th, onB_tr = sw.seq_of("xstar") + sw.seq_of("a_star_gap"), ta.x + ta.a
    o_off = occupancy(b.binding_dG(off_tr, off_th), conc_M)
    o_on = occupancy(b.binding_dG(onB_tr, onB_th), conc_M)
    return o_off, o_on, o_on / max(o_off, 1e-30)


def _paired_in(region, joined, md, offset):
    fc = RNA.fold_compound(joined, md)
    ss, e = fc.mfe()
    lo, hi = region
    sub = ss[offset + lo: offset + hi]
    return e, sum(1 for c in sub if c != "."), len(sub)


def truth_table(sw, cfg: PipelineConfig | None = None, verbose: bool = True):
    cfg = cfg or sw.cfg
    md = RNA.md(); md.temperature = cfg.temperature_c
    core = sw.core
    A, B = sw.triggerA, sw.triggerB
    region = sw.spans["primary_hairpin"]

    rows = []
    for label, trigs in (("switch alone", ()), ("+ A only", (A,)),
                         ("+ B only", (B,)), ("+ A + B", (A, B))):
        joined = "&".join(list(trigs) + [core])
        offset = sum(len(t) + 1 for t in trigs)
        e, paired, total = _paired_in(region, joined, md, offset)
        rows.append((label, e, paired, total))

    if verbose:
        closed = rows[0][2]
        print("primary (RBS/AUG) hairpin: %d nt; %d bp when shut" % (rows[0][3], closed))
        print()
        print("  %-14s %10s   %-18s %s" % ("condition", "MFE", "primary hairpin", "verdict"))
        print("  " + "-" * 62)
        for label, e, paired, total in rows:
            frac = paired / max(1, closed)
            verdict = "OPEN" if frac < 0.4 else ("partial" if frac < 0.8 else "shut")
            want = "OPEN" if label == "+ A + B" else "shut"
            mark = "ok" if verdict == want or (want == "shut" and verdict == "shut") else "<-- "
            print("  %-14s %10.2f   %2d/%2d bp %-9s %s %s"
                  % (label, e, paired, closed, "", verdict,
                     "" if verdict == want else "  <-- wanted %s" % want))
    return rows


def report(pair, cfg: PipelineConfig | None = None, reporter: str = ""):
    from .vista_switch import build
    cfg = cfg or PipelineConfig()
    b = get_backend(cfg)
    sw = build(pair, cfg, reporter)
    ta, tb = pair.triggerA, pair.triggerB

    print("=" * 64)
    print("Trigger A (gene 1, 5'->3' = k1-a-x-r1): %s" % ta.seq)
    print("   k1=%s a=%s x=%s r1=%s" % (ta.k1, ta.a, ta.x, ta.r1))
    print("Trigger B (gene 2, 5'->3' = k2-r2)    : %s" % tb.seq)
    print("   k2=%s r2=%s" % (tb.k2, tb.r2))
    print("hamming(revcomp(x), k2) = %d" % pair.hamming)
    print()
    print("switch core: %d nt" % len(sw.core))
    for n in ("r2star", "k2star", "r1_switch", "sec_loop", "r1star", "xstar",
              "a_star_gap", "k1star", "top", "ab_dom"):
        s, e = sw.spans[n]
        print("   %-12s %3d-%3d (%2d nt)  %s" % (n, s, e, e - s, sw.core[s:e]))
    print()

    # footprints contiguous?
    fp_A = sw.core[sw.spans["r1star"][0]:sw.spans["k1star"][1]]
    fp_B = sw.core[sw.spans["r2star"][0]:sw.spans["k2star"][1]]
    print("Trigger A footprint r1*|x*|a*|k1*  == revcomp(TriggerA)? %s"
          % (su.reverse_complement(fp_A) == ta.seq))
    print("Trigger B footprint r2*|k2*        == revcomp(TriggerB)? %s"
          % (su.reverse_complement(fp_B) == tb.seq))
    print()
    print("dG(B : switch) = %7.2f kcal/mol" % b.binding_dG(tb.seq, sw.core))
    print("dG(A : switch) = %7.2f kcal/mol" % b.binding_dG(ta.seq, sw.core))
    print("dG(A : B)      = %7.2f kcal/mol   (the x:k2 liability)" % b.binding_dG(ta.seq, tb.seq))
    up = b.unpaired_probabilities(sw.core)
    s, e = sw.spans["r2star"]
    print("P(unpaired) over r2* toehold  = %.3f   (Trigger B's landing pad)"
          % (sum(up[s:e]) / (e - s)))
    s, e = sw.spans["a_star_gap"]
    print("P(unpaired) over a* gap       = %.3f   (Kim's exposed 'a')" % (sum(up[s:e]) / (e - s)))
    s, e = sw.spans["masked_toehold"]
    print("P(unpaired) over masked r1*x* = %.3f   (should be LOW in OFF)"
          % (sum(up[s:e]) / (e - s)))
    print()
    truth_table(sw, cfg)
    print()
    o_off, o_on, ratio = and_ratio(sw, cfg)
    print("AND behaviour (nucleation occupancy at 10 nM -- what Kim actually measures):")
    print("   Trigger A can grab on, OFF        : %7.4f %%   (only a* exposed)" % (100 * o_off))
    print("   Trigger A can grab on, after B    : %7.4f %%   (x* + a* exposed)" % (100 * o_on))
    print("   AND ratio                         : %.0fx" % ratio)
    print()
    print("   (the MFE table above cannot show this: cofold has no concentration,")
    print("    so it always lets Trigger A bind -- see this module's docstring)")
    return sw


if __name__ == "__main__":
    from .target_scan import scan_both_orientations
    from . import examples
    cfg = PipelineConfig()
    pairs = scan_both_orientations(examples.GENE1, examples.GENE2, cfg)
    exact = [p for p in pairs if p.exact] or pairs
    print("candidates: %d (exact: %d)" % (len(pairs), sum(p.exact for p in pairs)))
    report(exact[0], cfg)
