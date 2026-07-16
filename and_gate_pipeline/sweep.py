"""Parameter sweeps for the AND/ON trade-off.

The governing equation (measured, see truth_table.py):

    toehold available to Trigger A, OFF      = |a|
    toehold available to Trigger A, after B  = |a| + Lx

Trigger B is complementary only to r2* and k2*, so it can displace exactly the
Lx-nt k2*:x* helix and never the r1:r1* helix -- whatever |r1| is.  So |r1|
sets the OFF lock and the post-nucleation binding energy, but NOT the toehold
Trigger A gets to grab.  Only |a| and Lx move that.

    AND ratio = occupancy(|a| + Lx) / occupancy(|a|)
    ON state  ~ occupancy(|a| + Lx)
    independence ~ dG(x:k2), which grows with Lx

so Lx is pulled in two directions at once and |a| sets the OFF floor.

    python -m and_gate_pipeline.sweep
"""

from __future__ import annotations

import math

from .config import PipelineConfig
from .kinetics import KineticParams, and_behaviour
from .target_scan import scan_both_orientations
from .thermo import get_backend
from .truth_table import _RT
from .vista_switch import build


def free_fraction(dG: float, C: float) -> float:
    """Fraction of each mRNA left unpaired if A and B can dimerise via x:k2."""
    Kd = math.exp(dG / _RT)
    disc = (2 * C + Kd) ** 2 - 4 * C * C
    x = ((2 * C + Kd) - math.sqrt(max(0.0, disc))) / 2
    return (C - x) / C


def evaluate(cfg: PipelineConfig, g1: str, g2: str, conc: float = 10e-9,
             kp: KineticParams | None = None):
    """Metrics for the best exact pair under this config, or None.

    Scored kinetically (see kinetics.py): P_fire is the fraction of transcripts
    that fire before being degraded.  Equilibrium occupancy is NOT used -- it
    cannot express Kim's AND, which lives in the displacement rate.
    """
    b = get_backend(cfg)
    kp = kp or KineticParams(trigger_conc_M=conc)
    pairs = [p for p in scan_both_orientations(g1, g2, cfg) if p.exact]
    if not pairs:
        return None
    p = pairs[0]
    sw = build(p, cfg)
    ta = p.triggerA
    k = and_behaviour(sw, cfg, kp)
    dg_ab = b.binding_dG(ta.x, p.triggerB.k2)
    s = sw.spans
    stem = b.mfe(sw.core[s["k2star"][0]:s["xstar"][1]])[1]
    return {
        "n_exact": len(pairs), "off": k["p_fire_off"], "on": k["p_fire_afterB"],
        "ratio": k["and_ratio"], "dg_ab": dg_ab,
        "free": free_fraction(dg_ab, conc), "stem_mfe": stem,
        "t_off": k["t_fire_off_s"], "t_on": k["t_fire_afterB_s"],
        "toehold_off": cfg.len_a, "toehold_onB": cfg.len_a + cfg.Lx,
    }


def sweep(g1: str, g2: str, conc: float = 10e-9):
    kp = KineticParams(trigger_conc_M=conc)
    print("AND/ON trade-off, scored kinetically (kinetics.py).")
    print("P_fire = fraction of transcripts that fire before decay "
          "(t1/2=%.0fs, [trigger]=%g nM)." % (kp.mrna_half_life_s, conc * 1e9))
    print("free = fraction of each mRNA left unpaired by the A:B duplex.")
    print("eff ON = free x P_fire(+B) -- a dimerised trigger cannot act.\n")

    hdr = ("%-4s %-4s %-6s %9s %9s %8s %8s %7s %8s"
           % ("Lx", "|a|", "exact", "leak OFF", "P_fire +B", "AND x", "dG(A:B)",
              "free", "eff ON"))
    for title, grid in (("sweep Lx  (|a| = 4)", [(lx, 4) for lx in (4, 6, 8, 10, 12)]),
                        ("sweep |a| (Lx = 6)", [(6, a) for a in (2, 3, 4, 5, 6)])):
        print(title); print(hdr); print("-" * len(hdr))
        for lx, a in grid:
            r1 = 30 - lx - a
            if r1 < 4:
                continue
            cfg = PipelineConfig(Lx=lx, len_a=a, len_k1=6, L_A=r1 + lx + a + 6,
                                 L_B=lx + 24, len_r2=24)
            m = evaluate(cfg, g1, g2, conc)
            if m is None:
                print("%-4d %-4d %-6s %8s" % (lx, a, "none", "-")); continue
            # a trigger that is dimerised with its partner cannot act, so the
            # ON state you actually get is free-fraction x occupancy
            eff = m["free"] * m["on"]
            print("%-4d %-4d %-6d %8.3f%% %8.3f%% %8.0f %8.1f %6.0f%% %7.2f%%"
                  % (lx, a, m["n_exact"], 100 * m["off"], 100 * m["on"],
                     m["ratio"], m["dg_ab"], 100 * m["free"], 100 * eff))
        print()


if __name__ == "__main__":
    from . import examples
    sweep(examples.GENE1, examples.GENE2)
