"""Kinetic model: does a transcript fire before it is degraded?

Why this module exists
----------------------
Two earlier metrics disagreed ~700x on the same design:

    equilibrium (accessibility-corrected):  32%  OFF leak  -> "it leaks"
    nucleation-only proxy:                0.047% OFF leak  -> "it holds"

Both were wrong in the same way: Kim 2019's AND is *kinetic*.  The exposed
toehold sets the strand-displacement RATE, not the equilibrium.  Equilibrium
asks "given infinite time, what fraction binds?" and answers 32%, because
Trigger A *can* eventually pay the cost of prising its binding site out of the
inhibitory hairpin.  A cell never grants infinite time: the transcript is
degraded first.  The mRNA lifetime is a kinetic filter that equilibrium cannot
see.

The model
---------
Three-step toehold-mediated strand displacement (Zhang & Winfree 2009):

    trigger + switch  <->  toehold duplex  ->  branch migration  ->  fired
                      k_on      k_off            k_bm

Steady state on the toehold intermediate gives

    k_eff = k_on * k_bm / (k_on * Kd_toe + k_bm)          [1/M/s]

    long/strong toehold : Kd_toe -> 0,  k_eff -> k_on      (saturated)
    short/weak toehold  : Kd_toe large, k_eff -> k_bm/Kd_toe
                          i.e. ~1 decade per nt, which is what Zhang & Winfree
                          measured -- but here it comes from the actual dG, so
                          it is sequence-specific rather than a length rule.

Crucially Kd_toe uses the *accessibility-corrected* energy

    dG_toe = dG_duplex(trigger : toehold) + opening_energy(toehold)

so a toehold that is buried inside the inhibitory hairpin is expensive to use
even though the duplex itself would be favourable.  This is what unifies the
two metrics above: accessibility enters the rate, and the rate competes with
degradation.

Then the transcript either fires or is degraded:

    k_obs  = k_eff * [trigger]                             [1/s]
    P_fire = k_obs / (k_obs + k_deg),   k_deg = ln2 / half_life

P_fire is the fraction of switch transcripts that produce protein -- the thing
an experiment actually measures.  The AND ratio is P_fire(+B) / P_fire(OFF).

Parameter honesty
-----------------
k_on and k_bm are order-of-magnitude literature values for DNA at 25 C, reused
for RNA at 37 C; k_bm in particular lumps a length-dependent random walk into
one constant.  Absolute P_fire values are therefore indicative.  The *ratio*
between two designs scored with the same constants is far more trustworthy than
either number alone -- rank designs with it, do not quote it as a yield.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .config import PipelineConfig
from .thermo import get_backend

_RT = 0.0019872 * 310.15          # kcal/mol at 37 C


@dataclass
class KineticParams:
    k_on: float = 3.0e6
    """Bimolecular hybridisation rate onto an exposed toehold [1/M/s].
    Zhang & Winfree 2009 measure ~3e6 for DNA; RNA is similar in magnitude."""

    k_bm: float = 1.0
    """Effective first-order rate for completing branch migration [1/s], once
    the toehold is engaged.  Lumps an N-step random walk (k_step/N^2) into one
    constant; ~1/s is conservative for an ~20-nt migration."""

    mrna_half_life_s: float = 300.0
    """Transcript half-life.  E. coli mRNA is typically 2-8 min."""

    trigger_conc_M: float = 10e-9
    """Free trigger concentration.  Should come from DE data per gene; 10 nM is
    a placeholder for a moderately expressed transcript."""


def displacement_rate(dG_toehold: float, kp: KineticParams) -> float:
    """k_eff [1/M/s] for a toehold of accessibility-corrected energy dG."""
    Kd = math.exp(dG_toehold / _RT)
    return kp.k_on * kp.k_bm / (kp.k_on * Kd + kp.k_bm)


def fire_probability(dG_toehold: float, kp: KineticParams) -> float:
    """Fraction of transcripts that fire before being degraded."""
    k_obs = displacement_rate(dG_toehold, kp) * kp.trigger_conc_M
    k_deg = math.log(2.0) / kp.mrna_half_life_s
    return k_obs / (k_obs + k_deg)


def time_to_fire_s(dG_toehold: float, kp: KineticParams) -> float:
    k_obs = displacement_rate(dG_toehold, kp) * kp.trigger_conc_M
    return float("inf") if k_obs <= 0 else 1.0 / k_obs


def toehold_dG(sw, cfg: PipelineConfig, state: str = "off") -> float:
    """Accessibility-corrected energy of the toehold Trigger A can nucleate on.

    state='off'  : only the a* gap is exposed (inhibitory hairpin shut)
    state='afterB': Trigger B has invaded k2*, so x* is released too -- scored
                    with B's own footprint held open (that is what B binding
                    means), not by assuming the result.
    """
    b = get_backend(cfg)
    s = sw.spans
    ta = sw.pair.triggerA
    if state == "off":
        site = list(range(*s["a_star_gap"]))
        trig = ta.a
        opening = b.opening_energy(sw.core, site)
    elif state == "afterB":
        site = list(range(s["xstar"][0], s["a_star_gap"][1]))
        trig = ta.x + ta.a
        given = list(range(s["r2star"][0], s["k2star"][1]))   # B bound here
        opening = b.opening_energy_conditioned(sw.core, site, given)
    else:
        raise ValueError(state)
    duplex = b.binding_dG(trig, sw.core[site[0]:site[-1] + 1])
    return duplex + opening


def and_behaviour(sw, cfg: PipelineConfig | None = None,
                  kp: KineticParams | None = None) -> dict:
    """P_fire without and with Trigger B, and the AND ratio."""
    cfg = cfg or sw.cfg
    kp = kp or KineticParams()
    dg_off = toehold_dG(sw, cfg, "off")
    dg_on = toehold_dG(sw, cfg, "afterB")
    p_off, p_on = fire_probability(dg_off, kp), fire_probability(dg_on, kp)
    return {
        "dG_toehold_off": dg_off, "dG_toehold_afterB": dg_on,
        "k_eff_off": displacement_rate(dg_off, kp),
        "k_eff_afterB": displacement_rate(dg_on, kp),
        "t_fire_off_s": time_to_fire_s(dg_off, kp),
        "t_fire_afterB_s": time_to_fire_s(dg_on, kp),
        "p_fire_off": p_off, "p_fire_afterB": p_on,
        "and_ratio": p_on / max(p_off, 1e-30),
        "half_life_s": kp.mrna_half_life_s,
    }


def report(sw, cfg: PipelineConfig | None = None, kp: KineticParams | None = None):
    cfg = cfg or sw.cfg
    kp = kp or KineticParams()
    r = and_behaviour(sw, cfg, kp)
    print("kinetic model  (k_on=%.0e /M/s, k_bm=%.1f /s, mRNA t1/2=%.0f s, "
          "[trigger]=%g nM)" % (kp.k_on, kp.k_bm, kp.mrna_half_life_s,
                                kp.trigger_conc_M * 1e9))
    print()
    print("  %-26s %12s %12s" % ("", "OFF", "after B"))
    print("  " + "-" * 52)
    print("  %-26s %12.2f %12.2f" % ("toehold dG (kcal/mol)",
                                     r["dG_toehold_off"], r["dG_toehold_afterB"]))
    print("  %-26s %12.2e %12.2e" % ("k_eff (1/M/s)", r["k_eff_off"], r["k_eff_afterB"]))
    print("  %-26s %12.0f %12.0f" % ("time to fire (s)",
                                     r["t_fire_off_s"], r["t_fire_afterB_s"]))
    print("  %-26s %11.4f%% %11.2f%%" % ("P(fire before decay)",
                                         100 * r["p_fire_off"], 100 * r["p_fire_afterB"]))
    print()
    print("  AND ratio = %.0fx" % r["and_ratio"])
    return r
