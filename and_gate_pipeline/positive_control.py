"""Positive control: a switch that is KNOWN to work must score as working.

Why this exists
---------------
Every number this project produces is model-dependent.  When the AND gate scores
badly there are two completely different explanations and, until now, no way to
tell them apart:

    (a) the design is bad          -- the gate really does not open
    (b) the measurement is bad     -- our scoring code cannot recognise an open
                                      switch even when it sees one

A positive control separates them.  We build a plain Green 2026 Series-A toehold
switch -- one hairpin, no inhibitory arm, the architecture the literature reports
working -- feed it its own cognate trigger, and push it through the SAME
primitives the AND gate uses: ``opening_energy``, ``displacement_rate``,
``spontaneous_rate``, ``fire_probability``.  If that reads OFF, the fault is in
(b) and no AND-gate number can be trusted.

This is deliberately NOT a separate scoring path.  It shares the code under test,
because a control that uses different code proves nothing about the code in use.

    python -m and_gate_pipeline.positive_control
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import sequence_utils as su
from .config import PipelineConfig
from .kinetics import (KineticParams, displacement_rate, fire_probability,
                       spontaneous_rate, time_to_fire_s)
from .thermo import get_backend
from .vista_switch import build_primary_module

# A 36-nt trigger with ordinary composition.  Nothing special about it: the
# point of the control is the ARCHITECTURE, not this sequence.
DEFAULT_TRIGGER = "AUGGCACGUUAACCGGAUUCCAUGCAUACAGGGCAU"


@dataclass
class ControlSwitch:
    """A plain Series-A switch, shaped like VistaAndSwitch so the same
    functions accept it."""
    cfg: PipelineConfig
    core: str
    spans: dict
    trigger: str
    toehold_len: int
    pair: object = None                       # unused; kept for API symmetry
    _extra: dict = field(default_factory=dict)

    def seq_of(self, name: str) -> str:
        s, e = self.spans[name]
        return self.core[s:e]


def build_control(cfg: PipelineConfig | None = None, trigger: str | None = None,
                  reporter: str = "") -> ControlSwitch:
    """Green 2026 Series-A switch for ``trigger``, built by VISTA's own code."""
    cfg = cfg or PipelineConfig()
    trig = su.to_rna(trigger or DEFAULT_TRIGGER)
    need = cfg.resolved_len_r1() + cfg.Lx + cfg.len_a + cfg.len_k1
    if len(trig) != need:
        raise ValueError(f"control trigger must be {need} nt for this config, "
                         f"got {len(trig)}")
    core, spans, toehold_len = build_primary_module(trig, cfg, reporter)

    # locate RBS + start codon exactly as the AND-gate builder does
    top = su.to_rna(cfg.hairpin_top)
    tp0 = spans["top"][0]
    r = top.find(su.to_rna(cfg.rbs_seq))
    if r >= 0:
        spans["rbs"] = (tp0 + r, tp0 + r + len(cfg.rbs_seq))
    a = top.rfind("AUG")
    if a >= 0:
        spans["start_codon"] = (tp0 + a, tp0 + a + 3)
    return ControlSwitch(cfg=cfg, core=core, spans=spans, trigger=trig,
                         toehold_len=toehold_len)


def score_control(sw: ControlSwitch, kp: KineticParams | None = None) -> dict:
    """Score the control with the AND gate's own primitives.

    A plain switch has two states, not four: no trigger, and trigger present.
    Its whole toehold is exposed -- there is no inhibitory hairpin masking it --
    so this measures the best case the architecture can produce.
    """
    cfg = sw.cfg
    kp = kp or KineticParams()
    b = get_backend(cfg)

    t0, t1 = sw.spans["toehold"]
    site = list(range(t0, t1))
    duplex = b.binding_dG(sw.trigger, sw.core[t0:t1])
    opening = b.opening_energy(sw.core, site)
    dg_toe = duplex + opening

    k_spont = spontaneous_rate(sw, cfg, kp)
    k_trig = displacement_rate(dg_toe, kp) * kp.conc("A")

    p_off = fire_probability_from_rate(k_spont, kp)
    p_on = fire_probability_from_rate(k_spont + k_trig, kp)
    return {
        "toehold_nt": t1 - t0,
        "dG_duplex": duplex, "dG_opening": opening, "dG_toehold": dg_toe,
        "k_spont": k_spont, "k_trigger": k_trig,
        "t_fire_s": time_to_fire_s(dg_toe, kp),
        "P_off": p_off, "P_on": p_on,
        "on_off": p_on / max(p_off, 1e-30),
    }


def fire_probability_from_rate(k_obs: float, kp: KineticParams) -> float:
    import math
    k_deg = math.log(2.0) / kp.mrna_half_life_s
    return k_obs / (k_obs + k_deg)


def check(cfg: PipelineConfig | None = None, kp: KineticParams | None = None,
          min_on: float = 0.5, max_off: float = 0.05) -> tuple[bool, dict]:
    """True when the control behaves like a working switch.

    ``min_on`` / ``max_off`` are deliberately loose.  This is not a performance
    benchmark -- it asks only whether the scoring code can tell ON from OFF at
    all.  A tight threshold here would just re-introduce a guessed number.
    """
    sw = build_control(cfg)
    r = score_control(sw, kp)
    ok = r["P_on"] >= min_on and r["P_off"] <= max_off
    return ok, r


def report(cfg: PipelineConfig | None = None, kp: KineticParams | None = None):
    cfg = cfg or PipelineConfig()
    kp = kp or KineticParams()
    sw = build_control(cfg)
    r = score_control(sw, kp)
    print("POSITIVE CONTROL -- plain Green Series-A switch, no inhibitory hairpin")
    print("scored with the same primitives as the AND gate\n")
    print(f"  switch length        {len(sw.core)} nt")
    print(f"  exposed toehold      {r['toehold_nt']} nt   (nothing masks it)")
    print(f"  dG duplex            {r['dG_duplex']:8.2f} kcal/mol")
    print(f"  dG opening cost      {r['dG_opening']:8.2f}")
    print(f"  dG toehold (net)     {r['dG_toehold']:8.2f}")
    print()
    print(f"  k_spont              {r['k_spont']:.3e} 1/s")
    print(f"  k_trigger            {r['k_trigger']:.3e} 1/s")
    print(f"  time to fire         {r['t_fire_s']:8.0f} s   "
          f"(mRNA t1/2 = {kp.mrna_half_life_s:.0f} s)")
    print()
    print(f"  P(fire) no trigger   {100 * r['P_off']:9.4f}%")
    print(f"  P(fire) + trigger    {100 * r['P_on']:9.4f}%")
    print(f"  ON/OFF               {r['on_off']:9.0f}x")
    ok = r["P_on"] >= 0.5 and r["P_off"] <= 0.05
    print()
    print("  VERDICT: " + ("PASS -- the scoring code recognises a working switch"
                           if ok else
                           "FAIL -- the scoring code cannot see a known-good "
                           "switch open.\n           No AND-gate number is "
                           "trustworthy until this passes."))
    return r


if __name__ == "__main__":
    r = report()
    raise SystemExit(0 if (r["P_on"] >= 0.5 and r["P_off"] <= 0.05) else 1)
