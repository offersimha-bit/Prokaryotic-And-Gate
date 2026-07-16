"""Stage 2 -- thermodynamic filtering / trigger accessibility.

A trigger is only worth building a switch around if its binding footprint is
*open* (single-stranded) inside its native transcript.  Following Toehold-VISTA
we evaluate the footprint in a series of flanking windows (+/-0..100 nt) and
summarise with:

* MFE of the flanked window,
* the specified ensemble defect against the fully-open reference (SED; lower ==
  more single-stranded), and the native ensemble defect (NED),
* the mean unpaired probability over the footprint itself (accessibility).

The +/-100 nt window is the pass/fail gate, matching the VISTA finding that
flanking structure dominates trigger accessibility.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import mean

from . import sequence_utils as su
from .config import PipelineConfig
from .thermo import ThermoBackend


@dataclass
class TriggerMetrics:
    name: str
    gene: str
    start: int
    end: int
    per_flank: dict = field(default_factory=dict)   # flank -> {mfe, sed, ned}
    accessibility: float = 0.0                       # at accessibility_flank
    open_sed: float = 0.0                             # at accessibility_flank
    passes: bool = False

    def summary(self) -> dict:
        return {
            "name": self.name, "start": self.start, "end": self.end,
            "accessibility": round(self.accessibility, 4),
            "open_sed": round(self.open_sed, 4),
            "mfe_100": round(self.per_flank.get(100, {}).get("mfe", float("nan")), 2),
            "passes": self.passes,
        }


def _flanked(gene: str, start: int, end: int, flank: int):
    lo = max(0, start - flank)
    hi = min(len(gene), end + flank)
    return gene[lo:hi], start - lo, end - lo


def evaluate_trigger(name: str, gene: str, start: int, end: int,
                     cfg: PipelineConfig, backend: ThermoBackend) -> TriggerMetrics:
    gene = su.to_rna(gene)
    tm = TriggerMetrics(name=name, gene=gene, start=start, end=end)
    for flank in cfg.flanking_lengths:
        sub, off_s, off_e = _flanked(gene, start, end, flank)
        if not sub:
            continue
        struct, energy = backend.mfe(sub)
        open_sed = backend.ensemble_defect(sub, "." * len(sub))
        ned = backend.ensemble_defect(sub, struct)
        up = backend.unpaired_probabilities(sub)
        local_acc = mean(up[off_s:off_e]) if off_e > off_s else 0.0
        tm.per_flank[flank] = {
            "mfe": energy, "sed": open_sed, "ned": ned,
            "accessibility": local_acc,
        }
    gate_flank = cfg.accessibility_flank
    if gate_flank not in tm.per_flank:
        gate_flank = max(tm.per_flank) if tm.per_flank else None
    if gate_flank is not None:
        tm.accessibility = tm.per_flank[gate_flank]["accessibility"]
        tm.open_sed = tm.per_flank[gate_flank]["sed"]
    tm.passes = (tm.accessibility >= cfg.min_accessibility
                 and tm.open_sed <= cfg.max_trigger_sed)
    return tm


def evaluate_pair_triggers(pair, cfg: PipelineConfig, backend: ThermoBackend):
    """Evaluate both triggers of a :class:`TriggerPair`."""
    ta, tb = pair.triggerA, pair.triggerB
    a_start = ta.pos_x - cfg.resolved_len_r1()
    a_end = ta.pos_x + cfg.Lx + cfg.len_a + cfg.len_k1
    b_start = tb.pos_k2 - cfg.resolved_len_r2()
    b_end = tb.pos_k2 + cfg.Lx
    tmA = evaluate_trigger("TriggerA", ta.gene, a_start, a_end, cfg, backend)
    tmB = evaluate_trigger("TriggerB", tb.gene, b_start, b_end, cfg, backend)
    return tmA, tmB
