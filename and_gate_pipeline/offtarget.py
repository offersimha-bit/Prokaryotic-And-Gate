"""Off-target risk -- transcriptome-wide sliding-window complementarity scan
(Section 7D).

A trigger can mis-activate (or be sequestered) if it is highly complementary to
a non-target transcript.  We slide a window across each supplied transcript and
score the best complementarity between the trigger and the window's reverse
complement (i.e. how well the trigger could hybridise there).  Hits at or above
``offtarget_max_identity`` are flagged; a hit against a gene tagged essential is
treated as disqualifying.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import sequence_utils as su
from .config import PipelineConfig


@dataclass
class OffTargetHit:
    trigger: str
    transcript: str
    position: int
    identity: float
    essential: bool


def _best_complementarity(trigger_rc: str, transcript: str, window: int):
    """Best fractional identity between ``trigger_rc`` and any ``window``-length
    stretch of ``transcript`` (both 5'->3').  ``trigger_rc`` is the reverse
    complement of the trigger, so identity here == complementarity to trigger."""
    best_id, best_pos = 0.0, -1
    n = len(transcript)
    L = min(window, len(trigger_rc))
    probe = trigger_rc[:L]
    if L == 0 or n < L:
        return best_id, best_pos
    for pos in range(0, n - L + 1):
        w = transcript[pos:pos + L]
        matches = sum(1 for x, y in zip(probe, w) if x == y)
        frac = matches / L
        if frac > best_id:
            best_id, best_pos = frac, pos
            if frac == 1.0:
                break
    return best_id, best_pos


def scan_offtargets(trigger_seq: str, transcriptome: dict, cfg: PipelineConfig,
                    essential: set | None = None,
                    exclude: set | None = None) -> list[OffTargetHit]:
    """``transcriptome``: {name: sequence}.  ``essential``: names whose hit is
    disqualifying.  ``exclude``: names to skip (e.g. the trigger's own gene)."""
    essential = essential or set()
    exclude = exclude or set()
    window = cfg.offtarget_window or cfg.Lx
    trig = su.to_rna(trigger_seq)
    trig_rc = su.reverse_complement(trig)
    hits: list[OffTargetHit] = []
    for name, seq in transcriptome.items():
        if name in exclude:
            continue
        ident, pos = _best_complementarity(trig_rc, su.to_rna(seq), window)
        if ident >= cfg.offtarget_max_identity:
            hits.append(OffTargetHit(trigger=trig, transcript=name, position=pos,
                                     identity=round(ident, 3),
                                     essential=name in essential))
    hits.sort(key=lambda h: (-h.essential, -h.identity))
    return hits
