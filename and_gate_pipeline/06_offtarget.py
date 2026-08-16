"""Off-target risk -- how else could the switch's recognition domain bind?

Two scans live here.

``scan_offtargets`` (legacy, section 7D) asks whether the TRIGGER could be
sequestered: it looks for windows the trigger is complementary to, scoring
FRACTIONAL IDENTITY against ``trigger_rc``.  Two limits are worth knowing before
reading its output.  It probes ``trigger_rc[:window]`` only, so with the default
6-nt window 30 of a 36-nt trigger's nucleotides are never examined; and 85% of
6 nt rounds up to an exact 6/6 match, so it fires on ~8% of random 400-nt
transcripts as a function of length rather than of binding.  It is kept because
``05_scoring_legacy.py`` still calls it.

``scan_offtargets_energetic`` asks the other, and for an AND gate more
dangerous, question: what else could bind the SWITCH?  A transcript resembling
the trigger competes for the recognition domain and opens the hairpin without
the real input -- a false ON, which destroys the logic rather than merely
weakening it.  (Trigger sequestration remains the legacy scan's question; the
two are not interchangeable and a full treatment wants both.)

Identity is not what decides whether a duplex forms -- free energy is.  So we
compare each candidate off-target duplex against the INTENDED one:

    hit(g) = exp(-(dG_offtarget - dG_cognate) / RT)

An off-target with the same dG as the real target scores 1; one that is 3
kcal/mol weaker scores ~0.008.  This is a physical transform with no fitted
weights.  Hits are then weighted by transcript abundance (FPKM), because a
perfect off-target in a gene that is not expressed cannot sequester anything:

    load = sum over transcripts of  FPKM(t) * hit(strongest window in t)

Folding every window transcriptome-wide is not affordable, so a candidate is
only folded if it shares an exact k-mer with the trigger (``offtarget_seed_k``).
The seed is a filter on what to fold, never on what to report.
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field

from . import sequence_utils as su
from .config import PipelineConfig

RT_KCAL_PER_K = 0.00198720425864083


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


# --------------------------------------------------------------------------- #
# energetic, expression-weighted scan                                         #
# --------------------------------------------------------------------------- #
@dataclass
class EnergeticHit:
    transcript: str
    position: int
    sequence: str
    dG: float                    # duplex dG against the recognition domain
    relative_strength: float     # vs the cognate duplex; 1.0 == just as good
    fpkm: float | None
    weighted: float              # fpkm * relative_strength
    essential: bool

    def summary(self) -> dict:
        return {"transcript": self.transcript, "position": self.position,
                "dG": round(self.dG, 2),
                "relative_strength": float(f"{self.relative_strength:.4g}"),
                "fpkm": self.fpkm,
                "weighted": float(f"{self.weighted:.4g}"),
                "essential": self.essential}


@dataclass
class OffTargetReport:
    trigger: str
    cognate_dG: float
    hits: list = field(default_factory=list)
    load: float | None = None          # sum of weighted hits; None if unweighted
    strongest_dG: float | None = None
    strongest: str | None = None
    margin: float | None = None        # strongest_dG - cognate_dG; higher is safer
    essential_hits: int = 0
    status: str = "evaluated"

    def summary(self) -> dict:
        return {"cognate_dG": round(self.cognate_dG, 2),
                "n_hits": len(self.hits),
                "load": None if self.load is None else float(f"{self.load:.4g}"),
                "strongest": self.strongest,
                "strongest_dG": (None if self.strongest_dG is None
                                 else round(self.strongest_dG, 2)),
                "margin": None if self.margin is None else round(self.margin, 2),
                "essential_hits": self.essential_hits,
                "status": self.status}


def relative_binding_hit(offtarget_dG: float, cognate_dG: float,
                         temperature_c: float = 37.0) -> float:
    """Boltzmann affinity of an off-target duplex relative to the cognate one.

    The cognate duplex is the reference, so an off-target with the same dG has
    strength exactly 1.  This is a physical transform, not a fitted score
    weight: it converts an energy difference into the ratio of equilibrium
    constants at the given temperature.  The exponent is clamped so a
    pathological dG cannot overflow the float.
    """
    rt = RT_KCAL_PER_K * (273.15 + temperature_c)
    exponent = -(offtarget_dG - cognate_dG) / rt
    return math.exp(max(-600.0, min(600.0, exponent)))


class SeedIndex:
    """Candidate-directed k-mer index over the background transcriptome.

    Only k-mers that actually occur in a queried trigger are stored, so the
    index costs what the queries need rather than what the transcriptome is.
    Positions are voted on by shared k-mers; the best-voted windows are the ones
    worth folding.
    """

    def __init__(self, transcriptome: dict, k: int):
        self.records = {name: su.to_rna(seq) for name, seq in transcriptome.items()}
        self.k = int(k)
        self.index: dict[str, list] = defaultdict(list)
        self.indexed: set[str] = set()

    def prime(self, query: str) -> None:
        k = self.k
        wanted = {query[i:i + k] for i in range(max(0, len(query) - k + 1))}
        missing = wanted - self.indexed
        if not missing:
            return
        for name, seq in self.records.items():
            for pos in range(max(0, len(seq) - k + 1)):
                kmer = seq[pos:pos + k]
                if kmer in missing:
                    self.index[kmer].append((name, pos))
        self.indexed |= missing

    def candidate_windows(self, query: str, limit: int, exclude=None,
                          own: tuple | None = None) -> list:
        """Best-voted ``len(query)``-nt windows, as ``(name, start, seq, votes)``.

        ``own`` is ``(record, start, end)`` of the trigger's intended site; any
        window overlapping it is skipped, or the trigger would be reported as
        its own strongest off-target.
        """
        self.prime(query)
        exclude = exclude or set()
        k, n = self.k, len(query)
        votes: dict[tuple, int] = defaultdict(int)
        for qpos in range(max(0, n - k + 1)):
            for name, bpos in self.index.get(query[qpos:qpos + k], []):
                if name in exclude:
                    continue
                start = bpos - qpos
                if start < 0 or start + n > len(self.records[name]):
                    continue
                if own and name == own[0] and start < own[2] and start + n > own[1]:
                    continue
                votes[(name, start)] += 1
        ranked = sorted(votes.items(), key=lambda kv: (-kv[1], kv[0]))
        return [(name, start, self.records[name][start:start + n], count)
                for (name, start), count in ranked[:limit]]


def scan_offtargets_energetic(trigger_seq: str, transcriptome: dict,
                              cfg: PipelineConfig, backend,
                              essential=None, exclude=None,
                              expression_fpkm: dict | None = None,
                              own_locus: tuple | None = None,
                              index: SeedIndex | None = None) -> OffTargetReport:
    """Score how well anything else in ``transcriptome`` competes with the
    intended target for the switch's recognition domain.

    ``expression_fpkm`` maps transcript name -> FPKM.  When supplied it must
    cover every transcript that gets scored: a missing entry raises rather than
    counting as zero abundance, because silently treating an unmeasured
    transcript as harmless is exactly the error this weighting exists to avoid
    (disable via ``cfg.offtarget_require_expression``).

    Pass ``index`` to reuse one seed index across many triggers.
    """
    essential = set(essential or ())
    exclude = set(exclude or ())
    trig = su.to_rna(trigger_seq)
    recognition = su.reverse_complement(trig)
    cognate = backend.binding_dG(recognition, trig)
    report = OffTargetReport(trigger=trig, cognate_dG=cognate)

    if not transcriptome:
        report.status = "not_evaluated_no_background"
        return report

    idx = index if index is not None else SeedIndex(transcriptome,
                                                    cfg.offtarget_seed_k)
    windows = idx.candidate_windows(trig, cfg.offtarget_max_hits,
                                    exclude=exclude, own=own_locus)
    if not windows:
        report.status = "no_window_shared_a_seed"
        report.load = 0.0 if expression_fpkm else None
        return report

    # one entry per transcript, its strongest window: overlapping windows of the
    # same transcript are the same physical competitor and must not be summed
    best: dict[str, tuple] = {}
    for name, start, seq, _votes in windows:
        dG = backend.binding_dG(recognition, seq)
        prev = best.get(name)
        if prev is None or dG < prev[0]:
            best[name] = (dG, start, seq)

    for name, (dG, start, seq) in best.items():
        strength = relative_binding_hit(dG, cognate, cfg.temperature_c)
        fpkm = None
        if expression_fpkm is not None:
            if name in expression_fpkm:
                fpkm = float(expression_fpkm[name])
            elif cfg.offtarget_require_expression:
                raise ValueError(
                    f"expression table has no FPKM for off-target record {name!r}; "
                    "supply it or set cfg.offtarget_require_expression = False")
            else:
                fpkm = 0.0
        report.hits.append(EnergeticHit(
            transcript=name, position=start, sequence=seq, dG=dG,
            relative_strength=strength, fpkm=fpkm,
            weighted=(fpkm if fpkm is not None else 1.0) * strength,
            essential=name in essential))

    report.hits.sort(key=lambda h: (-h.essential, -h.weighted, h.dG))
    report.essential_hits = sum(1 for h in report.hits if h.essential)
    if expression_fpkm is not None:
        report.load = sum(h.weighted for h in report.hits)
        report.status = "evaluated_expression_weighted"
    else:
        report.status = "evaluated_unweighted"
    strongest = min(report.hits, key=lambda h: h.dG)
    report.strongest_dG = strongest.dG
    report.strongest = f"{strongest.transcript}:{strongest.position}"
    # positive margin == every competitor binds more weakly than the real target
    report.margin = strongest.dG - cognate
    return report
