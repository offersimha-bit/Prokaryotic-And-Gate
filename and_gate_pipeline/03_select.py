"""STAGE 3 -- candidate selection: which pairs are worth the expensive scoring.

Stage 1 finds every connector coincidence; stage 2 measures how open each
trigger is in its own transcript.  Neither of those decides *which* candidates
get built and kinetically scored.  That decision used to be made twice, badly:

* the scanner took ``top_k = 8`` by mean accessibility -- with no diversity
  rule, so eight near-identical windows shifted by a nucleotide could fill the
  entire shortlist;
* ``pipeline._prescore()`` collapsed accessibility, Hamming and two pass/fail
  bonuses into one scalar and took the top 40.

This module replaces both with the scheme the team already uses in the NOT-gate
pipeline: score each candidate on independent criteria, peel Pareto fronts, and
take greedily subject to a maximum mutual overlap.

Criteria taken from the trigger-selection script
------------------------------------------------
These are the four scanner criteria that are genuine properties of the trigger
itself, and they keep the scanner's most important design decision: each is
converted to a 0-1 quality against a **fixed reference scale**, never min-maxed
against the rest of this run's pool.  A quality of 0.8 means the same thing
today as it will after changing the gene set or the shortlist size.

    (1) accessibility     mean P(unpaired) over the trigger footprint, measured
                          in gene context by stage 2 (already absolute)
    (2) openness          fraction of the functional region left unpaired when
                          the trigger folds alone (already absolute)
    (3) fold stability    p = exp(-(E_mfe - E_ens)/kT)          (already absolute)
    (5) sticking          longest unintended reverse-complement stretch between
                          Trigger A and Trigger B, connector masked out
    (6) substitution      longest unintended identical stretch, same masking

Criteria 4 (switch binding) and 7 (AND specificity) are NOT here.  Both were
computed against ``revcomp(r1 + r2 + a + k1)`` -- a switch in which r1 and r2 sit
next to each other, which is not the molecule this pipeline builds.  The
question they were trying to answer (does one trigger alone already fire the
gate?) is answered properly, and kinetically, in stage 5.

Masking note
------------
The scanner always masked the x/k2 connector before measuring cross-talk,
because in a split-trigger design that duplex is intended.  Here it is
parasitic: ``x = revcomp(k2)`` exists so k2* can pair x* inside the inhibitory
stem, and it makes the two triggers complementary as a side effect.  So the
default follows ``cfg.crosstalk_mask_connector`` (False), i.e. the liability is
measured rather than hidden -- while ``sticking``/``substitution`` still report
the masked value separately so the two effects stay distinguishable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import sequence_utils as su
from .config import PipelineConfig


# --------------------------------------------------------------------------- #
@dataclass
class CandidateQuality:
    """Per-criterion qualities in [0, 1], all against fixed reference scales."""
    accessibility: float = 0.0       # (1)
    openness: float = 0.0            # (2)
    fold_stability: float = 0.0      # (3)
    sticking: float = 0.0            # (5) 1 == no unintended complementarity
    substitution: float = 0.0        # (6) 1 == no unintended identity

    # raw values kept for reporting -- a quality with no raw number behind it
    # is not auditable
    raw: dict = field(default_factory=dict)

    def objectives(self) -> tuple[float, ...]:
        """The axes the Pareto front is computed on (all maximised).

        Sticking and substitution are combined into one "cross-talk" axis: they
        measure the same failure mode from two directions, and keeping them
        separate would inflate the front with candidates that are merely
        unusual on one of them.
        """
        return (self.accessibility, self.openness, self.fold_stability,
                min(self.sticking, self.substitution))

    def scalar(self, weights: dict | None = None) -> float:
        """Within-front tie-breaker only -- never the ranking axis itself."""
        w = weights or {"accessibility": 1.0, "openness": 1.0,
                        "fold_stability": 1.5, "crosstalk": 1.0}
        total = sum(w.values())
        return (w["accessibility"] * self.accessibility
                + w["openness"] * self.openness
                + w["fold_stability"] * self.fold_stability
                + w["crosstalk"] * min(self.sticking, self.substitution)) / total


# --------------------------------------------------------------------------- #
# criterion measurement                                                        #
# --------------------------------------------------------------------------- #
def _region_openness(structure: str, start: int, length: int) -> float:
    region = structure[start:start + length]
    return region.count(".") / len(region) if region else 0.0


def evaluate_quality(pair, tmA, tmB, backend, cfg: PipelineConfig
                     ) -> CandidateQuality:
    """Measure the four selection criteria for one trigger pair.

    ``tmA`` / ``tmB`` are stage-2 :class:`TriggerMetrics`; they already carry the
    in-context accessibility, so it is not recomputed here.
    """
    ta, tb = pair.triggerA, pair.triggerB
    seq_a, seq_b = ta.seq, tb.seq
    q = CandidateQuality()

    # (1) accessibility in native context -- straight from stage 2
    q.accessibility = 0.5 * (tmA.accessibility + tmB.accessibility)

    # (2) does the trigger hide its own functional region when folded alone?
    #     Trigger A's functional region is x + r1 (what pairs the switch's
    #     toehold); Trigger B's is r2 (its landing pad on r2*).
    struct_a, _ = backend.mfe(seq_a)
    struct_b, _ = backend.mfe(seq_b)
    a_fn_start = len(ta.k1) + len(ta.a)
    open_a = _region_openness(struct_a, a_fn_start, len(ta.x) + len(ta.r1))
    open_b = _region_openness(struct_b, len(tb.k2), len(tb.r2))
    q.openness = 0.5 * (open_a + open_b)

    # (3) how reliably does each trigger adopt its own MFE fold?
    p_a = backend.structure_probability(seq_a, cfg.temperature_c)
    p_b = backend.structure_probability(seq_b, cfg.temperature_c)
    q.fold_stability = 0.5 * (p_a + p_b)

    # (5)+(6) cross-talk between the two triggers
    if cfg.crosstalk_mask_connector:
        cmp_a = su.mask_region(seq_a, a_fn_start, len(ta.x))
        cmp_b = su.mask_region(seq_b, 0, len(tb.k2))
    else:
        cmp_a, cmp_b = seq_a, seq_b
    stick_nt = su.max_revcomp_match(cmp_a, cmp_b)
    subst_nt = su.max_identity_match(cmp_a, cmp_b)
    ref = max(1.0, float(cfg.unintended_match_nt_ref))
    q.sticking = max(0.0, 1.0 - stick_nt / ref)
    q.substitution = max(0.0, 1.0 - subst_nt / ref)

    # the same two numbers with the connector masked, for reporting only
    masked_a = su.mask_region(seq_a, a_fn_start, len(ta.x))
    masked_b = su.mask_region(seq_b, 0, len(tb.k2))

    q.raw = {
        "acc_A": tmA.accessibility, "acc_B": tmB.accessibility,
        "open_A": open_a, "open_B": open_b,
        "p_fold_A": p_a, "p_fold_B": p_b,
        "stick_nt": stick_nt, "subst_nt": subst_nt,
        "stick_nt_masked": su.max_revcomp_match(masked_a, masked_b),
        "subst_nt_masked": su.max_identity_match(masked_a, masked_b),
        "hamming": pair.hamming, "exact": pair.exact,
    }
    return q


# --------------------------------------------------------------------------- #
# Pareto fronts                                                                #
# --------------------------------------------------------------------------- #
def _dominates(p: tuple, q: tuple) -> bool:
    """p dominates q: at least as good on every axis, strictly better on one."""
    return all(a >= b for a, b in zip(p, q)) and any(a > b for a, b in zip(p, q))


def pareto_fronts(points: list[tuple]) -> list[int]:
    """Assign a 0-based front index to each point (all objectives maximised).

    Straightforward peeling.  The NOT-gate pipeline uses an O(n log n) sweep,
    but that only applies in two dimensions; here there are four axes and the
    pool reaching this stage is small (hundreds, not millions), so clarity wins.
    """
    n = len(points)
    front = [-1] * n
    remaining = set(range(n))
    level = 0
    while remaining:
        current = [i for i in remaining
                   if not any(_dominates(points[j], points[i])
                              for j in remaining if j != i)]
        if not current:                       # cycle guard; cannot normally happen
            current = list(remaining)
        for i in current:
            front[i] = level
        remaining -= set(current)
        level += 1
    return front


# --------------------------------------------------------------------------- #
# diversity                                                                    #
# --------------------------------------------------------------------------- #
def _overlap_fraction(w1: tuple[int, int], w2: tuple[int, int]) -> float:
    lo = max(w1[0], w2[0])
    hi = min(w1[1], w2[1])
    if hi <= lo:
        return 0.0
    shorter = min(w1[1] - w1[0], w2[1] - w2[0])
    return (hi - lo) / max(1, shorter)


def _too_similar(pair, chosen, max_overlap: float) -> bool:
    """True when ``pair`` reuses essentially the same windows as one already
    chosen.  Compared per source gene, so two pairs that share Trigger A but
    read completely different Trigger Bs are still kept."""
    for other in chosen:
        same_a = pair.meta.get("gene_a_name") == other.meta.get("gene_a_name")
        same_b = pair.meta.get("gene_b_name") == other.meta.get("gene_b_name")
        if not (same_a and same_b):
            continue
        ov_a = _overlap_fraction(pair.triggerA.window, other.triggerA.window)
        ov_b = _overlap_fraction(pair.triggerB.window, other.triggerB.window)
        if ov_a > max_overlap and ov_b > max_overlap:
            return True
    return False


# --------------------------------------------------------------------------- #
# the selection itself                                                         #
# --------------------------------------------------------------------------- #
def select(scored: list[tuple], cfg: PipelineConfig, k: int,
           max_overlap: float | None = None, progress=None) -> list[tuple]:
    """Pick ``k`` candidates from ``scored`` = [(pair, tmA, tmB, quality), ...].

    Peels Pareto fronts over the four criteria; inside each front sorts by the
    weighted scalar (tie-break only) and takes greedily while rejecting
    candidates that overlap an already-chosen one by more than ``max_overlap``.
    Continues into the next front if a front does not yield enough -- the same
    rule as ``choose_trigger.select_top_triggers`` in the NOT-gate pipeline.

    If the overlap rule cannot fill the quota (a small pool of genuinely similar
    windows), it relaxes rather than returning fewer than requested.
    """
    if not scored:
        return []
    max_overlap = (cfg.select_max_overlap if max_overlap is None
                   else max_overlap)
    fronts = pareto_fronts([q.objectives() for _p, _a, _b, q in scored])
    order = sorted(range(len(scored)),
                   key=lambda i: (fronts[i], -scored[i][3].scalar()))

    chosen: list[tuple] = []
    chosen_pairs: list = []
    for i in order:
        if len(chosen) >= k:
            break
        pair = scored[i][0]
        if _too_similar(pair, chosen_pairs, max_overlap):
            continue
        chosen.append(scored[i])
        chosen_pairs.append(pair)

    if len(chosen) < k:                       # relax the diversity rule
        picked = {id(c[0]) for c in chosen}
        for i in order:
            if len(chosen) >= k:
                break
            if id(scored[i][0]) not in picked:
                chosen.append(scored[i])
                picked.add(id(scored[i][0]))

    if progress:
        n_front0 = sum(1 for f in fronts if f == 0)
        progress(f"[select] {len(scored)} candidates -> {max(fronts) + 1} Pareto "
                 f"fronts ({n_front0} on front 1); kept {len(chosen)} "
                 f"(overlap cap {max_overlap:.0%})")
    return chosen
