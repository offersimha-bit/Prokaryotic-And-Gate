"""STAGE 7 -- rank the scored designs without inventing weights.

A design is good along two axes that genuinely trade off:

    P_11           absolute ON-state output      -- how much protein
    logic_margin   ON over the worst OFF state   -- how selective

Collapsing them into one number requires a weight, and a weight here would be a
guess.  Non-dominated (Pareto) sorting avoids that: a design belongs to front 1
if no other design beats it on BOTH axes.  Front 2 is what remains after front 1
is removed, and so on, so EVERY design gets a rank -- not just the optimal set --
and the user picks where on the front to sit rather than accepting our
preference.

A scalar is still offered for convenience, but only as an explicit user
preference (``pick``), never as a hidden default.

    python -m and_gate_pipeline.rank
"""

from __future__ import annotations

from dataclasses import dataclass, field


def _dominates(p: tuple, q: tuple) -> bool:
    """True when p is at least as good as q everywhere and better somewhere."""
    return all(a >= b for a, b in zip(p, q)) and any(a > b for a, b in zip(p, q))


def pareto_fronts(points: list[tuple]) -> list[int]:
    """Front index (0 = optimal) for every point.  All objectives maximised.

    O(n^2) in the general case, which is fine at pipeline scale: stage 7 sees
    tens of designs, not thousands, because stages 2-3 have already cut the
    field.  The 2-D O(n log n) sweep is not used because ``rank_axes`` is
    configurable and may hold more than two objectives.
    """
    n = len(points)
    fronts = [-1] * n
    remaining = set(range(n))
    current = 0
    while remaining:
        layer = [i for i in remaining
                 if not any(_dominates(points[j], points[i])
                            for j in remaining if j != i)]
        if not layer:                      # mutually dominating duplicates
            layer = list(remaining)
        for i in layer:
            fronts[i] = current
        remaining -= set(layer)
        current += 1
    return fronts


def _normalise(points: list[tuple]) -> list[tuple]:
    """Min-max each axis onto [0, 1] so axes of different units compare.

    Used ONLY for geometry (knee, preference blend), never for the front itself:
    domination is scale-invariant, so the front never depends on this.
    """
    if not points:
        return []
    cols = list(zip(*points))
    out = []
    for col in cols:
        lo, hi = min(col), max(col)
        span = hi - lo
        out.append([0.5 if span == 0 else (v - lo) / span for v in col])
    return list(zip(*out))


def knee_index(points: list[tuple], front: list[int]) -> int | None:
    """Design on front 0 closest to the ideal corner -- the best compromise.

    Reported, not imposed: it is one defensible pick among the front, and the
    front is the real answer.
    """
    idx = [i for i, f in enumerate(front) if f == 0]
    if not idx:
        return None
    norm = _normalise(points)
    best, best_d = None, float("inf")
    for i in idx:
        d = sum((1.0 - v) ** 2 for v in norm[i])
        if d < best_d:
            best, best_d = i, d
    return best


def pick(points: list[tuple], front: list[int], preference: float = 0.5) -> int | None:
    """Choose one design from front 0 by an EXPLICIT preference in [0, 1].

    0.0 = maximise the last axis (selectivity), 1.0 = maximise the first
    (output).  There is no default worth defending, which is why the caller has
    to name it; 0.5 simply means "no preference stated".
    """
    idx = [i for i, f in enumerate(front) if f == 0]
    if not idx:
        return None
    norm = _normalise(points)
    w = max(0.0, min(1.0, preference))
    best, best_s = None, -float("inf")
    for i in idx:
        v = norm[i]
        s = w * v[0] + (1.0 - w) * v[-1]
        if s > best_s:
            best, best_s = i, s
    return best


@dataclass
class Ranking:
    axes: tuple
    points: list[tuple]
    fronts: list[int]
    order: list[int]                        # indices, best first
    knee: int | None = None
    notes: list[str] = field(default_factory=list)

    def front_sizes(self) -> dict:
        out: dict[int, int] = {}
        for f in self.fronts:
            out[f] = out.get(f, 0) + 1
        return dict(sorted(out.items()))


def rank_designs(results, cfg, preference: float | None = None) -> Ranking:
    """Assign Pareto fronts to scored designs and order them.

    ``results`` are DesignResult objects whose ``score.metrics`` holds the
    four-state output.  Within a front, ties are broken by the first axis --
    an arbitrary but declared choice, since the whole point of a front is that
    its members are not comparable.
    """
    axes = tuple(getattr(cfg, "rank_axes", ("P_11", "logic_margin")))
    points = []
    for r in results:
        m = getattr(r.score, "metrics", None)
        if m is None:                       # legacy ScoreCard -- one axis only
            points.append((float(getattr(r.score, "total", 0.0)),))
        else:
            points.append(tuple(float(m.get(a, 0.0)) for a in axes))

    fronts = pareto_fronts(points)
    order = sorted(range(len(points)), key=lambda i: (fronts[i], -points[i][0]))
    kn = knee_index(points, fronts)
    rk = Ranking(axes=axes, points=points, fronts=fronts, order=order, knee=kn)

    n_floored = sum(1 for r in results
                    if getattr(r.score, "metrics", {}).get("leak_floor_active"))
    if n_floored:
        rk.notes.append(
            f"leak floor active on {n_floored}/{len(results)} designs: the "
            f"modelled OFF state fell below cfg.leak_floor, so their ratios are "
            f"capped by the floor rather than by the design")

    # A front per design means the objectives induce a TOTAL order -- nothing
    # trades off, so the Pareto machinery has added no information and the
    # ranking is really a sort on one axis.  Say so rather than let a degenerate
    # front be mistaken for a genuine trade-off analysis.
    if len(points) > 2 and len(set(fronts)) == len(fronts):
        why = ""
        if n_floored == len(results) and len(axes) == 2:
            why = (" -- every OFF state is below the floor, so logic_margin is "
                   "exactly P_11/leak_floor and the two axes are proportional. "
                   "Selectivity is not the limiting factor for these designs; "
                   "output is")
        rk.notes.append("axes induce a total order: no design trades one "
                        "objective for another" + why)
    if preference is not None:
        p = pick(points, fronts, preference)
        if p is not None:
            rk.order = [p] + [i for i in order if i != p]
            rk.notes.append(f"ordered by explicit preference={preference:.2f}")
    return rk


def report(rk: Ranking, results, top: int = 10):
    print(f"Pareto ranking over {rk.axes} -- both maximised, no weights")
    print(f"fronts: {rk.front_sizes()}   (front 0 = non-dominated)")
    for n in rk.notes:
        print(f"note: {n}")
    print()
    hdr = "%-5s %-6s " % ("rank", "front") + " ".join("%14s" % a for a in rk.axes)
    print(hdr); print("-" * len(hdr))
    for rank, i in enumerate(rk.order[:top], 1):
        mark = "  <- knee" if i == rk.knee else ""
        print("%-5d %-6d " % (rank, rk.fronts[i])
              + " ".join("%14.6g" % v for v in rk.points[i]) + mark)
