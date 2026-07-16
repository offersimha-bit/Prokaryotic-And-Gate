"""Stage 1 -- target scanning and trigger definition.

For two genes G1 (source of Trigger A) and G2 (source of Trigger B) we look for
a length-``Lx`` segment ``x`` in G1 whose reverse complement occurs in G2 (that
occurrence is ``k2``).  Perfect reverse complements in natural genes are rare,
so when no exact match exists we fall back to the minimum Hamming-distance
window and design the switch against the *actual* trigger sequences (Section 6,
"Mismatch Handling").

Trigger A (5'->3'):  r1 | x | a | k1        (length L_A)
Trigger B (5'->3'):  r2 | k2                (length L_B)

The role swap (G1<->G2) is performed by the caller running :func:`scan_pair`
twice; :func:`scan_both_orientations` does both.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import sequence_utils as su
from .config import PipelineConfig


@dataclass
class TriggerA:
    gene: str
    pos_x: int              # 0-based start of x within the gene
    r1: str
    x: str
    a: str
    k1: str

    @property
    def seq(self) -> str:
        return self.r1 + self.x + self.a + self.k1


@dataclass
class TriggerB:
    gene: str
    pos_k2: int             # 0-based start of k2 within the gene
    r2: str
    k2: str

    @property
    def seq(self) -> str:
        return self.r2 + self.k2


@dataclass
class TriggerPair:
    orientation: str        # "G1->A,G2->B" or the swap
    gene_a: str
    gene_b: str
    triggerA: TriggerA
    triggerB: TriggerB
    hamming: int            # between reverse_complement(x) and k2
    exact: bool
    meta: dict = field(default_factory=dict)

    @property
    def hamming_fraction(self) -> float:
        return self.hamming / max(1, len(self.triggerA.x))


def _valid_x_positions(gene: str, cfg: PipelineConfig) -> range:
    """x positions in gene_a with room for r1 upstream and a+k1 downstream."""
    lr1 = cfg.resolved_len_r1()
    tail = cfg.Lx + cfg.len_a + cfg.len_k1
    lo = lr1
    hi = len(gene) - tail
    return range(lo, hi + 1) if hi >= lo else range(0)


def _valid_k2_positions(gene: str, cfg: PipelineConfig) -> range:
    """k2 positions in gene_b with room for r2 upstream."""
    lr2 = cfg.resolved_len_r2()
    lo = lr2
    hi = len(gene) - cfg.Lx
    return range(lo, hi + 1) if hi >= lo else range(0)


def _best_k2_match(rc_x: str, gene_b: str, positions: range):
    """Return (best_pos, best_hamming, exact_flag).

    First tries exact occurrences of ``rc_x`` (Hamming 0); if none fall inside
    the allowed ``positions``, scans for the minimum Hamming-distance window.
    """
    allowed = set(positions)
    # fast path: exact reverse-complement occurrences
    for idx in su.find_all(gene_b, rc_x):
        if idx in allowed:
            return idx, 0, True
    # slow path: minimum Hamming window
    best_pos, best_h = -1, len(rc_x) + 1
    Lx = len(rc_x)
    for j in positions:
        h = su.hamming(rc_x, gene_b[j:j + Lx])
        if h < best_h:
            best_pos, best_h = j, h
            if h == 0:
                break
    return best_pos, best_h, False


def _build_triggerA(gene: str, i: int, cfg: PipelineConfig) -> TriggerA:
    lr1 = cfg.resolved_len_r1()
    x_end = i + cfg.Lx
    a_end = x_end + cfg.len_a
    k1_end = a_end + cfg.len_k1
    return TriggerA(
        gene=gene, pos_x=i,
        r1=gene[i - lr1:i],
        x=gene[i:x_end],
        a=gene[x_end:a_end],
        k1=gene[a_end:k1_end],
    )


def _build_triggerB(gene: str, j: int, cfg: PipelineConfig) -> TriggerB:
    lr2 = cfg.resolved_len_r2()
    return TriggerB(
        gene=gene, pos_k2=j,
        r2=gene[j - lr2:j],
        k2=gene[j:j + cfg.Lx],
    )


def scan_pair(gene_a: str, gene_b: str, cfg: PipelineConfig,
              orientation: str = "G1->A,G2->B",
              max_candidates: int | None = None) -> list[TriggerPair]:
    """Find Trigger A/B pairs with gene_a supplying A and gene_b supplying B."""
    ga = su.to_rna(gene_a)
    gb = su.to_rna(gene_b)
    xs = _valid_x_positions(ga, cfg)
    k2_positions = _valid_k2_positions(gb, cfg)
    max_h = int(cfg.max_hamming_fraction * cfg.Lx)

    pairs: list[TriggerPair] = []
    for i in xs:
        x = ga[i:i + cfg.Lx]
        if "N" in x:
            continue
        rc_x = su.reverse_complement(x)
        j, h, exact = _best_k2_match(rc_x, gb, k2_positions)
        if j < 0 or h > max_h:
            continue
        pairs.append(TriggerPair(
            orientation=orientation,
            gene_a=gene_a, gene_b=gene_b,
            triggerA=_build_triggerA(ga, i, cfg),
            triggerB=_build_triggerB(gb, j, cfg),
            hamming=h, exact=exact,
        ))
    pairs.sort(key=lambda p: (p.hamming, p.triggerA.pos_x))
    if max_candidates is not None:
        pairs = pairs[:max_candidates]
    return pairs


def scan_both_orientations(gene1: str, gene2: str, cfg: PipelineConfig,
                           max_candidates: int | None = None
                           ) -> list[TriggerPair]:
    """Run the scan and its role-swapped counterpart (Section 1, "Iteration")."""
    out = scan_pair(gene1, gene2, cfg, "G1->A,G2->B", max_candidates)
    out += scan_pair(gene2, gene1, cfg, "G2->A,G1->B", max_candidates)
    return out
