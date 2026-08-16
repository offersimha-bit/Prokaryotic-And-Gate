"""STAGE 1 -- target scanning and trigger definition.

This module is the merge of two previously separate implementations:

* ``and_gate_pipeline/target_scan.py``   -- exact two-gene scan, min-Hamming
  fallback, both role orientations, correct ``k1-a-x-r1`` domain order.
* ``Triger finding/and_gate_trigger.py`` -- pooled multi-gene discovery, the
  hash-join on the connector, the FASTA reader, and the opt-in window filters.

What was taken from the trigger-selection script
------------------------------------------------
``read_fasta_records``      pooled FASTA input, one or many records per file
``load_genes``              any number of files, all records pooled
hash-join on the connector  ``find_matching_pairs``'s ``defaultdict`` idea, here
                            generalised to a k-mer index so BOTH the pooled scan
                            and the two-gene scan are O(N) instead of O(N*M*Lx)
two-distinct-records rule   a real AND gate needs two inputs, not one gene
                            sensing itself twice
opt-in window filters       GC ceiling and forbidden motifs, now **off by
                            default** and applied to the switch, not the trigger
                            (see the note on Type IIS sites below)

What was deliberately NOT taken
-------------------------------
``build_switch_target``     ``revcomp(r1 + r2 + a + k1)`` assumes r1 and r2 are
                            adjacent on one contiguous binding site.  In this
                            architecture r2 binds ``r2*`` on the inhibitory
                            hairpin while r1 binds ``r1*`` inside the primary
                            toehold, separated by ``k2*``, the secondary loop
                            and ``x*``.  Criteria 4 and 7 of the scanner score
                            that non-existent molecule and are dropped; the real
                            question they were asking is answered kinetically in
                            stage 5.
its accessibility measure   ``unpaired_probs(whole_gene)`` is a global fold of a
                            full transcript.  Accessibility is stage 2's job and
                            is done there with flanked windows (and should move
                            to ``RNA.pfl_fold_up`` -- see 02_filtering.py).
Type IIS on the trigger     the trigger is a natural sequence that is never
                            synthesised, so Golden Gate sites in it are not a
                            constructability problem.  That filter belongs to
                            the switch.

Domain order (the bug this merge closes)
----------------------------------------
The scanner sliced its window as ``r1 | x | a | k1`` in genomic order, while the
switch requires Trigger A to read ``k1 - a - x - r1`` 5'->3'.  The old
``interop.windows_to_pair()`` handed the scanner's slices straight to
``TriggerA``, whose ``.seq`` reassembled them in the other order -- producing a
sequence that does not occur in the gene, and a ``pos_x`` that made stage 2 fold
the wrong window.  There is now exactly one slicer (:func:`_build_triggerA`),
one coordinate convention, and :func:`verify_pair` asserts that both triggers
are literally substrings of their source genes.

    Trigger A (5'->3'):  k1 | a | x | r1     (length L_A)   -- k1+a UPSTREAM of x
    Trigger B (5'->3'):  k2 | r2             (length L_B)   -- r2 DOWNSTREAM of k2
"""

from __future__ import annotations

import os
from collections import defaultdict
from dataclasses import dataclass, field

from . import sequence_utils as su
from .config import PipelineConfig


# --------------------------------------------------------------------------- #
# trigger records                                                             #
# --------------------------------------------------------------------------- #
@dataclass
class TriggerA:
    """Primary trigger: one contiguous window of gene A, read ``k1-a-x-r1``.

    The order is forced by the switch.  Trigger A's binding site runs
    ``r1* | x* | a* | k1*`` 5'->3' along the transcript, and a contiguous
    antiparallel duplex requires the trigger to be its reverse complement, i.e.
    k1 first.  Building it the other way round gave dG(A:switch) = -4.3 instead
    of -30.8 kcal/mol.
    """
    gene: str
    pos_x: int              # 0-based start of x within the gene
    r1: str
    x: str
    a: str
    k1: str

    @property
    def seq(self) -> str:
        return self.k1 + self.a + self.x + self.r1

    @property
    def window(self) -> tuple[int, int]:
        """[start, end) of the whole trigger inside its gene."""
        start = self.pos_x - len(self.a) - len(self.k1)
        return start, start + len(self.seq)


@dataclass
class TriggerB:
    """Secondary trigger: one contiguous window of gene B, read ``k2-r2``.

    Binds the switch's ``r2* | k2*`` (5'->3'): r2 lands on the single-stranded
    r2* toehold, then k2 invades the k2*:x* helix and opens the inhibitory
    hairpin.
    """
    gene: str
    pos_k2: int             # 0-based start of k2 within the gene
    r2: str
    k2: str

    @property
    def seq(self) -> str:
        return self.k2 + self.r2

    @property
    def window(self) -> tuple[int, int]:
        return self.pos_k2, self.pos_k2 + len(self.seq)


@dataclass
class TriggerPair:
    orientation: str
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


class TriggerIntegrityError(ValueError):
    """Raised when a trigger is not a contiguous slice of its source gene.

    This is the invariant the whole project rests on: nothing is synthesised,
    both triggers are real stretches of real genes.  If it breaks, every
    downstream number describes a molecule that cannot exist.
    """


def verify_pair(pair: TriggerPair) -> TriggerPair:
    """Assert both triggers really occur, contiguously, in their source genes.

    Checks the sequence *and* the recorded coordinates, because stage 2 folds
    the trigger by coordinate, not by sequence -- a correct sequence at a wrong
    offset would silently fold the wrong window.
    """
    for name, trig, gene in (("A", pair.triggerA, pair.gene_a),
                             ("B", pair.triggerB, pair.gene_b)):
        gene = su.to_rna(gene)
        s, e = trig.window
        if s < 0 or e > len(gene):
            raise TriggerIntegrityError(
                f"Trigger {name}: window [{s}, {e}) falls outside a gene of "
                f"length {len(gene)}")
        if gene[s:e] != trig.seq:
            raise TriggerIntegrityError(
                f"Trigger {name} is not a contiguous slice of its gene.\n"
                f"  gene[{s}:{e}] = {gene[s:e]}\n"
                f"  trigger.seq   = {trig.seq}\n"
                "The domain order of the slicer and of the .seq property have "
                "diverged -- see this module's docstring.")
    return pair


# --------------------------------------------------------------------------- #
# pooled FASTA input (from the trigger-selection script)                      #
# --------------------------------------------------------------------------- #
def read_fasta_records(path: str) -> list[tuple[str, str]]:
    """Parse a FASTA file that may hold one or many gene records.

    Returns [(name, sequence), ...] in file order.  Taken from the scanner so
    input behaviour is identical to what the team is already used to.
    """
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].strip() or f"record{len(records) + 1}"
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def load_genes(paths) -> list[tuple[str, str]]:
    """Pool every record from one path, a list of paths, or a directory."""
    if isinstance(paths, str):
        paths = [paths]
    out: list[tuple[str, str]] = []
    for p in paths:
        if os.path.isdir(p):
            for fn in sorted(os.listdir(p)):
                if fn.lower().endswith((".fa", ".fasta", ".fna")):
                    out.extend(read_fasta_records(os.path.join(p, fn)))
        else:
            out.extend(read_fasta_records(p))
    return out


# --------------------------------------------------------------------------- #
# window filters (opt-in; from the scanner, defaults changed)                  #
# --------------------------------------------------------------------------- #
def _window_rejected(window: str, cfg: PipelineConfig) -> str | None:
    """Return a reason string when a candidate window should be skipped.

    Both filters are OFF by default.  ``scan_max_gc`` is a real biological
    preference (GC-rich windows fold hard); ``scan_forbid_motifs`` exists only
    so the old scanner behaviour is reproducible -- it should normally stay off,
    because the trigger is endogenous and never synthesised.
    """
    if "N" in window:
        return "ambiguous base"
    max_gc = getattr(cfg, "scan_max_gc", 0.0)
    if max_gc and su.gc_fraction(window) > max_gc:
        return "GC above scan_max_gc"
    if getattr(cfg, "scan_forbid_motifs", False):
        rna = su.to_rna(window)
        for m in list(cfg.forbidden_runs) + list(cfg.forbidden_motifs):
            if su.to_rna(m) in rna:
                return f"contains {m}"
    return None


# --------------------------------------------------------------------------- #
# slicing                                                                      #
# --------------------------------------------------------------------------- #
def _valid_x_positions(gene: str, cfg: PipelineConfig) -> range:
    """x positions with room for k1+a UPSTREAM and r1 DOWNSTREAM."""
    lo = cfg.len_k1 + cfg.len_a
    hi = len(gene) - cfg.Lx - cfg.resolved_len_r1()
    return range(lo, hi + 1) if hi >= lo else range(0)


def _valid_k2_positions(gene: str, cfg: PipelineConfig) -> range:
    """k2 positions with room for r2 DOWNSTREAM."""
    hi = len(gene) - cfg.Lx - cfg.resolved_len_r2()
    return range(0, hi + 1) if hi >= 0 else range(0)


def _build_triggerA(gene: str, i: int, cfg: PipelineConfig) -> TriggerA:
    lr1 = cfg.resolved_len_r1()
    k1_start = i - cfg.len_a - cfg.len_k1
    a_start = i - cfg.len_a
    x_end = i + cfg.Lx
    return TriggerA(gene=gene, pos_x=i,
                    k1=gene[k1_start:a_start],
                    a=gene[a_start:i],
                    x=gene[i:x_end],
                    r1=gene[x_end:x_end + lr1])


def _build_triggerB(gene: str, j: int, cfg: PipelineConfig) -> TriggerB:
    lr2 = cfg.resolved_len_r2()
    k2_end = j + cfg.Lx
    return TriggerB(gene=gene, pos_k2=j,
                    k2=gene[j:k2_end],
                    r2=gene[k2_end:k2_end + lr2])


# --------------------------------------------------------------------------- #
# the k-mer index -- the scanner's hash-join, generalised                     #
# --------------------------------------------------------------------------- #
def build_k2_index(gene: str, cfg: PipelineConfig) -> dict[str, list[int]]:
    """{k2 sequence -> [valid start positions]} for one gene.

    This is the scanner's ``defaultdict`` join, moved onto the k2 side so a
    single pass over gene A can look every connector up in O(1).  It replaces
    ``_best_k2_match``'s O(|G_A| * |G_B| * Lx) Hamming sweep.
    """
    idx: dict[str, list[int]] = defaultdict(list)
    for j in _valid_k2_positions(gene, cfg):
        idx[gene[j:j + cfg.Lx]].append(j)
    return idx


def _neighbours(seq: str, max_h: int):
    """All strings within Hamming distance ``max_h`` of ``seq``.

    Only used for the mismatch-tolerant fallback.  Size is
    C(L, h) * 3^h -- at Lx=6 and max_h=2 that is 135 lookups, far cheaper than
    sweeping the whole gene, and it degrades gracefully because ``max_h`` is
    bounded by ``cfg.max_hamming_fraction``.
    """
    yield seq, 0
    if max_h <= 0:
        return
    alphabet = "ACGU"
    n = len(seq)
    seen = {seq}
    frontier = [seq]
    for dist in range(1, max_h + 1):
        nxt = []
        for s in frontier:
            for i in range(n):
                for b in alphabet:
                    if b == s[i]:
                        continue
                    cand = s[:i] + b + s[i + 1:]
                    if cand not in seen:
                        seen.add(cand)
                        nxt.append(cand)
                        yield cand, dist
        frontier = nxt


def _lookup_k2(rc_x: str, index: dict, max_h: int):
    """Best (position, hamming, exact) for ``rc_x`` in a k2 index."""
    for cand, dist in _neighbours(rc_x, max_h):
        hits = index.get(cand)
        if hits:
            return hits[0], dist, dist == 0
    return -1, max_h + 1, False


# --------------------------------------------------------------------------- #
# two-gene scan                                                                #
# --------------------------------------------------------------------------- #
def scan_pair(gene_a: str, gene_b: str, cfg: PipelineConfig,
              orientation: str = "G1->A,G2->B",
              max_candidates: int | None = None,
              name_a: str = "", name_b: str = "") -> list[TriggerPair]:
    """Trigger pairs with gene_a supplying Trigger A and gene_b supplying B."""
    ga, gb = su.to_rna(gene_a), su.to_rna(gene_b)
    index = build_k2_index(gb, cfg)
    max_h = int(cfg.max_hamming_fraction * cfg.Lx)

    pairs: list[TriggerPair] = []
    for i in _valid_x_positions(ga, cfg):
        ta = _build_triggerA(ga, i, cfg)
        reason = _window_rejected(ta.seq, cfg)
        if reason:
            continue
        j, h, exact = _lookup_k2(su.reverse_complement(ta.x), index, max_h)
        if j < 0:
            continue
        tb = _build_triggerB(gb, j, cfg)
        if _window_rejected(tb.seq, cfg):
            continue
        pairs.append(TriggerPair(
            orientation=orientation, gene_a=ga, gene_b=gb,
            triggerA=ta, triggerB=tb, hamming=h, exact=exact,
            meta={"gene_a_name": name_a or "geneA",
                  "gene_b_name": name_b or "geneB"}))
    pairs.sort(key=lambda p: (p.hamming, p.triggerA.pos_x))
    if max_candidates is not None:
        pairs = pairs[:max_candidates]
    return [verify_pair(p) for p in pairs]


def scan_both_orientations(gene1: str, gene2: str, cfg: PipelineConfig,
                           max_candidates: int | None = None
                           ) -> list[TriggerPair]:
    """Run the scan and its role-swapped counterpart."""
    out = scan_pair(gene1, gene2, cfg, "G1->A,G2->B", max_candidates,
                    name_a="gene1", name_b="gene2")
    out += scan_pair(gene2, gene1, cfg, "G2->A,G1->B", max_candidates,
                     name_a="gene2", name_b="gene1")
    return out


# --------------------------------------------------------------------------- #
# pooled multi-gene discovery (the scanner's job, done here)                   #
# --------------------------------------------------------------------------- #
def scan_pool(genes: list[tuple[str, str]], cfg: PipelineConfig,
              max_pairs: int | None = None,
              require_exact: bool = False,
              max_per_gene_pair: int | None = None,
              progress=None) -> list[TriggerPair]:
    """Pool any number of gene records and find every workable trigger pair.

    ``genes`` is [(name, sequence), ...] -- one FASTA with many records, several
    files pooled, or a mix, exactly as the standalone scanner accepted.

    Every gene is indexed once on the k2 side, then every gene is swept once on
    the x side against all *other* genes' indices.  A pair is only kept when its
    two windows come from two different records: an AND gate needs two inputs,
    not one gene sensing itself twice.

    ``require_exact=True`` reproduces the scanner's hard rule that the connector
    must be a perfect reverse complement.  The default allows the min-Hamming
    fallback, whose cost then shows up honestly as mismatches inside the
    secondary stem.
    """
    prepared = [(name, su.to_rna(seq)) for name, seq in genes]
    indices = {name: build_k2_index(seq, cfg) for name, seq in prepared}
    max_h = 0 if require_exact else int(cfg.max_hamming_fraction * cfg.Lx)

    pairs: list[TriggerPair] = []
    for name_a, ga in prepared:
        for i in _valid_x_positions(ga, cfg):
            ta = _build_triggerA(ga, i, cfg)
            if _window_rejected(ta.seq, cfg):
                continue
            rc_x = su.reverse_complement(ta.x)
            per_pair: dict[str, int] = {}
            for name_b, gb in prepared:
                if name_b == name_a:
                    continue                      # two DIFFERENT gene records
                if max_per_gene_pair and per_pair.get(name_b, 0) >= max_per_gene_pair:
                    continue
                j, h, exact = _lookup_k2(rc_x, indices[name_b], max_h)
                if j < 0:
                    continue
                tb = _build_triggerB(gb, j, cfg)
                if _window_rejected(tb.seq, cfg):
                    continue
                per_pair[name_b] = per_pair.get(name_b, 0) + 1
                pairs.append(TriggerPair(
                    orientation=f"pool:{name_a}->A,{name_b}->B",
                    gene_a=ga, gene_b=gb, triggerA=ta, triggerB=tb,
                    hamming=h, exact=exact,
                    meta={"gene_a_name": name_a, "gene_b_name": name_b,
                          "pos_A": ta.window[0], "pos_B": tb.window[0]}))

    pairs.sort(key=lambda p: (p.hamming, p.meta.get("gene_a_name", ""),
                              p.triggerA.pos_x))
    if progress:
        n_exact = sum(p.exact for p in pairs)
        progress(f"[scan] {len(prepared)} gene records -> {len(pairs)} trigger "
                 f"pairs ({n_exact} exact connector matches)")
    if max_pairs is not None:
        pairs = pairs[:max_pairs]
    return [verify_pair(p) for p in pairs]


def scan_from_fasta(paths, cfg: PipelineConfig | None = None, **kw
                    ) -> list[TriggerPair]:
    """Convenience entry point: FASTA path(s) or a folder -> trigger pairs."""
    cfg = cfg or PipelineConfig()
    return scan_pool(load_genes(paths), cfg, **kw)
