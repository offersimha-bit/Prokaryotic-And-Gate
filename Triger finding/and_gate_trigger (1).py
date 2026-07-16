#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
split_and_gate_trigger_designer.py
==================================

Given genes in two FASTA files (each file may contain one or more gene
records), this searches for a pair of ENDOGENOUS split triggers - one found
inside a gene from file 1, one found inside a gene from file 2 - that could
work together as a split AND-gate toehold switch, following this layout:

    Trigger 1 = [R1] + ['*'] + [A] + [K1]      (one contiguous window in a gene from file 1)
    Trigger 2 = [R2] + [K2]                    (one contiguous window in a gene from file 2)

    * R1  -> the target-recognition part of Trigger 1 (should be accessible/unfolded on its own gene).
    * R2  -> the target-recognition part of Trigger 2 (should be accessible/unfolded on its own gene).
    * '*' -> the connector part of Trigger 1. Length = `star_len`.
    * K2  -> the connector part of Trigger 2. Length = `star_len` (must equal len('*')).
    * A   -> the toehold ('a' domain) of Trigger 1. Length = `A_len`.
    * K1  -> the stem key of Trigger 1's side of the switch. Length = `K1_len`.

Unlike a typical trigger-RNA design tool (e.g. MODesign, Pelea et al., eLife 2025), NOTHING here
is synthesized: R1/'*'/A/K1 and R2/K2 are each a single contiguous, naturally-occurring stretch of
sequence, sliced directly out of a real gene in order. Nothing is invented. The only way the two
halves can plausibly work together is if the '*' found in Trigger 1 happens to be the exact
reverse-complement of the K2 found in Trigger 2 - so that is treated as a hard, exact requirement:
every Trigger-1 window and every Trigger-2 window are compared, and only pairs where
'*' == reverse_complement(K2) are considered at all. This is a search for a coincidence in existing
biology, not a design task.

------------------------------------------------------------------------------
Once a '*'/K2-matching pair is found, it's scored on 7 criteria (same spirit as MODesign's
scoring, adapted for two independently-found half-triggers instead of one designed trigger):
    1. Target accessibility     -> R1/R2 should be regions that are open/unfolded on their own gene.
    2. Trigger self-folding     -> neither trigger should fold on itself and hide the parts it
                                    needs (R + toehold) - checked for *both* triggers.
    3. Structure stability      -> the Boltzmann probability that each trigger folds into the
                                    structure it's given (p = exp(-(E_struct - E_ensemble)/kT));
                                    higher is better => less chance of an "unwanted" fold.
    4. Switch binding strength  -> how strongly each trigger *alone* binds the toehold+stem of the
                                    switch (duplex energy). See criterion 7 for why "alone" matters.
    5. Cross-binding (sticking) -> whether the two triggers are complementary to each other
                                    (reverse-complement) beyond the intended '*'/K2 connector => unwanted binding.
    6. Cross-similarity (substitution) -> whether the two triggers share long identical
                                    stretches => could confuse the switch or affect its structure.
    (5) and (6) are checked using a multi-window scan, similar to `has_local_complementarity`
    in MODesign, but here both checks are done: reverse-complement (sticking) and identity (substitution),
    while always masking out the intended '*'/K2 connector so it isn't flagged as a problem.
    7. AND-gate specificity (anti-leak) -> per Kim & Moon (ACS Synth Biol, 2020) "Modulating
                                    responses of toehold switches by an inhibitory hairpin": a
                                    two-hairpin/two-trigger AND gate only behaves like a true AND
                                    gate when a single trigger alone cannot already bind the switch
                                    strongly - otherwise it leaks like a one-input switch (their
                                    Figure 1C/1D). This checks that binding the switch with only
                                    Trigger 1 or only Trigger 2 (criterion 4) is much weaker than
                                    the fully combined R1+R2+A+K1 complex (bind_full).

Each criterion is first converted to a 0-1 "quality" against a FIXED reference scale (not against
the other candidates found in this run), then combined into a 0-100 score (see evaluate_pair). This
matters: three of the criteria are already true, absolute probabilities/percentages by construction
(1: unpaired probability, 2: % unpaired positions, 3: Boltzmann probability), so their quality IS
the raw value - no reference choice needed. The other four (4, 5, 6, 7) are physical quantities
with no natural 0-1 scale (kcal/mol energies, nt counts), so converting them to a quality requires
picking a reference scale - done here with named constants (Params.dg_per_bp_ref,
Params.unintended_match_nt_ref) documented at their definition, not derived from a specific paper.
Earlier versions of this script instead min-max-normalized every criterion against only the other
candidates found in the same run - with a small candidate pool (e.g. top_k=8) this routinely puts
several criteria at exactly 0% or 100% by construction, which is statistically meaningless and not
comparable between runs/gene pairs/parameter choices. The current approach fixes that: a score of
80/100 means the same thing whether you ran this today or after changing top_k or the input genes.

Optional hard filters (all opt-in, disabled by default; see Params) mirror the ones MODesign
itself uses as outright pass/fail cutoffs rather than soft weights:
    - min_functional_openness = 80.0    would match MODesign's 80% openness cutoff.
    - min_structure_probability = 0.4   would match MODesign's fold-probability cutoff.
    - max_gc_content                    caps the GC fraction of an R1/R2 candidate window.
    - max_leak_fraction = 0.5           rejects candidates where a single trigger alone binds the
                                        switch at more than half the strength of the full complex.

Usage:
    python split_and_gate_trigger_designer.py                       # see input priority below
    python split_and_gate_trigger_designer.py g1.fa g2.fa            # 2+ FASTA files, any number of them
    python split_and_gate_trigger_designer.py g1.fa g2.fa g3.fa ...  # each file may hold ONE gene or MANY
                                                                       # gene records. ALL records from ALL
                                                                       # files given are pooled together, then
                                                                       # every Trigger-1 window found in the
                                                                       # pool is compared against every
                                                                       # Trigger-2 window in the pool - the
                                                                       # two windows of a pair are required to
                                                                       # come from two DIFFERENT gene records
                                                                       # (so e.g. one file with 5 genes works
                                                                       # too: the script picks whichever 2 of
                                                                       # those 5 genes gives the best pair).

With no command-line arguments, input is picked in this order:
    1. A folder named "genes" next to this script - drop any number of .fa/.fasta files there
       (each may hold one or many gene records) and every one of them is pooled automatically.
       This is the easiest way to add/remove genes without touching the command line at all.
    2. gene1.fa and gene2.fa in the current directory, if both exist.
    3. A small built-in demo pair.

For real results it's recommended to install one of:
    pip install ViennaRNA        # a Python package named RNA - recommended, easy
    # or NUPACK 4 (register first at https://nupack.org)  - what MODesign itself uses
If neither is installed, a simple FALLBACK is used instead (not very accurate; does not fold RNA properly).
"""

from __future__ import annotations

import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Tuple


# ======================================================================
# 1) Parameters - everything tunable, with numbers/switches matching the paper
# ======================================================================
@dataclass
class Params:
    # --- domain lengths (nt) ---
    R1_len: int = 11        # length of Trigger 1's target-recognition region
    R2_len: int = 30        # length of Trigger 2's target-recognition region
    star_len: int = 6        # length of the connector '*' in Trigger 1 (== length of K2 in Trigger 2)
    A_len: int = 4           # length of Trigger 1's toehold 'a'
    K1_len: int = 15         # length of Trigger 1's stem key K1

    # --- scanning ---
    scan_step: int = 1
    top_k: int = 8          # how many '*'/K2-matching pairs to fully fold & score (most-accessible first)

    # --- folding conditions ---
    temperature_c: float = 37.0   # (MODesign uses 30C; can be changed here)
    material: str = "rna"

    # --- weights for scoring (higher => more preferred) ---
    # Final combined score is reported on a fixed 0-100 reference scale (see rank_candidates),
    # so these weights only set each criterion's *relative* share of that 100 points.
    w_access: float = 1.0        # (1) target accessibility
    w_openness: float = 1.0      # (2) trigger self-folding openness
    w_robustness: float = 1.5    # (3) structure stability (p) - MODesign's own composite score
                                  # uses P^3 (cubed), i.e. it weighs fold-reliability much more
                                  # heavily than a linear term; raised from 1.0 to reflect that.
    w_binding: float = 1.0       # (4) switch binding strength
    w_stick: float = 1.0         # (5) cross-binding (reverse-complement)
    w_substitute: float = 1.0    # (6) cross-similarity (identity)
    w_specificity: float = 1.5   # (7) AND-gate specificity / anti-leak (Kim & Moon, ACS Synth Biol 2020):
                                  # a single trigger alone must bind the switch much more weakly than
                                  # both triggers together, or the "AND gate" leaks like a one-input switch.

    # --- window sizes for scanning/reporting similarity between triggers (as in MODesign) ---
    sim_report_windows: Tuple[int, ...] = (4, 6, 8, 10, 12)

    # --- hard filter thresholds (opt-in; 0 => disabled) ---
    min_functional_openness: float = 0.0     # %: minimum required openness of the trigger's functional region
                                              # (MODesign itself hard-discards below 80% - set to 80.0 to match)
    min_structure_probability: float = 0.0   # 0-1: minimum required fold-stability probability p
                                              # (MODesign hard-discards below 0.4 - set to 0.4 to match)
    max_gc_content: float = 0.0              # 0-1: maximum allowed GC fraction of an R1/R2 candidate window
                                              # (0 => disabled)
    max_unintended_match: int = 0            # nt: maximum allowed unintended cross-match length
    max_leak_fraction: float = 0.0           # 0-1: max allowed (single-trigger binding / full binding)
                                              # strength ratio to the switch (0 => disabled). E.g. 0.5
                                              # would reject candidates where either trigger alone binds
                                              # the switch at more than half the strength of both together.

    # --- forbidden motifs ---
    # AAAA/CCCC/GGGG/UUUU-length homopolymers: generic synthesis/secondary-structure hazards (this
    #   threshold matches Green et al. 2014's own NUPACK design constraints, Cell 159:925-939 SI).
    # UUUUU specifically also approximates a Rho-independent (intrinsic) transcription terminator's
    #   U-tract in E. coli - though the real signal is a GC-rich hairpin immediately followed by the
    #   U-tract, not the U-tract alone; this plain substring check is a cheap proxy, not a full
    #   terminator predictor.
    # GGUCUC/GAGACC (BsaI), CGUCUC/GAGACG (BsmBI/Esp3I), GCUCUUC/GAAGAGC (SapI), GAAGAC/GUCUUC (BbsI)
    #   and their reverse complements: Type IIS restriction sites near-universally avoided inside
    #   cloned DNA parts in modern Golden Gate/MoClo synthetic-biology workflows.
    # NOT included by default (context-dependent, add them yourself if relevant to your workflow):
    #   classic 6-cutter sites (EcoRI/BamHI/HindIII etc - less relevant with Golden-Gate-style
    #   domestication), and Shine-Dalgarno-like (RBS) motifs, which only matter if this sequence
    #   ends up placed near a start codon in a larger transcript.
    forbidden_motifs: Tuple[str, ...] = (
        "AAAAA", "UUUUU", "GGGGG", "CCCCC",
        "GGUCUC", "GAGACC",     # BsaI + reverse complement
        "CGUCUC", "GAGACG",     # BsmBI/Esp3I + reverse complement
        "GCUCUUC", "GAAGAGC",   # SapI + reverse complement
        "GAAGAC", "GUCUUC",     # BbsI + reverse complement
    )

    # --- fixed reference scales used to convert physical quantities (criteria 4/5/6/7) into an
    # absolute 0-1 quality, independent of whatever else is in the candidate pool this run.
    # These are reasonable engineering heuristics, NOT values pinned down by a specific paper -
    # tune them if you have a better estimate for your system.
    dg_per_bp_ref: float = -1.5             # kcal/mol per base pair; rough typical RNA:RNA Watson-Crick
                                              # duplex stability at 37C (SantaLucia-style nearest-neighbor
                                              # params give roughly -1 to -3 kcal/mol/bp depending on
                                              # sequence). A trigger-vs-switch duplex energy this strong
                                              # *per possible base pair* counts as "perfect" (quality=1.0).
    unintended_match_nt_ref: float = 8.0     # nt; an unintended identical/complementary stretch this long
                                              # between the two triggers counts as "fully disqualifying"
                                              # (quality=0.0) for criteria 5/6. Chosen to line up with the
                                              # upper end of this script's own multi-window similarity scan
                                              # (sim_report_windows goes up to 12nt).

    # --- backend ---
    backend: str = "auto"   # 'auto' | 'nupack' | 'vienna' | 'fallback'


# ======================================================================
# 2) Sequence utilities (backend-agnostic)
# ======================================================================
_COMP_RNA = str.maketrans("ACGU", "UGCA")


def to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def reverse_complement(seq: str) -> str:
    # 'N' (unknown) stays 'N' after translation
    return to_rna(seq).translate(_COMP_RNA)[::-1]


def gc_content(seq: str) -> float:
    return sum(1 for b in seq.upper() if b in "GC") / len(seq) if seq else 0.0


def has_forbidden(seq: str, motifs: Tuple[str, ...]) -> bool:
    s = to_rna(seq)
    return any(m in s for m in motifs)


def _longest_common_substring(a: str, b: str) -> int:
    """Returns the length of the longest shared contiguous substring between a and b (DP). Ignores 'N'."""
    a, b = to_rna(a), to_rna(b)
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        ai = a[i - 1]
        for j in range(1, len(b) + 1):
            if ai != "N" and ai == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best = cur[j]
        prev = cur
    return best


def max_identity_match(a: str, b: str) -> int:
    """Finds the longest identical stretch between a and b => the 'substitution' score."""
    return _longest_common_substring(a, b)


def max_revcomp_match(a: str, b: str) -> int:
    """Finds the longest stretch where a's reverse complement (antiparallel) matches b => the 'sticking' score."""
    return _longest_common_substring(a, reverse_complement(b))


def window_hits(a: str, b: str, window: int, mode: str) -> bool:
    """As in MODesign: checks if there's any shared `window`-length substring between a and b, either
    identical (mode='id') or complementary (mode='rc'). Used only for reporting."""
    target = b if mode == "id" else reverse_complement(b)
    bwins = {target[i:i + window] for i in range(len(target) - window + 1)}
    for i in range(len(a) - window + 1):
        w = a[i:i + window]
        if "N" not in w and w in bwins:
            return True
    return False


def _region_openness(structure: str, start: int, length: int) -> float:
    """% of unpaired bases ('.') in the region [start:start+length] of a dot-bracket structure
    (similar to openess_mfe_structure in MODesign, but for a sub-region)."""
    region = structure[start:start + length]
    return 100.0 * region.count(".") / len(region) if region else 0.0


def _mask(seq: str, start: int, length: int) -> str:
    """Replaces a region with 'N' (used to hide the intended '*'/K2 connector so it's ignored in similarity checks)."""
    return seq[:start] + "N" * length + seq[start + length:]


def build_switch_target(r1: str, r2: str, a: str, k1: str) -> str:
    """The switch's target site = reverse complement of the mRNA's target region
    (R1+R2+A+K1). ('*'/K2 is the connector between the two triggers, it is not part of the switch)."""
    return reverse_complement(to_rna(r1) + to_rna(r2) + a + k1)


# ======================================================================
# 3) BACKEND layer - wrappers around NUPACK / ViennaRNA
# ======================================================================
# Common interface:
#   fold(seq)            -> (mfe_structure, mfe_energy, ensemble_energy)
#   unpaired_probs(seq)  -> [p_unpaired for each base]
#   duplex_energy(a, b)  -> energy
# plus helpers: mfe(), structure_probability().

class Backend:
    name = "base"

    def fold(self, seq: str) -> Tuple[str, float, float]:
        raise NotImplementedError

    def unpaired_probs(self, seq: str) -> List[float]:
        raise NotImplementedError

    def duplex_energy(self, a: str, b: str) -> float:
        raise NotImplementedError

    def mfe(self, seq: str) -> Tuple[str, float]:
        s, e, _ = self.fold(seq)
        return s, e

    def structure_probability(self, seq: str, temperature_c: float) -> float:
        """Boltzmann probability that a sequence folds into its own MFE structure (as in MODesign):
        p = exp(-(E_mfe - E_ensemble)/kT).  1.0 => structure is very stable; close to 0 => not very likely."""
        _, e_mfe, e_ens = self.fold(seq)
        kT = 0.001987 * (temperature_c + 273.15)
        return math.exp(-(e_mfe - e_ens) / kT) if kT else 0.0


# ---------------------------------------------------------------- ViennaRNA
class ViennaBackend(Backend):
    """ViennaRNA (import RNA).  Functions used:
        RNA.cvar.temperature, RNA.fold_compound, fc.mfe(), fc.pf() (for ensemble energy),
        fc.bpp() (base-pair probabilities => used for accessibility), RNA.duplexfold() (binding between two sequences)."""
    name = "vienna"

    def __init__(self, temperature_c: float = 37.0):
        import RNA  # type: ignore
        self.RNA = RNA
        RNA.cvar.temperature = temperature_c

    def fold(self, seq: str) -> Tuple[str, float, float]:
        seq = to_rna(seq)
        fc = self.RNA.fold_compound(seq)
        structure, mfe_e = fc.mfe()
        _, ens_e = fc.pf()   # for ensemble energy (partition function)
        return structure, float(mfe_e), float(ens_e)

    def unpaired_probs(self, seq: str) -> List[float]:
        seq = to_rna(seq)
        n = len(seq)
        fc = self.RNA.fold_compound(seq)
        fc.pf()
        bpp = fc.bpp()  # 1-based
        u = [1.0] * n
        for i in range(1, n + 1):
            paired = sum(bpp[i][j] for j in range(i + 1, n + 1)) + \
                     sum(bpp[j][i] for j in range(1, i))
            u[i - 1] = max(0.0, 1.0 - paired)
        return u

    def duplex_energy(self, a: str, b: str) -> float:
        return float(self.RNA.duplexfold(to_rna(a), to_rna(b)).energy)


# ------------------------------------------------------------------- NUPACK
class NupackBackend(Backend):
    """NUPACK 4 (import nupack) - matches the functions MODesign uses:
        nupack.Model, nupack.mfe, nupack.pfunc, nupack.energy,
        complex_analysis(compute=['pairs']) for pairing."""
    name = "nupack"

    def __init__(self, temperature_c: float = 37.0, material: str = "rna"):
        import nupack  # type: ignore
        self.nupack = nupack
        self.material = material
        self.model = nupack.Model(material=material, celsius=temperature_c)

    def fold(self, seq: str) -> Tuple[str, float, float]:
        seq = to_rna(seq)
        res = self.nupack.mfe(strands=[seq], model=self.model)
        struct, mfe_e = str(res[0].structure), float(res[0].energy)
        pf = self.nupack.pfunc(strands=[seq], model=self.model)
        ens_e = float(pf[1])   # for ensemble energy
        return struct, mfe_e, ens_e

    def unpaired_probs(self, seq: str) -> List[float]:
        seq = to_rna(seq)
        n = len(seq)
        try:
            strand = self.nupack.Strand(seq, name="s")
            cx = self.nupack.Complex([strand])
            result = self.nupack.complex_analysis([cx], model=self.model, compute=["pairs"])
            pairs = result[cx].pairs.to_array()  # (n+1)x(n+1); index n = "unpaired"
            return [float(pairs[i][n]) for i in range(n)]
        except Exception:
            return [0.5] * n

    def duplex_energy(self, a: str, b: str) -> float:
        res = self.nupack.mfe(strands=[to_rna(a), to_rna(b)], model=self.model)
        return float(res[0].energy)


# ----------------------------------------------------------------- FALLBACK
class FallbackBackend(Backend):
    """*** does not fold RNA properly! ***  a simple approximation (Nussinov) for when
    ViennaRNA/NUPACK are not installed. Use only for a quick check, not for real conclusions."""
    name = "fallback"
    _PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}

    def __init__(self, temperature_c: float = 37.0):
        self.temperature_c = temperature_c

    def _nussinov(self, seq: str, min_loop: int = 3) -> Tuple[List[Tuple[int, int]], int]:
        n = len(seq)
        dp = [[0] * n for _ in range(n)]
        for span in range(min_loop + 1, n):
            for i in range(0, n - span):
                j = i + span
                best = max(dp[i + 1][j], dp[i][j - 1])
                if (seq[i], seq[j]) in self._PAIRS and j - i > min_loop:
                    best = max(best, dp[i + 1][j - 1] + 1)
                for k in range(i + 1, j):
                    best = max(best, dp[i][k] + dp[k + 1][j])
                dp[i][j] = best
        pairs: List[Tuple[int, int]] = []
        stack = [(0, n - 1)]
        while stack:
            i, j = stack.pop()
            if i >= j:
                continue
            if dp[i][j] == dp[i + 1][j]:
                stack.append((i + 1, j))
            elif dp[i][j] == dp[i][j - 1]:
                stack.append((i, j - 1))
            elif (seq[i], seq[j]) in self._PAIRS and dp[i][j] == dp[i + 1][j - 1] + 1:
                pairs.append((i, j))
                stack.append((i + 1, j - 1))
            else:
                for k in range(i + 1, j):
                    if dp[i][j] == dp[i][k] + dp[k + 1][j]:
                        stack.append((i, k)); stack.append((k + 1, j)); break
        return pairs, dp[0][n - 1]

    def fold(self, seq: str) -> Tuple[str, float, float]:
        seq = to_rna(seq)
        pairs, npair = self._nussinov(seq)
        dot = ["."] * len(seq)
        for i, j in pairs:
            dot[i], dot[j] = "(", ")"
        mfe_e = -2.0 * npair
        return "".join(dot), mfe_e, mfe_e - 1.0   # rough ensemble energy estimate (approximate)

    def unpaired_probs(self, seq: str) -> List[float]:
        seq = to_rna(seq)
        pairs, _ = self._nussinov(seq)
        paired = {i for p in pairs for i in p}
        return [0.15 if i in paired else 0.85 for i in range(len(seq))]

    def duplex_energy(self, a: str, b: str) -> float:
        a, b = to_rna(a), to_rna(b)
        rc = reverse_complement(b)
        best = 0
        for off in range(-len(rc) + 1, len(a)):
            run = sum(1 for k in range(len(rc))
                      if 0 <= off + k < len(a) and a[off + k] == rc[k])
            best = max(best, run)
        return -2.0 * best


def select_backend(p: Params) -> Backend:
    order = [p.backend] if p.backend != "auto" else ["nupack", "vienna", "fallback"]
    errors = {}
    for name in order:
        try:
            if name == "nupack":
                return NupackBackend(p.temperature_c, p.material)
            if name == "vienna":
                return ViennaBackend(p.temperature_c)
            if name == "fallback":
                return FallbackBackend(p.temperature_c)
        except Exception as exc:
            errors[name] = str(exc)
            if p.backend != "auto":
                raise RuntimeError(f"Backend '{p.backend}' failed: {exc}") from exc
    raise RuntimeError(f"No backend available. Errors: {errors}")


# ======================================================================
# 4) Scanning genes for endogenous Trigger-1 / Trigger-2 windows
# ======================================================================
@dataclass
class Trigger1Window:
    """One contiguous candidate window for Trigger 1, sliced straight out of a gene: R1 + '*' + A + K1."""
    gene: str
    start: int
    seq: str             # full window = r1 + star + a + k1
    r1: str
    star: str
    a: str
    k1: str
    mean_unpaired: float  # accessibility of the R1 sub-region only


@dataclass
class Trigger2Window:
    """One contiguous candidate window for Trigger 2, sliced straight out of a gene: R2 + K2."""
    gene: str
    start: int
    seq: str             # full window = r2 + k2
    r2: str
    k2: str
    mean_unpaired: float  # accessibility of the R2 sub-region only


def scan_trigger1(genes: List[Tuple[str, str]], backend: Backend, p: Params) -> List[Trigger1Window]:
    """Slides a window of length R1_len+star_len+A_len+K1_len across every gene in `genes`
    (a list of (name, sequence) records) and returns every window found, cut into R1/'*'/A/K1."""
    win_len = p.R1_len + p.star_len + p.A_len + p.K1_len
    windows: List[Trigger1Window] = []
    for name, gene_seq in genes:
        gene = to_rna(gene_seq)
        if len(gene) < win_len:
            continue
        unpaired = backend.unpaired_probs(gene)
        for s in range(0, len(gene) - win_len + 1, p.scan_step):
            w = gene[s:s + win_len]
            if has_forbidden(w, p.forbidden_motifs):
                continue
            r1 = w[:p.R1_len]
            star = w[p.R1_len:p.R1_len + p.star_len]
            a = w[p.R1_len + p.star_len:p.R1_len + p.star_len + p.A_len]
            k1 = w[p.R1_len + p.star_len + p.A_len:]
            if p.max_gc_content > 0 and gc_content(r1) > p.max_gc_content:
                continue
            mu = sum(unpaired[s:s + p.R1_len]) / p.R1_len
            windows.append(Trigger1Window(name, s, w, r1, star, a, k1, mu))
    return windows


def scan_trigger2(genes: List[Tuple[str, str]], backend: Backend, p: Params) -> List[Trigger2Window]:
    """Slides a window of length R2_len+star_len across every gene in `genes`
    (a list of (name, sequence) records) and returns every window found, cut into R2/K2."""
    win_len = p.R2_len + p.star_len
    windows: List[Trigger2Window] = []
    for name, gene_seq in genes:
        gene = to_rna(gene_seq)
        if len(gene) < win_len:
            continue
        unpaired = backend.unpaired_probs(gene)
        for s in range(0, len(gene) - win_len + 1, p.scan_step):
            w = gene[s:s + win_len]
            if has_forbidden(w, p.forbidden_motifs):
                continue
            r2 = w[:p.R2_len]
            k2 = w[p.R2_len:]
            if p.max_gc_content > 0 and gc_content(r2) > p.max_gc_content:
                continue
            mu = sum(unpaired[s:s + p.R2_len]) / p.R2_len
            windows.append(Trigger2Window(name, s, w, r2, k2, mu))
    return windows


def find_matching_pairs(win1: List[Trigger1Window],
                         win2: List[Trigger2Window]) -> List[Tuple[Trigger1Window, Trigger2Window]]:
    """Pairs up Trigger-1 and Trigger-2 windows where '*' is the exact reverse-complement of K2 -
    the only way two independently-found endogenous half-triggers could plausibly bind each other.
    Uses a hash map on '*' so this is O(len(win1) + len(win2)) instead of comparing every pair.
    Both windows are scanned from the same pooled gene set, so a pair is only kept if its two
    windows come from two DIFFERENT gene records - a real AND gate needs two distinct inputs,
    not the same gene sensing itself twice."""
    by_star = defaultdict(list)
    for w1 in win1:
        by_star[w1.star].append(w1)
    pairs: List[Tuple[Trigger1Window, Trigger2Window]] = []
    for w2 in win2:
        needed_star = reverse_complement(w2.k2)
        for w1 in by_star.get(needed_star, []):
            if w1.gene == w2.gene:
                continue
            pairs.append((w1, w2))
    return pairs


# ======================================================================
# 5) Scoring a matched (Trigger 1, Trigger 2) pair
# ======================================================================
@dataclass
class Candidate:
    w1: Trigger1Window
    w2: Trigger2Window
    trigger1: str   # Trigger 1 = R1 + '*' + A + K1
    trigger2: str   # Trigger 2 = R2 + K2
    switch: str
    # scoring metrics
    access: float            # (1) target accessibility (R1+R2)/2                    [higher=better]
    open_trigger1: float     # (2) % openness of Trigger 1's functional region       [higher=better]
    open_trigger2: float     # (2) % openness of Trigger 2's functional region       [higher=better]
    prob_trigger1: float     # (3) structure stability probability, Trigger 1        [higher=better]
    prob_trigger2: float     # (3) structure stability probability, Trigger 2        [higher=better]
    bind_trigger1: float     # (4) binding energy Trigger 1<->switch                 [lower=better]
    bind_trigger2: float     # (4) binding energy Trigger 2<->switch                 [lower=better]
    stick_nt: int            # (5) longest unintended reverse-complement match        [lower=better]
    subst_nt: int            # (6) longest unintended identical match                 [lower=better]
    bind_full: float = 0.0   # (7) binding energy of the *fully combined* (R1+R2+A+K1) complex to the
                             #     switch - the ideal/reference strength when both triggers are present
    leak_margin: float = 0.0  # (7) bind_full minus the strongest single-trigger-alone binding
                              #     (bind_trigger1/bind_trigger2) - large & positive = low leak risk  [higher=better]
    score: float = 0.0
    # per-criterion quality, 0-1, against a FIXED reference scale (see module docstring) -
    # NOT relative to other candidates in this run, so these are comparable across runs.
    # Filled in by evaluate_pair(); used by report() to show how each criterion
    # contributed to the final 0-100 score.
    norm_access: float = 0.0
    norm_openness: float = 0.0
    norm_robustness: float = 0.0
    norm_binding: float = 0.0
    norm_stick: float = 0.0
    norm_substitute: float = 0.0
    norm_specificity: float = 0.0


def _binding_quality(energy: float, len_a: int, len_b: int, p: Params) -> float:
    """Converts a duplex free energy into an absolute 0-1 quality: energy per possible base pair
    (using the shorter of the two strands as the max possible bp count), relative to a fixed
    reference kcal/mol/bp (Params.dg_per_bp_ref). See that field's comment for what the reference means."""
    max_bp = max(1, min(len_a, len_b))
    q = energy / (p.dg_per_bp_ref * max_bp)
    return max(0.0, min(1.0, q))


def evaluate_pair(w1: Trigger1Window, w2: Trigger2Window, backend: Backend, p: Params) -> Candidate:
    trigger1 = w1.seq   # R1 + '*' + A + K1, exactly as found in the gene
    trigger2 = w2.seq   # R2 + K2, exactly as found in the gene
    switch = build_switch_target(w1.r1, w2.r2, w1.a, w1.k1)

    # (2) trigger self-folding: fold *both* triggers and check how open their functional regions are.
    trigger1_struct, _, _ = backend.fold(trigger1)
    trigger2_struct, _, _ = backend.fold(trigger2)
    # Trigger 1: functional region = R1 + toehold A.  Trigger 2: functional region = R2 (target region).
    open_trigger1_R1 = _region_openness(trigger1_struct, 0, len(w1.r1))
    a_start = len(w1.r1) + len(w1.star)
    open_trigger1_A = _region_openness(trigger1_struct, a_start, len(w1.a))
    open_trigger1 = 0.5 * (open_trigger1_R1 + open_trigger1_A)
    open_trigger2 = _region_openness(trigger2_struct, 0, len(w2.r2))

    # (3) structure stability
    prob_trigger1 = backend.structure_probability(trigger1, p.temperature_c)
    prob_trigger2 = backend.structure_probability(trigger2, p.temperature_c)

    # (4) switch binding strength - each trigger *alone* against the switch. Note this is a
    # single-strand binding energy, i.e. exactly the "one trigger alone" case that criterion (7)
    # below checks is NOT too strong (a strong single-trigger binding here means AND-gate leak).
    bind_trigger1 = backend.duplex_energy(trigger1, switch)
    bind_trigger2 = backend.duplex_energy(trigger2, switch)

    # (5)+(6) compare the two triggers to each other, masking the intended '*'/K2 connector
    trigger1_masked = _mask(trigger1, len(w1.r1), len(w1.star))              # mask '*'
    trigger2_masked = _mask(trigger2, len(w2.r2), len(w2.k2))                # mask K2
    stick_nt = max_revcomp_match(trigger1_masked, trigger2_masked)           # sticking (complementary)
    subst_nt = max_identity_match(trigger1_masked, trigger2_masked)          # substitution (identical)

    # (7) AND-gate specificity / anti-leak (Kim & Moon, ACS Synth Biol 2020): a single trigger alone
    # should bind the switch much more weakly than the fully combined R1+R2+A+K1 complex - otherwise
    # the "AND gate" leaks like a one-input switch (see their Figure 1C/1D).
    bind_full = backend.duplex_energy(w1.r1 + w2.r2 + w1.a + w1.k1, switch)
    leak_margin = min(bind_trigger1, bind_trigger2) - bind_full

    access = 0.5 * (w1.mean_unpaired + w2.mean_unpaired)

    # Convert every criterion to an absolute 0-1 quality against a FIXED reference scale (not
    # relative to any other candidate) - see the module docstring for why this matters.
    norm_access = access                                                     # already an absolute probability
    norm_openness = 0.01 * 0.5 * (open_trigger1 + open_trigger2)             # already an absolute %
    norm_robustness = 0.5 * (prob_trigger1 + prob_trigger2)                  # already an absolute probability
    norm_binding = 0.5 * (_binding_quality(bind_trigger1, len(trigger1), len(switch), p) +
                          _binding_quality(bind_trigger2, len(trigger2), len(switch), p))
    norm_stick = max(0.0, 1.0 - stick_nt / p.unintended_match_nt_ref)
    norm_substitute = max(0.0, 1.0 - subst_nt / p.unintended_match_nt_ref)
    norm_specificity = 0.0 if bind_full == 0 else max(0.0, min(1.0, 1.0 - min(bind_trigger1, bind_trigger2) / bind_full))

    max_possible = (p.w_access + p.w_openness + p.w_robustness +
                    p.w_binding + p.w_stick + p.w_substitute + p.w_specificity)
    raw = (p.w_access * norm_access + p.w_openness * norm_openness + p.w_robustness * norm_robustness +
           p.w_binding * norm_binding + p.w_stick * norm_stick + p.w_substitute * norm_substitute +
           p.w_specificity * norm_specificity)
    score = 100.0 * raw / max_possible if max_possible > 0 else 0.0

    return Candidate(w1, w2, trigger1, trigger2, switch,
                     access=access, open_trigger1=open_trigger1, open_trigger2=open_trigger2,
                     prob_trigger1=prob_trigger1, prob_trigger2=prob_trigger2,
                     bind_trigger1=bind_trigger1, bind_trigger2=bind_trigger2,
                     stick_nt=stick_nt, subst_nt=subst_nt,
                     bind_full=bind_full, leak_margin=leak_margin,
                     score=score,
                     norm_access=norm_access, norm_openness=norm_openness, norm_robustness=norm_robustness,
                     norm_binding=norm_binding, norm_stick=norm_stick, norm_substitute=norm_substitute,
                     norm_specificity=norm_specificity)


def _passes_hard_filters(c: Candidate, p: Params) -> bool:
    if p.min_functional_openness > 0:
        if c.open_trigger1 < p.min_functional_openness or c.open_trigger2 < p.min_functional_openness:
            return False
    if p.min_structure_probability > 0:
        if c.prob_trigger1 < p.min_structure_probability or c.prob_trigger2 < p.min_structure_probability:
            return False
    if p.max_unintended_match > 0:
        if max(c.stick_nt, c.subst_nt) >= p.max_unintended_match:
            return False
    if p.max_leak_fraction > 0 and c.bind_full != 0:
        leak_fraction = min(c.bind_trigger1, c.bind_trigger2) / c.bind_full
        if leak_fraction > p.max_leak_fraction:
            return False
    return True


def rank_candidates(cands: List[Candidate], p: Params) -> Tuple[List[Candidate], int]:
    """Filters out candidates that fail the hard filters, then sorts by score. The score itself
    is already computed per-candidate in evaluate_pair() against a fixed reference scale, so no
    normalization against the rest of the pool happens here (see module docstring)."""
    kept = [c for c in cands if _passes_hard_filters(c, p)]
    n_filtered = len(cands) - len(kept)
    kept.sort(key=lambda c: c.score, reverse=True)
    return kept, n_filtered


# ======================================================================
# 6) Main pipeline
# ======================================================================
def design(genes: List[Tuple[str, str]], p: Params):
    """`genes` is a pool of (name, sequence) records - it can come from one FASTA file with many
    records, several files pooled together, or a mix. Trigger 1 and Trigger 2 candidate windows are
    both scanned across the *entire* pool; find_matching_pairs() then picks whichever two distinct
    gene records give the best-matching pair."""
    backend = select_backend(p)
    win1 = scan_trigger1(genes, backend, p)
    win2 = scan_trigger2(genes, backend, p)
    pairs = find_matching_pairs(win1, win2)
    n_matches = len(pairs)
    # fold/score is the expensive step, so only do it for the most-accessible matches
    pairs.sort(key=lambda wp: 0.5 * (wp[0].mean_unpaired + wp[1].mean_unpaired), reverse=True)
    pairs = pairs[: p.top_k]
    cands = [evaluate_pair(w1, w2, backend, p) for w1, w2 in pairs]
    ranked, n_filtered = rank_candidates(cands, p)
    return backend, ranked, n_filtered, len(win1), len(win2), n_matches


# ======================================================================
# 7) Input/output
# ======================================================================
def read_fasta_records(path: str) -> List[Tuple[str, str]]:
    """Parses a FASTA file that may contain one or more gene records.
    Returns a list of (name, sequence) tuples, in file order."""
    records: List[Tuple[str, str]] = []
    name = None
    seq_chunks: List[str] = []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq_chunks)))
            name = line[1:].strip() or f"record{len(records) + 1}"
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if name is not None:
        records.append((name, "".join(seq_chunks)))
    return records


def report(backend, cands, n_filtered, n_win1, n_win2, n_matches, p: Params, top: int = 3) -> None:
    print("=" * 74)
    print(f"Backend: {backend.name}"
          + ("   (*** FALLBACK - not very accurate! install ViennaRNA/NUPACK ***)"
             if backend.name == "fallback" else ""))
    print(f"T={p.temperature_c}C | star_len={p.star_len} R1={p.R1_len} R2={p.R2_len} "
          f"A={p.A_len} K1={p.K1_len}")
    print(f"Scanned {n_win1} Trigger-1 candidate windows and "
          f"{n_win2} Trigger-2 candidate windows across the pooled gene set")
    print(f"Found {n_matches} exact '*' == reverse_complement(K2) matches"
          + (f"; folded/scored the top {len(cands) + n_filtered}" if n_matches > 0 else ""))
    if n_filtered:
        print(f"Filtered out by hard filters: {n_filtered} candidates")
    print("=" * 74)
    if not cands:
        print("No matching endogenous trigger pairs found.")
        print("Try: a smaller star_len (fewer nt need to match by chance), different genes, "
              "or relaxing the hard filters (min_functional_openness / max_unintended_match / max_leak_fraction).")
        return
    for rank, c in enumerate(cands[:top], 1):
        print(f"\n#{rank}  score={c.score:.1f}/100")
        print(f"  Trigger 1 found in '{c.w1.gene}' @{c.w1.start} (R1 accessibility {c.w1.mean_unpaired:.2f})")
        print(f"     R1({p.R1_len}nt) = {c.w1.r1}")
        print(f"     *({p.star_len}nt)  = {c.w1.star}")
        print(f"     A({p.A_len}nt)  = {c.w1.a}")
        print(f"     K1({p.K1_len}nt) = {c.w1.k1}")
        print(f"     full sequence   = {c.trigger1}")
        print(f"  Trigger 2 found in '{c.w2.gene}' @{c.w2.start} (R2 accessibility {c.w2.mean_unpaired:.2f})")
        print(f"     R2({p.R2_len}nt) = {c.w2.r2}")
        print(f"     K2({p.star_len}nt) = {c.w2.k2}")
        print(f"     full sequence   = {c.trigger2}")
        max_possible = (p.w_access + p.w_openness + p.w_robustness +
                        p.w_binding + p.w_stick + p.w_substitute + p.w_specificity)

        def _pts(weight: float, norm: float) -> float:
            return 100.0 * weight * norm / max_possible if max_possible > 0 else 0.0

        def _max_pts(weight: float) -> float:
            return 100.0 * weight / max_possible if max_possible > 0 else 0.0

        def _row(label: str, raw: str, norm: float, weight: float) -> str:
            return (f"       {label:<28}{raw:<48}"
                    f"quality={norm * 100:3.0f}%  ->  {_pts(weight, norm):5.1f} / {_max_pts(weight):4.1f} pts  (weight={weight})")

        print("     score breakdown  (quality% is absolute, against a fixed reference scale - "
              "comparable across runs, NOT just this run's candidates):")
        print(_row("(1) accessibility", f"raw={c.access:.2f}", c.norm_access, p.w_access))
        print(_row("(2) self-fold openness", f"raw: trigger1={c.open_trigger1:.0f}% trigger2={c.open_trigger2:.0f}%",
                    c.norm_openness, p.w_openness))
        print(_row("(3) structure stability (p)", f"raw: trigger1={c.prob_trigger1:.2f} trigger2={c.prob_trigger2:.2f}",
                    c.norm_robustness, p.w_robustness))
        print(_row("(4) switch binding dG", f"raw: trigger1={c.bind_trigger1:.1f} trigger2={c.bind_trigger2:.1f}",
                    c.norm_binding, p.w_binding))
        print(_row("(5) sticking (RC nt)", f"raw={c.stick_nt}nt", c.norm_stick, p.w_stick))
        print(_row("(6) substitution (nt)", f"raw={c.subst_nt}nt", c.norm_substitute, p.w_substitute))
        print(_row("(7) AND-gate specificity", f"raw: alone1={c.bind_trigger1:.1f} alone2={c.bind_trigger2:.1f} full={c.bind_full:.1f}",
                    c.norm_specificity, p.w_specificity))
        print(f"       {'-' * 96}")
        print(f"       {'TOTAL':<28}{'':<48}{'':16}->  {c.score:5.1f} / 100.0 pts")
        # multi-window similarity report between triggers (masking the intended '*'/K2 connector)
        trigger1_m = _mask(c.trigger1, len(c.w1.r1), len(c.w1.star))
        trigger2_m = _mask(c.trigger2, len(c.w2.r2), len(c.w2.k2))
        hits = [f"{w}:{('RC' if window_hits(trigger1_m, trigger2_m, w, 'rc') else '')}"
                f"{('ID' if window_hits(trigger1_m, trigger2_m, w, 'id') else '')}" or f"{w}:-"
                for w in p.sim_report_windows]
        print(f"     multi-window similarity between triggers (window:sticking/substitution): {' '.join(hits)}")


_DEMO_GENE1 = ("AUGAGCAAAGGUGAAGAACUGUUUACCGGCGUGGUGCCGAUUCUGGUGGAACUGGAUGGCGAU"
               "GUGAACGGCCAUAAAUUUAGCGUGAGCGGCGAAGGCGAAGGCGAUGCGACCUAUGGCAAACUG")
_DEMO_GENE2 = ("AUGGUGAGCAAGGGCGAGGAGCUGUUCACCGGGGUGGUGCCCAUCCUGGUCGAGCUGGACGGC"
               "GACGUAAACGGCCACAAGUUCAGCGUGUCCGGCGAGGGCGAGGGCGAUGCCACCUACGGCAAG")


def _genes_folder() -> str:
    """A folder named 'genes' next to this script file. Drop any FASTA files in there and they're
    picked up automatically, with no need to pass file paths on the command line.
    Falls back to the current working directory when run from a context with no `__file__`
    (e.g. a Jupyter/VS Code interactive cell instead of `python script.py`)."""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        base_dir = os.getcwd()
    return os.path.join(base_dir, "genes")


def _load_genes_folder(folder: str) -> List[Tuple[str, str]]:
    genes: List[Tuple[str, str]] = []
    for fname in sorted(os.listdir(folder)):
        if fname.lower().endswith((".fa", ".fasta")):
            genes.extend(read_fasta_records(os.path.join(folder, fname)))
    return genes


def main(argv: List[str]) -> int:
    p = Params()
    genes_folder = _genes_folder()
    folder_genes = _load_genes_folder(genes_folder) if os.path.isdir(genes_folder) else []
    # Only treat argv entries as gene files if they actually look like FASTA files on disk -
    # this also makes `main` safe to run from a Jupyter/VS Code interactive cell, where sys.argv
    # holds kernel launch flags (e.g. "-f", a connection.json path) rather than real arguments.
    cli_paths = [a for a in argv[1:] if a.lower().endswith((".fa", ".fasta", ".fna")) and os.path.isfile(a)]
    if cli_paths:
        genes: List[Tuple[str, str]] = []
        for path in cli_paths:
            genes.extend(read_fasta_records(path))
        print(f"(loaded {len(genes)} gene record(s) from {len(cli_paths)} file(s))\n")
    elif folder_genes:
        genes = folder_genes
        print(f"(loaded {len(genes)} gene record(s) from '{genes_folder}')\n")
    elif os.path.exists("gene1.fa") and os.path.exists("gene2.fa"):
        genes = read_fasta_records("gene1.fa") + read_fasta_records("gene2.fa")
        print(f"(loaded gene1.fa + gene2.fa: {len(genes)} gene record(s) total)\n")
    else:
        genes = [("demo_gene1", _DEMO_GENE1), ("demo_gene2", _DEMO_GENE2)]
        print("(no FASTA file provided - running with a demo example)\n")
    backend, cands, n_filtered, n_win1, n_win2, n_matches = design(genes, p)
    report(backend, cands, n_filtered, n_win1, n_win2, n_matches, p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
