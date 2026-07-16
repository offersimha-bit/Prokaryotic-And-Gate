"""Pure-sequence helpers: alphabet conversion, reverse complement, translation,
Hamming distance, and the restricted-sequence checks from Section 6.

Everything works in either DNA (T) or RNA (U) alphabets; internally the thermo
layer uses RNA (U), matching the Toehold-VISTA convention.
"""

from __future__ import annotations

from typing import Iterable

_COMPLEMENT = {
    "A": "T", "T": "A", "U": "A", "C": "G", "G": "C", "N": "N",
}

_CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 'S',
    'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L',
    'CUA': 'L', 'CUG': 'L', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAU': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {"UAA", "UAG", "UGA"}


def clean(seq: str) -> str:
    """Uppercase and strip whitespace; leaves alphabet (T vs U) untouched."""
    return "".join(seq.split()).upper()


def to_rna(seq: str) -> str:
    return clean(seq).replace("T", "U")


def to_dna(seq: str) -> str:
    return clean(seq).replace("U", "T")


def _is_rna(s: str) -> bool:
    """RNA unless the sequence actually contains T (DNA).  A fragment with
    neither T nor U (e.g. a U-less k-mer) defaults to RNA, matching this
    pipeline's RNA-internal convention -- otherwise its complement would leak a
    spurious T."""
    return "T" not in s


def reverse_complement(seq: str) -> str:
    """Reverse complement.  Output alphabet is RNA unless the input contains T.

    A mixed/unknown base raises ``KeyError`` so callers fail loudly rather than
    silently designing against a corrupt target.
    """
    s = clean(seq)
    comp = "".join(_COMPLEMENT[b] for b in reversed(s))
    return comp.replace("T", "U") if _is_rna(s) else comp


def complement(seq: str) -> str:
    """Element-wise (non-reversed) complement; RNA output unless input has T."""
    s = clean(seq)
    comp = "".join(_COMPLEMENT[b] for b in s)
    return comp.replace("T", "U") if _is_rna(s) else comp


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError(f"Hamming distance needs equal lengths ({len(a)} vs {len(b)})")
    return sum(1 for x, y in zip(clean(a), clean(b)) if x != y)


def translate(seq: str) -> str:
    """Translate an RNA (or DNA) sequence in frame 0 to single-letter aa,
    'X' for incomplete/unknown codons and '*' for stops."""
    rna = to_rna(seq)
    out = []
    for i in range(0, len(rna) - 2, 3):
        out.append(_CODON_TABLE.get(rna[i:i + 3], "X"))
    return "".join(out)


def has_inframe_stop(seq: str, start: int = 0) -> bool:
    """True if there is an in-frame stop codon at or after ``start`` (frame set
    by ``start``), excluding a stop that sits exactly at the final codon."""
    rna = to_rna(seq)
    codons = [rna[i:i + 3] for i in range(start, len(rna) - 2, 3)]
    # ignore a legitimate terminal stop codon
    body = codons[:-1] if codons and codons[-1] in STOP_CODONS else codons
    return any(c in STOP_CODONS for c in body)


def find_all(seq: str, motif: str) -> list[int]:
    """All 0-based start indices where ``motif`` occurs (overlapping)."""
    s, m = clean(seq), clean(motif)
    idx, out = s.find(m), []
    while idx != -1:
        out.append(idx)
        idx = s.find(m, idx + 1)
    return out


def count_aug_after_rbs(switch_rna: str, rbs_seq: str) -> int:
    """Number of AUG codons located downstream of the RBS.  The design should
    present exactly one (the intended start codon)."""
    s = to_rna(switch_rna)
    rbs = to_rna(rbs_seq)
    rbs_idx = s.find(rbs)
    if rbs_idx == -1:
        return len(find_all(s, "AUG"))
    return len([i for i in find_all(s, "AUG") if i > rbs_idx])


def has_forbidden_run(seq: str, runs: Iterable[str]) -> list[str]:
    """Return the list of forbidden homopolymer runs present in ``seq``."""
    s = to_rna(seq)
    return [r for r in runs if to_rna(r) in s]


def gc_fraction(seq: str) -> float:
    s = clean(seq)
    if not s:
        return 0.0
    return (s.count("G") + s.count("C")) / len(s)
