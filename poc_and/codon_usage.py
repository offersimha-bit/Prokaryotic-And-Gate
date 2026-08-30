"""
E. coli codon usage, and the codon metrics Green 2026 VISTA reports.

The table is the one VISTA itself uses, already in this repo at
``and_gate_pipeline/data/ecoli_codon_usage_table.csv`` (64 rows: Codon, Amino
acid, Fraction, Frequency/Thousand, Number). "Fraction" is the share of that
amino acid's codons which are this one, in E. coli -- so 1.0 would mean the
only choice and 0.05 means the organism almost never uses it.

Two uses here:

  * a TIE-BREAK when several synonymous codons are equally good on the
    objective we actually care about (fewest edits, or most divergence). This
    is free -- it never sacrifices the primary objective, it just stops the
    choice being made by list order, which is what it was before.

  * a REPORTED metric, so a recoded window that quietly filled up with codons
    E. coli dislikes is visible rather than silent.

Worth recording, because it changes a working assumption: the codon-max mCherry
that failed at the bench has a mean usage fraction of 0.471, BETTER than the
original's 0.384. So poor codon usage is unlikely to be why it failed, and we
should not design around that hypothesis. Our variants sit at 0.383-0.389 --
essentially the original, since only ~4% of the gene changes.
"""

import csv
import os


_TABLE_PATHS = [
    # the copy VISTA ships with, and the pipeline's own copy
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 "and_gate_pipeline", "data", "ecoli_codon_usage_table.csv"),
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 "vista", "toehold-VISTA", "ecoli_codon_usage_table.csv"),
]

_FRACTIONS = None


def fractions():
    """{codon: fraction} for E. coli, loaded once."""
    global _FRACTIONS
    if _FRACTIONS is None:
        for path in _TABLE_PATHS:
            if os.path.exists(path):
                with open(path, encoding="utf-8") as fh:
                    _FRACTIONS = {row["Codon"].upper().replace("U", "T"):
                                  float(row["Fraction"])
                                  for row in csv.DictReader(fh)}
                break
        else:
            # Not fatal: usage becomes a no-op tie-break and the metric reads 0.
            _FRACTIONS = {}
    return _FRACTIONS


def codon_fraction(codon):
    """Usage fraction of one codon, 0.0 if the table is missing."""
    return fractions().get(codon.upper().replace("U", "T"), 0.0)


def codon_span(start, end):
    """
    Widen a (start, end) window to whole codons.

    Windows in this pipeline routinely begin mid-codon -- k2 is 6-10 nt and the
    frame does not care. Any time DNA and amino acids are shown together they
    must come from the SAME codon-aligned span, or a perfectly good synonymous
    recoding looks like a broken translation.
    """
    return start - (start % 3), end + ((3 - end % 3) % 3)


def mean_fraction(cds, start=None, end=None):
    """
    Mean usage fraction over a window, codon-aligned.

    This is VISTA's "average codon fraction across the target binding region".
    """
    lo, hi = codon_span(0 if start is None else start,
                        len(cds) if end is None else end)
    hi = min(hi, len(cds) - (len(cds) % 3))
    codons = [cds[i:i + 3] for i in range(lo, hi, 3)]
    if not codons:
        return 0.0
    return sum(codon_fraction(c) for c in codons) / len(codons)


def first_two_codons_fraction(cds, start):
    """
    VISTA's "average of the first two codons" of the target region.

    Translation initiation is sensitive to the first few codons in particular,
    so VISTA reports them separately from the region average.
    """
    lo = start - (start % 3)
    codons = [cds[lo:lo + 3], cds[lo + 3:lo + 6]]
    codons = [c for c in codons if len(c) == 3]
    if not codons:
        return 0.0
    return sum(codon_fraction(c) for c in codons) / len(codons)


def report(cds, windows):
    """
    Usage summary for a gene and a set of named windows.

    `windows` is {name: (start, end)}. Returns {name: mean_fraction} plus
    "whole_gene".
    """
    out = {"whole_gene": mean_fraction(cds)}
    for name, (start, end) in windows.items():
        out[name] = mean_fraction(cds, start, end)
    return out
