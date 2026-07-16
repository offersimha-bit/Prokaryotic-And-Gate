"""Small bundled example inputs for a self-contained demo run.

GENE1 and GENE2 contain a planted, *imperfect* reverse-complement pair so the
demo exercises both the exact-match and the minimum-Hamming fallback paths, and
a tiny mock transcriptome for the off-target scan.  Real runs should supply
genuine gene sequences.
"""

# 12-nt core planted in GENE1 as "x"
_X = "GCAUACGGAUCA"

GENE1 = (
    "AUGGCACGUUAACCGGAUUCCAUGCAUACAGG"      # 5' context (provides r1)
    + _X                                     # x
    + "UGACAUGGCA"                            # a + k1 context
    + "CCGUUAACGGAUUCCGAUUACGCAUGGCACGUUAAUACGGACAU"
)

# GENE2 carries an approximate reverse complement of _X (one mismatch) so the
# Hamming fallback is used.
_K2_APPROX = "UGAUCCGUAUGG"   # revcomp(_X) = UGAUCCGUAUGC ; last base differs
GENE2 = (
    "GGCAUUAACGGGAUUCCAUUACGGCACAUUGGCAUAA"   # 5' context (provides r2)
    + _K2_APPROX                              # k2 (approx)
    + "CGGUUAACCGGAUUCCAUGCAUUACGGCACAUUAAGGCAU"
)

# mock transcriptome for the off-target scan
TRANSCRIPTOME = {
    "gene1": GENE1,
    "gene2": GENE2,
    "housekeeping_rpoB": "AUGGCUAGCUAGCUAGCAUCGAUCGUAGCUAGCUAGCAUCGAUCGAUCGUAG",
    "essential_ftsZ": "AUGUUCGAUCCGAUCGAUCGAUUUCCGGAUUAACGGCAUUAACCGGAUUCCU",
}
ESSENTIAL = {"essential_ftsZ"}
