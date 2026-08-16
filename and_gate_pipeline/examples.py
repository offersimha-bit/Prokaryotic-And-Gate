"""Bundled inputs for a self-contained run.

The default gene pair is REAL: the *E. coli* K-12 MG1655 threonine-operon genes
``thrA`` and ``thrC``, shipped as FASTA under ``data/genes/``.  Everything --
demo runs, sweeps, the truth table -- should be measured on these unless a run
explicitly supplies its own genes.

Why the switch away from the synthetic pair matters
---------------------------------------------------
The old ``GENE1``/``GENE2`` were 98 and 89 nt.  Stage 2's RNAplfold window is
240 nt, so on sequences that short the window never restricts anything and the
"local" accessibility silently degenerated to a global fold -- the pipeline
could not exercise the very model stage 2 exists to apply.  thrA (2463 nt) and
thrC (1287 nt) are longer than the window, so local folding is genuinely local
and a trigger's flanking context is real rather than an artefact of a stub.

They are also a defensible biological pair rather than an arbitrary one: both
sit in the same operon, are co-transcribed, and are therefore co-expressed --
which is the precondition an AND gate on endogenous transcripts needs.  Note
that co-expression is exactly what makes them a convenient TEST pair and a poor
DISCRIMINATING pair: a real sensor wants two genes that are high together only
in the disease state.  These are for exercising the pipeline, not a design
proposal.

The synthetic pair is kept as ``SYNTHETIC_GENE1``/``SYNTHETIC_GENE2`` because
several tests rely on its planted, deliberately imperfect reverse complement to
exercise the minimum-Hamming fallback path.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

_GENE_DIR = Path(__file__).resolve().parent / "data" / "genes"


def read_fasta(path) -> str:
    """Sequence of a single-record FASTA, uppercased, non-letters stripped."""
    lines = Path(path).read_text().splitlines()
    return "".join(ln.strip() for ln in lines if not ln.startswith(">")).upper()


@lru_cache(maxsize=None)
def load_gene(name: str) -> str:
    """Load a bundled gene by bare name, e.g. ``load_gene("thrA")``.

    Returned as DNA exactly as stored; callers convert with
    ``sequence_utils.to_rna`` the same way they would for user input, so the
    bundled genes travel the identical code path as a real FASTA.
    """
    path = _GENE_DIR / f"{name}.fa"
    if not path.exists():
        available = sorted(p.stem for p in _GENE_DIR.glob("*.fa"))
        raise FileNotFoundError(f"no bundled gene {name!r}; have {available}")
    return read_fasta(path)


def gene_path(name: str) -> Path:
    return _GENE_DIR / f"{name}.fa"


# --------------------------------------------------------------------------- #
# the default pair                                                            #
# --------------------------------------------------------------------------- #
GENE1_NAME = "thrA"     # aspartate kinase / homoserine dehydrogenase 1, 2463 nt
GENE2_NAME = "thrC"     # threonine synthase, 1287 nt

GENE1 = load_gene(GENE1_NAME)
GENE2 = load_gene(GENE2_NAME)


# --------------------------------------------------------------------------- #
# synthetic pair -- retained for the tests that need a planted near-match     #
# --------------------------------------------------------------------------- #
_X = "GCAUACGGAUCA"

SYNTHETIC_GENE1 = (
    "AUGGCACGUUAACCGGAUUCCAUGCAUACAGG"      # 5' context (provides r1)
    + _X                                     # x
    + "UGACAUGGCA"                            # a + k1 context
    + "CCGUUAACGGAUUCCGAUUACGCAUGGCACGUUAAUACGGACAU"
)

# carries an approximate reverse complement of _X (one mismatch) so the
# minimum-Hamming fallback is exercised.
_K2_APPROX = "UGAUCCGUAUGG"   # revcomp(_X) = UGAUCCGUAUGC ; last base differs
SYNTHETIC_GENE2 = (
    "GGCAUUAACGGGAUUCCAUUACGGCACAUUGGCAUAA"   # 5' context (provides r2)
    + _K2_APPROX                              # k2 (approx)
    + "CGGUUAACCGGAUUCCAUGCAUUACGGCACAUUAAGGCAU"
)


# --------------------------------------------------------------------------- #
# mock transcriptome for the off-target scan                                  #
# --------------------------------------------------------------------------- #
# Still a stub: four sequences cannot stand in for an E. coli transcriptome, and
# the energetic off-target scan additionally wants an FPKM table it does not yet
# have.  Off-target numbers from this input are illustrative, not measurements.
TRANSCRIPTOME = {
    GENE1_NAME: GENE1,
    GENE2_NAME: GENE2,
    "housekeeping_rpoB": "AUGGCUAGCUAGCUAGCAUCGAUCGUAGCUAGCUAGCAUCGAUCGAUCGUAG",
    "essential_ftsZ": "AUGUUCGAUCCGAUCGAUCGAUUUCCGGAUUAACGGCAUUAACCGGAUUCCU",
}
ESSENTIAL = {"essential_ftsZ"}
