"""
The folding engine, behind one small interface.

ViennaRNA is the default and is all you need to run the pipeline on Windows.
NUPACK is optional: if it cannot be imported the NUPACK columns come back as
None and the report prints "n/a" instead of crashing. To fill those columns,
re-run the same script under WSL where NUPACK is installed.

Why ViennaRNA is the reference: the numbers in Toehold_Candidates29.7.pdf were
produced with ViennaRNA's default model at 37 C. RNA.fold() and RNA.cofold()
reproduce all ten of the PDF's energies exactly, so anything we want to compare
against the document has to come from the same engine. NUPACK differs by 3-5
kcal/mol on the two-strand complex, which is why it is a second column rather
than a replacement.

Multi-strand note: ViennaRNA 2.7 handles three strands via fold_compound with
two '&' separators, so the four-state AND table does not require NUPACK.
"""

import math
import sys

try:
    import RNA
except ImportError as exc:                       # pragma: no cover - setup help
    # ViennaRNA is a compiled extension, so it is installed per interpreter.
    # The usual cause of landing here is running a DIFFERENT Python from the
    # one it was installed into -- most often the free-threaded build
    # (python3.14t.exe), which keeps its packages in a separate directory and
    # is registered as the default launcher on this machine.
    import site
    try:
        user_site = site.getusersitepackages()
    except Exception:
        user_site = "(unavailable)"
    raise ImportError(
        "ViennaRNA ('RNA') is not available to this interpreter.\n"
        "\n"
        "  running   : %s\n"
        "  python    : %s\n"
        "  user site : %s\n"
        "\n"
        "ViennaRNA is a compiled extension and is installed per interpreter,\n"
        "so a build that did not get it will not see another build's copy.\n"
        "The free-threaded build (python3.14t.exe) is a common culprit: it is\n"
        "often the default launcher but keeps packages in a separate folder.\n"
        "\n"
        "Fix: point your editor at the interpreter that has ViennaRNA, or\n"
        "install it into this one with:\n"
        "    \"%s\" -m pip install ViennaRNA\n"
        % (sys.executable, sys.version.split()[0], user_site, sys.executable)
    ) from exc

try:
    import nupack
    HAVE_NUPACK = True
except Exception:
    nupack = None
    HAVE_NUPACK = False


GAS_CONSTANT = 0.0019872          # kcal / (mol K)


def kt(temperature_c=37.0):
    """RT in kcal/mol at the given temperature."""
    return GAS_CONSTANT * (temperature_c + 273.15)


# ---------------------------------------------------------------------------
# ViennaRNA
# ---------------------------------------------------------------------------

def _joined(strands):
    """['AAA', 'CCC'] -> 'AAA&CCC'. A single string is passed through."""
    if isinstance(strands, str):
        return strands
    return "&".join(strands)


def _compound(strands, temperature_c=37.0):
    md = RNA.md()
    md.temperature = temperature_c
    return RNA.fold_compound(_joined(strands), md)


def mfe(strands, temperature_c=37.0):
    """
    Minimum free energy structure and its energy.

    `strands` is a sequence string, or a list of strings for a complex.
    Returns (structure, dG).
    """
    fc = _compound(strands, temperature_c)
    structure, energy = fc.mfe()
    return structure, energy


def ensemble(strands, temperature_c=37.0):
    """
    Partition-function quantities for one molecule or complex.

    Returns a dict with:
      mfe_structure, mfe_dG      the single most stable structure
      ensemble_dG                free energy of the whole ensemble
      p_mfe                      how much of the ensemble the MFE accounts for
      centroid_structure         the structure closest to the ensemble average
      centroid_dG, p_centroid
      diversity                  mean base-pair distance in the ensemble
      unpaired                   per-position probability of being unpaired

    p_mfe matters more than it looks: for these switches it is often only a few
    percent, which means conclusions drawn from the MFE picture alone can be
    unrepresentative of what the molecule actually does.
    """
    fc = _compound(strands, temperature_c)
    mfe_structure, mfe_dG = fc.mfe()
    fc.exp_params_rescale(mfe_dG)          # keeps the pf numerically stable
    _, ensemble_dG = fc.pf()
    centroid_structure, _ = fc.centroid()

    return {
        "mfe_structure": mfe_structure,
        "mfe_dG": mfe_dG,
        "ensemble_dG": ensemble_dG,
        "p_mfe": fc.pr_structure(mfe_structure),
        "centroid_structure": centroid_structure,
        "centroid_dG": fc.eval_structure(centroid_structure),
        "p_centroid": fc.pr_structure(centroid_structure),
        "diversity": fc.mean_bp_distance(),
        "unpaired": _unpaired_from(fc, len(_joined(strands).replace("&", ""))),
    }


def _unpaired_from(fold_compound, n):
    """
    Per-position probability that a base is NOT paired.

    fc.bpp() is a 1-indexed upper-triangular matrix of pair probabilities, so
    a base's unpaired probability is 1 minus everything it takes part in.
    """
    bpp = fold_compound.bpp()
    unpaired = [1.0] * n
    for i in range(1, n + 1):
        row = bpp[i]
        for j in range(i + 1, n + 1):
            p = row[j]
            if p:
                unpaired[i - 1] -= p
                unpaired[j - 1] -= p
    return unpaired


def unpaired_probs(strands, temperature_c=37.0):
    """Per-position unpaired probability. Convenience wrapper."""
    fc = _compound(strands, temperature_c)
    _, mfe_dG = fc.mfe()
    fc.exp_params_rescale(mfe_dG)
    fc.pf()
    return _unpaired_from(fc, len(_joined(strands).replace("&", "")))


def mean_unpaired(unpaired, span):
    """Average unpaired probability over a (start, end) span, as a percentage."""
    start, end = span
    if end <= start:
        return float("nan")
    return 100.0 * sum(unpaired[start:end]) / (end - start)


def duplex_dG(seq_a, seq_b, temperature_c=37.0):
    """
    Energy of the best duplex between two separate strands.

    Used for the trigger-vs-trigger cross-binding check, where we care about
    the intermolecular helix only and not about either strand's own folding.
    """
    md = RNA.md()
    md.temperature = temperature_c
    return RNA.duplexfold(seq_a, seq_b).energy


def dissociation_constant(dG, temperature_c=37.0):
    """
    Kd in molar, from a binding free energy.

    Kd is a concentration: the concentration at which half the molecules are
    paired up. Comparing it against the actual concentration in the cell is
    what tells you whether two strands will really find each other.
    """
    return math.exp(dG / kt(temperature_c))


def bound_fraction(dG, total_conc_m, temperature_c=37.0):
    """
    Fraction of A tied up as A:B when both start at total_conc_m.

    Solves A + B <-> AB at equilibrium for equal starting concentrations.
    """
    kd = dissociation_constant(dG, temperature_c)
    c = total_conc_m
    if c <= 0:
        return 0.0
    # complex concentration from the quadratic, then express as a fraction
    disc = (2.0 * c + kd) ** 2 - 4.0 * c * c
    complex_conc = (2.0 * c + kd - math.sqrt(disc)) / 2.0
    return complex_conc / c


# ---------------------------------------------------------------------------
# NUPACK -- optional second opinion
# ---------------------------------------------------------------------------

def nupack_mfe(strands, temperature_c=37.0):
    """
    MFE energy from NUPACK, or None if NUPACK is not installed.

    Kept deliberately thin: this exists to give the report a second engine
    column, not to become an alternative code path.
    """
    if not HAVE_NUPACK:
        return None
    model = nupack.Model(material="rna", celsius=temperature_c)
    seqs = [strands] if isinstance(strands, str) else list(strands)
    result = nupack.mfe(strands=seqs, model=model)
    return float(result[0].energy)


def engine_report():
    """One line saying which engines are live, for the run header."""
    parts = ["ViennaRNA %s" % RNA.__version__]
    parts.append("NUPACK %s" % nupack.__version__ if HAVE_NUPACK else "NUPACK not installed")
    return ", ".join(parts)
