"""
Every metric reported in Toehold_Candidates29.7.pdf, in one place.

Stage 0 uses this to check our tooling against the document; stage 4 will use
the same functions on the AND-gate designs, so the two are guaranteed to be
measuring the same things.

Each function's docstring says exactly what is computed and whether it
reproduces the PDF. See METRICS.md in this folder for the plain-language
version and the reproduction status of every row.

Reproduction status, briefly:
  exact       11 of the PDF's rows, to the last decimal it prints
  advisory     2 rows (the RBS pair) -- the PDF does not state its averaging
               window, and no window we tried matches exactly
  unresolved   1 row ("region just after the start codon open - ON")
  not wired    1 row (the VISTA similarity score)
"""

# Make relative imports work when this file is run on its own (the Run button in
# Visual Studio executes it as a plain script, with no package context). Runs
# only in that case; a normal "import poc_and.x" skips it entirely.
if __package__ in (None, ""):
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    import poc_and  # noqa: F401  -- makes the parent package real
    __package__ = "poc_and"

from .folding import RNA   # via folding, so the setup check runs first

from . import candidates as cd
from . import codon_usage
from . import folding


STOP_CODONS = ("UAA", "UAG", "UGA")
START_CODON = "AUG"


# ---------------------------------------------------------------------------
# Energies
# ---------------------------------------------------------------------------

def dG_switch(switch):
    """"dG switch alone". Reproduces the PDF exactly."""
    return RNA.fold(switch)[1]


def dG_trigger(trigger):
    """The trigger's own folding energy. The PDF does not print this, but it
    is needed for the corrected margin below."""
    return RNA.fold(trigger)[1]


def dG_complex(trigger, switch):
    """"dG switch + trigger complex". Reproduces the PDF exactly."""
    return RNA.cofold(trigger + "&" + switch)[1]


def dG_margin(trigger, switch):
    """
    "dG margin (net drive to bind)", as the PDF computes it.

    complex - switch. Note this leaves out the trigger's own folding, so it
    overstates the drive by 1.3 to 6.2 kcal/mol across the five candidates.
    Reproduces the PDF exactly.
    """
    return dG_complex(trigger, switch) - dG_switch(switch)


def dG_margin_corrected(trigger, switch):
    """
    The margin with the trigger's self-folding added back.

    complex - switch - trigger. This is the honest net drive: the trigger has
    to be unfolded before it can bind, and that costs energy.

    Worth knowing: the PDF itself uses THIS quantity as the denominator of its
    off-target percentage, while printing the uncorrected one in the dG margin
    row. The two rows in the document are not on the same footing.
    """
    return (dG_complex(trigger, switch)
            - dG_switch(switch)
            - dG_trigger(trigger))


# ---------------------------------------------------------------------------
# Accessibility
# ---------------------------------------------------------------------------

def toehold_open_off(switch, spans, temperature_c=37.0):
    """
    "Toehold open - OFF (trigger landing site)".

    Mean probability that a toehold base is unpaired, with the switch folded on
    its own. Reproduces the PDF exactly.
    """
    unpaired = folding.unpaired_probs(switch, temperature_c)
    return folding.mean_unpaired(unpaired, spans["toehold"])


def rbs_accessibility(switch, trigger, spans, temperature_c=37.0):
    """
    The two RBS rows: "RBS hidden - OFF" and "RBS exposed - ON".

    Both are the mean probability that an RBS base is UNPAIRED -- in the OFF
    state and in the switch+trigger complex respectively.

    The OFF row is mislabelled in the PDF. A high number there means the RBS is
    accessible, not protected. That is expected for this architecture, because
    the RBS sits in a hairpin LOOP and is unpaired by design; repression comes
    from the start codon being sequestered. The practical consequence is that
    this row barely changes between OFF and ON and so cannot rank candidates.

    Advisory, not exact: the PDF does not say which window it averaged over.
    We use the full 11-nt AACAGAGGAGA, which lands within ~0.6 points.
    """
    off = folding.unpaired_probs(switch, temperature_c)
    on = folding.unpaired_probs([trigger, switch], temperature_c)[len(trigger):]
    return (folding.mean_unpaired(off, spans["rbs"]),
            folding.mean_unpaired(on, spans["rbs"]))


def aug_accessibility(switch, trigger, spans, temperature_c=37.0):
    """
    "Start codon (AUG) exposed - ON", plus the OFF value the PDF omits.

    Mean probability that a start-codon base is unpaired. The ON value
    reproduces the PDF exactly.

    This is the row that actually discriminates: the start codon moves 13-28
    points between OFF and ON, where the RBS moves 0-7. Rank on this one.
    """
    off = folding.unpaired_probs(switch, temperature_c)
    on = folding.unpaired_probs([trigger, switch], temperature_c)[len(trigger):]
    return (folding.mean_unpaired(off, spans["aug"]),
            folding.mean_unpaired(on, spans["aug"]))


def aug_bulge_paired_off(switch, spans):
    """
    "AUG bulge in OFF (MFE structure)".

    How many of the start codon's three bases are paired in the single most
    stable structure. The PDF calls 0/3 "clean" and 1/3 "acceptable".
    Reproduces the PDF exactly (0,0,0,1,1).

    Read with care: this looks at the MFE structure only, and the MFE accounts
    for as little as 2% of the ensemble for these molecules. The ensemble
    columns in our report exist to put that in context.
    """
    structure, _ = RNA.fold(switch)
    pair_table = RNA.ptable(structure)
    start, end = spans["aug"]
    return sum(1 for i in range(start, end) if pair_table[i + 1] != 0)


def trigger_unfolded(trigger, temperature_c=37.0):
    """
    "Trigger unfolded on its own".

    Mean probability that a trigger base is unpaired when the trigger is folded
    by itself. A trigger tangled in its own structure has to pay to open before
    it can invade. Reproduces the PDF exactly.
    """
    unpaired = folding.unpaired_probs(trigger, temperature_c)
    return 100.0 * sum(unpaired) / len(unpaired)


def trigger_binding(trigger, switch, temperature_c=37.0):
    """
    "Trigger binding (at equilibrium)".

    Mean probability that a trigger base is paired TO THE SWITCH in the
    two-strand ensemble -- i.e. how completely the trigger lands on its target
    rather than staying free or folding on itself. Reproduces the PDF exactly.
    """
    n_trigger = len(trigger)
    total = n_trigger + len(switch)
    fc = RNA.fold_compound(trigger + "&" + switch)
    _, mfe_energy = fc.mfe()
    fc.exp_params_rescale(mfe_energy)
    fc.pf()
    bpp = fc.bpp()

    paired_to_switch = 0.0
    for i in range(1, total + 1):
        row = bpp[i]
        for j in range(i + 1, total + 1):
            p = row[j]
            if not p:
                continue
            # count the pair once for each partner that lives on the trigger
            if i <= n_trigger < j:
                paired_to_switch += p
            elif j <= n_trigger:
                pass          # both ends on the trigger: self-structure, not binding
    return 100.0 * paired_to_switch / n_trigger


def after_start_open_on(switch, trigger, spans, temperature_c=37.0):
    """
    "Region just after the start codon open - ON".  UNRESOLVED.

    We could not reproduce the PDF's values (3.4, 0.6, 1.7, 1.8, 4.6) from its
    description. Definitions tried and rejected:

      * mean unpaired probability of a window after the AUG, ON state
        -> comes out 88-95%, far too high
      * mean PAIRED probability of the same window (lengths 3 to 21)
        -> best fit at length 6, total error 16 points across five candidates
      * joint probability that the whole window is unpaired, via a constrained
        partition function (lengths 3 to 15)
        -> no length reproduces the pattern
      * the same on the b_pre* span, the linker span, and in the OFF state
        -> none match

    Returned here as the mean paired probability over the 6 nt following the
    start codon, which is the closest of the candidates tried, clearly marked
    so it is never mistaken for a reproduced value. To settle it we need the
    definition from whoever wrote the original validation script.
    """
    on = folding.unpaired_probs([trigger, switch], temperature_c)[len(trigger):]
    start = spans["aug"][1]
    end = min(start + 6, len(switch))
    return 100.0 - folding.mean_unpaired(on, (start, end))


# ---------------------------------------------------------------------------
# Off-target
# ---------------------------------------------------------------------------

def offtarget(trigger, reporter_rna):
    """
    "Off-target vs. engineered mCherry (duplex energy)" and its percentage.

    The best duplex between the TRIGGER and the engineered mCherry reporter --
    asking whether the trigger would rather stick to the reporter transcript
    than to its switch. Reproduces the PDF exactly (-23.3, -26.5, -25.0,
    -26.5, -24.5).

    The percentage the PDF prints alongside is this energy over the CORRECTED
    margin (complex - switch - trigger), not the margin it prints in the table.
    That reproduces exactly too: 52, 52, 49, 49, 43.

    Returns (duplex_dG, percent_of_intended). Percent needs the caller to pass
    the corrected margin, so it is computed in pdf_table() below.
    """
    return RNA.duplexfold(trigger, reporter_rna).energy


# ---------------------------------------------------------------------------
# Sequence integrity
# ---------------------------------------------------------------------------

def rbs_present(switch):
    """"RBS present". Is the Shine-Dalgarno sequence actually there."""
    return cd.RBS in switch


def reading_frame_intact(switch, spans):
    """
    "Reading frame intact".

    The distance from the start codon to the end of the construct must be a
    multiple of three, or the reporter downstream is translated out of frame.
    """
    aug_start = spans["aug"][0]
    return (len(switch) - aug_start) % 3 == 0


def unwanted_codons(switch, spans):
    """
    "Unwanted start / stop codons".

    Walking in frame from the real start codon to the end of the construct:
    any further AUG could start a competing product, and any stop codon would
    truncate the reporter before it begins. The PDF reports 0 / 0 for all five,
    which we reproduce.
    """
    aug_start = spans["aug"][0]
    tail = switch[aug_start:]
    codons = [tail[i:i + 3] for i in range(0, len(tail) - 2, 3)]
    extra_starts = sum(1 for c in codons[1:] if c == START_CODON)
    stops = sum(1 for c in codons if c in STOP_CODONS)
    return extra_starts, stops


# ---------------------------------------------------------------------------
# Ensemble defects and codon usage -- the extra panel Green 2026 VISTA reports
# ---------------------------------------------------------------------------

def ensemble_defects(sequence, temperature_c=37.0):
    """
    SED and NED for one molecule. ViennaRNA computes both directly.

    SED (specified ensemble defect, against the fully-unpaired reference) is
    the average number of nucleotides per base that are NOT single-stranded, so
    it reads as "how far from open is this molecule". Lower = more accessible.
    For a trigger that has to invade a hairpin, low SED is good.

    NED (native ensemble defect, against the molecule's own MFE structure) says
    how representative that MFE structure is of the whole ensemble. Lower = the
    MFE is a fair picture. This is the quantitative version of the warning that
    keeps recurring in this project: candidate 4's MFE carries only 2.2% of the
    ensemble, so single-structure claims about it are weak.

    Both are normalised to 0..1 by ViennaRNA, so they compare across lengths.
    """
    fc = RNA.fold_compound(sequence)
    mfe_structure, mfe_energy = fc.mfe()
    fc.exp_params_rescale(mfe_energy)
    fc.pf()
    return {
        "SED": fc.ensemble_defect("." * len(sequence)),
        "NED": fc.ensemble_defect(mfe_structure),
    }


def codon_metrics(gene_dna, region_start, region_end):
    """
    VISTA's codon-usage panel for a trigger's binding region.

    Three numbers, following compute_codon_fractions in the VISTA notebook:
      region_mean   average E. coli usage fraction across the region
      first_two     average over the first two codons, which matter
                    disproportionately for translation initiation
      whole_gene    the same average over the entire transcript, for context

    Usage is the share of an amino acid's codons that are this one in E. coli,
    so 0.5 is unremarkable and 0.05 means the organism almost never uses it.

    Context for reading these: the codon-max mCherry that FAILED at the bench
    scores 0.471, better than the working original's 0.384. Poor codon usage is
    therefore unlikely to explain that failure, and a low number here is a flag
    to investigate rather than a verdict.
    """
    return {
        "codon_region_mean": codon_usage.mean_fraction(
            gene_dna, region_start, region_end),
        "codon_first_two": codon_usage.first_two_codons_fraction(
            gene_dna, region_start),
        "codon_whole_gene": codon_usage.mean_fraction(gene_dna),
    }


# ---------------------------------------------------------------------------
# The whole table for one candidate
# ---------------------------------------------------------------------------

def pdf_table(cand_id, reporter_rna=None, temperature_c=37.0, trigger_seq=None):
    """
    Compute every PDF row for one candidate.

    reporter_rna is the engineered mCherry the off-target check runs against.
    Pass None to skip that row. trigger_seq overrides the printed trigger, so
    the same table can be produced for the real-mCherry version.
    """
    c = cd.CANDIDATES[cand_id]
    switch = c["switch"]
    trigger = trigger_seq if trigger_seq is not None else c["trigger"]
    spans = cd.domains(cand_id)

    g_switch = dG_switch(switch)
    g_trigger = dG_trigger(trigger)
    g_complex = dG_complex(trigger, switch)
    margin = g_complex - g_switch
    margin_corrected = margin - g_trigger

    rbs_off, rbs_on = rbs_accessibility(switch, trigger, spans, temperature_c)
    aug_off, aug_on = aug_accessibility(switch, trigger, spans, temperature_c)
    extra_starts, stops = unwanted_codons(switch, spans)

    row = {
        "cand": cand_id,
        "family": c["family"],
        "trigger_len": len(trigger),
        "k1_len": c["k1_len"],

        "dG_switch": g_switch,
        "dG_trigger": g_trigger,
        "dG_complex": g_complex,
        "dG_margin": margin,
        "dG_margin_corrected": margin_corrected,
        "dG_margin_per_nt": margin / len(trigger),

        "toehold_open_off": toehold_open_off(switch, spans, temperature_c),
        "rbs_off": rbs_off,
        "rbs_on": rbs_on,
        "aug_off": aug_off,
        "aug_on": aug_on,
        "aug_gap_on_minus_off": aug_on - aug_off,
        "aug_bulge_paired_off": aug_bulge_paired_off(switch, spans),

        "trigger_binding": trigger_binding(trigger, switch, temperature_c),
        "trigger_unfolded": trigger_unfolded(trigger, temperature_c),

        "rbs_present": rbs_present(switch),
        "reading_frame_intact": reading_frame_intact(switch, spans),
        "unwanted_starts": extra_starts,
        "unwanted_stops": stops,

        # unresolved -- see the docstring on after_start_open_on
        "after_start_on_UNRESOLVED": after_start_open_on(
            switch, trigger, spans, temperature_c),
    }

    switch_defects = ensemble_defects(switch, temperature_c)
    trigger_defects = ensemble_defects(trigger, temperature_c)
    row["switch_SED"] = switch_defects["SED"]
    row["switch_NED"] = switch_defects["NED"]
    row["trigger_SED"] = trigger_defects["SED"]
    row["trigger_NED"] = trigger_defects["NED"]

    if reporter_rna:
        ot = offtarget(trigger, reporter_rna)
        row["offtarget_dG"] = ot
        row["offtarget_pct_of_intended"] = 100.0 * abs(ot) / abs(margin_corrected)
    else:
        row["offtarget_dG"] = None
        row["offtarget_pct_of_intended"] = None

    return row


# Pressing Run on this file alone prints the metric table for every candidate.
if __name__ == "__main__":
    for _cid in sorted(cd.CANDIDATES):
        _row = pdf_table(_cid)
        print("cand %d  dG_switch %7.2f  AUG on %5.1f%%  SED %.3f  NED %.3f"
              % (_cid, _row["dG_switch"], _row["aug_on"],
                 _row["switch_SED"], _row["switch_NED"]))
