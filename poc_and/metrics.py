"""
Stage 4 -- score each AND switch across all four input states, and rank.

THE FOUR STATES
---------------
    00   switch alone                    should be OFF
    10   switch + trigger A              should be OFF  (B has not unmasked x*)
    01   switch + trigger B              should be OFF  (no A to fire it)
    11   switch + trigger A + trigger B  should be ON

The ON state is a THREE-strand complex. ViennaRNA 2.7 handles that directly
(`fold_compound("A&B&switch")`), verified, so the whole table runs on Windows
and NUPACK stays a second opinion rather than a requirement.


THE HEADLINE NUMBER
-------------------
How much better trigger A can open the switch once trigger B has done its job:

    dG_open(A)    = G(A+switch)   - G(switch)   - G(A)
    dG_open(A|B)  = G(A+B+switch) - G(B+switch) - G(A)
    ddG_AND       = dG_open(A|B)  - dG_open(A)
                  = (G11 - G01) - (G10 - G00)          <- G(A) cancels

More negative means the inhibitory hairpin is doing its job: A binds much
better with B present than without. Reported in kcal/mol and as the equivalent
fold-change, `exp(-ddG/RT)` -- the same number, since dG is logarithmic, which
answers "should it be a difference or a ratio" with "those are the same thing".

`ddG_AND` needs no length normalisation: it is a difference between two states
of the SAME molecule, so length cancels. That is a large part of why it is the
headline rather than any raw dG, which as stage 0 showed tracks trigger length
(r = -0.872) more than design quality.


THE READOUT
-----------
Energies say how tightly things bind; they do not say whether a ribosome can
start. For that we use START-CODON ACCESSIBILITY, which stage 0 established is
the only accessibility row that actually discriminates (it moves 13-28 points
between OFF and ON where the RBS moves 0-7). Separation is reported
threshold-free: state 11 against the WORST of the three OFF states, so the
answer does not depend on where anyone draws a line.


WHAT THIS STAGE DELIBERATELY DOES NOT DO
----------------------------------------
The triggers here are their binding windows plus a configurable margin of
native context, not the whole 711-nt transcripts: a three-strand fold of
711+711+160 nt is not tractable, and would mostly measure long-RNA noise.
Whole-gene interactions are covered separately, and already were -- stage 2's
`ortho` and `a_kill` both scan the entire gene. Read the two together.
"""

# Make relative imports work when this file is run on its own.
if __package__ in (None, ""):
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    import poc_and  # noqa: F401
    __package__ = "poc_and"

import csv
import math
import os

from .folding import RNA

from . import candidates as cd
from . import codon_usage
from . import folding
from . import pdf_metrics


STATES = ("00", "10", "01", "11")


def trigger_with_context(gene_dna, start, end, context_nt):
    """
    A trigger as the cell actually presents it: its window plus flanking gene.

    A trigger is not a free 30-mer -- it sits inside a transcript, and the
    flanking sequence can fold back over it. Taking some context in is the
    cheap half of that correction; `trigger_unfolded` in stage 0 is the other.
    """
    lo = max(0, start - context_nt)
    hi = min(len(gene_dna), end + context_nt)
    return cd.to_rna(gene_dna[lo:hi]), start - lo


def state_strands(built, trigger_a, trigger_b):
    """
    The strands present in each state, switch always LAST.

    Keeping the switch last means its coordinates are offset by exactly the
    combined length of the triggers ahead of it, which is what
    `switch_offset` below relies on.
    """
    switch = built["sequence"]
    return {
        "00": [switch],
        "10": [trigger_a, switch],
        "01": [trigger_b, switch],
        "11": [trigger_a, trigger_b, switch],
    }


def score_states(built, trigger_a, trigger_b, temperature_c=37.0):
    """
    Fold every state and read the switch's accessibility inside each.

    Returns per-state energy plus the accessibility of the start codon, the
    RBS, and the three toehold sub-domains.
    """
    spans = built["spans"]
    strands = state_strands(built, trigger_a, trigger_b)
    out = {}
    for state, parts in strands.items():
        offset = sum(len(p) for p in parts[:-1])
        ensemble = folding.ensemble(parts, temperature_c)
        unpaired = ensemble["unpaired"][offset:]

        def openness(name):
            return folding.mean_unpaired(unpaired, spans[name])

        out[state] = {
            "dG": ensemble["mfe_dG"],
            "ensemble_dG": ensemble["ensemble_dG"],
            "p_mfe": ensemble["p_mfe"],
            "aug_open": openness("primary_aug"),
            "rbs_open": openness("primary_rbs"),
            "a_star_open": openness("a_star"),
            "x_star_open": openness("x_star"),
            "r2_star_open": openness("r2_star"),
        }
    return out


def and_discrimination(states, temperature_c=37.0):
    """
    ddG_AND and the equivalent fold-change.

    G(A) cancels out of the difference, so this needs only the four state
    energies -- which also means it is insensitive to how the triggers were
    trimmed, a useful property given the context caveat above.
    """
    g = {s: states[s]["dG"] for s in STATES}
    open_a = g["10"] - g["00"]
    open_a_given_b = g["11"] - g["01"]
    ddg = open_a_given_b - open_a
    rt = folding.kt(temperature_c)
    return {
        "dG_open_A": open_a,
        "dG_open_A_given_B": open_a_given_b,
        "ddG_AND": ddg,
        "AND_fold_change": math.exp(-ddg / rt) if ddg > -60 else float("inf"),
    }


def separation(states):
    """
    How far state 11 stands above the worst OFF state, threshold-free.

    A confusion matrix needs a cutoff and its answer moves when the cutoff
    moves; the gap between two predicted states does not. Reported both as
    percentage points and as a ratio.
    """
    on = states["11"]["aug_open"]
    worst_off = max(states[s]["aug_open"] for s in ("00", "10", "01"))
    return {
        "aug_on_11": on,
        "aug_worst_off": worst_off,
        "separation_points": on - worst_off,
        "separation_ratio": (on / worst_off) if worst_off > 0 else float("inf"),
        "leakiest_off_state": max(("00", "10", "01"),
                                  key=lambda s: states[s]["aug_open"]),
    }


def confusion(states, threshold):
    """
    The predicted state table scored against an ideal AND gate.

    Only meaningful WITH a threshold, and we have no wet-lab data to set one,
    so this is a presentation aid rather than a result. State 11 is the single
    positive; 00, 10 and 01 are the three negatives. The threshold-free
    separation above is the number to trust.
    """
    calls = {s: states[s]["aug_open"] >= threshold for s in STATES}
    return {
        "threshold": threshold,
        "TP": int(calls["11"]),
        "FN": int(not calls["11"]),
        "FP": sum(int(calls[s]) for s in ("00", "10", "01")),
        "TN": sum(int(not calls[s]) for s in ("00", "10", "01")),
        "calls": "".join("1" if calls[s] else "0" for s in STATES),
    }


def cross_binding(trigger_a, trigger_b, built, trigger_conc_m, temperature_c):
    """
    Do the two inputs prefer each other over the switch?

    Two margins, both of which must be positive:
        risk_B  trigger B's landing site must out-compete trigger A for it
        risk_A  the switch must out-compete trigger B for trigger A
    Plus the absolute sequestered fraction, which is what the m bound in
    stage 2 exists to control.

    IMPORTANT: pass the BARE trigger windows here, never the context-extended
    ones. duplexfold returns the single best duplex anywhere between two
    sequences, so adding 25 nt of unrelated flanking gene to each side finds a
    strong duplex regardless of the design -- measured on candidate 5, the same
    pair goes from -7.0 kcal/mol (0.1% sequestered) at zero context to -34.9
    (100%) at 25 nt. The first number describes the design; the second
    describes the fact that any two long RNAs stick together somewhere. Whole-
    transcript interactions are covered by stage 2's whole-gene scans instead.
    """
    spans = built["spans"]
    switch = built["sequence"]
    # Trigger B's LANDING SITE: the free r2* toehold plus the k2* half of the
    # stem it invades. Not the hairpin's 5' arm -- that is k2* + r1copy, and
    # r1copy is set by trigger A. Naming this "arm" was misleading and got
    # corrected; r2* sits OUTSIDE the stem, which is what makes it a toehold.
    trigB_site = switch[spans["r2_star"][0]:spans["k2_star"][1]]
    footprint = switch[spans["r1_star"][0]:spans["a_star"][1]]

    ab = folding.duplex_dG(trigger_a, trigger_b, temperature_c)
    b_on_site = folding.duplex_dG(trigger_b, trigB_site, temperature_c)
    a_on_switch = folding.duplex_dG(trigger_a, footprint, temperature_c)
    return {
        "trigger_A_B_dG": ab,
        "sequestered_pct": 100.0 * folding.bound_fraction(
            ab, trigger_conc_m, temperature_c),
        "risk_B": ab - b_on_site,
        "risk_A": ab - a_on_switch,
    }


def mechanism(states, built):
    """
    Does the designed mechanism actually run?

    The whole architecture rests on one step: trigger B invades k2* and RELEASES
    x*, growing trigger A's nucleation toehold from |a| to |a| + Lx. That step is
    directly visible as the change in x* accessibility between state 00 and
    state 01, and it is a far more diagnostic number than any complex energy.

    `nucleation_nt` is the expected number of unpaired nucleotides in the region
    trigger A must nucleate on (a* plus x*). It is the quantity the AND gate is
    actually built to switch, and the one the kinetic model consumes.
    """
    len_a = built["len_a"]
    len_x = built["Lx"]

    def nucleation(state):
        return (states[state]["a_star_open"] / 100.0 * len_a
                + states[state]["x_star_open"] / 100.0 * len_x)

    return {
        "x_star_release_pts": states["01"]["x_star_open"] - states["00"]["x_star_open"],
        "nucleation_nt_00": nucleation("00"),
        "nucleation_nt_01": nucleation("01"),
        "nucleation_gain_nt": nucleation("01") - nucleation("00"),
    }


def structural(built, temperature_c=37.0):
    """Ensemble quality of the OFF-state switch on its own."""
    sequence = built["sequence"]
    defects = pdf_metrics.ensemble_defects(sequence, temperature_c)
    fold = built["fold"]
    return {
        "SED": defects["SED"],
        "NED": defects["NED"],
        "p_mfe_off": fold["p_mfe"],
        "intended_agreement_pct": fold["intended_agreement_pct"],
        "mfe_dG_off": fold["mfe_dG"],
        "ensemble_dG_off": fold["ensemble_dG"],
    }


# ---------------------------------------------------------------------------
# One construct
# ---------------------------------------------------------------------------

def score_build(built, original_gene, temperature_c, trigger_conc_m,
                context_nt, threshold, with_nupack):
    """Every stage-4 number for one built AND switch."""
    design = built["design"]
    cand_id = built["cand"]

    real = cd.trigger_real(cand_id, original_gene)
    a_lo = real["gene_pos"]
    a_hi = a_lo + len(real["sequence"])
    trigger_a, _ = trigger_with_context(original_gene, a_lo, a_hi, context_nt)

    b_lo = design["k2_span"][0]
    b_hi = design["r2_span"][1]
    trigger_b, _ = trigger_with_context(design["variant"], b_lo, b_hi, context_nt)

    states = score_states(built, trigger_a, trigger_b, temperature_c)
    row = {
        "cand": cand_id,
        "len_a": built["len_a"],
        "Lx": built["Lx"],
        "m": built["m"],
        "bulge_size": built["bulge_size"],
        "bulge_index": built["bulge_index"],
        "longest_duplex_bp": built["checks"]["bulge"]["longest_bp"],
        "total_nt": len(built["sequence"]),
        "trigger_A_nt": len(trigger_a),
        "trigger_B_nt": len(trigger_b),
    }
    for state in STATES:
        for key, value in states[state].items():
            row["%s_%s" % (state, key)] = value

    row.update(and_discrimination(states, temperature_c))
    row.update(separation(states))
    row.update(confusion(states, threshold))
    row.update(mechanism(states, built))
    # bare windows, deliberately -- see the note in cross_binding
    bare_a = cd.to_rna(original_gene[a_lo:a_hi])
    bare_b = cd.to_rna(design["variant"][b_lo:b_hi])
    row.update(cross_binding(bare_a, bare_b, built,
                             trigger_conc_m, temperature_c))
    row.update(structural(built, temperature_c))

    # carried through from stage 2, reported but not scored
    row["total_edits"] = design["total_edits"]
    row["gene_changed_pct"] = 100.0 * design["total_edits"] / len(design["variant"])
    row["usage_whole_gene"] = design.get("usage_whole_gene", 0.0)
    row["orthogonality_margin"] = design.get("orthogonality_margin", 0.0)
    row["a_kill"] = design.get("a_kill", 0.0)

    # per-nucleotide companion, so trigger-length effects stay visible
    row["dG_open_A_per_nt"] = row["dG_open_A"] / len(trigger_a)

    if with_nupack:
        row["nupack_dG_11"] = folding.nupack_mfe(
            [trigger_a, trigger_b, built["sequence"]], temperature_c)
        row["nupack_dG_00"] = folding.nupack_mfe(
            built["sequence"], temperature_c)
    else:
        row["nupack_dG_11"] = None
        row["nupack_dG_00"] = None

    row["_built"] = built
    return row


def rank(rows):
    """
    Order the designs.

    Equilibrium first, per the steer that the professor wants equilibrium
    energies. Cross-binding is a gate rather than a score -- a design whose two
    inputs sequester each other cannot work no matter how good its energies
    look. No pass/fail threshold beyond that: the table is ranked and the cut is
    yours to make.
    """
    return sorted(rows, key=lambda r: (
        # --- gates: a design failing either of these cannot work at all -----
        r["sequestered_pct"] > 10.0,
        not (r["risk_A"] > 0 and r["risk_B"] > 0),

        # --- the mechanism, and it IS an equilibrium quantity ---------------
        # Ranking on ddG_AND alone is misleading here, and the mechanism table
        # shows why: candidates 1 and 3 post the most negative ddG_AND of the
        # set while trigger B releases x* by under 1 point -- their gate does
        # not actually run. Nucleation gain measures the one step the whole
        # architecture depends on (B growing trigger A's landing site from |a|
        # to |a|+Lx), it is computed from the same equilibrium ensembles, and
        # unlike ddG_AND it is not confounded by state 10, where trigger A
        # thermodynamically displaces the hairpin no matter what we do.
        -r["nucleation_gain_nt"],
        -r["x_star_release_pts"],

        # --- then the energy headline, then the readout ---------------------
        r["ddG_AND"],
        -r["separation_points"],
    ))


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------

def run(config, results=None):
    """Stage 4 entry point."""
    print()
    print("STAGE 4 -- four-state equilibrium scoring")

    builds = (results or {}).get(3)
    if not builds:
        print("  (stage 3 has not run in this session -- running it now)")
        from . import build
        builds = build.run(config, results)

    mcherry = cd.read_mcherry(config["mcherry_file"])
    original = mcherry["original"]
    context_nt = config.get("trigger_context_nt", 25)
    threshold = config.get("on_threshold_pct", 80.0)
    with_nupack = folding.HAVE_NUPACK and config.get("engine", "auto") != "vienna"

    print("  triggers: binding window +/- %d nt of native context" % context_nt)
    print("  engines : ViennaRNA%s" % (" + NUPACK" if with_nupack else " only"))
    print("  scoring %d built construct(s)..." % len(builds))

    rows = []
    for built in builds:
        rows.append(score_build(built, original, config["temperature_C"],
                                config["trigger_conc_M"], context_nt,
                                threshold, with_nupack))

    rows = rank(rows)
    path = write_outputs(rows, config["out_dir"])
    _report(rows, threshold)
    print()
    print("wrote %s" % path)
    return rows


def _report(rows, threshold):
    print()
    print("=" * 108)
    print("FOUR-STATE TABLE -- start-codon accessibility per state (the ON readout)")
    print("=" * 108)
    print("%5s %6s %5s %4s %8s %8s %8s %8s %10s %10s" %
          ("cand", "bulge", "|a|", "Lx", "00 OFF", "10 OFF", "01 OFF", "11 ON",
           "separation", "ratio"))
    for r in rows:
        print("%5d %6s %5d %4d %7.1f%% %7.1f%% %7.1f%% %7.1f%% %9.1f  %9.2fx" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            r["len_a"], r["Lx"],
            r["00_aug_open"], r["10_aug_open"], r["01_aug_open"], r["11_aug_open"],
            r["separation_points"], r["separation_ratio"]))

    print()
    print("=" * 108)
    print("EQUILIBRIUM RANKING -- ddG_AND is the headline")
    print("=" * 108)
    print("%5s %6s %10s %12s %12s %10s %9s %9s %9s" %
          ("cand", "bulge", "ddG_AND", "fold-change", "dG_open(A)",
           "dG(A|B)", "seq'd%", "risk_A", "risk_B"))
    for r in rows:
        fold = r["AND_fold_change"]
        fold_text = "%.3g" % fold if fold < 1e6 else "%.1e" % fold
        print("%5d %6s %10.2f %12s %12.2f %10.2f %8.1f%% %9.1f %9.1f" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            r["ddG_AND"], fold_text, r["dG_open_A"], r["dG_open_A_given_B"],
            r["sequestered_pct"], r["risk_A"], r["risk_B"]))
    print()
    print("ddG_AND = (G11 - G01) - (G10 - G00); more negative = trigger A binds")
    print("          much better once B has unmasked x*. Needs no length")
    print("          normalisation -- it is a difference within one molecule.")
    print("risk_A / risk_B must both be positive: the switch has to out-compete")
    print("          the other trigger. Negative means the inputs prefer each other.")

    print()
    print("Confusion matrix at an ON threshold of %.0f%% start-codon accessibility."
          % threshold)
    print("We have no wet-lab data to set that threshold, so this is a")
    print("presentation aid -- the separation column above is threshold-free and")
    print("is the number to trust.")
    print()
    print("%5s %6s %10s %5s %5s %5s %5s" %
          ("cand", "bulge", "00/10/01/11", "TP", "TN", "FP", "FN"))
    for r in rows:
        print("%5d %6s %10s %5d %5d %5d %5d" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            r["calls"], r["TP"], r["TN"], r["FP"], r["FN"]))

    print()
    print("=" * 108)
    print("MECHANISM CHECK -- does trigger B actually release x*?")
    print("=" * 108)
    print("The entire architecture rests on one step: B invades k2* and frees")
    print("x*, growing trigger A's nucleation toehold from |a| to |a|+Lx.")
    print("That step is directly visible below, and it is more diagnostic than")
    print("any complex energy.")
    print()
    print("%5s %6s %4s %4s %10s %10s %12s %12s %11s" %
          ("cand", "bulge", "|a|", "Lx", "x* 00", "x* 01", "release",
           "nucleation", "gain"))
    for r in rows:
        print("%5d %6s %4d %4d %9.1f%% %9.1f%% %11.1f %6.1f->%-5.1f %10.1f" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            r["len_a"], r["Lx"], r["00_x_star_open"], r["01_x_star_open"],
            r["x_star_release_pts"], r["nucleation_nt_00"],
            r["nucleation_nt_01"], r["nucleation_gain_nt"]))
    print()
    print("nucleation = expected unpaired nt in (a* + x*), the site trigger A")
    print("             has to land on. The gate works by SWITCHING this number.")

    print()
    print("!" * 108)
    print("READ THIS BEFORE RANKING ON THE EQUILIBRIUM NUMBERS")
    print("!" * 108)
    print("Equilibrium cannot show AND behaviour for this architecture, and the")
    print("table above is where that becomes visible rather than a caveat.")
    print()
    print("In state 10 trigger A displaces the inhibitory hairpin COMPLETELY")
    print("(inspected on candidate 5: k2* 0% paired, x* 100% paired to A). It")
    print("has to: A forms ~38 bp with the switch while the hairpin holds only")
    print("~24 bp, so A wins on energy alone. And it must win -- a hairpin")
    print("strong enough to resist A at equilibrium could never be opened by A")
    print("even WITH B.")
    print()
    print("So the second input buys a KINETIC barrier, not a thermodynamic one:")
    print("without B, A has only |a| nucleotides to nucleate on, which slows it")
    print("by orders of magnitude against mRNA degradation. That is exactly the")
    print("~700x equilibrium-vs-kinetic disagreement recorded in HANDOFF.md.")
    print()
    print("What the equilibrium numbers here ARE good for: the OFF-state")
    print("structure, whether B releases x* at all, cross-binding, and")
    print("orthogonality. What they CANNOT settle is 10 versus 11. That needs")
    print("the kinetic model (stage 6).")

    print()
    print("Reported, not scored:")
    print("%5s %6s %10s %10s %9s %8s %8s %10s" %
          ("cand", "bulge", "RBS 00", "RBS 11", "SED", "NED", "P(MFE)",
           "gene chg"))
    for r in rows:
        print("%5d %6s %9.1f%% %9.1f%% %8.3f %8.3f %7.2f%% %9.1f%%" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            r["00_rbs_open"], r["11_rbs_open"], r["SED"], r["NED"],
            100 * r["p_mfe_off"], r["gene_changed_pct"]))
    print("RBS barely moves between states -- expected, it sits in a loop by")
    print("design. That is why the start codon is the readout, not the RBS.")


COLUMNS = (["cand", "len_a", "Lx", "m", "bulge_size", "bulge_index",
            "longest_duplex_bp", "total_nt", "trigger_A_nt", "trigger_B_nt"]
           + ["%s_%s" % (s, k) for s in STATES
              for k in ("dG", "ensemble_dG", "p_mfe", "aug_open", "rbs_open",
                        "a_star_open", "x_star_open", "r2_star_open")]
           + ["dG_open_A", "dG_open_A_given_B", "ddG_AND", "AND_fold_change",
              "dG_open_A_per_nt",
              "aug_on_11", "aug_worst_off", "separation_points",
              "separation_ratio", "leakiest_off_state",
              "threshold", "calls", "TP", "TN", "FP", "FN",
              "x_star_release_pts", "nucleation_nt_00", "nucleation_nt_01",
              "nucleation_gain_nt",
              "trigger_A_B_dG", "sequestered_pct", "risk_A", "risk_B",
              "SED", "NED", "p_mfe_off", "intended_agreement_pct",
              "mfe_dG_off", "ensemble_dG_off",
              "nupack_dG_00", "nupack_dG_11",
              "total_edits", "gene_changed_pct", "usage_whole_gene",
              "orthogonality_margin", "a_kill"])


def write_outputs(rows, out_dir):
    """One row per construct, every number, for the Excel deliverable."""
    os.makedirs(out_dir, exist_ok=True)
    from .find_second_trigger import _open_for_write
    fh, path = _open_for_write(os.path.join(out_dir, "stage4_scores.csv"),
                               newline="", encoding="utf-8")
    with fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    return path


# Pressing Run on this file alone does stage 4, using main.py's CONFIG.
if __name__ == "__main__":
    from poc_and.main import CONFIG
    run(dict(CONFIG))
