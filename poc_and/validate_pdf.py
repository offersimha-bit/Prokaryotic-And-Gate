"""
Stage 0 -- reproduce Toehold_Candidates29.7.pdf, then report what it got wrong.

Run this first, every time. It is the regression test for the whole pipeline:
if the folding engine, its model settings or the stored sequences ever drift,
this stage fails loudly before any design work is built on top of them.

Three parts:

  1. REPRODUCTION. Recompute every row the PDF prints and compare value by
     value. Eleven of them match exactly. Two (the RBS pair) are advisory
     because the PDF does not state its averaging window. One is unresolved.

  2. CORRECTIONS. Three rows are computed correctly but read wrongly, and they
     change how the candidates rank. Each is reprinted with its evidence.

  3. TRIGGER VERSIONS. Where the printed triggers differ from the gene that
     will actually be transcribed.

The metric definitions themselves live in pdf_metrics.py, and METRICS.md in
this folder explains each one in plain language.
"""

import csv
import os
import re

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
from . import folding
from . import pdf_metrics


# Rows the pipeline depends on. These reproduce the PDF exactly, so a mismatch
# means the engine or its model settings have drifted -- hard failure.
STRICT_CHECKS = {
    "dG_switch": 0.05,
    "dG_complex": 0.05,
    "dG_margin": 0.05,
    "toehold_open_off": 0.1,
    "aug_on": 0.1,
    "trigger_binding": 0.1,
    "trigger_unfolded": 0.1,
    "offtarget_dG": 0.05,
    "offtarget_pct_of_intended": 1.0,
    "aug_bulge_paired_off": 0.01,
}

# Reported with the delta, but not asserted. The PDF does not say which window
# it averaged the RBS over; the 11-nt AACAGAGGAGA lands within ~0.6 and the
# 8-nt AGAGGAGA within ~0.8, neither exact. This is also the row that turns out
# to be mislabelled (section 2a) and is excluded from ranking, so forcing
# agreement would mean fitting to a number we do not use.
ADVISORY_CHECKS = {
    "rbs_off": 1.0,
    "rbs_on": 1.0,
}

# The PDF prints these as "~52% of intended binding strength" etc.
PDF_OFFTARGET_PCT = {1: 52.0, 2: 52.0, 3: 49.0, 4: 49.0, 5: 43.0}


def _pdf_value(cand_id, key):
    """The PDF's own value for a metric, normalised to a number."""
    pdf = cd.CANDIDATES[cand_id]["pdf"]
    if key == "offtarget_pct_of_intended":
        return PDF_OFFTARGET_PCT[cand_id]
    if key == "aug_bulge_paired_off":
        # stored as text, e.g. "clean (0/3 paired)" or "1/3 paired (acceptable)"
        match = re.search(r"(\d)\s*/\s*3", pdf["aug_bulge_off"])
        return float(match.group(1)) if match else float("nan")
    return pdf[key]


def recompute(cand_id, reporter_rna, temperature_c=37.0):
    """Recompute every PDF row for one candidate. Thin wrapper on pdf_metrics."""
    return pdf_metrics.pdf_table(cand_id, reporter_rna, temperature_c)


def _print_check_block(rows, checks, title):
    """Print one comparison table and return the values that missed."""
    misses = []
    print()
    print(title)
    line = "%-26s" % "metric"
    for r in rows:
        line += "  cand %d      " % r["cand"]
    print(line)
    print("%-26s" % "" + "   ours    PDF ok " * len(rows))
    print("-" * 100)
    for key, tol in checks.items():
        line = "%-26s" % key
        for r in rows:
            ours = r[key]
            theirs = _pdf_value(r["cand"], key)
            ok = abs(ours - theirs) <= tol
            if not ok:
                misses.append((r["cand"], key, ours, theirs))
            line += " %6.2f %6.1f %s " % (ours, theirs, "OK" if ok else "XX")
        print(line)
    return misses


def check_against_pdf(rows):
    """Compare every recomputed row with the PDF."""
    print()
    print("=" * 100)
    print("1. REPRODUCTION -- do we get the same numbers the PDF reports?")
    print("=" * 100)

    failures = _print_check_block(
        rows, STRICT_CHECKS,
        "1a. Rows the pipeline depends on (must match exactly):")

    advisory = _print_check_block(
        rows, ADVISORY_CHECKS,
        "1b. RBS rows (reported, not asserted -- the PDF's window is unknown):")

    print()
    print("1c. Sequence-integrity rows (PDF reports yes / yes / 0 / 0 for all five):")
    print("%-26s %14s %16s %16s" %
          ("cand", "RBS present", "frame intact", "extra starts/stops"))
    for r in rows:
        print("%-26d %14s %16s %16s" % (
            r["cand"],
            "yes" if r["rbs_present"] else "NO",
            "yes" if r["reading_frame_intact"] else "NO",
            "%d / %d" % (r["unwanted_starts"], r["unwanted_stops"])))

    print()
    print("1d. Unresolved row -- 'Region just after the start codon open - ON':")
    print("    The PDF prints 3.4 / 0.6 / 1.7 / 1.8 / 4.6 and we cannot")
    print("    reproduce that from its description. Roughly twenty definitions")
    print("    were tried (mean unpaired, mean paired, joint-unpaired via a")
    print("    constrained partition function, over windows of 3-21 nt, on the")
    print("    b_pre* and linker spans, in both states). None match. Our")
    print("    closest is the mean paired probability over the 6 nt after the")
    print("    start codon, shown here ONLY so the gap is visible:")
    line = "      ours:"
    for r in rows:
        line += " %5.1f" % r["after_start_on_UNRESOLVED"]
    print(line)
    print("      PDF :" + "".join(" %5.1f" % cd.CANDIDATES[r["cand"]]["pdf"]["after_start_on"]
                                  for r in rows))
    print("    This row is NOT used anywhere. To settle it we need the")
    print("    definition from whoever wrote the original validation script.")

    print()
    print("1e. Not wired yet: the VISTA similarity score (PDF column 1).")

    print("-" * 100)
    n_strict = len(STRICT_CHECKS) * len(rows)
    if failures:
        print("FAILED: %d of %d strict values no longer reproduce the PDF:"
              % (len(failures), n_strict))
        for cand, key, ours, theirs in failures:
            print("   cand %d  %-24s ours %.2f  PDF %.2f" % (cand, key, ours, theirs))
    else:
        print("All %d strict values reproduce the PDF exactly." % n_strict)

    if advisory:
        worst = max(abs(o - t) for _, _, o, t in advisory)
        print("RBS rows differ by up to %.2f points on %d value(s) -- expected, see above."
              % (worst, len(advisory)))
    return failures


def report_rbs_correction(rows):
    """The row labelled "RBS hidden - OFF" is the UNPAIRED probability."""
    print()
    print("=" * 100)
    print("2a. CORRECTION -- 'RBS hidden - OFF' is the UNPAIRED probability")
    print("=" * 100)
    print("The PDF reads a high number here as good protection. It is the")
    print("opposite: it is how ACCESSIBLE the RBS is. Compare OFF with ON --")
    print("the RBS hardly moves, so this row carries no ranking information.")
    print()
    print("%-6s %12s %12s %10s   %12s %12s %10s" %
          ("cand", "RBS off", "RBS on", "change", "AUG off", "AUG on", "change"))
    for r in rows:
        print("%-6d %11.1f%% %11.1f%% %9.1f   %11.1f%% %11.1f%% %9.1f" % (
            r["cand"], r["rbs_off"], r["rbs_on"], r["rbs_on"] - r["rbs_off"],
            r["aug_off"], r["aug_on"], r["aug_gap_on_minus_off"]))
    print()
    print("The start codon is the row that actually discriminates: it moves")
    print("~13-28 points between states where the RBS moves ~0-7.")


def report_main_stem_correction():
    """The PDF's "fraction of the main stem that pairs correctly" is arithmetic."""
    print()
    print("=" * 100)
    print("2b. CORRECTION -- 'main stem 50/67/75%' is (n-3)/n, not folding quality")
    print("=" * 100)
    print("%-6s %8s %10s %14s %14s %10s" %
          ("cand", "k1_len", "(n-3)/n", "main stem", "upper stem", "PDF says"))

    pdf_pct = {1: "50%", 2: "50%", 3: "67%", 4: "75%", 5: "67%"}
    for cand_id in sorted(cd.CANDIDATES):
        switch = cd.CANDIDATES[cand_id]["switch"]
        k1_len = cd.CANDIDATES[cand_id]["k1_len"]
        spans = cd.domains(cand_id)
        structure, _ = RNA.fold(switch)
        pair_table = RNA.ptable(structure)

        main_formed, main_total = _count_helix(
            pair_table, spans["b_pre"], spans["b_pre_star"])
        up_formed, up_total = _count_helix(
            pair_table, spans["upper_stem"], spans["upper_stem_star"])

        print("%-6d %8d %10s %10d/%-3d %10d/%-3d %10s" % (
            cand_id, k1_len, "%.0f%%" % (100.0 * (k1_len - 3) / k1_len),
            main_formed, main_total, up_formed, up_total, pdf_pct[cand_id]))
    print()
    print("Every intended base pair forms in every candidate. The PDF's rising")
    print("percentage is the 3-nt bulge sitting in the denominator, so its")
    print("conclusion that a longer main stem folds more cleanly is not")
    print("supported by its own data.")


def _count_helix(pair_table, span_5, span_3):
    """
    How many of an intended helix's base pairs actually form in the MFE fold.

    The two arms run antiparallel: the first base of the 5' arm pairs with the
    LAST base of the 3' arm. pair_table is 1-indexed, 0 meaning unpaired.
    """
    a_start, a_end = span_5
    b_start, b_end = span_3
    total = a_end - a_start
    formed = 0
    for offset in range(total):
        i = a_start + offset
        expected_partner = b_end - 1 - offset
        if pair_table[i + 1] - 1 == expected_partner:
            formed += 1
    return formed, total


def report_length_correction(rows):
    """The PDF's dG margin ranking is largely just trigger length."""
    print()
    print("=" * 100)
    print("2c. CORRECTION -- the dG margin ranking tracks trigger length")
    print("=" * 100)
    print("%-6s %10s %12s %12s %16s %12s" %
          ("cand", "trig len", "margin", "per nt", "margin-dG(trig)", "offtarget%"))
    for r in rows:
        pct = r["offtarget_pct_of_intended"]
        print("%-6d %10d %12.1f %12.2f %16.1f %11.1f%%" % (
            r["cand"], r["trigger_len"], r["dG_margin"],
            r["dG_margin_per_nt"], r["dG_margin_corrected"],
            pct if pct is not None else float("nan")))

    by_raw = [r["cand"] for r in sorted(rows, key=lambda r: r["dG_margin"])]
    by_per_nt = [r["cand"] for r in sorted(rows, key=lambda r: r["dG_margin_per_nt"])]
    print()
    print("ranking by raw margin    : %s" % by_raw)
    print("ranking by margin per nt : %s" % by_per_nt)
    print()
    print("The 'margin-dG(trig)' column adds back the trigger's own folding")
    print("energy, which the PDF's margin row leaves out -- worth 1.3 to 6.2")
    print("kcal/mol. Note the PDF uses THAT corrected value as the denominator")
    print("of its off-target percentage, so its own two rows disagree on which")
    print("margin is the real one. Our offtarget%% column reproduces the PDF's,")
    print("which is how we know.")


def report_trigger_versions(rows, mcherry_path, versions):
    """Where the printed triggers differ from the gene that will be made."""
    print()
    print("=" * 100)
    print("3. TRIGGER VERSIONS -- printed vs real mCherry")
    print("=" * 100)
    original = versions["original"]
    out = {}
    print("%-6s %10s %8s   %s" % ("cand", "gene pos", "diffs", "difference"))
    for r in rows:
        cand_id = r["cand"]
        real = cd.trigger_real(cand_id, original)
        out[cand_id] = real
        printed = cd.CANDIDATES[cand_id]["trigger"]
        diffs = [(i, printed[i], real["sequence"][i])
                 for i in range(len(printed)) if printed[i] != real["sequence"][i]]
        desc = ", ".join("pos %d %s->%s" % (i, a, b) for i, a, b in diffs) or "identical"
        print("%-6d %10d %8d   %s" % (cand_id, real["gene_pos"], real["n_diff"], desc))
    return out


def write_csv(rows, out_dir):
    """Save the recomputed table so later stages need not redo stage 0."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "stage0_validation.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    return path


def run(config):
    """Stage 0 entry point. Returns the recomputed rows."""
    cand_ids = config["candidates"]
    temperature_c = config["temperature_C"]

    print()
    print("STAGE 0 -- validating our tooling against the PDF")
    print("engines: %s" % folding.engine_report())

    print("reading mCherry from %s" % config["mcherry_file"])
    versions = cd.read_mcherry(config["mcherry_file"], verbose=True)
    reporter_rna = cd.to_rna(versions["pdf_reference"])

    rows = [recompute(cid, reporter_rna, temperature_c) for cid in cand_ids]

    failures = check_against_pdf(rows)
    report_rbs_correction(rows)
    report_main_stem_correction()
    report_length_correction(rows)
    report_trigger_versions(rows, config["mcherry_file"], versions)

    path = write_csv(rows, config["out_dir"])
    print()
    print("wrote %s" % path)

    if failures:
        raise SystemExit(
            "Stage 0 failed: %d value(s) no longer reproduce the PDF. "
            "Fix this before trusting anything downstream." % len(failures))
    return rows


# Pressing Run on this file alone does stage 0, using main.py's CONFIG.
if __name__ == "__main__":
    from poc_and.main import CONFIG
    run(dict(CONFIG))
