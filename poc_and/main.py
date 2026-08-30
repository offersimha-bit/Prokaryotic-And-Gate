"""
poc_and -- two-input AND toehold switch pipeline.

HOW TO RUN
----------
Open this file in Visual Studio and press Run. That is all. Everything the
pipeline needs is in the CONFIG block below -- edit it, run again.

You can also run it from a terminal:

    python -m poc_and.main
    python -m poc_and.main --stages 0        (just the validation stage)
    python -m poc_and.main --explain 5       (walk through one candidate's design)

Each stage writes its results into out_dir and reads the previous stage's file
if it is already there, so re-running one stage does not redo the others.
"""

import os
import sys

# Make the package importable when this file is run directly (the Run button in
# Visual Studio executes it as a plain script, not as part of a package).
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)                 # ...\Prokaryotic-And-Gate
_PROJECT = os.path.dirname(_REPO)              # ...\Prokaryotic And Gate
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)


# ===========================================================================
# CONFIG -- edit this block, then Run
# ===========================================================================

CONFIG = {
    # ---- inputs ----------------------------------------------------------
    "mcherry_file": os.path.join(_PROJECT, "mCherry.txt"),

    # The reporter CDS that follows the linker. Leave "" to compute everything
    # switch-only, exactly as the PDF did. Paste the real GFP coding sequence
    # here when we have it and every stage picks it up automatically.
    "reporter_cds": "",

    "candidates": [1, 2, 3, 4, 5],

    # "printed" = the triggers as written in the PDF
    # "real"    = the same regions taken from mCherry.txt (4 of 5 differ by 1 nt)
    # "both"    = compute every metric twice and report side by side
    "trigger_version": "both",

    # ---- design sweep ----------------------------------------------------
    "a_range": range(3, 8),          # |a|, the toehold left exposed in OFF
    "Lx_range": range(6, 11),        # length of the k2*:x* stem
    "m_max": 7,                      # A:B complementarity -- the hard bound
    "Lx_minus_m_max": 3,             # unmatched bp that must fray for B to finish
    "bulge_sizes": [0, 1, 2, 3],     # bulge in r1copy, against RNase III
    "secondary_loop_len": 11,        # loop of the inhibitory hairpin
    "len_r2": 23,                    # trigger B's toehold, so |k2|+|r2| ~ 30 nt
    "allow_wobble_in_k2": True,      # G-U still pairs, so it widens the search
    "scan_break_offsets": True,      # let the broken block sit mid-stem, not only at the end
    "stage2_thermo_top_n": 40,       # whole-gene controls are the expensive part
    "stage3_designs_per_candidate": 5,  # let the KINETICS choose, not just stage 2

    # ---- physical --------------------------------------------------------
    "temperature_C": 37.0,
    "trigger_context_nt": 25,        # native gene flanking each trigger window
    "on_threshold_pct": 80.0,        # only for the confusion matrix, not the ranking
    "trigger_conc_M": 10e-9,         # ~10 nM, for sequestration and kinetics
    "mrna_half_life_s": 300.0,       # E. coli mRNA, 2-8 min typical
    "engine": "auto",                # "vienna" | "both"  (auto adds NUPACK if present)

    # ---- what to run, and where output goes ------------------------------
    # 5 runs last on purpose: the report needs stages 4 and 6 to have run
    "stages": [0, 2, 3, 4, 6, 5],
    "out_dir": os.path.join(_HERE, "output"),
}

# ===========================================================================


def _stage_0(config, results):
    from poc_and import validate_pdf
    return validate_pdf.run(config)


def _stage_2(config, results):
    from poc_and import find_second_trigger
    return find_second_trigger.run(config)


def _stage_3(config, results):
    from poc_and import build
    return build.run(config, results)


def _stage_4(config, results):
    from poc_and import metrics
    return metrics.run(config, results)


def _stage_6(config, results):
    from poc_and import appendix_kinetics
    return appendix_kinetics.run(config, results)


def _stage_5(config, results):
    from poc_and import report
    return report.run(config, results)


def _not_built(number, name):
    def runner(config, results):
        print()
        print("STAGE %d (%s) -- not built yet." % (number, name))
        print("We are adding stages one at a time; this one comes next.")
        return None
    return runner


STAGES = {
    0: ("validate the PDF", _stage_0),
    1: ("candidate data", _not_built(1, "candidate data")),
    2: ("find the second trigger", _stage_2),
    3: ("build the AND switches", _stage_3),
    4: ("metrics", _stage_4),
    5: ("Word and Excel report", _stage_5),
    6: ("kinetics -- decides 10 vs 11", _stage_6),
}


def parse_args(argv, config):
    """Tiny optional command-line override. Everything defaults to CONFIG."""
    args = list(argv)
    while args:
        flag = args.pop(0)
        if flag == "--stages":
            config["stages"] = [int(x) for x in args.pop(0).split(",")]
        elif flag == "--candidates":
            config["candidates"] = [int(x) for x in args.pop(0).split(",")]
        elif flag == "--out":
            config["out_dir"] = args.pop(0)
        elif flag == "--explain":
            config["explain"] = int(args.pop(0))
        elif flag == "--pick":
            # |a|,Lx,m -- trace a specific design rather than stage 2's winner
            config["explain_pick"] = [int(x) for x in args.pop(0).split(",")]
        elif flag in ("-h", "--help"):
            print(__doc__)
            raise SystemExit(0)
        else:
            raise SystemExit("unknown option: %s" % flag)
    return config


def main(argv=None):
    config = dict(CONFIG)
    config = parse_args(argv if argv is not None else sys.argv[1:], config)

    os.makedirs(config["out_dir"], exist_ok=True)

    # --explain prints the construction for one candidate and stops
    if config.get("explain"):
        from poc_and import explain as explain_mod
        explain_mod.explain(config["explain"], config,
                            config.get("explain_pick"))
        return None

    print("=" * 78)
    print("poc_and -- two-input AND toehold switch pipeline")
    print("=" * 78)
    print("candidates : %s" % config["candidates"])
    print("stages     : %s" % config["stages"])
    print("output     : %s" % config["out_dir"])
    print("reporter   : %s" % (("%d nt" % len(config["reporter_cds"]))
                               if config["reporter_cds"] else "none (switch-only, as in the PDF)"))

    results = {}
    for number in config["stages"]:
        if number not in STAGES:
            raise SystemExit("no such stage: %s" % number)
        name, runner = STAGES[number]
        results[number] = runner(config, results)

    print()
    print("=" * 78)
    print("done.")
    return results


if __name__ == "__main__":
    main()
