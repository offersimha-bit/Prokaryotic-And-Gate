"""
poc_and -- two-input AND toehold switch pipeline.

HOW TO RUN
----------
Open this file in Visual Studio and press Run. That is all. Everything the
pipeline needs is in the CONFIG block below -- edit it, run again.

You can also run it from a terminal:

    python -m poc_and.main
    python -m poc_and.main --stages 0        (just the validation stage)

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

    # ---- physical --------------------------------------------------------
    "temperature_C": 37.0,
    "trigger_conc_M": 10e-9,         # ~10 nM, for the sequestration calculation
    "engine": "auto",                # "vienna" | "both"  (auto adds NUPACK if present)

    # ---- what to run, and where output goes ------------------------------
    "stages": [0],                   # add 1,2,3,4,5 as we build them; 6 = kinetics
    "out_dir": os.path.join(_HERE, "output"),
}

# ===========================================================================


def _stage_0(config):
    from poc_and import validate_pdf
    return validate_pdf.run(config)


def _not_built(number, name):
    def runner(config):
        print()
        print("STAGE %d (%s) -- not built yet." % (number, name))
        print("We are adding stages one at a time; this one comes next.")
        return None
    return runner


STAGES = {
    0: ("validate the PDF", _stage_0),
    1: ("candidate data", _not_built(1, "candidate data")),
    2: ("find the second trigger", _not_built(2, "find the second trigger")),
    3: ("build the AND switches", _not_built(3, "build the AND switches")),
    4: ("metrics", _not_built(4, "metrics")),
    5: ("Word and Excel report", _not_built(5, "Word and Excel report")),
    6: ("kinetics appendix", _not_built(6, "kinetics appendix")),
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
        results[number] = runner(config)

    print()
    print("=" * 78)
    print("done.")
    return results


if __name__ == "__main__":
    main()
