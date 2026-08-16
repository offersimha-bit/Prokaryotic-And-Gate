#!/usr/bin/env python3
"""Run the whole AND-gate pipeline on the real genes, end to end.

    press Run in VS Code, or:      python main.py
    with options:                  python main.py --help

WHY THIS FILE EXISTS
--------------------
The stage modules are named ``01_target_scan.py`` ... ``07_rank.py`` so the
directory listing reads as the pipeline.  A Python module name cannot start with
a digit, so ``_loader.py`` loads them by path and registers stable aliases, and
they use package-relative imports.  The consequence: opening any stage file and
pressing Run gives

    ImportError: attempted relative import with no known parent package

That is expected.  This file sits OUTSIDE the package, so it has no relative
imports and the Run button works on it.  It is the one file to run.

Equivalent from a terminal:  python -m and_gate_pipeline run --demo
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

# Run button friendliness: VS Code launches with cwd set to whatever folder is
# open, and the package only imports from the repo root.  Pin it here rather
# than relying on the caller getting it right.
ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from and_gate_pipeline import examples                       # noqa: E402
from and_gate_pipeline.config import PipelineConfig          # noqa: E402
from and_gate_pipeline.pipeline import run_pipeline          # noqa: E402
from and_gate_pipeline.thermo import get_backend             # noqa: E402


DEFAULT_OUT = ROOT / "and_gate_results"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="AND-gate design pipeline, real genes end to end.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--gene1", default=examples.GENE1_NAME,
                   help="bundled gene name, or a path to a FASTA file")
    p.add_argument("--gene2", default=examples.GENE2_NAME,
                   help="bundled gene name, or a path to a FASTA file")
    p.add_argument("--out", default=str(DEFAULT_OUT), help="output directory")
    p.add_argument("--stage2-max", type=int, default=None,
                   help="cap on pairs measured in stage 2 (cost control). "
                        "Lower it for a fast smoke run, e.g. 40")
    p.add_argument("--top", type=int, default=40,
                   help="how many shortlisted designs to score fully")
    p.add_argument("--Lx", type=int, default=None,
                   help="connector length; the sweep found 7 optimal while the "
                        "config default is still 6")
    p.add_argument("--no-viz", action="store_true",
                   help="skip arc plots (they need networkx)")
    p.add_argument("--quick", action="store_true",
                   help="smoke run: stage2 cap 40, top 8, no plots")
    return p.parse_args(argv)


def _load(name_or_path: str) -> tuple[str, str]:
    """Accept either a bundled gene name or a path to a FASTA file."""
    path = Path(name_or_path)
    if path.exists():
        return path.stem, examples.read_fasta(path)
    return name_or_path, examples.load_gene(name_or_path)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.quick:
        args.stage2_max, args.top, args.no_viz = 40, 8, True

    name1, gene1 = _load(args.gene1)
    name2, gene2 = _load(args.gene2)

    overrides = {}
    if args.stage2_max is not None:
        overrides["stage2_max_candidates"] = args.stage2_max
    if args.Lx is not None:
        overrides["Lx"] = args.Lx
    cfg = PipelineConfig(**overrides)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print("AND-gate pipeline")
    print("=" * 72)
    print(f"  input      {name1}  {len(gene1):>5} nt")
    print(f"             {name2}  {len(gene2):>5} nt")
    print(f"  engine     {get_backend(cfg).name}")
    print(f"  geometry   Lx={cfg.Lx}  |a|={cfg.len_a}  L_A={cfg.L_A}  L_B={cfg.L_B}")
    print(f"  stage-2 cap {cfg.stage2_max_candidates}   full-score top {args.top}")
    print(f"  output     {out_dir}")
    print("-" * 72)

    started = time.time()
    out = run_pipeline(
        gene1, gene2, cfg,
        transcriptome=examples.TRANSCRIPTOME,
        essential_genes=examples.ESSENTIAL,
        max_full_score=args.top,
        out_dir=str(out_dir),
        # {} skips the whole-gene arc plots (the expensive ones); passing None
        # would make run_pipeline substitute the genes back in
        viz_genes={} if args.no_viz else {name1: gene1, name2: gene2},
        progress=print,
    )
    elapsed = time.time() - started

    print("-" * 72)
    if not out.results:
        print("no designs produced -- see the stage messages above")
        return 1

    print(f"{len(out.results)} designs scored in {elapsed:.0f}s "
          f"on {out.backend}; {out.n_candidates} pairs scanned")
    print()
    head = out.results[0].row()
    cols = [c for c in ("rank", "logic_margin", "P_11", "P_10", "P_00",
                        "on_off", "total") if c in head]
    print("  " + "".join(f"{c:>16}" for c in cols))
    for r in out.results[:10]:
        row = r.row()
        cells = []
        for c in cols:
            v = row[c]
            cells.append(f"{v:>16.4g}" if isinstance(v, float) else f"{v:>16}")
        print("  " + "".join(cells))

    print()
    for f in sorted(out_dir.iterdir()):
        print(f"  wrote  {f.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
