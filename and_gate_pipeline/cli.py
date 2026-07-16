"""Command-line entry point.

    python -m and_gate_pipeline --demo --out results
    python -m and_gate_pipeline --gene1 @g1.fasta --gene2 @g2.fasta --out results
    python -m and_gate_pipeline --gene1 ACGU... --gene2 ACGU... --Lx 12 --LA 36 --LB 30
"""

from __future__ import annotations

import argparse
import os
import sys

from .config import PipelineConfig
from .pipeline import run_pipeline


def _read_seq(arg: str) -> str:
    """A raw sequence, or ``@path`` to a FASTA/plain file (first record)."""
    if arg.startswith("@"):
        path = arg[1:]
        with open(path) as fh:
            lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith(">")]
        return "".join(lines)
    return arg


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="and_gate_pipeline",
        description="Design two-input RNA AND-gate toehold switches.")
    p.add_argument("--demo", action="store_true",
                   help="run on the bundled example genes")
    p.add_argument("--gene1", help="Gene 1 sequence or @file")
    p.add_argument("--gene2", help="Gene 2 sequence or @file")
    p.add_argument("--reporter", default="",
                   help="downstream reporter CDS (optional)")
    p.add_argument("--out", default="and_gate_results", help="output directory")
    p.add_argument("--no-viz", action="store_true", help="skip arc plots")
    p.add_argument("--no-optimize", action="store_true",
                   help="skip restricted-sequence optimisation")
    # tunables
    p.add_argument("--Lx", type=int)
    p.add_argument("--LA", type=int)
    p.add_argument("--LB", type=int)
    p.add_argument("--len-a", type=int)
    p.add_argument("--len-k1", type=int)
    p.add_argument("--len-r2", type=int)
    p.add_argument("--secondary-loop", type=int)
    p.add_argument("--off-mfe-target", type=float)
    p.add_argument("--top-n", type=int)
    p.add_argument("--max-full-score", type=int, default=40)
    p.add_argument("--no-nupack", action="store_true",
                   help="force the ViennaRNA backend even if NUPACK is present")
    return p


def config_from_args(args) -> PipelineConfig:
    cfg = PipelineConfig()
    for attr, val in (("Lx", args.Lx), ("L_A", args.LA), ("L_B", args.LB),
                      ("len_a", args.len_a), ("len_k1", args.len_k1),
                      ("len_r2", args.len_r2),
                      ("secondary_loop_len", args.secondary_loop),
                      ("off_state_mfe_target", args.off_mfe_target),
                      ("top_n", args.top_n)):
        if val is not None:
            setattr(cfg, attr, val)
    if args.no_nupack:
        cfg.prefer_nupack = False
    return cfg


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    cfg = config_from_args(args)

    transcriptome = essential = None
    if args.demo:
        from . import examples
        gene1, gene2 = examples.GENE1, examples.GENE2
        transcriptome, essential = examples.TRANSCRIPTOME, examples.ESSENTIAL
        reporter = args.reporter
    else:
        if not (args.gene1 and args.gene2):
            print("error: provide --gene1 and --gene2 (or use --demo)",
                  file=sys.stderr)
            return 2
        gene1, gene2 = _read_seq(args.gene1), _read_seq(args.gene2)
        reporter = args.reporter

    out = run_pipeline(
        gene1, gene2, cfg, reporter=reporter,
        transcriptome=transcriptome, essential_genes=essential,
        max_full_score=args.max_full_score,
        optimize=not args.no_optimize,
        out_dir=None if args.no_viz else args.out)

    if args.no_viz and args.out:
        os.makedirs(args.out, exist_ok=True)
        out.to_csv(os.path.join(args.out, "and_gate_designs_ranked.csv"))

    print()
    print(str(out.integrity))
    print(f"\nbackend={out.backend}  candidates={out.n_candidates}  "
          f"scored={out.n_scored}")
    for r in out.results[:cfg.top_n]:
        print(f"  rank {r.rank:2d}  score={r.score.total:6.3f}  "
              f"{r.pair.orientation:14s}  hamming={r.pair.hamming}  "
              f"flags={','.join(r.score.flags) or 'none'}")
    if args.out:
        print(f"\nwrote results to: {os.path.abspath(args.out)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
