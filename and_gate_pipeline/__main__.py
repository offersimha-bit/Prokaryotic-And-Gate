"""Package entry point.

Because the stage modules are named ``05_truth_table.py`` and friends, they can
no longer be reached with ``python -m and_gate_pipeline.truth_table`` (Python
will not import a module whose name starts with a digit).  Those entry points
are subcommands here instead:

    python -m and_gate_pipeline map                 # which file is which stage
    python -m and_gate_pipeline run --demo --out results
    python -m and_gate_pipeline scan  genes/        # stage 1 only, pooled FASTA
    python -m and_gate_pipeline truth-table         # stage 5, one design
    python -m and_gate_pipeline sweep               # stage 5, L_x / |a| sweeps
    python -m and_gate_pipeline audit               # spec audit

``run`` is the default, so the old invocation still works:

    python -m and_gate_pipeline --demo --out results
"""

from __future__ import annotations

import sys


def _cmd_map(argv) -> int:
    from . import STAGE_MAP
    width = max(len(f) for _s, f, _a in STAGE_MAP)
    current = None
    for stage, filename, alias in STAGE_MAP:
        if stage != current:
            print(f"\n[{stage}]")
            current = stage
        print(f"  {filename.ljust(width)}   ->  .{alias}")
    unavailable = getattr(sys.modules[__package__], "__unavailable__", {})
    if unavailable:
        print("\nnot loaded (optional dependency missing):")
        for alias, exc in unavailable.items():
            print(f"  .{alias}: {type(exc).__name__}: {exc}")
    return 0


def _cmd_scan(argv) -> int:
    """Stage 1 in isolation: pooled FASTA in, trigger pairs out."""
    import argparse
    from .config import PipelineConfig
    from .target_scan import scan_from_fasta

    p = argparse.ArgumentParser(prog="and_gate_pipeline scan")
    p.add_argument("paths", nargs="+", help="FASTA files or a folder of them")
    p.add_argument("--Lx", type=int)
    p.add_argument("--exact", action="store_true",
                   help="require a perfect connector match (no Hamming fallback)")
    p.add_argument("--top", type=int, default=10)
    args = p.parse_args(argv)

    cfg = PipelineConfig()
    if args.Lx:
        cfg.Lx = args.Lx
    pairs = scan_from_fasta(args.paths, cfg, require_exact=args.exact,
                            progress=print)
    for k, pr in enumerate(pairs[:args.top], 1):
        a, b = pr.triggerA, pr.triggerB
        print(f"\n#{k}  hamming={pr.hamming}  {pr.orientation}")
        print(f"   A @{a.window[0]}-{a.window[1]} in "
              f"{pr.meta.get('gene_a_name')}: {a.seq}")
        print(f"      k1={a.k1} a={a.a} x={a.x} r1={a.r1}")
        print(f"   B @{b.window[0]}-{b.window[1]} in "
              f"{pr.meta.get('gene_b_name')}: {b.seq}")
        print(f"      k2={b.k2} r2={b.r2}")
    return 0


def _run_module_main(alias: str) -> int:
    """Run a stage module's ``__main__``-style demo block.

    The numbered modules kept their original ``if __name__ == '__main__':``
    demos; those no longer fire, so the block was factored into ``main()`` in
    each of them.  This dispatches to it and reports clearly when a module could
    not be loaded at all (missing ViennaRNA, matplotlib, ...).
    """
    pkg = sys.modules[__package__]
    unavailable = getattr(pkg, "__unavailable__", {})
    if alias in unavailable:
        exc = unavailable[alias]
        print(f"cannot run '{alias}': {type(exc).__name__}: {exc}",
              file=sys.stderr)
        return 1
    module = getattr(pkg, alias, None)
    if module is None or not hasattr(module, "main"):
        print(f"module '{alias}' has no main()", file=sys.stderr)
        return 1
    return module.main() or 0


def _cmd_run(argv) -> int:
    from .cli import main
    return main(argv)


_COMMANDS = {
    "map": _cmd_map,
    "scan": _cmd_scan,
    "run": _cmd_run,
    "truth-table": lambda argv: _run_module_main("truth_table"),
    "sweep": lambda argv: _run_module_main("sweep"),
    "audit": lambda argv: _run_module_main("spec_audit"),
}


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in _COMMANDS:
        return _COMMANDS[argv[0]](argv[1:])
    return _cmd_run(argv)          # default: the design run


if __name__ == "__main__":
    raise SystemExit(main())
