"""Importer for the stage-numbered module files.

Why this exists
---------------
The modules are named after their pipeline stage (``01_target_scan.py``,
``02_filtering.py``, ...) so the folder listing *is* the pipeline.  Python
cannot import those names directly -- ``import 01_target_scan`` is a syntax
error, because an identifier may not start with a digit.

So each numbered file is loaded by path and registered in ``sys.modules`` under
a plain alias, in dependency order.  From that point on everything behaves
exactly as before:

    from . import sequence_utils as su          # unchanged
    from .target_scan import TriggerPair        # unchanged
    from and_gate_pipeline.kinetics import ...  # unchanged

No import statement anywhere in the package or the tests had to change: the
numbers are for humans reading the directory, the aliases are for Python.

Order matters.  A module's dependencies must already be registered when its body
executes, which is why ``_MODULES`` is a list and not a dict comprehension over
the directory.
"""

from __future__ import annotations

import importlib.util
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG = __name__.rsplit(".", 1)[0]

#: (filename on disk, import alias) -- in dependency order.
#: Unnumbered files are cross-cutting: they are used by several stages and do
#: not belong to any single one.
_MODULES: list[tuple[str, str]] = [
    # -- cross-cutting foundations -------------------------------------- #
    ("sequence_utils.py",        "sequence_utils"),
    ("config.py",                "config"),
    ("thermo.py",                "thermo"),
    ("constraints.py",           "constraints"),
    ("examples.py",              "examples"),
    # -- stage 1: find triggers ----------------------------------------- #
    ("01_target_scan.py",        "target_scan"),
    # -- stage 2: accessibility in gene context ------------------------- #
    ("02_filtering.py",          "filtering"),
    # -- stage 3: candidate selection ----------------------------------- #
    ("03_select.py",             "select"),
    # -- stage 4: build the switch -------------------------------------- #
    # legacy first: 04_build_switch.py reuses its _weaken_arm()
    ("04_build_switch_legacy.py", "architecture"),
    ("04_build_switch.py",       "vista_switch"),
    ("04_optimize.py",           "optimize"),
    # -- stage 6 (loaded early): 05_scoring_legacy.py imports offtarget -- #
    ("06_offtarget.py",          "offtarget"),
    # -- stage 5: score ------------------------------------------------- #
    ("05_kinetics.py",           "kinetics"),
    ("05_scoring_legacy.py",     "scoring"),
    # -- stage 6: output ------------------------------------------------ #
    ("06_visualize.py",          "visualize"),
    # -- orchestration and analysis tools ------------------------------- #
    ("pipeline.py",              "pipeline"),
    ("05_truth_table.py",        "truth_table"),
    ("05_sweep.py",              "sweep"),
    ("spec_audit.py",            "spec_audit"),
    ("cli.py",                   "cli"),
]

#: alias -> filename, for error messages and for the CLI's ``--map`` listing.
FILE_OF = {alias: fname for fname, alias in _MODULES}

#: alias -> stage number (None for cross-cutting modules).
STAGE_OF = {alias: (int(fname[:2]) if fname[:2].isdigit() else None)
            for fname, alias in _MODULES}


def _load_one(filename: str, alias: str):
    full_name = f"{_PKG}.{alias}"
    if full_name in sys.modules:
        return sys.modules[full_name]

    path = os.path.join(_HERE, filename)
    if not os.path.exists(path):
        raise ImportError(
            f"{_PKG}: expected {filename!r} (alias {alias!r}) in {_HERE}. "
            "If you renamed or moved a module, update _loader._MODULES so the "
            "stage numbering and the import aliases stay in sync.")

    spec = importlib.util.spec_from_file_location(full_name, path)
    module = importlib.util.module_from_spec(spec)
    # Register BEFORE exec so circular-ish `from . import x` resolves, and set
    # it as an attribute of the package so `from . import x` finds it via
    # getattr without the finder ever looking for an `x.py` on disk.
    sys.modules[full_name] = module
    pkg = sys.modules.get(_PKG)
    if pkg is not None:
        setattr(pkg, alias, module)
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(full_name, None)
        raise
    return module


def load_all(strict: bool = False) -> dict:
    """Load every numbered module and register its alias.

    ``strict=False`` (the default) skips modules whose optional dependencies are
    missing -- ``06_visualize.py`` needs matplotlib and networkx, and
    ``05_truth_table.py`` imports ViennaRNA directly -- so a stock install can
    still run the core pipeline.  Failures are recorded in the returned dict.
    """
    loaded, failed = {}, {}
    for filename, alias in _MODULES:
        try:
            loaded[alias] = _load_one(filename, alias)
        except Exception as exc:
            if strict:
                raise
            failed[alias] = exc
    if failed:
        pkg = sys.modules.get(_PKG)
        if pkg is not None:
            pkg.__unavailable__ = failed
    return loaded


def stage_map() -> list[tuple[str, str, str]]:
    """[(stage label, filename, alias), ...] sorted by pipeline stage.

    ``_MODULES`` is in *dependency* order (which is not stage order -- e.g.
    ``06_offtarget.py`` must load before ``05_scoring_legacy.py``), so this
    re-sorts for human consumption.
    """
    rows = [((STAGE_OF[alias] or 99), filename, alias)
            for filename, alias in _MODULES]
    rows.sort(key=lambda r: (r[0], r[1]))
    return [(f"stage {n}" if n != 99 else "shared / orchestration", f, a)
            for n, f, a in rows]
