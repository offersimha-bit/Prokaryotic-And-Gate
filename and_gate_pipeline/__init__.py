"""
and_gate_pipeline
=================

Computational design pipeline for two-input RNA AND-gate toehold switches.

Implements the specification distilled from:

* Kim et al. (2019) -- inhibitory-hairpin logic and the short spacer 'a' that
  couples two sequential strand-invasion events into an AND gate.
* Green et al. (2026) / Toehold-VISTA -- the Series-A primary stem (18 bp,
  6 nt invasion), the SED/NED accessibility features, and the +/-100 nt
  flanking-region emphasis.

Layout: the file names carry the stage number
---------------------------------------------
    01_target_scan.py          stage 1  find trigger pairs (pooled or two-gene)
    02_filtering.py            stage 2  accessibility of each trigger in context
    03_select.py               stage 3  Pareto + diversity shortlist
    04_build_switch.py         stage 4  Kim hairpin + VISTA Series-A module
    04_build_switch_legacy.py  stage 4  superseded builder (kept for comparison)
    04_optimize.py             stage 4  restricted-sequence repair
    05_kinetics.py             stage 5  displacement rate vs mRNA decay
    05_scoring_legacy.py       stage 5  superseded hand-weighted scorer
    05_truth_table.py          stage 5  four-condition check for one design
    05_sweep.py                stage 5  L_x and |a| parameter sweeps
    06_offtarget.py            stage 6  transcriptome complementarity scan
    06_visualize.py            stage 6  arc plots

    config.py sequence_utils.py thermo.py constraints.py examples.py
                               shared    used by several stages, no stage number
    pipeline.py                          orchestrator across all stages
    _loader.py                           makes the numbered files importable

Because a Python identifier cannot start with a digit, the numbered files are
loaded by path and registered under stable aliases (see :mod:`._loader`).  Every
import in the package and the tests keeps working unchanged::

    from and_gate_pipeline.target_scan import scan_pool     # 01_target_scan.py
    from and_gate_pipeline.kinetics import and_behaviour    # 05_kinetics.py

All thermodynamics go through one backend abstraction (:mod:`.thermo`) that
prefers NUPACK 4 (matching the VISTA reference model) and falls back to
ViennaRNA when NUPACK is absent.
"""

from . import _loader as _loader

_loader.load_all()

from .config import PipelineConfig                       # noqa: E402
from .thermo import get_backend, ThermoBackend           # noqa: E402

#: Which file implements which stage -- ``python -m and_gate_pipeline map``.
STAGE_MAP = _loader.stage_map()

__all__ = ["PipelineConfig", "get_backend", "ThermoBackend", "STAGE_MAP"]
__version__ = "0.2.0"
