"""
and_gate_pipeline
=================

Computational design pipeline for two-input RNA AND-gate toehold switches.

The pipeline implements the specification distilled from:

* Kim et al. (2019) -- logic of inhibitory hairpins and the short spacer 'a'
  that couples two sequential strand-invasion events into an AND gate.
* Green et al. (2026) / Toehold-VISTA -- the Series-A primary-stem architecture
  (18 bp stem, 6 nt invasion), the SED/NED accessibility features, and the
  +/-100 nt flanking-region emphasis.

Stages
------
1. Target scanning and trigger definition          (:mod:`.target_scan`)
2. Thermodynamic filtering (MFE + accessibility)    (:mod:`.filtering`)
3. AND-gate toehold-switch architecture             (:mod:`.architecture`)
4. System constraints / equations                   (:mod:`.constraints`)
5-7. Multi-stage scoring and ranking                 (:mod:`.scoring`)
     off-target scan                                 (:mod:`.offtarget`)
     sequence optimisation                           (:mod:`.optimize`)
     arc-plot visualisation                          (:mod:`.visualize`)

All thermodynamics go through a single backend abstraction
(:mod:`.thermo`) that prefers NUPACK 4 (matching the VISTA reference model)
and transparently falls back to ViennaRNA when NUPACK is not installed.
"""

from .config import PipelineConfig
from .thermo import get_backend, ThermoBackend

__all__ = ["PipelineConfig", "get_backend", "ThermoBackend"]
__version__ = "0.1.0"
