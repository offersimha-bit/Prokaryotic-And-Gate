"""Bridge to the standalone pooled trigger scanner
(``Triger finding/and_gate_trigger.py``).

That script solves a problem this pipeline does not: it pools **many** gene
records from any number of FASTA files and hunts for the lucky coincidence
where a Trigger-1 window's connector ``*`` is the exact reverse complement of
some Trigger-2 window's ``K2`` in a *different* gene.  This pipeline, by
contrast, takes exactly two genes but then builds and scores the real
two-hairpin switch.  Chaining them gives: pooled discovery -> full design.

Its trigger grammar is identical to ours, so the mapping is 1:1::

    Trigger1Window(r1, star, a, k1)  ->  TriggerA(r1, x=star, a, k1)
    Trigger2Window(r2, k2)           ->  TriggerB(r2, k2)
                                     ->  TriggerPair(hamming=0, exact=True)

The scanner file is imported **unmodified** and is never asked to choose its own
folding engine -- see :class:`_ScannerBackendShim`.

Usage
-----
    from and_gate_pipeline.interop import load_genes, run_from_scanner
    genes = load_genes(["genes/amr.fasta"])
    out = run_from_scanner(genes, out_dir="results")
    print(out.results[0].switch.core)
"""

from __future__ import annotations

import importlib.util
import os
import sys

from . import sequence_utils as su
from .config import PipelineConfig
from .pipeline import run_pipeline, PipelineOutput
from .target_scan import TriggerA, TriggerB, TriggerPair
from .thermo import ThermoBackend, get_backend

_DEFAULT_SCRIPT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "Triger finding", "and_gate_trigger.py")


# --------------------------------------------------------------------------- #
# loading the scanner without touching it                                     #
# --------------------------------------------------------------------------- #
_MODULE_NAME = "_and_gate_trigger_scanner"
_scanner_cache: dict = {}


def load_scanner(path: str | None = None):
    """Import the standalone scanner as a module, unmodified.

    Importing under a non-``__main__`` name means its ``main()`` never runs, so
    importing has no side effects.
    """
    path = os.path.abspath(path or _DEFAULT_SCRIPT)
    if path in _scanner_cache:
        return _scanner_cache[path]
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"trigger scanner not found at {path!r}; pass an explicit path")
    spec = importlib.util.spec_from_file_location(_MODULE_NAME, path)
    mod = importlib.util.module_from_spec(spec)
    # Must be in sys.modules *before* exec_module: the scanner defines
    # @dataclass classes, and dataclasses resolves type hints via
    # sys.modules[cls.__module__], which would otherwise be None.
    sys.modules[_MODULE_NAME] = mod
    try:
        spec.loader.exec_module(mod)                 # type: ignore[union-attr]
    except Exception:
        sys.modules.pop(_MODULE_NAME, None)
        raise
    _scanner_cache[path] = mod
    return mod


# --------------------------------------------------------------------------- #
# backend injection -- the reason the scanner is safe to use from here        #
# --------------------------------------------------------------------------- #
class _ScannerBackendShim:
    """Presents this pipeline's verified :class:`ThermoBackend` through the
    scanner's ``Backend`` duck-type (``fold`` / ``unpaired_probs`` /
    ``duplex_energy`` / ``mfe`` / ``structure_probability``).

    Why this exists
    ---------------
    The scanner's own ``select_backend()`` tries NUPACK first, and its NUPACK
    ``unpaired_probs()`` assumes the pair matrix is (n+1)x(n+1) with an unpaired
    column at index n.  NUPACK 4.1 actually returns an n x n matrix with
    P(unpaired) on the **diagonal**, so that lookup raises IndexError, which the
    scanner swallows (``except Exception: return [0.5] * n``) and silently
    reports every base as 50% unpaired.  That would corrupt accessibility *and*
    the accessibility-ordered top-k pre-selection.

    We never call the scanner's ``select_backend()``.  Its ``scan_trigger1`` /
    ``scan_trigger2`` take the backend as an argument and only ever call
    ``unpaired_probs``, so injecting this shim routes that call to the
    pipeline's engine (whose NUPACK pair-matrix handling is fixed and verified
    against ViennaRNA).  The scanner file itself needs no edit.
    """

    name = "pipeline-backend"

    def __init__(self, backend: ThermoBackend, temperature_c: float = 37.0):
        self._b = backend
        self.temperature_c = temperature_c

    def fold(self, seq: str):
        struct, energy = self._b.mfe(seq)
        return struct, float(energy), float(self._b.ensemble_energy(seq))

    def unpaired_probs(self, seq: str):
        return self._b.unpaired_probabilities(seq)

    def duplex_energy(self, a: str, b: str) -> float:
        return float(self._b.complex_mfe([a, b]))

    def mfe(self, seq: str):
        return self._b.mfe(seq)

    def structure_probability(self, seq: str, temperature_c: float) -> float:
        return self._b.structure_probability(seq, temperature_c)


# --------------------------------------------------------------------------- #
# parameter <-> config translation                                            #
# --------------------------------------------------------------------------- #
def params_from_config(mod, cfg: PipelineConfig):
    """Build the scanner's ``Params`` from a :class:`PipelineConfig`."""
    p = mod.Params()
    p.R1_len = cfg.resolved_len_r1()
    p.star_len = cfg.Lx
    p.A_len = cfg.len_a
    p.K1_len = cfg.len_k1
    p.R2_len = cfg.resolved_len_r2()
    p.temperature_c = cfg.temperature_c
    p.material = cfg.material
    return p


def config_from_params(p, **overrides) -> PipelineConfig:
    """Build a :class:`PipelineConfig` matching the scanner's ``Params``.

    Keeps the Section-4 equations satisfied:
    ``|r1|+|x|+|a|+|k1| == L_A`` and ``|k2|+|r2| == L_B``.
    """
    cfg = PipelineConfig(
        Lx=p.star_len,
        len_a=p.A_len,
        len_k1=p.K1_len,
        L_A=p.R1_len + p.star_len + p.A_len + p.K1_len,
        L_B=p.R2_len + p.star_len,
        len_r2=p.R2_len,
        temperature_c=p.temperature_c,
    )
    for k, v in overrides.items():
        setattr(cfg, k, v)
    return cfg


# --------------------------------------------------------------------------- #
# window -> TriggerPair conversion                                            #
# --------------------------------------------------------------------------- #
def windows_to_pair(w1, w2, seq_by_name: dict) -> TriggerPair:
    """Convert one scanner (Trigger1Window, Trigger2Window) match into a
    :class:`TriggerPair`.

    The scanner stores the gene *name*; this pipeline needs the gene *sequence*
    (it re-folds the trigger in its native context), so the name is resolved
    through ``seq_by_name`` and kept in ``pair.meta``.
    """
    gene_a = su.to_rna(seq_by_name[w1.gene])
    gene_b = su.to_rna(seq_by_name[w2.gene])
    ta = TriggerA(gene=gene_a, pos_x=w1.start + len(w1.r1),
                  r1=w1.r1, x=w1.star, a=w1.a, k1=w1.k1)
    tb = TriggerB(gene=gene_b, pos_k2=w2.start + len(w2.r2),
                  r2=w2.r2, k2=w2.k2)
    h = su.hamming(su.reverse_complement(w1.star), w2.k2)
    return TriggerPair(
        orientation="scanner:G(A)->A,G(B)->B",
        gene_a=gene_a, gene_b=gene_b,
        triggerA=ta, triggerB=tb,
        hamming=h, exact=(h == 0),
        meta={"gene_a_name": w1.gene, "gene_b_name": w2.gene,
              "scanner_start_A": w1.start, "scanner_start_B": w2.start,
              "scanner_access_A": w1.mean_unpaired,
              "scanner_access_B": w2.mean_unpaired},
    )


def scan_with_scanner(genes: list[tuple[str, str]],
                      cfg: PipelineConfig | None = None,
                      scanner_path: str | None = None,
                      max_pairs: int = 40,
                      progress=print) -> list[TriggerPair]:
    """Run the pooled scanner over ``genes`` [(name, seq), ...] using this
    pipeline's folding engine, and return :class:`TriggerPair` candidates
    ordered by mean trigger accessibility (best first)."""
    cfg = cfg or PipelineConfig()
    mod = load_scanner(scanner_path)
    backend = get_backend(cfg)
    shim = _ScannerBackendShim(backend, cfg.temperature_c)
    p = params_from_config(mod, cfg)

    win1 = mod.scan_trigger1(genes, shim, p)
    win2 = mod.scan_trigger2(genes, shim, p)
    matches = mod.find_matching_pairs(win1, win2)
    progress(f"[scanner] {len(win1)} Trigger-A windows, {len(win2)} Trigger-B "
             f"windows, {len(matches)} exact x==revcomp(k2) matches "
             f"(engine={backend.name})")

    # most-accessible first, then hand the best over to the full designer
    matches.sort(key=lambda wp: 0.5 * (wp[0].mean_unpaired + wp[1].mean_unpaired),
                 reverse=True)
    seq_by_name = dict(genes)
    return [windows_to_pair(w1, w2, seq_by_name)
            for w1, w2 in matches[:max_pairs]]


def load_genes(paths) -> list[tuple[str, str]]:
    """Read FASTA records from one or more files, pooled (reuses the scanner's
    own parser so behaviour matches it exactly)."""
    mod = load_scanner()
    genes: list[tuple[str, str]] = []
    for path in ([paths] if isinstance(paths, str) else paths):
        genes.extend(mod.read_fasta_records(path))
    return genes


# --------------------------------------------------------------------------- #
# end-to-end: pooled scan -> full switch design                               #
# --------------------------------------------------------------------------- #
def run_from_scanner(genes: list[tuple[str, str]],
                     cfg: PipelineConfig | None = None,
                     scanner_path: str | None = None,
                     max_pairs: int = 40,
                     progress=print, **pipeline_kwargs) -> PipelineOutput:
    """Pooled multi-gene discovery (scanner) -> full AND-gate design (pipeline).

    Any extra keyword goes straight to :func:`run_pipeline` (``reporter``,
    ``out_dir``, ``transcriptome``, ``max_full_score``, ...).
    """
    cfg = cfg or PipelineConfig()
    pairs = scan_with_scanner(genes, cfg, scanner_path, max_pairs, progress)
    if not pairs:
        progress("[scanner] no exact connector matches found; "
                 "try a shorter Lx or more genes")
    viz_genes = {name: seq for name, seq in genes}
    return run_pipeline(cfg=cfg, pairs=pairs, viz_genes=viz_genes,
                        progress=progress, **pipeline_kwargs)
