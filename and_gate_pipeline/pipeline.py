"""End-to-end orchestrator: genes in, ranked AND-gate switch designs out.

    run_pipeline(gene1, gene2, cfg) ->
        1. logical-integrity check                (constraints)
        2. target scan + trigger definition       (target_scan)  [both orientations]
        3. trigger thermodynamic filtering         (filtering)    [Stage 2]
        4. cheap pre-ranking + top-K selection
        5. build switch + optimise restricted seqs (architecture, optimize)
        6. full hierarchical scoring               (scoring)      [Stage 7]
        7. rank, write CSV, optional arc plots
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field

from . import sequence_utils as su
from .architecture import build_switch, DesignedSwitch
from .config import PipelineConfig
from .constraints import validate_config, IntegrityReport
from .filtering import evaluate_pair_triggers, TriggerMetrics
from .optimize import optimize_switch
from .scoring import DesignScorer, ScoreCard
from .target_scan import scan_both_orientations, TriggerPair
from .thermo import get_backend

def _find_codon_table() -> str | None:
    """Locate the E. coli codon-usage table.

    Prefers the copy vendored inside the package (so the repo is self-contained
    and portable); falls back to the sibling VISTA checkout if present.  Returns
    None rather than a bad path, so a missing table is visible instead of
    silently degrading the translation-efficiency score.
    """
    here = os.path.dirname(__file__)
    candidates = [
        os.path.join(here, "data", "ecoli_codon_usage_table.csv"),
        os.path.join(os.path.dirname(here), "vista", "toehold-VISTA",
                     "ecoli_codon_usage_table.csv"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


_DEFAULT_CODON = _find_codon_table()


@dataclass
class DesignResult:
    pair: TriggerPair
    switch: DesignedSwitch
    tmA: TriggerMetrics
    tmB: TriggerMetrics
    score: ScoreCard
    rank: int = 0

    def row(self) -> dict:
        p = self.pair
        base = {
            "rank": self.rank,
            "orientation": p.orientation,
            "gene_A": p.gene_a, "gene_B": p.gene_b,
            "x": p.triggerA.x, "k2": p.triggerB.k2,
            "hamming": p.hamming, "exact": p.exact,
            "triggerA": p.triggerA.seq, "triggerB": p.triggerB.seq,
            "switch_core": self.switch.core,
            "switch_full": self.switch.full,
            "off_structure": self.switch.off_structure,
        }
        base.update(self.score.as_row())
        return base


@dataclass
class PipelineOutput:
    integrity: IntegrityReport
    results: list[DesignResult] = field(default_factory=list)
    n_candidates: int = 0
    n_scored: int = 0
    backend: str = ""

    def to_csv(self, path: str) -> str:
        if not self.results:
            open(path, "w").close()
            return path
        rows = [r.row() for r in self.results]
        fields = list(rows[0].keys())
        for r in rows:                       # union of detail keys
            for k in r:
                if k not in fields:
                    fields.append(k)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        return path


def _prescore(tmA: TriggerMetrics, tmB: TriggerMetrics, pair: TriggerPair) -> float:
    """Cheap ranking to pick which candidates get the expensive full score:
    reward accessible triggers and low Hamming (good complementary core)."""
    return (tmA.accessibility + tmB.accessibility
            - 0.05 * pair.hamming
            + (0.2 if tmA.passes else 0) + (0.2 if tmB.passes else 0))


def run_pipeline(gene1: str | None = None, gene2: str | None = None,
                 cfg: PipelineConfig | None = None,
                 reporter: str = "",
                 transcriptome: dict | None = None,
                 essential_genes: set | None = None,
                 expression: dict | None = None,
                 codon_table_path: str | None = None,
                 max_candidates_per_orientation: int | None = None,
                 max_full_score: int = 40,
                 optimize: bool = True,
                 out_dir: str | None = None,
                 pairs: list | None = None,
                 viz_genes: dict | None = None,
                 progress=print) -> PipelineOutput:
    """Design AND-gate switches from two genes.

    ``pairs``: pre-computed :class:`TriggerPair` list.  When supplied, the
    built-in two-gene scan is skipped and these pairs are built/scored instead
    -- this is how :mod:`.interop` feeds in candidates found by the standalone
    pooled multi-gene scanner.  ``viz_genes``: {name: sequence} to arc-plot
    (defaults to gene1/gene2).
    """
    cfg = cfg or PipelineConfig()
    backend = get_backend(cfg)

    # 1. logical integrity ------------------------------------------------ #
    integrity = validate_config(cfg)
    out = PipelineOutput(integrity=integrity, backend=backend.name)
    if not integrity.ok:
        progress("Logical-integrity check FAILED; aborting.\n" + str(integrity))
        return out

    # 2. target scan ------------------------------------------------------ #
    if pairs is None:
        if gene1 is None or gene2 is None:
            raise ValueError("supply gene1 and gene2, or a pre-computed `pairs` list")
        pairs = scan_both_orientations(gene1, gene2, cfg,
                                       max_candidates_per_orientation)
    out.n_candidates = len(pairs)
    progress(f"[scan] {len(pairs)} candidate trigger pairs "
             f"(exact matches: {sum(p.exact for p in pairs)})")
    if not pairs:
        return out

    # 3. trigger filtering + 4. pre-rank --------------------------------- #
    scored_pre = []
    for p in pairs:
        try:
            tmA, tmB = evaluate_pair_triggers(p, cfg, backend)
        except Exception as ex:                      # pragma: no cover
            progress(f"  filtering skipped a pair: {ex}")
            continue
        scored_pre.append((_prescore(tmA, tmB, p), p, tmA, tmB))
    scored_pre.sort(key=lambda t: -t[0])
    shortlist = scored_pre[:max_full_score]
    progress(f"[filter] scoring top {len(shortlist)} of {len(scored_pre)} "
             f"candidates in full")

    # 5-6. build, optimise, full score ----------------------------------- #
    codon_path = codon_table_path or _DEFAULT_CODON
    if codon_path is None:
        progress("  [warn] no codon-usage table found; translation-efficiency "
                 "scoring will use a neutral default")
    scorer = DesignScorer(
        cfg, backend, codon_table_path=codon_path,
        transcriptome=transcriptome, essential_genes=essential_genes,
        expression=expression)

    results: list[DesignResult] = []
    for k, (_pre, p, tmA, tmB) in enumerate(shortlist, 1):
        try:
            sw = build_switch(p, cfg, reporter=reporter)
            if optimize:
                sw, _rep = optimize_switch(sw, cfg, backend)
            sc = scorer.score(sw, tmA, tmB)
        except Exception as ex:                      # pragma: no cover
            progress(f"  [{k}] scoring failed: {ex}")
            continue
        results.append(DesignResult(pair=p, switch=sw, tmA=tmA, tmB=tmB, score=sc))

    results.sort(key=lambda r: -r.score.total)
    for i, r in enumerate(results, 1):
        r.rank = i
    out.results = results
    out.n_scored = len(results)
    progress(f"[score] {len(results)} designs scored; "
             f"top total = {results[0].score.total:.3f}" if results else "no designs")

    # 7. outputs ---------------------------------------------------------- #
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        out.to_csv(os.path.join(out_dir, "and_gate_designs_ranked.csv"))
        _write_final_designs(out, cfg, os.path.join(out_dir, "final_designs.txt"))
        if viz_genes is None:
            viz_genes = {}
            if gene1 is not None:
                viz_genes["gene1"] = gene1
            if gene2 is not None:
                viz_genes["gene2"] = gene2
        _emit_visuals(out, viz_genes, cfg, backend, out_dir, progress)
    return out


def _write_final_designs(out: PipelineOutput, cfg: PipelineConfig, path: str):
    top = out.results[:cfg.top_n]
    lines = ["# AND-gate toehold-switch designs (ranked)",
             f"# backend: {out.backend}   candidates: {out.n_candidates}   "
             f"scored: {out.n_scored}", ""]
    for r in top:
        p = r.pair
        lines += [
            f">>> rank {r.rank}   score={r.score.total:.3f}   "
            f"orientation={p.orientation}   hamming={p.hamming}",
            f"    TriggerA ({p.gene_a}): {p.triggerA.seq}",
            f"    TriggerB ({p.gene_b}): {p.triggerB.seq}",
            f"    switch  : {r.switch.core}",
            f"    OFF-struct: {r.switch.off_structure}",
            f"    B={r.score.triggerB_activation:.3f} "
            f"INT={r.score.intermediate_state:.3f} "
            f"ON={r.score.triggerA_on_state:.3f} PEN={r.score.penalty:.3f}",
            f"    flags: {', '.join(r.score.flags) or 'none'}",
            "",
        ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))


def _emit_visuals(out, viz_genes, cfg, backend, out_dir, progress):
    try:
        from .visualize import arc_plot, export_pair_fraction_csv
    except Exception as ex:                          # pragma: no cover
        progress(f"  visualisation unavailable: {ex}")
        return
    viz = os.path.join(out_dir, "viz")
    os.makedirs(viz, exist_ok=True)
    # arc plots of the target genes (trigger accessibility context)
    for name, seq in (viz_genes or {}).items():
        try:
            arc_plot(seq, os.path.join(viz, f"{name}_arcs.png"), backend,
                     title=f"{name} base-pairing")
            export_pair_fraction_csv(seq, os.path.join(viz, f"{name}_pair_fraction.csv"),
                                     backend)
        except Exception as ex:                      # pragma: no cover
            progress(f"  arc plot for {name} failed: {ex}")
    # arc plots of the top designs' switches with domain shading
    for r in out.results[:min(3, len(out.results))]:
        sw = r.switch
        regions = {k: v for k, v in sw.domains.spans.items()
                   if k in ("sec_k2star", "sec_xstar", "prim_k1star", "prim_loop")}
        try:
            arc_plot(sw.core, os.path.join(viz, f"rank{r.rank}_switch_arcs.png"),
                     backend, title=f"rank {r.rank} switch (OFF state)",
                     regions=regions)
        except Exception as ex:                      # pragma: no cover
            progress(f"  arc plot for rank {r.rank} failed: {ex}")
