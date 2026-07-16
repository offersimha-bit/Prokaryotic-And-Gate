"""Base-pairing arc-plot visualisation (integrates the networkx workflow).

Reproduces the VISTA ``pair_probability_arc_diagram`` style: pairwise
probabilities are turned into a networkx graph, then drawn as arcs along a
linear backbone -- high-probability pairs above the axis, low-probability pairs
below -- coloured by probability.  Also exports the pairwise probabilities in
the VISTA ``pair_fraction.csv`` layout so the original notebook can consume
them unchanged.

Use it to show base-pairing across a full RNA target (trigger accessibility) or
across a designed switch, optionally shading domain/trigger regions.
"""

from __future__ import annotations

import csv

from . import sequence_utils as su
from .thermo import ThermoBackend


def export_pair_fraction_csv(seq: str, path: str, backend: ThermoBackend,
                             threshold: float = 0.01) -> str:
    """Write pairwise probabilities in the VISTA ``pair_fraction.csv`` layout:
    columns (no header) = StrandA, i, StrandB, j, Probability; unpaired rows use
    j = -1."""
    seq = su.to_rna(seq)
    pairs = backend.pair_probabilities(seq, threshold)
    unpaired = backend.unpaired_probabilities(seq)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        for i, j, p in pairs:
            w.writerow([0, i, 0, j, f"{p:.6f}"])
        for i, up in enumerate(unpaired, start=1):
            w.writerow([0, i, 0, -1, f"{up:.6f}"])
    return path


def arc_plot(seq: str, path: str, backend: ThermoBackend,
             threshold: float = 0.02, title: str | None = None,
             regions: dict | None = None, figsize=(24, 8)):
    """Render an arc diagram of base-pairing probabilities to ``path`` (PNG).

    ``regions``: optional {label: (start, end)} 0-based spans shaded on the
    backbone (e.g. trigger footprint, RBS)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import networkx as nx

    seq = su.to_rna(seq)
    n = len(seq)
    pairs = backend.pair_probabilities(seq, threshold)

    G = nx.Graph()
    G.add_nodes_from(range(1, n + 1))
    for i, j, p in pairs:
        G.add_edge(i, j, weight=p)

    fig, ax = plt.subplots(figsize=figsize)
    if pairs:
        cmap = plt.cm.viridis.reversed()
        probs = [p for _, _, p in pairs]
        norm = plt.Normalize(vmin=min(probs), vmax=max(probs))
        for i, j, p in pairs:
            span = abs(j - i)
            height = span / n
            xs = np.linspace(i, j, 100)
            sign = 1.0 if p > 0.5 else -1.0
            ys = sign * height * np.sin(np.linspace(0, np.pi, 100))
            ax.plot(xs, ys, color=cmap(norm(p)),
                    linewidth=1.0 + 2.5 * p, alpha=0.85)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.02)
        cbar.set_label("pairing probability")

    ax.axhline(0, color="black", linewidth=0.8, zorder=1)

    if regions:
        colors = plt.cm.tab10.colors
        for k, (label, (s, e)) in enumerate(regions.items()):
            ax.axvspan(s + 1, e, ymin=0.47, ymax=0.53,
                       color=colors[k % len(colors)], alpha=0.6)
            ax.text((s + 1 + e) / 2, 0, label, ha="center", va="center",
                    fontsize=8, rotation=90)

    ax.set_xlim(0, n + 1)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel("nucleotide position")
    if title:
        ax.set_title(title)
    ax.text(0.01, 0.95, "above axis: P>0.5   below: P<=0.5",
            transform=ax.transAxes, fontsize=8, color="gray")
    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path
