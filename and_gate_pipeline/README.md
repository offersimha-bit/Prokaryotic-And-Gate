# AND-gate toehold-switch design pipeline

A Python pipeline that designs and ranks **two-input RNA AND-gate toehold
switches** from two natural genes. Trigger B opens an upstream *inhibitory*
stem; Trigger A then opens the downstream *main* (Series-A) stem — so output is
produced only when **both** triggers are present (Kim 2019 sequential-hairpin
logic; Green 2026 / Toehold-VISTA Series-A architecture and accessibility
scoring).

It implements the full specification: target scanning → thermodynamic
filtering → AND-gate construction → constraint checking → multi-stage scoring
and ranking → off-target scanning → arc-plot visualisation.

---

## Setting up (for teammates)

This repo uses a **git submodule** (`vista/` — the AlexGreenLab Toehold-VISTA
reference implementation, pinned to the exact commit this pipeline was built
against). Submodules are *not* fetched by a plain `git clone`, so use:

```bash
git clone --recurse-submodules https://github.com/offersimha-bit/Prokaryotic-And-Gate.git
cd Prokaryotic-And-Gate

# already cloned without --recurse-submodules? run this once:
git submodule update --init --recursive
```

Then create an environment and install the pinned dependencies:

```bash
python -m venv .venv && source .venv/bin/activate      # Linux/WSL
python -m pip install -r and_gate_pipeline/requirements.txt
python -m and_gate_pipeline --demo --out results       # verify
python -m and_gate_pipeline.tests.test_pipeline        # 14/14 should pass
```

NUPACK is optional and must be installed separately (see below); without it the
pipeline runs on ViennaRNA automatically. **The pipeline does not need the
`vista/` submodule at runtime** — the one file it used from there (the E. coli
codon-usage table) is vendored at `and_gate_pipeline/data/`. The submodule is
pinned for reproducibility and reference (notebooks, PLS-DA model params).

## Requirements

| Package | Role | Verified version |
|---|---|---|
| **ViennaRNA** (`import RNA`) | folding, partition function, base-pair probabilities | ✅ 2.7.2 |
| **NUPACK 4** (`import nupack`) | preferred engine (matches the VISTA reference model) | ✅ 4.1.0.1 — optional; used automatically when present |
| numpy, pandas, matplotlib, networkx | scoring tables + arc plots | ✅ pinned in `requirements.txt` |

> **Running under WSL.** NUPACK is installed in the Linux virtual environment
> `.venv` (Python 3.12), which lives one level **above** this repo (it is not
> committed). Run the pipeline from the repo root through it:
>
> ```bash
> wsl bash -c "cd '/mnt/c/Users/Dell/OneDrive - mail.tau.ac.il/IGEM/Toehold/Prokaryotic And Gate/Prokaryotic-And-Gate' \
>   && ../.venv/bin/python -m and_gate_pipeline --demo --out results"
> ```
>
> With NUPACK present the pipeline uses it by default; add `--no-nupack` to force
> ViennaRNA. If you clone this repo somewhere without that `.venv`, create one and
> `pip install -r and_gate_pipeline/requirements.txt` (NUPACK must be installed
> separately — see below); the pipeline then runs on ViennaRNA alone.
>
> **NUPACK vs ViennaRNA (verified).** On identical inputs the two engines agree
> to <1 kcal/mol on OFF-state / stem MFE and to <0.01 on SED/NED/accessibility;
> per-design total scores differ by ~0.08 on average and rankings are ~89%
> concordant. The one systematic difference is the multi-strand ON-state
> complex MFE (~5 kcal/mol, NUPACK more negative), reflecting different
> strand-association models. NUPACK's `defect()` already returns the *normalised*
> ensemble defect, and its pair matrix carries P(unpaired) on the diagonal — both
> handled in `thermo.py`. Constraint-conditioned accessibility (the AND-mechanism
> intermediate/ON sub-scores) has no NUPACK analysis hook, so it is delegated to
> ViennaRNA, whose unconstrained accessibility matches NUPACK to <0.01.

> **NUPACK note.** NUPACK is not pip-installable (licensed manual download).
> When it is present the pipeline uses it with the exact VISTA model
> (`Model(material='rna', ensemble='stacking', celsius=T, sodium, magnesium)`).
> When it is absent, every quantity — MFE, **SED** (specified ensemble defect),
> **NED** (native ensemble defect), accessibility, binding ΔG — is computed from
> ViennaRNA's partition function and base-pair-probability matrix instead. The
> ensemble-defect definition is identical (expected number of incorrectly
> paired nucleotides relative to a reference structure, normalised by length).
> No code change is needed to switch engines.

## Quick start

```bash
# self-contained demo on bundled example genes
python -m and_gate_pipeline --demo --out results

# your own genes (raw sequence or @fasta)
python -m and_gate_pipeline --gene1 @geneA.fasta --gene2 @geneB.fasta \
    --Lx 12 --LA 36 --LB 30 --secondary-loop 11 --out results

# force the ViennaRNA backend, skip plots
python -m and_gate_pipeline --demo --no-nupack --no-viz --out results
```

Outputs in `results/`:

* `and_gate_designs_ranked.csv` — every scored design, all sub-scores + details
* `final_designs.txt` — human-readable top-N with sequences and OFF-state structure
* `viz/*_arcs.png`, `viz/*_pair_fraction.csv` — arc diagrams of the target genes
  and the top switches (VISTA `pair_fraction.csv` layout)

## Pooled multi-gene discovery (`interop.py`)

The standalone scanner at `Triger finding/and_gate_trigger.py` does something
this pipeline doesn't: it pools **many** gene records and hunts for the
coincidence where one gene's connector `x` is the exact reverse complement of
another gene's `k2`. This pipeline takes only two genes but builds the real
switch. `interop.py` chains them — **pooled discovery → full design** — without
editing the scanner:

```python
from and_gate_pipeline.interop import load_scanner, load_genes, run_from_scanner, config_from_params

mod   = load_scanner()                       # imported unmodified; its main() never runs
cfg   = config_from_params(mod.Params())     # Lx=6, r1=11, a=4, k1=15, r2=30
genes = load_genes(["genes/amr.fasta"])      # any number of FASTA files, pooled
out   = run_from_scanner(genes, cfg=cfg, out_dir="results")

best = out.results[0]
print(best.score.quality_percent, best.switch.core)
print(best.pair.meta["gene_a_name"], best.pair.meta["gene_b_name"])
```

> **Why this is safe to run on NUPACK.** The scanner's own `select_backend()`
> tries NUPACK first, and its NUPACK `unpaired_probs()` assumes an
> `(n+1)×(n+1)` pair matrix with an unpaired column — but NUPACK 4.1 returns
> `n×n` with P(unpaired) on the **diagonal**. The lookup raises `IndexError`,
> which the scanner swallows (`except Exception: return [0.5]*n`), silently
> reporting **every base as 50 % unpaired** — which also scrambles its
> accessibility-ordered top-k pre-selection. `interop` never calls
> `select_backend()`: it injects `_ScannerBackendShim`, wrapping this pipeline's
> verified engine, into the scanner's `scan_trigger1`/`scan_trigger2` (which
> take `backend` as an argument). Covered by
> `test_interop_shim_is_not_the_broken_backend`.
>
> ⚠️ Running `python and_gate_trigger.py` **standalone** still hits the bug.
> Until it's fixed upstream, set `backend: str = "vienna"` in its `Params`.

Ported from that scanner into this pipeline's scoring (Section 7D):
**cross-trigger crosstalk** (`crosstalk_stick_nt` / `crosstalk_subst_nt` — do the
two triggers hybridise to or mimic each other beyond the intended connector,
with the connector masked out), **Type IIS sites** (BsaI/BsmBI/SapI/BbsI —
Golden Gate hazards), and an **absolute 0–100 `quality_percent`** with a
per-criterion `breakdown()`, measured against fixed reference scales so it stays
comparable across runs.

## Library use

```python
from and_gate_pipeline import PipelineConfig
from and_gate_pipeline.pipeline import run_pipeline

cfg = PipelineConfig(Lx=12, L_A=36, L_B=30, secondary_loop_len=11)
out = run_pipeline(gene1, gene2, cfg, reporter=gfp_cds, out_dir="results")
best = out.results[0]
print(best.switch.core, best.score.total)
```

---

## How the stages map to the spec

| Stage | Module | What it does |
|---|---|---|
| 1. Target scan / triggers | `target_scan.py` | finds length-`Lx` reverse-complement cores in G1/G2 (exact, else **minimum-Hamming** fallback); builds Trigger A = `r1·x·a·k1` and Trigger B = `r2·k2`; runs both role-swap orientations |
| 2. Thermo filtering | `filtering.py` | MFE + SED + NED + accessibility of each trigger over ±0/10/25/50/100 nt windows; ±100 nt gate |
| 3. Architecture | `architecture.py` | builds the two-hairpin switch (secondary inhibitory stem + spacer `a*` + Series-A 18 bp primary stem, RBS loop, AUG bulge) and the intended OFF-state structure |
| 4. Constraints | `constraints.py` | checks the Section-4 equations and logical integrity |
| 5–7. Scoring | `scoring.py` | hierarchical score: (A) Trigger-B activation, (B) post-B intermediate, (C) Trigger-A ON state, (D) penalties |
| — Off-target | `offtarget.py` | transcriptome-wide sliding-window complementarity scan; essential-gene hits disqualify |
| 6. Optimisation | `optimize.py` | repairs forbidden runs / in-frame stops / extra AUGs and nudges OFF-MFE toward −54.25 |
| Visualisation | `visualize.py` | networkx + matplotlib arc plots and `pair_fraction.csv` export |
| Engine | `thermo.py` | NUPACK↔ViennaRNA backend abstraction |

## Tunable variables (Section 5.2)

All live on `PipelineConfig`: `secondary_loop_len`, `secondary_arm_gc_bias`
(r1-clamp strength), `L_A`, `L_B`, `len_r2`, plus `Lx`, `len_a`, `len_k1`,
`primary_stem_len`, `off_state_mfe_target`, scoring `weights`, and the physical
model (`temperature_c`, `sodium`, `magnesium`). CLI flags cover the common ones.

## Scoring summary (Section 7)

* **7A Trigger-B activation** — target-region accessibility (G2), switch toehold
  availability, and Trigger-B : inhibitory-stem binding ΔG; optional
  expression/encounter weighting.
* **7B Intermediate** — with Trigger B's toehold held open (ViennaRNA hard
  constraint), re-measures Trigger-A-site accessibility (must *increase*) and
  intermediate-complex stability (kinetic-trap check).
* **7C Trigger-A / ON** — ON-state MFE of the ternary complex, RBS/AUG
  liberation when both triggers are bound, and codon-usage translational
  efficiency of the first codons after the start.
* **7D Penalties** — leakage vs −54.25, secondary-must-be-stronger-than-primary,
  weak spacer `a`, forbidden runs / in-frame stops / spurious AUGs, off-target
  hits, and a structural half-life proxy.

## Design decisions & caveats

* **Exact-complement construction (Section 6).** Every switch domain is the
  exact reverse complement of the *actual* trigger domain, so each trigger is
  captured with perfect complementarity even when `x` and `k2` match only
  approximately. The residual `x`/`k2` mismatch appears as `hamming` mismatches
  inside the secondary stem — the true biophysical cost of using natural genes.
* **What is asserted vs measured.** The OFF-state *lock* is built explicitly per
  Section 3 (5′ arm `k2*·r1`, 3′ arm `r1*·x*`, 3-nt junction bulge) and scored
  with SED against that intended structure. Trigger-binding steps are **not**
  assumed — they are evaluated with real cofold thermodynamics and
  constraint-conditioned accessibility, so the numbers stay physical regardless
  of annotation. The sequential (two-contact) opening is exactly the Kim-2019
  AND mechanism.
* **Performance.** The off-target scan is an O(N·L) sliding window; full scoring
  is applied only to the top `--max-full-score` pre-ranked candidates.
* The primary loop reuses a conserved RBS/AUG context so the reading frame and
  single start codon are valid by construction; the optimiser only edits the
  free primary-stem body.

## Tests

```bash
python -m pytest and_gate_pipeline/tests -q
```
```

## Sources

* Kim, J. et al. (2019) *Modulating responses of toehold switches by an
  inhibitory hairpin* — inhibitory-hairpin logic and the short spacer `a`.
* Green, A. et al. (2026) / **Toehold-VISTA** — Series-A stem (18 bp / 6 nt
  invasion), SED/NED accessibility features, ±100 nt flanking emphasis,
  −54.25 kcal/mol low-leak OFF-state target.
