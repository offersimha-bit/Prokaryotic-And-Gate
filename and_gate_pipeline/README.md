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
python -m and_gate_pipeline map                        # stage <-> file <-> alias
python -m and_gate_pipeline.tests.test_pipeline        # 40/40 should pass
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

## Pooled multi-gene discovery (now part of STAGE 1)

The standalone scanner at `Triger finding/and_gate_trigger.py` used to be the
only way to pool many gene records and hunt for the coincidence where one gene's
connector `x` is the exact reverse complement of another gene's `k2`. That
capability now lives inside `01_target_scan.py`, and `interop.py` is gone.

```python
from and_gate_pipeline import PipelineConfig
from and_gate_pipeline.target_scan import scan_from_fasta
from and_gate_pipeline.pipeline import run_pipeline

cfg   = PipelineConfig()
pairs = scan_from_fasta("genes/", cfg)          # any number of FASTA files
out   = run_pipeline(cfg=cfg, pairs=pairs, out_dir="results")
```

### Why the bridge was removed rather than fixed

`interop.windows_to_pair()` handed the scanner's slices straight into
`TriggerA`. The scanner cut its window in genomic order `r1 | x | a | k1`, while
`TriggerA.seq` reassembles `k1 + a + x + r1` — the order the switch requires.
The result was a "trigger" that **does not occur in the gene**, at a `pos_x`
that made stage 2 fold a different window than the one selected:

```
scanner window (contiguous in gene): GUGAAGAACUGUUUACCGGCGUGGUGCCGAUUCUGG
TriggerA.seq as the pipeline used it: UGGUGCCGAUUCUGGGGCGUUUACCGUGAAGAACUG
is TriggerA.seq present in gene?      False
scanner A window span in gene:        10 -> 46
filtering.py folded span:              2 -> 38
```

There is now exactly one slicer, one coordinate convention, and
`target_scan.verify_pair()` asserts on every returned pair that both triggers
are literally substrings of their source genes at the recorded offsets.
`test_trigger_is_a_contiguous_slice_of_its_gene` and
`test_verify_pair_rejects_a_scrambled_trigger` keep it that way.

### What was taken from the scanner, and what was not

| Taken into `01_target_scan.py` | Why |
|---|---|
| pooled FASTA input (`read_fasta_records`, `load_genes`) | the pipeline only handled two genes |
| the hash-join on the connector | generalised to a k-mer index: **~170× faster** on 4 kb × 4 kb (0.07 s vs ~12 s) |
| the two-distinct-records rule | an AND gate needs two inputs, not one gene sensing itself twice |
| GC / motif window filters | kept, but **off by default** — see below |

| Taken into `03_select.py` | Why |
|---|---|
| criteria 1, 2, 3, 5, 6 | genuine properties of the trigger itself |
| the **fixed reference scale** | the scanner's best decision: a quality of 0.8 means the same thing across runs, unlike min-max over a pool of 8 |
| the per-criterion breakdown | scoring you cannot audit is scoring you cannot trust |

| Deliberately **not** taken | Why |
|---|---|
| `build_switch_target` and criteria 4 & 7 | they assume the switch is `revcomp(r1+r2+a+k1)`, i.e. r1 adjacent to r2. In this architecture r2 binds `r2*` on the inhibitory hairpin while r1 binds `r1*` in the primary toehold, separated by `k2*`, an 11-nt loop and `x*`. Those criteria scored a molecule that is never built; the question they asked is answered kinetically in stage 5 |
| `unpaired_probs(whole_gene)` | a global fold of a full transcript. Accessibility is stage 2's job, on flanked windows |
| Type IIS filtering of the **trigger** | the trigger is endogenous and never synthesised; Golden Gate sites only matter for the switch. `cfg.scan_forbid_motifs` defaults to `False` |
| the NUPACK `unpaired_probs` path | its `(n+1)×(n+1)` assumption raises `IndexError`, swallowed by a bare `except` that returns `[0.5]*n`. `thermo.py` reads the diagonal correctly instead |

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

## Repository layout — the file names carry the stage number

```
and_gate_pipeline/
  01_target_scan.py            STAGE 1  find trigger pairs (pooled or two-gene)
  02_filtering.py              STAGE 2  accessibility of each trigger in context
  03_select.py                 STAGE 3  Pareto + diversity shortlist
  04_build_switch.py           STAGE 4  Kim hairpin + VISTA Series-A module
  04_build_switch_legacy.py    STAGE 4  superseded builder (kept for comparison)
  04_optimize.py               STAGE 4  restricted-sequence repair
  05_kinetics.py               STAGE 5  displacement rate vs mRNA decay
  05_scoring_legacy.py         STAGE 5  superseded hand-weighted scorer
  05_truth_table.py            STAGE 5  four-condition check for one design
  05_sweep.py                  STAGE 5  L_x and |a| sweeps
  06_offtarget.py              STAGE 6  transcriptome complementarity scan
  06_visualize.py              STAGE 6  arc plots

  config.py  thermo.py  sequence_utils.py  constraints.py  examples.py
                               shared — used by several stages, no stage number
  pipeline.py                  orchestrator across all stages
  cli.py  spec_audit.py        entry point and spec checker
  _loader.py                   makes the numbered files importable
```

A Python identifier cannot start with a digit, so `import 01_target_scan` is a
syntax error. `_loader.py` therefore loads each numbered file by path and
registers it under a stable alias, in dependency order. **No import statement in
the package or the tests changed** — the numbers are for humans reading the
directory, the aliases are for Python:

```python
from and_gate_pipeline.target_scan import scan_from_fasta   # 01_target_scan.py
from and_gate_pipeline.kinetics import and_behaviour        # 05_kinetics.py
```

`python -m and_gate_pipeline map` prints the file ↔ alias ↔ stage table.

## Command line

Because the stage modules are numbered, `python -m and_gate_pipeline.truth_table`
no longer resolves. Those entry points are subcommands:

```bash
python -m and_gate_pipeline map                        # file <-> stage <-> alias
python -m and_gate_pipeline run --demo --out results   # full design run
python -m and_gate_pipeline scan genes/ --exact        # STAGE 1 only, pooled FASTA
python -m and_gate_pipeline truth-table                # STAGE 5, one design
python -m and_gate_pipeline sweep                      # STAGE 5, parameter sweeps
python -m and_gate_pipeline audit                      # spec audit
python -m and_gate_pipeline --demo --out results       # `run` is the default
```

## What each stage does

| Stage | File | What it does |
|---|---|---|
| 1. Find triggers | `01_target_scan.py` | pooled multi-gene **or** two-gene scan; k-mer index join on `x == revcomp(k2)`; min-Hamming fallback; both role orientations; asserts every trigger is a contiguous slice of its gene |
| 2. Accessibility | `02_filtering.py` | MFE + SED + NED + accessibility of each trigger over ±0/10/25/50/100 nt windows in its native transcript |
| 3. Selection | `03_select.py` | four absolute criteria (accessibility, self-fold openness, fold stability, cross-talk) → Pareto fronts → greedy pick with ≤50 % window overlap |
| 4. Build | `04_build_switch.py` | Kim 2019 inhibitory hairpin prepended to VISTA's own Series-A builder |
| 4. Repair | `04_optimize.py` | forbidden runs / in-frame stops / extra AUGs; nudges OFF-MFE toward −54.25 |
| 5. Score | `05_kinetics.py` | accessibility-corrected ΔG → displacement rate → P(fire before decay) |
| 6. Off-target | `06_offtarget.py` | transcriptome-wide sliding-window complementarity; essential-gene hits disqualify |
| 6. Visualise | `06_visualize.py` | networkx + matplotlib arc plots, `pair_fraction.csv` export |
| — | `thermo.py` | NUPACK ↔ ViennaRNA behind one interface |

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
python -m pytest and_gate_pipeline/tests -q      # or, without pytest installed:
python -m and_gate_pipeline.tests.test_pipeline  # 40/40
```

## Sources

* Kim, J. et al. (2019) *Modulating responses of toehold switches by an
  inhibitory hairpin* — inhibitory-hairpin logic and the short spacer `a`.
* Green, A. et al. (2026) / **Toehold-VISTA** — Series-A stem (18 bp / 6 nt
  invasion), SED/NED accessibility features, ±100 nt flanking emphasis,
  −54.25 kcal/mol low-leak OFF-state target.
