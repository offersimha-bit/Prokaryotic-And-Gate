# HANDOFF — AND-gate design pipeline

State as of commit `13e7ab1`. Written for someone opening this repo cold.

---

## 1. What we are building

A **two-input RNA AND gate** on **endogenous** genes. Feed the pipeline a gene
list (e.g. from differential expression), and it looks for two genes that
together — and only together — switch on a reporter.

Architecture follows **Kim 2019**: an upstream *inhibitory* hairpin masks the
binding site of the trigger that opens the downstream *primary* hairpin (which
carries the RBS and AUG). Trigger B unmasks; Trigger A then fires. The primary
hairpin is **Green 2026 / VISTA Series-A** geometry (18-nt arms, 15 bp, 3-nt
AUG bulge, 6-nt invasion, 30-nt toehold).

Both triggers are **real contiguous stretches of real genes** — nothing is
synthesised. That is the hard constraint the whole design bends around.

---

## 2. ⚠️ Read this before running anything: the repo has two architectures

| | **Track A — what the CLI runs** | **Track B — what we validated** |
|---|---|---|
| builder | `architecture.build_switch` | `vista_switch.build` |
| scoring | `scoring.DesignScorer` (hand-invented 0–100 weights) | `kinetics.py` (rate vs mRNA decay) |
| imported by | `pipeline.py`, `cli.py`, `optimize.py`, `spec_audit.py` | `sweep.py`, `truth_table.py` **only** |
| status | **architecture we disproved** | correct, but **not wired into the pipeline** |

`python -m and_gate_pipeline --demo` exercises **Track A**. Every conclusion in
§4 comes from **Track B**. Wiring Track B into `pipeline.py` is the top of the
to-do list (§6).

---

## 3. How to run it

NUPACK only works inside the **WSL** virtualenv, which sits one level *above*
this repo and is not committed.

```bash
wsl bash -c "cd '/mnt/c/Users/Dell/OneDrive - mail.tau.ac.il/IGEM/Toehold/Prokaryotic And Gate/Prokaryotic-And-Gate' \
  && ../.venv/bin/python -m and_gate_pipeline --demo --out results"
```

| command | what it does |
|---|---|
| `-m and_gate_pipeline --demo` | full pipeline, demo genes → ranked CSV + arc plots (**Track A**) |
| `-m and_gate_pipeline.tests.test_pipeline` | 33 tests, both engines (pytest optional) |
| `-m and_gate_pipeline.spec_audit` | 44 spec demands → PASS / NOTE / DECISION / FAIL |
| `-m and_gate_pipeline.sweep` | L_x and \|a\| sweeps, kinetically scored (**Track B**) |
| `-m and_gate_pipeline.truth_table` | 4-condition truth table for one design (**Track B**) |

Fresh clone needs `--recurse-submodules` (VISTA is pinned as a submodule at
`14eef32` and *is* imported at runtime). Deps pinned in
`and_gate_pipeline/requirements.txt`; NUPACK must be installed separately.

---

## 4. What we measured

All on **one synthetic demo gene pair** — never yet on real genes.

**Two spec errors, both fixed in Track B:**

- **Trigger A domain order.** The spec said `r1-x-a-k1` (5'→3'). The switch's
  own footprint requires **`k1-a-x-r1`** — the spec order is reversed and
  cannot lay down on its own binding site. Fixing it: ΔG(A:switch)
  **−4.31 → −30.82 kcal/mol**.
- **Trigger B had no toehold.** The spec put `k2*` inside the stem's 5' arm
  with nothing single-stranded upstream. Adding the **`r2*`** toehold is what
  lets B bind at all.

**Governing equation** (measured, not assumed):

```
toehold available to Trigger A  =  |a|         in OFF
                                =  |a| + Lx    after Trigger B binds
```

B can only ever displace the `Lx`-nt `k2*:x*` helix — never `r1:r1*`. So `|r1|`
sets the OFF lock and the post-nucleation energy, but **not** the toehold A gets
to grab. Only `|a|` and `Lx` move that.

**Parameter findings:**

| finding | evidence |
|---|---|
| **L_x = 7 optimal** | effective ON 71%; agreed by two independent metrics |
| **\|a\| = 4 is a knife-edge** | \|a\|=3 → dead switch; \|a\|=5 → 85% leak. Reproduces Kim's own choice **from physics**, not by fitting |
| **Secondary-stem mismatches are a dead lever** | one mismatch *anywhere* kills the AND — the helix that hides the toehold is the helix that holds it shut. Keep `secondary_arm_gc_bias = 0` |
| **±1 nt sensitivity** | rate is ~1 decade per nucleotide. Synthesis error or a mis-annotated TSS is not cosmetic |

**Kinetic result** (L_x=7, \|a\|=4):

| | OFF | after B |
|---|---|---|
| toehold ΔG | −5.34 | −11.16 kcal/mol |
| time to fire | 17,177 s | 35 s (mRNA t½ = 300 s) |
| P(fire before decay) | **2.46%** | **92.58%** |

**AND ratio 38×.** Inside the published range for toehold switches.

---

## 5. Why scoring is kinetic, not equilibrium

Two metrics disagreed **~700×** on the same design's OFF leak (equilibrium 32%
vs kinetic 2.5%). Equilibrium asks *"given infinite time, what fraction binds?"*
— and Trigger A **can** eventually prise its site out of the inhibitory hairpin.
A cell never grants infinite time: the transcript is degraded first. **mRNA
lifetime is a kinetic filter equilibrium cannot see.**

`kinetics.py` implements Zhang & Winfree (2009):

```
k_eff  = k_on·k_bm / (k_on·Kd_toe + k_bm)
P_fire = k_obs/(k_obs + k_deg),   k_obs = k_eff·[trigger],  k_deg = ln2/t½
```

`Kd_toe` uses the **accessibility-corrected** energy
`ΔG_duplex + opening_energy(toehold)` — which is what unifies the two failed
metrics: accessibility enters the *rate*, and the rate races degradation.

> **Caveat that must travel with any number from this model.** `k_on`/`k_bm`
> are order-of-magnitude **DNA** values reused for RNA at 37 °C, and `k_bm`
> lumps a length-dependent random walk into one constant. **Absolute P_fire is
> indicative only — rank designs by the ratio between them; never quote P_fire
> as a yield.**

---

## 6. What is not done

1. **Wire Track B into `pipeline.py`** (§2). Until then the CLI builds the
   disproved architecture.
2. **Use VISTA's ML model for ranking.** `all_trained_model_params.pkl` holds a
   trained PLS-DA model; we currently use VISTA's *scanner* but not its model.
   Needs the scaler fix below.
3. **Run on real genes.** Everything above is one synthetic pair.
4. **Step 3 sweep:** `|r1|`, `|r2|`, secondary loop size/strength — unblocked
   now that the OFF lock is scorable.
5. **Migrate `truth_table.py` to `kinetics.py`.** It still uses the superseded
   nucleation-only `occupancy()`, so **the repo currently prints two different
   AND ratios for the same design**: `truth_table` says 78×, `kinetics`/`sweep`
   say 38×. The 38× is the one to trust (§5). Not a disagreement about the
   design — just two modules on two metrics.
5. **Skip-x variant:** drop the `x = revcomp(k2)` constraint entirely and let
   Trigger A loop out over `k2`. Removes the gene–gene coincidence requirement
   at the cost of a large internal loop. Untested.

**Known bug in VISTA's own code (not ours):** `rank_new_designs()` recomputes
the scaler from whatever pool you pass it, while the trained
`scaler_mean`/`scaler_scale` sit unread in the pkl. Its rankings are therefore
pool-relative. Use the saved scaler.

**Known bug in `Triger finding/and_gate_trigger.py` (teammate's file):** its
NUPACK `unpaired_probs()` hits an IndexError that a bare `except` swallows,
silently reporting every base as 50% unpaired — and it prefers NUPACK by
default, so it broke the moment NUPACK was installed. Standalone workaround:
set `backend="vienna"` in its `Params`. Our `interop.py` is immune (it injects
our verified engine instead of calling `select_backend()`).

---

## 7. Module map

| module | role |
|---|---|
| `target_scan.py` | find `x`/`k2` pairs across two genes; min-Hamming fallback; both orientations |
| `filtering.py` | fold triggers in native context ±{0,10,25,50,100} nt; SED/NED/accessibility gate |
| `vista_switch.py` | **Track B** builder — Kim masking + VISTA's own switch builder |
| `architecture.py` | Track A builder (superseded) |
| `kinetics.py` | **the scoring model** — displacement rate vs mRNA decay |
| `thermo.py` | NUPACK ↔ ViennaRNA behind one interface |
| `scoring.py` | Track A scoring, hand-invented weights (superseded by `kinetics.py`) |
| `sweep.py` / `truth_table.py` | parameter sweeps and the 4-condition AND check |
| `interop.py` | bridge to the pooled multi-gene trigger scanner |
| `spec_audit.py` | every spec clause checked against a real design |
| `offtarget.py`, `optimize.py`, `visualize.py`, `constraints.py` | off-target scan, sequence repair, arc plots, length equations |

Engines cross-checked: <1 kcal/mol on stem MFE, ~89% rank concordance. The one
real divergence is multi-strand ON-state complex MFE (~5 kcal/mol).
