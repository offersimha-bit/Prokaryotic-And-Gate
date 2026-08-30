# What every metric in the candidates PDF actually measures

Companion to `Toehold_Candidates29.7.pdf`. For each row the document reports:
what it means in plain terms, exactly how it is computed, whether our code
reproduces it, and what to watch out for when reading it.

Implementation: [`pdf_metrics.py`](pdf_metrics.py). Checked by
[`validate_pdf.py`](validate_pdf.py) (stage 0).

Engine throughout: **ViennaRNA 2.x, default model, 37 °C**. This is not a
preference — it is what the PDF used, confirmed by the fact that
`RNA.fold()` and `RNA.cofold()` return its energies to the last decimal.
NUPACK disagrees by 3–5 kcal/mol on the two-strand complex, which is why it
belongs in a second column rather than as a substitute.

---

## Reproduction status

| # | PDF row | status |
|---|---|---|
| 1 | Similarity score | **not wired** (VISTA model) |
| 2 | ΔG switch alone | exact |
| 3 | ΔG switch + trigger complex | exact |
| 4 | ΔG margin | exact |
| 5 | Toehold open — OFF | exact |
| 6 | RBS hidden — OFF | advisory (±0.6) + **mislabelled** |
| 7 | RBS exposed — ON | advisory (±0.5) |
| 8 | Start codon exposed — ON | exact |
| 9 | AUG bulge in OFF (MFE structure) | exact |
| 10 | Region just after start codon open — ON | **unresolved** |
| 11 | Trigger binding (at equilibrium) | exact |
| 12 | RBS present | exact |
| 13 | Reading frame intact | exact |
| 14 | Unwanted start / stop codons | exact |
| 15 | Trigger unfolded on its own | exact |
| 16 | Off-target vs engineered mCherry | exact (energy **and** percentage) |

50 of 50 strict values reproduce. Details on the three exceptions below.

---

## The energy rows

### ΔG switch alone
**Meaning.** How tightly the switch folds up on itself with no trigger around.
More negative = a more stable OFF hairpin.
**Computed.** `RNA.fold(switch)` → the minimum free energy.
**Watch out.** This is the single most stable structure, not the average one.
For candidate 4 that structure accounts for only **2.2%** of the ensemble, so
"the" OFF structure is a much looser idea than a single picture suggests.

### ΔG switch + trigger complex
**Meaning.** How stable the switch and trigger are once bound together.
**Computed.** `RNA.cofold(trigger + "&" + switch)`.

### ΔG margin (net drive to bind)
**Meaning.** The document's measure of how strongly the trigger wants to bind.
**Computed.** `complex − switch`.
**Watch out — two problems.**
1. It leaves out the trigger's own folding energy. The trigger has to come
   unfolded before it can invade, and that costs 1.3–6.2 kcal/mol here. The
   honest figure is `complex − switch − trigger`, which we report as
   **ΔG margin corrected**.
2. **The ranking it produces is mostly just trigger length.** The five triggers
   are 35–41 nt and a longer one binds harder for free. Correlation between
   length and margin is **r = −0.872**. Per nucleotide the five span only 12%
   where the raw column spans 25%, and the top two swap places:

   | ranking | order |
   |---|---|
   | by raw margin | 4, 5, 3, 2, 1 |
   | by margin per nt | **5**, **4**, 3, 2, 1 |

   So always read the per-nucleotide column alongside the raw one.

There is a third oddity worth knowing: the PDF uses the **corrected** margin as
the denominator of its off-target percentage (row 16) while printing the
**uncorrected** one in this row. Its own two rows disagree about which margin is
the real one. We know this because using the corrected margin reproduces the
percentages exactly and using the printed one does not.

---

## The accessibility rows

These all rest on one idea: fold the molecule, then ask for each base how
likely it is to be **unpaired** across the whole ensemble of structures — not
just in the single best one.

### Toehold open — OFF
**Meaning.** How exposed the trigger's landing site is before the trigger
arrives. If the toehold is buried in the switch's own structure the trigger has
nothing to grab.
**Computed.** Mean unpaired probability over the toehold span, switch folded
alone.
**Reading it.** All five sit at 50–59%, which the PDF fairly calls modest.

### RBS hidden — OFF  ← mislabelled
**Meaning as computed.** The mean probability the RBS is **unpaired**, i.e. how
**accessible** it is. That is the opposite of "hidden".
**What the numbers say once read correctly:**

| cand | RBS unpaired OFF | RBS unpaired ON | change |
|---|---|---|---|
| 1 | 87.6% | 94.8% | +7.2 |
| 2 | 97.1% | 97.6% | +0.5 |
| 3 | 97.1% | 97.6% | +0.5 |
| 4 | 94.6% | 97.0% | +2.4 |
| 5 | 96.9% | 95.2% | **−1.7** |

The RBS is highly accessible in **both** states and barely moves — candidate
5's actually goes down. **This row cannot rank candidates.**

**Is that a design failure? No.** In this architecture the RBS deliberately sits
in an 18-nt hairpin **loop**, and loops are unpaired by construction. Repression
comes from the start codon being sequestered in the stem, plus the loop being
too constrained for a ribosome to engage productively — not from the RBS being
base-paired. The number is doing what the architecture intends; the *label* and
the verdicts built on it ("strong RBS-off", "excellent RBS protection") read it
backwards.

**Why advisory rather than exact.** The PDF never states which window it
averaged. The 11-nt `AACAGAGGAGA` lands within ~0.6 points, the 8-nt
`AGAGGAGA` within ~0.8, neither exact. Since the row is excluded from ranking
anyway, forcing agreement would mean fitting to a number we do not use.

### Start codon (AUG) exposed — ON
**Meaning.** How free the start codon is once the trigger has bound — i.e.
whether a ribosome can actually start.
**Computed.** Mean unpaired probability of the three AUG bases in the
switch+trigger complex.
**This is the row that matters.** Unlike the RBS it moves sharply between
states:

| cand | AUG off | AUG on | gap |
|---|---|---|---|
| 1 | 79.2% | 98.5% | +19.3 |
| 2 | 73.1% | 99.0% | +25.9 |
| 3 | 67.0% | 95.4% | **+28.3** |
| 4 | 70.0% | 96.0% | +26.0 |
| 5 | 77.7% | 90.6% | **+12.9** |

Candidate 5 has the weakest ON/OFF discrimination of the five on the only
accessibility metric that discriminates — the opposite of what the RBS row
implied about it.

### AUG bulge in OFF (MFE structure)
**Meaning.** In the single most stable OFF structure, how many of the start
codon's three bases are paired. The PDF calls 0/3 "clean" and 1/3 "acceptable".
**Computed.** Fold the switch, read the pair table, count.
**Watch out.** This looks at the MFE structure alone, and the MFE can be a
small slice of the ensemble (2.2% for candidate 4). The ensemble columns we add
exist to put exactly this kind of single-structure claim in context.

### Region just after the start codon open — ON  ← unresolved
The PDF prints 3.4 / 0.6 / 1.7 / 1.8 / 4.6 and **we cannot reproduce it**.
Definitions tried and rejected:

- mean unpaired probability of a window after the AUG in the ON state → 88–95%,
  far too high;
- mean **paired** probability of the same window, lengths 3 to 21 → best fit at
  length 6, still 16 points of total error across five candidates;
- joint probability that the whole window is unpaired, via a constrained
  partition function, lengths 3 to 15 → no length reproduces the pattern;
- the same on the `b_pre*` span, on the linker span, and in the OFF state →
  none match.

The code returns the closest candidate (mean paired probability over the 6 nt
after the start codon), clearly flagged, and **uses it nowhere**. Settling it
needs the definition from whoever wrote the original validation script.

---

## The binding rows

### Trigger binding (at equilibrium)
**Meaning.** How completely the trigger lands on its switch rather than
drifting free or folding on itself.
**Computed.** In the two-strand ensemble, the mean probability that a trigger
base is paired **to the switch specifically** (pairs where both partners lie on
the trigger are self-structure and are not counted).
**Reading it.** 98–99.8% for all five — high, and not very discriminating.

### Trigger unfolded on its own
**Meaning.** How much of the trigger is free rather than tangled in its own
structure. A trigger folded up on itself has to pay to open before it can
invade.
**Computed.** Mean unpaired probability of the trigger folded alone.
**Reading it.** 60–71%. Candidate 4's trigger is the most self-structured
(60.6%), which partly offsets its otherwise strong numbers.

---

## The off-target row

### Off-target vs engineered mCherry (duplex energy)
**Meaning.** Would the trigger rather stick to the reporter transcript than to
its own switch? The reporter here is the engineered (codon-max) mCherry.
**Computed.** `RNA.duplexfold(trigger, engineered_mCherry)` — the best duplex
between the two, ignoring each one's internal folding.
**The percentage.** `|off-target ΔG| / |corrected margin|`, where corrected
margin is `complex − switch − trigger`. This reproduces the PDF's 52 / 52 / 49 /
49 / 43 exactly.
**Reading it.** 43–52% of the intended binding strength is not nothing. The PDF
says as much — "worth double-checking through another method before ruling it
out completely."

---

## The sequence-integrity rows

Cheap checks, all passing for all five.

| row | what it verifies |
|---|---|
| RBS present | the Shine–Dalgarno sequence is actually in the construct |
| Reading frame intact | distance from the start codon to the end divides by 3, so the reporter downstream is in frame |
| Unwanted start / stop codons | reading in frame from the real AUG: no second AUG that could start a competing product, no stop codon that would truncate the reporter |

---

## Metrics we add beyond the PDF

| metric | why |
|---|---|
| ΔG margin **corrected** | adds back the trigger's self-folding, which the PDF's margin row omits |
| ΔG margin **per nucleotide** | strips out the trigger-length effect that dominates the raw ranking |
| AUG unpaired **OFF**, and the ON−OFF gap | the PDF prints only the ON value, so its discriminating power is invisible |
| **P(MFE)** — how much of the ensemble the best structure accounts for | as low as 2.2%, which decides how much weight any MFE-based row deserves |
| **Centroid** structure and its probability | the structure closest to the ensemble average. For candidate 4 the centroid leaves the toehold completely unpaired, suggesting the toehold–CDS spillover the PDF worries about is an MFE artefact |
| **Ensemble ΔG** and diversity | how broad the structural ensemble is |
| **NUPACK column** for each energy | a second engine as a cross-check, not a replacement |

---

## Three things to carry into any reading of the document

1. **"RBS hidden — OFF" is accessibility, not protection.** High is exposed.
   The row cannot rank candidates, and the per-candidate verdicts resting on it
   are inverted.
2. **"Main stem pairs 50 / 67 / 75%" is arithmetic, not folding quality.** It is
   exactly (n−3)/n for n = 6/9/12, because the 3-nt bulge sits in the
   denominator. Folding each switch and counting shows **all five form 100% of
   their intended base pairs** — main stem and upper stem alike. The conclusion
   that a longer main stem folds more cleanly is not supported by the data.
3. **The ΔG margin ranking largely tracks trigger length** (r = −0.872). Read
   the per-nucleotide column with it.
