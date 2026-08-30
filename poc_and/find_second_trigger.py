r"""
Stage 2 -- find the second trigger, and design the mCherry variant that carries it.

WHAT THIS STAGE IS FOR
----------------------
The AND gate needs a second input. Input A is the original mCherry (the gene
that already works at our bench). Input B has to be a DIFFERENT transcript, so
we can express them separately and read all four states 00 / 01 / 10 / 11.

Input B is the original mCherry with a small number of SYNONYMOUS edits -- same
protein, different nucleotides -- concentrated in two windows and leaving the
rest of the gene untouched. The codon-max variant (244 edits, 34% of the gene)
failed at the bench, so staying close to the working sequence is a hard goal.
Our variants change roughly 4% of it.


THE GEOMETRY
------------
    5'- r2* - k2* - r1copy - [loop] - r1* - x* - a* - k1* - [primary hairpin] -3'
         \______ added ______/         \____ existing toehold, untouched ____/

    TRIGGER A  5'- k1 - a - x - r1 -3'      (a stretch of the ORIGINAL mCherry)
    TRIGGER B  5'- k2 - r2 -3'              (a stretch of the VARIANT mCherry)

In OFF the inhibitory hairpin pairs `k2*` with `x*` and `r1copy` with `r1*`,
leaving only |a| nucleotides of trigger A's landing site exposed -- too short to
fire. Trigger B lands on `r2*`, invades `k2*`, releases `x*`, and trigger A's
toehold grows from |a| to |a| + Lx. Only then can A fire.


THE THREE KNOBS, AND WHY m IS NOT Lx
------------------------------------
    |a|   toehold left exposed in OFF        (Kim tested 4, 7, 10)
    Lx    length of the k2*:x* stem          (the OFF lock)
    m     how many of those Lx positions trigger B actually matches

Because `k2* = x` and `k2 = x*`, trigger A and trigger B are reverse complements
over the matched positions -- so they can bind EACH OTHER instead of the switch,
and then neither does its job. That coupling is governed by m, not Lx:

    m = 6  ->  Kd ~ 43 uM   ->   0% of trigger A tied up
    m = 7  ->  Kd ~ 386 nM  ->   2.5%
    m = 8  ->  Kd ~ 7.9 nM  ->  42%     <- crosses cellular mRNA (~10 nM)

So m is capped (config m_max) while Lx stays free for a tighter OFF lock.


TWO THINGS THAT ARE EASY TO GET WRONG, AND WERE
-----------------------------------------------
1. The Lx - m positions trigger B must NOT match are broken DELIBERATELY, and
   the break must avoid the G-U wobble partner as well as the Watson-Crick one.
   Leaving them merely unconstrained failed twice over: the gene matched them by
   accident (candidate 4's nominal m=6 design had a real A:B duplex of -17.5
   kcal/mol and 99% of trigger A tied up, where m=6 should give -6.2 and 0%),
   and simply picking a different base often produced a G-U pair, which pairs.

2. Sequestration and orthogonality are scored from the ACTUAL sequences of the
   whole genes, never from the nominal m and never from a single window. Both
   inputs are 711-nt transcripts and can interact anywhere.
"""

import csv
import os

# Make relative imports work when this file is run on its own (the Run button in
# Visual Studio executes it as a plain script, with no package context). Runs
# only in that case; a normal "import poc_and.x" skips it entirely.
if __package__ in (None, ""):
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    import poc_and  # noqa: F401  -- makes the parent package real
    __package__ = "poc_and"

from .folding import RNA   # via folding, so the setup check runs first

from . import candidates as cd
from . import codon_usage
from . import folding


# ---------------------------------------------------------------------------
# Genetic code (bacterial table 11), hard-coded so this module needs nothing
# beyond ViennaRNA.
# ---------------------------------------------------------------------------

CODON_TO_AA = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_TO_CODONS = {}
for _codon, _aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)


# Which bases can pair with a given base, in DNA spelling. G-U (here G-T) is a
# real wobble pair and still holds a helix together, so it has to be handled
# explicitly in BOTH directions: allowed when we want pairing, forbidden when we
# are trying to break it.
PAIRS_WITH_WOBBLE = {"A": "T", "C": "G", "G": "CT", "T": "AG"}
PAIRS_WATSON_CRICK = {"A": "T", "C": "G", "G": "C", "T": "A"}


def translate(dna):
    """Translate a DNA coding sequence. Trailing partial codons are ignored."""
    return "".join(CODON_TO_AA[dna[i:i + 3]]
                   for i in range(0, len(dna) - 2, 3))


def synonymous_choices(cds):
    """For each codon, every codon encoding the same amino acid."""
    return [AA_TO_CODONS[CODON_TO_AA[cds[i:i + 3]]]
            for i in range(0, len(cds) - 2, 3)]


def _best(usable, key):
    """
    Pick the best codon on `key`, breaking ties by E. coli codon usage.

    The tie-break is the point: before it, ties were settled by whatever order
    the codon happened to sit in the table, which is arbitrary. It never costs
    anything on the primary objective.
    """
    return max(usable, key=lambda c: (key(c), codon_usage.codon_fraction(c)))


# ---------------------------------------------------------------------------
# What trigger B has to read, position by position
# ---------------------------------------------------------------------------

def k2_pairing_sets(x_rna, allow_wobble):
    """
    For each position of k2, the set of bases that PAIR with the facing base of
    k2* (which is x itself).

    k2 runs antiparallel to k2*, so k2 position i faces x position Lx-1-i.
    With wobble allowed a G in x accepts C or T, and a T accepts A or G, which
    roughly doubles the number of gene positions able to host k2.
    """
    table = PAIRS_WITH_WOBBLE if allow_wobble else PAIRS_WATSON_CRICK
    x_dna = cd.to_dna(x_rna)
    n = len(x_dna)
    return [table[x_dna[n - 1 - i]] for i in range(n)]


def break_block(len_x, m, offset):
    """
    Which k2 positions are deliberately broken.

    They form a contiguous block of length Lx - m starting at `offset`, so
    offset=0 puts them at the 5' end of k2 (far from r2) and larger offsets move
    them into the middle.

    The trade-off is worth measuring rather than assuming: trigger B invades
    from the r2 side, so a break in the middle stalls invasion earlier than one
    at the far end -- but it also splits the trigger-A:trigger-B duplex into two
    short helices, which suppresses sequestration better than one truncation.
    """
    return set(range(offset, offset + (len_x - m)))


def build_constraints(x_rna, m, break_offset, allow_wobble):
    """
    Per-position requirement for k2: a string of allowed bases, or None for the
    positions we are going to break (left free during the search, broken after).
    """
    pairing = k2_pairing_sets(x_rna, allow_wobble)
    broken = break_block(len(x_rna), m, break_offset)
    return [None if i in broken else pairing[i] for i in range(len(x_rna))]


# ---------------------------------------------------------------------------
# Searching the gene
# ---------------------------------------------------------------------------

def achievable_bases(cds, choices):
    """
    Pre-filter index: which bases each gene position could hold synonymously.

    Checking this first throws out most positions in one cheap pass, so the
    exact per-codon check only runs on the few that survive. It is a
    necessary-but-not-sufficient test (a codon must satisfy all of its
    constrained offsets at once), which is exactly what a pre-filter should be.
    """
    out = [""] * len(cds)
    for ci, options in enumerate(choices):
        lo = 3 * ci
        for k in range(3):
            out[lo + k] = "".join(sorted({c[k] for c in options}))
    return out


def install_sites(cds, choices, constraints, achievable, forbid_span=None):
    """
    Every gene position where `constraints` can be satisfied synonymously.

    `constraints[i]` is a string of allowed bases, or None for "any".

    `forbid_span` keeps the site clear of trigger A's own window. Those bases
    would otherwise have two contradictory jobs -- be recoded so variant B
    cannot fire the primary hairpin, AND spell k2 -- and since k2 is the reverse
    complement of part of trigger A, putting them close together would make
    variant B fold a hairpin over its own trigger.

    Returns [(n_edits, position)], cheapest first.
    """
    hits = []
    n = len(constraints)
    for pos in range(len(cds) - n + 1):
        if forbid_span and pos < forbid_span[1] and forbid_span[0] < pos + n:
            continue

        ok = True
        for i, allowed in enumerate(constraints):
            if allowed is None:
                continue
            reachable = achievable[pos + i]
            if not any(b in reachable for b in allowed):
                ok = False
                break
        if not ok:
            continue

        total_edits = 0
        for ci in range(pos // 3, (pos + n - 1) // 3 + 1):
            lo = 3 * ci
            if ci >= len(choices):
                ok = False
                break
            required = {j - lo: constraints[j - pos]
                        for j in range(max(lo, pos), min(lo + 3, pos + n))
                        if constraints[j - pos] is not None}
            usable = [c for c in choices[ci]
                      if all(c[k] in v for k, v in required.items())]
            if not usable:
                ok = False
                break
            original = cds[lo:lo + 3]
            total_edits += min(
                sum(1 for a, b in zip(c, original) if a != b) for c in usable)
        if ok:
            hits.append((total_edits, pos))
    return sorted(hits)


def apply_constraints(cds, choices, constraints, pos):
    """
    Write the constrained bases into the gene at `pos`.

    Cheapest synonymous codon wins, ties broken by codon usage.
    Returns (new_cds, locked_positions).
    """
    out = list(cds)
    n = len(constraints)
    locked = set()
    for ci in range(pos // 3, (pos + n - 1) // 3 + 1):
        lo = 3 * ci
        required = {j - lo: constraints[j - pos]
                    for j in range(max(lo, pos), min(lo + 3, pos + n))
                    if constraints[j - pos] is not None}
        if not required:
            continue
        usable = [c for c in choices[ci]
                  if all(c[k] in v for k, v in required.items())]
        original = cds[lo:lo + 3]
        best = _best(usable,
                     lambda c: -sum(1 for a, b in zip(c, original) if a != b))
        for k in range(3):
            out[lo + k] = best[k]
        locked.update(lo + k for k in required)
    return "".join(out), locked


def break_pairing(cds, choices, pos, x_rna, broken_indices, locked):
    """
    Force the broken k2 positions NOT to pair with k2* -- wobble included.

    A base merely different from the Watson-Crick partner is not enough: on the
    first version of this code, two of candidate 5's three broken positions came
    out as G-U wobbles, which still pair. So the forbidden set is the FULL
    pairing set and the base is chosen from outside it.

    Returns (new_cds, n_broken, n_still_pairing).
    """
    pairing = k2_pairing_sets(x_rna, allow_wobble=True)
    out = list(cds)
    broken = still = 0
    for ci in sorted({(pos + i) // 3 for i in broken_indices}):
        lo = 3 * ci
        if ci >= len(choices):
            continue
        original = cds[lo:lo + 3]
        targets = {(pos + i) - lo: pairing[i] for i in broken_indices
                   if lo <= pos + i < lo + 3}
        usable = []
        for c in choices[ci]:
            good = True
            for k in range(3):
                gp = lo + k
                if (gp in locked or k not in targets) and c[k] != original[k]:
                    good = False
                    break
            if good:
                usable.append(c)
        if not usable:
            for k, forbidden in targets.items():
                still += 1 if original[k] in forbidden else 0
                broken += 0 if original[k] in forbidden else 1
            continue

        def n_unpaired(codon):
            return sum(1 for k, forbidden in targets.items()
                       if codon[k] not in forbidden)

        best = _best(usable, n_unpaired)
        for k in range(3):
            out[lo + k] = best[k]
        for k, forbidden in targets.items():
            if best[k] in forbidden:
                still += 1
            else:
                broken += 1
    return "".join(out), broken, still


def recode_region(cds, choices, start, end, locked=frozenset()):
    """
    Change as many bases as possible inside [start, end) without changing the
    protein, without touching anything outside the window, and without
    disturbing a locked base.

    This is what buys orthogonality: the more the variant's window differs from
    the original, the less able the original is to act as trigger B (or the
    variant as trigger A). Codons are independent given the constraints, so the
    most divergent synonymous codon for each is optimal -- no search needed.
    """
    out = list(cds)
    for ci in range(start // 3, (end - 1) // 3 + 1):
        lo = 3 * ci
        if ci >= len(choices):
            break
        original = cds[lo:lo + 3]
        usable = []
        for c in choices[ci]:
            good = True
            for k in range(3):
                gp = lo + k
                if (gp < start or gp >= end or gp in locked) and c[k] != original[k]:
                    good = False
                    break
            if good:
                usable.append(c)
        if not usable:
            continue

        def divergence(codon):
            return sum(1 for k in range(3)
                       if start <= lo + k < end and lo + k not in locked
                       and codon[k] != original[k])

        best = _best(usable, divergence)
        for k in range(3):
            out[lo + k] = best[k]
    return "".join(out)


def hamming(a, b, start=None, end=None):
    """
    Count differing positions. Every edit number in this stage comes from here.

    Bookkeeping the counts as passes ran drifted the moment a pass was added:
    the break pass changed bases without counting them, so candidates 4 and 5
    reported 24 and 28 edits while actually differing by 25 and 29. Deriving
    them from the sequences cannot drift.
    """
    lo = 0 if start is None else start
    hi = len(a) if end is None else end
    return sum(1 for i in range(lo, hi) if a[i] != b[i])


# ---------------------------------------------------------------------------
# One design
# ---------------------------------------------------------------------------

def build_variant(cand_id, original, choices, len_a, len_x, m, break_offset,
                  len_r2, install_pos, a_window, allow_wobble):
    """
    Produce the variant-B gene for one point of the sweep. Four passes:

      1. the m positions of k2 trigger B must match -- installed
      2. the Lx - m positions it must NOT match -- broken, wobble included
      3. r2, recoded hard: this is the orthogonality work
      4. the trigger-A window, recoded hard so variant B cannot fire the primary

    Everything outside those windows is byte-identical to the working original.
    """
    trig = cd.trigger_domains(cand_id, len_a, len_x)
    x_rna = trig["x"]
    constraints = build_constraints(x_rna, m, break_offset, allow_wobble)

    variant, locked = apply_constraints(original, choices, constraints, install_pos)

    broken_indices = break_block(len_x, m, break_offset)
    n_broken = n_still = 0
    if broken_indices:
        variant, n_broken, n_still = break_pairing(
            variant, choices, install_pos, x_rna, broken_indices, locked)
        locked = set(locked) | {install_pos + i for i in broken_indices}

    r2_start = install_pos + len_x
    r2_end = min(r2_start + len_r2, len(original))
    variant = recode_region(variant, choices, r2_start, r2_end, locked)
    variant = recode_region(variant, choices, a_window[0], a_window[1])

    return {
        "variant": variant,
        "k2_span": (install_pos, install_pos + len_x),
        "r2_span": (r2_start, r2_end),
        "a_window": a_window,
        "x": x_rna,
        "k2_required": cd.revcomp(x_rna),
        "unmatched_broken": n_broken,
        "unmatched_still_pairing": n_still,
        # every edit number derived from the sequences, never bookkept
        "k2_edits": hamming(original, variant, install_pos, install_pos + len_x),
        "r2_edits": hamming(original, variant, r2_start, r2_end),
        "trigA_window_edits": hamming(original, variant, a_window[0], a_window[1]),
        "total_edits": hamming(original, variant),
        "k2_window_nt": len_x,
        "r2_window_nt": r2_end - r2_start,
        "trigA_window_nt": a_window[1] - a_window[0],
    }


def cheap_score(design, cand_id, original, trigger_conc_m, temperature_c):
    """
    Trigger-to-trigger sequestration, from the real sequences of both triggers.

    Cheap enough to run on every design, and it is the hard gate -- if the two
    inputs stick to each other, nothing else about the design matters.
    """
    k2_start, _ = design["k2_span"]
    _, r2_end = design["r2_span"]
    trigger_b = cd.to_rna(design["variant"][k2_start:r2_end])
    trigger_a = cd.trigger_real(cand_id, original)["sequence"]

    ab = folding.duplex_dG(trigger_a, trigger_b, temperature_c)
    design["trigger_B"] = trigger_b
    design["trigger_A_B_dG"] = ab
    design["sequestered_pct"] = 100.0 * folding.bound_fraction(
        ab, trigger_conc_m, temperature_c)
    return design


def full_score(design, cand_id, original, temperature_c):
    """
    The two orthogonality controls, scanned across the WHOLE gene.

    Comparing only the corresponding window is not a control: both inputs are
    711-nt transcripts and can interact anywhere.

    Note the A-control targets the switch's toehold FOOTPRINT, not the whole
    switch. Against a whole switch any long RNA finds a -70 kcal/mol duplex
    somewhere, which says nothing; the footprint is the only place a trigger can
    productively bind.
    """
    variant = design["variant"]
    r2_start, r2_end = design["r2_span"]

    r2_rna = cd.to_rna(variant[r2_start:r2_end])
    trigB_site = cd.revcomp(r2_rna) + design["x"]

    b_on_site = folding.duplex_dG(design["trigger_B"], trigB_site, temperature_c)
    original_on_site = folding.duplex_dG(cd.to_rna(original), trigB_site, temperature_c)

    switch = cd.CANDIDATES[cand_id]["switch"]
    f_start, f_end = cd.domains(cand_id)["footprint"]
    footprint = switch[f_start:f_end]
    variant_on_foot = folding.duplex_dG(cd.to_rna(variant), footprint, temperature_c)
    original_on_foot = folding.duplex_dG(cd.to_rna(original), footprint, temperature_c)

    design.update({
        "trigB_site_r2star_k2star": trigB_site,
        "B_on_site_dG": b_on_site,
        "original_on_site_dG": original_on_site,
        "orthogonality_margin": original_on_site - b_on_site,
        "variant_on_footprint_dG": variant_on_foot,
        "original_on_footprint_dG": original_on_foot,
        "a_kill": variant_on_foot - original_on_foot,
        "usage_whole_gene": codon_usage.mean_fraction(variant),
        "usage_r2_window": codon_usage.mean_fraction(variant, r2_start, r2_end),
        "usage_r2_original": codon_usage.mean_fraction(original, r2_start, r2_end),
    })
    return design


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------

def sweep_candidate(cand_id, original, choices, achievable, config):
    """Every viable point of the (|a|, Lx, m, break offset, site) sweep."""
    real = cd.trigger_real(cand_id, original)
    a_window = (real["gene_pos"], real["gene_pos"] + len(real["sequence"]))

    m_max = config["m_max"]
    gap_max = config["Lx_minus_m_max"]
    len_r2 = config["len_r2"]
    allow_wobble = config.get("allow_wobble_in_k2", True)
    scan_offsets = config.get("scan_break_offsets", True)

    trigger = cd.CANDIDATES[cand_id]["trigger"]
    k1_len = cd.CANDIDATES[cand_id]["k1_len"]

    designs = []
    seen = set()
    for len_a in config["a_range"]:
        for len_x in config["Lx_range"]:
            if k1_len + len_a + len_x > len(trigger):
                continue
            trig = cd.trigger_domains(cand_id, len_a, len_x)
            x_rna = trig["x"]
            for m in range(min(len_x, m_max), max(0, len_x - gap_max) - 1, -1):
                if m <= 0:
                    continue
                gap = len_x - m
                # the broken block of length `gap` can start anywhere that
                # keeps it inside k2, so offsets run 0..m
                if gap and scan_offsets:
                    offsets = range(0, m + 1)
                else:
                    offsets = [0]
                for offset in offsets:
                    constraints = build_constraints(x_rna, m, offset, allow_wobble)
                    widened = (a_window[0] - (len_x + len_r2), a_window[1])
                    sites = install_sites(original, choices, constraints,
                                          achievable, forbid_span=widened)
                    if not sites:
                        continue
                    _, pos = sites[0]
                    if pos + len_x + len_r2 > len(original):
                        continue

                    d = build_variant(cand_id, original, choices, len_a, len_x,
                                      m, offset, len_r2, pos, a_window,
                                      allow_wobble)
                    key = (d["variant"], pos)
                    if key in seen:
                        continue
                    seen.add(key)
                    d.update({
                        "cand": cand_id,
                        "len_a": len_a,
                        "Lx": len_x,
                        "m": m,
                        "Lx_minus_m": gap,
                        "break_offset": offset,
                        "len_r1": len(trig["r1"]),
                        "install_pos": pos,
                        "n_install_sites": len(sites),
                        "wobble_allowed": allow_wobble,
                    })
                    designs.append(d)
    return designs


def run(config):
    """Stage 2 entry point."""
    print()
    print("STAGE 2 -- find the second trigger and design the variant-B gene")

    mcherry = cd.read_mcherry(config["mcherry_file"])
    original = mcherry["original"]
    choices = synonymous_choices(original)
    achievable = achievable_bases(original, choices)
    reference_protein = translate(original)
    top_n = config.get("stage2_thermo_top_n", 40)

    print("sweep: |a| in %s, Lx in %s, m <= %d, Lx-m <= %d, |r2| = %d"
          % (list(config["a_range"]), list(config["Lx_range"]),
             config["m_max"], config["Lx_minus_m_max"], config["len_r2"]))
    print("       G-U wobble allowed in k2 install: %s | break offsets scanned: %s"
          % (config.get("allow_wobble_in_k2", True),
             config.get("scan_break_offsets", True)))

    all_designs = []
    for cand_id in config["candidates"]:
        designs = sweep_candidate(cand_id, original, choices, achievable, config)

        for d in designs:
            if translate(d["variant"]) != reference_protein:
                raise SystemExit(
                    "cand %d (a=%d Lx=%d m=%d off=%d): recoded variant changes "
                    "the protein -- that is a bug, not a design choice."
                    % (cand_id, d["len_a"], d["Lx"], d["m"], d["break_offset"]))
            cheap_score(d, cand_id, original,
                        config["trigger_conc_M"], config["temperature_C"])

        # cheap gate first: if the triggers sequester each other nothing else
        # matters, and the whole-gene controls below are the expensive part
        designs.sort(key=lambda d: (d["sequestered_pct"], -d["Lx"],
                                    d["total_edits"]))
        for d in designs[:top_n]:
            full_score(d, cand_id, original, config["temperature_C"])
        scored = [d for d in designs[:top_n] if "a_kill" in d]
        scored.sort(key=lambda d: (
            d["sequestered_pct"] > 10.0,
            -d["orthogonality_margin"],
            -d["a_kill"],
            # prefer designs where every intended break actually took. Some
            # cannot: Met and Trp have a single codon each, so a position
            # landing there is unbreakable. Sequestration already measures the
            # consequence, so this only settles ties.
            d["unmatched_still_pairing"],
            -d["Lx"],
            d["total_edits"],
        ))
        all_designs.extend(scored)
        _report_candidate(cand_id, designs, scored)

    csv_path, fasta_path = write_outputs(all_designs, config["out_dir"])
    print()
    print("wrote %s" % csv_path)
    print("wrote %s" % fasta_path)
    _report_summary(all_designs)
    return all_designs


COLUMNS = ["cand", "len_a", "Lx", "m", "Lx_minus_m", "break_offset", "len_r1",
           "install_pos", "n_install_sites", "wobble_allowed",
           "unmatched_broken", "unmatched_still_pairing",
           "k2_edits", "k2_window_nt", "r2_edits", "r2_window_nt",
           "trigA_window_edits", "trigA_window_nt", "total_edits",
           "trigger_A_B_dG", "sequestered_pct",
           "B_on_site_dG", "original_on_site_dG", "orthogonality_margin",
           "variant_on_footprint_dG", "original_on_footprint_dG", "a_kill",
           "usage_whole_gene", "usage_r2_original", "usage_r2_window",
           "x", "k2_required", "trigger_B", "trigB_site_r2star_k2star"]


def _open_for_write(path, **kwargs):
    """
    Open a file for writing, stepping aside if it is locked.

    Excel takes an exclusive lock on an open .csv, so leaving the previous
    results open in a window would otherwise throw away the whole run at the
    last step. Falls back to path_1, path_2, ... and says so.
    """
    try:
        return open(path, "w", **kwargs), path
    except PermissionError:
        stem, ext = os.path.splitext(path)
        for n in range(1, 100):
            alt = "%s_%d%s" % (stem, n, ext)
            try:
                handle = open(alt, "w", **kwargs)
            except PermissionError:
                continue
            print("  NOTE: %s is locked (open in Excel?) -- wrote %s instead"
                  % (os.path.basename(path), os.path.basename(alt)))
            return handle, alt
        raise


def write_outputs(all_designs, out_dir):
    """One CSV of every scored design, plus a FASTA of the best per candidate."""
    os.makedirs(out_dir, exist_ok=True)
    fh, csv_path = _open_for_write(
        os.path.join(out_dir, "stage2_second_trigger.csv"),
        newline="", encoding="utf-8")
    with fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for d in all_designs:
            writer.writerow(d)

    best = {}
    for d in all_designs:
        best.setdefault(d["cand"], d)
    fh, fasta_path = _open_for_write(
        os.path.join(out_dir, "stage2_variantB_genes.fasta"), encoding="utf-8")
    with fh:
        for cand_id in sorted(best):
            d = best[cand_id]
            fh.write(">mCherry_variantB_cand%d  a=%d Lx=%d m=%d break_offset=%d "
                     "edits=%d\n" % (cand_id, d["len_a"], d["Lx"], d["m"],
                                     d["break_offset"], d["total_edits"]))
            for i in range(0, len(d["variant"]), 60):
                fh.write(d["variant"][i:i + 60] + "\n")
    return csv_path, fasta_path


def _report_candidate(cand_id, designs, scored):
    print()
    print("-" * 104)
    print("Candidate %d -- %d designs swept, %d scored in full"
          % (cand_id, len(designs), len(scored)))
    if not scored:
        print("  none viable.")
        return
    print("%4s %4s %3s %5s %5s %7s %8s %10s %9s %8s" %
          ("|a|", "Lx", "m", "brk", "|r1|", "edits", "seq'd%", "B on site",
           "ortho", "a_kill"))
    for d in scored[:6]:
        print("%4d %4d %3d %5d %5d %7d %7.1f%% %10.1f %9.1f %8.1f" % (
            d["len_a"], d["Lx"], d["m"], d["break_offset"], d["len_r1"],
            d["total_edits"], d["sequestered_pct"], d["B_on_site_dG"],
            d["orthogonality_margin"], d["a_kill"]))
    b = scored[0]
    print("  best: |a|=%d Lx=%d m=%d break@%d, k2 at gene %d" %
          (b["len_a"], b["Lx"], b["m"], b["break_offset"], b["install_pos"]))
    print("        edits %d of %d nt (%.1f%%): k2 %d/%d, r2 %d/%d, trigA window %d/%d"
          % (b["total_edits"], len(b["variant"]),
             100.0 * b["total_edits"] / len(b["variant"]),
             b["k2_edits"], b["k2_window_nt"], b["r2_edits"], b["r2_window_nt"],
             b["trigA_window_edits"], b["trigA_window_nt"]))
    print("        broken positions: %d truly unpaired, %d still pairing"
          % (b["unmatched_broken"], b["unmatched_still_pairing"]))


def _report_summary(all_designs):
    print()
    print("=" * 104)
    print("SUMMARY -- best design per candidate")
    print("=" * 104)
    print("%5s %4s %4s %3s %5s %7s %8s %9s %8s %11s %8s" %
          ("cand", "|a|", "Lx", "m", "brk", "edits", "seq'd%", "ortho",
           "a_kill", "gene chg", "usage"))
    seen = set()
    for d in all_designs:
        if d["cand"] in seen:
            continue
        seen.add(d["cand"])
        print("%5d %4d %4d %3d %5d %7d %7.1f%% %9.1f %8.1f %10.1f%% %8.3f" % (
            d["cand"], d["len_a"], d["Lx"], d["m"], d["break_offset"],
            d["total_edits"], d["sequestered_pct"], d["orthogonality_margin"],
            d["a_kill"], 100.0 * d["total_edits"] / len(d["variant"]),
            d["usage_whole_gene"]))
    print()
    print("ortho  = dG(original : B-site) - dG(trigger B : B-site), whole-gene scan.")
    print("         Positive means trigger B out-competes the original mCherry")
    print("         for its own landing site, so state 10 reads OFF.")
    print("a_kill = dG(variant : toehold footprint) - dG(original : footprint).")
    print("         Positive means the recoded variant grips the primary toehold")
    print("         WORSE than the original, so state 01 reads OFF.")
    print("usage  = mean E. coli codon fraction of the variant (original 0.384).")


# Pressing Run on this file alone does stage 2, using main.py's CONFIG.
if __name__ == "__main__":
    from poc_and.main import CONFIG
    run(dict(CONFIG))
