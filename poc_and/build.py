"""
Stage 3 -- build the AND switch by prepending the inhibitory hairpin.

WHAT GETS BUILT
---------------
    5'- GGG - r2* - k2* - r1copy - loop - [ the ENTIRE original switch, from
                                            position 3 onward, byte-identical ]

Everything from `r1*` rightward -- the toehold, the primary stem, the RBS loop,
the AUG bulge, the linker -- is the sequence straight out of
Toehold_Candidates29.7.pdf. We only add sequence in front of it. That is what
makes these constructs a controlled comparison against the single-input
candidates already tested in the lab, and it is why the primary hairpin is never
rebuilt here.

Of the four added domains, only ONE is free sequence:

    r2*      = revcomp(r2)         fixed by the variant gene stage 2 designed
    k2*      = x                   fixed by trigger A
    r1copy   = trigger A's own r1  fixed by trigger A (see below)
    loop     free -- ours to choose, and it must not look like an RBS


WHY r1copy IS TRIGGER A'S OWN r1, NOT THE PERFECT COMPLEMENT OF r1*
-------------------------------------------------------------------
`r1*` carries G-U wobbles, because the original NUPACK design ran with
wobble_mutations enabled. A perfectly Watson-Crick `r1copy` would therefore grip
`r1*` HARDER than trigger A itself does -- measured at -22.10 vs -21.80 kcal/mol,
so trigger A would be climbing 0.30 kcal/mol uphill to displace the very clamp
that is supposed to release it.

Copying trigger A's own r1, wobbles included, makes the exchange isoenergetic:
the clamp and the trigger form the identical helix, so displacement is driven
purely by the toehold, which is textbook strand displacement. The effect is
small but it is free.


THE BULGE
---------
The inhibitory stem is Lx + |r1| base pairs -- 22 to 30 bp of continuous duplex.
E. coli RNase III cleaves double-stranded RNA and generally wants ~20 bp or more,
so a stem that long is a plausible substrate in a way Kim's ~11 bp masked region
never was. A single unpaired base inserted into `r1copy` splits it into two
shorter helices.

It is an INSERTION into the added strand, so `r1*`, the toehold and the main
trigger are all untouched -- which is the constraint that matters.

The bulge also biases displacement in trigger A's favour, which was not the
original motive: trigger A's own r1 has no bulge, so when A displaces r1copy it
forms the unbulged helix and the exchange becomes downhill rather than neutral.
That is bought with exactly the same amount of OFF lock, so it is a knife-edge of
the same kind as |a| and gets swept rather than assumed.
"""

# Make relative imports work when this file is run on its own (the Run button in
# Visual Studio executes it as a plain script, with no package context).
if __package__ in (None, ""):
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    import poc_and  # noqa: F401
    __package__ = "poc_and"

import csv
import os

from .folding import RNA

from . import candidates as cd
from . import folding


# A loop for the inhibitory hairpin. Free sequence, but it sits in the 5'UTR
# upstream of the real RBS, so the one thing it must not do is look like a
# second ribosome binding site. Checked by _utr_checks below rather than
# assumed.
LOOP_SCAFFOLD = "GAAACAGAACGAAUCAGACUUCGGAUCAG"

# Shine-Dalgarno-like motifs we refuse to create in the added 5' region.
SD_MOTIFS = ("AGGAGG", "AAGGAG", "AGGAGA", "GGAGGU", "AAGGAGG")

WATSON_CRICK = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "U"), ("U", "G")}


def fit(scaffold, n):
    """Trim or tile a scaffold to exactly n nucleotides."""
    if n <= len(scaffold):
        return scaffold[:n]
    return (scaffold * (n // len(scaffold) + 1))[:n]


# ---------------------------------------------------------------------------
# Where to put the bulge
# ---------------------------------------------------------------------------

def bulge_positions(len_x, len_r1, r1copy, r1_star, max_helix=15):
    """
    Legal insertion points inside r1copy, best first.

    Three rules, in the order they matter:

      1. Neither resulting helix may exceed `max_helix` bp, or the split has not
         actually solved the RNase III problem. Inserting at r1copy index j
         splits the stem into (Lx + j) and (|r1| - j).
      2. Stay at least 2 bp clear of the k2*:x* junction. A bulge there frays
         the very helix trigger B has to invade, and the OFF lock depends on it.
      3. Prefer a spot whose neighbouring pairs are Watson-Crick. Stacking a
         bulge on top of an existing G-U wobble concentrates two
         destabilisations in one place and the helix opens there, which costs
         more OFF lock than the RNase III benefit is worth.
    """
    lo = max(2, len_r1 - max_helix)
    hi = min(len_r1 - 1, max_helix - len_x)
    out = []
    for j in range(lo, hi + 1):
        # r1copy[j] pairs with r1_star[len_r1 - 1 - j]
        wobble_neighbours = 0
        for k in (j - 1, j):
            if 0 <= k < len_r1:
                partner = r1_star[len_r1 - 1 - k]
                if (r1copy[k], partner) in WOBBLE:
                    wobble_neighbours += 1
        balance = abs((len_x + j) - (len_r1 - j))
        out.append((wobble_neighbours, balance, j))
    out.sort()
    return [j for _, _, j in out]


# ---------------------------------------------------------------------------
# Building one switch
# ---------------------------------------------------------------------------

def build_and_switch(cand_id, design, loop_len=11, bulge_size=0,
                     bulge_index=None, max_helix=15):
    """
    Assemble the AND switch for one stage-2 design.

    Returns the sequence, the span of every domain, and the OFF structure the
    design INTENDS -- written from the spans, not folded, so that comparing it
    against the actual fold is a real test rather than a tautology.
    """
    switch = cd.CANDIDATES[cand_id]["switch"]
    spans = cd.domains(cand_id)
    len_a, len_x = design["len_a"], design["Lx"]

    trig = cd.trigger_domains(cand_id, len_a, len_x)
    r1_trigger = trig["r1"]                       # trigger A's own r1, wobbles kept
    len_r1 = len(r1_trigger)

    # r1* is the 5'-most part of the existing toehold
    toehold_start, toehold_end = spans["toehold"]
    r1_star = switch[toehold_start:toehold_start + len_r1]
    x_star = switch[toehold_start + len_r1:toehold_start + len_r1 + len_x]

    r2_star = cd.revcomp(cd.to_rna(design["variant"][
        design["r2_span"][0]:design["r2_span"][1]]))
    k2_star = trig["x"]

    # optional bulge, inserted into the ADDED strand only
    r1copy = r1_trigger
    used_index = None
    if bulge_size:
        options = bulge_positions(len_x, len_r1, r1_trigger, r1_star, max_helix)
        if bulge_index is not None and bulge_index in options:
            used_index = bulge_index
        elif options:
            used_index = options[0]
        if used_index is not None:
            r1copy = (r1_trigger[:used_index] + "A" * bulge_size
                      + r1_trigger[used_index:])

    loop = fit(LOOP_SCAFFOLD, loop_len)
    leader = cd.LEADER

    added = r2_star + k2_star + r1copy + loop
    sequence = leader + added + switch[len(leader):]

    # --- domain spans in the finished construct ---------------------------
    pos = len(leader)
    layout = []
    for name, seq in (("r2_star", r2_star), ("k2_star", k2_star),
                      ("r1copy", r1copy), ("loop", loop)):
        layout.append((name, pos, pos + len(seq)))
        pos += len(seq)
    shift = len(added)
    domain_spans = {name: (s, e) for name, s, e in layout}
    domain_spans["leader"] = (0, len(leader))
    for name, (s, e) in spans.items():
        if name == "leader":
            continue
        domain_spans["primary_" + name] = (s + shift, e + shift)
    # the three parts of the old toehold, now individually meaningful
    domain_spans["r1_star"] = (toehold_start + shift,
                               toehold_start + len_r1 + shift)
    domain_spans["x_star"] = (toehold_start + len_r1 + shift,
                              toehold_start + len_r1 + len_x + shift)
    domain_spans["a_star"] = (toehold_start + len_r1 + len_x + shift,
                              toehold_end + shift)

    # The GGG cap belongs at the 5' END of the finished construct, not stranded
    # in the middle. The PDF added it to each switch as a transcription cap, so
    # prepending our hairpin must MOVE it, not leave a copy sitting upstream of
    # r1* where it would be a random GGG inside the molecule. Asserted rather
    # than assumed, because it is invisible in the output either way.
    assert sequence.startswith(leader), "the cap is not at the 5' end"
    assert sequence[len(leader) + len(added):] == switch[len(leader):], (
        "everything from r1* onward must be the PDF switch byte-for-byte")

    structure = intended_structure(cand_id, design, sequence, domain_spans,
                                   len_r1, len_x, bulge_size, used_index)

    return {
        "cand": cand_id,
        "sequence": sequence,
        "intended_structure": structure,
        "spans": domain_spans,
        "r2_star": r2_star,
        "k2_star": k2_star,
        "r1copy": r1copy,
        "r1_star": r1_star,
        "x_star": x_star,
        "loop": loop,
        "bulge_size": bulge_size,
        "bulge_index": used_index,
        "added_nt": len(added),
        "len_a": len_a,
        "Lx": len_x,
        "m": design["m"],
        "len_r1": len_r1,
        "stem_bp": len_x + len_r1,
        "design": design,
    }


def intended_structure(cand_id, design, sequence, spans, len_r1, len_x,
                       bulge_size, bulge_index):
    """
    The OFF-state dot-bracket the design intends, assembled from the spans.

    Written independently of any folding, so that diffing it against the real
    MFE fold in stage 4 actually tests whether the molecule does what we asked.
    """
    n = len(sequence)
    out = ["."] * n

    def mark(span, char):
        for i in range(*span):
            out[i] = char

    # inhibitory hairpin: k2* pairs x*, r1copy pairs r1*
    mark(spans["k2_star"], "(")
    mark(spans["x_star"], ")")
    r1c_start, r1c_end = spans["r1copy"]
    mark((r1c_start, r1c_end), "(")
    if bulge_size and bulge_index is not None:
        for i in range(r1c_start + bulge_index,
                       r1c_start + bulge_index + bulge_size):
            out[i] = "."                       # the bulge is unpaired
    mark(spans["r1_star"], ")")

    # primary hairpin, exactly as the original candidate intends it
    mark(spans["primary_b_pre"], "(")
    mark(spans["primary_asc_bulge"], ".")
    mark(spans["primary_upper_stem"], "(")
    mark(spans["primary_loop"], ".")
    mark(spans["primary_upper_stem_star"], ")")
    mark(spans["primary_aug"], ".")
    mark(spans["primary_b_pre_star"], ")")
    mark(spans["primary_linker"], ".")

    return "".join(out)


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_helices(built):
    """
    Assert every helix the construct depends on can actually form.

    This is the guard against silent orientation errors -- the class of mistake
    that already cost us one full rebuild of the search. Arms are compared
    antiparallel: the first base of the 5' arm faces the LAST base of the 3' arm.
    """
    problems = []

    def helix(name, arm5, arm3, allow_wobble=True):
        rev = arm3[::-1]
        wc = wob = mismatch = 0
        for a, b in zip(arm5, rev):
            if (a, b) in WATSON_CRICK:
                wc += 1
            elif (a, b) in WOBBLE:
                wob += 1
            else:
                mismatch += 1
        if len(arm5) != len(arm3):
            problems.append("%s: arms differ in length (%d vs %d)"
                            % (name, len(arm5), len(arm3)))
        if mismatch:
            problems.append("%s: %d mismatch(es)" % (name, mismatch))
        if wob and not allow_wobble:
            problems.append("%s: %d wobble(s) where none allowed" % (name, wob))
        return {"wc": wc, "wobble": wob, "mismatch": mismatch}

    # k2*:x* is the OFF lock. Mismatches are forbidden, but WOBBLES ARE NOT A
    # DEFECT here -- they are trigger A's own, inherited by copying its x.
    #
    # This is the same trade-off as r1copy, with the same answer. Measured on
    # candidate 5: a perfectly Watson-Crick k2* binds x* at -6.50 where trigger
    # A's own x binds it at -5.80, so a "flawless" lock would make trigger A
    # climb 0.70 kcal/mol uphill to displace it. Copying trigger A's x makes
    # the exchange isoenergetic, and since trigger A re-forms the identical
    # helix the RELATIVE lock is not weakened at all.
    #
    # So the invariant worth asserting is not "no wobbles" but "k2* IS trigger
    # A's x". An earlier version of this check demanded zero wobbles and
    # reported a correct design as broken.
    lock = helix("k2*:x*", built["k2_star"], built["x_star"], allow_wobble=True)
    expected_k2_star = cd.trigger_domains(
        built["cand"], built["len_a"], built["Lx"])["x"]
    if built["k2_star"] != expected_k2_star:
        problems.append(
            "k2* is %s but should be trigger A's own x (%s), or trigger A has "
            "to displace a helix stronger than the one it forms"
            % (built["k2_star"], expected_k2_star))

    # r1copy:r1* legitimately carries trigger A's wobbles; that is the point.
    #
    # The bulge is removed before comparing the arms, because a bulge is NOT a
    # mismatch and lumping them together would be wrong in both directions. A
    # mismatch is two bases opposite each other that fail to pair; a bulge is an
    # insertion with nothing opposite it at all, and the arms pair perfectly
    # around it. So "mismatch 0" here is a statement about the arms, and says
    # nothing about whether a bulge is present -- which is why the bulge is
    # reported separately below rather than left implicit.
    r1copy = built["r1copy"]
    bulge_size = built["bulge_size"]
    bulge_index = built["bulge_index"]
    if bulge_size and bulge_index is not None:
        r1copy = r1copy[:bulge_index] + r1copy[bulge_index + bulge_size:]
    clamp = helix("r1copy:r1*", r1copy, built["r1_star"])

    if bulge_size and bulge_index is not None:
        # the bulge splits one long duplex into two, which is its whole purpose
        upper = built["Lx"] + bulge_index
        lower = built["len_r1"] - bulge_index
        bulge = {"size": bulge_size, "index": bulge_index,
                 "helix_upper_bp": upper, "helix_lower_bp": lower,
                 "longest_bp": max(upper, lower)}
    else:
        bulge = {"size": 0, "index": None,
                 "helix_upper_bp": built["stem_bp"], "helix_lower_bp": 0,
                 "longest_bp": built["stem_bp"]}

    return {"k2_x_lock": lock, "r1_clamp": clamp, "bulge": bulge,
            "problems": problems}


def utr_checks(built):
    """
    The only sequence rules that apply to what we added.

    It sits in the 5'UTR, upstream of the real RBS, so stop codons are
    irrelevant -- nothing is being translated there yet. What DOES matter is
    that we have not accidentally created a second ribosome binding site, and
    that the real RBS still appears exactly once.
    """
    added = built["sequence"][:built["added_nt"] + len(cd.LEADER)]
    notes = []
    for motif in SD_MOTIFS:
        if motif in added:
            notes.append("added region contains an SD-like motif: %s" % motif)
    if cd.RBS in added:
        notes.append("added region contains the RBS sequence itself")
    if built["sequence"].count(cd.RBS) != 1:
        notes.append("RBS appears %d times, expected once"
                     % built["sequence"].count(cd.RBS))
    return notes


def fold_report(built, temperature_c=37.0):
    """Fold the finished construct and compare against what we intended."""
    sequence = built["sequence"]
    ensemble = folding.ensemble(sequence, temperature_c)
    intended = built["intended_structure"]
    actual = ensemble["mfe_structure"]

    agree = sum(1 for a, b in zip(intended, actual) if a == b)
    unpaired = folding.unpaired_probs(sequence, temperature_c)
    spans = built["spans"]
    return {
        "mfe_dG": ensemble["mfe_dG"],
        "ensemble_dG": ensemble["ensemble_dG"],
        "p_mfe": ensemble["p_mfe"],
        "mfe_structure": actual,
        "centroid_structure": ensemble["centroid_structure"],
        "intended_agreement_pct": 100.0 * agree / len(intended),
        "a_star_open_off": folding.mean_unpaired(unpaired, spans["a_star"]),
        "x_star_open_off": folding.mean_unpaired(unpaired, spans["x_star"]),
        "r2_star_open_off": folding.mean_unpaired(unpaired, spans["r2_star"]),
        "rbs_open_off": folding.mean_unpaired(unpaired, spans["primary_rbs"]),
        "aug_open_off": folding.mean_unpaired(unpaired, spans["primary_aug"]),
    }


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------

def run(config, results=None):
    """
    Build one AND switch per candidate, sweeping the bulge.

    Takes stage 2's designs from `results` if stage 2 ran in this invocation,
    and re-runs the search otherwise so this file works standalone.
    """
    print()
    print("STAGE 3 -- build the AND switches")

    designs = (results or {}).get(2)
    if not designs:
        print("  (stage 2 has not run in this session -- running it now)")
        from . import find_second_trigger
        designs = find_second_trigger.run(config)

    # Carry the top K stage-2 designs per candidate, not just the winner.
    #
    # Stage 2 ranks on orthogonality, a-kill and sequestration. It cannot see
    # ddG_toe -- the kinetic gate that stage 6 showed is what actually decides
    # whether a design works -- because that needs the built construct. Passing
    # only stage 2's winner would lock in the (|a|, Lx, m) choice before the
    # deciding criterion is ever computed, from a field of ~40 scored designs.
    # Carrying several lets the kinetics choose among real alternatives.
    per_candidate = max(1, config.get("stage3_designs_per_candidate", 3))
    shortlist = {}
    for d in designs:
        shortlist.setdefault(d["cand"], []).append(d)
    for cand_id in shortlist:
        shortlist[cand_id] = shortlist[cand_id][:per_candidate]

    loop_len = config.get("secondary_loop_len", 11)
    built_all = []
    print("  carrying the top %d stage-2 design(s) per candidate" % per_candidate)
    for cand_id in sorted(shortlist):
      for design in shortlist[cand_id]:
        for bulge in config.get("bulge_sizes", [0]):
            built = build_and_switch(cand_id, design, loop_len=loop_len,
                                     bulge_size=bulge)
            if bulge and built["bulge_index"] is None:
                continue          # no legal position for a bulge this size
            built["checks"] = check_helices(built)
            built["utr_notes"] = utr_checks(built)
            built["fold"] = fold_report(built, config["temperature_C"])
            built_all.append(built)
        _report_candidate(cand_id, [b for b in built_all if b["cand"] == cand_id])

    path = write_outputs(built_all, config["out_dir"])
    print()
    print("wrote %s" % path)
    _report_summary(built_all)
    return built_all


def _report_candidate(cand_id, builds):
    if not builds:
        return
    b0 = builds[0]
    print()
    print("-" * 100)
    print("Candidate %d -- |a|=%d Lx=%d m=%d, stem %d bp, added %d nt, total %d nt"
          % (cand_id, b0["len_a"], b0["Lx"], b0["m"], b0["stem_bp"],
             b0["added_nt"], len(b0["sequence"])))
    lock = b0["checks"]["k2_x_lock"]
    clamp = b0["checks"]["r1_clamp"]
    print("  k2*:x* lock  %d WC, %d wobble, %d mismatch   (0 mismatches required;"
          " wobbles are trigger A's own)"
          % (lock["wc"], lock["wobble"], lock["mismatch"]))
    print("  r1copy:r1*   %d WC, %d wobble, %d mismatch   (arms only -- the bulge"
          " is counted separately)"
          % (clamp["wc"], clamp["wobble"], clamp["mismatch"]))
    if b0["checks"]["problems"]:
        for p in b0["checks"]["problems"]:
            print("  PROBLEM: %s" % p)
    if b0["utr_notes"]:
        for note in b0["utr_notes"]:
            print("  5'UTR note: %s" % note)
    print()
    print("  Bulge SWEEP -- all sizes are built; none is chosen yet. A bulge is an")
    print("  unpaired insertion in r1copy (not a mismatch), and it splits the one")
    print("  long stem into two shorter helices, which is the RNase III point.")
    print()
    print("  %8s %9s %9s %9s %8s %10s %10s %10s" %
          ("bulge", "helices", "longest", "MFE dG", "P(MFE)", "a* open",
           "x* open", "r2* open"))
    for b in builds:
        f = b["fold"]
        bl = b["checks"]["bulge"]
        if bl["size"]:
            tag = "%d nt @%d" % (bl["size"], bl["index"])
            helices = "%d + %d" % (bl["helix_upper_bp"], bl["helix_lower_bp"])
        else:
            tag = "none"
            helices = "%d" % bl["helix_upper_bp"]
        print("  %8s %9s %8dbp %9.1f %7.2f%% %9.1f%% %9.1f%% %9.1f%%" % (
            tag, helices, bl["longest_bp"], f["mfe_dG"], 100 * f["p_mfe"],
            f["a_star_open_off"], f["x_star_open_off"], f["r2_star_open_off"]))
    print("  longest = the longest continuous duplex; RNase III wants ~20bp+,")
    print("  so getting this under ~15 is the reason the bulge exists.")


def _report_summary(built_all):
    print()
    print("=" * 100)
    print("SUMMARY -- the OFF state we care about")
    print("=" * 100)
    print("a* is the toehold trigger A can reach BEFORE trigger B arrives; it")
    print("should be open. x* is the part the inhibitory stem hides; it should")
    print("be SHUT. r2* is trigger B's landing site; it should be open.")
    print()
    print("Showing the NO-BULGE build for each candidate. The bulged variants are")
    print("in the CSV and in the per-candidate sweeps above; which to use is a")
    print("stage-4 decision, not settled here.")
    print()
    print("%5s %7s %9s %10s %10s %10s %10s %9s" %
          ("cand", "bulge", "stem bp", "a* open", "x* shut", "r2* open",
           "RBS open", "P(MFE)"))
    for b in built_all:
        if b["bulge_size"]:
            continue
        f = b["fold"]
        print("%5d %7s %9d %9.1f%% %9.1f%% %9.1f%% %9.1f%% %8.2f%%" % (
            b["cand"], "none", b["stem_bp"], f["a_star_open_off"],
            100 - f["x_star_open_off"], f["r2_star_open_off"],
            f["rbs_open_off"], 100 * f["p_mfe"]))


COLUMNS = ["cand", "len_a", "Lx", "m", "len_r1", "stem_bp", "bulge_size",
           "bulge_index", "added_nt", "total_nt", "mfe_dG", "ensemble_dG",
           "p_mfe", "intended_agreement_pct", "a_star_open_off",
           "x_star_open_off", "r2_star_open_off", "rbs_open_off",
           "aug_open_off", "lock_mismatch", "clamp_wobble",
           "helix_upper_bp", "helix_lower_bp", "longest_bp", "utr_notes",
           "sequence", "intended_structure", "mfe_structure"]


def write_outputs(built_all, out_dir):
    """One row per built construct, with its sequence and both structures."""
    os.makedirs(out_dir, exist_ok=True)
    from .find_second_trigger import _open_for_write
    fh, path = _open_for_write(os.path.join(out_dir, "stage3_and_switches.csv"),
                               newline="", encoding="utf-8")
    with fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for b in built_all:
            row = dict(b)
            row.update(b["fold"])
            row["total_nt"] = len(b["sequence"])
            row["lock_mismatch"] = b["checks"]["k2_x_lock"]["mismatch"]
            row["clamp_wobble"] = b["checks"]["r1_clamp"]["wobble"]
            row.update(b["checks"]["bulge"])
            row["utr_notes"] = "; ".join(b["utr_notes"])
            writer.writerow(row)
    return path


# Pressing Run on this file alone does stage 3, using main.py's CONFIG.
if __name__ == "__main__":
    from poc_and.main import CONFIG
    run(dict(CONFIG))
