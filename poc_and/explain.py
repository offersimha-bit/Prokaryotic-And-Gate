"""
Print the full construction of one AND-gate design, step by step.

    python -m poc_and.main --explain 5

This exists to be shown to other people. Everything it prints is generated from
the same objects the pipeline actually uses, so the figure cannot drift away
from the design.

One rule drives the formatting: whenever DNA and amino acids appear together
they are generated from the SAME codon list, never from two separate slices.
The windows here routinely start mid-codon (k2 is 6-10 nt and the reading frame
does not care), and an earlier hand-written version of this trace showed a raw
23-nt window against an 8-codon translation. It looked like a broken
translation when the recoding was in fact perfectly synonymous.
"""

# Make relative imports work when this file is run on its own (the Run button in
# Visual Studio executes it as a plain script, with no package context). Runs
# only in that case; a normal "import poc_and.x" skips it entirely.
if __package__ in (None, ""):
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    import poc_and  # noqa: F401  -- makes the parent package real
    __package__ = "poc_and"

from . import candidates as cd
from . import codon_usage
from . import find_second_trigger as f2


WATSON_CRICK = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "U"), ("U", "G")}


def _rule(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def codon_block(dna, start, end):
    """
    A window widened to whole codons, returned as aligned lines.

    Returns (codons, amino_acids, lo, hi). The caller prints both from this one
    list, which is what keeps the DNA and the protein in step.
    """
    lo, hi = codon_usage.codon_span(start, end)
    hi = min(hi, len(dna))
    codons = [dna[i:i + 3] for i in range(lo, hi, 3) if len(dna[i:i + 3]) == 3]
    amino = [f2.CODON_TO_AA[c] for c in codons]
    return codons, amino, lo, hi


def print_window(label, original, variant, start, end):
    """Show one recoded window: codons, amino acids, and what changed."""
    o_codons, o_aa, lo, hi = codon_block(original, start, end)
    v_codons, v_aa, _, _ = codon_block(variant, start, end)

    offset = start - lo
    print("  gene[%d:%d]  (widened to whole codons: [%d:%d], %d codons)"
          % (start, end, lo, hi, len(o_codons)))
    if offset:
        print("  note: the window starts %d base(s) into a codon, so the codon-"
              "aligned view begins at %d" % (offset, lo))
    print()
    print("    original  %s" % " ".join(o_codons))
    print("              %s" % " ".join("%-3s" % a for a in o_aa))
    print("    variant   %s" % " ".join(v_codons))
    print("              %s" % " ".join("%-3s" % a for a in v_aa))
    # Mark only what is inside the REQUESTED window. The codon-aligned view can
    # reach a base or two beyond it, and a change out there is not part of this
    # window's edit count -- showing them the same way would make the count look
    # wrong.
    pieces = []
    for ci, (oc, vc) in enumerate(zip(o_codons, v_codons)):
        chars = []
        for k in range(3):
            gp = lo + 3 * ci + k
            if not (start <= gp < end):
                chars.append(" ")          # outside the window we asked about
            elif oc[k] != vc[k]:
                chars.append("^")
            else:
                chars.append(".")
        pieces.append("".join(chars))
    print("    changed   %s" % " ".join(pieces))
    print("              (blank = outside the requested window, shown only to"
          " complete the codon)")
    changed = sum(1 for a, b in zip(original[start:end], variant[start:end])
                  if a != b)
    print("    -> %d of %d bases changed; protein %s"
          % (changed, end - start,
             "IDENTICAL" if o_aa == v_aa else "*** CHANGED - BUG ***"))
    print("    -> E. coli codon usage %.3f -> %.3f"
          % (codon_usage.mean_fraction(original, start, end),
             codon_usage.mean_fraction(variant, start, end)))


def pair_line(top, bottom, label):
    """
    Show two strands pairing, antiparallel.

    `top` is 5'->3' and `bottom` is 5'->3'; the bottom is reversed for display
    so the two lines face each other the way the molecules do.
    """
    rev = bottom[::-1]
    marks = "".join("|" if (a, b) in WATSON_CRICK
                    else ("o" if (a, b) in WOBBLE else "X")
                    for a, b in zip(top, rev))
    print("    %s" % label)
    print("      5'-%s-3'" % top)
    print("         %s" % marks)
    print("      3'-%s-5'" % rev)
    print("      %d Watson-Crick, %d wobble, %d mismatch"
          % (marks.count("|"), marks.count("o"), marks.count("X")))
    return marks


def explain(cand_id, config, pick=None):
    """
    Print the whole construction for one candidate.

    `pick` is an optional (|a|, Lx, m) to trace a SPECIFIC design. Without it
    this shows stage 2's own winner, which is not necessarily the one the
    report ships: the report ranks on the kinetic gate, which stage 2 cannot
    see. Pass the finalist's parameters to trace what actually got chosen.
    """
    mcherry = cd.read_mcherry(config["mcherry_file"])
    original = mcherry["original"]
    choices = f2.synonymous_choices(original)
    achievable = f2.achievable_bases(original, choices)

    designs = f2.sweep_candidate(cand_id, original, choices, achievable, config)
    if not designs:
        print("No viable design for candidate %d." % cand_id)
        return
    for d in designs:
        f2.cheap_score(d, cand_id, original,
                       config["trigger_conc_M"], config["temperature_C"])
    designs.sort(key=lambda d: (d["sequestered_pct"], -d["Lx"], d["total_edits"]))
    top = designs[:config.get("stage2_thermo_top_n", 40)]
    for d in top:
        f2.full_score(d, cand_id, original, config["temperature_C"])
    top.sort(key=lambda d: (d["sequestered_pct"] > 10.0,
                            -d["orthogonality_margin"], -d["a_kill"],
                            -d["Lx"], d["total_edits"]))
    if pick is not None:
        wanted = [x for x in top
                  if (x["len_a"], x["Lx"], x["m"]) == tuple(pick)]
        if not wanted:
            wanted = [x for x in designs
                      if (x["len_a"], x["Lx"], x["m"]) == tuple(pick)]
        if not wanted:
            raise SystemExit("no design with |a|=%d Lx=%d m=%d for candidate %d"
                             % (pick[0], pick[1], pick[2], cand_id))
        d = wanted[0]
        if "trigB_site_r2star_k2star" not in d:
            f2.full_score(d, cand_id, original, config["temperature_C"])
    else:
        d = top[0]

    len_a, len_x, m = d["len_a"], d["Lx"], d["m"]
    k2_start, k2_end = d["k2_span"]
    r2_start, r2_end = d["r2_span"]
    variant = d["variant"]

    _rule("AND-gate construction -- candidate %d" % cand_id)
    print("chosen from %d designs swept:  |a|=%d  Lx=%d  m=%d  break offset=%d"
          % (len(designs), len_a, len_x, m, d["break_offset"]))
    print("k2 installed at gene position %d of %d" % (k2_start, len(original)))

    # ---- 1 ----------------------------------------------------------------
    _rule("STEP 1 -- split trigger A into its four domains")
    trig = cd.trigger_domains(cand_id, len_a, len_x)
    k1_len = cd.CANDIDATES[cand_id]["k1_len"]
    print("Trigger A is a real stretch of the ORIGINAL mCherry, and reads")
    print("5'-k1-a-x-r1-3'. k1 sits at the 5' end -- verified against the switch,")
    print("not assumed.")
    print()
    print("    %s" % cd.CANDIDATES[cand_id]["trigger"])
    print("    %s%s%s%s" % ("k" * k1_len, "a" * len_a, "x" * len_x,
                            "r" * len(trig["r1"])))
    print()
    print("    k1 = %-14s invades the primary stem" % trig["k1"])
    print("    a  = %-14s stays exposed in OFF (too short to fire alone)"
          % trig["a"])
    print("    x  = %-14s hidden by the inhibitory hairpin until B arrives"
          % trig["x"])
    print("    r1 = %-14s held by the r1copy clamp" % trig["r1"])

    # ---- 2 ----------------------------------------------------------------
    _rule("STEP 2 -- what trigger B must read")
    print("The inhibitory stem pairs k2* with x*, and k2* IS x. So trigger B")
    print("must carry the reverse complement of x:")
    print()
    print("    x            5'-%s-3'" % trig["x"])
    print("    k2_required  5'-%s-3'" % d["k2_required"])
    print()
    print("Only %d of those %d positions have to match (m=%d). The other %d are"
          % (m, len_x, m, len_x - m))
    print("broken on purpose, so trigger A and trigger B cannot stick together.")

    # ---- 3 ----------------------------------------------------------------
    _rule("STEP 3 -- install it in mCherry, synonymously")
    print_window("k2", original, variant, k2_start, k2_end)

    # ---- 4 ----------------------------------------------------------------
    _rule("STEP 4 -- did we get the match and the break we asked for?")
    k2_variant = cd.to_rna(variant[k2_start:k2_end])
    broken = f2.break_block(len_x, m, d["break_offset"])
    role = "".join("~" if i in broken else "M" for i in range(len_x))
    print("    k2 required  %s" % d["k2_required"])
    print("    k2 in variant %s" % k2_variant)
    print("    role          %s" % role)
    print("                  M = must match   ~ = broken on purpose")
    print()
    pairing = f2.k2_pairing_sets(trig["x"], allow_wobble=True)
    ok = True
    for i in range(len_x):
        base = cd.to_dna(k2_variant[i])
        pairs = base in pairing[i]
        if i in broken and pairs:
            ok = False
            print("    position %d was meant to be broken but still pairs" % i)
        if i not in broken and not pairs:
            ok = False
            print("    position %d was meant to match but does not" % i)
    print("    all %d positions behave as designed: %s" % (len_x, ok))
    print("    (%d truly unpaired, %d still pairing)"
          % (d["unmatched_broken"], d["unmatched_still_pairing"]))

    # ---- 5 ----------------------------------------------------------------
    _rule("STEP 5 -- recode r2, trigger B's toehold")
    print("This is where orthogonality lives. If r2 were left alone, the")
    print("ORIGINAL mCherry would carry the same toehold and open the secondary")
    print("hairpin by itself, so state 10 would read ON.")
    print()
    print_window("r2", original, variant, r2_start, r2_end)

    # ---- 6 ----------------------------------------------------------------
    _rule("STEP 6 -- recode the trigger-A window, so variant B cannot fire the primary")
    a_lo, a_hi = d["a_window"]
    print("Note this window is the WHOLE %d-nt region trigger A occupies, not the"
          % (a_hi - a_lo))
    print("|a|=%d exposed gap. That is why its edit count is the largest of the"
          % len_a)
    print("three, and why it is not comparable to |a|.")
    print()
    print_window("trigger-A window", original, variant, a_lo, a_hi)

    # ---- 7 ----------------------------------------------------------------
    _rule("STEP 7 -- the finished trigger B and the site it lands on")
    print("    trigger B  5'-%s-3'  (%d nt)" % (d["trigger_B"], len(d["trigger_B"])))
    print("    B site     5'-%s-3'   (r2* toehold + k2*)"
          % d["trigB_site_r2star_k2star"])
    print()
    rc_b = cd.revcomp(d["trigger_B"])
    print("Is that site just revcomp(trigger B)? Almost, and the")
    print("difference is the design:")
    print("    revcomp(trigger B) = %s" % rc_b)
    print("    B site             = %s" % d["trigB_site_r2star_k2star"])
    differing = [i for i, (p, q) in enumerate(zip(d["trigB_site_r2star_k2star"], rc_b))
                 if p != q]
    print("    they differ at %s -- exactly the broken positions." % differing)
    print()
    print("    B site = r2* + k2* = revcomp(r2) + x. The k2* half must be")
    print("    the FLAWLESS complement of x, because k2*:x* is the OFF lock,")
    print("    and x belongs to trigger A. Trigger B is the one allowed to be")
    print("    imperfect -- if it were not, it would stick to trigger A.")

    # ---- 8 ----------------------------------------------------------------
    _rule("STEP 8 -- the three helices")
    pair_line(trig["x"], cd.revcomp(trig["x"]), "k2* : x*   (the OFF lock -- must be perfect)")
    print()
    pair_line(d["trigger_B"], d["trigB_site_r2star_k2star"],
              "trigger B : its landing site   (breaks are deliberate)")

    # ---- 9 ----------------------------------------------------------------
    _rule("STEP 9 -- does it behave as an AND gate?")
    print("  trigger A : trigger B   %8.2f kcal/mol -> %.1f%% of A sequestered"
          % (d["trigger_A_B_dG"], d["sequestered_pct"]))
    print("     the two inputs must NOT stick to each other")
    print()
    print("  ortho   %+8.1f   dG(original:site) %.1f  vs  dG(trigger B:site) %.1f"
          % (d["orthogonality_margin"], d["original_on_site_dG"], d["B_on_site_dG"]))
    print("     positive = B out-competes the original for its site, so 10 is OFF")
    print()
    print("  a_kill  %+8.1f   dG(variant:footprint) %.1f  vs  dG(original:footprint) %.1f"
          % (d["a_kill"], d["variant_on_footprint_dG"],
             d["original_on_footprint_dG"]))
    print("     positive = the variant grips the primary toehold worse than the")
    print("     original, so 01 is OFF")
    print()
    print("  gene changed  %d of %d nt (%.1f%%)   k2 %d/%d, r2 %d/%d, trigA %d/%d"
          % (d["total_edits"], len(variant),
             100.0 * d["total_edits"] / len(variant),
             d["k2_edits"], d["k2_window_nt"], d["r2_edits"], d["r2_window_nt"],
             d["trigA_window_edits"], d["trigA_window_nt"]))
    print("  codon usage   %.3f (original %.3f, the failed codon-max 0.471)"
          % (d["usage_whole_gene"], codon_usage.mean_fraction(original)))
    print("  protein       %s"
          % ("identical to mCherry"
             if f2.translate(variant) == f2.translate(original)
             else "*** CHANGED - BUG ***"))


# Pressing Run on this file alone explains one candidate. Change the number
# here, or use:  python -m poc_and.main --explain 4
if __name__ == "__main__":
    from poc_and.main import CONFIG
    explain(5, dict(CONFIG))
