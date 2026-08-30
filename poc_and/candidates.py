"""
Stage 1 -- the five toehold-switch candidates from Toehold_Candidates29.7.pdf.

Nothing is computed in this file. It only holds sequences and the numbers the
PDF reports, so every other stage reads its inputs from one place.


Architecture of each switch (verified against the PDF sequences)
----------------------------------------------------------------

  5'- GGG - toehold - b_pre - bulge(3) - upper_stem(6) - LOOP(18)
          - upper_stem*(6) - AUG(3) - b_pre* - LINKER(21) -3'

Which parts pair with which:

  * upper_stem  pairs with  upper_stem*      (6 bp, closes the 18-nt RBS loop)
  * b_pre       pairs with  b_pre*           (the "main stem")
  * bulge(3) and AUG(3) are the two unpaired 3-nt segments, opposite each
    other. Only the second one is the real start codon.

The trigger covers  toehold + b_pre + bulge  and reads 5'->3' as

      k1 - a - x - r1          with   k1 = b_pre + bulge

so k1 sits at the trigger's 5' END. Verified on candidate 1: the switch's
main-stem segment is UGUAUG and revcomp(trigger[0:6]) = UGUACG -- a match with
one G-U wobble -- while the trigger's 3' end does not match at all.

Note on the triggers: the PDF's printed triggers pair perfectly with their
switches (0 mismatches) but contain 4-7 G-U wobbles, which come from
wobble_mutations=True in the NUPACK design script. Four of the five are also
one synonymous base away from the mCherry in mCherry.txt. Both versions are
carried through the pipeline -- see trigger_real() below.
"""

LINKER = "AACCUGGCGGCAGCGCAAAAG"   # Green 2014 linker, identical in all five
RBS = "AACAGAGGAGA"                # identical in all five, sits inside the loop
LOOP = "AGACAAGAACAGAGGAGA"        # 18 nt, identical in all five
LEADER = "GGG"                     # transcription leader
UPPER_STEM_LEN = 6                 # upper stem, closes the RBS loop
BULGE_LEN = 3                      # both the ascending bulge and the AUG


# ---------------------------------------------------------------------------
# The candidates.
#
#   switch  : full switch sequence as printed in the PDF (RNA alphabet)
#   trigger : the trigger as printed in the PDF
#   k1_len  : length of k1 = b_pre + bulge. The PDF calls this the "main stem"
#             and reports 6 / 9 / 12.
#   pdf     : the values the PDF reports, used by stage 0 to check our tooling.
# ---------------------------------------------------------------------------

CANDIDATES = {
    1: {
        "switch": ("GGGGGAGGUGAUGUUUGAUUUGAUGUUGAUGUUGUAUGCCGCGGAGACAAGAACAGAGG"
                   "AGACCGUGGAUGACAAACCUGGCGGCAGCGCAAAAG"),
        "trigger": "CGUACAACGUCAACAUCAAGUUGGACAUCACCUCC",
        "k1_len": 6,
        "family": "A",
        "pdf": {
            "similarity": 6.91,
            "dG_switch": -18.3, "dG_complex": -66.8, "dG_margin": -48.5,
            "toehold_open_off": 50.0,
            "rbs_off": 87.8, "rbs_on": 94.6, "aug_on": 98.5,
            "aug_bulge_off": "clean (0/3 paired)",
            "after_start_on": 3.4,
            "trigger_binding": 99.2, "trigger_unfolded": 67.3,
            "offtarget_dG": -23.3,
        },
    },
    2: {
        "switch": ("GGGUGUGGGAGGUGGUGUUUAGUUUGAUGUUGAUGUUGUAUGCCGCCCAGACAAGAACA"
                   "GAGGAGAGGGCGGAUGACAAACCUGGCGGCAGCGCAAAAG"),
        "trigger": "CGUACAACGUCAACAUCAAGUUGGACAUCACCUCCCACA",
        "k1_len": 6,
        "family": "A",
        "pdf": {
            "similarity": 7.62,
            "dG_switch": -22.1, "dG_complex": -76.7, "dG_margin": -54.6,
            "toehold_open_off": 59.0,
            "rbs_off": 96.5, "rbs_on": 97.1, "aug_on": 99.0,
            "aug_bulge_off": "clean (0/3 paired)",
            "after_start_on": 0.6,
            "trigger_binding": 98.0, "trigger_unfolded": 70.7,
            "offtarget_dG": -26.5,
        },
    },
    3: {
        "switch": ("GGGGGGGGUGAUGUUUAAUUUGAUGUUGAUGUUGUAUGUGCUCCCGCAGACAAGAACAG"
                   "AGGAGAGCGGGAAUGUAUACAAACCUGGCGGCAGCGCAAAAG"),
        "trigger": "GCGCGUACAACGUCAACAUCAAGUUGGACAUCACCUCC",
        "k1_len": 9,
        "family": "A",
        "pdf": {
            "similarity": 6.98,
            "dG_switch": -21.3, "dG_complex": -76.4, "dG_margin": -55.1,
            "toehold_open_off": 50.1,
            "rbs_off": 96.7, "rbs_on": 97.1, "aug_on": 95.4,
            "aug_bulge_off": "clean (0/3 paired)",
            "after_start_on": 1.7,
            "trigger_binding": 99.8, "trigger_unfolded": 65.8,
            "offtarget_dG": -25.0,
        },
    },
    4: {
        "switch": ("GGGUGUUCGUAUUGUUUCAUGAUGGUGUAGUCCUCGUUGUGGGAGCGACCAGACAAGAA"
                   "CAGAGGAGAGGUCGUAUGUAUAACGAGAACCUGGCGGCAGCGCAAAAG"),
        "trigger": "UCCCACAACGAGGACUACACCAUCGUGGAACAGUACGAACG",
        "k1_len": 12,
        "family": "A",
        "pdf": {
            "similarity": 7.46,
            "dG_switch": -27.2, "dG_complex": -87.6, "dG_margin": -60.4,
            "toehold_open_off": 52.1,
            "rbs_off": 94.7, "rbs_on": 96.6, "aug_on": 96.0,
            "aug_bulge_off": "1/3 paired (acceptable)",
            "after_start_on": 1.8,
            "trigger_binding": 98.0, "trigger_unfolded": 60.6,
            "offtarget_dG": -26.5,
        },
    },
    5: {
        "switch": ("GGGCUAUAGUCUUUUUCUGUAUUAUGGGGUCGUCGGAGGGGUCCUGCAGACAAGAACAG"
                   "AGGAGAGCAGGAAUGUUUCGAAACCUGGCGGCAGCGCAAAAG"),
        "trigger": "CCCCUCCGACGGCCCCGUAAUGCAGAAGAAGACUAUGG",
        "k1_len": 9,
        "family": "B",
        "pdf": {
            "similarity": 8.33,
            "dG_switch": -24.9, "dG_complex": -83.7, "dG_margin": -58.8,
            "toehold_open_off": 51.9,
            "rbs_off": 96.5, "rbs_on": 95.6, "aug_on": 90.6,
            "aug_bulge_off": "1/3 paired (acceptable)",
            "after_start_on": 4.6,
            "trigger_binding": 99.7, "trigger_unfolded": 69.4,
            "offtarget_dG": -24.5,
        },
    },
}


# ---------------------------------------------------------------------------
# Small sequence helpers (RNA alphabet, U not T)
# ---------------------------------------------------------------------------

_COMPLEMENT_RNA = str.maketrans("ACGU", "UGCA")
_COMPLEMENT_DNA = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    """Reverse complement of an RNA sequence, Watson-Crick only."""
    return seq.translate(_COMPLEMENT_RNA)[::-1]


def revcomp_dna(seq):
    """Reverse complement of a DNA sequence."""
    return seq.translate(_COMPLEMENT_DNA)[::-1]


def to_rna(seq):
    return seq.upper().replace("T", "U")


def to_dna(seq):
    return seq.upper().replace("U", "T")


# ---------------------------------------------------------------------------
# Where each part of a switch starts and ends.
#
# Returned as a dict of (start, end) half-open index pairs into the switch
# string, so switch[start:end] gives that part. Every stage uses this instead
# of re-deriving offsets by hand -- that is where the earlier orientation bug
# came from.
# ---------------------------------------------------------------------------

def domains(cand_id):
    """Index spans for every named part of candidate cand_id's switch."""
    c = CANDIDATES[cand_id]
    sw = c["switch"]
    trig = c["trigger"]
    k1_len = c["k1_len"]

    leader_end = len(LEADER)                    # 3
    footprint_end = leader_end + len(trig)      # end of the trigger footprint
    toehold_end = footprint_end - k1_len        # toehold stops where k1* starts

    linker_start = sw.find(LINKER)
    if linker_start < 0:
        raise ValueError("linker not found in candidate %s" % cand_id)
    rbs_start = sw.find(RBS)
    if rbs_start < 0:
        raise ValueError("RBS not found in candidate %s" % cand_id)

    mid = footprint_end                          # start of the "middle" region
    return {
        "leader": (0, leader_end),
        "toehold": (leader_end, toehold_end),
        "k1_star": (toehold_end, footprint_end),          # b_pre + ascending bulge
        "b_pre": (toehold_end, footprint_end - BULGE_LEN),
        "asc_bulge": (footprint_end - BULGE_LEN, footprint_end),
        "upper_stem": (mid, mid + UPPER_STEM_LEN),
        "loop": (mid + UPPER_STEM_LEN, mid + UPPER_STEM_LEN + len(LOOP)),
        "upper_stem_star": (mid + 24, mid + 30),
        "aug": (mid + 30, mid + 33),
        "b_pre_star": (mid + 33, linker_start),
        "linker": (linker_start, len(sw)),
        "footprint": (leader_end, footprint_end),
        "rbs": (rbs_start, rbs_start + len(RBS)),
    }


def part(cand_id, name):
    """The actual subsequence of a named part."""
    start, end = domains(cand_id)[name]
    return CANDIDATES[cand_id]["switch"][start:end]


def trigger_domains(cand_id, len_a, len_x, trigger_seq=None):
    """
    Split trigger A into its four domains.

    Trigger A reads 5'->3' as  k1 - a - x - r1.  len_a and len_x are the design
    knobs |a| and Lx; |r1| is whatever is left over.
    """
    c = CANDIDATES[cand_id]
    trig = trigger_seq if trigger_seq is not None else c["trigger"]
    k1_len = c["k1_len"]
    if k1_len + len_a + len_x > len(trig):
        raise ValueError(
            "|k1|=%d + |a|=%d + Lx=%d exceeds trigger length %d"
            % (k1_len, len_a, len_x, len(trig)))
    i = k1_len
    j = i + len_a
    k = j + len_x
    return {"k1": trig[:i], "a": trig[i:j], "x": trig[j:k], "r1": trig[k:]}


# ---------------------------------------------------------------------------
# mCherry
# ---------------------------------------------------------------------------

# Which words in a FASTA title mean which version. Checked in this order, so
# put the more specific patterns first.
_MCHERRY_LABELS = [
    ("codon_max", ("maximize", "maximiz", "max the nucle", "maximum")),
    ("pdf_reference", ("supplied", "reference", "end of the pdf", "engineered")),
    ("original", ("original", "wild", "wt")),
]


def read_mcherry(path, verbose=False):
    """
    Read the three mCherry versions from mCherry.txt.

    The file is FASTA-like: each entry is a title line starting with '>' and
    the sequence on the following line(s). Entries are told apart BY THEIR
    TITLE TEXT, not by their position in the file -- see _MCHERRY_LABELS above.
    Relying on order would silently give wrong answers if the file were ever
    reordered or a fourth sequence added, and every downstream result depends
    on getting 'original' right.

    Returns a dict with keys 'original', 'codon_max', 'pdf_reference', each an
    uppercase DNA string. The PDF reference carries a 5'UTR and a 3' tail; the
    other two are bare CDS.

    Set verbose=True to print which title was matched to which key.
    """
    import re

    with open(path, encoding="utf-8") as fh:
        raw = fh.read()

    # Split on lines beginning with '>'; the first line of each block is the
    # title and everything after it, joined, is the sequence.
    entries = []
    for block in re.split(r"^>+\s*", raw, flags=re.M):
        if not block.strip():
            continue
        lines = block.strip().split("\n")
        title = lines[0].strip()
        sequence = "".join(lines[1:])
        sequence = re.sub(r"[^A-Za-z]", "", sequence).upper().replace("U", "T")
        if sequence:
            entries.append((title, sequence))

    found = {}
    unmatched = []
    for title, sequence in entries:
        lowered = title.lower()
        for key, patterns in _MCHERRY_LABELS:
            if key in found:
                continue
            if any(p in lowered for p in patterns):
                found[key] = sequence
                if verbose:
                    print("  mCherry '%s' <- \"%s\"" % (key, title[:60]))
                break
        else:
            unmatched.append(title)

    missing = [k for k, _ in _MCHERRY_LABELS if k not in found]
    if missing:
        raise ValueError(
            "could not identify %s in %s.\n"
            "  titles found: %s\n"
            "  Titles are matched on keywords (see _MCHERRY_LABELS in "
            "candidates.py). Either rename the title in the file or add the "
            "wording used to that list."
            % (", ".join(missing), path, [t for t, _ in entries]))
    if unmatched and verbose:
        print("  (ignored unlabelled entries: %s)" % unmatched)

    return found


def trigger_real(cand_id, mcherry_original):
    """
    The trigger as the real gene actually transcribes it.

    The PDF's printed triggers for candidates 1, 2, 3 and 5 differ from the
    original mCherry by exactly one synonymous base; candidate 4 matches
    exactly. This finds the best ungapped match and returns the gene's own
    sequence there, plus where it sits and how many bases differ.
    """
    printed_dna = to_dna(CANDIDATES[cand_id]["trigger"])
    n = len(printed_dna)
    best_score = -1
    best_pos = -1
    for i in range(len(mcherry_original) - n + 1):
        window = mcherry_original[i:i + n]
        score = sum(1 for a, b in zip(window, printed_dna) if a == b)
        if score > best_score:
            best_score = score
            best_pos = i
    return {
        "sequence": to_rna(mcherry_original[best_pos:best_pos + n]),
        "gene_pos": best_pos,
        "n_diff": n - best_score,
    }
