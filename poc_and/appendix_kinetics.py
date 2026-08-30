"""
Stage 6 -- the kinetic model. For this architecture it is not an appendix.

WHY THIS STAGE MATTERS MORE THAN PLANNED
----------------------------------------
Stage 4 established, by inspecting the folded complex, that trigger A displaces
the inhibitory hairpin COMPLETELY on its own: ~38 bp of A:switch against ~24 bp
of hairpin, so A wins on energy alone. It has to -- a hairpin strong enough to
resist A at equilibrium could never be opened by A even with B. The consequence
is that state 10 and state 11 are indistinguishable at equilibrium, and the
start-codon separation is 0.0 points for every design.

So the AND behaviour is entirely KINETIC. Without trigger B, trigger A has only
|a| nucleotides to nucleate on; with B it has |a| + Lx. Displacement rate falls
roughly a decade per nucleotide of toehold, so that difference is the gate --
and it only counts because the transcript is degraded on a timescale the slow
route cannot beat. mRNA lifetime is a filter equilibrium cannot see.

This is the same conclusion HANDOFF.md reached from a different direction (two
metrics disagreeing ~700x on OFF leak), reproduced here independently.


THE MODEL
---------
Three-step toehold-mediated strand displacement (Zhang & Winfree 2009):

    trigger + switch <-> toehold duplex -> branch migration -> fired
                    k_on      k_off           k_bm

Steady state on the intermediate gives

    k_eff  = k_on * k_bm / (k_on * Kd_toe + k_bm)      [1/M/s]
    k_obs  = k_eff * [trigger]                         [1/s]
    P_fire = k_obs / (k_obs + k_deg),  k_deg = ln2 / half_life

with the toehold energy corrected for how buried it is:

    dG_toe = dG_duplex(trigger : site) + opening_energy(site)

`opening_energy` is the work needed to make the site single-stranded, computed
here from a constrained partition function rather than a per-base
approximation. That correction is the whole point: a toehold locked inside the
inhibitory hairpin is expensive to use even when the duplex itself is
favourable, and that expense is what state 10 pays and state 11 does not.

Reimplemented rather than imported. `and_gate_pipeline/05_kinetics.py` has the
same three equations and the same constants, but its entry points take that
pipeline's own switch and config objects, which assume VISTA Series-A geometry
(30-nt toehold, 18-nt stem, 11-nt loop) -- not ours. The physics below is
twelve lines; the geometry assumptions were the risk.


PARAMETER HONESTY -- READ BEFORE QUOTING ANY NUMBER
---------------------------------------------------
k_on and k_bm are order-of-magnitude literature values for DNA at 25 C, reused
for RNA at 37 C, and k_bm lumps a length-dependent random walk into a single
constant. **Absolute P_fire is indicative only. Never quote it as a yield.**
The RATIO between two states scored with identical constants is far more
trustworthy than either number, which is why the AND ratio is the output and
P_fire is shown only to make the ratio auditable. The sensitivity sweep at the
end exists to show how much the ranking depends on parameters we guessed.
"""

# Make relative imports work when this file is run on its own.
if __package__ in (None, ""):
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    import poc_and  # noqa: F401
    __package__ = "poc_and"

import csv
import math
import os
from dataclasses import dataclass

from .folding import RNA

from . import candidates as cd
from . import folding
from . import metrics as stage4


@dataclass
class KineticParams:
    """
    Rate constants, carried over unchanged from `and_gate_pipeline/05_kinetics.py`
    so numbers stay comparable with the earlier work.
    """

    k_on: float = 3.0e6
    """Hybridisation onto an exposed toehold [1/M/s]. Zhang & Winfree 2009
    measure ~3e6 for DNA; RNA is similar in magnitude."""

    k_bm: float = 1.0
    """Completing branch migration once the toehold is engaged [1/s]. Lumps an
    N-step random walk into one constant; ~1/s is conservative for ~20 nt."""

    mrna_half_life_s: float = 300.0
    """E. coli mRNA is typically 2-8 minutes."""

    trigger_conc_M: float = 10e-9
    """Free trigger concentration; 10 nM is a moderately expressed transcript."""

    k_ribosome: float = 0.1
    """Ribosome loading onto a fully accessible start codon [1/s]. With
    equilibrium accessibility this gives the spontaneous (trigger-independent)
    firing rate -- the floor no amount of gating can get under."""


RT_37 = 0.0019872 * 310.15          # kcal/mol


# ---------------------------------------------------------------------------
# The three equations
# ---------------------------------------------------------------------------

def displacement_rate(dG_toehold, kp):
    """k_eff [1/M/s] for a toehold of accessibility-corrected energy dG."""
    kd = math.exp(dG_toehold / RT_37)
    return kp.k_on * kp.k_bm / (kp.k_on * kd + kp.k_bm)


def fire_probability(dG_toehold, kp, conc_m=None):
    """Fraction of transcripts that fire before being degraded."""
    conc = kp.trigger_conc_M if conc_m is None else conc_m
    k_obs = displacement_rate(dG_toehold, kp) * conc
    k_deg = math.log(2.0) / kp.mrna_half_life_s
    return k_obs / (k_obs + k_deg)


def time_to_fire_s(dG_toehold, kp, conc_m=None):
    conc = kp.trigger_conc_M if conc_m is None else conc_m
    k_obs = displacement_rate(dG_toehold, kp) * conc
    return float("inf") if k_obs <= 0 else 1.0 / k_obs


def spontaneous_probability(aug_open_fraction, kp):
    """
    Leak with no trigger at all: a ribosome finding the start codon anyway.

    This is the floor. If it is comparable to the gated ON probability, the
    gate cannot be measured no matter how good its ratio looks.
    """
    k_leak = kp.k_ribosome * aug_open_fraction
    k_deg = math.log(2.0) / kp.mrna_half_life_s
    return k_leak / (k_leak + k_deg)


# ---------------------------------------------------------------------------
# Accessibility, done properly
# ---------------------------------------------------------------------------

def opening_energy(strands, span, temperature_c=37.0):
    """
    Work to hold `span` fully single-stranded, in kcal/mol (>= 0).

    Computed as the free-energy difference between the constrained ensemble
    (every base in the span forced unpaired) and the free one. This is the
    joint probability of the whole site being open, which is what a trigger
    actually needs -- not the average of per-base probabilities, which is a
    much easier bar and would flatter a site that is never open all at once.
    """
    joined = "&".join(strands) if not isinstance(strands, str) else strands
    md = RNA.md()
    md.temperature = temperature_c

    free = RNA.fold_compound(joined, md)
    _, mfe_energy = free.mfe()
    free.exp_params_rescale(mfe_energy)
    _, f_free = free.pf()

    held = RNA.fold_compound(joined, md)
    held.exp_params_rescale(mfe_energy)
    for i in range(*span):
        held.hc_add_up(i + 1)          # 1-indexed
    _, f_open = held.pf()

    return max(0.0, f_open - f_free)


WATSON_CRICK = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
WOBBLE = {("G", "U"), ("U", "G")}


def _assert_pairs(arm5, arm3, label):
    """
    Refuse to score a duplex whose domains are in the wrong order.

    Two arms that are supposed to pair must do so with zero mismatches when
    laid antiparallel. Getting a domain order backwards produces a plausible
    but meaningless energy rather than an error, which is exactly how the
    a+x / x+a mistake survived until the numbers looked odd. Wobbles are
    allowed -- they are trigger A's own, inherited from the original design.
    """
    rev = arm3[::-1]
    mismatches = sum(1 for p, q in zip(arm5, rev)
                     if (p, q) not in WATSON_CRICK and (p, q) not in WOBBLE)
    if mismatches:
        raise AssertionError(
            "%s: %d mismatch(es) -- the domains are almost certainly in the "
            "wrong order.\n  5'-%s-3'\n  3'-%s-5'"
            % (label, mismatches, arm5, rev))


# ---------------------------------------------------------------------------
# One construct
# ---------------------------------------------------------------------------

def score_kinetics(built, original_gene, kp, context_nt, temperature_c=37.0):
    """
    The four-state kinetic table for one AND switch.

    The two routes that matter:

      A alone (state 10)  nucleate on a* only, with the switch in state 00
      A given B (state 11) nucleate on a* + x*, with the switch in state 01

    Both use the accessibility the switch actually has in that state, which is
    where the gate lives: the site is the same, its availability is not.
    """
    design = built["design"]
    cand_id = built["cand"]
    spans = built["spans"]
    switch = built["sequence"]

    real = cd.trigger_real(cand_id, original_gene)
    a_lo, a_hi = real["gene_pos"], real["gene_pos"] + len(real["sequence"])
    trigger_a, _ = stage4.trigger_with_context(original_gene, a_lo, a_hi, context_nt)
    b_lo, b_hi = design["k2_span"][0], design["r2_span"][1]
    trigger_b, _ = stage4.trigger_with_context(design["variant"], b_lo, b_hi, context_nt)

    trig = cd.trigger_domains(cand_id, built["len_a"], built["Lx"])
    a_span = spans["a_star"]
    ax_span = (spans["x_star"][0], spans["a_star"][1])     # x* and a* are adjacent

    # --- the two nucleation routes trigger A could use --------------------
    #
    # Score BOTH routes in BOTH states and take the better one, because that
    # is what the molecule does -- A nucleates wherever it can go fastest.
    #
    # Scoring one route per state (a* in OFF, a*+x* with B) was wrong, and
    # wrong in a way that inverted the answer: holding 13 nt open always costs
    # more than holding 5 nt open, so the larger toehold looked worse however
    # well trigger B performed. Candidate 5 came out with ddG_toe = +3.17, the
    # gate apparently running backwards. Taking the minimum over routes fixes
    # it structurally -- state 01 always has at least the options state 00 has,
    # so the gate can never appear to run the wrong way as an artefact.
    # Domain ORDER matters and is easy to get wrong. Trigger A reads
    # 5'-k1-a-x-r1-3' and the switch site reads 5'-x*-a*-3', so the facing
    # trigger segment is "a then x" -- NOT "x then a". Coding it the other way
    # round gave 6 mismatches and -7.60 kcal/mol where the correct order gives
    # 0 mismatches and -18.80, which is an 11 kcal/mol error straight into the
    # rate. Asserted below so it cannot recur silently.
    nucleation_segment = trig["a"] + trig["x"]
    site_ax = switch[ax_span[0]:ax_span[1]]
    _assert_pairs(nucleation_segment, site_ax, "trigger a+x : x*a* site")

    dup_a = folding.duplex_dG(trig["a"], switch[a_span[0]:a_span[1]], temperature_c)
    dup_ax = folding.duplex_dG(nucleation_segment, site_ax, temperature_c)

    b_offset = len(trigger_b)
    with_b = [trigger_b, switch]

    def shifted(span):
        return (span[0] + b_offset, span[1] + b_offset)

    open_a_off = opening_energy(switch, a_span, temperature_c)
    open_ax_off = opening_energy(switch, ax_span, temperature_c)
    open_a_b = opening_energy(with_b, shifted(a_span), temperature_c)
    open_ax_b = opening_energy(with_b, shifted(ax_span), temperature_c)

    route_a_off = dup_a + open_a_off
    route_ax_off = dup_ax + open_ax_off
    route_a_b = dup_a + open_a_b
    route_ax_b = dup_ax + open_ax_b

    dG_toe_alone = min(route_a_off, route_ax_off)
    dG_toe_given_b = min(route_a_b, route_ax_b)
    route_alone = "a*" if route_a_off <= route_ax_off else "a*+x*"
    route_given_b = "a*" if route_a_b <= route_ax_b else "a*+x*"
    off_open = open_a_off if route_alone == "a*" else open_ax_off
    given_b_open = open_a_b if route_given_b == "a*" else open_ax_b

    # --- the gate itself: displacement only -------------------------------
    # This is the part the model supports. Both states use identical rate
    # constants, so the RATIO between them is meaningful even though neither
    # absolute number is.
    p_disp_10 = fire_probability(dG_toe_alone, kp)
    p_disp_11 = fire_probability(dG_toe_given_b, kp)

    # --- the trigger-independent floor: reported, NOT folded into the gate --
    #
    # Kept separate on purpose. k_ribosome applied to raw start-codon
    # accessibility assumes every unpaired AUG initiates translation, and the
    # candidates PDF explicitly warns against exactly that: an unpaired AUG
    # inside a constrained loop is an UPPER BOUND on ribosome access, not a
    # rate. With AUG ~77% accessible in OFF -- a property inherited from the
    # original single-input candidate, not caused by anything we added -- this
    # term returns a 97% spontaneous leak, which swamped every AND ratio to
    # exactly 1.0 when it was combined in.
    #
    # It cannot be right: the single-input candidates demonstrably work at the
    # bench. So k_ribosome needs calibrating against their measured ON/OFF
    # before this term means anything, and until then it is shown beside the
    # gate rather than multiplied into it.
    aug_00 = built.get("_aug_00", 0.0) / 100.0
    aug_01 = built.get("_aug_01", 0.0) / 100.0
    p_00 = spontaneous_probability(aug_00, kp)
    p_01 = spontaneous_probability(aug_01, kp)

    p_10 = p_disp_10
    p_11 = p_disp_11
    worst_off = p_disp_10
    return {
        "cand": cand_id,
        "len_a": built["len_a"],
        "Lx": built["Lx"],
        "m": built["m"],
        "bulge_size": built["bulge_size"],
        "dup_a_dG": dup_a,
        "dup_ax_dG": dup_ax,
        "opening_off_kcal": off_open,
        "opening_given_b_kcal": given_b_open,
        "route_alone": route_alone,
        "route_given_b": route_given_b,
        "open_a_off": open_a_off,
        "open_ax_off": open_ax_off,
        "open_a_b": open_a_b,
        "open_ax_b": open_ax_b,
        "dG_toe_A_alone": dG_toe_alone,
        "dG_toe_A_given_B": dG_toe_given_b,
        "ddG_toe": dG_toe_given_b - dG_toe_alone,
        "t_fire_alone_s": time_to_fire_s(dG_toe_alone, kp),
        "t_fire_given_b_s": time_to_fire_s(dG_toe_given_b, kp),
        "P_00": p_00,
        "P_10": p_10,
        "P_01": p_01,
        "P_11": p_11,
        "P_disp_10": p_disp_10,
        "P_disp_11": p_disp_11,
        "AND_ratio": (p_disp_11 / worst_off) if worst_off > 0 else float("inf"),
        "leak_floor_00": p_00,
        "leak_floor_01": p_01,
        "floor_dominates": p_00 > p_disp_11,
        "_aug_00_fraction": aug_00,
        "leak_state": max(("00", "01", "10"),
                          key=lambda s: {"00": p_00, "01": p_01, "10": p_10}[s]),
        "_built": built,
        "_trigger_a": trigger_a,
        "_trigger_b": trigger_b,
        "_a_span": a_span,
        "_ax_span": ax_span,
    }


def sensitivity(row, kp, half_lives=(60.0, 300.0, 600.0),
                concentrations=(1e-9, 10e-9, 100e-9)):
    """
    How much does the AND ratio depend on the constants we guessed?

    k_on, k_bm, the half-life and the trigger concentration are all estimates.
    A ranking that survives this sweep is worth acting on; one that does not
    should be reported as undetermined rather than dressed up.
    """
    aug_00_fraction = row.get("_aug_00_fraction", 0.0)
    ratios = []
    for half_life in half_lives:
        for conc in concentrations:
            kp2 = KineticParams(k_on=kp.k_on, k_bm=kp.k_bm,
                                mrna_half_life_s=half_life,
                                trigger_conc_M=conc,
                                k_ribosome=kp.k_ribosome)
            p_alone = fire_probability(row["dG_toe_A_alone"], kp2, conc)
            p_given = fire_probability(row["dG_toe_A_given_B"], kp2, conc)
            # displacement only, matching the headline metric
            ratios.append((p_given / p_alone) if p_alone > 0 else float("inf"))
    finite = [r for r in ratios if r != float("inf")]
    return {
        "AND_ratio_min": min(finite) if finite else float("inf"),
        "AND_ratio_max": max(finite) if finite else float("inf"),
        "AND_ratio_median": sorted(finite)[len(finite) // 2] if finite else float("inf"),
        "n_conditions": len(ratios),
    }


# ---------------------------------------------------------------------------
# Stage entry point
# ---------------------------------------------------------------------------

def run(config, results=None):
    """Stage 6 entry point."""
    print()
    print("STAGE 6 -- kinetic scoring (the stage that decides 10 vs 11)")

    stage4_rows = (results or {}).get(4)
    if not stage4_rows:
        print("  (stage 4 has not run in this session -- running it now)")
        stage4_rows = stage4.run(config, results)

    mcherry = cd.read_mcherry(config["mcherry_file"])
    original = mcherry["original"]
    kp = KineticParams(
        mrna_half_life_s=config.get("mrna_half_life_s", 300.0),
        trigger_conc_M=config["trigger_conc_M"],
    )
    context_nt = config.get("trigger_context_nt", 25)

    print("  k_on %.0e /M/s, k_bm %.1f /s, mRNA half-life %.0f s, [trigger] %.0f nM"
          % (kp.k_on, kp.k_bm, kp.mrna_half_life_s, kp.trigger_conc_M * 1e9))

    rows = []
    for r4 in stage4_rows:
        built = r4["_built"]
        built["_aug_00"] = r4["00_aug_open"]
        built["_aug_01"] = r4["01_aug_open"]
        row = score_kinetics(built, original, kp, context_nt,
                             config["temperature_C"])
        row.update(sensitivity(row, kp))
        row["nucleation_gain_nt"] = r4["nucleation_gain_nt"]
        row["x_star_release_pts"] = r4["x_star_release_pts"]
        row["sequestered_pct"] = r4["sequestered_pct"]
        rows.append(row)

    rows.sort(key=lambda r: (r["sequestered_pct"] > 10.0, -r["AND_ratio"]))
    path = write_outputs(rows, config["out_dir"])
    _report(rows, kp)
    print()
    print("wrote %s" % path)
    return rows


def _fmt(value):
    if value == float("inf"):
        return "inf"
    if value >= 1000:
        return "%.1e" % value
    return "%.1f" % value


def _report(rows, kp):
    print()
    print("=" * 110)
    print("KINETIC TOEHOLD -- what trigger A actually has to grab, per state")
    print("=" * 110)
    print("%5s %6s %4s %4s %9s %9s %8s %9s %8s %9s" %
          ("cand", "bulge", "|a|", "Lx", "route", "dG_toe", "route",
           "dG_toe", "ddG_toe", "fold"))
    print("%5s %6s %4s %4s %9s %9s %8s %9s %8s %9s" %
          ("", "", "", "", "-- A alone --", "", "-- given B --", "", "", "faster"))
    for r in rows:
        fold = math.exp(-r["ddG_toe"] / RT_37)
        print("%5d %6s %4d %4d %9s %9.2f %8s %9.2f %8.2f %9s" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            r["len_a"], r["Lx"], r["route_alone"], r["dG_toe_A_alone"],
            r["route_given_b"], r["dG_toe_A_given_B"], r["ddG_toe"],
            _fmt(fold)))
    print()
    print("open OFF / open +B = work to expose the nucleation site, from a")
    print("constrained partition function. ddG_toe is the whole gate: how much")
    print("cheaper trigger A's landing becomes once B has done its job.")

    print()
    print("=" * 110)
    print("FOUR-STATE KINETIC TABLE -- fraction firing before degradation")
    print("=" * 110)
    print("Displacement only. The spontaneous floor is shown separately below,")
    print("because as parameterised it is not calibrated and swamps everything.")
    print()
    print("%5s %6s %13s %13s %12s %11s %12s" %
          ("cand", "bulge", "10 LEAK", "11 ON", "AND ratio", "t_fire ON",
           "floor 00"))
    for r in rows:
        print("%5d %6s %13.2e %13.2e %12s %10s %12.2e %s" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            r["P_disp_10"], r["P_disp_11"], _fmt(r["AND_ratio"]),
            ("%.0fs" % r["t_fire_given_b_s"]) if r["t_fire_given_b_s"] < 1e6
            else "never",
            r["leak_floor_00"],
            "<- floor above ON" if r["floor_dominates"] else ""))

    print()
    print("=" * 110)
    print("ROBUSTNESS -- AND ratio across %d parameter combinations" % rows[0]["n_conditions"])
    print("=" * 110)
    print("mRNA half-life 60/300/600 s x [trigger] 1/10/100 nM. A ranking that")
    print("survives this is worth acting on; one that does not is undetermined.")
    print()
    print("%5s %6s %12s %12s %12s %10s %12s" %
          ("cand", "bulge", "min", "median", "max", "gain(nt)", "release"))
    for r in rows:
        print("%5d %6s %12s %12s %12s %9.1f %11.1f" % (
            r["cand"], ("%dnt" % r["bulge_size"]) if r["bulge_size"] else "none",
            _fmt(r["AND_ratio_min"]), _fmt(r["AND_ratio_median"]),
            _fmt(r["AND_ratio_max"]), r["nucleation_gain_nt"],
            r["x_star_release_pts"]))

    print()
    print("!" * 110)
    print("P_fire values are INDICATIVE ONLY -- never quote one as a yield.")
    print("k_on and k_bm are order-of-magnitude DNA constants at 25 C reused for")
    print("RNA at 37 C, and k_bm lumps a length-dependent random walk into one")
    print("number. The RATIO between two states scored with identical constants")
    print("is the trustworthy output; the probabilities are shown so the ratio")
    print("can be audited, not so they can be reported.")
    print("!" * 110)


COLUMNS = ["cand", "len_a", "Lx", "m", "bulge_size",
           "dup_a_dG", "dup_ax_dG", "opening_off_kcal", "opening_given_b_kcal",
           "route_alone", "route_given_b", "open_a_off", "open_ax_off",
           "open_a_b", "open_ax_b",
           "dG_toe_A_alone", "dG_toe_A_given_B", "ddG_toe",
           "t_fire_alone_s", "t_fire_given_b_s",
           "P_disp_10", "P_disp_11", "AND_ratio",
           "leak_floor_00", "leak_floor_01", "floor_dominates", "leak_state",
           "AND_ratio_min", "AND_ratio_median", "AND_ratio_max", "n_conditions",
           "nucleation_gain_nt", "x_star_release_pts", "sequestered_pct"]


def write_outputs(rows, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    from .find_second_trigger import _open_for_write
    fh, path = _open_for_write(os.path.join(out_dir, "stage6_kinetics.csv"),
                               newline="", encoding="utf-8")
    with fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    return path


# Pressing Run on this file alone does stage 6, using main.py's CONFIG.
if __name__ == "__main__":
    from poc_and.main import CONFIG
    run(dict(CONFIG))
