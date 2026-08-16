"""Central, fully-tunable configuration for the AND-gate design pipeline.

Every knob the specification calls "tunable" (Section 5.2) lives here, along
with the biophysical thresholds pulled from the source literature.  A single
``PipelineConfig`` instance is threaded through every stage so a design run is
reproducible from one object.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Sequence


@dataclass
class PipelineConfig:
    # ------------------------------------------------------------------ #
    # Stage 1 -- target scanning and trigger geometry                    #
    # ------------------------------------------------------------------ #
    Lx: int = 6
    """Length of the complementary core segment (|x| == |k2|).

    6, not 12.  Lx sets three things at once:
      * findability  -- 4^6 = 4k (exact matches occur between real genes)
                        vs 4^12 = 16.7M (essentially never)
      * trigger independence -- x:k2 is also the A:B duplex.  At Lx=6,
        dG ~ -6.9 (Kd ~ 15 uM) so ~100% of both mRNAs stay free at cellular
        concentrations; at Lx=12, dG ~ -19.4 (Kd ~ 0.02 pM) so ~0% stay free
        and the two transcripts would sequester each other completely.
      * the k2*:k2 helix that Trigger B invades -- 6 nt, matching Green 2026's
        6-nt invasion.
    """

    L_A: int = 36
    """Total length of Trigger A (|r1| + |x| + |a| + |k1|)."""

    L_B: int = 30
    """Total length of Trigger B (|r2| + |k2|)."""

    len_a: int = 4
    """Spacer 'a' between x and k1 on Trigger A.  ~4 nt per Kim 2019."""

    len_k1: int = 6
    """k1 = the 3' invasion domain that nucleates opening of the primary stem
    (6 nt invasion, Green 2026 Series A)."""

    len_r1: int = 8
    """r1 = upstream portion of Trigger A that seeds the secondary stem arm.
    Derived quantity: len_r1 == L_A - Lx - len_a - len_k1 unless overridden."""

    len_r2: int | None = None
    """r2 = upstream portion of Trigger B.  If ``None`` it is derived as
    L_B - Lx.  Exposed as a tunable per Section 5.2."""

    max_hamming_fraction: float = 0.34
    """Reject a G1/G2 segment pair whose minimum Hamming distance exceeds this
    fraction of Lx.  Perfect reverse complements score 0."""

    # --- stage-1 window filters (ported from the trigger-selection script,
    #     with the defaults deliberately changed) ------------------------ #
    scan_max_gc: float = 0.0
    """Maximum GC fraction of a candidate trigger window (0 == disabled).

    Off by default: a GC-rich window is a soft liability that stage 2 already
    measures as poor accessibility, so hard-rejecting it discards real biology
    before it has been weighed."""

    scan_forbid_motifs: bool = False
    """Apply ``forbidden_runs`` / ``forbidden_motifs`` to the *trigger* windows.

    The standalone scanner did this by default, which is wrong for this project:
    Type IIS sites (BsaI / BsmBI / SapI / BbsI) are Golden Gate constructability
    hazards, and the trigger is a natural sequence that is never synthesised --
    only the switch is cloned.  Leave False; the same motifs are still enforced
    on the switch in stage 4.  Set True only to reproduce old scanner runs."""

    # ------------------------------------------------------------------ #
    # Stage 2 -- thermodynamic filtering / accessibility                 #
    # ------------------------------------------------------------------ #
    flanking_lengths: Sequence[int] = (0, 10, 25, 50, 100)
    """Windows (nt each side) around a trigger binding site for MFE / SED /
    NED, matching the Toehold-VISTA feature set (+/-100 nt)."""

    accessibility_flank: int = 100
    """Flank used for the pass/fail accessibility gate in stage 2."""

    min_accessibility: float = 0.5
    """Minimum mean unpaired probability over a trigger for it to be 'open'.
    Higher == more single-stranded == more accessible."""

    max_trigger_sed: float = 0.5
    """Maximum normalised specified-ensemble-defect (vs. fully-open structure)
    allowed for a trigger to pass the accessibility gate."""

    enforce_performance_gates: bool = False
    """Whether ``min_accessibility`` / ``max_trigger_sed`` DISCARD a candidate.

    False by default, deliberately.  These are performance thresholds, not
    validity checks, and at this stage we want to see the whole distribution
    rather than throw away real biology against a number we have not calibrated.
    Accessibility is not ignored: it is a Pareto objective in stage 3, so a
    poorly accessible candidate loses on merit instead of being silently
    compensated by a bonus.  Set True to hard-filter (matching the design
    review's recommendation) once the thresholds are calibrated."""

    # -- local (RNAplfold) accessibility ------------------------------------ #
    use_local_accessibility: bool = True
    """Measure trigger accessibility with RNAplfold (ViennaRNA ``probs_window``)
    in addition to the global fold of the flanked window.

    Folding a whole mRNA globally is not trustworthy at kilobase length: the MFE
    structure of a 3 kb transcript is dominated by long-range pairs the ribosome
    never sees.  RNAplfold restricts pairing to a sliding window, which is the
    standard accessibility model (and what the NOT-gate pipeline already uses).

    The measurements are recorded but do NOT yet drive the stage-2 gate --
    changing the gate is a separate task (worst-track aggregation)."""

    local_window_size: int = 240
    """RNAplfold sliding-window width (nt).  240/200 are the values validated
    against ``RNA.pfl_fold_up`` on the 189-trigger mCherry panel."""

    local_max_bp_span: int = 200
    """Maximum base-pair span allowed inside the RNAplfold window."""

    stage2_max_candidates: int | None = 400
    """Cap on how many scanned pairs reach stage 2.  None == no cap.

    Stage 2 costs ~1.2 s per pair (five flanked folds per trigger, plus the
    RNAplfold tracks), and a real gene pair produces thousands of candidates:
    thrA x thrC yields 3680, which is over an hour before stage 3 even starts.
    The old synthetic demo produced 113 and so never exposed this.

    Candidates are ordered by (exact connector match, then lowest Hamming)
    before the cap.  That ordering is cheap but weak -- it is a triage rule, not
    a quality ranking, and among the 1331 exact matches thrA/thrC produces it is
    effectively arbitrary.  A real cheap prefilter (GC, homopolymer, connector
    duplex dG) belongs here; until then, raise the cap rather than trust that
    the dropped candidates were the bad ones."""

    accessibility_metric: str = "global"
    """Which accessibility number stage 3 ranks on: ``"global"`` or
    ``"local_worst"``.

    ``global``      the historical one -- a per-base mean from folding the
                    +/-``accessibility_flank`` window as a whole.
    ``local_worst`` the per-nucleotide worst case across every flank, measured
                    with RNAplfold (see ``use_local_accessibility``).

    Default is ``global`` so this switch changes nothing on its own; flip it to
    compare rankings.  ``local_worst`` is strictly the more pessimistic of the
    two, so absolute values are NOT comparable between the settings -- only
    rankings are.  Note ``min_accessibility`` was calibrated against ``global``
    and would reject far more under ``local_worst``; it is inert while
    ``enforce_performance_gates`` is False, which is the default."""

    binding_seed_len: int = 8
    """Length of the nucleation seed scored at each trigger end.

    Reported as the JOINT probability that all ``binding_seed_len`` nucleotides
    are simultaneously unpaired -- not the mean of per-base probabilities.  A
    toehold nucleates on a contiguous open stretch, so the joint quantity is the
    physically relevant one; the mean is optimistic because it lets an open base
    at one end compensate for a paired base at the other."""

    # ------------------------------------------------------------------ #
    # Stage 3 -- candidate selection (Pareto + diversity)                #
    # ------------------------------------------------------------------ #
    select_max_overlap: float = 0.5
    """Two shortlisted candidates from the same gene pair may not overlap by
    more than this fraction on BOTH triggers.  Without it, the shortlist fills
    with the same site shifted by a nucleotide -- the failure mode the old
    ``top_k = 8`` had.  0.5 matches the NOT-gate pipeline's own rule."""

    select_use_pareto: bool = True
    """Shortlist by Pareto front over the stage-3 criteria instead of by a
    single scalar.  Set False to fall back to the legacy ``_prescore`` scalar."""

    # ------------------------------------------------------------------ #
    # Stage 4 -- toehold-switch architecture                             #
    # ------------------------------------------------------------------ #
    primary_stem_len: int = 18
    """Primary (downstream) stem length -- Green 2026 Series A."""

    primary_invasion: int = 6
    """Toehold invasion depth into the primary stem."""

    bulge_len: int = 3
    """AUG-style bulge size in both stems."""

    secondary_loop_len: int = 11
    """Secondary (upstream/inhibitory) stem loop size.  TUNABLE (Section 5.2).
    Must NOT contain an RBS and r1 must not enter the loop."""

    rbs_seq: str = "AGAGGAGA"
    """Shine-Dalgarno / RBS presented in the primary loop."""

    ribosome_footprint_3p: int = 15
    """How far past the start codon the 30S initiation complex must find
    unstructured RNA.  The leak rate uses P(RBS -> AUG+this is simultaneously
    unpaired), NOT the accessibility of the RBS alone: in a toehold switch the
    RBS sits in the loop and is single-stranded by design, so measuring it alone
    reports a locked switch as ~80% open."""

    include_trigger_opening_cost: bool = True
    """Charge a trigger for opening its OWN structure before it can pair.

    The RNAup balance has three terms -- binding + switch-opening +
    TRIGGER-opening -- and we were scoring only the first two, which flattered
    every candidate.  GROOT's RNAup decomposition includes it as well.

    Measured on the isolated trigger, so it is a LOWER bound: in vivo the trigger
    sits inside its transcript and the true cost is larger.  Stage 2 measures
    that context separately (and does not yet feed it into the rate).
    Set False to reproduce runs made before this term existed."""

    # Conserved TSgen2 top element (RBS + start context), from Toehold-VISTA.
    hairpin_top: str = "GUUAUAGUUAUGAACAGAGGAGACAUAACAUGAAC"
    linker_suffix: str = "AACCUGGCGGCAGCGCAAAAG"
    d_domain: str = "AAC"

    # --- detuning the inhibitory hairpin (spec tunable 2) ---------------- #
    # Both act ONLY on the switch's own copy of r1 on the 5' arm; r1* is
    # Trigger A's binding site and is never touched.
    secondary_arm_gc_bias: float = 0.0
    """Fraction of the r1:r1* clamp deliberately mismatched.  0 = perfect clamp
    (strongest OFF lock); 1 = every clamp position broken (weakest).  The open
    question this knob exists to answer: Trigger B can only displace the Lx-nt
    k2*:x* helix, so the toehold it frees is capped at |a|+Lx -- unless a
    detuned r1:r1* helix melts further once B has invaded.  Raising this trades
    OFF lock for extra unmasking; sweep it, don't guess."""

    secondary_bulge_len: int = 0
    """Length of a deliberate non-pairing bulge placed immediately 3' of k2* on
    the 5' arm (0 = none).  Section 3 of the spec asked for 3."""

    # ------------------------------------------------------------------ #
    # Stage 4/5 -- relative-stability requirements                       #
    # ------------------------------------------------------------------ #
    off_state_mfe_target: float = -54.25
    """Empirical low-leak OFF-state MFE target (kcal/mol)."""

    off_state_mfe_tolerance: float = 6.0
    """OFF-state MFE within +/- this of the target incurs no leakage penalty."""

    require_secondary_stronger: bool = True
    """Enforce that the inhibitory (secondary) stem is more stable (lower MFE)
    than the primary stem, so the AND 'lock' holds."""

    # ------------------------------------------------------------------ #
    # Stage 7 -- scoring weights (hierarchical)                          #
    # ------------------------------------------------------------------ #
    weights: dict = field(default_factory=lambda: {
        "triggerB_activation": 1.0,   # 7A
        "intermediate_state": 0.8,    # 7B
        "triggerA_on_state": 1.0,     # 7C
        "penalties": 1.0,             # 7D
    })

    expression_weighting: bool = True
    """Multiply Trigger-B activation score by relative transcript abundance
    (encounter probability, spec 7A) when DE data is supplied as
    {gene_name: abundance}.  Neutral (1.0) when no data is given."""

    expression_weight_range: tuple = (0.5, 2.0)
    """Clamp on the abundance multiplier (abundance / median) so a single
    highly-expressed transcript cannot dominate the ranking."""

    # ------------------------------------------------------------------ #
    # Stage 6 -- sequence-optimisation constraints                       #
    # ------------------------------------------------------------------ #
    forbidden_runs: Sequence[str] = ("AAAA", "CCCC", "GGGG", "UUUU")
    forbid_inframe_stops: bool = True
    single_aug_after_rbs: bool = True

    # Type IIS restriction sites (+ reverse complements) that break Golden
    # Gate / MoClo assembly if they occur inside a part.  Checked on the switch
    # in RNA space.  BsaI, BsmBI/Esp3I, SapI, BbsI.
    forbidden_motifs: Sequence[str] = (
        "GGUCUC", "GAGACC",     # BsaI
        "CGUCUC", "GAGACG",     # BsmBI / Esp3I
        "GCUCUUC", "GAAGAGC",   # SapI
        "GAAGAC", "GUCUUC",     # BbsI
    )

    # --- cross-trigger crosstalk (Section 7D, ported from the trigger scanner)
    unintended_match_nt_ref: float = 8.0
    """An unintended identical / reverse-complementary stretch this long between
    Trigger A and Trigger B counts as fully disqualifying (quality 0).  The
    intended x/k2 connector is masked out before the comparison."""

    max_unintended_match: int = 0
    """Hard filter: reject when the longest unintended cross-match reaches this
    many nt (0 == disabled, score-only)."""

    crosstalk_mask_connector: bool = False
    """Whether to hide the x/k2 duplex before measuring Trigger-A/Trigger-B
    crosstalk.

    In the *split-trigger* design (the standalone scanner) the two halves are
    MEANT to hybridise through that connector, so it is masked out.  In THIS
    two-hairpin architecture the same duplex is parasitic: x = revcomp(k2) is
    required only so k2* can pair x* inside the inhibitory stem, and it makes
    the two triggers complementary to each other as an unwanted side effect.
    Default False so that liability is measured, not hidden."""

    # Off-target scan -- legacy identity scan (superseded, see below)
    offtarget_window: int | None = None      # defaults to Lx if None
    offtarget_max_identity: float = 0.85     # fraction identity that disqualifies

    # Off-target scan -- energetic, expression-weighted
    offtarget_energetic: bool = True
    """Score off-targets by duplex free energy relative to the cognate duplex
    instead of by fractional identity.

    Measured behaviour of the identity scan above, at its defaults: the window
    is ``Lx`` = 6 nt, and 85% of 6 nt rounds up to an exact 6/6 match, so it
    tests a single 6-mer -- and only the FIRST one, since it probes
    ``trigger_rc[:window]``.  For a 36-nt trigger, 30 nt are never examined.  It
    then fires on ~8% of random 400-nt transcripts, which is a function of
    transcript length rather than of binding.

    The energetic scan uses the full trigger footprint and compares free
    energies, so a GC-rich near-match and an AU-rich exact match are no longer
    scored the same."""

    offtarget_seed_k: int = 8
    """k-mer seed length for the candidate-directed background index.  A hit
    must share at least one exact k-mer with the trigger to be folded, which is
    what keeps a transcriptome-wide scan affordable."""

    offtarget_max_hits: int = 8
    """Maximum seed-ranked windows folded per trigger.  Windows are ranked by
    shared-k-mer votes, so this keeps the strongest candidates."""

    max_off_target_load: float | None = None
    """Optional hard gate on total expression-weighted off-target load.  None
    == measure and report only.  Leave None until calibrated: the load is a sum
    over the transcriptome and its scale depends on the background supplied."""

    offtarget_require_expression: bool = True
    """When an FPKM table is supplied it must cover every off-target record that
    is scored.  A missing entry raises instead of being treated as zero
    abundance -- silently scoring an unmeasured transcript as harmless is the
    failure mode this flag exists to prevent.  Set False only to reproduce
    unweighted runs."""

    # ------------------------------------------------------------------ #
    # Runtime                                                            #
    # ------------------------------------------------------------------ #
    temperature_c: float = 37.0
    sodium: float = 1.0
    magnesium: float = 0.0
    material: str = "rna"
    ensemble: str = "stacking"
    prefer_nupack: bool = True
    random_seed: int = 0
    top_n: int = 12

    # ------------------------------------------------------------------ #
    def resolved_len_r2(self) -> int:
        return self.len_r2 if self.len_r2 is not None else (self.L_B - self.Lx)

    def resolved_len_r1(self) -> int:
        derived = self.L_A - self.Lx - self.len_a - self.len_k1
        # len_r1 is an explicit knob but must stay consistent with L_A.
        return derived if derived >= 0 else self.len_r1

    def as_dict(self) -> dict:
        return asdict(self)
