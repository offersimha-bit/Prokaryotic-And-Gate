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
    Lx: int = 12
    """Length of the complementary core segment (|x| == |k2|)."""

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

    # ------------------------------------------------------------------ #
    # Stage 3 -- toehold-switch architecture                             #
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

    # Conserved TSgen2 top element (RBS + start context), from Toehold-VISTA.
    hairpin_top: str = "GUUAUAGUUAUGAACAGAGGAGACAUAACAUGAAC"
    linker_suffix: str = "AACCUGGCGGCAGCGCAAAAG"
    d_domain: str = "AAC"

    # Strength balance between the two upper arms of the secondary stem.
    # +1 favours a stronger r1*/(r1) clamp, -1 favours a weaker one. TUNABLE.
    secondary_arm_gc_bias: float = 0.0

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
    (encounter probability) when DE data is supplied."""

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

    # Off-target scan
    offtarget_window: int | None = None      # defaults to Lx if None
    offtarget_max_identity: float = 0.85     # fraction identity that disqualifies

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
