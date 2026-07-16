"""Stage 7 -- multi-stage, hierarchical scoring and ranking.

Each design is scored through four blocks that mirror the specification:

7A  Trigger-B activation      -- can B open the inhibitory stem?
7B  Post-activation state     -- is the Trigger-A site exposed & stable after B?
7C  Trigger-A activation / ON -- ON-state MFE, RBS liberation, translation
7D  Global penalties          -- leakage, off-target, restricted seqs, half-life

Sub-scores are normalised to ~[0, 1] (higher == better); penalties are
subtracted.  Trigger-binding steps use real cofold thermodynamics and
constraint-conditioned accessibility rather than assumed base-pairing.
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field

from . import sequence_utils as su
from .architecture import DesignedSwitch
from .config import PipelineConfig
from .filtering import TriggerMetrics
from .offtarget import scan_offtargets
from .optimize import check_restricted
from .thermo import ThermoBackend


# --------------------------------------------------------------------------- #
# small normalisation helpers                                                 #
# --------------------------------------------------------------------------- #
def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def _neg_to_01(dg: float, best: float, worst: float) -> float:
    """Map a free energy so that ``best`` (most negative) -> 1 and ``worst`` -> 0."""
    if best == worst:
        return 0.0
    return _clamp01((worst - dg) / (worst - best))


def _load_codon_fractions(path: str | None) -> dict:
    if path and os.path.exists(path):
        out = {}
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                codon = su.to_rna(row["Codon"])
                try:
                    out[codon] = float(row["Fraction"])
                except (KeyError, ValueError):
                    continue
        if out:
            return out
    return {}


# --------------------------------------------------------------------------- #
@dataclass
class ScoreCard:
    triggerB_activation: float = 0.0
    intermediate_state: float = 0.0
    triggerA_on_state: float = 0.0
    penalty: float = 0.0
    total: float = 0.0
    quality_percent: float = 0.0
    """Absolute 0-100 quality of the three positive criteria, measured against
    FIXED reference scales -- never min-maxed against the other candidates in
    this run.  A 70 today means the same as a 70 after changing the gene set or
    top-N, so it is comparable across runs.  ``total`` is this quality (as a
    0-1 weighted sum) minus penalties, and is what ranking uses."""

    details: dict = field(default_factory=dict)
    flags: list[str] = field(default_factory=list)

    def breakdown(self, weights: dict) -> list[tuple]:
        """Per-criterion (label, quality, weight, points, max_points) on the
        0-100 scale -- mirrors the standalone scanner's score breakdown."""
        pos = [
            ("(A) Trigger-B activation", self.triggerB_activation,
             weights["triggerB_activation"]),
            ("(B) intermediate: A site opens", self.intermediate_state,
             weights["intermediate_state"]),
            ("(C) Trigger-A / ON state", self.triggerA_on_state,
             weights["triggerA_on_state"]),
        ]
        wsum = sum(w for _l, _q, w in pos) or 1.0
        return [(label, q, w, 100.0 * w * q / wsum, 100.0 * w / wsum)
                for label, q, w in pos]

    def as_row(self) -> dict:
        row = {
            "score_total": round(self.total, 4),
            "quality_percent": round(self.quality_percent, 2),
            "score_triggerB": round(self.triggerB_activation, 4),
            "score_intermediate": round(self.intermediate_state, 4),
            "score_triggerA_on": round(self.triggerA_on_state, 4),
            "penalty": round(self.penalty, 4),
            "flags": ";".join(self.flags),
        }
        for k, v in self.details.items():
            row[k] = round(v, 4) if isinstance(v, float) else v
        return row


class DesignScorer:
    def __init__(self, cfg: PipelineConfig, backend: ThermoBackend,
                 codon_table_path: str | None = None,
                 transcriptome: dict | None = None,
                 essential_genes: set | None = None,
                 expression: dict | None = None):
        self.cfg = cfg
        self.backend = backend
        self.codon = _load_codon_fractions(codon_table_path)
        self.transcriptome = transcriptome or {}
        self.essential = essential_genes or set()
        self.expression = expression or {}

    # ------- 7A ---------------------------------------------------------- #
    def _score_triggerB(self, sw: DesignedSwitch, tmB: TriggerMetrics,
                        details: dict) -> float:
        b = self.backend
        core = sw.core
        k2_idx = sw.triggerB_toehold_indices()
        toehold_access = b.region_accessibility(core, k2_idx)

        s_k2, e_x = sw.domains.spans["sec_k2star"][0], sw.domains.spans["sec_xstar"][1]
        secondary_region = core[s_k2:e_x]
        dg_bind = b.binding_dG(sw.triggerB, secondary_region)

        target_access = tmB.accessibility
        bind_norm = _neg_to_01(dg_bind, best=-30.0, worst=0.0)

        encounter = 1.0
        if self.cfg.expression_weighting and self.expression:
            g = sw.pair.triggerB.gene
            # expression keyed by gene name isn't available here; use a supplied
            # relative-abundance multiplier if present, else neutral.
            encounter = float(self.expression.get("TriggerB", 1.0))

        details.update({
            "B_target_access": target_access,
            "B_toehold_access": toehold_access,
            "B_bind_dG": dg_bind,
            "B_encounter": encounter,
        })
        base = (0.4 * target_access + 0.3 * toehold_access + 0.3 * bind_norm)
        return _clamp01(base) * encounter

    # ------- 7B ---------------------------------------------------------- #
    def _score_intermediate(self, sw: DesignedSwitch, details: dict) -> float:
        b = self.backend
        core = sw.core
        a_idx = sw.triggerA_footprint_indices()
        # Trigger B has seized its toehold + r1 arm -> hold those unpaired
        forced = sw.triggerB_toehold_indices() + sw.domains.region("sec_r1")
        off_access = b.region_accessibility(core, a_idx)
        int_access = b.region_accessibility(core, a_idx, forced_unpaired=forced)
        delta = int_access - off_access

        int_mfe = b.complex_mfe([sw.triggerB, core])
        stability = _neg_to_01(int_mfe, best=self.cfg.off_state_mfe_target, worst=0.0)

        details.update({
            "int_A_access_off": off_access,
            "int_A_access_afterB": int_access,
            "int_A_access_gain": delta,
            "int_complex_mfe": int_mfe,
        })
        # reward exposure of the A site and a stable (non-trap) intermediate
        return _clamp01(0.6 * int_access + 0.25 * _clamp01(delta * 3)
                        + 0.15 * stability)

    # ------- 7C ---------------------------------------------------------- #
    def _score_on_state(self, sw: DesignedSwitch, details: dict) -> float:
        b = self.backend
        core = sw.core
        on_mfe = b.complex_mfe([sw.triggerA, sw.triggerB, core])
        rbs_idx = sw.primary_rbs_indices()
        forced = (sw.triggerA_footprint_indices()
                  + sw.triggerB_toehold_indices()
                  + sw.domains.region("sec_r1"))
        rbs_liberation = b.region_accessibility(core, rbs_idx, forced_unpaired=forced)
        transl = self._translation_efficiency(sw)

        on_norm = _neg_to_01(on_mfe, best=-120.0, worst=-20.0)
        details.update({
            "on_complex_mfe": on_mfe,
            "rbs_liberation": rbs_liberation,
            "translation_eff": transl,
        })
        return _clamp01(0.45 * on_norm + 0.35 * rbs_liberation + 0.20 * transl)

    def _translation_efficiency(self, sw: DesignedSwitch, n_codons: int = 10) -> float:
        if not self.codon:
            return 0.5
        rbs = su.to_rna(self.cfg.rbs_seq)
        idx = sw.full.find(rbs)
        augs = [i for i in su.find_all(sw.full, "AUG") if i > idx]
        if not augs:
            return 0.0
        start = augs[0]
        coding = sw.full[start:]
        fracs = []
        for c in range(3, 3 + 3 * n_codons, 3):  # skip AUG itself
            codon = coding[c:c + 3]
            if len(codon) < 3:
                break
            fracs.append(self.codon.get(codon, 0.0))
        return _clamp01(sum(fracs) / len(fracs)) if fracs else 0.5

    # ------- 7D ---------------------------------------------------------- #
    def _penalties(self, sw: DesignedSwitch, tmA: TriggerMetrics,
                   tmB: TriggerMetrics, details: dict, flags: list) -> float:
        b = self.backend
        cfg = self.cfg
        core = sw.core
        penalty = 0.0

        # OFF-state MFE / leakage --------------------------------------- #
        off_mfe = b.mfe(core)[1]
        details["off_state_mfe"] = off_mfe
        leak_gap = off_mfe - (cfg.off_state_mfe_target + cfg.off_state_mfe_tolerance)
        if leak_gap > 0:                      # less negative than allowed -> leaky
            penalty += 0.05 * leak_gap
            flags.append(f"leaky_off_mfe(+{leak_gap:.1f})")

        # relative stability: secondary must be stronger than primary ---- #
        s0 = sw.domains.spans["sec_k2star"][0]
        s1 = sw.domains.spans["sec_xstar"][1]
        p0 = sw.domains.spans["prim_k1star"][0]
        sec_mfe = b.mfe(core[s0:s1])[1]
        prim_mfe = b.mfe(core[p0:])[1]
        details["secondary_mfe"] = sec_mfe
        details["primary_mfe"] = prim_mfe
        if cfg.require_secondary_stronger and sec_mfe >= prim_mfe:
            penalty += 0.25
            flags.append("secondary_not_stronger")

        # spacer 'a' binding strength ----------------------------------- #
        a_seq = sw.pair.triggerA.a
        astar = sw.domain_seq("spacer_astar")
        a_dg = b.binding_dG(a_seq, astar) if a_seq and astar else 0.0
        details["spacer_a_dG"] = a_dg
        if a_dg > -3.0:
            penalty += 0.10
            flags.append("weak_spacer_a")

        # cross-trigger crosstalk (ported from the standalone scanner) ---- #
        # Do the two triggers stick to (reverse-complement) or mimic
        # (identity) each other *beyond* the intended x/k2 connector?  The
        # connector is masked so it is not flagged as a problem.
        ta, tb = sw.pair.triggerA, sw.pair.triggerB
        tA_masked = su.mask_region(ta.seq, len(ta.r1), len(ta.x))
        tB_masked = su.mask_region(tb.seq, len(tb.r2), len(tb.k2))
        stick_nt = su.max_revcomp_match(tA_masked, tB_masked)
        subst_nt = su.max_identity_match(tA_masked, tB_masked)
        ref = max(1.0, cfg.unintended_match_nt_ref)
        details["crosstalk_stick_nt"] = stick_nt
        details["crosstalk_subst_nt"] = subst_nt
        details["crosstalk_stick_quality"] = max(0.0, 1.0 - stick_nt / ref)
        details["crosstalk_subst_quality"] = max(0.0, 1.0 - subst_nt / ref)
        penalty += 0.20 * min(1.0, stick_nt / ref)
        penalty += 0.10 * min(1.0, subst_nt / ref)
        if cfg.max_unintended_match and max(stick_nt, subst_nt) >= cfg.max_unintended_match:
            penalty += 0.4
            flags.append(f"crosstalk({max(stick_nt, subst_nt)}nt)")

        # restricted sequences ------------------------------------------ #
        rep = check_restricted(sw, cfg)
        details["forbidden_runs"] = len(rep.forbidden_runs)
        details["aug_after_rbs"] = rep.aug_after_rbs
        if rep.forbidden_runs:
            penalty += 0.05 * len(rep.forbidden_runs)
            flags.append("forbidden_run")

        # Type IIS restriction sites (Golden Gate / MoClo hazard) -------- #
        motifs = su.has_forbidden_run(sw.full, cfg.forbidden_motifs)
        details["type2s_sites"] = len(motifs)
        if motifs:
            penalty += 0.10 * len(motifs)
            flags.append(f"type2s({','.join(motifs)})")
        if rep.inframe_stop:
            penalty += 0.30
            flags.append("inframe_stop")
        if rep.aug_after_rbs != 1:
            penalty += 0.20
            flags.append(f"aug_after_rbs={rep.aug_after_rbs}")

        # trigger accessibility gate (Stage 2 re-checked as a penalty) --- #
        if not tmA.passes:
            penalty += 0.15
            flags.append("triggerA_inaccessible")
        if not tmB.passes:
            penalty += 0.15
            flags.append("triggerB_inaccessible")

        # off-target risk ------------------------------------------------ #
        if self.transcriptome:
            hitsA = scan_offtargets(sw.triggerA, self.transcriptome, cfg,
                                    self.essential, exclude={sw.pair.triggerA.gene})
            hitsB = scan_offtargets(sw.triggerB, self.transcriptome, cfg,
                                    self.essential, exclude={sw.pair.triggerB.gene})
            n_ess = sum(1 for h in hitsA + hitsB if h.essential)
            n_hit = len(hitsA) + len(hitsB)
            details["offtarget_hits"] = n_hit
            details["offtarget_essential"] = n_ess
            if n_ess:
                penalty += 0.5 * n_ess
                flags.append(f"offtarget_essential({n_ess})")
            elif n_hit:
                penalty += 0.05 * n_hit
                flags.append(f"offtarget({n_hit})")

        # temporal stability / half-life proxy -------------------------- #
        mean_paired = 1.0 - (sum(b.unpaired_probabilities(core)) / len(core))
        details["mean_paired_frac"] = mean_paired
        if mean_paired < 0.35:               # too unstructured -> short-lived
            penalty += 0.10
            flags.append("low_structure_halflife")

        return penalty

    # ------- driver ------------------------------------------------------ #
    def score(self, sw: DesignedSwitch, tmA: TriggerMetrics,
              tmB: TriggerMetrics) -> ScoreCard:
        details: dict = {}
        flags: list[str] = []
        w = self.cfg.weights

        sB = self._score_triggerB(sw, tmB, details)
        sInt = self._score_intermediate(sw, details)
        sOn = self._score_on_state(sw, details)
        pen = self._penalties(sw, tmA, tmB, details, flags)

        positive = (w["triggerB_activation"] * sB
                    + w["intermediate_state"] * sInt
                    + w["triggerA_on_state"] * sOn)
        total = positive - w["penalties"] * pen

        # Absolute 0-100 quality: every sub-score is already measured against a
        # fixed reference scale, so dividing by the weight sum keeps it
        # comparable across runs / gene sets / top-N choices.
        wsum = (w["triggerB_activation"] + w["intermediate_state"]
                + w["triggerA_on_state"]) or 1.0
        quality = 100.0 * positive / wsum

        return ScoreCard(triggerB_activation=sB, intermediate_state=sInt,
                         triggerA_on_state=sOn, penalty=pen, total=total,
                         quality_percent=quality, details=details, flags=flags)
