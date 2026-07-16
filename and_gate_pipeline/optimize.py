"""Stage 6 -- sequence-optimisation constraints.

Only a small part of the switch is free to redesign: the primary-stem body
(``prim_asc_body``, minus the fixed 3-nt CAU that codes the AUG start), and the
loop scaffolds.  Everything else is dictated by the trigger sequences.  The
optimiser makes greedy single-base substitutions in the free body positions to:

* remove forbidden homopolymer runs (AAAA/CCCC/GGGG/UUUU),
* remove in-frame stop codons downstream of the start codon,
* keep exactly one AUG after the RBS,
* nudge the OFF-state MFE toward the -54.25 kcal/mol low-leak target,

re-deriving the reverse-complementary descending stem after every change so the
primary stem stays paired.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from . import sequence_utils as su
from .architecture import build_switch, DesignedSwitch
from .config import PipelineConfig
from .thermo import ThermoBackend


@dataclass
class RestrictedReport:
    forbidden_runs: list[str]
    inframe_stop: bool
    aug_after_rbs: int
    ok: bool


def check_restricted(sw: DesignedSwitch, cfg: PipelineConfig) -> RestrictedReport:
    runs = su.has_forbidden_run(sw.core, cfg.forbidden_runs)
    # Count starts only within the switch module (5'UTR + start codon); a real
    # reporter CDS legitimately contains internal Met (AUG) codons downstream.
    aug_count = su.count_aug_after_rbs(sw.core, cfg.rbs_seq)
    # reading frame starts at the AUG that opens the descending stem
    start = _start_codon_index(sw, cfg)
    stop = su.has_inframe_stop(sw.full, start) if start is not None else True
    ok = (not runs) and (not stop) and (aug_count == 1)
    return RestrictedReport(forbidden_runs=runs, inframe_stop=stop,
                            aug_after_rbs=aug_count, ok=ok)


def _start_codon_index(sw: DesignedSwitch, cfg: PipelineConfig):
    rbs_idx = sw.full.find(su.to_rna(cfg.rbs_seq))
    augs = [i for i in su.find_all(sw.full, "AUG") if i > rbs_idx]
    return augs[0] if augs else None


def _score_violation(sw, cfg, backend, mfe_weight=0.02):
    rep = check_restricted(sw, cfg)
    v = 0
    v += 5 * len(rep.forbidden_runs)
    v += 10 if rep.inframe_stop else 0
    v += 8 * abs(rep.aug_after_rbs - 1)
    mfe = backend.mfe(sw.core)[1]
    v += mfe_weight * abs(mfe - cfg.off_state_mfe_target)
    return v, rep


def optimize_switch(sw: DesignedSwitch, cfg: PipelineConfig,
                    backend: ThermoBackend, max_iter: int = 120,
                    seed: int | None = None) -> tuple[DesignedSwitch, RestrictedReport]:
    """Greedy repair of the free primary-stem body.  Returns the best switch and
    its restricted-sequence report (may still be imperfect if a violation lives
    in a fixed, trigger-derived region such as k1)."""
    rng = random.Random(cfg.random_seed if seed is None else seed)
    best = sw
    best_v, best_rep = _score_violation(best, cfg, backend)
    if best_v == 0:
        return best, best_rep

    body_lo, body_hi = sw.domains.spans["prim_asc_body"]
    # the trailing 3 nt (CAU -> AUG start) are fixed; leave them alone
    free_positions = list(range(0, (body_hi - body_lo) - 3))
    seed_body = sw.domain_seq("prim_asc_body")

    def rebuild(free_body: str) -> DesignedSwitch:
        # builder replaces the trailing 3 nt with CAU (-> AUG start codon)
        return build_switch(sw.pair, cfg, reporter=su.to_dna(sw.reporter),
                            asc_body_override=free_body + "CAU")

    body = list(seed_body[:-3])
    for _ in range(max_iter):
        if not free_positions:
            break
        pos = rng.choice(free_positions)
        original = body[pos]
        for base in "ACGU":
            if base == original:
                continue
            trial = body[:]
            trial[pos] = base
            cand = rebuild("".join(trial))
            v, rep = _score_violation(cand, cfg, backend)
            if v < best_v:
                best, best_v, best_rep = cand, v, rep
                body[pos] = base
                break
        if best_v == 0:
            break
    return best, best_rep
