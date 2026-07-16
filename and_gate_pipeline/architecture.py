"""Stage 3 -- AND-gate toehold-switch architecture builder.

The switch chains two hairpins in one transcript (5'->3'):

    [ secondary / inhibitory stem ] -- [ spacer a* ] -- [ primary / main stem ] -- linker -- reporter

Trigger B opens the upstream *secondary* stem; Trigger A then opens the
downstream *primary* stem, giving AND logic (Kim 2019).  The primary stem is a
Green-2026 Series-A hairpin (18 bp stem, ~6 nt invasion, conserved RBS/AUG).

Interpretation decisions (documented so a bench scientist can re-tune them)
--------------------------------------------------------------------------
* Every switch-side domain is built as the *exact reverse complement of the
  actual trigger domain* (Section 6, "Mismatch Handling"), so each trigger is
  captured with perfect complementarity even when x and k2 are only an
  approximate reverse-complement pair.  The residual x/k2 mismatch shows up as
  ``hamming`` mismatches inside the secondary stem -- the true biophysical cost.
* Secondary stem (Section 3):
      5' arm (base->loop):  k2*  then r1          (k2* = reverse_complement(k2))
      3' arm (loop->base):  r1*  then x*          (r1* , x* = reverse complements)
  r1 pairs r1*; k2* pairs x*; a 3-nt junction bulge (the "AUG-style" bulge) is
  left unpaired between the two sub-stems, exactly as specified.
* The OFF-state *lock* is what we encode as an explicit intended structure and
  score with SED.  Trigger-binding steps are scored with real cofold
  thermodynamics (see :mod:`.scoring`) rather than by assuming a binding
  orientation, so the numbers stay physical regardless of annotation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import sequence_utils as su
from .config import PipelineConfig
from .target_scan import TriggerPair


# A neutral, low-structure secondary-loop scaffold (no RBS, no AUG).  It is
# padded/trimmed to ``secondary_loop_len``.
_SECONDARY_LOOP_SCAFFOLD = "GAAACAGAACGAAUCAGACUUCGGAUCAG"
# Primary-stem body (fills the ascending arm to 18 nt).  The final 3 nt are
# fixed to CAU so the paired descending base reads AUG (the start codon).
_ASC_BODY_SEED = "GACUGACUGCUCACGUACCAU"


@dataclass
class DomainMap:
    """0-based [start, end) spans of every named domain within ``core``."""
    spans: dict = field(default_factory=dict)

    def add(self, name: str, start: int, length: int) -> int:
        self.spans[name] = (start, start + length)
        return start + length

    def region(self, *names: str) -> list[int]:
        idx: list[int] = []
        for n in names:
            if n in self.spans:
                s, e = self.spans[n]
                idx.extend(range(s, e))
        return idx


@dataclass
class DesignedSwitch:
    pair: TriggerPair
    cfg: PipelineConfig
    core: str                 # secondary + spacer + primary hairpin (no reporter)
    off_structure: str        # intended OFF-state dot-bracket over ``core``
    full: str                 # core + linker + reporter (translated construct)
    domains: DomainMap
    reporter: str = ""
    notes: list[str] = field(default_factory=list)

    # convenient sub-sequences ------------------------------------------- #
    @property
    def triggerA(self) -> str:
        return self.pair.triggerA.seq

    @property
    def triggerB(self) -> str:
        return self.pair.triggerB.seq

    def domain_seq(self, name: str) -> str:
        s, e = self.domains.spans[name]
        return self.core[s:e]

    def triggerA_footprint_indices(self) -> list[int]:
        """Positions on the core that are complementary to Trigger A domains."""
        return self.domains.region("sec_r1star", "sec_xstar", "spacer_astar",
                                    "prim_k1star")

    def triggerB_toehold_indices(self) -> list[int]:
        return self.domains.region("sec_k2star")

    def primary_rbs_indices(self) -> list[int]:
        return self.domains.region("prim_loop", "prim_descending")


def _fit_loop(scaffold: str, length: int) -> str:
    if length <= len(scaffold):
        return scaffold[:length]
    reps = (length // len(scaffold)) + 1
    return (scaffold * reps)[:length]


def build_switch(pair: TriggerPair, cfg: PipelineConfig,
                 reporter: str = "",
                 asc_body_override: str | None = None) -> DesignedSwitch:
    ta, tb = pair.triggerA, pair.triggerB
    r1, x, a, k1 = ta.r1, ta.x, ta.a, ta.k1
    k2 = tb.k2
    B = cfg.bulge_len

    rc = su.reverse_complement
    k2star = rc(k2)               # captures Trigger B (5' arm base)
    r1star = rc(r1)               # r1 clamp partner (3' arm, loop-proximal)
    xstar = rc(x)                 # complementary to x (3' arm, base)
    astar = rc(a)                 # spacer-a binding site (inter-stem)
    k1star = rc(k1)               # primary invasion toehold (ascending base)

    lr1 = len(r1)
    Lx = len(x)
    loop = _fit_loop(_SECONDARY_LOOP_SCAFFOLD, cfg.secondary_loop_len)

    # ---- secondary (inhibitory) stem ----------------------------------- #
    # 5' arm : k2*  | r1[:B] (bulge) | r1[B:]
    # 3' arm : r1*[: -B] | r1*[-B:] (bulge) | x*
    sec_5arm = k2star + r1
    sec_3arm = r1star + xstar
    secondary = sec_5arm + loop + sec_3arm

    # ---- spacer a* between the stems ----------------------------------- #
    # ---- primary (Series-A) hairpin ------------------------------------ #
    body_len = cfg.primary_stem_len - cfg.len_k1
    seed = su.to_rna(asc_body_override) if asc_body_override else _ASC_BODY_SEED
    asc_body = _fit_loop(seed, body_len)
    # guarantee the descending base spells the AUG start codon
    asc_body = asc_body[:-3] + su.reverse_complement("AUG")  # ...CAU
    ascending_paired = k1star + asc_body                      # 18 nt paired core
    aug_bulge = "GGA"[:B] if B <= 3 else ("GGA" + "A" * (B - 3))
    descending = su.reverse_complement(ascending_paired)      # 18 nt
    prim_loop = _primary_loop(cfg)

    primary = k1star + aug_bulge + asc_body + prim_loop + descending

    core = secondary + astar + primary

    # ---- intended OFF-state structure over the core -------------------- #
    off = _off_structure(cfg, lr1, Lx, len(loop), len(astar),
                         len(aug_bulge), len(asc_body), len(prim_loop))
    assert len(off) == len(core), (len(off), len(core))

    # ---- domain position map ------------------------------------------- #
    dm = DomainMap()
    p = 0
    p = dm.add("sec_k2star", p, len(k2star))
    p = dm.add("sec_r1", p, lr1)
    p = dm.add("sec_loop", p, len(loop))
    p = dm.add("sec_r1star", p, len(r1star))
    p = dm.add("sec_xstar", p, len(xstar))
    p = dm.add("spacer_astar", p, len(astar))
    p = dm.add("prim_k1star", p, len(k1star))
    p = dm.add("prim_aug_bulge", p, len(aug_bulge))
    p = dm.add("prim_asc_body", p, len(asc_body))
    p = dm.add("prim_loop", p, len(prim_loop))
    p = dm.add("prim_descending", p, len(descending))

    reporter_rna = su.to_rna(reporter) if reporter else ""
    if reporter_rna:
        reporter_rna = reporter_rna[:len(reporter_rna) - len(reporter_rna) % 3]
    full = core + su.to_rna(cfg.linker_suffix) + reporter_rna

    return DesignedSwitch(pair=pair, cfg=cfg, core=core, off_structure=off,
                          full=full, domains=dm, reporter=reporter_rna)


def _primary_loop(cfg: PipelineConfig) -> str:
    """RBS-bearing 5'UTR loop.  No AUG after the RBS (the only AUG downstream of
    the RBS must be the start codon in the descending stem)."""
    pre, post = "GGACUU", "UACA"
    loop = su.to_rna(pre + cfg.rbs_seq + post)
    return loop


def _off_structure(cfg, lr1, Lx, loop_len, astar_len,
                   aug_bulge_len, asc_body_len, prim_loop_len) -> str:
    B = cfg.bulge_len
    # secondary stem
    sec = (
        "(" * Lx                    # k2*  (pairs x*)
        + "." * B                   # r1[:B] bulge
        + "(" * (lr1 - B)           # r1[B:] (pairs r1*)
        + "." * loop_len            # secondary loop
        + ")" * (lr1 - B)           # r1*
        + "." * B                   # r1* bulge
        + ")" * Lx                  # x*
    )
    spacer = "." * astar_len        # a* single-stranded (Trigger A spacer site)
    # primary stem: k1*(6) [bulge] asc_body(12) loop descending(18)
    prim = (
        "(" * cfg.len_k1
        + "." * aug_bulge_len
        + "(" * asc_body_len
        + "." * prim_loop_len
        + ")" * (cfg.len_k1 + asc_body_len)
    )
    return sec + spacer + prim
