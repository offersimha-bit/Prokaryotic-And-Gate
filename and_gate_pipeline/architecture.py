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
# Retained only for ``asc_body_override``; the default primary stem now takes
# its upper arm from the conserved element (see ``_primary_parts``).
_ASC_BODY_SEED = "GACUGACUGCUCACGUACCAU"

# Watson-Crick + G:U wobble partners -- a base X pairs with P if X in _PAIRS[P].
_PAIRS = {"A": {"U"}, "U": {"A", "G"}, "G": {"C", "U"}, "C": {"G"}}


def _non_pairing_base(partner: str) -> str:
    """A base that cannot pair with ``partner`` (Watson-Crick or G:U wobble).

    Used to *build* a bulge into the sequence.  A bulge is a property of the
    sequence, not of the annotation: if the two arms are exact reverse
    complements the helix simply closes through the intended gap.
    """
    for b in "AGCU":
        if b not in _PAIRS.get(partner, set()) and partner not in _PAIRS.get(b, set()):
            return b
    return "A"                                          # pragma: no cover


def _weaken_arm(r1: str, bulge_len: int, gc_bias: float) -> tuple[str, list[int]]:
    """Return the switch's internal copy of r1 with mismatches introduced, plus
    the 0-based offsets (within r1) that were made non-pairing.

    Two effects, both acting **only on the switch's own r1 copy** on the 5' arm.
    ``r1*`` on the 3' arm is Trigger A's binding site and is never touched, so
    Trigger A keeps binding it with perfect complementarity (spec section 6).

    * ``bulge_len``: the first ``bulge_len`` nt after k2* are made non-pairing
      -> the section-3 junction bulge actually forms.
    * ``gc_bias`` (tunable 2, "binding strength of the upper segments of the
      secondary stem"): fraction of the remaining r1/r1* clamp additionally
      mismatched.  0.0 = perfect clamp (strongest); 1.0 = every clamp position
      broken (weakest).  Positions are spread evenly so the clamp weakens
      smoothly rather than losing one contiguous block.
    """
    seq = list(r1)
    broken: list[int] = []
    for i in range(min(bulge_len, len(r1))):
        seq[i] = _non_pairing_base(su.complement(r1[i]))
        broken.append(i)

    clamp = list(range(bulge_len, len(r1)))
    n_extra = int(round(max(0.0, min(1.0, gc_bias)) * len(clamp)))
    if n_extra and clamp:
        step = len(clamp) / n_extra
        for k in range(n_extra):
            i = clamp[min(len(clamp) - 1, int(k * step))]
            seq[i] = _non_pairing_base(su.complement(r1[i]))
            broken.append(i)
    return "".join(seq), sorted(set(broken))


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
    # 5' arm : k2*  | r1_sw   (switch's own r1 copy -- mismatched to form the
    #                          section-3 junction bulge, and to tune the clamp)
    # 3' arm : r1*  | x*      (r1* is Trigger A's site -- never modified)
    r1_sw, broken = _weaken_arm(r1, B, cfg.secondary_arm_gc_bias)
    sec_5arm = k2star + r1_sw
    sec_3arm = r1star + xstar
    secondary = sec_5arm + loop + sec_3arm

    # ---- spacer a* between the stems ----------------------------------- #
    # ---- primary (Series-A) hairpin ------------------------------------ #
    ascending_paired, aug_bulge, prim_loop, d_tail = _primary_parts(
        cfg, k1star, asc_body_override)
    asc_body = ascending_paired[len(k1star):]
    descending = su.reverse_complement(ascending_paired)       # 18 nt

    primary = k1star + aug_bulge + asc_body + prim_loop + descending + d_tail

    core = secondary + astar + primary

    # ---- intended OFF-state structure over the core -------------------- #
    off = _off_structure(cfg, lr1, Lx, len(loop), len(astar),
                         len(aug_bulge), len(asc_body), len(prim_loop),
                         broken, len(d_tail))
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
    if d_tail:
        p = dm.add("prim_d_domain", p, len(d_tail))

    reporter_rna = su.to_rna(reporter) if reporter else ""
    if reporter_rna:
        reporter_rna = reporter_rna[:len(reporter_rna) - len(reporter_rna) % 3]
    full = core + su.to_rna(cfg.linker_suffix) + reporter_rna

    return DesignedSwitch(pair=pair, cfg=cfg, core=core, off_structure=off,
                          full=full, domains=dm, reporter=reporter_rna)


def _primary_parts(cfg: PipelineConfig, k1star: str,
                   asc_body_override: str | None = None):
    """Series-A primary hairpin, taken from the conserved Toehold-VISTA element.

    ``cfg.hairpin_top`` (GUUAUAGUUAUGAACAGAGGAGACAUAACAUGAAC) decomposes exactly
    as VISTA uses it::

        [:12]   GUUAUAGUUAUG          conserved upper ascending stem
        [12:-3] AACAGAGGAGACAUAACAUG  the loop -- RBS at +3, start codon at its 3' end
        [-3:]   AAC                   the d domain

    so ``k1*(6) + hairpin_top[:12]`` gives an 18-bp stem invaded 6 nt by k1 --
    the Series-A geometry -- while the RBS/AUG/reading frame come from Green's
    validated element rather than from filler of our own.  Returns
    (ascending_paired, aug_bulge, loop, d_tail).
    """
    top = su.to_rna(cfg.hairpin_top)
    body_len = cfg.primary_stem_len - cfg.len_k1
    B = cfg.bulge_len

    if asc_body_override is not None:
        # optimiser path: caller supplies the upper arm; keep the CAU->AUG tail
        body = _fit_loop(su.to_rna(asc_body_override), body_len)
        body = body[:-3] + su.reverse_complement("AUG")
        loop = su.to_rna(_LEGACY_LOOP_PRE + cfg.rbs_seq + _LEGACY_LOOP_POST)
        d_tail = ""
    else:
        body = top[:body_len]                    # conserved upper stem
        loop = top[body_len:-3]                  # conserved loop: RBS ... AUG
        d_tail = top[-3:]                        # == cfg.d_domain

    aug_bulge = ("GGA" + "A" * max(0, B - 3))[:B]   # one-sided bulge after k1*
    return k1star + body, aug_bulge, loop, d_tail


_LEGACY_LOOP_PRE, _LEGACY_LOOP_POST = "GGACUU", "UACA"


def _off_structure(cfg, lr1, Lx, loop_len, astar_len,
                   aug_bulge_len, asc_body_len, prim_loop_len,
                   broken=(), d_tail_len=0) -> str:
    """Intended OFF-state dot-bracket.

    ``broken`` lists 0-based offsets within r1 that were made non-pairing, so
    the annotation matches the sequence we actually built (bulge + any clamp
    weakening from tunable 2) instead of asserting pairs that cannot form.
    """
    broken = set(broken)
    five = "".join("." if i in broken else "(" for i in range(lr1))
    three = "".join("." if (lr1 - 1 - i) in broken else ")" for i in range(lr1))
    sec = ("(" * Lx + five + "." * loop_len + three + ")" * Lx)
    spacer = "." * astar_len
    prim = ("(" * cfg.len_k1
            + "." * aug_bulge_len
            + "(" * asc_body_len
            + "." * prim_loop_len
            + ")" * (cfg.len_k1 + asc_body_len)
            + "." * d_tail_len)
    return sec + spacer + prim
