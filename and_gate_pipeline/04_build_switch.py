"""AND-gate switch built on VISTA's own toehold-switch construction.

The insight that makes this small: with the spec's own numbers,

    |r1| + |x| + |a|  =  20 + 6 + 4  =  30   <- Green 2026's toehold
    |k1|              =              6       <- Green 2026's invasion
    L_A               =             36       <- Green 2026's trigger

Trigger A's binding site on the switch is ``r1* | x* | a* | k1*`` (5'->3'),
which is exactly ``reverse_complement(TriggerA)``.  So feeding Trigger A into
VISTA's own builder as its 36-nt "target" reproduces our primary module
verbatim -- toehold(30) + k1*(6) + the conserved Series-A hairpin.  We do not
re-implement it; hand-rolling that element got the loop wrong three times.

Our contribution is only the Kim 2019 inhibitory hairpin, prepended:

    5'-[r2*][k2*][r1]-loop-[ r1* x* | a* ][k1*][conserved hairpin]-linker-CDS-3'
        \___/ \________ masks 26 of the 30-nt toehold ________/
        toehold                                    \__/ exposed gap = Kim's 'a'
        for B

OFF : the inhibitory hairpin pairs (k2*+r1) with (r1*+x*), so only a* (4 nt) of
      Trigger A's 30-nt toehold is exposed -- too short to fire.  Kim's a=4.
+B  : B binds the r2* toehold and invades k2*, which displaces the k2*:x*
      helix and releases x* -- and ONLY x*.  B has no complementarity to
      r1_switch, so the r1*:r1_switch helix stays closed.  Trigger A's
      NUCLEATION toehold therefore grows from |a| to |a|+Lx, never to the full
      30 nt.  This ceiling was measured, not assumed.
+A+B: Trigger A nucleates on x*+a* (|a|+Lx nt) and then branch-migrates through
      r1*, displacing the switch's own r1 copy, so the final duplex does cover
      all 30 nt -- but the RATE is set by the nucleation toehold, not by the
      final duplex size.  k1 then invades the 6-bp stem base and the RBS/AUG
      hairpin opens.  From that point it is a plain TSgen2 switch, which is why
      VISTA's trained model applies to it in-domain.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import sequence_utils as su
from .architecture import _weaken_arm          # detunes only the switch's r1 copy
from .config import PipelineConfig
from .target_scan import TriggerPair

_SECONDARY_LOOP_SCAFFOLD = "GAAACAGAACGAAUCAGACUUCGGAUCAG"


def _fit(scaffold: str, n: int) -> str:
    if n <= len(scaffold):
        return scaffold[:n]
    return (scaffold * (n // len(scaffold) + 1))[:n]


@dataclass
class VistaAndSwitch:
    pair: TriggerPair
    cfg: PipelineConfig
    core: str
    spans: dict
    primary_only: str          # the TSgen2 switch alone (what VISTA's model scores)
    reporter: str = ""
    notes: list = field(default_factory=list)
    broken: list = field(default_factory=list)
    """Offsets inside r1_switch deliberately mismatched (bulge / clamp detune),
    so the intended structure can annotate them as unpaired instead of asserting
    pairs that cannot form."""

    @property
    def full(self) -> str:
        """The translated construct.  ``core`` already ends with linker+CDS from
        VISTA's builder, so this is the same string -- exposed under the name the
        rest of the pipeline uses."""
        return self.core

    @property
    def off_structure(self) -> str:
        """Intended OFF-state dot-bracket over ``core``.

        Written from the spans rather than folded, so that comparing it to the
        actual fold (SED) tests whether the design we MEANT was built.  That
        comparison is what caught a bulge that had been annotated but never put
        into the sequence.
        """
        s, n = self.spans, len(self.core)
        db = ["."] * n

        def fill(name, ch):
            if name in s:
                for i in range(*s[name]):
                    db[i] = ch

        # inhibitory hairpin:  k2*( r1_sw( loop )r1* )x*
        fill("k2star", "(")
        fill("r1_switch", "(")
        fill("r1star", ")")
        fill("xstar", ")")
        r1s = s["r1_switch"][0]
        r1e = s["r1star"][1]
        for off in self.broken:                 # mismatched -> unpaired both sides
            if r1s + off < n:
                db[r1s + off] = "."
            mirror = r1e - 1 - off
            if s["r1star"][0] <= mirror < r1e:
                db[mirror] = "."

        # primary hairpin:  k1*( top0-3( top6-12( RBSloop )top23-29 )top32-35 )k1
        fill("k1star", "(")
        fill("ab_dom", ")")
        if "top" in s:
            t0 = s["top"][0]
            body = self.cfg.primary_stem_len - self.cfg.len_k1     # 12
            for i in range(t0, t0 + 3):
                db[i] = "("
            for i in range(t0 + 3, t0 + 6):
                db[i] = "."                      # 3x3 internal loop, ascending
            for i in range(t0 + 6, t0 + body):
                db[i] = "("
            for i in range(t0 + body, s["top"][1] - 12):
                db[i] = "."                      # RBS loop
            for i in range(s["top"][1] - 12, s["top"][1] - 6):
                db[i] = ")"
            for i in range(s["top"][1] - 6, s["top"][1] - 3):
                db[i] = "."                      # the AUG bulge (start codon)
            for i in range(s["top"][1] - 3, s["top"][1]):
                db[i] = ")"
        return "".join(db)

    def seq_of(self, name: str) -> str:
        s, e = self.spans[name]
        return self.core[s:e]

    @property
    def triggerA(self) -> str:
        return self.pair.triggerA.seq

    @property
    def triggerB(self) -> str:
        return self.pair.triggerB.seq


def build_primary_module(target: str, cfg: PipelineConfig, reporter: str = "") -> tuple:
    """VISTA's toehold-switch construction, verbatim, for a 36-nt target.

    Mirrors Toehold_VISTA.ipynb::return_design_object:
        hairpin = toehold + b*_dom + a*_dom + hairpin_top + a_dom + b_dom
                  + d_dom + a*_dom + rc(hairpin_top[-3:]) + suffix + output[:30]
    Returns (sequence, spans, toehold_len).
    """
    rc = su.reverse_complement
    target = su.to_rna(target)
    toehold_len = cfg.resolved_len_r1() + cfg.Lx + cfg.len_a      # == 30

    target_comp = rc(target)                    # r1* x* a* k1*
    toehold = target_comp[:toehold_len]         # r1* + x* + a*
    a_dom, b_dom = target[:3], target[3:6]      # first 6 nt of the target == k1
    a_star, b_star = rc(a_dom), rc(b_dom)
    top = su.to_rna(cfg.hairpin_top)
    d_dom = su.to_rna(cfg.d_domain)

    parts = [("toehold", toehold), ("k1star", b_star + a_star), ("top", top),
             ("ab_dom", a_dom + b_dom), ("d_dom", d_dom), ("a_star", a_star),
             ("top_tail_rc", rc(top[-3:]))]
    seq, spans, p = "", {}, 0
    for name, s in parts:
        spans[name] = (p, p + len(s)); p += len(s); seq += s
    rep = su.to_rna(reporter)
    if rep.startswith("AUG"):
        rep = rep[3:]                            # the switch supplies the start codon
    rep = rep[:len(rep) - len(rep) % 3]
    tail = su.to_rna(cfg.linker_suffix) + rep
    spans["linker_cds"] = (p, p + len(tail)); seq += tail
    return seq, spans, toehold_len


def build(pair: TriggerPair, cfg: PipelineConfig, reporter: str = "") -> VistaAndSwitch:
    """Kim inhibitory hairpin + VISTA primary module."""
    rc = su.reverse_complement
    ta, tb = pair.triggerA, pair.triggerB

    primary, pspans, toehold_len = build_primary_module(ta.seq, cfg, reporter)

    # --- the inhibitory hairpin (our only addition) --------------------- #
    r2star = rc(tb.r2)                 # single-stranded 5' toehold for Trigger B
    k2star = rc(tb.k2)                 # Lx nt; pairs x* -> Trigger B invades here
    # The switch's OWN copy of r1 (pairs r1*).  This is the only part of the
    # inhibitory hairpin we are free to detune: r1* is Trigger A's binding site
    # and must stay an exact reverse complement of the real trigger.
    r1_sw, broken = _weaken_arm(ta.r1, cfg.secondary_bulge_len,
                                cfg.secondary_arm_gc_bias)
    loop = _fit(_SECONDARY_LOOP_SCAFFOLD, cfg.secondary_loop_len)
    secondary = r2star + k2star + r1_sw + loop

    core = secondary + primary
    off = len(secondary)

    spans, p = {}, 0
    for name, s in (("r2star", r2star), ("k2star", k2star),
                    ("r1_switch", r1_sw), ("sec_loop", loop)):
        spans[name] = (p, p + len(s)); p += len(s)
    for name, (s, e) in pspans.items():
        spans[name] = (s + off, e + off)

    # sub-spans inside the 30-nt toehold: r1* | x* | a*
    t0 = spans["toehold"][0]
    lr1 = cfg.resolved_len_r1()
    spans["r1star"] = (t0, t0 + lr1)
    spans["xstar"] = (t0 + lr1, t0 + lr1 + cfg.Lx)
    spans["a_star_gap"] = (t0 + lr1 + cfg.Lx, t0 + toehold_len)
    # the primary hairpin proper = k1* .. ab_dom
    spans["primary_hairpin"] = (spans["k1star"][0], spans["ab_dom"][1])
    # what the inhibitory hairpin masks
    spans["masked_toehold"] = (t0, t0 + lr1 + cfg.Lx)

    # RBS and start codon, located inside the conserved element rather than
    # hard-coded, so they follow cfg.hairpin_top if it is ever changed.
    top_seq = su.to_rna(cfg.hairpin_top)
    tp0 = spans["top"][0]
    rbs_off = top_seq.find(su.to_rna(cfg.rbs_seq))
    if rbs_off >= 0:
        spans["rbs"] = (tp0 + rbs_off, tp0 + rbs_off + len(cfg.rbs_seq))
    aug_off = top_seq.rfind("AUG")          # the start codon is the LAST AUG
    if aug_off >= 0:
        spans["start_codon"] = (tp0 + aug_off, tp0 + aug_off + 3)

    return VistaAndSwitch(pair=pair, cfg=cfg, core=core, spans=spans,
                          primary_only=primary, reporter=su.to_rna(reporter),
                          broken=list(broken))
