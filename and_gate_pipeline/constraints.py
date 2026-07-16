"""Stage 4 -- system constraints and logical-integrity checks.

These are the *closed-form* checks (lengths and equations from Section 4).
Thermodynamic requirements (secondary stem stronger than primary, strong 'a'
site, trigger accessibility) are evaluated later with the folding engine and
live in :mod:`.scoring`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import PipelineConfig


@dataclass
class IntegrityReport:
    ok: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    derived: dict = field(default_factory=dict)

    def __str__(self) -> str:
        lines = [f"logical integrity: {'OK' if self.ok else 'CONTRADICTION'}"]
        for e in self.errors:
            lines.append(f"  ERROR   {e}")
        for w in self.warnings:
            lines.append(f"  warning {w}")
        for k, v in self.derived.items():
            lines.append(f"  {k} = {v}")
        return "\n".join(lines)


def validate_config(cfg: PipelineConfig) -> IntegrityReport:
    """Verify the characterisation is well defined and free of contradictions
    (Section 5.1)."""
    errors: list[str] = []
    warnings: list[str] = []

    lr1 = cfg.resolved_len_r1()
    lr2 = cfg.resolved_len_r2()

    # --- integer / non-negative lengths ---------------------------------- #
    lengths = {
        "Lx": cfg.Lx, "len_a": cfg.len_a, "len_k1": cfg.len_k1,
        "len_r1": lr1, "len_r2": lr2, "L_A": cfg.L_A, "L_B": cfg.L_B,
        "primary_stem_len": cfg.primary_stem_len,
        "secondary_loop_len": cfg.secondary_loop_len, "bulge_len": cfg.bulge_len,
    }
    for name, val in lengths.items():
        if not isinstance(val, int):
            errors.append(f"{name} must be an integer (got {val!r})")
        elif val < 0:
            errors.append(f"{name} must be non-negative (got {val})")

    # --- Section 4 equations --------------------------------------------- #
    # |a| + |x| + |k1| + |r1| == L_A
    lhs_a = cfg.len_a + cfg.Lx + cfg.len_k1 + lr1
    if lhs_a != cfg.L_A:
        errors.append(
            f"|a|+|x|+|k1|+|r1| = {lhs_a} != L_A = {cfg.L_A}")
    # |k2| + |r2| == L_B  (|k2| == |x| == Lx)
    lhs_b = cfg.Lx + lr2
    if lhs_b != cfg.L_B:
        errors.append(f"|k2|+|r2| = {lhs_b} != L_B = {cfg.L_B}")

    # --- structural feasibility ------------------------------------------ #
    if lr1 <= cfg.bulge_len:
        errors.append(
            f"len_r1 ({lr1}) must exceed bulge_len ({cfg.bulge_len}) so the "
            f"r1/r1* clamp survives the junction bulge")
    if cfg.len_k1 > cfg.primary_stem_len:
        errors.append(
            f"len_k1 ({cfg.len_k1}) cannot exceed primary_stem_len "
            f"({cfg.primary_stem_len})")
    if cfg.len_k1 < cfg.primary_invasion:
        warnings.append(
            f"len_k1 ({cfg.len_k1}) < primary_invasion ({cfg.primary_invasion}); "
            f"k1 will not fully seat the {cfg.primary_invasion}-nt invasion toehold")

    # RBS must not fall inside the secondary loop, and r1 must stay out of it.
    if cfg.rbs_seq in "N" * cfg.secondary_loop_len:
        warnings.append("secondary loop may be long enough to harbour an RBS")

    if cfg.Lx <= 0:
        errors.append("Lx must be positive")

    derived = {"len_r1": lr1, "len_r2": lr2,
               "trigger_A_footprint": lr1 + cfg.Lx + cfg.len_a + cfg.len_k1,
               "trigger_B_footprint": lr2 + cfg.Lx}

    return IntegrityReport(ok=not errors, errors=errors,
                           warnings=warnings, derived=derived)
