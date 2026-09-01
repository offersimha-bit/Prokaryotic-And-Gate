"""Figure 1 for the A0_AND specification - Kim 2019 style architecture schematic."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Arc, FancyArrowPatch

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")

C = dict(r2="#2A9D8F", x="#E76F51", k2="#9B5DE5", mainpre="#43AA8B",
         k1="#4361EE", bulge="#98A2AD", aug="#D00000", rbs="#F9C74F",
         sloop="#8ECAE6", link="#CED4DA", gfp="#1B4332", cap="#6C757D")
NTH = 1.55          # vertical units per nt
W = 11.0            # arm width
GAP = 5.5           # gap between the two arms of a hairpin


def box(ax, x, y, w, h, color, label, sub=None, fs=8, tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=0.5",
                                fc=color, ec="#22252A", lw=1.0, zorder=3))
    t = label if sub is None else label + "\n" + sub
    ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=fs,
            color=tc, zorder=4, fontweight="bold", linespacing=1.05)


def rungs(ax, xl, xr, y0, y1):
    n = max(2, int(round((y1 - y0) / NTH)))
    for i in range(n):
        y = y0 + (i + 0.5) * (y1 - y0) / n
        ax.plot([xl, xr], [y, y], color="#5C636A", lw=0.8, zorder=2)


def _fs(lab):
    """shrink the font for long domain names so they stay inside the box"""
    return 7.6 if len(lab) <= 11 else 6.6


def hairpin(ax, cx, y0, arm5, arm3, loop_label, loop_col, note=None, gap=GAP):
    """arm5: bottom->top [(label, nt, color, in_internal_loop)]; arm3: top->bottom.

    Elements flagged `in_internal_loop` are unpaired and are drawn *between* the
    two arms - which is what a symmetric internal loop actually is.
    """
    xl, xr = cx - gap / 2 - W, cx + gap / 2
    bw = min(W * 0.62, gap / 2 - 0.6)
    y, spans5 = y0, []
    for lab, nt, col, loop in arm5:
        h = nt * NTH
        if loop:
            box(ax, xl + W + 0.3, y, bw, h, col, lab, fs=6.8)
        else:
            box(ax, xl, y, W, h, col, lab, sub=str(nt) + " nt", fs=_fs(lab))
            spans5.append((y, y + h))
        y += h
    top5 = y
    y, spans3 = y0, []
    for lab, nt, col, loop in reversed(arm3):
        h = nt * NTH
        if loop:
            box(ax, xr - 0.3 - bw, y, bw, h, col, lab, fs=6.8)
        else:
            box(ax, xr, y, W, h, col, lab, sub=str(nt) + " nt", fs=_fs(lab))
            spans3.append((y, y + h))
        y += h
    top3 = y
    for (a, b) in spans5:
        for (c, d) in spans3:
            lo, hi = max(a, c), min(b, d)
            if hi - lo > 0.05:
                rungs(ax, xl + W, xr, lo, hi)
    top = max(top5, top3)
    ax.add_patch(Arc((cx, top), gap + W, 11.0, theta1=0, theta2=180,
                     lw=11, color=loop_col, zorder=3))
    ax.text(cx, top + 7.6, loop_label, ha="center", va="center", fontsize=8.5,
            fontweight="bold", color="#22252A", zorder=5, linespacing=1.2)
    if note:
        ax.text(cx, y0 - 3.2, note, ha="center", va="top", fontsize=8,
                style="italic", color="#495057", linespacing=1.25)
    return xl, xr + W, top


fig = plt.figure(figsize=(13.4, 11.2))
ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
ax.set_xlim(0, 152)
ax.set_ylim(-30, 84)
ax.axis("off")

# ---------------- Panel A : OFF-state switch ----------------
ax.text(1, 76, "A", fontsize=18, fontweight="bold")
ax.text(7.5, 76, "AND switch RNA  —  OFF state (no trigger bound)",
        fontsize=13, fontweight="bold", va="center")

Y0 = 30
MAIN_GAP = 14.0
sec_l, sec_r, sec_top = hairpin(
    ax, 52, Y0,
    arm5=[("Secondary_pre*", 9, C["x"], False), ("k2*", 9, C["k2"], False)],
    arm3=[("SecondaryZ", 9, C["k2"], False), ("x*", 9, C["x"], False)],
    loop_label="Secondary loop\n(no RBS)", loop_col=C["sloop"],
    note="secondary (inhibitory) hairpin\n18 nt per arm = 18 bp (no bulge)")
main_l, main_r, main_top = hairpin(
    ax, 87, Y0, gap=MAIN_GAP,
    arm5=[("Main_pre*", 9, C["mainpre"], False),
          ("oppos.\nbulge\n3", 3, C["bulge"], True),
          ("k1*", 6, C["k1"], False)],
    arm3=[("MainZ", 6, C["k1"], False), ("AUG", 3, C["aug"], True),
          ("Main_pre", 9, C["mainpre"], False)],
    loop_label="RBS loop (18 nt)", loop_col=C["rbs"],
    note="main (switch) hairpin\n18 nt per arm = 15 bp + 3-nt bulge")
ax.annotate("3 × 3 internal loop\n(unpaired)",
            xy=(87 + MAIN_GAP / 2 + W - 1, Y0 + 13.5 * NTH),
            xytext=(87 + MAIN_GAP / 2 + W + 6, Y0 + 9.0 * NTH),
            fontsize=7.2, ha="left", va="center", color="#495057", style="italic",
            linespacing=1.25, zorder=6,
            arrowprops=dict(arrowstyle="->", color="#868E96", lw=1.1))

# upper-stem placeholders: label the secondary one on the left, the main one on the right
ax.plot([52 - GAP / 2 - W, 52 + GAP / 2 + W], [sec_top, sec_top],
        ls=(0, (2.5, 2.5)), color="#495057", lw=1.6, zorder=6)
ax.text(52 - GAP / 2 - W - 1.6, sec_top, "SUSA / SUSD = 0 nt\n(placeholder)",
        fontsize=7.6, va="center", ha="right", color="#495057", linespacing=1.25)
ax.plot([87 - MAIN_GAP / 2 - W, 87 + MAIN_GAP / 2 + W], [main_top, main_top],
        ls=(0, (2.5, 2.5)), color="#495057", lw=1.6, zorder=6)
ax.text(87 + MAIN_GAP / 2 + W + 1.6, main_top, "MUSA / MUSD = 0 nt\n(placeholder)",
        fontsize=7.6, va="center", color="#495057", linespacing=1.25)

box(ax, 4, Y0 - 2.4, 7.0, 4.8, C["cap"], "GGG", fs=8)
box(ax, 11, Y0 - 2.4, 22.5, 4.8, C["r2"], "r2*   trigger-B toehold", sub="free", fs=8)
ax.plot([33.5, sec_l], [Y0, Y0], color="#22252A", lw=1.8, zorder=1)
ax.text(1.8, Y0 + 4.4, "5'", fontsize=12, fontweight="bold")

ax.plot([sec_r, main_l], [Y0, Y0], color="#22252A", lw=2.4, zorder=1)
mid = (sec_r + main_l) / 2
ax.annotate("a = 0 nt", xy=(mid, Y0 - 0.4), xytext=(mid, Y0 - 15.0),
            ha="center", fontsize=11.5, fontweight="bold", color="#D00000",
            arrowprops=dict(arrowstyle="->", color="#D00000", lw=1.8))
ax.text(mid, Y0 - 18.0, "the two hairpins are directly adjacent.\n"
                        "With no exposed spacer, trigger A has\nno landing site in state 10.",
        ha="center", va="top", fontsize=8, color="#D00000", style="italic",
        linespacing=1.35)

box(ax, main_r, Y0 - 2.4, 17, 4.8, C["link"], "LINKER", sub="21 nt", fs=8, tc="#22252A")
box(ax, main_r + 17, Y0 - 2.4, 21, 4.8, C["gfp"], "GFP CDS", sub="759 nt", fs=8)
ax.text(main_r + 39.4, Y0 + 4.4, "3'", fontsize=12, fontweight="bold")

# ---------------- Panel B : triggers ----------------
ax.text(1, 3, "B", fontsize=18, fontweight="bold")
ax.text(7.5, 3, "Trigger RNAs (5'→3').  A trigger binds antiparallel: its 5' end pairs "
                "with the 3'-most element of its footprint.",
        fontsize=11, fontweight="bold", va="center")


def strand(ax, x0, y, segs, name, total):
    ax.text(x0 - 1.8, y + 2.5, "5'", fontsize=11, fontweight="bold", ha="right")
    x = x0
    for lab, nt, col in segs:
        w = max(nt * 1.30, 6.0)
        tc = "#22252A" if col == "#F1F3F5" else "white"
        box(ax, x, y, w, 5.0, col, lab, sub=str(nt), fs=7.6, tc=tc)
        x += w
    ax.text(x + 1.6, y + 2.5, "3'", fontsize=11, fontweight="bold")
    ax.text(x0 - 5.6, y + 2.5, name, fontsize=10, fontweight="bold",
            ha="right", va="center", linespacing=1.25)
    ax.text(x + 6.0, y + 2.5, str(total) + " nt", fontsize=9.5, va="center",
            color="#495057", fontweight="bold")


strand(ax, 33, -8, [("k1", 6, C["k1"]), ("oppos.\nbulge", 3, C["bulge"]),
                    ("Main_pre", 9, C["mainpre"]), ("a", 0, "#F1F3F5"),
                    ("x", 9, C["x"])],
       "Trigger A\n(gene A, main)", 27)
strand(ax, 33, -16, [("k2", 9, C["k2"]), ("Secondary_pre = x*", 9, C["x"]),
                     ("r2", 32, C["r2"])],
       "Trigger B\n(gene B, secondary)", 50)
ax.text(33, -20.5, "The two triggers share the x domain — trigger A carries x, trigger B carries x*.  "
                   "That overlap is what couples them;\nit is also what lets them bind each other "
                   "instead of the switch.",
        fontsize=8.4, color="#7B2D26", style="italic", linespacing=1.4, va="top")

out = os.path.join(OUTDIR, "fig1_architecture.png")
fig.savefig(out, dpi=220, bbox_inches="tight", facecolor="white")
print("wrote", out)

# ================= Figure 2 : mechanism / truth table =================
fig2 = plt.figure(figsize=(13.0, 3.5))
ax2 = fig2.add_axes([0.01, 0.02, 0.98, 0.96])
ax2.set_xlim(0, 152)
ax2.set_ylim(0, 34)
ax2.axis("off")
steps = [("00", "no trigger", "both hairpins closed;\nx* is paired", "OFF", "#ADB5BD"),
         ("10", "trigger A only", "x* still paired — trigger A\nhas nothing to nucleate on", "OFF", "#ADB5BD"),
         ("01", "trigger B only", "secondary hairpin opens,\nx* is freed; AUG still locked", "OFF", "#F9C74F"),
         ("11", "A  +  B", "A nucleates on the freed x*,\ndisplaces MainZ + Main_pre", "ON", "#2D6A4F")]
x = 2
for i, (st, inp, txt, out_s, col) in enumerate(steps):
    ax2.add_patch(FancyBboxPatch((x, 3), 33.5, 26,
                                 boxstyle="round,pad=0,rounding_size=1.2",
                                 fc=("#EBF7F0" if out_s == "ON" else "#F8F9FA"),
                                 ec=col, lw=2.6))
    ax2.text(x + 2.2, 24.6, "state " + st, fontsize=12, fontweight="bold", color="#22252A")
    ax2.text(x + 31.3, 24.6, out_s, fontsize=13, fontweight="bold", ha="right",
             color=("#2D6A4F" if out_s == "ON" else "#868E96"))
    ax2.text(x + 2.2, 19.4, inp, fontsize=9.5, color="#495057", fontweight="bold")
    ax2.text(x + 2.2, 11.5, txt, fontsize=9, va="center", color="#343A40", linespacing=1.45)
    if i < 3:
        ax2.add_patch(FancyArrowPatch((x + 34.2, 16), (x + 36.8, 16),
                                      arrowstyle="-|>", mutation_scale=16,
                                      lw=1.8, color="#868E96"))
    x += 37.5
out2 = os.path.join(OUTDIR, "fig2_truthtable.png")
fig2.savefig(out2, dpi=220, bbox_inches="tight", facecolor="white")
print("wrote", out2)
