# -*- coding: utf-8 -*-
"""Figure 4 - the RBS loop: three published choices, and what sets the SD-AUG spacing."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
DESIGNED = "#8ECAE6"
RBSCOL = "#F9C74F"
SDCOL = "#E76F51"
AUG = "#D00000"
STEM = "#4361EE"

RBS = "AACAGAGGAGA"
SD = "AGGAGA"
LOOPS = [("Green 2014 / VISTA", "AACAGAGGAGA", 11),
         ("Kim 2019", "CAAGAACAGAGGAGA", 15),
         ("ours (from the supervisor)", "AGACAAGAACAGAGGAGA", 18)]

fig = plt.figure(figsize=(13.4, 6.4))

# ---------------- Panel A : the three loops, aligned on their 3' end -------
axA = fig.add_axes([0.02, 0.50, 0.96, 0.44])
axA.set_xlim(0, 100); axA.set_ylim(0, 26); axA.axis("off")
axA.text(0, 24.2, "A", fontsize=15, fontweight="bold")
axA.text(4.2, 24.2, "The RBS loop in the three designs — aligned on their 3' end",
         fontsize=11.5, fontweight="bold", va="center")
axA.text(4.2, 21.0, "Every loop ends with the same 11-nt RBS. They differ only by a prefix, "
                    "and each is the previous one with more prefix added.",
         fontsize=8.8, color="#495057", va="center", style="italic")

CW = 2.55                      # cell width
RIGHT = 88.0                   # 3' end, common to all rows
for row, (name, seq, n) in enumerate(LOOPS):
    y = 13.2 - row * 5.4
    x0 = RIGHT - len(seq) * CW
    axA.text(x0 - 1.6, y + 1.55, name, fontsize=9, fontweight="bold", ha="right", va="center")
    axA.text(x0 - 1.6, y - 0.15, f"{n} nt", fontsize=8, color="#6C757D", ha="right", va="center")
    r = seq.find(RBS)
    sd = seq.find(SD)
    for i, ch in enumerate(seq):
        col = DESIGNED if i < r else RBSCOL
        if sd >= 0 and sd <= i < sd + len(SD):
            col = SDCOL
        axA.add_patch(FancyBboxPatch((x0 + i * CW, y), CW * 0.9, 3.1,
                                     boxstyle="round,pad=0,rounding_size=0.25",
                                     fc=col, ec="#22252A", lw=0.7))
        axA.text(x0 + i * CW + CW * 0.45, y + 1.55, ch, fontsize=9,
                 ha="center", va="center", fontweight="bold",
                 color=("white" if col == SDCOL else "#22252A"))
    if r > 0:
        axA.text(x0 + r * CW / 2, y + 3.9, f"{r} nt free / designed",
                 fontsize=7.6, ha="center", color="#0B7285", style="italic")
axA.text(RIGHT - 11 * CW / 2, 13.2 + 3.9, "fixed RBS, 11 nt — 3'-flush in all three",
         fontsize=7.8, ha="center", color="#8A6D0B", style="italic")
axA.add_patch(FancyBboxPatch((90.0, 1.0), 9.5, 16.0, boxstyle="round,pad=0,rounding_size=0.5",
                             fc="#F8F9FA", ec="#CED4DA", lw=1.0))
for i, (lab, col) in enumerate((("SD core\nAGGAGA", SDCOL), ("rest of\nthe RBS", RBSCOL),
                                ("designed\nflank (N)", DESIGNED))):
    axA.add_patch(FancyBboxPatch((91.2, 12.6 - i * 5.0), 2.1, 2.1,
                                 boxstyle="round,pad=0,rounding_size=0.25",
                                 fc=col, ec="#22252A", lw=0.7))
    axA.text(94.0, 13.65 - i * 5.0, lab, fontsize=7.4, va="center", linespacing=1.25)

# ---------------- Panel B : what sets the SD -> AUG spacing ---------------
axB = fig.add_axes([0.02, 0.02, 0.96, 0.44])
axB.set_xlim(0, 100); axB.set_ylim(0, 26); axB.axis("off")
axB.text(0, 24.2, "B", fontsize=15, fontweight="bold")
axB.text(4.2, 24.2, "Because the RBS is 3'-flush, the spacing from the SD to the start codon "
                    "is exactly len_k1",
         fontsize=11.5, fontweight="bold", va="center")

def layout(y, len_k1, len_mainpre, tag, colour, note):
    x = 12.0
    parts = [("RBS loop", 18, RBSCOL, "#22252A"),
             (f"MainZ = {len_k1}", len_k1, STEM, "white"),
             ("AUG", 3, AUG, "white"),
             (f"Main_pre = {len_mainpre}", len_mainpre, "#43AA8B", "white")]
    SC = 1.55
    axB.text(x - 1.6, y + 1.6, tag, fontsize=9.5, fontweight="bold", ha="right", va="center",
             color=colour)
    sd_end = None
    for lab, n, c, tc in parts:
        w = n * SC
        axB.add_patch(FancyBboxPatch((x, y), w, 3.2, boxstyle="round,pad=0,rounding_size=0.3",
                                     fc=c, ec="#22252A", lw=0.9))
        axB.text(x + w / 2, y + 1.6, lab, fontsize=8, ha="center", va="center",
                 color=tc, fontweight="bold")
        if lab == "RBS loop":
            axB.add_patch(FancyBboxPatch((x + w - 6 * SC, y), 6 * SC, 3.2,
                                         boxstyle="round,pad=0,rounding_size=0.3",
                                         fc=SDCOL, ec="#22252A", lw=0.9))
            axB.text(x + w - 3 * SC, y + 1.6, "AGGAGA", fontsize=7.4, ha="center",
                     va="center", color="white", fontweight="bold")
            sd_end = x + w
        x += w
    aug_x = sd_end + len_k1 * SC
    axB.annotate("", xy=(aug_x, y - 1.1), xytext=(sd_end, y - 1.1),
                 arrowprops=dict(arrowstyle="<->", lw=1.4, color=colour))
    axB.text((sd_end + aug_x) / 2, y - 2.9, f"spacing = {len_k1} nt", fontsize=8.6,
             ha="center", color=colour, fontweight="bold")
    axB.text(x + 2.0, y + 1.6, note, fontsize=8.2, va="center", color="#343A40",
             linespacing=1.3)

layout(14.5, 6, 12, "len_main_pre = 12", "#2D6A4F",
       "spacing 6 nt — the value used by Green,\nby Kim, and by all five of our tested switches")
layout(4.5, 9, 9, "len_main_pre = 9", "#A4133C",
       "spacing 9 nt — no published precedent\nin this scaffold")

out = os.path.join(OUTDIR, "fig4_loops.png")
fig.savefig(out, dpi=220, bbox_inches="tight", facecolor="white")
print("wrote", out)
