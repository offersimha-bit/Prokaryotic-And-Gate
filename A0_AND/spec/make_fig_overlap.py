"""Figure 3: the overlap trade-off - |x| vs trigger-trigger sequestration."""
import os, math, random, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import RNA

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures")
RT = 0.001987204258 * 310.15
CONC = 10e-9
random.seed(7)


def rc(s):
    return s.translate(str.maketrans("ACGU", "UGCA"))[::-1]


def seq_frac(dG, c=CONC):
    Kd = math.exp(dG / RT)
    return (((2 * c + Kd) - math.sqrt((2 * c + Kd) ** 2 - 4 * c * c)) / 2) / c


def spread_mismatch(partner, n):
    p = list(partner)
    L = len(p)
    for i in range(n):
        q = min(max(round((i + 0.5) * L / n) - 1, 0), L - 1)
        p[q] = [b for b in "ACGU" if b != p[q]][0]
    return "".join(p)


Ls = list(range(4, 19))
perfect, mism, runs = [], [], []
for L in Ls:
    nmm = L // 4
    dp, dm = [], []
    for _ in range(500):
        x = "".join(random.choice("ACGU") for _ in range(L))
        dp.append(RNA.duplexfold(x, rc(x)).energy)
        dm.append(RNA.duplexfold(x, spread_mismatch(rc(x), nmm)).energy)
    perfect.append(seq_frac(statistics.mean(dp)) * 100)
    mism.append(seq_frac(statistics.mean(dm)) * 100)
    runs.append(nmm)

fig, ax = plt.subplots(figsize=(9.0, 4.6))
ax.axvspan(5, 11, color="#8ECAE6", alpha=0.22, zorder=0)
ax.text(8, 88, "range tested\nby Kim 2019\n(* = 5, 8, 11)", ha="center", fontsize=8.5,
        color="#1B6E8C", fontweight="bold", linespacing=1.3)

ax.plot(Ls, perfect, "o-", color="#D00000", lw=2.4, ms=6,
        label="perfect complementarity over x")
ax.plot(Ls, mism, "s-", color="#2D6A4F", lw=2.4, ms=6,
        label="with mismatches at a 1 : 4 ratio, evenly spread")
ax.axhline(10, ls="--", lw=1.2, color="#868E96")
ax.text(18.4, 12, "10 %", fontsize=8, color="#868E96", ha="right")

ax.set_xlabel("|x|  —  length of the overlap between the two triggers (nt)", fontsize=10)
ax.set_ylabel("trigger A sequestered by trigger B  (%)", fontsize=10)
ax.set_title("Why the overlap has to be swept, not maximised",
             fontsize=11.5, fontweight="bold", pad=26)
ax.text(0.5, 1.015, "ViennaRNA duplex energies, 37 °C, both transcripts at 10 nM\n"
                    "mean over 500 random sequences per point",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=8.5,
        color="#495057", linespacing=1.35)
ax.set_xticks(Ls)
ax.set_ylim(-4, 104)
ax.set_xlim(3.5, 18.5)
ax.grid(alpha=0.25, ls=":")
ax.legend(fontsize=9, loc="center right", framealpha=0.95)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout()
out = os.path.join(OUTDIR, "fig3_overlap.png")
fig.savefig(out, dpi=220, facecolor="white", bbox_inches="tight")
print("wrote", out)
for L, p, m, n in zip(Ls, perfect, mism, runs):
    print(f"  |x|={L:2d}  perfect {p:6.1f}%   {n} mm -> {m:6.1f}%")
