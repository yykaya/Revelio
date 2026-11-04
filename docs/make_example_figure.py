#!/usr/bin/env python3
"""Draw docs/example_sites.png from the example run (requires matplotlib >= 3.6).

  ./revelio -i examples/col0_chr1_3631-4630.fa --model weights.tsv --header \
      -o examples/col0_chr1_3631-4630.model.tsv
  python3 docs/make_example_figure.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FASTA = os.path.join(ROOT, "examples", "col0_chr1_3631-4630.fa")
CALLS = os.path.join(ROOT, "examples", "col0_chr1_3631-4630.model.tsv")
OUT = os.path.join(ROOT, "docs", "example_sites.png")

ZOOM = (37, 72)  # 1-based inclusive window shown in panel a

COLORS = {"CG": "#2a78d6", "CHG": "#eb6834", "CHH": "#1baf7a"}
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COMP = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    with open(path) as f:
        return "".join(l.strip() for l in f if not l.startswith(">")).upper()


def read_calls(path):
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        r["pos"] = int(r["pos_1based"])
        r["prob"] = float(r["prob"])
    return rows


def panel_zoom(ax, seq, calls):
    lo, hi = ZOOM
    sites = {(r["pos"], r["strand"]): r["context"] for r in calls if lo <= r["pos"] <= hi}
    y_top, y_bot = 0.56, 0.30
    for i, pos in enumerate(range(lo, hi + 1)):
        x = i
        for strand, y, base in (("+", y_top, seq[pos - 1]), ("-", y_bot, seq[pos - 1].translate(COMP))):
            ctx = sites.get((pos, strand))
            if ctx:
                ax.add_patch(FancyBboxPatch((x - 0.42, y - 0.11), 0.84, 0.22,
                                            boxstyle="round,pad=0,rounding_size=0.12",
                                            fc=COLORS[ctx], ec="none"))
            ax.text(x, y, base, ha="center", va="center", family="monospace", fontsize=11,
                    color="white" if ctx else INK, fontweight="bold" if ctx else "normal")
        ctx_p, ctx_m = sites.get((pos, "+")), sites.get((pos, "-"))
        # context names rotated so neighbouring sites never collide
        if ctx_p:
            ax.text(x, y_top + 0.14, ctx_p, ha="center", va="bottom", rotation=90, fontsize=7.5, color=INK2)
        if ctx_m:
            ax.text(x, y_bot - 0.14, ctx_m, ha="center", va="top", rotation=90, fontsize=7.5, color=INK2)
    n = hi - lo + 1
    ax.text(-1.0, y_top, "+ 5'", ha="right", va="center", fontsize=9, color=INK2)
    ax.text(n, y_top, "3'", ha="left", va="center", fontsize=9, color=INK2)
    ax.text(-1.0, y_bot, "− 3'", ha="right", va="center", fontsize=9, color=INK2)
    ax.text(n, y_bot, "5'", ha="left", va="center", fontsize=9, color=INK2)
    for i in range(0, n, 5):
        ax.text(i, 0.93, str(lo + i), ha="center", va="bottom", fontsize=7, color=INK2)
        ax.plot([i, i], [0.89, 0.92], color=INK2, lw=0.6)
    ax.set_xlim(-3.2, n + 1.2)
    ax.set_ylim(-0.02, 1.0)
    ax.axis("off")
    ax.set_title(f"a  Both strands are scanned (positions {lo}–{hi})", loc="left",
                 fontsize=11, color=INK, fontweight="bold")


def panel_track(ax, seq, calls):
    lo, hi = ZOOM
    ax.add_patch(Rectangle((lo - 0.5, -1.05), hi - lo + 1, 2.1, fc=GRID, ec="none", alpha=0.6, zorder=0))
    for ctx in ("CHH", "CHG", "CG"):  # draw CG last so it sits on top
        pts = [r for r in calls if r["context"] == ctx]
        xs = [r["pos"] for r in pts]
        ys = [r["prob"] if r["strand"] == "+" else -r["prob"] for r in pts]
        ax.vlines(xs, 0, ys, color=COLORS[ctx], lw=1.0, zorder=2)
        n_plus = sum(r["strand"] == "+" for r in pts)
        ax.plot(xs, ys, "o", ms=3.6, color=COLORS[ctx], mec=SURFACE, mew=0.6, zorder=3,
                label=f"{ctx}  ({n_plus} + / {len(pts) - n_plus} −)")
    ax.axhline(0, color=INK2, lw=0.8, zorder=1)
    ax.set_xlim(0, len(seq) + 1)
    ax.set_ylim(-1.05, 1.05)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticklabels(["1.0", "0.5", "0", "0.5", "1.0"])
    ax.text(-0.065, 0.75, "+ strand", transform=ax.transAxes, rotation=90, ha="center", va="center",
            fontsize=9, color=INK2)
    ax.text(-0.065, 0.25, "− strand", transform=ax.transAxes, rotation=90, ha="center", va="center",
            fontsize=9, color=INK2)
    ax.set_xlabel("Position in example region (bp; region = TAIR10 Chr1:3,631–4,630)", color=INK2, fontsize=9)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    leg = ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=9,
                    title="Context (sites)", title_fontsize=9, alignment="left", handletextpad=0.3)
    for t in leg.get_texts():
        t.set_color(INK)
    ax.set_title("b  Prior probability per site (--model weights.tsv); grey band = panel a",
                 loc="left", fontsize=11, color=INK, fontweight="bold")


def main():
    seq = read_fasta(FASTA)
    calls = read_calls(CALLS)
    fig = plt.figure(figsize=(11, 6.8), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.25, 2.1], hspace=0.30,
                          left=0.08, right=0.80, top=0.93, bottom=0.09)
    panel_zoom(fig.add_subplot(gs[0]), seq, calls)
    ax = fig.add_subplot(gs[1])
    ax.set_facecolor(SURFACE)
    panel_track(ax, seq, calls)
    fig.savefig(OUT, dpi=180, facecolor=SURFACE)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
