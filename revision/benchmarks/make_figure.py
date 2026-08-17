#!/usr/bin/env python
"""Supplementary scaling figure for the PoolParty manuscript.

Three panels sharing a log x-axis:
  A  wall-clock generation time vs library size (log-log)
  B  throughput vs library size (log-log)  -- the plateau IS the linearity claim
  C  peak resident memory vs library size (log x, linear y)

Reads paper_benchmark.json, writes figS_scaling.pdf (vector, for print) and
figS_scaling.png (300 dpi, for preview).

A and B are the same measurement expressed two ways: B = N / A. They are both
shown because readers use them differently -- A answers "how long will my run
take", B answers "does cost per sequence degrade with size". Carrying B in the
figure means the throughput table can be dropped from the supplement.

Deliberately minimal: no in-plot annotations, no fitted exponents, no per-point
values, no special markers. Every interpretive statement belongs in the caption.
Only integer decades are plotted, plus each example's own library size.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent

# Categorical slots 1-3 of the validated palette. Validated as a set for
# all-pairs separation on a light surface: worst CVD dE 9.2, worst
# normal-vision dE 24.0. Aqua sits below 3:1 contrast, so the relief rule
# applies -- discharged by the companion table, which carries the values in text.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2 = "#0b0b0b", "#52514e"

# Presentation order: largest library first. Colours stay bound to the example,
# not to its position in this list -- reordering the series must not repaint
# them, or a reader comparing against an earlier draft sees the wrong lines.
SERIES = [
    ("dms",      "GB1 deep mutational scan", AQUA),
    ("mpra",     "MPRA regulatory grammar",  BLUE),
    ("spliceai", "SpliceAI surrogate",       ORANGE),
]

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9,
    "axes.linewidth": 0.6,
    "axes.edgecolor": INK2,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "pdf.fonttype": 42,      # TrueType, not Type 3 -- journal requirement
    "ps.fonttype": 42,
    "savefig.bbox": "tight",
})


def decade_points(records):
    """Integer decades only, plus the example's own library size."""
    biggest = max(r["n"] for r in records)
    keep = []
    for r in sorted(records, key=lambda z: z["n"]):
        lg = np.log10(r["n"])
        if abs(lg - round(lg)) < 1e-9 or r["n"] == biggest:
            keep.append(r)
    return keep


def render_table(recs):
    """Table S1 as an image, for showing around before it is set in LaTeX.

    Drawn as text on a bare axes rather than with plt.table: booktabs-style
    horizontal rules only, right-aligned numerics, no vertical lines.
    """
    rows = []
    for key, label, _ in SERIES:
        r = next(x for x in recs if x["example"] == key and x["native"])
        # The SpliceAI library is two matched pools of this size; that is stated
        # in the caption rather than crowding the Example column here.
        rows.append([
            label,
            f"{r['n']:,}",
            f"{r['gen_mean']:.2f} ± {r['gen_sd']:.2f}",
            f"{r['n'] / r['gen_mean']:,.0f}",
            f"{r['peak_mb_mean']:.0f}",
            f"{r['build_mean']:.2f}",
        ])

    head = ["Example", "Library size", "Generation\ntime (s)",
            "Throughput\n(seq s$^{-1}$)", "Peak\nmemory (MB)", "DAG\nconstruction (s)"]
    # One x per column, used for BOTH the header and the values. Centring the
    # header on a span while right-aligning values to that span's edge leaves
    # the two visibly offset -- measured at 30-160 px on the rendered PNG, and
    # worsening left to right. A single shared anchor makes that impossible.
    centres = [0.00, 0.300, 0.455, 0.610, 0.755, 0.905]
    aligns = ["left"] + ["center"] * 5

    fig, ax = plt.subplots(figsize=(7.4, 1.55))
    ax.axis("off")
    y_head, y0, dy = 0.66, 0.44, 0.155
    y_top, y_mid = 0.97, 0.575
    y_bot = y0 - 2 * dy - 0.085

    for x, h, al in zip(centres, head, aligns):
        ax.text(x, y_head, h, ha=al, va="bottom", fontsize=7.8, color=INK,
                fontweight="bold", linespacing=1.25, transform=ax.transAxes)
    for i, row in enumerate(rows):
        y = y0 - i * dy
        for x, cell, al in zip(centres, row, aligns):
            ax.text(x, y, cell, ha=al, va="center", fontsize=8, color=INK,
                    transform=ax.transAxes)

    for y, lw in ((y_top, 1.1), (y_mid, 0.7), (y_bot, 1.1)):
        ax.plot([0, 1], [y, y], transform=ax.transAxes, color=INK, lw=lw,
                clip_on=False, zorder=5)
    ax.set_ylim(y_bot - 0.06, y_top + 0.03)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = HERE / f"tableS1.{ext}"
        fig.savefig(out, dpi=300)
        print(f"wrote {out}")
    plt.close(fig)


def main():
    recs = json.load(open(HERE / "paper_benchmark.json"))
    render_table(recs)
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(9.6, 3.0))

    for key, label, color in SERIES:
        r = decade_points([x for x in recs if x["example"] == key])
        n = np.array([x["n"] for x in r], float)
        t = np.array([x["gen_mean"] for x in r], float)
        tsd = np.array([x["gen_sd"] for x in r], float)
        m = np.array([x["peak_mb_mean"] for x in r], float)
        msd = np.array([x["peak_mb_sd"] for x in r], float)

        style = dict(color=color, lw=1.6, marker="o", ms=4.4, mfc=color,
                     mec="white", mew=0.7, capsize=2, elinewidth=0.9, zorder=3)
        axA.errorbar(n, t, yerr=tsd, label=label, **style)
        # Throughput SD propagated from the timing SD: d(N/T) = N/T^2 * dT
        axB.errorbar(n, n / t, yerr=n * tsd / t**2, **style)
        axC.errorbar(n, m, yerr=msd, **style)

    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlabel("sequences generated")
    axA.set_ylabel("generation time (s)")
    axA.set_title("A", loc="left", color=INK, fontweight="bold")

    axB.set_xscale("log"); axB.set_yscale("log")
    axB.set_xlabel("sequences generated")
    axB.set_ylabel("throughput (sequences s$^{-1}$)")
    axB.set_title("B", loc="left", color=INK, fontweight="bold")

    axC.set_xscale("log")
    axC.set_xlabel("sequences generated")
    axC.set_ylabel("peak memory (MB)")
    axC.set_title("C", loc="left", color=INK, fontweight="bold")
    axC.set_ylim(0, 540)

    for ax in (axA, axB, axC):
        ax.grid(True, which="major", axis="both", lw=0.4, color="#e6e5e0", zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    handles, labels = axA.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.05), handlelength=1.8, columnspacing=1.6)
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    for ext in ("pdf", "png"):
        out = HERE / f"figS_scaling.{ext}"
        fig.savefig(out, dpi=300)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
