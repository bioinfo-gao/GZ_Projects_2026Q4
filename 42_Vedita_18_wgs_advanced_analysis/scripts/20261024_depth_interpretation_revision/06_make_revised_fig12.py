#!/usr/bin/env python3
"""Render the corrected P42 Figure 12 in a fresh directory.

This figure explicitly separates the approximate WGS data-yield context from the
restricted candidate-comparison denominator.  It does not impose a retrospective
depth pass/fail gate or make an absence claim from comparator zero counts.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/gao/.codex/skills/scientific-figure-making/assets")
from pubfig import PALETTE, assert_no_text_collisions, finalize_figure, soften_for_fill  # noqa: E402


def finish(fig: plt.Figure, out_dir: Path) -> None:
    fig.tight_layout(pad=1.4)
    assert_no_text_collisions(fig)
    finalize_figure(fig, out_dir / "P42_Fig12_evidence_boundary_summary", formats=["png", "pdf"], dpi=300, close=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"Refusing to overwrite {args.out_dir}")
    args.out_dir.mkdir(parents=True)

    fig, ax = plt.subplots(figsize=(15.2, 7.0))
    ax.set_axis_off()
    ax.text(0.03, 0.94, "P42 Fig12. Interpretation boundary for clone/doxycycline comparison", ha="left", va="center", transform=ax.transAxes, fontsize=15, fontweight="bold")
    boxes = [
        (0.03, 0.58, 0.21, 0.24, "Client FISH\nQualitative chromosome-fragment\nabnormality reported\n(no raw image/table supplied)", PALETTE["highlight"]),
        (0.29, 0.58, 0.21, 0.24, "WGS data-yield context\nApproximately 25–30× for\nmost samples\n(separate from analysis denominator)", soften_for_fill(PALETTE["blue_main"], 0.55)),
        (0.55, 0.58, 0.21, 0.24, "Candidate-comparison context\nPrimary, non-duplicate,\nMAPQ≥20 evidence denominator\n(not a service-yield metric)", soften_for_fill(PALETTE["green_3"], 0.55)),
        (0.81, 0.58, 0.16, 0.24, "Mapping evidence\nC1–C4 remain\ncandidate observations\nfrom the current WGS", soften_for_fill(PALETTE["red_strong"], 0.55)),
    ]
    for x0, y0, width, height, label, color in boxes:
        ax.add_patch(plt.Rectangle((x0, y0), width, height, facecolor=color, edgecolor="black", linewidth=1.1, transform=ax.transAxes))
        ax.text(x0 + width / 2, y0 + height / 2, label, ha="center", va="center", transform=ax.transAxes, fontsize=9.5)
    for start in (0.24, 0.50, 0.76):
        ax.annotate("", xy=(start + 0.04, 0.70), xytext=(start, 0.70), xycoords=ax.transAxes,
                    arrowprops={"arrowstyle": "->", "lw": 1.4, "color": PALETTE["ink"]})
    ax.text(0.50, 0.34, "Clone/doxycycline specificity is not established from the current WGS comparison alone", ha="center", va="center", transform=ax.transAxes, fontsize=15, fontweight="bold", color=PALETTE["red_strong"])
    ax.text(0.50, 0.20, "C2 has internally consistent split-read structure; exact Manta/DELLY concordance is 0/20 and rescue did not change C1–C4 counts.", ha="center", va="center", transform=ax.transAxes, fontsize=10)
    ax.text(0.50, 0.10, "Comparator zeros are reported as not observed in the available evidence context; no retrospective numeric depth gate or direct allele-fraction cutoff is applied.", ha="center", va="center", transform=ax.transAxes, fontsize=9.5)
    finish(fig, args.out_dir)


if __name__ == "__main__":
    main()
