#!/usr/bin/env python3
"""Create a delivery-grade P41 §7.3 to P42 assay-follow-up schematic."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


COLORS = {
    "primary": "#FDB462", "second": "#B3DE69", "external": "#BEBADA",
    "line": "#003C30", "text": "#272727", "background": "#FFFFFF",
}


def read_tsv(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selection", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    args = ap.parse_args()
    rows = sorted(read_tsv(args.selection), key=lambda x: x["candidate_label"])
    assert [r["candidate_label"] for r in rows] == ["C1", "C2", "C3", "C4"]
    assert len({COLORS["primary"], COLORS["second"], COLORS["external"]}) == 3
    assert all(c.lower() not in {"#000000", "#272727"} for c in (COLORS["primary"], COLORS["second"], COLORS["external"]))

    fig = plt.figure(figsize=(15, 8), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.1, 1.0])
    ax = fig.add_subplot(gs[0])
    values = [int(r["raw_dual_endpoint_qname_count"]) for r in rows]
    fills = [COLORS["primary"] if r["candidate_label"] == "C2" else COLORS["second"] for r in rows]
    bars = ax.bar([r["candidate_label"] for r in rows], values, color=fills, edgecolor=COLORS["text"], linewidth=1.0)
    ax.set_ylabel("Raw dual-endpoint QNAMEs")
    ax.set_ylim(0, max(values) + 1.2)
    ax.set_title("P42 assay-design priorities derived from P41 Section 7.3", loc="left", weight="bold")
    for bar, row, val in zip(bars, rows, values):
        label = "P41 EV0623 carry-forward" if row["candidate_label"] == "C2" else "P42 second-tier"
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.08, f"{val}\n{label}", ha="center", va="bottom", fontsize=10, color=COLORS["text"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25, linewidth=0.7)
    ax.set_axisbelow(True)

    flow = fig.add_subplot(gs[1])
    flow.set_axis_off()
    boxes = [
        (0.02, "P41 Section 7.3\nEV0623: S05 candidate\n3 split reads / 0 pairs", COLORS["primary"]),
        (0.27, "P42 crosswalk\nEV0623 = C2 geometry\n(raw WGS evidence)", COLORS["primary"]),
        (0.52, "In-silico assay design\norientation-aware template\n+ whole-reference audit", COLORS["second"]),
        (0.77, "External validation pending\nindependent DNA, empirical PCR\nnew local WGS / raw FISH", COLORS["external"]),
    ]
    for x, text, color in boxes:
        flow.add_patch(FancyBboxPatch((x, 0.30), 0.19, 0.43, boxstyle="round,pad=0.02,rounding_size=0.02",
                                      linewidth=1.2, edgecolor=COLORS["text"], facecolor=color, transform=flow.transAxes))
        flow.text(x + 0.095, 0.515, text, ha="center", va="center", fontsize=10.2, color=COLORS["text"], transform=flow.transAxes)
    for x in (0.215, 0.465, 0.715):
        flow.annotate("", xy=(x + 0.04, 0.515), xytext=(x, 0.515), xycoords="axes fraction",
                      arrowprops={"arrowstyle": "->", "color": COLORS["line"], "lw": 1.8})
    flow.text(0.02, 0.08, "In-silico design does not validate a structural variant; a negative PCR/WGS result requires measured assay sensitivity.",
              transform=flow.transAxes, fontsize=10.2, color=COLORS["text"])
    fig.savefig(args.out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(args.out_prefix.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")


if __name__ == "__main__":
    main()
