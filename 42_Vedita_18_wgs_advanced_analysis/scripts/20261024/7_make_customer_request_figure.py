#!/usr/bin/env python3
"""Render the client-question completion figure from new P42 answer tables only."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker import NullLocator

sys.path.insert(0, "/home/gao/.codex/skills/scientific-figure-making/assets")
from pubfig import (  # noqa: E402
    FigureStyle,
    apply_publication_style,
    assert_fill_not_near_black,
    assert_label_ceiling,
    assert_no_palette_collision,
    assert_no_text_collisions,
    assert_no_yellow_text,
    assert_visible_on_white,
    legend_outside,
)


COLORS = {
    "timeline": "#80B1D3",
    "c1": "#BEBADA",
    "c2": "#FB8072",
    "c3": "#FCCDE5",
    "c4": "#BC80BD",
    "anchor_only": "#B3DE69",
    "no_exact": "#999999",
    "ink": "#272727",
    "line": "#003C30",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    args = parse_args()
    in_dir = args.input_dir.resolve()
    out = args.out_prefix.resolve()
    png = out.with_suffix(".png")
    pdf = out.with_suffix(".pdf")
    if png.exists() or pdf.exists():
        raise FileExistsError(f"Refusing to overwrite an existing figure: {out}")

    c2 = read_tsv(in_dir / "c2_junction_resolution_summary.tsv")
    rescue = read_tsv(in_dir / "read_rescue_candidate_impact.tsv")
    caller = read_tsv(in_dir / "candidate_manta_delly_concordance.tsv")
    if len(c2) != 1 or c2[0]["candidate_label"] != "C2":
        raise ValueError("Expected one C2 junction-resolution row")
    c2_row = c2[0]
    if (c2_row["independent_qname_count"], c2_row["split_read_end_count"],
            c2_row["query_overlap_bp_distribution"], c2_row["unmapped_query_gap_bp_distribution"]) != (
                "3", "5", "2bp:5", "0;0;0;0;0"):
        raise ValueError("Unexpected C2 base-pair resolution summary")

    rescue = sorted(rescue, key=lambda row: row["candidate_label"])
    if [row["candidate_label"] for row in rescue] != ["C1", "C2", "C3", "C4"]:
        raise ValueError("Expected C1-C4 rescue-impact rows")
    before = [int(row["raw_dual_endpoint_qname_count_before_rescue"]) for row in rescue]
    after = [int(row["raw_dual_endpoint_qname_count_after_rescue"]) for row in rescue]
    overlaps = [int(row["rescued_unmapped_mate_qname_overlap"]) for row in rescue]
    if before != after or any(overlaps):
        raise ValueError("Figure is specified only for the verified zero-impact C1-C4 rescue result")

    status_counts: dict[str, int] = {}
    manta_exact = delly_exact = both_exact = 0
    for row in caller:
        status = row["caller_concordance_status"]
        status_counts[status] = status_counts.get(status, 0) + 1
        manta_count = int(row["manta_exact_two_end_record_count"])
        delly_count = int(row["delly_exact_two_end_record_count"])
        manta_exact += manta_count
        delly_exact += delly_count
        both_exact += int(manta_count > 0 and delly_count > 0)
    if len(caller) != 20 or status_counts != {
        "anchor_only_or_nonconcordant": 5, "no_exact_Manta_or_DELLY_match_in_audited_sources": 15
    } or (manta_exact, delly_exact, both_exact) != (0, 0, 0):
        raise ValueError("Unexpected Manta/DELLY concordance aggregate")

    timeline_fills = [COLORS["timeline"]]
    candidate_fills = [COLORS[k] for k in ("c1", "c2", "c3", "c4")]
    caller_fills = [COLORS["anchor_only"], COLORS["no_exact"]]
    assert_no_palette_collision({
        "P41-to-P42 timeline": timeline_fills,
        "candidate support bars": candidate_fills,
        "caller-concordance bar": caller_fills,
    })
    all_fills = timeline_fills + candidate_fills + caller_fills
    assert_visible_on_white(all_fills, role="fill")
    assert_fill_not_near_black(all_fills)
    assert_visible_on_white([COLORS["ink"], COLORS["line"]], role="line")
    assert_no_yellow_text([COLORS["ink"], COLORS["line"]])

    apply_publication_style(FigureStyle(font_size=14, axes_linewidth=1.7))
    fig = plt.figure(figsize=(15.5, 9.3), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.08, 1.0], hspace=0.10, wspace=0.20)
    timeline = fig.add_subplot(gs[0, :], label="timeline")
    support = fig.add_subplot(gs[1, 0], label="support")
    callers = fig.add_subplot(gs[1, 1], label="callers")

    timeline.set_axis_off()
    # ``set_axis_off`` leaves default tick Text objects visible to geometry
    # inspection in some Matplotlib releases; hide both axis containers too.
    timeline.xaxis.set_visible(False)
    timeline.yaxis.set_visible(False)
    timeline.xaxis.set_major_locator(NullLocator())
    timeline.xaxis.set_minor_locator(NullLocator())
    timeline.yaxis.set_major_locator(NullLocator())
    timeline.yaxis.set_minor_locator(NullLocator())
    timeline.set_title("P42 answers the P41 Section 7.3 customer questions", loc="left", pad=8, weight="bold")
    boxes = [
        (0.035, "P41 Section 7.3\nEV0623 in S05\n3 split reads / 0 pairs"),
        (0.365, "P42 C2 crosswalk\nchr11:114,607,191 ↔\nchr21:20,648,242"),
        (0.695, "P42 read-end resolution\n3 distinct QNAMEs; 5 read ends\n2-bp AT overlap; 0-bp gap"),
    ]
    for x, text in boxes:
        patch = FancyBboxPatch(
            (x, 0.36), 0.255, 0.36, boxstyle="round,pad=0.025,rounding_size=0.025",
            linewidth=1.2, edgecolor=COLORS["ink"], facecolor=COLORS["timeline"],
            transform=timeline.transAxes,
        )
        timeline.add_patch(patch)
        timeline.text(x + 0.1275, 0.54, text, ha="center", va="center", fontsize=12.0,
                      color=COLORS["ink"], transform=timeline.transAxes)
    for x in (0.300, 0.630):
        timeline.annotate("", xy=(x + 0.052, 0.54), xytext=(x, 0.54), xycoords="axes fraction",
                          arrowprops={"arrowstyle": "->", "lw": 1.8, "color": COLORS["line"]})
    timeline.text(
        0.035, 0.13,
        "Interpretation boundary: the five split-read alignment representations are internally consistent with a 2-bp microhomology; "
        "repeat-context placement ambiguity and lack of independent validation remain.",
        transform=timeline.transAxes, fontsize=10.6, color=COLORS["ink"], wrap=True,
    )

    labels = [row["candidate_label"] for row in rescue]
    bars = support.bar(labels, after, color=candidate_fills, edgecolor=COLORS["ink"], linewidth=1.1)
    support.set_ylim(0, 3.85)
    support.set_yticks([0, 1, 2, 3])
    support.set_ylabel("Raw dual-endpoint QNAMEs")
    support.set_title("Read rescue: C1-C4 counts were unchanged\n(5 unmapped mates retained across the locus; 0 overlap C1-C4)",
                      loc="left", weight="bold", fontsize=13)
    support.grid(axis="y", alpha=0.22, linewidth=0.7)
    support.set_axisbelow(True)
    for bar, value in zip(bars, after):
        support.text(bar.get_x() + bar.get_width() / 2, value + 0.10,
                     f"{value}\nbefore = after", ha="center", va="bottom", fontsize=10.3, color=COLORS["ink"])

    anchor_n = status_counts["anchor_only_or_nonconcordant"]
    no_exact_n = status_counts["no_exact_Manta_or_DELLY_match_in_audited_sources"]
    callers.barh([0], [anchor_n], color=COLORS["anchor_only"], edgecolor=COLORS["ink"], linewidth=1.0,
                 label="Anchor-only / nonconcordant")
    callers.barh([0], [no_exact_n], left=[anchor_n], color=COLORS["no_exact"], edgecolor=COLORS["ink"], linewidth=1.0,
                 label="No exact Manta or DELLY match")
    callers.set_xlim(0, 20)
    callers.set_ylim(-0.7, 0.7)
    callers.set_xticks([0, 5, 10, 15, 20])
    callers.set_yticks([])
    callers.set_xlabel("C1-C4 × DLD1-sample comparisons (n = 20)\nAnchor-only: C4 DELLY chr11 INS records without the chr3 partner")
    callers.set_title("Caller concordance: exact two-end BND matches are 0/20\nManta 0/20; DELLY 0/20; both callers 0/20",
                      loc="left", weight="bold", fontsize=13)
    callers.text(anchor_n / 2, 0, "5\nanchor-only", ha="center", va="center", fontsize=11, color=COLORS["ink"])
    callers.text(anchor_n + no_exact_n / 2, 0, "15\nno exact match", ha="center", va="center", fontsize=11, color=COLORS["ink"])

    for axis in (timeline, support, callers):
        assert_label_ceiling(axis, max_annotations=6)
    assert_no_text_collisions(fig)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
