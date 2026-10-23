#!/usr/bin/env python3
"""Create P42 external-validation audit figures from stage-7 canonical tables.

Every panel is descriptive.  In particular, caller support is not labelled as
an independent-molecule count or a clone-specific biological event.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/home/gao/.codex/skills/scientific-figure-making/assets")
from pubfig import (  # noqa: E402
    PALETTE,
    FigureStyle,
    apply_publication_style,
    assert_fill_not_near_black,
    assert_label_ceiling,
    assert_no_palette_collision,
    assert_no_text_collisions,
    assert_shared_palette,
    assert_visible_on_white,
    finalize_figure,
    legend_outside,
    make_heatmap,
    panel_crops,
    soften_for_fill,
)


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_manifest(path: Path, payload):
    def json_default(value):
        if hasattr(value, "item"):
            return value.item()
        raise TypeError(f"Unsupported manifest value: {type(value).__name__}")
    path.write_text(json.dumps(payload, indent=2, default=json_default) + "\n")


def bar_labels(ax, bars, values):
    for bar, value in zip(bars, values):
        ax.annotate(str(value), (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                    ha="center", va="bottom", xytext=(0, 4), textcoords="offset points",
                    fontsize=13, color=PALETTE["ink"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage7", type=Path, required=True)
    parser.add_argument("--sample-table", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing figure directory: {args.out_dir}")
    args.out_dir.mkdir(parents=True)

    events = read_tsv(args.stage7 / "normalized_external_sv_events.tsv")
    summary = read_tsv(args.stage7 / "dld1_external_caller_candidate_summary.tsv")
    matrix = read_tsv(args.stage7 / "dld1_external_caller_matrix.tsv")
    raw_summary = read_tsv(args.stage7 / "full_evidence_cluster_summary_by_sample.tsv")
    with (args.stage7 / "audit.json").open() as handle:
        audit = json.load(handle)
    sample_rows = read_csv(args.sample_table)

    apply_publication_style(FigureStyle(font_size=16, axes_linewidth=2.2))

    # Figure 1: source counts and DLD1 caller-support patterns.  Its panels
    # encode different variables and intentionally use disjoint fill families.
    source_order = ["joint_18sample_targeted", "dld1_multisample_germline", "dld1_multisample_somatic"]
    source_labels = ["GRIDSS\nall 18", "SvABA\ngermline", "SvABA\nsomatic"]
    source_counts = [sum(row["source"] == source for row in events) for source in source_order]
    source_colors = [soften_for_fill(PALETTE["blue_main"], 0.55),
                     soften_for_fill(PALETTE["blue_secondary"], 0.55),
                     soften_for_fill(PALETTE["teal"], 0.55)]
    pattern_order = ["case_and_control_observed", "case_only_observed", "control_only_observed", "not_called_in_dld1_samples"]
    pattern_labels = ["Case + control", "Case only", "Control only", "No DLD1 support"]
    pattern_colors = [soften_for_fill(PALETTE["red_strong"], 0.55),
                      soften_for_fill(PALETTE["red_2"], 0.62),
                      soften_for_fill(PALETTE["magenta"], 0.52), PALETTE["highlight"]]
    assert_no_palette_collision({"Fig1 source fills": source_colors, "Fig1 support-pattern fills": pattern_colors})
    assert_visible_on_white(source_colors + pattern_colors, role="fill")
    assert_fill_not_near_black(source_colors + pattern_colors)
    fig, (ax1, ax2, ax_legend) = plt.subplots(1, 3, figsize=(18.5, 6.2),
                                               gridspec_kw={"width_ratios": [0.95, 1.25, 0.72]})
    bars = ax1.bar(source_labels, source_counts, color=source_colors, edgecolor="black", linewidth=1.1)
    bar_labels(ax1, bars, source_counts)
    ax1.set_ylabel("Normalized caller events")
    ax1.set_title("A. External-caller event inventory", loc="left", fontweight="bold")
    ax1.set_ylim(0, max(source_counts) * 1.18 + 1)
    ax1.tick_params(axis="x", labelsize=11)
    assert_label_ceiling(ax1, max_annotations=6)

    x = np.arange(len(source_order))
    bottoms = np.zeros(len(source_order), dtype=int)
    totals = []
    for pattern, label, color in zip(pattern_order, pattern_labels, pattern_colors):
        values = np.array([sum(row["source"] == source and row["caller_support_pattern"] == pattern for row in summary)
                           for source in source_order])
        ax2.bar(x, values, bottom=bottoms, color=color, edgecolor="black", linewidth=0.7, label=label)
        bottoms += values
    totals = list(bottoms)
    for xpos, total in zip(x, totals):
        ax2.annotate(str(total), (xpos, total), ha="center", va="bottom", xytext=(0, 4),
                     textcoords="offset points", fontsize=13, color=PALETTE["ink"])
    ax2.set_xticks(x, source_labels)
    ax2.set_ylabel("Caller-defined events")
    ax2.set_title("B. DLD1 caller-support pattern", loc="left", fontweight="bold")
    ax2.set_ylim(0, max(totals) * 1.18 + 1)
    handles, labels = ax2.get_legend_handles_labels()
    ax_legend.set_axis_off()
    ax_legend.legend(handles, labels, loc="center", frameon=False, title="DLD1 caller-support pattern")
    assert_label_ceiling(ax2, max_annotations=6)
    assert_no_text_collisions(fig)
    panel_crops(fig)
    finalize_figure(fig, args.out_dir / "P42_Fig1_external_caller_audit", formats=["png", "pdf"], dpi=300, close=True)

    # Figure 2: every GRIDSS event in the stage-7 DLD1 matrix; no top-N
    # truncation.  Rows are sorted only for visual readability.
    gridss = [row for row in matrix if row["caller"] == "GRIDSS" and row["source"] == "joint_18sample_targeted"]
    event_ids = sorted({row["caller_event_id"] for row in gridss})
    sample_order = ["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl", "S04_DLD1_c3_dox", "S05_DLD1_c4_dox"]
    support = defaultdict(dict)
    for row in gridss:
        support[row["caller_event_id"]][row["sample_id"]] = float(row["caller_support"] or 0)
    event_ids.sort(key=lambda event: (-sum(support[event].get(sample, 0) for sample in sample_order), event))
    heat = np.array([[support[event].get(sample, 0) for sample in sample_order] for event in event_ids], dtype=float)
    fig, ax = plt.subplots(figsize=(8.2, 13.0))
    make_heatmap(ax, heat, x_labels=["S01\nbaseline", "S02\ncontrol", "S03\ncontrol", "S04\n+dox", "S05\n+dox"],
                 y_labels=None, cmap="Blues", cbar_label="GRIDSS VF (allocated caller support)", vmin=0, vmax=max(1, heat.max()))
    ax.set_ylabel(f"GRIDSS events (n={len(event_ids)})")
    ax.set_title("P42 Fig2. DLD1 GRIDSS support\n(all 278 events; no top-N truncation)",
                 loc="left", fontweight="bold", fontsize=14)
    assert_no_text_collisions(fig)
    panel_crops(fig)
    finalize_figure(fig, args.out_dir / "P42_Fig2_DLD1_GRIDSS_support_matrix", formats=["png", "pdf"], dpi=300, close=True)

    # Figure 3: raw evidence has a distinct evidence unit from the caller
    # matrix.  The same cell-line palette is deliberately reused in both
    # panels because its semantics are identical (sample lineage).
    line_by_sample = {row["sample_id"]: row["cell_line"] for row in sample_rows}
    group_colors = {
        "DLD1-136": soften_for_fill(PALETTE["blue_main"], 0.55),
        "GM12878": soften_for_fill(PALETTE["green_3"], 0.55),
        "Mutu I": soften_for_fill(PALETTE["red_strong"], 0.55),
        "DG75": soften_for_fill(PALETTE["violet"], 0.55),
    }
    assert_shared_palette({"Fig3 all clusters": group_colors, "Fig3 clusters >=3 fragments": group_colors})
    assert_visible_on_white(list(group_colors.values()), role="fill")
    assert_fill_not_near_black(list(group_colors.values()))
    raw_summary.sort(key=lambda row: int(row["sample_id"].split("_", 1)[0][1:]))
    labels = [row["sample_id"].split("_", 1)[0] for row in raw_summary]
    colors = [group_colors[line_by_sample[row["sample_id"]]] for row in raw_summary]
    all_clusters = [int(row["cluster_count"]) for row in raw_summary]
    strong_clusters = [int(row["clusters_ge3_fragments"]) for row in raw_summary]
    fig, (ax1, ax2, ax_legend) = plt.subplots(3, 1, figsize=(15.5, 10.2), sharex=True,
                                               gridspec_kw={"height_ratios": [1, 1, 0.17]})
    for ax, values, title in [
        (ax1, all_clusters, "A. Full-evidence junction clusters"),
        (ax2, strong_clusters, "B. Clusters with at least 3 independent fragments"),
    ]:
        ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.7)
        ax.set_ylabel("Cluster count")
        ax.set_title(title, loc="left", fontweight="bold")
    # ``sharex`` suppresses the upper labels by design, but categorical bars
    # do not guarantee visible labels on the bottom axes.  Set them explicitly
    # so every bar remains traceable to its sample.
    ax2.set_xticks(np.arange(len(labels)))
    ax2.set_xticklabels(labels, fontsize=12)
    ax2.tick_params(axis="x", which="both", labelbottom=True)
    for tick_label in ax2.get_xticklabels():
        tick_label.set_visible(True)
    # The explicit S01–S18 tick labels already identify this axis; omit a
    # redundant axis title so it cannot collide with the cell-line legend.
    handles = [plt.Rectangle((0, 0), 1, 1, color=color, ec="black") for color in group_colors.values()]
    ax_legend.set_axis_off()
    ax_legend.legend(handles, list(group_colors), ncol=4, loc="center", frameon=False, title="Cell line")
    if [tick.get_text() for tick in ax2.get_xticklabels()] != labels or not all(
        tick.get_visible() for tick in ax2.get_xticklabels()
    ):
        raise AssertionError("Figure 3 bottom-axis sample labels are absent or hidden")
    assert_no_text_collisions(fig)
    panel_crops(fig)
    finalize_figure(fig, args.out_dir / "P42_Fig3_full_evidence_cluster_audit", formats=["png", "pdf"], dpi=300, close=True)

    write_manifest(args.out_dir / "figure_manifest.json", {
        "source_stage": str(args.stage7),
        "figure_1": {"event_counts": dict(zip(source_order, source_counts)), "dld1_pattern_totals": dict(zip(source_order, totals))},
        "figure_2": {"gridss_events_rendered": len(event_ids), "dld1_samples_rendered": sample_order},
        "figure_3": {"samples_rendered": labels, "raw_cluster_table_rows": len(raw_summary)},
        "interpretation_boundary": "Figures report caller fields and raw-cluster counts. They are not final breakpoint or clone-specific classifications.",
        "contig_all_placement_audit": audit["contig_all_placement_summary"],
    })


if __name__ == "__main__":
    main()
