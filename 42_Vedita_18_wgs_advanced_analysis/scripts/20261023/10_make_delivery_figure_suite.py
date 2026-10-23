#!/usr/bin/env python3
"""Generate the P42 delivery figure suite from canonical audited tables.

All plots preserve the distinction between raw mapping evidence, caller fields,
and the pre-registered high-confidence gate.  No panel interprets a caller
field as an independent molecule count or treats lack of a call as absence.
"""
from __future__ import annotations

import argparse
import csv
import gzip
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
    PALETTE, FigureStyle, apply_publication_style, assert_fill_not_near_black,
    assert_no_palette_collision, assert_no_text_collisions, assert_shared_palette,
    assert_visible_on_white, finalize_figure, make_heatmap, panel_crops, soften_for_fill,
)


DLD1 = ["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl", "S04_DLD1_c3_dox", "S05_DLD1_c4_dox"]
CASES = {"S04_DLD1_c3_dox", "S05_DLD1_c4_dox"}
LINE_COLORS = {
    "DLD1-136": soften_for_fill(PALETTE["blue_main"], 0.55),
    "GM12878": soften_for_fill(PALETTE["green_3"], 0.55),
    "Mutu I": soften_for_fill(PALETTE["red_strong"], 0.55),
    "DG75": soften_for_fill(PALETTE["violet"], 0.55),
}


def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def numeric(value, default=0.0):
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def sid_short(sample_id):
    return sample_id.split("_", 1)[0]


def finish(fig, out: Path, stem: str):
    fig.tight_layout(pad=2.0)
    assert_no_text_collisions(fig)
    panel_crops(fig)
    finalize_figure(fig, out / stem, formats=["png", "pdf"], dpi=300, close=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.exists():
        raise SystemExit(f"Refusing to overwrite figure directory: {args.out_dir}")
    args.out_dir.mkdir(parents=True)
    root = args.project / "analysis/20261023_fullcohort_repeat_resolve"
    samples = read_csv(args.project / "scripts/0_sample_table.csv")
    order = [row["sample_id"] for row in samples]
    line_by_sample = {row["sample_id"]: row["cell_line"] for row in samples}
    colors = [LINE_COLORS[line_by_sample[sample]] for sample in order]
    assert_shared_palette({"cohort depth": LINE_COLORS, "raw clusters": LINE_COLORS})
    assert_visible_on_white(list(LINE_COLORS.values()), role="fill")
    assert_fill_not_near_black(list(LINE_COLORS.values()))

    depths = {row["sample_id"]: numeric(row["canonical_genome_mean_depth"]) for row in read_tsv(root / "03_coverage/canonical_genome_mean_depth.tsv")}
    normalized = read_tsv(root / "03_coverage/discovery_1kb_normalized.tsv")
    inventory = {row["sample_id"]: row for row in read_tsv(root / "01_full_window_inventory/inventory_audit_summary.tsv")}
    cluster_summary = {row["sample_id"]: row for row in read_tsv(root / "07_external_sv_dld1_normalization_v1/full_evidence_cluster_summary_by_sample.tsv")}
    assembly = read_tsv(root / "05_local_assembly/assemblies_v1/assembly_run_status.tsv")
    contigs = read_tsv(root / "05_local_assembly/contig_placements_v1/contig_placement_summary.tsv")
    rescue = read_tsv(root / "05_local_assembly/inputs_v1/unmapped_mate_rescue/unmapped_mate_recovery_status.tsv")
    external = read_tsv(root / "07_external_sv_dld1_normalization_v1/normalized_external_sv_events.tsv")
    gridss = read_tsv(root / "07_external_sv_dld1_normalization_v1/dld1_external_caller_matrix.tsv")
    fixed = read_tsv(root / "09_fixed_candidate_regenotyping_v1/fixed_candidate_dld1_matrix.tsv")
    fixed_summary = read_tsv(root / "09_fixed_candidate_regenotyping_v1/fixed_candidate_summary.tsv")
    with (root / "09_fixed_candidate_regenotyping_v1/audit.json").open() as handle:
        regenotype_audit = json.load(handle)

    # Fig 01 — 18-sample whole-genome depth and cohort design.
    fig, ax = plt.subplots(figsize=(15.6, 6.1))
    bars = ax.bar([sid_short(s) for s in order], [depths[s] for s in order], color=colors, edgecolor="black", linewidth=0.7)
    for bar, depth in zip(bars, [depths[s] for s in order]):
        ax.annotate(f"{depth:.1f}", (bar.get_x() + bar.get_width()/2, bar.get_height()), ha="center", va="bottom", fontsize=8, xytext=(0, 3), textcoords="offset points")
    ax.set_ylabel("Canonical genome mean depth (×)")
    ax.set_title("P42 Fig01. Cohort-wide standard coverage", loc="left", fontweight="bold")
    ax.set_ylim(0, max(depths.values()) * 1.18)
    handles = [plt.Rectangle((0, 0), 1, 1, color=color, ec="black") for color in LINE_COLORS.values()]
    ax.legend(handles, list(LINE_COLORS), ncol=4, frameon=False, title="Cell line", loc="upper left")
    finish(fig, args.out_dir, "P42_Fig01_cohort_coverage")

    # Fig 02 — DLD1 chr11 normalized depth trajectories, cases dashed.
    fig, ax = plt.subplots(figsize=(14.2, 6.4))
    dld1_style = {DLD1[0]: (PALETTE["blue_main"], "-"), DLD1[1]: (PALETTE["blue_secondary"], "-"), DLD1[2]: (PALETTE["teal"], "-"), DLD1[3]: (PALETTE["red_strong"], "--"), DLD1[4]: (PALETTE["red_2"], "--")}
    assert_no_palette_collision({"DLD1 control tracks": [dld1_style[s][0] for s in DLD1[:3]], "DLD1 +dox tracks": [dld1_style[s][0] for s in DLD1[3:]]})
    for sample in DLD1:
        rows = [row for row in normalized if row["sample_id"] == sample and row["track"] == "standard_primary_mapq20"]
        rows.sort(key=lambda row: numeric(row["start_1based"]))
        x = [numeric(row["start_1based"]) / 1e6 for row in rows]
        y = [numeric(row["depth_over_genome_mean"]) for row in rows]
        color, style = dld1_style[sample]
        ax.plot(x, y, color=color, linestyle=style, linewidth=1.5, label=(sid_short(sample) + (" +dox" if sample in CASES else " control")))
    ax.axhline(1, color=PALETTE["ink"], linewidth=0.8, linestyle=":")
    ax.set_xlabel("chr11 coordinate (Mb)")
    ax.set_ylabel("Depth / sample genome mean")
    ax.set_title("P42 Fig02. DLD1 chr11 discovery-window coverage", loc="left", fontweight="bold")
    ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.18), frameon=False)
    finish(fig, args.out_dir, "P42_Fig02_DLD1_chr11_depth")

    # Fig 03 — record class and mapping-quality audit; MAPQ0 is not stacked because it overlaps classes.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15.6, 9.0), sharex=True)
    class_colors = [soften_for_fill(PALETTE["blue_main"], 0.55), soften_for_fill(PALETTE["magenta"], 0.55), soften_for_fill(PALETTE["teal"], 0.55)]
    assert_visible_on_white(class_colors, role="fill")
    x = np.arange(len(order)); bottom = np.zeros(len(order))
    for field, label, color in zip(["records_primary", "records_secondary", "records_supplementary"], ["Primary", "Secondary", "Supplementary"], class_colors):
        values = np.array([numeric(inventory[s][field]) for s in order])
        ax1.bar(x, values, bottom=bottom, color=color, edgecolor="black", linewidth=0.35, label=label)
        bottom += values
    ax1.set_ylabel("Records in full window")
    ax1.set_title("A. Alignment classes retained in the raw ledger", loc="left", fontweight="bold")
    ax1.legend(ncol=3, frameon=False, loc="upper right")
    mapq0_pct = [100 * numeric(inventory[s]["records_mapq0"]) / max(1, numeric(inventory[s]["records_total"])) for s in order]
    ax2.bar(x, mapq0_pct, color=colors, edgecolor="black", linewidth=0.7)
    ax2.set_xticks(x, [sid_short(s) for s in order])
    ax2.set_ylabel("MAPQ 0 records (%)")
    ax2.set_title("B. MAPQ 0 is audited separately, not excluded from discovery", loc="left", fontweight="bold")
    finish(fig, args.out_dir, "P42_Fig03_alignment_inventory")

    # Fig 04 — full-evidence cluster audit, using the same lineage colors as Fig01.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15.6, 9.0), sharex=True)
    all_clusters = [numeric(cluster_summary[s]["cluster_count"]) for s in order]
    strong_clusters = [numeric(cluster_summary[s]["clusters_ge3_fragments"]) for s in order]
    for ax, values, title in [(ax1, all_clusters, "A. Full-evidence junction clusters"), (ax2, strong_clusters, "B. Clusters with ≥3 independent fragments")]:
        ax.bar(np.arange(len(order)), values, color=colors, edgecolor="black", linewidth=0.65)
        ax.set_ylabel("Cluster count")
        ax.set_title(title, loc="left", fontweight="bold")
    ax2.set_xticks(np.arange(len(order)), [sid_short(s) for s in order])
    finish(fig, args.out_dir, "P42_Fig04_full_evidence_clusters")

    # Fig 05 — all 59 assemblies and five rescued unmapped mates.
    by_sample_assembly = Counter(row["sample_id"] for row in assembly if row["status"] == "completed")
    templates = defaultdict(float)
    for row in assembly:
        templates[row["sample_id"]] += numeric(row["paired_templates"])
    rescue_by_sample = Counter(row["sample_id"] for row in rescue if row["recovery_status"] == "recovered_unmapped_sequence")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16.0, 6.4))
    ax1.bar(np.arange(len(order)), [by_sample_assembly[s] for s in order], color=colors, edgecolor="black", linewidth=0.7)
    ax1.set_xticks(np.arange(len(order)))
    ax1.tick_params(axis="x", labelbottom=False)
    ax1.set_xlabel("Samples in canonical order S01–S18")
    ax1.set_ylabel("Completed local assemblies")
    ax1.set_title("A. Local-assembly coverage (59/59 loci completed)", loc="left", fontweight="bold")
    ax2.scatter([templates[s] for s in order], [rescue_by_sample[s] for s in order], s=95, c=colors, edgecolor="black", linewidth=0.7)
    ax2.set_xlabel("Paired templates supplied to local assemblies")
    ax2.set_ylabel("Recovered unmapped mate read ends")
    ax2.set_title("B. Targeted sequence rescue (5 samples; 1 read end each)", loc="left", fontweight="bold")
    finish(fig, args.out_dir, "P42_Fig05_local_assembly_audit")

    # Fig 06 — all-placement contig ambiguity, no primary-only interpretation.
    contig_by_sample = defaultdict(lambda: Counter())
    for row in contigs:
        sample = row["contig_id"].split("_L", 1)[0]
        loci = int(row["distinct_loci"] or 0)
        contig_by_sample[sample]["Unmapped"] += int(loci == 0)
        contig_by_sample[sample]["Unique locus"] += int(loci == 1)
        contig_by_sample[sample]["Alternative loci"] += int(loci > 1)
    fig, ax = plt.subplots(figsize=(15.6, 6.2))
    contig_colors = [PALETTE["grey_soft"], soften_for_fill(PALETTE["green_3"], 0.55), soften_for_fill(PALETTE["red_2"], 0.55)]
    assert_no_palette_collision({"contig placement categories": contig_colors})
    bottom = np.zeros(len(order))
    for label, color in zip(["Unmapped", "Unique locus", "Alternative loci"], contig_colors):
        vals = np.array([contig_by_sample[s][label] for s in order])
        ax.bar(np.arange(len(order)), vals, bottom=bottom, color=color, edgecolor="black", linewidth=0.45, label=label)
        bottom += vals
    ax.set_xticks(np.arange(len(order)), [sid_short(s) for s in order])
    ax.set_ylabel("Assembled contigs")
    ax.set_title("P42 Fig06. All-placement local-assembly audit", loc="left", fontweight="bold")
    ax.legend(ncol=3, frameon=False, loc="upper right")
    finish(fig, args.out_dir, "P42_Fig06_contig_placement_ambiguity")

    # Fig 07 — external caller inventory, with event geometry retained separately.
    source_order = ["joint_18sample_targeted", "dld1_multisample_germline", "dld1_multisample_somatic"]
    source_labels = ["GRIDSS\nall 18", "SvABA\ngermline", "SvABA\nsomatic"]
    event_counts = [sum(row["source"] == source for row in external) for source in source_order]
    geometry = regenotype_audit["filter_impact"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.0))
    caller_colors = [soften_for_fill(PALETTE["blue_main"], 0.55), soften_for_fill(PALETTE["blue_secondary"], 0.55), soften_for_fill(PALETTE["teal"], 0.55)]
    bars = ax1.bar(source_labels, event_counts, color=caller_colors, edgecolor="black", linewidth=0.7)
    for bar, value in zip(bars, event_counts): ax1.annotate(str(value), (bar.get_x()+bar.get_width()/2, value), ha="center", va="bottom", xytext=(0, 4), textcoords="offset points")
    ax1.set_ylabel("Normalized caller BND events")
    ax1.set_title("A. External-caller inventory", loc="left", fontweight="bold")
    geom_labels = ["All caller\nevents", "Two-endpoint\ngeometry", "Anchor-only /\nunparsed"]
    geom_values = [geometry["candidate_events_total"], geometry["two_endpoint_bnd_events"], geometry["not_two_endpoint_evaluable_events"]]
    geom_colors = [soften_for_fill(PALETTE["blue_main"], 0.55), soften_for_fill(PALETTE["green_3"], 0.55), PALETTE["highlight"]]
    bars = ax2.bar(geom_labels, geom_values, color=geom_colors, edgecolor="black", linewidth=0.7)
    for bar, value in zip(bars, geom_values): ax2.annotate(str(value), (bar.get_x()+bar.get_width()/2, value), ha="center", va="bottom", xytext=(0, 4), textcoords="offset points")
    ax2.set_ylabel("Caller BND events")
    ax2.set_title("B. Fixed-candidate geometry audit", loc="left", fontweight="bold")
    finish(fig, args.out_dir, "P42_Fig07_external_caller_and_geometry")

    # Fig 08 — all GRIDSS event-by-sample caller support (no top-N truncation).
    grid_rows = [row for row in gridss if row["caller"] == "GRIDSS" and row["source"] == "joint_18sample_targeted"]
    event_ids = sorted({row["caller_event_id"] for row in grid_rows})
    support = defaultdict(dict)
    for row in grid_rows: support[row["caller_event_id"]][row["sample_id"]] = numeric(row["caller_support"])
    event_ids.sort(key=lambda event: (-sum(support[event][sample] for sample in DLD1), event))
    heat = np.array([[support[event][sample] for sample in DLD1] for event in event_ids])
    fig, ax = plt.subplots(figsize=(8.4, 13.2))
    make_heatmap(ax, heat, ["S01\nbaseline", "S02\ncontrol", "S03\ncontrol", "S04\n+dox", "S05\n+dox"], None, cmap="Blues", cbar_label="GRIDSS VF (caller field)", vmin=0, vmax=max(1, heat.max()))
    ax.set_ylabel(f"GRIDSS events (n={len(event_ids)})")
    ax.set_title("P42 Fig08. DLD1 GRIDSS support\n(all events; caller field, not template count)", loc="left", fontweight="bold", fontsize=14)
    finish(fig, args.out_dir, "P42_Fig08_DLD1_GRIDSS_support_matrix")

    # Fig 09 — every case-only raw candidate and all five samples, no ranking/truncation.
    case_only = [row for row in fixed_summary if row["classification"] == "case_only_raw_candidate_requires_sensitivity_and_assembly"]
    case_ids = [row["candidate_id"] for row in case_only]
    fixed_lookup = {(row["candidate_id"], row["sample_id"]): row for row in fixed}
    raw_heat = np.array([[numeric(fixed_lookup[(candidate, sample)]["raw_dual_endpoint_fragment_count"]) for sample in DLD1] for candidate in case_ids])
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    make_heatmap(ax, raw_heat, ["S01", "S02", "S03", "S04", "S05"], [f"C{i+1}" for i in range(len(case_ids))], cmap="Reds", cbar_label="Raw dual-endpoint QNAME fragments", annotate=True, fmt="{:.0f}", vmin=0, vmax=max(1, raw_heat.max()))
    ax.set_xlabel("DLD1 samples (S01–S03 control; S04–S05 +dox)")
    ax.set_title("P42 Fig09. All case-only raw-read candidates (n=4)", loc="left", fontweight="bold")
    finish(fig, args.out_dir, "P42_Fig09_case_only_raw_candidates")

    # Fig 10 — exact BND geometry for all four case-only candidates.
    fig, ax = plt.subplots(figsize=(14.8, 5.4))
    y = np.arange(len(case_only))
    for idx, row in enumerate(case_only):
        endpoints = row["fixed_endpoints"].split(";")
        chr11 = next(point for point in endpoints if point.startswith("chr11:"))
        partner = next(point for point in endpoints if point != chr11)
        anchor_mb = numeric(chr11.split(":")[1]) / 1e6
        ax.scatter(anchor_mb, y[idx], s=130, color=PALETTE["red_strong"], edgecolor="black", zorder=3)
        case_counts = row["case_dual_endpoint_fragments"].replace("S04_DLD1_c3_dox:", "S04=").replace(";S05_DLD1_c4_dox:", ", S05=")
        ax.text(0.53, y[idx], f"{partner}  |  {case_counts}", transform=ax.get_yaxis_transform(), va="center", fontsize=10)
    ax.set_yticks(y, [f"C{i+1}" for i in range(len(case_only))])
    ax.invert_yaxis()
    ax.set_xlabel("chr11 anchor coordinate (Mb)")
    ax.set_title("P42 Fig10. Case-only candidate geometries and raw support", loc="left", fontweight="bold")
    finish(fig, args.out_dir, "P42_Fig10_case_only_breakpoint_schematics")

    # Fig 11 — filter impact is quantified rather than treating unmatched rows as discarded evidence.
    raw_by_sample = regenotype_audit["filter_impact"]["input_raw_evidence_rows_by_sample"]
    matched_by_sample = Counter()
    with gzip.open(root / "09_fixed_candidate_regenotyping_v1/raw_evidence_match_audit.tsv.gz", "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            matched_by_sample[row["sample_id"]] += 1
    fig, ax = plt.subplots(figsize=(12.8, 6.1))
    sample_x = np.arange(len(DLD1)); incoming = np.array([raw_by_sample[s] for s in DLD1]); matched = np.array([matched_by_sample[s] for s in DLD1])
    ax.bar(sample_x, incoming - matched, color=PALETTE["grey_soft"], edgecolor="black", linewidth=0.5, label="Raw rows not linked to a fixed caller geometry")
    ax.bar(sample_x, matched, bottom=incoming-matched, color=soften_for_fill(PALETTE["blue_main"], 0.55), edgecolor="black", linewidth=0.5, label="Anchor-context matches")
    ax.set_xticks(sample_x, [sid_short(s) for s in DLD1])
    ax.set_ylabel("Ungated raw-evidence rows")
    ax.set_title("P42 Fig11. Fixed-candidate re-genotyping filter-impact audit", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="upper right")
    finish(fig, args.out_dir, "P42_Fig11_raw_evidence_filter_impact")

    # Fig 12 — concise evidence boundary.  The FISH statement is qualitative
    # provenance from the plan; no FISH image/table was supplied for quantitation.
    fig, ax = plt.subplots(figsize=(15.2, 6.7))
    ax.set_axis_off()
    boxes = [
        (0.03, 0.57, 0.21, 0.25, "Client FISH\nQualitative chromosome-fragment\nabnormality reported\n(no raw image/table supplied)", PALETTE["highlight"]),
        (0.29, 0.57, 0.21, 0.25, "Coverage\nS04/S05 depth pattern is\na screening observation\nnot a breakpoint", soften_for_fill(PALETTE["blue_main"], 0.55)),
        (0.55, 0.57, 0.21, 0.25, "External callers\n282 BND events; 94 two-endpoint\ngeometries; same-BAM\nalgorithmic evidence", soften_for_fill(PALETTE["green_3"], 0.55)),
        (0.81, 0.57, 0.16, 0.25, "Raw re-genotyping\n4 case-only patterns\nnone high confidence", soften_for_fill(PALETTE["red_strong"], 0.55)),
    ]
    for x0, y0, width, height, text, color in boxes:
        rect = plt.Rectangle((x0, y0), width, height, facecolor=color, edgecolor="black", linewidth=1.1, transform=ax.transAxes)
        ax.add_patch(rect); ax.text(x0+width/2, y0+height/2, text, ha="center", va="center", transform=ax.transAxes, fontsize=10)
    for start in [0.24, 0.50, 0.76]: ax.annotate("", xy=(start+0.04, 0.695), xytext=(start, 0.695), xycoords=ax.transAxes, arrowprops={"arrowstyle": "->", "lw": 1.4, "color": PALETTE["ink"]})
    ax.text(0.50, 0.26, "High-confidence Dox/clone-specific conclusion: NOT MET", ha="center", va="center", transform=ax.transAxes, fontsize=17, fontweight="bold", color=PALETTE["red_strong"])
    ax.text(0.50, 0.12, "Reasons: controls are below the ≥60× depth needed for a 3-fragment / 5% sensitivity context, and no sequence-level contig junction is demonstrated.", ha="center", va="center", transform=ax.transAxes, fontsize=10)
    finish(fig, args.out_dir, "P42_Fig12_evidence_boundary_summary")

    manifest = {
        "figure_count": 12,
        "source_tables": [
            "01_full_window_inventory/inventory_audit_summary.tsv", "03_coverage/canonical_genome_mean_depth.tsv",
            "03_coverage/discovery_1kb_normalized.tsv", "05_local_assembly/assemblies_v1/assembly_run_status.tsv",
            "05_local_assembly/contig_placements_v1/contig_placement_summary.tsv",
            "07_external_sv_dld1_normalization_v1/*.tsv", "09_fixed_candidate_regenotyping_v1/*.tsv",
        ],
        "case_only_candidate_count": len(case_only),
        "high_confidence_dox_candidate_count": 0,
        "interpretation_boundary": "All caller, raw-mapping and contig-placement evidence remains algorithmic/mapping evidence from the same WGS data; no panel establishes an independently verified unique breakpoint.",
        "fisher_note": "Client FISH is a qualitative reported observation; no source image/table was supplied for quantification in this figure suite.",
    }
    (args.out_dir / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
