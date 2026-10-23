#!/usr/bin/env python3
"""Re-genotype fixed external-caller BND candidates against ungated raw evidence.

This script deliberately keeps two evidence layers apart.  The fixed universe
is every stage-7 caller BND, while the raw universe is the complete previously
extracted chr11 evidence ledger, not caller QNAMEs.  Raw primary, secondary,
supplementary, duplicate, QC-fail and MAPQ0 records stay auditable.  QNAME is
used only to report a non-duplicated independent-fragment count for an event
within one sample; it never deletes the underlying alignment-context rows.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


DLD1_CONTROLS = ("S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl")
DLD1_CASES = ("S04_DLD1_c3_dox", "S05_DLD1_c4_dox")
DLD1_SAMPLES = DLD1_CONTROLS + DLD1_CASES
TOLERANCE_BP = 50
STANDARD_TRACK = "standard_primary_mapq20"
REMOTE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
LOCAL_RE = re.compile(r"^([^:;]+):(\d+):[^>]*>(.*)$")
CONTIG_ORIGIN_RE = re.compile(r"^(.+)_L\d+_(chr[^_]+)_(\d+)_(\d+)\|")


def write_tsv(path: Path, fields, rows) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def as_int(value: str | int | None) -> int | None:
    try:
        return int(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def as_float(value: str | float | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def endpoint_string(endpoint) -> str:
    return f"{endpoint[0]}:{endpoint[1]}"


def parse_bnd_event(row: dict) -> tuple[tuple[str, int] | None, tuple[str, int] | None, str]:
    """Return canonical two-endpoint BND geometry from stage-7 endpoint text.

    A linked VCF BND is represented twice; parsing both records avoids treating
    its mate representation as two candidates.  Single breakends are retained
    but cannot be tested as a two-endpoint event.
    """
    pairs = set()
    anchors = set()
    for token in row["endpoints"].split(";"):
        match = LOCAL_RE.match(token)
        if not match:
            continue
        local = (match.group(1), int(match.group(2)))
        anchors.add(local)
        remote_match = REMOTE_RE.search(match.group(3))
        if remote_match:
            remote = (remote_match.group(1), int(remote_match.group(2)))
            pairs.add(tuple(sorted((local, remote))))
    if len(pairs) == 1:
        first, second = next(iter(pairs))
        return first, second, "two_endpoint_bnd"
    chr11_anchors = sorted(endpoint for endpoint in anchors if endpoint[0] == "chr11")
    if chr11_anchors and not pairs:
        return chr11_anchors[0], None, "anchor_only_breakend"
    if len(pairs) > 1:
        return None, None, "multi_geometry_unparsed"
    return None, None, "unparsed_bnd"


def parse_sa_destinations(tag: str):
    for entry in tag.rstrip(";").split(";"):
        fields = entry.split(",")
        if len(fields) >= 2:
            pos = as_int(fields[1])
            if pos is not None:
                yield fields[0], pos


def row_anchor_positions(row: dict):
    if row["evidence_type"] == "softclip":
        pos = as_int(row["breakend_1based"])
        if pos is not None:
            yield row["chrom"], pos
        return
    for key in ("start_1based", "end_1based"):
        pos = as_int(row[key])
        if pos is not None:
            yield row["chrom"], pos


def near(first: tuple[str, int], second: tuple[str, int]) -> bool:
    return first[0] == second[0] and abs(first[1] - second[1]) <= TOLERANCE_BP


def context_matches(row: dict, candidate: dict) -> bool:
    return any(near(anchor, endpoint) for anchor in row_anchor_positions(row) for endpoint in candidate["endpoints"])


def pair_matches(local: tuple[str, int], remote: tuple[str, int], candidate: dict) -> bool:
    if len(candidate["endpoints"]) != 2:
        return False
    first, second = candidate["endpoints"]
    return (near(local, first) and near(remote, second)) or (near(local, second) and near(remote, first))


def exact_dual_mode(row: dict, candidate: dict) -> str:
    if len(candidate["endpoints"]) != 2:
        return ""
    local_positions = list(row_anchor_positions(row))
    if row["evidence_type"] == "sa_split":
        for local in local_positions:
            for destination in parse_sa_destinations(row["sa_tag"]):
                if pair_matches(local, destination, candidate):
                    return "dual_endpoint_sa"
    if row["evidence_type"] == "discordant_pair":
        mate_pos = as_int(row["mate_start_1based"])
        if mate_pos is not None:
            remote = (row["mate_chrom"], mate_pos)
            for local in local_positions:
                if pair_matches(local, remote, candidate):
                    return "dual_endpoint_discordant_pair"
    return ""


def standard_eligible(row: dict) -> bool:
    return (
        row["record_class"] == "primary"
        and row["is_duplicate"] == "0"
        and row["is_qcfail"] == "0"
        and as_int(row["mapq"]) is not None
        and as_int(row["mapq"]) >= 20
    )


def read_stage7_events(path: Path):
    candidates = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["svtypes"] != "BND":
                continue
            first, second, geometry_status = parse_bnd_event(row)
            endpoints = tuple(endpoint for endpoint in (first, second) if endpoint is not None)
            candidate_id = f"{row['caller']}|{row['source']}|{row['caller_event_id']}"
            candidates.append({
                "candidate_id": candidate_id,
                "caller": row["caller"],
                "source": row["source"],
                "caller_event_id": row["caller_event_id"],
                "filters": row["filters"],
                "vcf_record_count": row["vcf_record_count"],
                "endpoints": endpoints,
                "geometry_status": geometry_status,
                "caller_endpoint_text": row["endpoints"],
            })
    return candidates


def make_candidate_index(candidates):
    by_bin = defaultdict(set)
    for candidate in candidates:
        for chrom, pos in candidate["endpoints"]:
            if chrom != "chr11":
                continue
            for bin_id in range((pos - TOLERANCE_BP) // 100, (pos + TOLERANCE_BP) // 100 + 1):
                by_bin[bin_id].add(candidate["candidate_id"])
    return by_bin


def load_raw_evidence(raw_dir: Path):
    rows_by_sample = {}
    input_counts = {}
    for sample in DLD1_SAMPLES:
        path = raw_dir / sample / "junction_evidence_records.tsv.gz"
        if not path.is_file():
            raise SystemExit(f"Required raw evidence ledger missing: {path}")
        with gzip.open(path, "rt", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows_by_sample[sample] = rows
        input_counts[sample] = len(rows)
    return rows_by_sample, input_counts


def load_depths(path: Path):
    by_sample = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["track"] != STANDARD_TRACK:
                continue
            depth = as_float(row["mean_depth"])
            if depth is not None:
                by_sample[row["sample_id"]].append((as_int(row["start_1based"]), as_int(row["end_1based"]), depth))
    return by_sample


def depth_at(depths, sample: str, endpoints):
    values = []
    for chrom, position in endpoints:
        if chrom != "chr11":
            continue
        for start, end, value in depths.get(sample, []):
            if start is not None and end is not None and start <= position <= end:
                values.append(value)
                break
    return min(values) if values else None


def load_contig_context(path: Path):
    by_sample = defaultdict(list)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            origin = CONTIG_ORIGIN_RE.match(row["contig_id"])
            if not origin:
                continue
            sample, chrom, start, end = origin.groups()
            by_sample[sample].append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "contig_id": row["contig_id"],
                "distinct_loci": as_int(row["distinct_loci"]) or 0,
            })
    return by_sample


def assembly_context(contigs, sample: str, endpoints):
    matched = []
    for contig in contigs.get(sample, []):
        for chrom, position in endpoints:
            if chrom == contig["chrom"] and contig["start"] - TOLERANCE_BP <= position <= contig["end"] + TOLERANCE_BP:
                matched.append(contig)
                break
    return {
        "assembly_context_contig_count": len(matched),
        "assembly_context_has_alternative_locus": int(any(row["distinct_loci"] > 1 for row in matched)),
        "assembly_junction_status": "not_tested_by_locus_level_contig_summary",
    }


def summarize_template_state(state: dict):
    return {
        "raw_context_record_count": state["context_records"],
        "raw_context_fragment_count": int(bool(state["contexts"])),
        "raw_dual_endpoint_fragment_count": int(bool(state["dual_modes"])),
        "raw_standard_eligible_dual_fragment_count": int(state["standard_dual"]),
        "raw_dual_modes": ",".join(sorted(state["dual_modes"])),
        "raw_context_evidence_types": ",".join(sorted(state["evidence_types"])),
        "raw_context_record_classes": ",".join(sorted(state["record_classes"])),
        "raw_context_read_ends": ",".join(sorted(state["read_ends"])),
        "raw_context_max_mapq": max(state["mapqs"]) if state["mapqs"] else "",
    }


def classify(candidate: dict, rows: list[dict]):
    if candidate["geometry_status"] != "two_endpoint_bnd":
        return "not_evaluable_as_two_endpoint_bnd", "caller event has no uniquely parsed two-endpoint geometry"
    by_sample = {row["sample_id"]: row for row in rows}
    cases = [by_sample[sample] for sample in DLD1_CASES]
    controls = [by_sample[sample] for sample in DLD1_CONTROLS]
    case_dual = [int(row["raw_dual_endpoint_fragment_count"]) for row in cases]
    control_dual = [int(row["raw_dual_endpoint_fragment_count"]) for row in controls]
    case_standard = [int(row["raw_standard_eligible_dual_fragment_count"]) for row in cases]
    control_depth_ok = all((as_float(row["standard_mean_depth_at_endpoint"]) or 0) >= 60 for row in controls)
    if any(case_dual) and any(control_dual):
        return "shared_raw_dual_endpoint_evidence", "matching raw dual-endpoint fragments occur in at least one case and control"
    if any(case_dual) and not any(control_dual):
        if all(value >= 3 for value in case_standard) and control_depth_ok:
            return "dox_pattern_raw_candidate_high_confidence_gate_not_met", "raw count/depth pattern present, but no sequence-level assembly junction assessment"
        return "case_only_raw_candidate_requires_sensitivity_and_assembly", "case-only raw dual-endpoint evidence; control sensitivity and/or assembly gate is not met"
    if any(control_dual):
        return "control_only_raw_dual_endpoint_evidence", "matching raw dual-endpoint fragments occur only in controls"
    if any(int(row["raw_context_fragment_count"]) for row in rows):
        return "anchor_context_without_matching_dual_endpoint", "raw anchor context exists but no same-QNAME two-endpoint match within 50 bp"
    return "no_matching_raw_dual_endpoint_evidence", "no raw anchor-context row matched this fixed caller geometry"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage7-events", type=Path, required=True)
    parser.add_argument("--raw-evidence-dir", type=Path, required=True)
    parser.add_argument("--coverage", type=Path, required=True)
    parser.add_argument("--contigs", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=False)

    candidates = read_stage7_events(args.stage7_events)
    if len(candidates) != 282:
        raise SystemExit(f"Expected all 282 stage-7 BND caller events, observed {len(candidates)}")
    by_id = {candidate["candidate_id"]: candidate for candidate in candidates}
    candidate_index = make_candidate_index(candidates)
    raw_rows, raw_input_counts = load_raw_evidence(args.raw_evidence_dir)
    depths = load_depths(args.coverage)
    contigs = load_contig_context(args.contigs)

    states = defaultdict(lambda: {
        "context_records": 0, "contexts": set(), "dual_modes": set(), "standard_dual": False,
        "evidence_types": set(), "record_classes": set(), "read_ends": set(), "mapqs": [],
    })
    matched_record_rows = []
    match_counts = Counter()
    for sample, rows in raw_rows.items():
        for row_number, row in enumerate(rows, start=2):
            possible = set()
            for chrom, position in row_anchor_positions(row):
                if chrom == "chr11":
                    possible.update(candidate_index.get(position // 100, set()))
            for candidate_id in sorted(possible):
                candidate = by_id[candidate_id]
                if not context_matches(row, candidate):
                    continue
                mode = exact_dual_mode(row, candidate)
                key = (candidate_id, sample, row["qname"])
                state = states[key]
                state["context_records"] += 1
                state["contexts"].add(row["evidence_type"])
                state["evidence_types"].add(row["evidence_type"])
                state["record_classes"].add(row["record_class"])
                state["read_ends"].add(row["read_end"])
                if as_int(row["mapq"]) is not None:
                    state["mapqs"].append(as_int(row["mapq"]))
                if mode:
                    state["dual_modes"].add(mode)
                    state["standard_dual"] = state["standard_dual"] or standard_eligible(row)
                    match_counts[f"{sample}_dual_{mode}"] += 1
                else:
                    match_counts[f"{sample}_anchor_context_{row['evidence_type']}"] += 1
                matched_record_rows.append({
                    "candidate_id": candidate_id, "sample_id": sample, "raw_evidence_row_number": row_number,
                    "qname": row["qname"], "read_end": row["read_end"], "match_mode": mode or "anchor_context_only",
                    "evidence_type": row["evidence_type"], "record_class": row["record_class"], "flag": row["flag"],
                    "is_duplicate": row["is_duplicate"], "is_qcfail": row["is_qcfail"], "mapq": row["mapq"],
                    "chrom": row["chrom"], "start_1based": row["start_1based"], "end_1based": row["end_1based"],
                    "breakend_1based": row["breakend_1based"], "mate_chrom": row["mate_chrom"],
                    "mate_start_1based": row["mate_start_1based"], "sa_tag": row["sa_tag"],
                })

    qname_ledger_rows = []
    matrix_rows = []
    matrix_by_candidate = defaultdict(list)
    for candidate in candidates:
        endpoints_text = ";".join(endpoint_string(endpoint) for endpoint in candidate["endpoints"])
        for sample in DLD1_SAMPLES:
            qname_states = [
                (qname, state) for (candidate_id, state_sample, qname), state in states.items()
                if candidate_id == candidate["candidate_id"] and state_sample == sample
            ]
            for qname, state in sorted(qname_states):
                qname_ledger_rows.append({
                    "candidate_id": candidate["candidate_id"], "sample_id": sample, "qname": qname,
                    **summarize_template_state(state),
                })
            context_qnames = [state for _, state in qname_states if state["contexts"]]
            dual_qnames = [state for _, state in qname_states if state["dual_modes"]]
            standard_dual_qnames = [state for _, state in qname_states if state["dual_modes"] and state["standard_dual"]]
            raw_classes = sorted({klass for state in context_qnames for klass in state["record_classes"]})
            raw_modes = sorted({mode for state in dual_qnames for mode in state["dual_modes"]})
            depth = depth_at(depths, sample, candidate["endpoints"])
            assembly = assembly_context(contigs, sample, candidate["endpoints"])
            row = {
                "candidate_id": candidate["candidate_id"], "caller": candidate["caller"], "source": candidate["source"],
                "caller_event_id": candidate["caller_event_id"], "filters": candidate["filters"],
                "geometry_status": candidate["geometry_status"], "fixed_endpoints": endpoints_text,
                "dld1_role": "case" if sample in DLD1_CASES else "control", "sample_id": sample,
                "raw_context_fragment_count": len(context_qnames),
                "raw_dual_endpoint_fragment_count": len(dual_qnames),
                "raw_standard_eligible_dual_fragment_count": len(standard_dual_qnames),
                "raw_dual_modes": ",".join(raw_modes), "raw_context_record_classes": ",".join(raw_classes),
                "standard_mean_depth_at_endpoint": "" if depth is None else f"{depth:.6f}",
                "three_fragment_fraction_context": "" if not depth else f"{3.0 / depth:.6f}",
                **assembly,
            }
            matrix_rows.append(row)
            matrix_by_candidate[candidate["candidate_id"]].append(row)

    summary_rows = []
    classification_counts = Counter()
    for candidate in candidates:
        rows = matrix_by_candidate[candidate["candidate_id"]]
        classification, rationale = classify(candidate, rows)
        classification_counts[classification] += 1
        by_sample = {row["sample_id"]: row for row in rows}
        summary_rows.append({
            "candidate_id": candidate["candidate_id"], "caller": candidate["caller"], "source": candidate["source"],
            "caller_event_id": candidate["caller_event_id"], "filters": candidate["filters"],
            "geometry_status": candidate["geometry_status"],
            "fixed_endpoints": ";".join(endpoint_string(endpoint) for endpoint in candidate["endpoints"]),
            "case_dual_endpoint_fragments": ";".join(
                f"{sample}:{by_sample[sample]['raw_dual_endpoint_fragment_count']}" for sample in DLD1_CASES
            ),
            "control_dual_endpoint_fragments": ";".join(
                f"{sample}:{by_sample[sample]['raw_dual_endpoint_fragment_count']}" for sample in DLD1_CONTROLS
            ),
            "case_standard_eligible_dual_fragments": ";".join(
                f"{sample}:{by_sample[sample]['raw_standard_eligible_dual_fragment_count']}" for sample in DLD1_CASES
            ),
            "control_standard_mean_depth": ";".join(
                f"{sample}:{by_sample[sample]['standard_mean_depth_at_endpoint'] or 'not_available'}" for sample in DLD1_CONTROLS
            ),
            "assembly_junction_gate": "not_assessed_by_locus_level_contig_summary",
            "classification": classification, "classification_rationale": rationale,
            "interpretation_boundary": "raw QNAME support is mapping evidence from the same BAM used by callers; it is not independent biological validation or a unique-locus proof",
        })

    candidate_fields = [
        "candidate_id", "caller", "source", "caller_event_id", "filters", "geometry_status", "fixed_endpoints",
        "dld1_role", "sample_id", "raw_context_fragment_count", "raw_dual_endpoint_fragment_count",
        "raw_standard_eligible_dual_fragment_count", "raw_dual_modes", "raw_context_record_classes",
        "standard_mean_depth_at_endpoint", "three_fragment_fraction_context", "assembly_context_contig_count",
        "assembly_context_has_alternative_locus", "assembly_junction_status",
    ]
    qname_fields = [
        "candidate_id", "sample_id", "qname", "raw_context_record_count", "raw_context_fragment_count",
        "raw_dual_endpoint_fragment_count", "raw_standard_eligible_dual_fragment_count", "raw_dual_modes",
        "raw_context_evidence_types", "raw_context_record_classes", "raw_context_read_ends", "raw_context_max_mapq",
    ]
    match_fields = [
        "candidate_id", "sample_id", "raw_evidence_row_number", "qname", "read_end", "match_mode", "evidence_type",
        "record_class", "flag", "is_duplicate", "is_qcfail", "mapq", "chrom", "start_1based", "end_1based",
        "breakend_1based", "mate_chrom", "mate_start_1based", "sa_tag",
    ]
    summary_fields = [
        "candidate_id", "caller", "source", "caller_event_id", "filters", "geometry_status", "fixed_endpoints",
        "case_dual_endpoint_fragments", "control_dual_endpoint_fragments", "case_standard_eligible_dual_fragments",
        "control_standard_mean_depth", "assembly_junction_gate", "classification", "classification_rationale",
        "interpretation_boundary",
    ]
    write_tsv(args.out_dir / "fixed_candidate_dld1_matrix.tsv", candidate_fields, matrix_rows)
    write_tsv(args.out_dir / "raw_qname_evidence_ledger.tsv", qname_fields, qname_ledger_rows)
    with gzip.open(args.out_dir / "raw_evidence_match_audit.tsv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=match_fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(matched_record_rows)
    write_tsv(args.out_dir / "fixed_candidate_summary.tsv", summary_fields, summary_rows)
    audit = {
        "script": "9_regenotype_fixed_candidates.py",
        "research_question": "For every stage-7 external-caller BND, quantify DLD1 raw dual-endpoint mapping evidence using the ungated raw ledger.",
        "candidate_universe": "All 282 normalized stage-7 BND caller events; this does not replace or close the ungated full-evidence discovery universe.",
        "raw_input_universe": "All previously extracted chr11 evidence rows for S01-S05, including primary/secondary/supplementary, duplicate, QC-fail and MAPQ0 context.",
        "independent_fragment_unit": "unique (candidate_id, sample_id, QNAME); multiple records/read ends from a QNAME are retained in raw_evidence_match_audit.tsv.gz and collapsed only for the reported fragment count.",
        "two_endpoint_rule": "SA destination or discordant mate must match the two caller endpoints in either order within 50 bp; a soft clip or supplementary anchor alone is context, not dual-endpoint support.",
        "standard_subset_rule": "Primary, MAPQ>=20, non-duplicate, non-QC-fail matched record; full raw evidence remains separately reported.",
        "assembly_rule": "Locus-overlapping contigs are context only. No sequence-level contig junction is inferred from the placement summary.",
        "filter_impact": {
            "input_raw_evidence_rows_by_sample": raw_input_counts,
            "input_raw_evidence_rows_total": sum(raw_input_counts.values()),
            "matched_anchor_context_rows": len(matched_record_rows),
            "raw_rows_not_linked_to_fixed_caller_geometry": sum(raw_input_counts.values()) - len(matched_record_rows),
            "candidate_events_total": len(candidates),
            "two_endpoint_bnd_events": sum(row["geometry_status"] == "two_endpoint_bnd" for row in candidates),
            "not_two_endpoint_evaluable_events": sum(row["geometry_status"] != "two_endpoint_bnd" for row in candidates),
            "candidate_sample_matrix_rows": len(matrix_rows),
            "matched_qname_rows": len(qname_ledger_rows),
        },
        "match_counts_by_sample_and_mode": dict(sorted(match_counts.items())),
        "classification_counts": dict(sorted(classification_counts.items())),
        "interpretation_limit": "Re-genotyping and external calls use the same BAM data. Results evaluate mapping/caller concordance, not independent biological confirmation or a unique genomic locus.",
    }
    (args.out_dir / "audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    (args.out_dir / ".complete").write_text("completed\n")


if __name__ == "__main__":
    main()
