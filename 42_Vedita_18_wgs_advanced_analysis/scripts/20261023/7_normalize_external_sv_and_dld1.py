#!/usr/bin/env python3
"""Normalize P42 external caller outputs without promoting calls to biological truth.

Evidence model
--------------
Research question: are caller-defined structural candidates preferentially
observed in DLD1 +dox cases (S04/S05) compared with S01/S02/S03?

Observation unit: one caller event, represented by all of its linked VCF
breakend records.  Per-sample caller fields are reduced by maximum across
mate records only to avoid double-counting a two-record BND representation.
This is *not* an independent-molecule count.

Candidate universe: all target-related GRIDSS events and all target-related
SvABA SV events from the documented v3/v5 runs.  The code does not close the
raw-evidence universe: full-window clusters and assembly all-placement data
are summarized separately and remain required for final classification.

Allowed conclusion: caller-support pattern only.  A sample with no support in
one caller is written ``not_called_by_this_caller``, never "event absent".
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import pysam


DLD1_CASES = ("S04_DLD1_c3_dox", "S05_DLD1_c4_dox")
DLD1_CONTROLS = ("S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl")
DLD1_SAMPLES = DLD1_CONTROLS + DLD1_CASES
NUMERIC_GRIDSS = ("VF", "SR", "RP", "REF", "AF", "QUAL")
NUMERIC_SVABA = ("AD", "DP", "SR", "DR", "LO", "LR", "GQ")


def scalar(value):
    if value is None:
        return ""
    if isinstance(value, tuple):
        return ",".join("" if x is None else str(x) for x in value)
    return str(value)


def float_or_none(value):
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def sample_id(vcf_sample_name: str) -> str:
    """GRIDSS uses sample IDs; SvABA v5 uses complete BAM paths."""
    name = Path(vcf_sample_name).name
    return name.split(".markdup", 1)[0]


def vcf_event_id(record, caller: str) -> str:
    if caller == "GRIDSS":
        event = record.info.get("EVENT")
        if event:
            return str(event)
    mate_ids = record.info.get("MATEID") or ()
    if isinstance(mate_ids, str):
        mate_ids = (mate_ids,)
    labels = [str(record.id)] + [str(x) for x in mate_ids]
    # SvABA linked breakend IDs commonly end in :1 / :2.  Preserve the full
    # lexical group when another caller uses a different ID convention.
    roots = [label.rsplit(":", 1)[0] if label.rsplit(":", 1)[-1].isdigit() else label for label in labels]
    return min(roots)


def endpoint(record) -> str:
    alt = record.alts[0] if record.alts else "."
    return f"{record.contig}:{record.pos}:{record.ref}>{alt}"


def target_related(record) -> bool:
    alts = "|".join(record.alts or ())
    return record.contig == "chr11" or "chr11:" in alts or "]chr11" in alts or "[chr11" in alts


def numeric_max(values):
    parsed = [float_or_none(v) for v in values]
    parsed = [v for v in parsed if v is not None and not math.isnan(v)]
    return max(parsed) if parsed else ""


def parse_events(path: Path, caller: str, source: str):
    vf = pysam.VariantFile(str(path))
    sample_map = {name: sample_id(name) for name in vf.header.samples}
    groups = defaultdict(list)
    input_records = 0
    for record in vf:
        input_records += 1
        svtype = str(record.info.get("SVTYPE", ""))
        if svtype not in {"BND", "DEL", "DUP", "INV"} or not target_related(record):
            continue
        groups[vcf_event_id(record, caller)].append(record.copy())
    vf.close()
    return groups, sample_map, input_records


def summarize_event(event_id, records, caller, source, sample_map):
    svtypes = sorted({str(rec.info.get("SVTYPE", "")) for rec in records})
    filters = sorted({flt for rec in records for flt in rec.filter.keys()}) or ["."]
    endpoints = sorted({endpoint(rec) for rec in records})
    info = {
        "caller": caller,
        "source": source,
        "caller_event_id": event_id,
        "svtypes": ",".join(svtypes),
        "filters": ";".join(filters),
        "vcf_record_count": len(records),
        "endpoints": ";".join(endpoints),
    }
    per_sample = {}
    for vcf_name, sid in sample_map.items():
        rows = [rec.samples[vcf_name] for rec in records]
        if caller == "GRIDSS":
            metrics = {field: numeric_max([row.get(field) for row in rows]) for field in NUMERIC_GRIDSS}
            support = float_or_none(metrics["VF"])
            per_sample[sid] = {
                **metrics,
                "GT": "",
                "caller_support": support if support is not None else "",
                "caller_support_field": "GRIDSS_VF_allocated_fragment_count",
                "call_status": "caller_support_observed" if support and support > 0 else "not_called_by_this_caller",
                "detection_limit_3_support_fraction": "",
            }
        else:
            metrics = {field: numeric_max([row.get(field) for row in rows]) for field in NUMERIC_SVABA}
            gt = "/".join("." if x is None else str(x) for x in rows[0].get("GT", ())) if rows else ""
            support = float_or_none(metrics["AD"])
            dp = float_or_none(metrics["DP"])
            detection_limit = 3.0 / dp if dp and dp > 0 else ""
            per_sample[sid] = {
                **metrics,
                "GT": gt,
                "caller_support": support if support is not None else "",
                "caller_support_field": "SvABA_AD_caller_support_count",
                "call_status": "caller_support_observed" if support and support > 0 else "not_called_by_this_caller",
                "detection_limit_3_support_fraction": detection_limit,
            }
    return info, per_sample


def read_cluster_summary(path: Path):
    by_sample = defaultdict(lambda: {"cluster_count": 0, "clusters_ge3_fragments": 0, "max_fragment_count": 0})
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = row["sample_id"]
            fragments = int(row["independent_fragment_count"])
            by_sample[sample]["cluster_count"] += 1
            by_sample[sample]["clusters_ge3_fragments"] += int(fragments >= 3)
            by_sample[sample]["max_fragment_count"] = max(by_sample[sample]["max_fragment_count"], fragments)
    return by_sample


def read_contig_summary(path: Path):
    rows = list(csv.DictReader(path.open(newline=""), delimiter="\t"))
    return {
        "contigs": len(rows),
        "contigs_with_alternative_loci": sum(int(row["distinct_loci"]) > 1 for row in rows),
        "contigs_with_primary_mapq60": sum(row["primary_mapq"] == "60" for row in rows),
    }


def write_tsv(path: Path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gridss-vcf", type=Path, required=True)
    parser.add_argument("--svaba-germline-vcf", type=Path, required=True)
    parser.add_argument("--svaba-somatic-vcf", type=Path, required=True)
    parser.add_argument("--cluster-tsv", type=Path, required=True)
    parser.add_argument("--contig-summary", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=False)

    input_specs = [
        (args.gridss_vcf, "GRIDSS", "joint_18sample_targeted"),
        (args.svaba_germline_vcf, "SvABA", "dld1_multisample_germline"),
        (args.svaba_somatic_vcf, "SvABA", "dld1_multisample_somatic"),
    ]
    all_events, matrix_rows = [], []
    count_by_source = Counter()
    sample_sets = {}
    for path, caller, source in input_specs:
        groups, sample_map, raw_records = parse_events(path, caller, source)
        sample_sets[source] = sorted(set(sample_map.values()))
        count_by_source[f"{source}_input_vcf_records"] = raw_records
        for event_id, records in sorted(groups.items()):
            event, per_sample = summarize_event(event_id, records, caller, source, sample_map)
            all_events.append(event)
            count_by_source[f"{source}_normalized_events"] += 1
            for sid in DLD1_SAMPLES:
                values = per_sample.get(sid)
                if values is None:
                    values = {field: "" for field in NUMERIC_GRIDSS + NUMERIC_SVABA}
                    values.update({"GT": "", "caller_support": "", "caller_support_field": "",
                                   "call_status": "sample_not_in_this_caller", "detection_limit_3_support_fraction": ""})
                matrix_rows.append({**event, "sample_id": sid, "dld1_role": "case" if sid in DLD1_CASES else "control", **values})

    svaba_expected = set(DLD1_SAMPLES)
    for source in ("dld1_multisample_germline", "dld1_multisample_somatic"):
        if set(sample_sets[source]) != svaba_expected:
            raise SystemExit(f"{source} sample audit failed: {sample_sets[source]}")

    summary_rows = []
    grouped_matrix = defaultdict(list)
    for row in matrix_rows:
        grouped_matrix[(row["caller"], row["source"], row["caller_event_id"])].append(row)
    for _, rows in sorted(grouped_matrix.items()):
        case_support = sum(float_or_none(row["caller_support"]) not in (None, 0) for row in rows if row["dld1_role"] == "case")
        control_support = sum(float_or_none(row["caller_support"]) not in (None, 0) for row in rows if row["dld1_role"] == "control")
        pattern = (
            "case_and_control_observed" if case_support and control_support else
            "case_only_observed" if case_support else
            "control_only_observed" if control_support else
            "not_called_in_dld1_samples"
        )
        base = rows[0]
        summary_rows.append({
            "caller": base["caller"], "source": base["source"], "caller_event_id": base["caller_event_id"],
            "svtypes": base["svtypes"], "filters": base["filters"], "endpoints": base["endpoints"],
            "case_samples_with_caller_support": case_support,
            "control_samples_with_caller_support": control_support,
            "caller_support_pattern": pattern,
            "interpretation_limit": "caller pattern only; requires full-evidence, alternative-placement, and independent-template audit before clone-specific interpretation",
        })

    cluster_summary = read_cluster_summary(args.cluster_tsv)
    cluster_rows = [{"sample_id": sample, **cluster_summary[sample]} for sample in sorted(cluster_summary)]
    contig_summary = read_contig_summary(args.contig_summary)

    event_fields = ["caller", "source", "caller_event_id", "svtypes", "filters", "vcf_record_count", "endpoints"]
    matrix_fields = event_fields + ["sample_id", "dld1_role"] + list(NUMERIC_GRIDSS) + list(NUMERIC_SVABA) + [
        "GT", "caller_support", "caller_support_field", "call_status", "detection_limit_3_support_fraction"
    ]
    write_tsv(args.out_dir / "normalized_external_sv_events.tsv", event_fields, all_events)
    write_tsv(args.out_dir / "dld1_external_caller_matrix.tsv", matrix_fields, matrix_rows)
    write_tsv(args.out_dir / "dld1_external_caller_candidate_summary.tsv", [
        "caller", "source", "caller_event_id", "svtypes", "filters", "endpoints",
        "case_samples_with_caller_support", "control_samples_with_caller_support", "caller_support_pattern", "interpretation_limit"
    ], summary_rows)
    write_tsv(args.out_dir / "full_evidence_cluster_summary_by_sample.tsv", [
        "sample_id", "cluster_count", "clusters_ge3_fragments", "max_fragment_count"
    ], cluster_rows)
    audit = {
        "script": "7_normalize_external_sv_and_dld1.py",
        "research_question": "Caller-defined SV support across S01-S05; no caller result is treated as a biological conclusion.",
        "candidate_universe": "Target-related GRIDSS/SvABA BND/DEL/DUP/INV events only; raw evidence remains a separate non-closed universe.",
        "event_unit": "Linked VCF breakend records; mate records are collapsed per caller event to prevent representation double-counting.",
        "support_unit_limit": "GRIDSS VF and SvABA AD are caller fields, not independent-template counts.",
        "negative_limit": "not_called_by_this_caller is not event absence; SvABA DP and 3/DP are reported only as caller-depth sensitivity context.",
        "input_sample_sets": sample_sets,
        "counts": dict(sorted(count_by_source.items())),
        "contig_all_placement_summary": contig_summary,
    }
    (args.out_dir / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    (args.out_dir / ".complete").write_text("completed\n")


if __name__ == "__main__":
    main()
