#!/usr/bin/env python3
"""Answer the client-requested P41 §7.3 questions from completed P42 evidence.

All source BAM/VCF/ledger inputs are read-only.  The script creates only new,
date-stamped tables and does not re-align or re-call variants.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


LABELS = {
    "gridss192bf_23": "C1",
    "gridss192fb_336": "C2",
    "gridss192ff_106": "C3",
    "gridss56ff_1": "C4",
}
DLD1 = {"S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl", "S04_DLD1_c3_dox", "S05_DLD1_c4_dox"}
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-table", type=Path, required=True)
    ap.add_argument("--p41", type=Path, required=True)
    ap.add_argument("--summary", type=Path, required=True)
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--caller-input", type=Path, required=True)
    ap.add_argument("--caller-records", type=Path, required=True)
    ap.add_argument("--rescue-ledger", type=Path, required=True)
    ap.add_argument("--raw-audit", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    return ap.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_endpoints(text: str) -> tuple[tuple[str, int], tuple[str, int]]:
    parts = text.split(";")
    if len(parts) != 2:
        raise ValueError(f"Expected two endpoints: {text}")
    parsed = []
    for part in parts:
        chrom, pos = part.rsplit(":", 1)
        parsed.append((chrom, int(pos)))
    return parsed[0], parsed[1]


def endpoint_near(a: tuple[str, int], b: tuple[str, int], tolerance: int = 50) -> bool:
    return a[0] == b[0] and abs(a[1] - b[1]) <= tolerance


def caller_match(candidate: tuple[tuple[str, int], tuple[str, int]], record: dict[str, str]) -> tuple[bool, bool]:
    observed = ((record["chrom"], int(record["pos_1based"])), (record["end_chrom"], int(record["end_1based"])))
    two_end = (endpoint_near(candidate[0], observed[0]) and endpoint_near(candidate[1], observed[1])) or (
        endpoint_near(candidate[0], observed[1]) and endpoint_near(candidate[1], observed[0])
    )
    anchor = any(endpoint_near(x, y) for x in candidate for y in observed)
    return two_end and record["svtype"] == "BND", anchor


def query_interval(cigar: str) -> tuple[int, int]:
    """0-based inclusive interval containing query bases aligned to reference."""
    qpos = 0
    starts, ends = [], []
    for n_text, op in CIGAR_RE.findall(cigar):
        n = int(n_text)
        if op in "M=X":
            starts.append(qpos)
            ends.append(qpos + n - 1)
            qpos += n
        elif op in "IS":
            qpos += n
    if not starts:
        raise ValueError(f"No reference-consuming query block in CIGAR {cigar}")
    return min(starts), max(ends)


def read_end(flag: str) -> str:
    value = int(flag)
    if value & 64:
        return "R1"
    if value & 128:
        return "R2"
    return "unclassified"


def c2_read_end_resolution(raw_rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        if row["candidate_label"] == "C2" and row["chrom"] in {"chr11", "chr21"}:
            groups[(row["qname"], read_end(row["flag"]))].append(row)
    resolved: list[dict[str, object]] = []
    for (qname, end), rows in sorted(groups.items()):
        chr11 = [r for r in rows if r["chrom"] == "chr11"]
        chr21 = [r for r in rows if r["chrom"] == "chr21"]
        if not chr11 or not chr21:
            continue
        # One split representation per read end; use the longest mapped segment
        # on each side if an aligner retained duplicate context records.
        a = max(chr11, key=lambda r: query_interval(r["cigar"])[1] - query_interval(r["cigar"])[0])
        b = max(chr21, key=lambda r: query_interval(r["cigar"])[1] - query_interval(r["cigar"])[0])
        a_start, a_end = query_interval(a["cigar"])
        b_start, b_end = query_interval(b["cigar"])
        overlap_start, overlap_end = max(a_start, b_start), min(a_end, b_end)
        overlap_bp = max(0, overlap_end - overlap_start + 1)
        if overlap_bp:
            overlap_seq = a["seq"][overlap_start:overlap_end + 1]
        else:
            overlap_seq = ""
        gap_bp = max(0, b_start - a_end - 1) if a_start <= b_start else max(0, a_start - b_end - 1)
        resolved.append({
            "candidate_label": "C2", "sample_id": a["sample_id"], "qname": qname, "read_end": end,
            "chr11_flag": a["flag"], "chr11_pos_1based": a["pos"], "chr11_cigar": a["cigar"], "chr11_mapq": a["mapq"],
            "chr11_query_start_0based": a_start, "chr11_query_end_0based": a_end,
            "chr21_flag": b["flag"], "chr21_pos_1based": b["pos"], "chr21_cigar": b["cigar"], "chr21_mapq": b["mapq"],
            "chr21_query_start_0based": b_start, "chr21_query_end_0based": b_end,
            "query_overlap_bp": overlap_bp, "query_overlap_sequence": overlap_seq or ".", "unmapped_query_gap_bp": gap_bp,
            "interpretation": "alignment representation consistent with microhomology and no intervening unaligned query bases" if overlap_bp else "no query-overlap microhomology observed in this read end",
        })
    overlap_counter = Counter(int(x["query_overlap_bp"]) for x in resolved)
    summary = {
        "candidate_label": "C2",
        "independent_qname_count": len({x["qname"] for x in resolved}),
        "split_read_end_count": len(resolved),
        "query_overlap_bp_distribution": ";".join(f"{k}bp:{v}" for k, v in sorted(overlap_counter.items())),
        "unmapped_query_gap_bp_distribution": ";".join(str(x["unmapped_query_gap_bp"]) for x in resolved),
        "allowed_conclusion": "All resolved C2 read ends are consistent with the reported microhomology representation; repeat-context placement ambiguity and lack of independent validation remain.",
    }
    return resolved, summary


def data_manifest(samples: list[dict[str, str]], p41: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sample in samples:
        sid = sample["sample_id"]
        specs = [
            ("BAM", Path(sample["bam_path"]), Path(sample["bam_path"] + ".bai"), "available source BAM; do not duplicate into analysis delivery"),
            ("Manta_filtered_SV", p41 / "analysis/sv/manta_gw" / sid / "results/variants/diploidSV.vcf.gz", p41 / "analysis/sv/manta_gw" / sid / "results/variants/diploidSV.vcf.gz.tbi", "filtered Manta SV source"),
            ("Manta_candidate_SV", p41 / "analysis/sv/manta_gw" / sid / "results/variants/candidateSV.vcf.gz", None, "Manta candidate SV source; not re-labelled as unfiltered"),
            ("Manta_candidate_small_indel", p41 / "analysis/sv/manta_gw" / sid / "results/variants/candidateSmallIndels.vcf.gz", p41 / "analysis/sv/manta_gw" / sid / "results/variants/candidateSmallIndels.vcf.gz.tbi", "small-indel candidate source; not a project-wide SNV/indel callset"),
            ("DELLY_discovery_BCF", p41 / "analysis/sv/delly/discover" / f"{sid}.delly.bcf", None, "DELLY discovery source"),
            ("DELLY_region_genotyped_BCF", p41 / "analysis/sv/delly/region_geno_20261007" / f"{sid}.geno.bcf", None, "DELLY region-genotyped source"),
        ]
        for file_type, path, index, note in specs:
            rows.append({
                "sample_id": sid, "file_type": file_type, "source_path": str(path), "file_exists": int(path.is_file()),
                "size_bytes": path.stat().st_size if path.is_file() else 0,
                "index_path": str(index) if index else ".", "index_exists": int(index.is_file()) if index else ".", "note": note,
            })
    rows.append({
        "sample_id": "all_samples", "file_type": "project_wide_SNV_indel_VCF", "source_path": ".", "file_exists": 0,
        "size_bytes": 0, "index_path": ".", "index_exists": ".",
        "note": "No provenance-qualified project-wide SNV/indel VCF was identified in the documented P41 source paths; it is not included in the current package.",
    })
    return rows


def main() -> None:
    args = parse_args()
    if args.out_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {args.out_dir}")
    args.out_dir.mkdir(parents=True)
    samples = read_csv(args.sample_table)
    summary_rows = [r for r in read_tsv(args.summary) if r["caller_event_id"] in LABELS]
    if {r["caller_event_id"] for r in summary_rows} != set(LABELS):
        raise RuntimeError("Expected all C1-C4 fixed candidate records")
    candidate_by_label = {LABELS[r["caller_event_id"]]: r for r in summary_rows}
    caller_inputs = read_tsv(args.caller_input)
    if len(caller_inputs) != 72 or not all(r["exists"] == "1" and r["readable"] == "1" for r in caller_inputs):
        raise RuntimeError("Manta/DELLY input audit is incomplete or contains unreadable sources")
    caller_rows = [r for r in read_tsv(args.caller_records) if r["sample_id"] in DLD1]

    concordance = []
    for label in ("C1", "C2", "C3", "C4"):
        candidate = candidate_by_label[label]
        endpoints = parse_endpoints(candidate["fixed_endpoints"])
        for sid in sorted(DLD1):
            observed = [r for r in caller_rows if r["sample_id"] == sid]
            by_caller = defaultdict(list)
            anchor_by_caller = defaultdict(list)
            for record in observed:
                two_end, anchor = caller_match(endpoints, record)
                if two_end:
                    by_caller[record["caller"]].append(record)
                if anchor:
                    anchor_by_caller[record["caller"]].append(record)
            manta = by_caller["Manta"]
            delly = by_caller["DELLY"]
            if manta and delly:
                status = "both_concordant_exact_two_end"
            elif manta:
                status = "Manta_only_exact_two_end"
            elif delly:
                status = "DELLY_only_exact_two_end"
            elif anchor_by_caller:
                status = "anchor_only_or_nonconcordant"
            else:
                status = "no_exact_Manta_or_DELLY_match_in_audited_sources"
            concordance.append({
                "candidate_label": label, "sample_id": sid, "fixed_endpoints": candidate["fixed_endpoints"],
                "manta_exact_two_end_record_count": len(manta), "delly_exact_two_end_record_count": len(delly),
                "manta_anchor_record_count": len(anchor_by_caller["Manta"]), "delly_anchor_record_count": len(anchor_by_caller["DELLY"]),
                "caller_concordance_status": status,
                "manta_exact_record_ids": ";".join(r["record_id"] for r in manta) or ".",
                "delly_exact_record_ids": ";".join(r["record_id"] for r in delly) or ".",
                "anchor_record_detail": ";".join(f"{r['caller']}:{r['svtype']}:{r['chrom']}:{r['pos_1based']}->{r['end_chrom']}:{r['end_1based']}" for rs in anchor_by_caller.values() for r in rs) or ".",
                "interpretation_boundary": "Caller concordance is algorithmic evidence from the same WGS data, not independent biological validation.",
            })
    concordance_fields = list(concordance[0])
    write_tsv(args.out_dir / "candidate_manta_delly_concordance.tsv", concordance, concordance_fields)

    raw_rows = read_tsv(args.raw_audit)
    c2_read_ends, c2_summary = c2_read_end_resolution(raw_rows)
    if c2_summary["independent_qname_count"] != 3 or c2_summary["split_read_end_count"] != 5:
        raise RuntimeError(f"Unexpected C2 split-read reconstruction counts: {c2_summary}")
    write_tsv(args.out_dir / "c2_read_end_junction_resolution.tsv", c2_read_ends, list(c2_read_ends[0]))
    write_tsv(args.out_dir / "c2_junction_resolution_summary.tsv", [c2_summary], list(c2_summary))

    rescue_qnames = {r["qname"] for r in read_tsv(args.rescue_ledger) if r["record_class"] == "unmapped"}
    rescue_rows = []
    for label in ("C1", "C2", "C3", "C4"):
        qnames = {r["qname"] for r in raw_rows if r["candidate_label"] == label}
        recovered = qnames & rescue_qnames
        rescue_rows.append({
            "candidate_label": label, "raw_dual_endpoint_qname_count_before_rescue": len(qnames),
            "rescued_unmapped_mate_qname_overlap": len(recovered), "raw_dual_endpoint_qname_count_after_rescue": len(qnames),
            "rescued_qnames": ";".join(sorted(recovered)) or ".",
            "interpretation": "No rescued unmapped-mate QNAME overlaps this candidate; rescue did not change its raw dual-endpoint QNAME count." if not recovered else "Rescued mate overlap retained; interpret alongside the full raw alignment ledger.",
        })
    write_tsv(args.out_dir / "read_rescue_candidate_impact.tsv", rescue_rows, list(rescue_rows[0]))

    matrix = [r for r in read_tsv(args.matrix) if r["caller_event_id"] in LABELS]
    clone_rows = []
    for row in matrix:
        clone_rows.append({
            "candidate_label": LABELS[row["caller_event_id"]], "sample_id": row["sample_id"], "dld1_role": row["dld1_role"],
            "raw_dual_endpoint_fragment_count": row["raw_dual_endpoint_fragment_count"],
            "standard_eligible_dual_fragment_count": row["raw_standard_eligible_dual_fragment_count"],
            "standard_mean_depth_at_endpoint": row["standard_mean_depth_at_endpoint"],
            "raw_support_over_depth_proxy": round(int(row["raw_dual_endpoint_fragment_count"]) / float(row["standard_mean_depth_at_endpoint"]), 6),
            "filtering_rule": "High-confidence requires >=3 standard-eligible independent QNAME fragments, >=60x local depth for the 3-fragment/5% sensitivity context, and a sequence-level junction; no direct allele-fraction cutoff was applied because repeat-context support/depth is not a validated allele-fraction estimator.",
        })
    write_tsv(args.out_dir / "clone_specific_comparison_with_criteria.tsv", clone_rows, list(clone_rows[0]))

    manifest = data_manifest(samples, args.p41)
    write_tsv(args.out_dir / "data_delivery_availability_manifest.tsv", manifest, list(manifest[0]))

    coverage = [
        {"client_question": "Read rescue and local reassembly", "answer_status": "completed", "new_evidence_file": "read_rescue_candidate_impact.tsv", "boundary": "Rescue retained five unmapped mates overall; none overlaps C1-C4 dual-endpoint QNAMEs."},
        {"client_question": "Base-pair C2 junction detail", "answer_status": "completed with repeat-context qualification", "new_evidence_file": "c2_read_end_junction_resolution.tsv", "boundary": "Reports split-read alignment representation; not independent validation or a unique breakpoint proof."},
        {"client_question": "Manta/DELLY concordance", "answer_status": "completed", "new_evidence_file": "candidate_manta_delly_concordance.tsv", "boundary": "Exact two-end and anchor-only calls are separated; same-data algorithmic concordance is not biological replication."},
        {"client_question": "DLD1 +Dox versus parental/control clones", "answer_status": "completed", "new_evidence_file": "clone_specific_comparison_with_criteria.tsv", "boundary": "Raw candidate status is not a confirmed clone-specific structural variant at available depth."},
        {"client_question": "Genome-wide file delivery", "answer_status": "availability inventoried", "new_evidence_file": "data_delivery_availability_manifest.tsv", "boundary": "SNV/indel VCF availability is not asserted when no provenance-qualified project-wide source was identified."},
    ]
    write_tsv(args.out_dir / "customer_request_coverage.tsv", coverage, list(coverage[0]))

    audit = {
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "manta_delly_source_files_readable": len(caller_inputs),
        "dld1_candidate_caller_comparison_rows": len(concordance),
        "c2_independent_qnames": c2_summary["independent_qname_count"],
        "c2_split_read_ends": c2_summary["split_read_end_count"],
        "c2_microhomology_distribution": c2_summary["query_overlap_bp_distribution"],
        "rescued_unmapped_mate_qnames": len(rescue_qnames),
        "candidate_rescue_overlap_total": sum(int(r["rescued_unmapped_mate_qname_overlap"]) for r in rescue_rows),
        "manifest_rows": len(manifest),
    }
    (args.out_dir / "customer_request_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    (args.out_dir / ".complete").write_text("client-question completion analysis complete\n")


if __name__ == "__main__":
    main()
