#!/usr/bin/env python3
"""Run a whole-reference exact-primer specificity audit and select in-silico pairs."""
from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(x: str) -> float:
    return float(x)


def expected_mapping(pair: dict, context: dict, side: str) -> tuple[str, int, int, str]:
    """Return the exact reference locus/strand that the primer must hit."""
    if side == "left":
        start = int(pair["left_template_start_0based"])
        end = int(pair["left_template_end_0based"])
        local_pos = int(context["template_local_position"])
        return (
            context["template_local_chrom"],
            local_pos - 399 + start,
            local_pos - 399 + end - 1,
            "plus",
        )
    start = int(pair["right_template_start_0based"]) - int(context["junction_index_0based"])
    end = int(pair["right_template_end_0based"]) - int(context["junction_index_0based"])
    remote_pos = int(context["template_remote_position"])
    if context["remote_orientation"] == "forward_from_remote_position":
        return context["template_remote_chrom"], remote_pos + start, remote_pos + end - 1, "minus"
    return context["template_remote_chrom"], remote_pos - end + 1, remote_pos - start, "plus"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--unchecked-pairs", type=Path, required=True)
    ap.add_argument("--template-audit", type=Path, required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--makeblastdb", required=True)
    ap.add_argument("--blastn", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    pairs = read_tsv(args.unchecked_pairs)
    contexts = {r["candidate_label"]: r for r in read_tsv(args.template_audit)}
    if set(contexts) != {"C1", "C2", "C3", "C4"}:
        raise RuntimeError("Expected exactly C1-C4 template contexts")
    primer_queries = []
    seen = set()
    for pair in pairs:
        context = contexts[pair["candidate_label"]]
        for side in ("left", "right"):
            seq = pair[f"{side}_primer"]
            key = (pair["candidate_label"], side, seq)
            if key not in seen:
                seen.add(key)
                chrom, start, end, strand = expected_mapping(pair, context, side)
                primer_queries.append({"query_id": "|".join(key), "candidate_label": pair["candidate_label"], "side": side,
                                       "primer_sequence": seq, "expected_chrom": chrom, "expected_start": start,
                                       "expected_end": end, "expected_strand": strand})
    blast_out = args.out_dir / "primer_blastn_short.tsv"
    if primer_queries:
        db_dir = args.out_dir / "grch38_blast_db"
        db_dir.mkdir(exist_ok=False)
        db_prefix = db_dir / "Homo_sapiens_assembly38"
        subprocess.run([
            args.makeblastdb, "-in", args.reference, "-dbtype", "nucl", "-parse_seqids",
            "-out", str(db_prefix), "-title", "P42_GRCh38_assay_specificity_audit",
        ], check=True)
        fasta = args.out_dir / "primer_queries.fasta"
        with fasta.open("w") as handle:
            for r in primer_queries:
                handle.write(f">{r['query_id']}\n{r['primer_sequence']}\n")
        subprocess.run([
            args.blastn, "-task", "blastn-short", "-query", str(fasta), "-db", str(db_prefix),
            "-word_size", "7", "-dust", "no", "-evalue", "1000", "-max_target_seqs", "100000", "-max_hsps", "100000",
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand", "-out", str(blast_out),
        ], check=True)
    else:
        blast_out.write_text("")

    hits = defaultdict(list)
    if blast_out.stat().st_size:
        with blast_out.open() as handle:
            for line in handle:
                qid, sid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, strand = line.rstrip("\n").split("\t")
                hits[qid].append({"subject_chrom": sid, "pident": float(pident), "alignment_length": int(length),
                                  "mismatch": int(mismatch), "gapopen": int(gapopen), "subject_start": int(sstart),
                                  "subject_end": int(send), "subject_strand": strand})

    specificity = []
    primer_status = {}
    for query in primer_queries:
        exact = [h for h in hits[query["query_id"]] if h["pident"] == 100.0 and h["alignment_length"] == len(query["primer_sequence"])
                 and h["mismatch"] == 0 and h["gapopen"] == 0]
        expected = [h for h in exact if h["subject_chrom"] == query["expected_chrom"]
                    and min(h["subject_start"], h["subject_end"]) == int(query["expected_start"])
                    and max(h["subject_start"], h["subject_end"]) == int(query["expected_end"])
                    and h["subject_strand"] == query["expected_strand"]]
        status = "unique_expected_locus" if len(exact) == 1 and len(expected) == 1 else "nonunique_or_expected_locus_not_unique"
        primer_status[query["query_id"]] = status
        specificity.append({**query, "exact_full_length_hit_count": len(exact), "expected_chrom_exact_hit_count": len(expected),
                            "specificity_status": status,
                            "exact_hit_loci": ";".join(f"{h['subject_chrom']}:{h['subject_start']}-{h['subject_end']}({h['subject_strand']})" for h in exact)})

    selections = []
    for label in ("C1", "C2", "C3", "C4"):
        context = contexts[label]
        cands = [p for p in pairs if p["candidate_label"] == label]
        chosen = None
        for p in cands:
            lk = f"{label}|left|{p['left_primer']}"
            rk = f"{label}|right|{p['right_primer']}"
            if primer_status.get(lk) == "unique_expected_locus" and primer_status.get(rk) == "unique_expected_locus":
                chosen = p
                break
        concordant_n = int(context["raw_junction_concordant_qname_count"])
        if concordant_n == 0:
            status = "no_in_silico_template_supported_by_raw_junction_sequence"
            reason = "No reference-derived template offset had an exact 24-mer spanning match in the recovered supporting read sequences."
        elif chosen is None:
            status = "no_pair_passed_whole_reference_exact_specificity"
            reason = "Unchecked candidates exist, but no pair had one full-length exact hit for each primer at the expected GRCh38 loci."
        elif label == "C2":
            status = "primary_in_silico_assay_ready_not_biologically_validated"
            reason = "C2 is the P41 EV0623 carry-forward; the selected pair is raw-template concordant and passes the stated in-silico specificity criterion."
        else:
            status = "second_tier_in_silico_assay_ready_not_biologically_validated"
            reason = "A second-tier candidate has a raw-template-concordant pair passing the stated in-silico specificity criterion."
        output = {
            "candidate_label": label, "caller_event_id": context["caller_event_id"], "fixed_endpoints": context["fixed_endpoints"],
            "raw_dual_endpoint_qname_count": context["raw_dual_endpoint_qname_count"],
            "raw_junction_concordant_qname_count": concordant_n, "in_silico_design_status": status,
            "interpretation": reason,
            "external_requirement": "PCR requires independently prepared DNA and empirical sensitivity; WGS comparison requires new measured local depth if a negative control conclusion is sought.",
        }
        if chosen:
            output.update({
                "selected_unchecked_pair_rank": chosen["unchecked_pair_rank"], "left_primer": chosen["left_primer"], "right_primer": chosen["right_primer"],
                "expected_junction_product_bp": chosen["expected_junction_product_bp"], "pair_penalty": chosen["pair_penalty"],
                "left_specificity_status": primer_status[f"{label}|left|{chosen['left_primer']}"],
                "right_specificity_status": primer_status[f"{label}|right|{chosen['right_primer']}"],
            })
        selections.append(output)

    specificity_fields = [
        "query_id", "candidate_label", "side", "primer_sequence", "expected_chrom", "expected_start", "expected_end", "expected_strand",
        "exact_full_length_hit_count", "expected_chrom_exact_hit_count", "specificity_status", "exact_hit_loci",
    ]
    write_tsv(args.out_dir / "primer_specificity.tsv", specificity, specificity_fields)
    fields = [
        "candidate_label", "caller_event_id", "fixed_endpoints", "raw_dual_endpoint_qname_count", "raw_junction_concordant_qname_count",
        "in_silico_design_status", "selected_unchecked_pair_rank", "left_primer", "right_primer", "expected_junction_product_bp", "pair_penalty",
        "left_specificity_status", "right_specificity_status", "interpretation", "external_requirement",
    ]
    write_tsv(args.out_dir / "assay_design_selection.tsv", selections, fields)
    (args.out_dir / ".specificity_complete").write_text("whole-reference primer specificity audit complete\n")


if __name__ == "__main__":
    main()
