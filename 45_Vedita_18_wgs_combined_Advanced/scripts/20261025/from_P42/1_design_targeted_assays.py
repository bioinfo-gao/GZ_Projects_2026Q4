#!/usr/bin/env python3
"""Build orientation-aware virtual BND templates and unchecked PCR-primer pairs.

This is an in-silico design step. It neither runs PCR nor converts WGS support
into independent biological confirmation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import primer3


LABELS = {
    "gridss192bf_23": "C1",
    "gridss192fb_336": "C2",
    "gridss192ff_106": "C3",
    "gridss56ff_1": "C4",
}
CASE_ONLY = "case_only_raw_candidate_requires_sensitivity_and_assembly"
BASES = set("ACGT")


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_tsv(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fetch_fasta(samtools: str, fasta: str, chrom: str, start: int, end: int) -> str:
    if start < 1:
        raise ValueError(f"Invalid reference interval {chrom}:{start}-{end}")
    out = subprocess.check_output([samtools, "faidx", fasta, f"{chrom}:{start}-{end}"], text=True)
    seq = "".join(line.strip() for line in out.splitlines() if not line.startswith(">"))
    expected = end - start + 1
    if len(seq) != expected:
        raise RuntimeError(f"Reference extraction length mismatch for {chrom}:{start}-{end}: {len(seq)} != {expected}")
    return seq.upper()


def parse_endpoint_record(text: str) -> dict:
    match = re.fullmatch(r"([^:;]+):(\d+):([ACGTN]+)>(.+)", text)
    if not match:
        raise ValueError(f"Cannot parse normalized BND endpoint: {text}")
    chrom, pos, ref, alt = match.groups()
    # Both BND placement forms are valid: local REF can precede the remote
    # location (``C]chr:pos]``) or follow it (``]chr:pos]C``).  Their
    # orientation is represented by the bracket character, not by assuming
    # REF is the first ALT character.
    if not (alt.startswith(ref) or alt.endswith(ref)):
        raise ValueError(f"Unsupported local-side BND representation: {text}")
    bracket = "[" if "[" in alt else "]" if "]" in alt else None
    if bracket is None:
        raise ValueError(f"No BND bracket in ALT: {text}")
    remote = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
    if remote is None:
        raise ValueError(f"Cannot parse remote endpoint from ALT: {text}")
    return {
        "local_chrom": chrom,
        "local_pos": int(pos),
        "ref": ref,
        "alt": alt,
        "remote_chrom": remote.group(1),
        "remote_pos": int(remote.group(2)),
        "bracket": bracket,
    }


def normalized_events(path: Path) -> dict[str, list[dict]]:
    result: dict[str, list[dict]] = {}
    for row in read_tsv(path):
        if row["caller_event_id"] not in LABELS:
            continue
        endpoints = [parse_endpoint_record(x) for x in row["endpoints"].split(";")]
        result[row["caller_event_id"]] = endpoints
    missing = set(LABELS).difference(result)
    if missing:
        raise RuntimeError(f"Missing normalized endpoint records: {sorted(missing)}")
    return result


def bam_paths(samples: Path) -> dict[str, str]:
    with samples.open() as handle:
        return {row["sample_id"]: row["bam_path"] for row in csv.DictReader(handle)}


def supporting_qnames(ledger: Path) -> dict[str, dict[str, set[str]]]:
    out: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in read_tsv(ledger):
        event_id = row["candidate_id"].rsplit("|", 1)[-1]
        if event_id in LABELS and int(row["raw_dual_endpoint_fragment_count"]) > 0:
            out[event_id][row["sample_id"]].add(row["qname"])
    return out


def supporting_sam_records(samtools: str, bam: str, endpoints: list[dict], qnames: set[str]) -> list[dict]:
    records: dict[tuple[str, int, str, str], dict] = {}
    for ep in endpoints:
        for chrom, pos in ((ep["local_chrom"], ep["local_pos"]), (ep["remote_chrom"], ep["remote_pos"])):
            region = f"{chrom}:{max(1, pos - 300)}-{pos + 300}"
            text = subprocess.check_output([samtools, "view", bam, region], text=True)
            for line in text.splitlines():
                fields = line.split("\t")
                if len(fields) < 11 or fields[0] not in qnames:
                    continue
                key = (fields[0], int(fields[1]), fields[2], fields[3])
                records[key] = {
                    "qname": fields[0], "flag": int(fields[1]), "chrom": fields[2],
                    "pos": int(fields[3]), "mapq": int(fields[4]), "cigar": fields[5],
                    "seq": fields[9].upper(),
                    "sa": next((x[5:] for x in fields[11:] if x.startswith("SA:Z:")), ""),
                }
    if not records:
        raise RuntimeError("No supporting SAM records recovered from the source BAM")
    return sorted(records.values(), key=lambda x: (x["qname"], x["chrom"], x["pos"], x["flag"]))


def crossing_support(template: str, junction: int, records: list[dict]) -> tuple[set[str], int]:
    """Count QNAMEs containing exact 24-mers spanning the proposed junction."""
    kmers = {template[start:start + 24] for start in range(junction - 16, junction - 7)}
    kmers = {k for k in kmers if len(k) == 24 and set(k).issubset(BASES)}
    qnames: set[str] = set()
    hit_count = 0
    for rec in records:
        orientations = (rec["seq"], revcomp(rec["seq"]))
        if any(kmer in seq for kmer in kmers for seq in orientations):
            qnames.add(rec["qname"])
            hit_count += 1
    return qnames, hit_count


def choose_template(samtools: str, fasta: str, endpoint: dict, records: list[dict]) -> dict:
    best = None
    for local_shift, remote_shift in itertools.product(range(-4, 5), repeat=2):
        local_pos = endpoint["local_pos"] + local_shift
        remote_pos = endpoint["remote_pos"] + remote_shift
        left = fetch_fasta(samtools, fasta, endpoint["local_chrom"], local_pos - 399, local_pos)
        # A BND ALT describes breakend geometry, but raw split-read sequence is
        # the decisive evidence for which reference-strand sequence crosses the
        # observed molecule.  Evaluate both remote orientations instead of
        # enforcing an ALT-string assumption.  This prevents a representation
        # convention from silently turning observed split-read support into a
        # false zero-concordance result.
        orientations = [
            ("forward_from_remote_position", fetch_fasta(samtools, fasta, endpoint["remote_chrom"], remote_pos, remote_pos + 399)),
            ("reverse_complement_to_remote_position", revcomp(fetch_fasta(samtools, fasta, endpoint["remote_chrom"], remote_pos - 399, remote_pos))),
        ]
        for remote_orientation, right in orientations:
            template = left + right
            qnames, hits = crossing_support(template, len(left), records)
            score = (
                len(qnames), hits, -abs(local_shift) - abs(remote_shift),
                -abs(local_shift), -abs(remote_shift),
                int(remote_orientation == "forward_from_remote_position"),
            )
            candidate = {
                "template": template, "junction_index": len(left), "local_shift": local_shift,
                "remote_shift": remote_shift, "raw_concordant_qnames": qnames,
                "raw_concordant_record_count": hits, "remote_orientation": remote_orientation,
                "orientation_selection_basis": "maximum exact 24-mer support in recovered raw split-read sequences",
                "score": score,
            }
            if best is None or candidate["score"] > best["score"]:
                best = candidate
    assert best is not None
    return best


def primer_metrics(seq: str) -> dict | None:
    if len(seq) < 19 or len(seq) > 25 or not set(seq).issubset(BASES):
        return None
    gc = 100.0 * (seq.count("G") + seq.count("C")) / len(seq)
    tm = primer3.calc_tm(seq)
    hairpin = primer3.calc_hairpin(seq).tm
    homodimer = primer3.calc_homodimer(seq).tm
    if not (58.0 <= tm <= 62.5 and 40.0 <= gc <= 60.0):
        return None
    if re.search(r"([ACGT])\1{4,}", seq) or seq.endswith("GGGG") or seq.endswith("CCCC"):
        return None
    if hairpin > 45.0 or homodimer > 45.0:
        return None
    return {"seq": seq, "tm": round(tm, 3), "gc_percent": round(gc, 3),
            "hairpin_tm": round(hairpin, 3), "homodimer_tm": round(homodimer, 3)}


def design_pairs(template: str, junction: int, raw_records: list[dict]) -> list[dict]:
    left: dict[str, dict] = {}
    right: dict[str, dict] = {}
    for length in range(20, 25):
        for end in range(junction - 25, junction - 181, -1):
            start = end - length
            if start < 0:
                continue
            metric = primer_metrics(template[start:end])
            if metric:
                metric.update({"start": start, "end": end, "side": "left"})
                left.setdefault(metric["seq"], metric)
        for start in range(junction + 25, min(len(template) - length, junction + 181)):
            end = start + length
            metric = primer_metrics(revcomp(template[start:end]))
            if metric:
                metric.update({"start": start, "end": end, "side": "right"})
                right.setdefault(metric["seq"], metric)
    def single_rank(m):
        return (abs(m["tm"] - 60.0) + abs(m["gc_percent"] - 50.0) / 10.0 + max(m["hairpin_tm"], m["homodimer_tm"]) / 100.0, m["start"])
    lefts = sorted(left.values(), key=single_rank)[:8]
    rights = sorted(right.values(), key=single_rank)[:8]
    raw_sequences = [(r["qname"], s) for r in raw_records for s in (r["seq"], revcomp(r["seq"]))]
    pairs = []
    for l, r in itertools.product(lefts, rights):
        if l["seq"] == r["seq"]:
            continue
        product_bp = r["end"] - l["start"]
        if not 100 <= product_bp <= 360:
            continue
        heterodimer_tm = primer3.calc_heterodimer(l["seq"], r["seq"]).tm
        if heterodimer_tm > 45.0:
            continue
        lraw = sorted({q for q, s in raw_sequences if l["seq"] in s})
        rraw = sorted({q for q, s in raw_sequences if r["seq"] in s})
        penalty = (abs(l["tm"] - 60) + abs(r["tm"] - 60) + abs(product_bp - 220) / 100.0)
        pairs.append({
            "left": l, "right": r, "product_bp": product_bp, "pair_penalty": round(penalty, 4),
            "heterodimer_tm": round(heterodimer_tm, 3),
            "left_raw_exact_qnames": ";".join(lraw), "right_raw_exact_qnames": ";".join(rraw),
        })
    return sorted(pairs, key=lambda x: (x["pair_penalty"], -len(x["left_raw_exact_qnames"]), -len(x["right_raw_exact_qnames"])))[:25]


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, type=Path)
    ap.add_argument("--ledger", required=True, type=Path)
    ap.add_argument("--normalized", required=True, type=Path)
    ap.add_argument("--samples", required=True, type=Path)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--samtools", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=False)

    case_rows = [r for r in read_tsv(args.summary) if r["classification"] == CASE_ONLY]
    if {r["caller_event_id"] for r in case_rows} != set(LABELS):
        raise RuntimeError("The C1-C4 fixed-candidate set differs from the recorded four case-only candidates")
    event_map = normalized_events(args.normalized)
    bam_map = bam_paths(args.samples)
    qnames = supporting_qnames(args.ledger)

    contexts, raw_audit, pair_rows, fasta_records = [], [], [], []
    for row in sorted(case_rows, key=lambda x: LABELS[x["caller_event_id"]]):
        event_id = row["caller_event_id"]
        label = LABELS[event_id]
        endpoints = event_map[event_id]
        endpoint = next((x for x in endpoints if x["local_chrom"] == "chr11"), endpoints[0])
        event_qnames = qnames[event_id]
        records = []
        for sample, names in event_qnames.items():
            recs = supporting_sam_records(args.samtools, bam_map[sample], endpoints, names)
            for rec in recs:
                rec["candidate_label"] = label
                rec["sample_id"] = sample
            records.extend(recs)
        template_info = choose_template(args.samtools, args.reference, endpoint, records)
        template = template_info["template"]
        concordant = template_info["raw_concordant_qnames"]
        fasta_records.append((label, template, template_info))
        raw_support_qnames = sorted({r["qname"] for r in records})
        contexts.append({
            "candidate_label": label, "caller_event_id": event_id,
            "fixed_endpoints": row["fixed_endpoints"], "normalized_bnd_records": ";".join(
                f"{x['local_chrom']}:{x['local_pos']}:{x['ref']}>{x['alt']}" for x in endpoints),
            "template_local_chrom": endpoint["local_chrom"], "template_local_position": endpoint["local_pos"] + template_info["local_shift"],
            "template_remote_chrom": endpoint["remote_chrom"], "template_remote_position": endpoint["remote_pos"] + template_info["remote_shift"],
            "vcf_remote_bracket": endpoint["bracket"], "remote_orientation": template_info["remote_orientation"],
            "orientation_selection_basis": template_info["orientation_selection_basis"], "template_length_bp": len(template),
            "junction_index_0based": template_info["junction_index"], "local_coordinate_shift_bp": template_info["local_shift"],
            "remote_coordinate_shift_bp": template_info["remote_shift"],
            "raw_dual_endpoint_qnames": ";".join(raw_support_qnames),
            "raw_dual_endpoint_qname_count": len(raw_support_qnames),
            "raw_junction_concordant_qnames": ";".join(sorted(concordant)),
            "raw_junction_concordant_qname_count": len(concordant),
            "raw_junction_concordant_record_count": template_info["raw_concordant_record_count"],
            "template_sha256": hashlib.sha256(template.encode()).hexdigest(),
            "classification": row["classification"],
            "design_boundary": "in-silico template only; not PCR evidence or biological confirmation",
        })
        for rec in records:
            raw_audit.append({"candidate_label": label, "caller_event_id": event_id, **rec})
        if concordant:
            for i, pair in enumerate(design_pairs(template, template_info["junction_index"], records), start=1):
                pair_rows.append({
                    "candidate_label": label, "caller_event_id": event_id, "unchecked_pair_rank": i,
                    "left_primer": pair["left"]["seq"], "left_tm": pair["left"]["tm"], "left_gc_percent": pair["left"]["gc_percent"],
                    "left_hairpin_tm": pair["left"]["hairpin_tm"], "left_homodimer_tm": pair["left"]["homodimer_tm"],
                    "left_template_start_0based": pair["left"]["start"], "left_template_end_0based": pair["left"]["end"],
                    "right_primer": pair["right"]["seq"], "right_tm": pair["right"]["tm"], "right_gc_percent": pair["right"]["gc_percent"],
                    "right_hairpin_tm": pair["right"]["hairpin_tm"], "right_homodimer_tm": pair["right"]["homodimer_tm"],
                    "right_template_start_0based": pair["right"]["start"], "right_template_end_0based": pair["right"]["end"],
                    "expected_junction_product_bp": pair["product_bp"], "pair_penalty": pair["pair_penalty"], "pair_heterodimer_tm": pair["heterodimer_tm"],
                    "left_raw_exact_qnames": pair["left_raw_exact_qnames"], "right_raw_exact_qnames": pair["right_raw_exact_qnames"],
                    "unchecked_status": "requires whole-reference primer specificity audit",
                })

    with (args.out_dir / "virtual_junction_templates.fasta").open("w") as handle:
        for label, seq, info in fasta_records:
            handle.write(f">{label}|junction_index={info['junction_index']}|raw_concordant_qnames={len(info['raw_concordant_qnames'])}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")
    write_tsv(args.out_dir / "candidate_template_audit.tsv", contexts, list(contexts[0]))
    write_tsv(args.out_dir / "raw_support_alignment_audit.tsv", raw_audit, list(raw_audit[0]))
    fields = [
        "candidate_label", "caller_event_id", "unchecked_pair_rank", "left_primer", "left_tm", "left_gc_percent", "left_hairpin_tm", "left_homodimer_tm", "left_template_start_0based", "left_template_end_0based",
        "right_primer", "right_tm", "right_gc_percent", "right_hairpin_tm", "right_homodimer_tm", "right_template_start_0based", "right_template_end_0based", "expected_junction_product_bp", "pair_penalty", "pair_heterodimer_tm",
        "left_raw_exact_qnames", "right_raw_exact_qnames", "unchecked_status",
    ]
    write_tsv(args.out_dir / "primer_pairs_unchecked.tsv", pair_rows, fields)
    (args.out_dir / ".design_complete").write_text("orientation-aware template and unchecked primer design complete\n")


if __name__ == "__main__":
    main()
