"""Stage 13 — pre-delivery self-audit: every key number in the report is recomputed from the DELIVERED tables
(not from the analysis directories) and must appear in the report text; figures exist and were drawn with the
expected item counts. Exit 1 on any mismatch.
Usage: 13_self_audit.py <delivery_dir>
"""
import collections
import csv
import glob
import gzip
import os
import sys

D = sys.argv[1]
rep = open(glob.glob(f"{D}/*.md")[0]).read()
fails = []


def need(txt, why):
    if txt not in rep:
        fails.append(f"{why}: expected text not in report -> {txt!r}")


def tsv(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


T = tsv(f"{D}/tables/junction_table_target.tsv")
W = tsv(f"{D}/tables/junction_table_window_context.tsv")
need(f"{len(T)} distinct junctions have at least one end in the target region", "target total")
need(f"({len(W)} more lie only", "window context total")
cc = collections.Counter(r["classification"] for r in T)
need(f"Well supported {cc['Well_supported']}, Limited {cc['Limited']}, Ambiguous (repeat) {cc['Ambiguous_repeat']}, "
     f"Weak {cc['Weak']}, Not supported {cc['Not_supported']}", "class counts")
ns_before = sum(any(int(r[f"S{i:02d}_before"]) >= 1 for i in range(1, 11)) for r in T if r["classification"] == "Not_supported")
need(f"Of the {cc['Not_supported']} not supported, {ns_before} had at least one fragment", "not supported w/ before")
for i in range(1, 11):
    s = f"S{i:02d}"
    b = sum(int(r[f"{s}_before"]) >= 1 for r in T)
    a = sum(int(r[f"{s}_after"]) >= 1 for r in T)
    a2 = sum(int(r[f"{s}_after"]) >= 2 for r in T)
    need(f"| {s} | {b} | {a} | {a2} |", f"before/after row {s}")
cs = [r for r in T if r["DLD1_dox_comparison"].startswith("case_specific")]
need(f"DLD1: {len(cs)} target-region junctions are present in S04 or S05", "DLD1 case-specific")
gm = [r for r in T if r["GM12878_clone_specific"]]
need(f"GM12878: {len(gm)} junctions are present in one clone", "GM clone-specific")
assert max(int(r["max_after_frag_S01_S10"]) for r in gm) == 1 or fails.append("GM max fragments != 1")
bg = sum(int(r["background_S11_S18_samples_present"]) > 0 for r in T if int(r["max_after_frag_S01_S10"]) >= 1)
n1 = sum(1 for r in T if int(r["max_after_frag_S01_S10"]) >= 1)
need(f"{bg} of the {n1} events with at least one fragment", "background panel")
ge2 = [r for r in T if int(r["max_after_frag_S01_S10"]) >= 2]
need(f"After read rescue, {len(ge2)} have two or more independent fragments", "core conclusions >=2")
calls = collections.Counter(r["caller_support"] for r in T)
need(f"{calls.get('DELLY', 0)} called by DELLY only", "DELLY only")
need(f"{calls.get('Manta', 0)} by Manta only", "Manta only")
assert calls.get("Manta+DELLY", 0) == 0 or fails.append("report says 0 by both")
# key events
key = {r["junction_id"]: r for r in T}
for jid, s, bef, aft, af in (("P44J0332", "S04", 0, 4, "0.29"), ("P44J0115", "S05", 1, 2, "0.17"),
                             ("P44J0289", "S05", 3, 3, "0.27"), ("P44J0389", "S03", 1, 2, "0.09")):
    r = key[jid]
    if int(r[f"{s}_before"]) != bef or int(r[f"{s}_after"]) != aft:
        fails.append(f"{jid} {s} before/after {r[f'{s}_before']}/{r[f'{s}_after']} != {bef}/{aft}")
    if f"{float(r[f'{s}_AF']):.2f}" != af:
        fails.append(f"{jid} AF {r[f'{s}_AF']} != {af}")
    need(f"{bef} → {aft}", f"{jid} before->after")
    for b in ("S01", "S02", "S03"):
        if jid != "P44J0389" and r[f"{b}_state"] != "absent":
            fails.append(f"{jid} {b} state {r[f'{b}_state']} != absent")
# rescue table
for r in tsv(f"{D}/tables/rescue_source_summary.tsv"):
    s = r["sample"][:3]
    if int(s[1:]) <= 10:
        need(f"| {s} | {int(r['records_in_W']):,} | {int(r['qnames_total']):,} | {int(r['qnames_only_d']):,} |", f"rescue row {s}")
# depth boundary and BAF
for r in tsv(f"{D}/supporting_analyses/depth_boundary/boundary_estimates.tsv"):
    if r["mapq_min"] == "20" and r["sample"] in ("S04_DLD1_c3_dox", "S05_DLD1_c4_dox"):
        need(f"{int(r['boundary_bin_start']):,} | {int(r['ci95_lo']):,}–{int(r['ci95_hi']):,}", f"boundary {r['sample']}")
        need(f"{float(r['ratio_left']):.2f} / {float(r['ratio_right']):.2f}", f"boundary ratios {r['sample']}")
seg = {(r["sample"][:3], r["segment"]): r for r in tsv(f"{D}/supporting_analyses/chr11_BAF/baf_segments.tsv")}
need(f"{float(seg[('S04', 'distal_114.7-135.0Mb')]['frac_LOH_like']) * 100:.1f}% (S04)", "S04 LOH %")
need(f"{float(seg[('S05', 'distal_114.7-135.0Mb')]['frac_LOH_like']) * 100:.1f}% (S05)", "S05 LOH %")
rh = tsv(f"{D}/supporting_analyses/chr11_BAF/retained_homolog_S04_vs_S05.tsv")[0]
need(f"{int(rh['sites_LOH_in_both']):,} sites", "retained homolog n")
need(f"{float(rh['fraction_same']) * 100:.2f}%", "retained homolog fraction")
# figures
for f in ("P44_Fig1_locus_overview", "P44_Fig2_chr11_BAF", "P44_Fig3_rescue_sources", "P44_Fig4_before_after_and_classes"):
    for ext in ("png", "pdf"):
        if not os.path.exists(f"{D}/figures/{f}.{ext}"):
            fails.append(f"missing figure {f}.{ext}")
jv = sorted(glob.glob(f"{D}/figures/junction_views/*.png"))
ge2_all = [r for r in T + W if int(r["max_after_frag_S01_S10"]) >= 2]
if len(jv) != len(ge2_all):
    fails.append(f"junction views {len(jv)} != events with >=2 fragments {len(ge2_all)}")
need(f"({len(ge2_all)} junctions; PNG and PDF)", "junction view count in section 9")
# ---- round-2 additions (P1-P5)
import subprocess
ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"
paths = 0
for b in glob.glob(f"{D}/**/*.bam", recursive=True):
    paths += subprocess.run(f"{ENV}/samtools view --no-PG -H {b} | grep -c '/home/\\|/mnt/\\|/Work_bio'", shell=True,
                            capture_output=True, text=True).stdout.strip() not in ("0", "")
for v in glob.glob(f"{D}/**/*.vcf.gz", recursive=True):
    paths += subprocess.run(f"{ENV}/bcftools view --no-version -h {v} 2>/dev/null | grep -c '/home/\\|/mnt/\\|/Work_bio'",
                            shell=True, capture_output=True, text=True).stdout.strip() not in ("0", "")
if paths:
    fails.append(f"P1: {paths} BAM/VCF headers still contain local paths")
aud = tsv(f"{D}/methods/header_sanitization_audit.tsv")
if any(r["status"] != "ok" for r in aud):
    fails.append("P1: header_sanitization_audit has non-ok rows")
cov = tsv(f"{D}/tables/client_request_coverage.tsv")
for r in cov:
    for p in r["delivered_files"].split(";"):
        if p and not os.path.exists(f"{D}/{p}"):
            fails.append(f"P2: coverage names missing {p}")
led = tsv(f"{D}/tables/key_junction_read_evidence.tsv")
if {r["event_id"] for r in led} != {r["junction_id"] for r in T + W if int(r["max_after_frag_S01_S10"]) >= 2}:
    fails.append("P3: ledger events != events with >= 2 fragments")
pr = tsv(f"{D}/tables/primer_design.tsv")
if {r["event_id"] for r in pr} != {r["junction_id"] for r in T + W if int(r["max_after_frag_S01_S10"]) >= 2}:
    fails.append("P4: primer table events != events with >= 2 fragments")
inv = tsv(f"{D}/tables/data_inventory.tsv")
nfiles = sum(1 for f in glob.glob(f"{D}/**/*", recursive=True) if os.path.isfile(f) and not f.endswith("md5sums_delivery.txt"))
if len(inv) != nfiles:
    fails.append(f"P5: inventory rows {len(inv)} != files {nfiles}")

# banned recommendation (long-read platforms) and repudiation wording
for bad in ("long-read", "long read", "Nanopore", "PacBio", "incorrect", "mistake", "wrong", "erroneous"):
    if bad.lower() in rep.lower():
        fails.append(f"banned wording present: {bad}")
if fails:
    print("SELF-AUDIT FAILED")
    print("\n".join(fails))
    sys.exit(1)
print(f"SELF-AUDIT OK: {len(T)} target junctions, {len(jv)} junction views, all report numbers match delivered tables")
