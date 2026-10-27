"""Stage 13 (P45; extends P44 13_self_audit.py) — pre-delivery self-audit. Every key number in the report is recomputed
from the DELIVERED files (never from analysis directories) and must appear in the report text; exit 1 on any mismatch.
P45 additions: key events located by coordinates (IDs are not hard-coded), claim levels, four-caller counts, flank
events, G1 genome-wide SNV/indel counts, G2 genome-wide DLD1 counts, primer loci, DELLY conversion audit, no internal
project numbers, figure counts.
Usage: 13_self_audit.py <delivery_dir>
"""
import collections
import csv
import glob
import gzip
import os
import re
import subprocess
import sys

D = sys.argv[1]
rep = open(glob.glob(f"{D}/*.md")[0]).read()
fails = []
ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"


def need(txt, why):
    if txt not in rep:
        fails.append(f"{why}: expected text not in report -> {txt!r}")


def tsv(p):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"))


def c(n):
    return f"{n:,}"


T = tsv(f"{D}/tables/junction_table_target.tsv")
W = tsv(f"{D}/tables/junction_table_window_context.tsv")
need(f"{len(T)} distinct junctions have at least one end in the target region", "target total")
need(f"({c(len(W))} more lie only in the ± 50 kb flanks", "window context total")
cc = collections.Counter(r["classification"] for r in T)
need(f"Well supported {cc['Well_supported']}, Limited {cc['Limited']}, Ambiguous (repeat) {cc['Ambiguous_repeat']}, "
     f"Weak {cc['Weak']}, Not supported {cc['Not_supported']}", "class counts")
ns_before = sum(any(int(r[f"S{i:02d}_before"]) >= 1 for i in range(1, 11)) for r in T if r["classification"] == "Not_supported")
need(f"Of the {cc['Not_supported']} not supported, {ns_before} had at least one fragment", "not supported w/ before")
need(f"{cc['Not_supported']} of {len(T)} target-region candidates are better or equally explained", "key finding 2c")
need(f"{cc['Not_supported']} are explained at least as well by an unrearranged copy", "core conclusions NS")
for i in range(1, 11):
    s = f"S{i:02d}"
    b = sum(int(r[f"{s}_before"]) >= 1 for r in T)
    a = sum(int(r[f"{s}_after"]) >= 1 for r in T)
    a2 = sum(int(r[f"{s}_after"]) >= 2 for r in T)
    need(f"| {s} | {b} | {a} | {a2} |", f"before/after row {s}")
ge2_t = [r for r in T if int(r["max_after_frag_S01_S10"]) >= 2]
ge2_all = [r for r in T + W if int(r["max_after_frag_S01_S10"]) >= 2]
need(f"After read rescue, {len(ge2_t)} have two or more independent fragments", "core >=2")
need(f"In the ± 50 kb flanks, {len(ge2_all) - len(ge2_t)} further junctions have two or more fragments", "flank >=2")
need(f"All {len(ge2_all)} junctions in the window (target region ± 50 kb) with two or more fragments", "Figure 3 caption (summary matrix; 1027 wording)")
need(f"(one of {len(ge2_all)}; all in", "Figure 4 caption (junction views)")
need(f"({len(ge2_all)} junctions; PNG and PDF)", "section 9 junction views")
cs = [r for r in T if r["DLD1_dox_comparison"].startswith("case_specific")]
cnt = collections.Counter(r["DLD1_dox_case_samples"] for r in cs)
need(f"{len(cs)} junctions are present in S04 or S05 and in none of S01–S03 ({cnt['S04']} in S04, {cnt['S05']} in S05)", "DLD1 case-specific")
need(f"{sum(r['DLD1_dox_comparison'] == 'case_specific_underpowered' for r in cs)} of the {len(cs)} have at least one underpowered", "underpowered")
gm = [r for r in T if r["GM12878_clone_specific"]]
need(f"{len(gm)} junctions are present in one clone and absent from S06", "GM clone-specific")
if max(int(r["max_after_frag_S01_S10"]) for r in gm) != 1:
    fails.append("GM clone-specific max fragments != 1")
bg = sum(int(r["background_S11_S18_samples_present"]) > 0 for r in T if int(r["max_after_frag_S01_S10"]) >= 1)
n1 = sum(1 for r in T if int(r["max_after_frag_S01_S10"]) >= 1)
need(f"{bg} of the {n1} events with at least one fragment", "background panel")
if any(int(r["background_S11_S18_samples_present"]) > 0 for r in ge2_t):
    fails.append("a Section 7.4 event has background-panel support but the report says none does")
# four-caller counts
calls = collections.Counter(r["caller_support"] for r in T)
for k, txt in (("DELLY+GRIDSS+SvABA", "| DELLY + GRIDSS + SvABA | {} (J0786, S05) |"), ("Manta+GRIDSS", "| Manta + GRIDSS | {} |"),
               ("GRIDSS", "| GRIDSS only | {} |"), ("SvABA", "| SvABA only | {} |"), ("none", "| none | {} |")):
    need(txt.format(calls.get(k, 0)), f"caller row {k}")
both_md = sum(1 for r in T if "Manta" in r["caller_support"] and "DELLY" in r["caller_support"])
need(f"| Manta and DELLY together | {both_md} |", "Manta+DELLY")
need(f"{calls.get('GRIDSS', 0) + calls.get('Manta+GRIDSS', 0)} are called by GRIDSS alone or together with Manta", "GRIDSS count text")
# key events by coordinate
def find(rows, a_start):
    """the junction with >= 2 fragments whose breakpoint interval A starts at a_start (0-fragment junctions can share it)"""
    hit = [r for r in rows if r["breakpoint_interval_A"].startswith(str(a_start)) and int(r["max_after_frag_S01_S10"]) >= 2]
    if len(hit) != 1:
        fails.append(f"key event at {a_start}: {len(hit)} rows")
        return None
    return hit[0]


KEY = ((114605216, "S04", 0, 4, "0.29", "Ambiguous_repeat", "structurally_coherent_limited_candidate"),
       (114605971, "S05", 1, 2, "0.17", "Limited", "structurally_coherent_limited_candidate"),
       (114607190, "S05", 3, 3, "0.27", "Limited", "limited_junction_candidate"),
       (114609273, "S03", 1, 2, "0.09", "Limited", "limited_junction_candidate"))
for a0, s, bef, aft, af, cls, claim in KEY:
    r = find(T, a0)
    if not r:
        continue
    if (int(r[f"{s}_before"]), int(r[f"{s}_after"])) != (bef, aft):
        fails.append(f"{r['junction_id']} {s} before/after {r[f'{s}_before']}/{r[f'{s}_after']} != {bef}/{aft}")
    if f"{float(r[f'{s}_AF']):.2f}" != af:
        fails.append(f"{r['junction_id']} AF {r[f'{s}_AF']} != {af}")
    if r["classification"] != cls or r["claim_level"] != claim:
        fails.append(f"{r['junction_id']} class/claim {r['classification']}/{r['claim_level']} != {cls}/{claim}")
    need(f"| {r['junction_id']} |", f"{r['junction_id']} row in 7.4")
    need(f"{bef} → {aft}", f"{r['junction_id']} before->after")
    if s != "S03":
        refs = ", ".join(r[f"{b}_ref"] for b in ("S01", "S02", "S03"))
        if any(r[f"{b}_state"] != "absent" for b in ("S01", "S02", "S03")):
            fails.append(f"{r['junction_id']} baseline not all absent")
        need(f"({refs}", f"{r['junction_id']} baseline reference fragments")
need("absent (13, 12, 24 reference fragments)", "S04 baseline refs wording")
# flank events named in 7.8
for a0, samples in ((114590721, "S06,S07,S08,S09,S10"), (114590798, "S06,S07,S08,S09,S10")):
    r = find(W, a0)
    if r and (r["samples_present"] != samples or r["classification"] != "Well_supported" or r["claim_level"] != "well_supported_junction_candidate"):
        fails.append(f"flank event {a0}: {r['samples_present']} {r['classification']} {r['claim_level']}")
    if r:
        lo, hi = min(int(r[f"S{i:02d}_after"]) for i in range(6, 11)), max(int(r[f"S{i:02d}_after"]) for i in range(6, 11))
        if any(int(r[f"S{i:02d}_after"]) for i in range(1, 6)):
            fails.append(f"flank event {a0} present in DLD1")
        if (lo, hi) != (3, 10) and a0 == 114590721:
            fails.append(f"J1365 fragment range {lo}-{hi}")
for a0, s in ((114629853, "S05"), (114602525, "S04")):
    r = find(W, a0)
    if r and r["samples_present"] != s:
        fails.append(f"flank event {a0} samples {r['samples_present']} != {s}")
# coherence table agrees with claim levels
coh = tsv(f"{D}/tables/structural_coherence_tests.tsv")
coh_yes = {(r["junction_id"]) for r in coh if r["coherent"] == "1"}
claim_yes = {r["junction_id"] for r in T + W if r["claim_level"] == "structurally_coherent_limited_candidate"}
if coh_yes != claim_yes:
    fails.append(f"coherence tests {sorted(coh_yes)} != claim level {sorted(claim_yes)}")
if any(r["claim_level"] == "confirmed" for r in T + W):
    fails.append("claim level 'confirmed' present without an independent assay")
# rescue table
for r in tsv(f"{D}/tables/rescue_source_summary.tsv"):
    s = r["sample"][:3]
    if int(s[1:]) <= 10:
        need(f"| {s} | {c(int(r['records_in_W']))} | {c(int(r['qnames_total']))} | {c(int(r['qnames_only_d']))} |", f"rescue row {s}")
rs = [r for r in tsv(f"{D}/tables/rescue_source_summary.tsv") if int(r["sample"][1:3]) <= 10]
need(f"{c(min(int(r['qnames_total']) for r in rs))}–{c(max(int(r['qnames_total']) for r in rs))} read pairs were rescued", "rescue range")
need(f"{c(min(int(r['qnames_only_d']) for r in rs))}–{c(max(int(r['qnames_only_d']) for r in rs))} of them were placed elsewhere", "XA range")
# filter impact
fi = {r["item"]: int(r["count"]) for r in tsv(f"{D}/tables/evidence_filter_impact.tsv")}
need(f"{c(fi['orig_split_adapter_readthrough'])} such split alignments, {c(fi['orig_clip_adapter'])} adapter soft clips and "
     f"{c(fi['disc_low_complexity_mate'])} discordant pairs", "artefact counts")
# depth boundary and BAF
for r in tsv(f"{D}/supporting_analyses/depth_boundary/boundary_estimates.tsv"):
    if r["mapq_min"] == "20" and r["sample"] in ("S04_DLD1_c3_dox", "S05_DLD1_c4_dox"):
        need(f"{c(int(r['boundary_bin_start']))} | {c(int(r['ci95_lo']))}–{c(int(r['ci95_hi']))}", f"boundary {r['sample']}")
        need(f"{float(r['ratio_left']):.2f} / {float(r['ratio_right']):.2f}", f"boundary ratios {r['sample']}")
seg = {(r["sample"][:3], r["segment"]): r for r in tsv(f"{D}/supporting_analyses/chr11_BAF/baf_segments.tsv")}
need(f"{float(seg[('S04', 'distal_114.7-135.0Mb')]['frac_LOH_like']) * 100:.1f}% (S04)", "S04 LOH %")
need(f"{float(seg[('S05', 'distal_114.7-135.0Mb')]['frac_LOH_like']) * 100:.1f}% (S05)", "S05 LOH %")
rh = tsv(f"{D}/supporting_analyses/chr11_BAF/retained_homolog_S04_vs_S05.tsv")[0]
need(f"{c(int(rh['sites_LOH_in_both']))} sites", "retained homolog n")
need(f"{float(rh['fraction_same']) * 100:.2f}%", "retained homolog fraction")
# G1 genome-wide SNV/indel
G1 = f"{D}/genome_wide_SNV_indel_VCF"
n_raw = int(subprocess.run(f"{ENV}/bcftools index -n {G1}/S01_S10.genome.raw.vcf.gz", shell=True, capture_output=True, text=True).stdout)
n_q = int(subprocess.run(f"{ENV}/bcftools index -n {G1}/S01_S10.genome.QUAL30.vcf.gz", shell=True, capture_output=True, text=True).stdout)
n_11 = int(subprocess.run(f"{ENV}/bcftools view -H -r chr11 {G1}/S01_S10.genome.QUAL30.vcf.gz | wc -l", shell=True, capture_output=True, text=True).stdout)
need(f"produced {c(n_raw)} records after splitting multi-allelic sites ({c(n_q)} with QUAL ≥ 30)", "G1 counts")
need(f"({c(n_11)} at QUAL ≥ 30)", "G1 chr11")
smp = subprocess.run(f"{ENV}/bcftools query -l {G1}/S01_S10.genome.QUAL30.vcf.gz", shell=True, capture_output=True, text=True).stdout.split()
if len(smp) != 10:
    fails.append(f"G1 VCF has {len(smp)} samples")
# G2 genome-wide DLD1
G2 = {r["metric"]: r["value"] for r in tsv(f"{D}/genome_wide_DLD1_comparison/summary_counts.tsv")}
ev = tsv(f"{D}/genome_wide_DLD1_comparison/dld1_genome_wide_events.tsv")
dx = tsv(f"{D}/genome_wide_DLD1_comparison/dld1_dox_specific_events.tsv")
if int(G2["events_total"]) != len(ev) or len(dx) != sum(r["comparison_S04_S05_vs_S01_S03"].startswith("dox_specific") for r in ev):
    fails.append("G2 summary/table counts disagree")
need(f"{c(len(ev))} DELLY events", "G2 total in report")
for k in ("dox_specific", "dox_specific_underpowered", "shared_with_baseline", "baseline_only", "not_present"):
    need(f"| `{k}` | {c(int(G2.get('comparison_' + k, 0)))} |", f"G2 {k}")
for r in dx:
    if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific":
        need(f"| {r['event_id']} |", "G2 event row")
        s = r["delly_present_samples"]
        need(f"| {r[s + '_alt']} / {r[s + '_ref']} |", f"G2 counts {r['event_id']}")
        if min(int(r[f"S0{i}_ref"]) for i in (1, 2, 3)) < 11 or any(int(r[f"S0{i}_alt"]) for i in (1, 2, 3)):
            fails.append(f"G2 {r['event_id']} baseline claim (0 variant, >= 11 reference) violated")
n_dx = sum(1 for r in dx if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific")
n_m = sum(1 for r in dx if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific" and r["manta_two_end_samples"])
n_p = sum(1 for r in dx if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific" and r["delly_somatic_filter"] == "PASS")
need(f"{n_dx} structural variants elsewhere in the genome", "G2 core n")
need(f"{n_m} of the {n_dx} are also called by Manta", "G2 Manta n")
need(f"5 pass DELLY's somatic filter, {n_m} are also called by Manta" if n_p == 5 else "@@", "G2 PASS n")
if "@@" in rep:
    fails.append("placeholder left in report")
# primers
pr = tsv(f"{D}/tables/primer_design.tsv")
if {r["event_id"] for r in pr} != {r["junction_id"] for r in ge2_all}:
    fails.append("primer table events != events with >= 2 fragments")
sel = [r for r in pr if r["design_status"].startswith("pair_selected")]
need(f"For two of the {len(ge2_all)} junctions with two or more fragments" if len(sel) == 2 else "@@", "primer count")
need(f"For the other {len(ge2_all) - len(sel)} junctions", "primer others")
for r in sel:
    need(r["left_primer"], "primer seq"); need(r["right_primer"], "primer seq")
    need(f"{r['product_bp_on_junction_template']} bp", "product size")
pb = tsv(f"{D}/tables/primer_binding_sites.tsv")
for r in pb:
    if r["pair"] == "selected" and r["within_1kb_of_event_breakend"] != "1":
        fails.append(f"selected primer binds away from its event: {r}")
# P42-style audits and boundary
for f in ("methods/header_sanitization_audit.tsv", "methods/delly_conversion_audit.tsv"):
    if any(r["status"] != "ok" for r in tsv(f"{D}/{f}")):
        fails.append(f"{f} has non-ok rows")
paths = 0
for b in glob.glob(f"{D}/**/*.bam", recursive=True):
    paths += subprocess.run(f"{ENV}/samtools view --no-PG -H {b} | grep -c '/home/\\|/mnt/\\|/Work_bio'", shell=True,
                            capture_output=True, text=True).stdout.strip() not in ("0", "")
for v in glob.glob(f"{D}/**/*.vcf.gz", recursive=True):
    paths += subprocess.run(f"{ENV}/bcftools view --no-version -h {v} 2>/dev/null | grep -c '/home/\\|/mnt/\\|/Work_bio'",
                            shell=True, capture_output=True, text=True).stdout.strip() not in ("0", "")
if paths:
    fails.append(f"{paths} BAM/VCF headers still contain local paths")
for r in tsv(f"{D}/tables/client_request_coverage.tsv"):
    for p in r["delivered_files"].split(";"):
        if p and not os.path.exists(f"{D}/{p}"):
            fails.append(f"coverage names missing {p}")
led = tsv(f"{D}/tables/key_junction_read_evidence.tsv")
if {r["event_id"] for r in led} != {r["junction_id"] for r in ge2_all}:
    fails.append("ledger events != events with >= 2 fragments")
inv = tsv(f"{D}/tables/data_inventory.tsv")
nfiles = sum(1 for f in glob.glob(f"{D}/**/*", recursive=True) if os.path.isfile(f) and not f.endswith("md5sums_delivery.txt"))
if len(inv) != nfiles:
    fails.append(f"inventory rows {len(inv)} != files {nfiles}")
for f in ("Fig1_locus_overview", "Fig2_rescue_sources", "Fig3_key_events_summary", "Fig5_before_after_and_classes", "Fig6_chr11_BAF"):   # 1027: report-order numbers
    for ext in ("png", "pdf"):
        if not os.path.exists(f"{D}/figures/{f}.{ext}"):
            fails.append(f"missing figure {f}.{ext}")
jv = sorted(glob.glob(f"{D}/figures/junction_views/*.png"))
if len(jv) != len(ge2_all):
    fails.append(f"junction views {len(jv)} != events with >=2 fragments {len(ge2_all)}")
# client text boundary: no internal project numbers, banned wording
for f in glob.glob(f"{D}/**/*", recursive=True):
    if f.endswith((".md", ".tsv", ".json", ".fa")) and os.path.isfile(f):
        if re.search(r"\bP4[0-9]J?\d*\b|P4[0-9]_", open(f, errors="ignore").read()):
            fails.append(f"internal project number in {os.path.relpath(f, D)}")
for bad in ("long-read", "long read", "Nanopore", "PacBio", "incorrect", "mistake", "wrong", "erroneous"):
    if bad.lower() in rep.lower():
        fails.append(f"banned wording present: {bad}")
if fails:
    print("SELF-AUDIT FAILED")
    print("\n".join(fails))
    sys.exit(1)
print(f"SELF-AUDIT OK: {len(T)} target + {len(W)} flank junctions, {len(jv)} junction views, G2 {len(ev)} events; "
      "all report numbers match delivered tables")
