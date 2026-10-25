"""Stage 12c (P45; P42's DELLY conversion audit) — every DELLY VCF in the delivery folder against its source BCF.
For genome_wide_SV_VCF/delly/<S>.delly.{unfiltered,PASS}.vcf.gz the source is the P41 discovery BCF (all records /
FILTER=PASS records); for genome_wide_DLD1_comparison/dld1_S01_S05.delly_genotyped.vcf.gz it is the G2a BCF.
Compared: record count and SHA-256 of the record lines (bcftools view -H; INFO/FORMAT text identical by construction,
header excluded because 12a rewrites paths there). Writes <delivery>/methods/delly_conversion_audit.tsv; exit 1 on mismatch.
Usage: 12c_delly_conversion_audit.py <delivery_dir>
"""
import glob
import hashlib
import os
import subprocess
import sys

D = sys.argv[1]
BCF = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/bcftools"
P41 = "/home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/sv/delly/discover"
G2 = "/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced/analysis/20261025_p45/G2_dld1_genome_wide"


def digest(args):
    h, n = hashlib.sha256(), 0
    p = subprocess.Popen([BCF, "view", "-H"] + args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in p.stdout:
        h.update(line)
        n += 1
    if p.wait():
        raise SystemExit(f"bcftools failed: {args}")
    return n, h.hexdigest()


pairs = []
for v in sorted(glob.glob(f"{D}/genome_wide_SV_VCF/delly/*.delly.*.vcf.gz")):
    s, kind = os.path.basename(v).split(".delly.")[0], os.path.basename(v).split(".delly.")[1].split(".")[0]
    pairs.append((v, f"{P41}/{s}.delly.bcf", ["-f", "PASS"] if kind == "PASS" else []))
for v in glob.glob(f"{D}/genome_wide_DLD1_comparison/*.vcf.gz"):
    pairs.append((v, f"{G2}/dld1_S01_S05.somatic_tagged.bcf", []))
rows, bad = [], 0
for v, src, flt in pairs:
    n0, h0 = digest(flt + [src])
    n1, h1 = digest([v])
    ok = n0 == n1 and h0 == h1
    bad += not ok
    rows.append((os.path.relpath(v, D), os.path.basename(src) + (" (FILTER=PASS)" if flt else ""), n0, n1, h0, h1,
                 "ok" if ok else "FAIL"))
os.makedirs(f"{D}/methods", exist_ok=True)
with open(f"{D}/methods/delly_conversion_audit.tsv", "w") as fh:
    fh.write("delivered_vcf\tsource_bcf\trecords_source\trecords_delivered\tbody_sha256_source\tbody_sha256_delivered\tstatus\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(f"12c: {len(rows)} DELLY files, {bad} mismatches")
sys.exit(1 if bad or not rows else 0)
