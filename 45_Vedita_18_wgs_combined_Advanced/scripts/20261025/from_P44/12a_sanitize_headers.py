"""Stage 12a (P1) — remove local absolute paths from the headers of every BAM / VCF in the delivery folder.
Runs on the DELIVERY COPIES only (sources are never touched). Each absolute path (/home, /mnt, /Work_bio, /tmp,
/data, /scratch) is replaced by its file name; nothing else in the header changes. Data body is proven unchanged:
record count and SHA-256 of the record stream (`samtools view` / `bcftools view -H`) before == after.
Writes <delivery>/methods/header_sanitization_audit.tsv; exit 1 if any body changed or any path remains.
Usage: 12a_sanitize_headers.py <delivery_dir>
"""
import glob
import hashlib
import os
import re
import subprocess
import sys

D = sys.argv[1]
ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"
SAM, BCF = f"{ENV}/samtools", f"{ENV}/bcftools"
PATH_RE = re.compile(r"(?:/home|/mnt|/Work_bio|/tmp|/data|/scratch)(?:/[^\s\t'\",;:]+)+")


def sh(cmd, binary=False):
    p = subprocess.run(cmd, shell=True, check=True, capture_output=True, executable="/bin/bash")
    return p.stdout if binary else p.stdout.decode()


def body_digest(cmd):
    h = hashlib.sha256()
    n = 0
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, executable="/bin/bash")
    for line in p.stdout:
        h.update(line)
        n += 1
    p.wait()
    if p.returncode:
        raise SystemExit(f"failed: {cmd}")
    return n, h.hexdigest()


def clean(text):
    k = len(PATH_RE.findall(text))
    return PATH_RE.sub(lambda m: os.path.basename(m.group(0)), text), k


rows, fail = [], False
for bam in sorted(glob.glob(f"{D}/**/*.bam", recursive=True)):
    hdr = sh(f"{SAM} view --no-PG -H {bam}")
    new, k = clean(hdr)
    n0, h0 = body_digest(f"{SAM} view {bam}")
    if k:
        hf = bam + ".hdr.tmp"
        open(hf, "w").write(new)
        tmp = bam + ".tmp"
        sh(f"{SAM} reheader --no-PG {hf} {bam} > {tmp} && mv {tmp} {bam} && rm -f {hf} && {SAM} index {bam}")
    n1, h1 = body_digest(f"{SAM} view {bam}")
    left = len(PATH_RE.findall(sh(f"{SAM} view --no-PG -H {bam}")))
    ok = n0 == n1 and h0 == h1 and left == 0
    fail |= not ok
    rows.append((os.path.relpath(bam, D), "BAM", k, n0, n1, h0, h1, left, "ok" if ok else "FAIL"))

for vcf in sorted(glob.glob(f"{D}/**/*.vcf.gz", recursive=True)):
    hdr = sh(f"{BCF} view --no-version -h {vcf} 2>/dev/null")
    new, k = clean(hdr)
    n0, h0 = body_digest(f"{BCF} view -H {vcf} 2>/dev/null")
    if k:
        hf = vcf + ".hdr.tmp"
        open(hf, "w").write(new)
        tmp = vcf + ".tmp.vcf.gz"
        sh(f"{BCF} reheader -h {hf} -o {tmp} {vcf} && mv {tmp} {vcf} && rm -f {hf} {vcf}.tbi {vcf}.csi && {BCF} index -f -t {vcf}")
    n1, h1 = body_digest(f"{BCF} view -H {vcf} 2>/dev/null")
    left = len(PATH_RE.findall(sh(f"{BCF} view --no-version -h {vcf} 2>/dev/null")))
    ok = n0 == n1 and h0 == h1 and left == 0
    fail |= not ok
    rows.append((os.path.relpath(vcf, D), "VCF", k, n0, n1, h0, h1, left, "ok" if ok else "FAIL"))

os.makedirs(f"{D}/methods", exist_ok=True)
with open(f"{D}/methods/header_sanitization_audit.tsv", "w") as fh:
    fh.write("file\ttype\tpaths_removed\trecords_before\trecords_after\tbody_sha256_before\tbody_sha256_after\t"
             "paths_remaining\tstatus\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(f"header sanitisation: {len(rows)} files, {sum(r[2] for r in rows)} paths removed, "
      f"{sum(r[8] == 'FAIL' for r in rows)} failed")
sys.exit(1 if fail else 0)
