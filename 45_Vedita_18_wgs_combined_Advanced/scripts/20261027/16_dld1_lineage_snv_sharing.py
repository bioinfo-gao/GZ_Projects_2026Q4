"""Stage 16 — why the 18 distal-11q sites where S04 and S05 retained different alleles are all S04 ref / S05 alt.

Hypothesis: S02, S03 and S05 share a sub-lineage of the DLD1 population that S04 does not belong to; SNVs that arose
in that sub-lineage are heterozygous in S02, S03 and S05, subclonal in S01, absent from S04. Such SNVs lying on the
retained homolog pass the baseline-heterozygous filter of the BAF analysis (het in S01-S03) and appear as S04 ref /
S05 alt. S04-private SNVs can never pass that filter (absent from S02/S03), so the reverse combination is 0 by
construction.

Evidence unit: biallelic SNV record of the joint S01-S10 call (QUAL >= 30), allele depths per sample.
Kept: every DLD1 sample S01-S05 with DP >= 8; no alt read in any GM12878 sample (S06-S10), so the SNV is not a
germline variant shared with an unrelated line. Carrier = AF 0.25-0.75; non-carrier = 0 alt reads (DP >= 8).
Output: analysis/20261027_dld1_lineage/ (new directory).
"""
import collections
import os
import statistics as st
import subprocess

PROJ = "/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced"
BCF = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/bcftools"
VCF = f"{PROJ}/custom_research_report_20261026/genome_wide_SNV_indel_VCF/S01_S10.genome.QUAL30.vcf.gz"
BAF = f"{PROJ}/custom_research_report_20261026/supporting_analyses/chr11_BAF/baf_sites.tsv.gz"
OUT = f"{PROJ}/analysis/20261027_dld1_lineage"
CLONES = ["S02", "S03", "S04", "S05"]
os.makedirs(OUT, exist_ok=True)

# 1. the 18 discordant distal-11q sites, with allele depths in all ten samples
sites = []
for line in subprocess.run(["zcat", BAF], capture_output=True, text=True, check=True).stdout.splitlines():
    x = line.split("\t")
    if x[0] != "DLD1-136" or not (114_700_000 < int(x[1]) <= 135_000_000) or "" in x[5:7] or "NA" in x[5:7]:
        continue
    a, b = float(x[5]), float(x[6])
    if (a < 0.1 and b > 0.9) or (a > 0.9 and b < 0.1):
        sites.append((int(x[1]), a, b))
reg = f"{OUT}/_sites.tsv"
open(reg, "w").writelines(f"chr11\t{p}\n" for p, *_ in sites)
q = subprocess.run([BCF, "query", "-R", reg, "-f", "%POS\t%REF\t%ALT\t%QUAL[\t%AD]\n", VCF],
                   capture_output=True, text=True, check=True).stdout.splitlines()
os.remove(reg)
with open(f"{OUT}/distal11q_discordant_sites.tsv", "w") as fh:
    fh.write("pos\tref\talt\tqual\t" + "\t".join(f"S{i:02d}_AD" for i in range(1, 11))
             + "\tS04_retained\tS05_retained\n")
    for l in q:
        x = l.split("\t")
        p = int(x[0])
        s04, s05 = next((a, b) for pp, a, b in sites if pp == p)
        fh.write(l + f"\t{'alt' if s04 > 0.9 else 'ref'}\t{'alt' if s05 > 0.9 else 'ref'}\n")
assert len(q) == len(sites) == 18, (len(q), len(sites))

# 2. genome-wide sharing pattern of DLD1 clone SNVs absent from GM12878
PATTERNS = {"S02+S03+S05": (0, 1, 3), "S04 only": (2,), "S05 only": (3,), "S02+S03": (0, 1),
            "S02+S03+S04": (0, 1, 2), "S02+S04+S05": (0, 2, 3), "S03+S04+S05": (1, 2, 3), "S04+S05": (2, 3),
            "S02 only": (0,), "S03 only": (1,)}
rows = []
for c in [f"chr{i}" for i in range(1, 23)]:
    q = subprocess.run([BCF, "query", "-r", c, "-i", 'TYPE="snp"', "-f", "[%AD\t]\n", VCF],
                       capture_output=True, text=True, check=True).stdout.splitlines()
    cnt, s01 = collections.Counter(), collections.defaultdict(list)
    for l in q:
        try:
            ad = [tuple(map(int, a.split(",")))[:2] for a in l.rstrip("\t").split("\t")]
        except ValueError:
            continue
        if len(ad) != 10 or any(len(t) < 2 for t in ad) or any(sum(t) < 8 for t in ad[:5]) \
                or any(t[1] > 0 for t in ad[5:]):
            continue
        cl = ad[1:5]
        af = [v / (r + v) for r, v in cl]
        for name, carriers in PATTERNS.items():
            if all(0.25 <= af[i] <= 0.75 for i in carriers) and all(cl[i][1] == 0 for i in range(4) if i not in carriers):
                cnt[name] += 1
                s01[name].append(ad[0][1] / sum(ad[0]))
    for name in PATTERNS:
        rows.append((c, name, cnt[name], f"{st.median(s01[name]):.3f}" if s01[name] else "NA"))
    print(c, dict(cnt), flush=True)
with open(f"{OUT}/clone_snv_sharing_by_chrom.tsv", "w") as fh:
    fh.write("chrom\tcarrier_pattern\tn_snv\tS01_median_AF\n")
    fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
tot = collections.Counter()
for c, name, n, _ in rows:
    tot[name] += n
with open(f"{OUT}/clone_snv_sharing_autosomes.tsv", "w") as fh:
    fh.write("carrier_pattern\tn_snv_autosomes\n")
    fh.writelines(f"{k}\t{v}\n" for k, v in tot.most_common())
print("TOTAL", dict(tot.most_common()))
print("DONE", OUT)
