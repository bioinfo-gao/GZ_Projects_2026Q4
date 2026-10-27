"""Stage 15 — independent check of the two largest +Dox-specific DELLY events of Section 7.9
(DEL00003167 chr4 888 kb in S05; DUP00005675 chr8 327 kb in S04).

Question: does the copy number and allele balance of each interval change in the clone that carries the call?
Evidence layers, none of them DELLY's own counts:
  1. read depth (samtools bedcov, MAPQ>=20, deletions excluded) in bins, each sample normalised to its own flank
     median, then divided by the mean of S01-S03;
  2. allele fraction at SNPs heterozygous in all of S01-S03 (joint S01-S10 call, QUAL>=30);
  3. reads at the chr4 junction: split fragments (SA to the other breakend), reference-spanning fragments.
Output: analysis/20261027_sv_verification/*.tsv (new directory; nothing earlier is rewritten).
"""
import collections
import os
import re
import statistics as st
import subprocess
import sys

sys.path.insert(0, "/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced/scripts/20261025")
import p45lib as L  # noqa: E402

ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"
VCF = f"{L.PROJ}/custom_research_report_20261026/genome_wide_SNV_indel_VCF/S01_S10.genome.QUAL30.vcf.gz"
OUT = f"{L.PROJ}/analysis/20261027_sv_verification"
BAM = {r["sample_id"]: r["bam_path"] for r in L.samples()}
DLD1 = ["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl", "S04_DLD1_c3_dox", "S05_DLD1_c4_dox"]
BG = ["S06_GM12878_pop", "S11_MutuI_pop", "S15_DG75_pop"]
EVENTS = {  # id: chrom, start, end, carrier, coarse bin, padding
    "DEL00003167": ("chr4", 156284227, 157172208, "S05_DLD1_c4_dox", 50000, 2000000),
    "DUP00005675": ("chr8", 46747734, 47074827, "S04_DLD1_c3_dox", 25000, 1000000),
}
os.makedirs(OUT, exist_ok=True)


def bedcov(sample, bins):
    bed = f"{OUT}/_bins.bed"
    with open(bed, "w") as fh:
        fh.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in bins)
    res = subprocess.run([f"{ENV}/samtools", "bedcov", "-Q", "20", "-j", bed, BAM[sample]],
                         capture_output=True, text=True, check=True).stdout.split("\n")
    out = [int(l.split("\t")[3]) / (int(l.split("\t")[2]) - int(l.split("\t")[1])) for l in res if l]
    assert len(out) == len(bins), (sample, len(out), len(bins))
    return out


def depth_layer(eid):
    c, a, b, carrier, w, pad = EVENTS[eid]
    lo = max(a - pad, 0) // w * w
    bins = [(c, x, x + w) for x in range(lo, b + pad, w)]
    samples = DLD1 + BG
    d = {s: bedcov(s, bins) for s in samples}
    inside = [i for i, (_, s0, e0) in enumerate(bins) if s0 >= a and e0 <= b]
    flank = [i for i, (_, s0, e0) in enumerate(bins) if (e0 < a - w or s0 > b + w) and d["S01_DLD1_pop"][i] > 5]
    fl = {s: st.median(d[s][i] for i in flank) for s in samples}
    ratio = {s: st.median(d[s][i] for i in inside) / fl[s] for s in samples}
    base = st.mean(ratio[s] for s in DLD1[:3])
    with open(f"{OUT}/{eid}_depth_summary.tsv", "w") as fh:
        fh.write("sample\tflank_median_depth\tinside_median_depth\tinside_over_flank\tvs_mean_S01_S03\tn_inside_bins\n")
        for s in samples:
            fh.write(f"{s}\t{fl[s]:.1f}\t{fl[s] * ratio[s]:.1f}\t{ratio[s]:.3f}\t{ratio[s] / base:.3f}\t{len(inside)}\n")
    with open(f"{OUT}/{eid}_depth_bins.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\tinside_event\t" + "\t".join(f"{s}_norm" for s in samples)
                 + f"\t{carrier}_vs_S01_S03\n")
        for i, (cc, s0, e0) in enumerate(bins):
            n = {s: d[s][i] / fl[s] for s in samples}
            den = st.mean(n[s] for s in DLD1[:3])
            r = f"{n[carrier] / den:.3f}" if den > 0 else "NA"
            fh.write(f"{cc}\t{s0}\t{e0}\t{int(i in inside)}\t" + "\t".join(f"{n[s]:.3f}" for s in samples) + f"\t{r}\n")
    return ratio[carrier] / base


def fine_depth(eid, centre, half, w, tag):
    c, a, b, carrier, *_ = EVENTS[eid]
    bins = [(c, x, x + w) for x in range(centre - half, centre + half, w)]
    norm = {s: float(open(f"{OUT}/{eid}_depth_summary.tsv").read().split(s + "\t")[1].split("\t")[0]) for s in DLD1}
    d = {s: [x / norm[s] for x in bedcov(s, bins)] for s in DLD1}
    with open(f"{OUT}/{eid}_fine_depth_{tag}.tsv", "w") as fh:
        fh.write("chrom\tstart\tend\toffset_from_breakend\t" + "\t".join(f"{s}_norm" for s in DLD1)
                 + f"\t{carrier}_vs_S01_S03\n")
        for i, (cc, s0, e0) in enumerate(bins):
            den = st.mean(d[s][i] for s in DLD1[:3])
            r = f"{d[carrier][i] / den:.3f}" if den > 0 else "NA"
            fh.write(f"{cc}\t{s0}\t{e0}\t{s0 - centre}\t" + "\t".join(f"{d[s][i]:.3f}" for s in DLD1) + f"\t{r}\n")


def allele_layer(eid):
    c, a, b, carrier, w, pad = EVENTS[eid]
    q = subprocess.run([f"{ENV}/bcftools", "query", "-r", f"{c}:{max(a - pad, 1)}-{b + pad}", "-i", 'TYPE="snp"',
                        "-s", ",".join(DLD1), "-f", "%POS[\t%AD]\n", VCF], capture_output=True, text=True, check=True)
    dev = {"inside": collections.defaultdict(list), "flank": collections.defaultdict(list)}
    for line in q.stdout.splitlines():
        x = line.split("\t")
        try:
            ad = [tuple(map(int, s.split(",")))[:2] for s in x[1:]]
        except ValueError:
            continue
        if any(len(t) < 2 for t in ad) or not all(r + v >= 10 and 0.3 <= v / (r + v) <= 0.7 for r, v in ad[:3]):
            continue
        k = "inside" if a < int(x[0]) < b else "flank"
        for s, (r, v) in zip(DLD1, ad):
            if r + v >= 15:
                dev[k][s].append(abs(v / (r + v) - 0.5))
    with open(f"{OUT}/{eid}_allele_balance.tsv", "w") as fh:
        fh.write("region\tsample\tn_baseline_het_snps_DP15\tmedian_abs_AF_minus_0.5\tfrac_abs_dev_0.10_0.30\t"
                 "frac_abs_dev_gt_0.40\n")
        for k in ("inside", "flank"):
            for s in DLD1:
                v = dev[k][s]
                fh.write(f"{k}\t{s}\t{len(v)}\t{st.median(v):.3f}\t{sum(.1 <= x <= .3 for x in v) / len(v):.3f}\t"
                         f"{sum(x > .4 for x in v) / len(v):.3f}\n")
    return st.median(dev["inside"][carrier]), st.median(dev["flank"][carrier])


def junction_layer():
    """Reads at the two chr4 breakends: split fragments joining them, and reference-spanning fragments."""
    ends = [("A", 156284229, 156283900, 156284550), ("B", 157172208, 157171900, 157172550)]
    rows = []
    frag = []
    for s in DLD1:
        for lab, bp, r0, r1 in ends:
            sam = subprocess.run([f"{ENV}/samtools", "view", BAM[s], f"chr4:{r0}-{r1}"],
                                 capture_output=True, text=True, check=True).stdout.splitlines()
            split, ref, adapt = set(), set(), set()
            for l in sam:
                x = l.split("\t")
                fl = int(x[1])
                if fl & (4 | 1024):
                    continue
                m = re.search(r"SA:Z:(\S+)", l)
                other = 157172208 if lab == "A" else 156284229
                if m and any(p.split(",")[0] == "chr4" and abs(int(p.split(",")[1]) - other) < 600
                             for p in m.group(1).split(";") if p):
                    split.add(x[0])
                    if any(ad[:13] in x[9] or L.revcomp(ad[:13]) in x[9] for ad in L.ADAPTERS.values()):
                        adapt.add(x[0])
                    if s == "S05_DLD1_c4_dox":
                        frag.append((s, lab, x[0], fl, x[3], x[4], x[5], x[8], m.group(1)))
                if fl & 256 or int(x[4]) < 20:
                    continue
                pos, end = int(x[3]), int(x[3]) - 1 + sum(int(n) for n, o in re.findall(r"(\d+)([MDN=X])", x[5]))
                if pos <= bp - 20 and end >= bp + 20:
                    ref.add(x[0])
            rows.append((s, lab, bp, len(split), len(adapt), len(ref)))
    with open(f"{OUT}/DEL00003167_junction_reads.tsv", "w") as fh:
        fh.write("sample\tbreakend\tposition\tsplit_fragments_to_other_breakend\tof_which_contain_adapter\t"
                 "reference_spanning_fragments_MAPQ20\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    with open(f"{OUT}/DEL00003167_S05_split_records.tsv", "w") as fh:
        fh.write("sample\tbreakend_window\tread_name\tflag\tpos\tmapq\tcigar\ttlen\tSA\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in sorted(frag, key=lambda r: (r[2], r[1], r[3])))
    return rows


if __name__ == "__main__":
    for eid in EVENTS:
        print(eid, "carrier depth vs S01-S03 inside interval:", f"{depth_layer(eid):.3f}")
        print(eid, "carrier median |AF-0.5| inside / flank:", "%.3f / %.3f" % allele_layer(eid))
    fine_depth("DEL00003167", 156284229, 15000, 1000, "breakendA_1kb")
    fine_depth("DEL00003167", 157172208, 15000, 1000, "breakendB_1kb")
    fine_depth("DEL00003167", 157172208, 8000, 250, "breakendB_250bp")
    for r in junction_layer():
        print("junction", *r)
    os.remove(f"{OUT}/_bins.bed")
    print("DONE", OUT)
