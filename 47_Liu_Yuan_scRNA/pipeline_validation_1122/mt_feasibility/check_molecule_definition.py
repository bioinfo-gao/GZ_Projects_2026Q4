"""位点 12850：不同"分子"定义下的计数，对照 maegatk depth（27,450）与审计分子数（11,407）"""
import pysam, collections
cells = {l.strip() for l in open("cells.tsv")}
bam = pysam.AlignmentFile("pbmc_chrM.bam"); p0 = 12850 - 1
reads = 0; cbub = set(); cbub_start = set(); cbub_start_strand = set()
for r in bam.fetch("chrM", p0, p0 + 1):
    if r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
    if not (r.has_tag("CB") and r.has_tag("UB")) or r.get_tag("CB") not in cells: continue
    if p0 not in set(r.get_reference_positions()): continue
    reads += 1; k = (r.get_tag("CB"), r.get_tag("UB"))
    cbub.add(k); cbub_start.add(k + (r.reference_start,)); cbub_start_strand.add(k + (r.reference_start, r.is_reverse))
print(f"reads={reads}  (CB,UB)={len(cbub)}  (CB,UB,start)={len(cbub_start)}  (CB,UB,start,strand)={len(cbub_start_strand)}")
