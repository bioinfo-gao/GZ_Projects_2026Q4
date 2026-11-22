#!/usr/bin/env bash
# 2_build_starsolo_ref.sh —— P47 单细胞专用 STAR 索引（2026-11-22 建）
# 基因组：GENCODE v45 GRCh38 primary assembly（本机共享参考，见 /Work_bio/references/Structure_and_Readme.md）
# GTF：按 10x GRCh38-2024-A 公布的规则重现过滤（不是 10x 官方文件本身）：
#   - gene_type 与 transcript_type 都须在允许名单（protein_coding / lncRNA / IG_* / TR_* 及其 pseudogene）
#   - 去掉 tag "readthrough_transcript" 的转录本；去掉 *_PAR_Y 基因
#   为什么：全量 GTF 里重叠注释（readthrough 等）让 read 落到多个基因，STARsolo 会丢弃这些 read
#   注意：Mt_rRNA / Mt_tRNA 不在允许名单 ⇒ 不计数；chrM 序列仍完整在索引里，allele pileup 不受影响
# 比对器：STAR 2.7.11b（p47_align env）。sjdbOverhang 留默认 100（10x R2 常见 90–91 bp）。
# 线程：8（1122 实测 p46 同时占 ~40 线程，16 会使总量超 56 上限）。内存：p46 同时占约 77 GB ⇒ 用稀疏 SA（genomeSAsparseD 2，
#   建库/比对内存约减半，比对结果相同、稍慢），limitGenomeGenerateRAM 28 GB（超限 STAR 自行报错，不会 OOM 他人任务）。
set -euo pipefail
SRC=/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45
OUT=$SRC/starsolo_10x2024Afilter_star2711b
ENV=/mnt/ex_8T_SSD/conda_envs/p47_align
mkdir -p "$OUT"
GTF_OUT=$OUT/gencode.v45.primary_assembly.10x2024Afilter.gtf

if [ ! -s "$GTF_OUT" ]; then
python3 -I - "$SRC/gencode.v45.annotation.gtf" "$GTF_OUT" <<'PY'
import sys, re, collections
src, out = sys.argv[1:3]
ALLOW = {"protein_coding","protein_coding_LoF","lncRNA",
         "IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene",
         "IG_V_pseudogene","IG_J_pseudogene","IG_C_pseudogene",
         "TR_C_gene","TR_D_gene","TR_J_gene","TR_V_gene","TR_V_pseudogene","TR_J_pseudogene"}
def attr(s, k):
    m = re.search(k + r' "([^"]+)"', s); return m.group(1) if m else None
keep_tx, keep_gene = set(), set()
why = collections.Counter()
with open(src) as f:                      # 第一遍：决定保留哪些 transcript / gene
    for ln in f:
        if ln.startswith("#"): continue
        c = ln.rstrip("\n").split("\t")
        if c[2] != "transcript": continue
        a = c[8]; g = attr(a, "gene_id"); t = attr(a, "transcript_id")
        if g.endswith("_PAR_Y"): why["PAR_Y"] += 1; continue
        if attr(a, "gene_type") not in ALLOW: why["gene_type"] += 1; continue
        if attr(a, "transcript_type") not in ALLOW: why["transcript_type"] += 1; continue
        if 'tag "readthrough_transcript"' in a: why["readthrough"] += 1; continue
        keep_tx.add(t); keep_gene.add(g)
n_in = collections.Counter(); n_out = collections.Counter()
with open(src) as f, open(out, "w") as o:  # 第二遍：写出
    for ln in f:
        if ln.startswith("#"): o.write(ln); continue
        c = ln.split("\t"); n_in[c[2]] += 1
        a = c[8]
        ok = attr(a, "gene_id") in keep_gene if c[2] == "gene" else attr(a, "transcript_id") in keep_tx
        if ok: o.write(ln); n_out[c[2]] += 1
print("transcripts dropped by reason:", dict(why))
print("kept genes:", len(keep_gene), " kept transcripts:", len(keep_tx))
for k in ("gene", "transcript", "exon"): print(f"  {k}: {n_in[k]} -> {n_out[k]}")
PY
fi

if [ ! -s "$OUT/SA" ]; then
  "$ENV/bin/STAR" --runMode genomeGenerate --runThreadN 8 \
    --genomeDir "$OUT" --genomeFastaFiles "$SRC/GRCh38.primary_assembly.genome.fa" \
    --sjdbGTFfile "$GTF_OUT" --sjdbOverhang 100 \
    --genomeSAsparseD 2 --limitGenomeGenerateRAM 28000000000 --outFileNamePrefix "$OUT/build_"
fi
"$ENV/bin/STAR" --version
grep -E "versionGenome|sjdbGTFfile|sjdbOverhang" "$OUT/genomeParameters.txt"
echo REF_DONE
