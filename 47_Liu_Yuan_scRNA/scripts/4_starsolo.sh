#!/usr/bin/env bash
# 4_starsolo.sh —— 单样本 STARsolo 定量 + 带 CB/UB 标签的 BAM（2026-11-22 建）
# 用法: bash 4_starsolo.sh <sample> <chem v2|v3|v4> <outdir> <R1,R1,...> <R2,R2,...> [threads=8]
#   chem 由 3_detect_chemistry.py 判定；多 lane 用逗号拼接，R1/R2 顺序须一一对应。
# 参数取舍（STARsolo 文档中"复现 Cell Ranger"的设置）：
#   CB 1MM_multi_Nbase_pseudocounts / UMI 1MM_CR / MultiGeneUMI_CR / clipAdapterType CellRanger4 / EmptyDrops_CR
#   soloFeatures: Gene（外显子，= CR <7）+ GeneFull_Ex50pAS（含内含子，≈ CR ≥7 默认）+ SJ + Velocyto（spliced/unspliced，供轨迹）
#   soloMultiMappers: Unique（主矩阵）+ EM（多比对基因另出一份，做敏感性）
# BAM：保留全部 secondary alignment（STAR 默认 outSAMmultNmax -1，NH/HI 标签），未比对 read 留在 BAM（outSAMunmapped Within）
#   ⇒ 规则 §2A：secondary 不是独立分子，但是定位歧义（NUMT / 同源区）的直接证据，不在这一层丢。
set -euo pipefail
S=$1; CHEM=$2; OUT=$3; R1=$4; R2=$5; T=${6:-8}
ENV=/mnt/ex_8T_SSD/conda_envs/p47_align
IDX=/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/starsolo_10x2024Afilter_star2711b
WLD=/Work_bio/references/10x_whitelists_ngs_tools_1.8.6
case $CHEM in
  v2) WL=$WLD/10x_version2_whitelist.txt.gz; UMILEN=10 ;;
  v3) WL=$WLD/10x_version3_whitelist.txt.gz; UMILEN=12 ;;
  v4) WL=$WLD/10x_version4_whitelist.txt.gz; UMILEN=12 ;;
  *) echo "unknown chem $CHEM" >&2; exit 2 ;;
esac
mkdir -p "$OUT/$S"
WLTXT=$OUT/$S/whitelist_${CHEM}.txt
[ -s "$WLTXT" ] || zcat "$WL" > "$WLTXT"          # STAR 不读 .gz 白名单
[ -s "$IDX/SA" ] || { echo "index not built: $IDX" >&2; exit 3; }

"$ENV/bin/STAR" --runThreadN "$T" --genomeDir "$IDX" \
  --readFilesIn "$R2" "$R1" --readFilesCommand zcat \
  --outFileNamePrefix "$OUT/$S/" \
  --soloType CB_UMI_Simple --soloCBwhitelist "$WLTXT" \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen "$UMILEN" --soloBarcodeReadLength 0 \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --clipAdapterType CellRanger4 --outFilterScoreMin 30 \
  --soloCellFilter EmptyDrops_CR \
  --soloFeatures Gene GeneFull_Ex50pAS SJ Velocyto --soloMultiMappers Unique EM \
  --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --limitBAMsortRAM 20000000000 --soloOutFormatFeaturesGeneField3 "Gene Expression"
"$ENV/bin/samtools" index -@ 4 "$OUT/$S/Aligned.sortedByCoord.out.bam"
# 与 Cell Ranger 同格式：矩阵/条码/特征 gzip（scanpy read_10x_mtx / Seurat Read10X 直接可读）；无损、原位替换本次新产物
find "$OUT/$S/Solo.out" -type f \( -name '*.mtx' -o -name 'barcodes.tsv' -o -name 'features.tsv' \) -exec "$ENV/bin/pigz" -p 4 {} +
rm -f "$WLTXT"   # 本脚本自己刚解压的临时白名单（与 .gz 内容相同）
echo "STARSOLO_DONE $S"
