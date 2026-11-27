#!/usr/bin/env bash
# P46 Stage 3 v3b（2026-11-27）—— 同一口径回贴 pMSCV 单载体索引，使 Fig 6 两个 panel 方法对称。
# 另一作用：Cas9/EGFP 在 pLPC 里无同源序列，v2(双载体) 与 v3b(单载体) 结果应当一致 ⇒ 互为一致性检查。
#
# 为什么：pMSCV 8407-10390 与 pLPC 1-1973（ori+AmpR）99.1% 同源，双载体索引下骨架读段一律 MAPQ=0，
#   而深度用 `samtools depth -Q 20` ⇒ 全部丢弃 ⇒ pLPC 骨架恒为 0.00，与骨架是否真的存在无关。
#   实测：pLPC AmpR 区有 46–1,763 条读段，MAPQ≥20 的为 0 条。
# 口径：输入沿用 v2 已保存的候选读段（不重新从 BAM 提取，保证与 v2 可比）；v2 产物一律不改不删。
set -uo pipefail
source /Work_bio/gao/miniforge3/etc/profile.d/conda.sh; conda activate regular_bioinfo
PROJ=/home/gao/projects_2026Q4/46_park_wgs
V2=$PROJ/analysis/stage3_vector_v2_20261124
OUT=${OUT:-$PROJ/analysis/stage3_pMSCVonly_v3_20261127}
NPAR=${NPAR:-5}; THR=${THR:-5}          # 5×5=25 线程，远低于 28 核/56 线程上限
say(){ echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT"/{ref,bam,tables,logs,flags}

# ── 单载体参考：报告载体 pLPC 自己一条，没有同源竞争者 ──────────────────
if [[ ! -s "$OUT/ref/pMSCV_only.fa.ann" ]]; then
  samtools faidx "$PROJ/analysis/stage3_refs_20261121/vectors.fa" pMSCV_Cas9_EGFP > "$OUT/ref/pMSCV_only.fa"
  samtools faidx "$OUT/ref/pMSCV_only.fa"
  bwa-mem2 index "$OUT/ref/pMSCV_only.fa" > "$OUT/logs/index.log" 2>&1 || { say "FATAL index"; exit 1; }
fi
[[ -s "$OUT/ref/pMSCV_only.fa.ann" ]] || { say "FATAL 索引未生成"; exit 1; }

one(){
  local s=$1
  [[ -e "$OUT/flags/$s.done" ]] && { say "SKIP_CACHED $s"; return 0; }
  local cand="$V2/seq/$s.candidates.fa.gz"
  [[ -s "$cand" ]] || { say "FATAL 缺输入 $cand"; return 1; }
  zcat "$cand" \
  | bwa-mem2 mem -t "$THR" -v 1 "$OUT/ref/pMSCV_only.fa" - 2> "$OUT/logs/$s.bwa.log" \
  | samtools view -b -F 0x904 - 2> "$OUT/logs/$s.view.log" \
  | samtools sort -@ 2 -o "$OUT/bam/$s.pMSCVonly.bam" - 2> "$OUT/logs/$s.sort.log"
  local rc=("${PIPESTATUS[@]}")
  for i in "${!rc[@]}"; do [[ "${rc[$i]}" -eq 0 ]] || { say "FATAL $s pipe[$i] rc=${rc[$i]}"; return 1; }; done
  samtools index "$OUT/bam/$s.pMSCVonly.bam" || { say "FATAL index $s"; return 1; }
  # 与 v2 完全相同的深度口径，唯一差别是参考里没有同源竞争者
  samtools depth -a -Q 20 -J "$OUT/bam/$s.pMSCVonly.bam" > "$OUT/tables/$s.pMSCVonly_depth.txt"
  [[ -s "$OUT/tables/$s.pMSCVonly_depth.txt" ]] || { say "FATAL $s 深度为空"; return 1; }
  # 审计：总读段 / 贴上 / MAPQ>=20
  { echo -e "metric\tvalue"
    echo -e "mapped_primary\t$(samtools view -c "$OUT/bam/$s.pMSCVonly.bam")"
    echo -e "mapq_ge20\t$(samtools view -c -q 20 "$OUT/bam/$s.pMSCVonly.bam")"
  } > "$OUT/tables/$s.audit.tsv"
  touch "$OUT/flags/$s.done"; say "DONE $s"
}

mapfile -t S < <(python3 -c "import csv;[print(r['safe_name']) for r in csv.DictReader(open('$PROJ/scripts/0_sample_table.csv'))]")
say "样本数=${#S[@]} 并行=$NPAR 线程/样本=$THR"
pids=()
for s in "${S[@]}"; do
  one "$s" & pids+=($!)
  while (( $(jobs -rp | wc -l) >= NPAR )); do sleep 5; done
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
(( fail == 0 )) || { say "FATAL 有样本失败"; exit 1; }
n=$(ls "$OUT/flags"/*.done 2>/dev/null | wc -l)
[[ "$n" -eq "${#S[@]}" ]] || { say "FATAL 完成 $n/${#S[@]}"; exit 1; }
say "V3_MAP_ALL_DONE $n/${#S[@]}"
