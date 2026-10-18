#!/usr/bin/env bash
# P42 Step 3 —— 把项目 41 的 18 个 markdup BAM + 索引【复制】进 ~/Dropbox/Vedita_Bam（只准备，不分享、不发送）。
#
# 设计（每条都有理由）：
#  1) 复制而不是移动：源文件留在 项目 41 的 analysis/bam/，不动。
#  2) 先复制到与 Dropbox 同一文件系统的 staging（/Work_bio/dropbox/staging_...，在 Dropbox 同步目录之外），
#     源与副本分别算 MD5 并逐个比对，全部一致才 mv 进 ~/Dropbox/Vedita_Bam（同文件系统 ⇒ 瞬间重命名）。
#     这样 Dropbox 只会同步【已校验完整】的文件，不会同步写到一半的 BAM。
#  3) 只交付 *.markdup.bam + *.markdup.bam.bai（18 + 18）。不交付 *.sorted.bam（未标重复的中间产物，706 GB）。
#  4) 写 MD5.txt（相对文件名，客户用 md5sum -c 即可）与 README.txt（英文，不含我方本机路径）。
#  5) 线程：rsync 1 个 + md5sum 并行 4 个，远低于 28 核上限。
#  6) 不覆盖：若 ~/Dropbox/Vedita_Bam 已有同名文件则中止。
#
# 用法: bash 3_prepare_bam_for_dropbox.sh        # 自动进 tmux 会话 p42_bam_dropbox
set -euo pipefail

SRC=/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/analysis/bam
STAGE=/Work_bio/dropbox/staging_Vedita_Bam_20261014
DEST="$HOME/Dropbox/Vedita_Bam"
SAMTOOLS=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools
LOG=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis/logs
mkdir -p "$LOG"

if [[ -z "${TMUX:-}" ]]; then
  tmux new-session -d -s p42_bam_dropbox "bash '$0' 2>&1 | tee -a '$LOG/3_prepare_bam_for_dropbox_$(date +%Y%m%d_%H%M).log'"
  echo "已在 tmux 会话 p42_bam_dropbox 启动；查看: tmux capture-pane -t p42_bam_dropbox -p | tail -20"; exit 0
fi

say(){ echo "[$(date '+%F %T')] $*"; }
die(){ say "FATAL: $*"; exit 1; }

mapfile -t BAMS < <(ls "$SRC"/*.markdup.bam | sort)
[[ ${#BAMS[@]} -eq 18 ]] || die "markdup.bam 数量 ${#BAMS[@]} != 18"

say "STEP3_START  preflight"
for b in "${BAMS[@]}"; do
  sid=$(basename "$b" .markdup.bam)
  [[ -e "$SRC/$sid.markdup.done" ]] || die "$sid: 无 .markdup.done"
  [[ -s "$b.bai" ]] || die "$sid: 无 .bai"
  "$SAMTOOLS" quickcheck -v "$b" || die "$sid: samtools quickcheck 失败"
done
need=$(du -cb "${BAMS[@]}" "${BAMS[@]/%/.bai}" | tail -1 | cut -f1)
free=$(df -B1 --output=avail /Work_bio | tail -1)
say "需要 $((need/1024**3)) GiB，/Work_bio 可用 $((free/1024**3)) GiB"
[[ $free -gt $((need*12/10)) ]] || die "空间不足（要求 ≥ 1.2×）"

mkdir -p "$STAGE" "$DEST"
if compgen -G "$DEST/*" > /dev/null; then die "$DEST 非空，中止（不覆盖）"; fi

say "COPY_START -> $STAGE"
rsync -a --info=stats1 "${BAMS[@]}" "${BAMS[@]/%/.bai}" "$STAGE"/
say "COPY_DONE"

say "MD5_START  (源与副本各 4 并行)"
( cd "$SRC"   && ls *.markdup.bam *.markdup.bam.bai | sort | xargs -P 4 -n 1 md5sum | sort -k2 ) > "$STAGE/.src.md5"
( cd "$STAGE" && ls *.markdup.bam *.markdup.bam.bai | sort | xargs -P 4 -n 1 md5sum | sort -k2 ) > "$STAGE/.dst.md5"
cmp "$STAGE/.src.md5" "$STAGE/.dst.md5" || die "源与副本 MD5 不一致"
[[ $(wc -l < "$STAGE/.dst.md5") -eq 36 ]] || die "MD5 行数 != 36"
say "MD5_OK  36 个文件源与副本一致"
cp "$STAGE/.dst.md5" "$STAGE/MD5.txt"; rm -f "$STAGE/.src.md5" "$STAGE/.dst.md5"

cat > "$STAGE/README.txt" <<'EOT'
BAM files for the 18 WGS samples (Athenomics)

Contents
  <sample>.markdup.bam      coordinate-sorted alignments, duplicates marked (not removed)
  <sample>.markdup.bam.bai  index for the BAM file
  MD5.txt                   MD5 checksum of every file above (36 lines)
  Sample names: S01_DLD1_pop ... S18_DG75_c14 (S01-S05 DLD1, S06-S10 GM12878, S11-S14 Mutu I, S15-S18 DG75).

How the files were made
  Reference: GRCh38 (hg38 coordinates), GATK analysis set, including ALT, decoy and HLA contigs and chrEBV.
  Aligner: bwa-mem2 2.2.1 (mem -Y), reads not trimmed before alignment.
  Duplicates: marked with GATK MarkDuplicates 4.6.2.0.
  Read group: SM tag = sample name.

Checking the files after download
  Linux/macOS:  md5sum -c MD5.txt
  Windows:      certutil -hashfile <file> MD5   and compare with the matching line in MD5.txt
  The BAM files are large (18 files, about 845 GiB in total); please keep at least 900 GiB free.
EOT

say "MOVE_START -> $DEST (同文件系统重命名)"
mv "$STAGE"/*.markdup.bam "$STAGE"/*.markdup.bam.bai "$STAGE"/MD5.txt "$STAGE"/README.txt "$DEST"/
rmdir "$STAGE" 2>/dev/null || say "staging 目录非空，保留"
n=$(ls "$DEST" | wc -l)
say "DONE  $DEST 内 $n 个文件（预期 38 = 18 BAM + 18 bai + MD5.txt + README.txt）"
[[ $n -eq 38 ]] || die "文件数 != 38"
( cd "$DEST" && md5sum -c MD5.txt --quiet && say "FINAL md5sum -c MD5.txt 全部通过" )
