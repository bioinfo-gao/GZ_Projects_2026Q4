#!/usr/bin/env bash
# 把 WALKUP-21683 的 20 个 FASTQ 补齐到 Dropbox 交付目录（用户 2026-11-21 要求：必须分享给客户）
# - 只放 20 个 FASTQ；Reports/SampleSheet/delivery_manifest 含别家客户样本编号，不放
# - 前提：全量校验日志以 ALL_DONE 结尾（size+crc32c 与 bucket 一致）
# - 先拷到 Dropbox 目录外的暂存区（同盘 /Work_bio），拷完 cmp 大小再 mv 进去 ⇒ Dropbox 不会上传半截文件
# - 已存在的同名文件一律不覆盖（Dropbox 里原有 8 个，已核对 crc32c 与 bucket 一致）
set -euo pipefail
SRC=/srv/external_data/gao_projects/46_park_wgs_raw/WALKUP-21683
DBX=/Work_bio/dropbox/Dropbox_Data/Athenomics_Delivery/WALKUP-21683/Yeji_WGS
STG=/Work_bio/dropbox/_staging_WALKUP-21683
VLOG=$SRC/_final_verify_20261121_1131.log
mkdir -p "$STG"

while tmux has-session -t gcs_final_verify 2>/dev/null; do sleep 60; done
tail -1 "$VLOG" | grep -q '^ALL_DONE$' || { echo "全量校验未通过，不拷贝：$VLOG"; exit 1; }
echo "[$(date '+%F %T')] 校验通过，开始补齐 Dropbox"

for f in "$SRC"/fastq/*WGS0000*.fastq.gz; do
  n=$(basename "$f")
  if [[ -e "$DBX/$n" ]]; then
    [[ $(stat -c %s "$DBX/$n") == $(stat -c %s "$f") ]] && echo "已存在且大小一致，跳过 $n" || echo "⚠ 已存在但大小不同，未动 $n"
    continue
  fi
  echo "[$(date '+%F %T')] 拷贝 $n"
  cp "$f" "$STG/$n"
  [[ $(stat -c %s "$STG/$n") == $(stat -c %s "$f") ]] || { echo "FAIL 大小不符 $n"; exit 1; }
  mv -n "$STG/$n" "$DBX/$n"
done
N=$(ls "$DBX"/*.fastq.gz | wc -l)
echo "[$(date '+%F %T')] Dropbox 本地目录 FASTQ 数：$N（应为 20）"
[[ $N -eq 20 ]] && echo DROPBOX_LOCAL_DONE || exit 1
