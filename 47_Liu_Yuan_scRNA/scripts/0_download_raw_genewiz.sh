#!/usr/bin/env bash
# P47 原始数据下载：GENEWIZ 项目 30-1371703644（来源：客户 2026-11-24 转发的 GENEWIZ 10-01 通知邮件；远端文件时间戳 2026-11-19）
# 用法：bash 0_download_raw_genewiz.sh <sftp_user>
# 密码：若存在 $PASSFILE（用户本人创建，chmod 600，git 仓库外）⇒ 经环境变量 LFTP_PASSWORD 交给 lftp（--env-password），
#       脚本不打印、不记录；否则 lftp 交互提示由用户输入。agent 不读取该文件内容。
# 落盘：原始数据放 /srv/external_data（与 P46 同惯例，SSD 留给分析产物），项目 raw_data/ 内为软链接。
# 目标目录名固定为 _1124（首次下载日），任何一天重跑都续传同一目录。
# 抗中断：net:max-retries 0（无限重连）+ repeat --until-ok（mirror 失败自动再来），均在同一 lftp 进程内，无需重输密码；
#         进程真退出时由 6_download_watch.sh 自动重启（需 $PASSFILE）。
set -euo pipefail
USER_NAME="${1:?usage: $0 <sftp_user>}"
PROJ=/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA
RAW=/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw
DEST="$RAW/genewiz_30-1371703644_1124"
PASSFILE=/home/gao/projects_2026H1/psd/psd.txt
LOG="$PROJ/logs/download_genewiz_$(date +%Y%m%d_%H%M).log"
LISTING="$RAW/remote_listing_genewiz_$(date +%Y%m%d_%H%M).txt"
mkdir -p "$DEST" "$PROJ/logs"
OPEN="open -u $USER_NAME sftp://sftp.genewiz.com"
if [ -s "$PASSFILE" ]; then
  LFTP_PASSWORD=$(<"$PASSFILE"); export LFTP_PASSWORD
  OPEN="open --env-password -u $USER_NAME sftp://sftp.genewiz.com"
  echo "[$(date '+%F %T')] password source: file (non-interactive)" | tee -a "$LOG"
fi
echo "[$(date '+%F %T')] dest=$DEST" | tee -a "$LOG"
# --parallel=4 同时 4 个文件，--use-pget-n=6 每个大文件拆 6 段；总带宽受本机出口 ~13 MiB/s 限制（1124 实测）
lftp -c "
set sftp:auto-confirm yes;
set net:max-retries 0; set net:reconnect-interval-base 30; set net:reconnect-interval-max 600; set net:timeout 120;
$OPEN;
find -l . > $LISTING;
repeat --until-ok -d 120 mirror --continue --parallel=4 --use-pget-n=6 --verbose . $DEST
" > >(tee -a "$LOG") 2>&1
unset LFTP_PASSWORD
echo "[$(date '+%F %T')] remote listing: $LISTING ($(wc -l < "$LISTING") lines)" | tee -a "$LOG"
echo "[$(date '+%F %T')] mirror finished; local file count: $(find "$DEST" -type f | wc -l)" | tee -a "$LOG"
find "$DEST" -iname '*.md5' -type f | while read -r m; do
  (cd "$(dirname "$m")" && md5sum -c "$(basename "$m")") 2>&1 | tee -a "$LOG" || true
done
echo "[$(date '+%F %T')] ALL_DONE" | tee -a "$LOG"
