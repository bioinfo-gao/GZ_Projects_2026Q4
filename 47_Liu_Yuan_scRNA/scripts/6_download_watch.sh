#!/usr/bin/env bash
# P47 下载看门狗（user directive 2026-11-24：「万一中断请重启，我不可能睡觉时候来工作」⇒ 授权自动重启下载）。
# 每小时把状态写入 logs/download_watch.log；事件同时写 ALERTS + 发邮件（wd_sendmail.py，自带 SMTP）+ 广播 tty。
#   ① lftp 退出但最新下载日志无 ALL_DONE ⇒ 等 5 min 后自动重启下载（--continue 续传）
#   ② lftp 在、但连续 2 次巡检（~2 h）网卡接收 < 1 MiB/s ⇒ 判为卡死，kill lftp 并重启
#   ③ 下载 + md5 校验全部结束 ⇒ 完成汇总邮件，写 COMPLETE 标记，退出
#   重启最多 MAX_RESTART 次；超出 ⇒ 写 GAVEUP 标记、发信、退出（cron 守护见到标记不再拉起）
# 自动重启前提：用户本人创建的密码文件 $PASSFILE 存在（本脚本只检查其存在与非空，不读内容）；
#   不存在 ⇒ 只告警不重启（重启会停在密码提示上空等）。
# 速度口径：网卡 rx_bytes 增量（文件大小/磁盘块在 pget 分段下均不可信，见 1124 实测）。
set -uo pipefail
PROJ=/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA
DEST=/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124
DL_SESS=p47_download
DL_CMD="cd $PROJ && bash scripts/0_download_raw_genewiz.sh hyu3_med_miami; exec bash"
PASSFILE=/home/gao/projects_2026H1/psd/psd.txt
WLOG="$PROJ/logs/download_watch.log"
ALERTS="$PROJ/logs/download_watch_ALERTS.log"
COMPLETE="$PROJ/logs/download_COMPLETE"
GAVEUP="$PROJ/logs/download_GAVEUP"
MAILER=/home/gao/Utilities/Server/monitoring_alerting/wd_sendmail.py
POLL="${POLL:-3600}"; MAX_RESTART="${MAX_RESTART:-12}"
EXPECT_FILES=380; EXPECT_MD5=8; EXPECT_GIB=996.7
DEV=$(ip route get 52.5.58.88 | grep -o 'dev [^ ]*' | cut -d' ' -f2)
rx() { cat /sys/class/net/$DEV/statistics/rx_bytes; }
dllog() { ls -t "$PROJ"/logs/download_genewiz_*.log | head -1; }
log() { echo "[$(date '+%F %T')] $*" | tee -a "$WLOG"; }
alert() { # $1=subject $2=body
  echo "[$(date '+%F %T')] $1 | $2" >> "$ALERTS"
  printf '%s\n\n%s\n' "$2" "$(tail -20 "$WLOG")" | python3 "$MAILER" "$1" - || log "MAIL FAILED: $1"
  for t in $(ps -eo tty= | grep pts | sort -u); do echo "[P47 download] $1" > /dev/$t 2>/dev/null; done
}
restarts=$(grep -c 'RESTART #' "$WLOG" 2>/dev/null || true); restarts=${restarts:-0}
restart_dl() { # $1=reason
  if [ ! -s "$PASSFILE" ]; then
    [ -z "${nopass_alerted:-}" ] && alert "🚨[P47] 下载中断，无密码文件无法自动重启" "$1。需在 tmux $DL_SESS 手动重跑并输入密码。"
    nopass_alerted=1; return 1; fi
  if [ "$restarts" -ge "$MAX_RESTART" ]; then
    : > "$GAVEUP"; alert "🚨[P47] 下载已自动重启 ${restarts} 次仍失败，停止重试" "$1"; exit 1; fi
  restarts=$((restarts+1))
  tmux kill-session -t "$DL_SESS" 2>/dev/null; sleep 300
  tmux new-session -d -s "$DL_SESS" "$DL_CMD"
  log "RESTART #$restarts ($1)"
  alert "⚠[P47] 下载中断，已自动重启（第 $restarts 次）" "$1"
}
low=0
log "watch start; dev=$DEV poll=${POLL}s restarts_so_far=$restarts"
while true; do
  [ -f "$COMPLETE" ] && { log "COMPLETE marker present, exit"; exit 0; }
  a=$(rx); t0=$(date +%s); sleep 60; b=$(rx); t1=$(date +%s)
  rate=$(python3 -c "print(f'{($b-$a)/($t1-$t0)/2**20:.1f}')")
  nf=$(find "$DEST" -type f ! -name '*.lftp-pget-status' | wc -l)
  disk=$(du -s --block-size=1G "$DEST" | cut -f1)
  L=$(dllog)
  if pgrep -x lftp >/dev/null; then
    log "RUNNING rx=${rate} MiB/s files=${nf}/${EXPECT_FILES} on_disk~${disk}G/${EXPECT_GIB}G"
    if python3 -c "import sys; sys.exit(0 if $rate < 1 else 1)"; then low=$((low+1)); else low=0; fi
    if [ $low -ge 2 ]; then low=0; pkill -x lftp; sleep 10; restart_dl "lftp 在但连续 2 次巡检网卡接收 <1 MiB/s（当前 ${rate}），判为卡死"; fi
  elif pgrep -f 'md5sum -c' >/dev/null; then
    log "MD5 checking (lftp finished) files=${nf}"
  elif grep -q 'ALL_DONE' "$L"; then
    ok=$(grep -c ': OK$' "$L"); bad=$(grep -c ': FAILED' "$L")
    msg="下载结束。本地文件 ${nf}（远端清单 ${EXPECT_FILES}），md5 OK ${ok}/${EXPECT_MD5}，FAILED ${bad}。自动重启 ${restarts} 次。目录 ${DEST}"
    log "DONE $msg"; : > "$COMPLETE"
    if [ "$ok" -eq "$EXPECT_MD5" ] && [ "$bad" -eq 0 ] && [ "$nf" -eq "$EXPECT_FILES" ]; then
      alert "✅[P47] 原始数据下载完成，md5 全部通过" "$msg"
    else
      alert "⚠[P47] 下载结束但核对未全过" "$msg"
    fi
    exit 0
  else
    log "STOPPED lftp gone, no ALL_DONE in $L; files=${nf}/${EXPECT_FILES}"
    restart_dl "lftp 已退出且日志 $(basename "$L") 无 ALL_DONE（本地文件 ${nf}/${EXPECT_FILES}）"
    sleep 600; continue   # 给新进程 10 min 建连后再巡检
  fi
  sleep $((POLL-60))
done
