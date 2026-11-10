#!/usr/bin/env bash
set -uo pipefail
P=/home/gao/projects_2026Q4/52_Lijian_Wu_RNAseq_by_Codex
# Integrity verification and downstream have no Nextflow heartbeat; session and
# resource checks remain active throughout. Long-task inspection uses P/work.
exec env WD_MAILER=/home/gao/Utilities/Server/monitoring_alerting/wd_sendmail.py WD_PERSIST=1 WD_POLL=600 WD_LOG="$P/logs/watchdog.log" WD_ALERT_LOG="$P/logs/watchdog_ALERTS.log" WD_STALE_MIN=0 WD_LOW_LOAD=10 WD_LOW_LOAD_MIN=240 WD_MAX_TASK_H=5 bash /home/gao/Utilities/Server/monitoring_alerting/nextflow_watchdog.sh "p52_full:$P:$P/logs/master.log"
