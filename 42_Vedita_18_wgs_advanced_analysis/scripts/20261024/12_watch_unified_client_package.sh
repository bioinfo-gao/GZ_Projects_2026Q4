#!/usr/bin/env bash
# Local/TTY watchdog for the long client-package transfer; no external mail is sent.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
LOG_DIR="${PROJECT}/logs/20261024_client_package"
WATCHDOG=/home/gao/Utilities/Server/monitoring_alerting/nextflow_watchdog.sh
logs=("${LOG_DIR}"/client_package_*.log)
[[ -e "${logs[0]}" ]] || { echo "No package log exists yet" >&2; exit 1; }
latest=$(ls -1t "${logs[@]}" | head -n 1)
exec env WD_PERSIST=1 WD_POLL=600 WD_HEARTBEAT=6 WD_LOW_MEM_MB=2500 \
    WD_STALE_MIN=0 WD_LOW_LOAD=0 WD_MAX_TASK_H=0 WD_MAIL=0 WD_TTY=1 \
    WD_LOG="${LOG_DIR}/watchdog.log" WD_ALERT_LOG="${LOG_DIR}/watchdog_ALERTS.log" \
    bash "${WATCHDOG}" "p42_client_package_20261024:${PROJECT}:${latest}"
