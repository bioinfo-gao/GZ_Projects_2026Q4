#!/usr/bin/env bash
# Persistent watchdog for the long-lived P42 coverage / assembly tmux session.
set -uo pipefail
PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
WD=/home/gao/Utilities/Server/monitoring_alerting/nextflow_watchdog.sh
# This is a bash/mosdepth batch, not Nextflow: .nextflow.log has no meaning.
# Keep session, memory and load checks; effectively disable only STALLED.
exec env WD_PERSIST=1 WD_POLL=600 \
  WD_LOG="$PROJECT/logs/20261023_fullcohort_repeat_resolve/watchdog.log" \
  WD_ALERT_LOG="$PROJECT/logs/20261023_fullcohort_repeat_resolve/watchdog_ALERTS.log" \
  WD_STALE_MIN=0 WD_LOW_LOAD=10 WD_LOW_LOAD_MIN=240 WD_MAX_TASK_H=5 \
  bash "$WD" "p42_coverage_20261023:$PROJECT:$PROJECT/logs/20261023_fullcohort_repeat_resolve/03_coverage_driver.log"
