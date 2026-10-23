#!/usr/bin/env bash
# Persistent watchdog for the corrected all-five-sample P42 SvABA v5 stage.
set -uo pipefail
PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
WD=/home/gao/Utilities/Server/monitoring_alerting/nextflow_watchdog.sh
exec env WD_PERSIST=1 WD_POLL=600 \
  WD_LOG="$PROJECT/logs/20261023_fullcohort_repeat_resolve/watchdog_external_svaba_v5.log" \
  WD_ALERT_LOG="$PROJECT/logs/20261023_fullcohort_repeat_resolve/watchdog_external_svaba_v5_ALERTS.log" \
  WD_STALE_MIN=0 WD_LOW_LOAD=10 WD_LOW_LOAD_MIN=240 WD_MAX_TASK_H=12 \
  bash "$WD" "p42_external_svaba_v5_20261023:$PROJECT:$PROJECT/logs/20261023_fullcohort_repeat_resolve/06_external_svaba_dld1_v5_multisample.run.log"
