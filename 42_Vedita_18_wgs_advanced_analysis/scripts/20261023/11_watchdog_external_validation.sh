#!/usr/bin/env bash
# Persistent watchdog for P42 targeted external SV validation.
set -uo pipefail
PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
WD=/home/gao/Utilities/Server/monitoring_alerting/nextflow_watchdog.sh
exec env WD_PERSIST=1 WD_POLL=600 \
  WD_LOG="$PROJECT/logs/20261023_fullcohort_repeat_resolve/watchdog_external_validation.log" \
  WD_ALERT_LOG="$PROJECT/logs/20261023_fullcohort_repeat_resolve/watchdog_external_validation_ALERTS.log" \
  WD_STALE_MIN=0 WD_LOW_LOAD=10 WD_LOW_LOAD_MIN=240 WD_MAX_TASK_H=12 \
  bash "$WD" "p42_external_sv_validation_20261023:$PROJECT:$PROJECT/logs/20261023_fullcohort_repeat_resolve/06_external_sv_validation.run.log"
