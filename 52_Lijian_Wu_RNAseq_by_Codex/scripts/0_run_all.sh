#!/usr/bin/env bash
set -Eeuo pipefail
P=/home/gao/projects_2026Q4/52_Lijian_Wu_RNAseq_by_Codex
trap 'echo "MASTER_FAILED line=$LINENO $(date -Is)" >> "$P/logs/MASTER_STATUS"; echo "P52 master failed; inspect $P/logs/master.log" | /home/gao/Utilities/Server/monitoring_alerting/wd_sendmail.py "[P52] MASTER_FAILED" - || true' ERR
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$P/scripts"
echo "FASTQ_VERIFICATION_RUNNING $(date -Is)" >> "$P/logs/MASTER_STATUS"
python3 -u 1_verify_fastq.py
bash 2_nextflow.sh
echo "QUANTIFICATION_COMPLETE $(date -Is)" >> "$P/logs/MASTER_STATUS"
# Downstream is prepared by the active agent while quantification runs.
while [ ! -f "$P/logs/DOWNSTREAM_SCRIPTS_READY.marker" ]; do
 echo "$(date -Is) quantification complete; waiting for reviewed downstream scripts"
 sleep 60
done
bash 7_orchestrate_downstream.sh
echo "DOWNSTREAM_COMPLETE $(date -Is)" >> "$P/logs/MASTER_STATUS"
# Keep active lifecycle attached until report and delivery have been audited.
while [ ! -f "$P/logs/DELIVERY_COMPLETE.marker" ]; do sleep 60; done
