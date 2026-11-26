#!/usr/bin/env bash
# P47 round1 编排链：step13b 跑完 -> 14 -> 15 -> 16 -> 17 -> 18 -> 19
# 为什么写成链而不是我下次检查再启动：定时自检是汇报机制不是执行保证（P46 DELLY 空转 9.7h 的教训）。
# 为什么用一个常驻 session：阶段性 tmux 跑完即消失，会让 continuity_guard 每 10 分钟误报。
set -uo pipefail
P=/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA
L=/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/logs
C=/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125
PY=/mnt/ex_8T_SSD/conda_envs/sc_tut_main/bin/python
RS=/mnt/ex_8T_SSD/conda_envs/p47_r/bin/Rscript
mkdir -p "$L"
ts(){ date '+%Y-%m-%d %H:%M:%S'; }
note(){ echo "[$(ts)] $*"; }

note "CHAIN_START  等待 step 13b 的 ckpt_13_qc.h5ad"
# 等上游：既接受 ckpt 文件出现，也接受 13b 日志里的 PASS；若 13b 进程已死且无产物则失败退出
until [ -s "$C/ckpt_13_qc.h5ad" ] && grep -q STEP13B_PASS "$L/13b_qc_1125.log" 2>/dev/null; do
  if ! pgrep -f "13b_qc_mad_filter.py" >/dev/null 2>&1; then
    if [ -s "$C/ckpt_13_qc.h5ad" ] && grep -q STEP13B_PASS "$L/13b_qc_1125.log" 2>/dev/null; then break; fi
    note "CHAIN_FAIL  13b 进程已退出但没有产出 ckpt_13 / STEP13B_PASS"; exit 1
  fi
  sleep 30
done
note "step 13b 就绪：$(ls -la "$C/ckpt_13_qc.h5ad" | awk '{print $5" bytes"}')"

run(){  # run <名字> <日志名> <命令...>
  local name="$1" log="$2"; shift 2
  note "RUN_START  $name"
  if "$@" > "$L/$log" 2>&1; then
    note "RUN_OK     $name  ($(grep -cE 'PASS$' "$L/$log" 2>/dev/null || echo 0) PASS 标记)"
  else
    note "RUN_FAIL   $name  -> 见 $L/$log"
    tail -25 "$L/$log" | sed 's/^/    | /'
    return 1
  fi
}

run "14 归一化/HVG/整合/聚类" 14_cluster_1125.log  $PY -u "$P/scripts/14_norm_hvg_integrate_cluster.py" || exit 1
run "15 ambient decontX"      15_ambient_1125.log  $RS    "$P/scripts/15_ambient_decontx.R"            || note "15 失败，继续（ambient 非后续硬依赖）"
run "16 pseudobulk 导出"      16_pseudobulk_1125.log $PY -u "$P/scripts/16_pseudobulk_export.py"        || exit 1
run "17 描述性效应量"         17_contrasts_1125.log  $PY -u "$P/scripts/17_descriptive_contrasts.py"    || exit 1
run "18 自审"                 18_audit_1125.log      $PY -u "$P/scripts/18_self_audit.py"               || note "18 有 FAIL 条目，见日志"
run "19 图件"                 19_figures_1125.log    $PY -u "$P/scripts/19_figures.py"                  || note "19 失败，见日志"

note "CHAIN_DONE"
echo "=== 产物清单 ==="
ls -la "$P/output_results/round1_descriptive_1125/" | tail -40
echo "=== 图件 ==="
ls -la "$P/output_results/round1_descriptive_1125/figures_1125/" 2>/dev/null
echo "=== 自审 ==="
cat "$P/output_results/round1_descriptive_1125/18_self_audit.csv" 2>/dev/null
