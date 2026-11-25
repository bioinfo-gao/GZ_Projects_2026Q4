#!/usr/bin/env bash
# P46 —— DELLY 之后的自动续跑链（2026-11-25 建；教训见 memory feedback_never_defer_a_stage_to_my_next_check）
#
# 为什么存在：11-24 夜间 4a/Stage3 跑完后 DELLY 因"等磁盘空闲"被我留到次日手动启动 ⇒ 空转 9.7 h。
#   规则改为：任何"等 X 完成后再跑 Y"必须当场以 tmux + until 等待链启动，不得留给下一次自检。
# 本脚本等 DELLY 的 ALL_DONE，随后自动跑完三项可计算的复核，并把结果写成表供交付阶段引用。
# 异常分支：DELLY 的 tmux 消失而 ALL_DONE 仍缺 ⇒ FATAL 退出并留日志，不静默空转。
set -uo pipefail
source /Work_bio/gao/miniforge3/etc/profile.d/conda.sh; conda activate regular_bioinfo
PROJ=/home/gao/projects_2026Q4/46_park_wgs
DELLY_LOG=$PROJ/logs/stage4b_delly_pipelined_20261125.log
DELLY_OUT=$PROJ/analysis/stage4b_delly_n10_20261125
OUT=$PROJ/analysis/stage5_prep_$(date +%Y%m%d)
say(){ echo "[$(date '+%F %T')] $*"; }
mkdir -p "$OUT/tables"

say "CHAIN_START 等 DELLY ALL_DONE"
GRACE=0
until grep -q 'STAGE4B_DELLY_EXIT=0' "$DELLY_LOG" 2>/dev/null; do
  if ! tmux has-session -t p46_stage4b_delly 2>/dev/null; then
    GRACE=$((GRACE+1))
    if [[ $GRACE -ge 3 ]]; then
      say "FATAL DELLY session 已消失而 STAGE4B_DELLY_EXIT=0 仍不存在 ⇒ 异常退出，不启动后续"
      say "  查 $DELLY_LOG 与 $DELLY_OUT/logs/"
      exit 1
    fi
  else GRACE=0; fi
  sleep 120
done
say "DELLY 已完成，开始交付前复核"

# ---- 复核 5：DELLY 私有 SV 与脱靶候选 / 整合位点清单比对 ----
say "CHECK5_START 私有 SV × 脱靶候选 / 整合位点"
python3 -I "$PROJ/scripts/29_crosscheck_sv_vs_lists.py" "$DELLY_OUT/delly_joint_n10.bcf" "$OUT/tables" 2>&1
say "CHECK5_EXIT=$?"

say "CHAIN_DONE 可计算复核已全部完成；剩下的是需要人眼的部分（图件逐张看、报告撰写）"
