#!/usr/bin/env python3
"""P46 — 把本项目登记进跨 Agent 的活跃项目 SSOT + 本 Agent 的续跑清单。

两个文件、两个职责（不得互相顶替）：
  1) Agent_Assets/runtime/active_projects.json       —— 跨 Agent 的唯一 SSOT（状态在这里）
  2) Agent_Assets/runtime/claude_code_active_project_autostart.json —— 本 Agent 的会话续跑登记

只做「追加本项目条目」或「更新本项目条目」，一律保留其他项目的所有既有字段。
用法：
  python3 -I scripts/7_register_active_project.py \
      --status running --current-stage "..." --next-stage "..." --tmux p46_stage0_align
"""
import argparse
import datetime
import json
import os
import sys

RUNTIME = "/home/gao/Utilities/AI_and_Agent/Agent_Assets/runtime"
SSOT = os.path.join(RUNTIME, "active_projects.json")
AUTOSTART = os.path.join(RUNTIME, "claude_code_active_project_autostart.json")

PID = "P46_Yeji_Park_10sample_CRISPR_dual_cut_wgs"
PROOT = "/home/gao/projects_2026Q4/46_park_wgs"
SESSION_ID = "28f7fbf1-a55a-48f3-a128-c5570f180ea3"

COMPLETION = (
    "Complete only when custom_research_report_YYYYMMDD/ exists for all 10 samples and covers: "
    "(1) PPP1R2 target-site deletion status per clone with base-pair breakpoints and allele "
    "state (3 controls WT, Sample 1-4 single-allele, Sample 5-7 bi-allelic, cross-checked "
    "against the client's gel/Sanger expectation as the external calibration point); "
    "(2) genome-wide off-target evidence at the enumerated gRNA1/gRNA2 candidate sites, "
    "7 treated vs 3 controls; (3) vector screen against both P31 maps (Cas9 was RNP this "
    "round, so vector sequence can only be carried over); (4) private-variant burden "
    "comparison treated vs controls with filter-impact audits; plus counts/QC, per-sample "
    "VCFs, PDF figures, report sections 1-9, and the pre-delivery self-audit. "
    "User review precedes any sending; the agent never sends."
)

CONTINUATION = (
    "Automated continuity resumption for active P46 (Yeji Park, 10-sample RPE-1 CRISPR "
    "dual-cut WGS at ~45-65x raw). Resume only P46 from this saved session. Read "
    "/home/gao/Utilities/AI_and_Agent/Agent_Assets/runtime/active_projects.json and "
    "/home/gao/projects_2026Q4/46_park_wgs/docs/P46_research_plan_and_review_1028.md "
    "(sections 8.0-8.5 are the stage plan), then inspect tmux p46_* sessions and "
    "/home/gao/projects_2026Q4/46_park_wgs/logs/, and continue the recorded next stage "
    "until the completion definition is met. Full-workflow automatic completion is "
    "authorized; external sending is NOT. Never delete or overwrite result directories; "
    "BAMs live under /srv/external_data/gao_projects/46_park_wgs_bam/. Update active-project "
    "state only to running, complete, or genuinely blocked with evidence."
)


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def upsert(path, key, finder, builder, updater):
    with open(path) as fh:
        doc = json.load(fh)
    items = doc[key]
    for it in items:
        if finder(it):
            updater(it)
            break
    else:
        items.append(builder())
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(doc, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    os.replace(tmp, path)
    return doc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--status", required=True, choices=["running", "blocked", "complete"])
    ap.add_argument("--current-stage", required=True)
    ap.add_argument("--next-stage", required=True)
    ap.add_argument("--tmux", default="")
    ap.add_argument("--blocked-reason", default="")
    ap.add_argument("--blocking-question", default="")
    a = ap.parse_args()

    if a.status == "blocked" and not a.blocked_reason:
        sys.exit("blocked 必须同时给 --blocked-reason")

    def ssot_entry():
        return {
            "id": PID,
            "project_root": PROOT,
            "agent": "claude_code",
            "status": a.status,
            "current_stage": a.current_stage,
            "next_stage": a.next_stage,
            "expected_tmux_session": a.tmux,
            "blocked_reason": a.blocked_reason,
            "blocking_question": a.blocking_question,
            "last_verified_at": now(),
            "completion_definition": COMPLETION,
        }

    def ssot_update(it):
        it.update({
            "agent": "claude_code",
            "status": a.status,
            "current_stage": a.current_stage,
            "next_stage": a.next_stage,
            "expected_tmux_session": a.tmux,
            "blocked_reason": a.blocked_reason,
            "blocking_question": a.blocking_question,
            "last_verified_at": now(),
            "completion_definition": COMPLETION,
        })

    upsert(SSOT, "projects", lambda it: it.get("id") == PID, ssot_entry, ssot_update)

    def auto_entry():
        return {
            "active_project_id": PID,
            "project_root": PROOT,
            "status": "armed" if a.status == "running" else "disarmed",
            "usage_limit_retry_offsets_seconds": [10800, 14400, 18000],
            "reset_grace_seconds": 300,
            "fallback_stale_seconds": 10800,
            "min_transcript_idle_seconds": 1200,
            "claude_binary": "/home/gao/.local/bin/claude",
            "claude_session_id": SESSION_ID,
            "session_cwd": "/home/gao",
            "session_transcript": "/home/gao/.claude/projects/-home-gao/%s.jsonl" % SESSION_ID,
            "permission_mode": "auto",
            "continuation_prompt": CONTINUATION,
            "last_agent_heartbeat": "",
            "last_agent_heartbeat_source": "",
            "last_usage_limit_at": "",
            "last_launch_at": "",
            "last_worker_session": "",
            "last_worker_log": "",
            "last_worker_result": "",
            "last_trigger": "",
            "missing_tmux_observed_at": "",
            "usage_episode_started_at": "",
            "usage_reset_at": "",
            "usage_last_attempt_at": "",
        }

    def auto_update(it):
        it["status"] = "armed" if a.status == "running" else "disarmed"
        it["claude_session_id"] = SESSION_ID
        it["session_transcript"] = "/home/gao/.claude/projects/-home-gao/%s.jsonl" % SESSION_ID
        it["continuation_prompt"] = CONTINUATION

    upsert(AUTOSTART, "entries",
           lambda it: it.get("active_project_id") == PID, auto_entry, auto_update)

    # 回读验证：两个文件都必须能解析，且本项目条目确实在里面、状态与本次一致
    with open(SSOT) as fh:
        ss = json.load(fh)
    me = [p for p in ss["projects"] if p["id"] == PID]
    assert len(me) == 1, "SSOT 里 P46 条目数为 %d" % len(me)
    assert me[0]["status"] == a.status and me[0]["agent"] == "claude_code"
    with open(AUTOSTART) as fh:
        au = json.load(fh)
    mine = [e for e in au["entries"] if e["active_project_id"] == PID]
    assert len(mine) == 1, "autostart 里 P46 条目数为 %d" % len(mine)

    print("OK  SSOT 共 %d 个项目（P46 status=%s, tmux=%s）；autostart 共 %d 条（P46 %s）"
          % (len(ss["projects"]), me[0]["status"], me[0]["expected_tmux_session"] or "-",
             len(au["entries"]), mine[0]["status"]))


if __name__ == "__main__":
    sys.exit(main())
