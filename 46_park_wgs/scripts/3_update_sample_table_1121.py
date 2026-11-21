#!/usr/bin/env python3
"""P46 — 把测序交付的事实回填进样本表 SSOT（scripts/0_sample_table.csv）。

本轮新增的事实来源（两个，都不是推断）：
  1) Kuo-Chan Hung 2026-11-20 17:01 UTC 邮件 "Yeji Data Delivery and Metadata"
     正文表格：Sample Name ↔ Sample ID ↔ i7/i5 index ↔ Required reads ↔ Real Sequencing Reads。
     ⇒ 给出 sample_id、i7/i5、vendor_read_pairs。
  2) 本机 FASTQ 实测：文件名解码（flowcell/lane/库编号/样本ID）+ 文件大小。
     ⇒ 给出 fastq_r1 / fastq_r2 / R1_GiB / R2_GiB / lane。

⚠ 口径说明（为什么不在这里填 raw_read_pairs / raw_depth_x）：
  - `vendor_read_pairs` 是 **厂方报告值**，不是我方实测。两者必须分列，不能互相顶替。
  - `raw_read_pairs` 与 `raw_depth_x` 留空，待 Stage 0 的 flagstat / mosdepth 实测后回填，
    由 5_fill_sample_table_from_stage0.py 写入。退出码 0 不等于数字已核对。

⚠ 仅追加列、不改任何既有列的值（客户给定事实与分组裁定一律原样保留）。
"""
import csv
import os
import re
import sys

PROJ = "/home/gao/projects_2026Q4/46_park_wgs"
CSV = os.path.join(PROJ, "scripts", "0_sample_table.csv")
FASTQ_DIR = "/srv/external_data/gao_projects/46_park_wgs_raw/WALKUP-21683/fastq"

# --- 事实来源 1：Kuo-Chan 2026-11-20 邮件正文表格（逐行抄录，未改一个字符） ---
VENDOR = {
    # sample_name: (sample_id, i7, i5, required_reads_M, library_size_bp, real_sequencing_reads)
    "Ctr_1":    ("WGS000029", "ATTGAGGA", "CGAGATCT", 500, 300, 576728998),
    "Ctr_2":    ("WGS000030", "CAACCACA", "CGAGATCT", 500, 300, 544862727),
    "Ctr_3":    ("WGS000031", "GACTAGTA", "CGAGATCT", 500, 300, 497440989),
    "Sample 1": ("WGS000032", "CAATGGAA", "CGAGATCT", 500, 300, 567126162),
    "Sample 2": ("WGS000033", "CACTTCGA", "CGAGATCT", 500, 300, 560871311),
    "Sample 3": ("WGS000034", "CAGCGTTA", "CGAGATCT", 500, 300, 596138972),
    "Sample 4": ("WGS000035", "CATACCAA", "CGAGATCT", 500, 300, 535487760),
    "Sample 5": ("WGS000036", "CCAGTTCA", "CGAGATCT", 500, 300, 499253985),
    "Sample 6": ("WGS000037", "CCGAAGTA", "CGAGATCT", 500, 300, 675460326),
    "Sample 7": ("WGS000038", "CCGTGAGA", "CGAGATCT", 500, 300, 569098815),
}

NEW_COLS = ["safe_name", "sample_id", "i7_index", "i5_index", "vendor_required_reads_M",
            "vendor_library_size_bp", "vendor_read_pairs",
            "flowcell", "lane", "fastq_r1", "fastq_r2"]

# Illumina 文件名：<flowcell>_<lane>_<库编号>_<样本ID>_<_?>S<n>_L00<n>_R<n>_001.fastq.gz
PAT = re.compile(
    r"^(?P<flowcell>[A-Z0-9]+)_(?P<lanefield>\d+)_(?P<libno>\d+)_"
    r"(?P<sid>WGS\d+)_+S(?P<snum>\d+)_L(?P<lane>\d{3})_R(?P<read>[12])_001\.fastq\.gz$")


def index_fastqs():
    """扫描 FASTQ 目录，按 sample_id 归档 R1/R2。多 lane 时按 lane 分别保留。"""
    found = {}
    for fn in sorted(os.listdir(FASTQ_DIR)):
        m = PAT.match(fn)
        if not m:
            continue
        sid = m.group("sid")
        rec = found.setdefault(sid, {"flowcell": m.group("flowcell"), "lanes": set()})
        rec["lanes"].add(int(m.group("lane")))
        rec["R%s" % m.group("read")] = os.path.join(FASTQ_DIR, fn)
    return found


def gib(path):
    return round(os.path.getsize(path) / 1024 ** 3, 2)


def main():
    fq = index_fastqs()

    with open(CSV, newline="") as fh:
        rows = list(csv.DictReader(fh))

    fieldnames = list(rows[0].keys())
    # 新列插在 notes 之前，notes 永远留在最后一列（人读的自由文本）
    for c in NEW_COLS:
        if c not in fieldnames:
            fieldnames.insert(fieldnames.index("notes"), c)

    n_filled = 0
    for r in rows:
        name = r["sample_name"]
        assert name in VENDOR, "样本表里有厂方元数据未覆盖的样本：%r" % name
        sid, i7, i5, req_m, libsz, real_pairs = VENDOR[name]
        assert sid in fq, "厂方 sample_id %s 在 FASTQ 目录里找不到对应文件" % sid
        rec = fq[sid]
        assert "R1" in rec and "R2" in rec, "%s 缺 R1 或 R2" % sid
        assert len(rec["lanes"]) == 1, \
            "%s 跨 %d 个 lane，必须按 lane 汇总而不是按文件数：%s" % (sid, len(rec["lanes"]), sorted(rec["lanes"]))

        # safe_name：客户样本名含空格（"Sample 1"），不能直接做文件名/RG 的 SM 值。
        # 规则固定为「空格→下划线」，可逆、与客户命名一一对应，交付时再换回客户原名。
        r["safe_name"] = name.replace(" ", "_")
        r["sample_id"] = sid
        r["i7_index"] = i7
        r["i5_index"] = i5
        r["vendor_required_reads_M"] = req_m
        r["vendor_library_size_bp"] = libsz
        r["vendor_read_pairs"] = real_pairs
        r["flowcell"] = rec["flowcell"]
        r["lane"] = sorted(rec["lanes"])[0]
        r["fastq_r1"] = rec["R1"]
        r["fastq_r2"] = rec["R2"]
        r["R1_GiB"] = gib(rec["R1"])
        r["R2_GiB"] = gib(rec["R2"])
        n_filled += 1

    # --- 硬断言：10 个样本、20 个 FASTQ、sample_id 与文件一一对应、无重复 ---
    assert n_filled == 10, "填了 %d 个样本，应为 10" % n_filled
    sids = [r["sample_id"] for r in rows]
    assert len(set(sids)) == 10, "sample_id 有重复：%s" % sids
    paths = [r["fastq_r1"] for r in rows] + [r["fastq_r2"] for r in rows]
    assert len(set(paths)) == 20, "FASTQ 路径有重复，某个文件被两个样本认领"
    safe = [r["safe_name"] for r in rows]
    assert len(set(safe)) == 10, "safe_name 有重复：%s" % safe
    assert all(" " not in s for s in safe), "safe_name 仍含空格：%s" % safe
    # 分组列未被动过
    roles = [r["role"] for r in rows]
    assert sum(1 for x in roles if x.startswith("control")) == 3, "对照数不是 3"
    assert sum(1 for x in roles if x == "test") == 7, "处理样本数不是 7"

    with open(CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    total = sum(float(r["R1_GiB"]) + float(r["R2_GiB"]) for r in rows)
    print("OK  10 samples, 20 FASTQ, %.1f GiB total" % total)
    for r in rows:
        print("  %-9s %s  lane %s  %5.1f + %5.1f GiB  vendor %d pairs  [%s]"
              % (r["sample_name"], r["sample_id"], r["lane"],
                 float(r["R1_GiB"]), float(r["R2_GiB"]),
                 int(r["vendor_read_pairs"]), r["role"]))


if __name__ == "__main__":
    sys.exit(main())
