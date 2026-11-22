#!/usr/bin/env python3
"""3_detect_chemistry.py —— 从 R1 判 10x 化学（2026-11-22 建）

用法: python3 -I 3_detect_chemistry.py <R1.fastq.gz> [n_reads=200000]
判据: R1 前 16 bp 对 v2/v3/v4 白名单的精确匹配率 + R1 长度（v2: 26 bp = 16 CB + 10 UMI；v3/v4: 28 bp = 16 + 12）。
输出: 每个白名单的匹配率、R1 长度分布、判定。最高匹配率 < 70% 或与第二名差 < 20 个百分点 ⇒ 判 AMBIGUOUS，须人工确认。
v3 与 v4 白名单可能大量重叠 ⇒ 同时看重叠率；Flex/5'/多组学不在本脚本判定范围内，匹配率全低时报 UNKNOWN。
"""
import gzip, sys, collections, itertools
WL = "/Work_bio/references/10x_whitelists_ngs_tools_1.8.6"
LISTS = {"v2": f"{WL}/10x_version2_whitelist.txt.gz", "v3": f"{WL}/10x_version3_whitelist.txt.gz",
         "v4": f"{WL}/10x_version4_whitelist.txt.gz"}
r1 = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 200_000

bcs, lens = [], collections.Counter()
with gzip.open(r1, "rt") as f:
    for i, ln in enumerate(itertools.islice(f, 4 * n)):
        if i % 4 == 1:
            s = ln.strip(); lens[len(s)] += 1; bcs.append(s[:16])
uniq = set(bcs)
rate = {}
for k, p in LISTS.items():
    with gzip.open(p, "rt") as f:
        wl = {x.split()[0] for x in f}
    rate[k] = sum(b in wl for b in bcs) / len(bcs)
    print(f"{k}: whitelist={len(wl):,}  exact-match reads={rate[k]:.1%}")
print(f"reads sampled={len(bcs):,}  unique CB16={len(uniq):,}")
print("R1 length:", dict(sorted(lens.items(), key=lambda x: -x[1])[:4]))
best, second = sorted(rate, key=rate.get, reverse=True)[:2]
L = lens.most_common(1)[0][0]
if rate[best] < 0.70: call = "UNKNOWN (no 3' GEX whitelist >70%; check Flex/5'/multiome)"
elif rate[best] - rate[second] < 0.20: call = f"AMBIGUOUS {best} vs {second} — confirm with client/library prep"
else: call = best
if call in ("v3", "v4") and L < 28: call += f" BUT R1 length {L} < 28 — UMI may be truncated, check"
if call == "v2" and L < 26: call += f" BUT R1 length {L} < 26"
print("CALL:", call)
