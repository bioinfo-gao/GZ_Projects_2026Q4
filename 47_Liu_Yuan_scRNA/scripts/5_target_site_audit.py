#!/usr/bin/env python3
"""5_target_site_audit.py —— 靶位点逐 record 证据审计 + 分子层面计数（2026-11-22 建）

对应计划 §3 decision table 与规则 §2A：
  - data universe = BAM 中覆盖该位点的全部 record（不预先用 caller / 旧表封闭）
  - primary / secondary / supplementary 分开记账；secondary/supplementary 不计入独立分子，但进审计表
  - 每道过滤都记录"过滤前→后"与被滤原因；被滤 record 写进 excluded 表（可追溯 FLAG/CIGAR/MAPQ）
  - 分子 = (CB, UB)；分子内多条 read 取多数碱基，内部冲突的分子单列（不悄悄丢）
用法:
  python -I 5_target_site_audit.py --bam X.bam --sites sites.tsv --out prefix [--cells barcodes.tsv] [--minmapq 20] [--minbq 20]
  sites.tsv: chrom<TAB>pos(1-based)<TAB>ref<TAB>alt<TAB>label   （ref/alt 为参考正链碱基）
输出: <prefix>.site_audit.tsv / <prefix>.molecules.tsv / <prefix>.excluded_records.tsv.gz
"""
import argparse, collections, gzip, sys
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True); ap.add_argument("--sites", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--cells"); ap.add_argument("--minmapq", type=int, default=20); ap.add_argument("--minbq", type=int, default=20)
a = ap.parse_args()

cells = None
if a.cells:
    op = gzip.open if a.cells.endswith(".gz") else open
    with op(a.cells, "rt") as f: cells = {l.strip().split("\t")[0] for l in f if l.strip()}

def cls(r):
    return "supplementary" if r.is_supplementary else ("secondary" if r.is_secondary else "primary")

bam = pysam.AlignmentFile(a.bam)
audit_rows, mol_rows = [], []
exc = gzip.open(a.out + ".excluded_records.tsv.gz", "wt")
exc.write("label\tchrom\tpos\tread\tflag\tclass\tmapq\tcigar\tCB\tUB\tbase\tbq\treason\n")

with open(a.sites) as f:
    sites = [l.rstrip("\n").split("\t") for l in f if l.strip() and not l.startswith("#")]
for chrom, pos, ref, alt, label in sites:
    p0 = int(pos) - 1
    c = collections.Counter()                     # 逐级计数（过滤影响审计）
    mols = collections.defaultdict(collections.Counter)
    for r in bam.fetch(chrom, p0, p0 + 1):
        k = cls(r); c[f"rec_{k}"] += 1; c["rec_all"] += 1
        # 该 read 在位点上的碱基（缺失 / 剪接跨越 / 软剪切 都要分开）
        base, bq = None, None
        for qp, rp in r.get_aligned_pairs(matches_only=False):
            if rp == p0:
                if qp is None: base = "DEL_or_SKIP"
                else: base = r.query_sequence[qp] if r.query_sequence else "NOSEQ"; bq = r.query_qualities[qp] if r.query_qualities is not None else None
                break
        cb = r.get_tag("CB") if r.has_tag("CB") else None
        ub = r.get_tag("UB") if r.has_tag("UB") else None
        reason = None
        if k != "primary": reason = f"not_primary:{k}"
        elif r.mapping_quality < a.minmapq: reason = f"mapq<{a.minmapq}"
        elif base in (None, "DEL_or_SKIP", "NOSEQ"): reason = f"no_base:{base}"
        elif bq is not None and bq < a.minbq: reason = f"bq<{a.minbq}"
        elif cb is None: reason = "no_CB"
        elif ub is None: reason = "no_UB"
        elif cells is not None and cb not in cells: reason = "CB_not_in_QC_cells"
        if reason:
            c[f"excl_{reason.split(':')[0] if reason.startswith(('not_primary','no_base')) else reason}"] += 1
            exc.write(f"{label}\t{chrom}\t{pos}\t{r.query_name}\t{r.flag}\t{k}\t{r.mapping_quality}\t{r.cigarstring}\t{cb}\t{ub}\t{base}\t{bq}\t{reason}\n")
            continue
        c["kept_reads"] += 1
        mols[(cb, ub)][base] += 1
    # 分子层面：多数碱基；并列或少数碱基占比 >0 记为 conflict
    bc = collections.Counter(); n_conf = 0
    for (cb, ub), cnt in mols.items():
        top = cnt.most_common(); b, nb = top[0]
        tie = len(top) > 1 and top[1][1] == nb
        conflict = len(top) > 1
        n_conf += conflict
        call = "TIE" if tie else b
        bc[call] += 1
        mol_rows.append((label, chrom, pos, cb, ub, call, sum(cnt.values()), int(conflict), ";".join(f"{x}:{y}" for x, y in top)))
    n_mol = len(mols)
    audit_rows.append({
        "label": label, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
        **{k: c.get(k, 0) for k in ["rec_all", "rec_primary", "rec_secondary", "rec_supplementary"]},
        **{k: v for k, v in c.items() if k.startswith("excl_")},
        "kept_reads": c.get("kept_reads", 0), "molecules": n_mol, "molecules_conflict": n_conf,
        "cells_with_molecule": len({cb for cb, _ in mols}),
        "mol_ref": bc.get(ref, 0), "mol_alt": bc.get(alt, 0),
        "mol_other": n_mol - bc.get(ref, 0) - bc.get(alt, 0),
        "alt_fraction_mol": (bc.get(alt, 0) / (bc.get(ref, 0) + bc.get(alt, 0))) if (bc.get(ref, 0) + bc.get(alt, 0)) else "NA",
    })
exc.close()

keys = []
for d in audit_rows:
    for k in d:
        if k not in keys: keys.append(k)
with open(a.out + ".site_audit.tsv", "w") as o:
    o.write("\t".join(keys) + "\n")
    for d in audit_rows: o.write("\t".join(str(d.get(k, 0)) for k in keys) + "\n")
with open(a.out + ".molecules.tsv", "w") as o:
    o.write("label\tchrom\tpos\tCB\tUB\tcall\tn_reads\tconflict\tbase_counts\n")
    for row in mol_rows: o.write("\t".join(map(str, row)) + "\n")
# 自洽断言：每个位点 primary+secondary+supplementary = 全部 record；kept + 各排除项 = 全部 record
for d in audit_rows:
    assert d["rec_primary"] + d["rec_secondary"] + d["rec_supplementary"] == d["rec_all"], d
    assert d["kept_reads"] + sum(v for k, v in d.items() if k.startswith("excl_")) == d["rec_all"], d
print(f"sites={len(audit_rows)}  molecules={len(mol_rows)}  -> {a.out}.site_audit.tsv")
