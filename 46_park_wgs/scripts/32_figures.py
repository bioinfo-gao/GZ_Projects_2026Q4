"""P46 客户报告图件（2026-11-26）。

规则来源：`~/.claude/skills/scientific-figure-making/SKILL.md`（V21 十色系三色调为唯一色值来源）
+ `references/house-style.md`（字号/线宽/导出）+ `references/api.md`（断言）。
**交付只放 PDF**；PNG 仅作我方逐张查看，落 tmp/，不进交付目录。

每个数字都从源文件现读，不从内部文档抄（verify-every-number 规则）。
"""
import matplotlib; matplotlib.use("Agg")
import sys, csv, json, os, collections
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/gao/.claude/skills/scientific-figure-making/assets")
sys.path.insert(0, "/home/gao/Utilities/color_and_format/python/color_figures_with_preconfigured_sets")
from pubfig import (apply_publication_style, FigureStyle, finalize_figure, legend_outside,
                    assert_no_palette_collision, assert_visible_on_white, assert_points_drawn)
from palette_priority import select_palette

P = "/home/gao/projects_2026Q4/46_park_wgs"
QC = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/qc"
COV = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/coverage"
FIG = f"{P}/analysis/stage5_figures_20261126"
TMP = f"{FIG}/tmp_png_preview"
os.makedirs(TMP, exist_ok=True)
G = 3215721559
apply_publication_style(FigureStyle(font_size=15, axes_linewidth=2))

st = list(csv.DictReader(open(f"{P}/scripts/0_sample_table.csv")))
NAMES = [r["safe_name"] for r in st]
ROLE = {r["safe_name"]: ("Control" if r["role"].startswith("control") else "Edited") for r in st}
VENDOR = {r["safe_name"]: int(r["vendor_read_pairs"]) for r in st}

pal = select_palette(4)
COL = pal["fill"] if isinstance(pal, dict) else list(pal)
C_CTRL, C_EDIT, C_A, C_B = COL[0], COL[1], COL[2], COL[3]
assert_visible_on_white([C_CTRL, C_EDIT, C_A, C_B], role="fill")
GRP = {"Control": C_CTRL, "Edited": C_EDIT}


def fig_legend(fig, axes, top=0.76):
    """图例统一放【全图顶部】、子图整体下移 —— 避免与 set_title 同处 axes 顶端造成压字。
    skill 要求图例在绘图区外且禁用 loc='best'；fig.legend 满足前者，且显式给 loc。"""
    seen, H, L = set(), [], []
    for ax in np.atleast_1d(axes).ravel():
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l and l not in seen:
                seen.add(l); H.append(h); L.append(l)
        if ax.get_legend() is not None: ax.get_legend().remove()
    if H:
        fig.legend(H, L, loc="upper center", ncol=min(len(H), 4), frameon=False,
                   bbox_to_anchor=(0.5, 1.0))
    fig.subplots_adjust(top=top)

def tsv(p, d="\t"): return list(csv.DictReader(open(p), delimiter=d))
def bar_groups(ax, vals, ylab, title, annot=None):
    x = np.arange(len(NAMES))
    for g, c in GRP.items():
        idx = [i for i, s in enumerate(NAMES) if ROLE[s] == g]
        ax.bar(x[idx], [vals[NAMES[i]] for i in idx], color=c, edgecolor="black", linewidth=0.8, label=g)
    ax.set_xticks(x); ax.set_xticklabels(NAMES, rotation=45, ha="right")
    ax.set_ylabel(ylab); ax.set_title(title, loc="left", pad=34)
    if annot:
        for i, s in enumerate(NAMES): ax.text(i, vals[s], annot(vals[s]), ha="center", va="bottom", fontsize=9)

# ── Fig 1：测序深度双口径 + insert size ───────────────────────────────
raw, eff, ins, dup = {}, {}, {}, {}
for s in NAMES:
    v = {}
    for ln in open(f"{QC}/{s}.samtools_stats.txt"):
        p = ln.rstrip("\n").split("\t")
        if p[0] == "SN": v[p[1]] = float(p[2])
    ins[s] = v["insert size average:"]
    raw[s] = VENDOR[s] * 2 * 150 / G
    for ln in open(f"{COV}/{s}.mosdepth.summary.txt"):
        f = ln.split("\t")
        if f[0] == "total": eff[s] = float(f[3])
    fl = open(f"{QC}/{s}.flagstat.txt").read()
    pr = int([l for l in fl.split("\n") if l.endswith(" primary")][0].split()[0])
    du = int([l for l in fl.split("\n") if "primary duplicates" in l][0].split()[0])
    dup[s] = 100 * du / pr

fig, axes = plt.subplots(1, 3, figsize=(19, 5.4))
x = np.arange(len(NAMES)); w = 0.38
axes[0].bar(x - w/2, [raw[s] for s in NAMES], w, color=C_A, edgecolor="black", linewidth=.8, label="Raw output")
axes[0].bar(x + w/2, [eff[s] for s in NAMES], w, color=C_B, edgecolor="black", linewidth=.8, label="Usable (dedup + mate-overlap corrected)")
axes[0].axhline(45, color="#E41A1C", linestyle="--", linewidth=1.6)
axes[0].set_xlim(-0.8, 10.6)
axes[0].text(9.7, 45, "45x\nordered", color="#E41A1C", ha="left", va="center", fontsize=10)
axes[0].set_xticks(x); axes[0].set_xticklabels(NAMES, rotation=45, ha="right")
axes[0].set_ylabel("Depth (x)"); axes[0].set_title("A  Sequencing depth, two conventions", loc="left", pad=34)
for g, c in GRP.items():
    idx = [i for i, s in enumerate(NAMES) if ROLE[s] == g]
    axes[1].scatter(idx, [ins[NAMES[i]] for i in idx], s=95, color=c, edgecolor="black", linewidth=.8, label=g, zorder=3)
assert_points_drawn(axes[1], len(NAMES))
axes[1].axhline(150, color="#E41A1C", linestyle="--", linewidth=1.6)
axes[1].text(9.4, 151.5, "read length 150 bp", color="#E41A1C", fontsize=10, ha="right")
axes[1].set_xticks(x); axes[1].set_xticklabels(NAMES, rotation=45, ha="right")
axes[1].set_ylim(135, 196); axes[1].set_ylabel("Mean insert size (bp)")
axes[1].set_title("B  Library insert size", loc="left", pad=34)
bar_groups(axes[2], dup, "Duplicate rate (%)", "C  PCR/optical duplicate rate")
assert_no_palette_collision({"depth": [C_A, C_B], "group": [C_CTRL, C_EDIT]})
fig_legend(fig, axes)
finalize_figure(fig, f"{FIG}/Fig1_sequencing_depth_and_library", formats=["pdf"], dpi=300)
finalize_figure(fig, f"{TMP}/Fig1", formats=["png"], dpi=150); plt.close(fig)

# ── Fig 2：靶点 Δ202 ─────────────────────────────────────────────────
on = tsv(f"{P}/analysis/stage1_ontarget_20261123/tables/ontarget_summary.tsv")
frac = {r["sample"]: (float(r["deleted_allele_frac_CUT1"]) + float(r["deleted_allele_frac_CUT2"])) / 2 for r in on}
jr = {r["sample"]: int(r["junction_reads"]) for r in on}
fig, axes = plt.subplots(1, 2, figsize=(14, 5.4))
bar_groups(axes[0], frac, "Mean deleted-allele fraction", "A  Fraction of alleles carrying the 202 bp deletion",
           annot=lambda v: f"{v:.2f}")
for y, lab in [(0.5, "heterozygous (0.5)"), (1.0, "both alleles (1.0)")]:
    axes[0].axhline(y, color="#8C510A", linestyle=":", linewidth=1.4)
    axes[0].text(-0.35, y + .03, lab, ha="left", fontsize=9, color="#8C510A")
axes[0].annotate("Sample_5: one allele deleted,\nthe other carries a small indel at cut 1",
                 xy=(7, 0.75), xytext=(2.6, 0.86), fontsize=9, color="#8C510A",
                 arrowprops=dict(arrowstyle="->", color="#8C510A"))
axes[0].set_ylim(0, 1.18); legend_outside(axes[0])
FLOOR = 0.45
for g, c in GRP.items():
    idx = [i for i, s in enumerate(NAMES) if ROLE[s] == g]
    nz = [i for i in idx if jr[NAMES[i]] > 0]; zr = [i for i in idx if jr[NAMES[i]] == 0]
    if nz: axes[1].scatter(nz, [jr[NAMES[i]] for i in nz], s=90, color=c, edgecolor="black",
                           linewidth=.8, label=g, zorder=3)
    if zr: axes[1].scatter(zr, [FLOOR] * len(zr), s=90, facecolors="white", edgecolor=c,
                           linewidth=1.8, zorder=3, label="Zero reads" if g == "Control" else None)
assert_points_drawn(axes[1], len(NAMES))
axes[1].axhline(FLOOR * 1.35, color="#999999", linewidth=.9, linestyle="-")
axes[1].text(3.4, FLOOR * 1.02, "open symbols = 0 reads (no junction detected)", fontsize=9, color="#5A5A5A")
axes[1].set_yscale("log"); axes[1].set_xticks(x); axes[1].set_xticklabels(NAMES, rotation=45, ha="right")
axes[1].set_ylabel("Junction-spanning reads (log scale)")
axes[1].set_title("B  Reads crossing the deletion junction", loc="left", pad=34)
axes[1].annotate("Sample_6: amplicon-like contamination\n(genotype not readable at face value)",
                 xy=(8, jr["Sample_6"]), xytext=(4.2, 1600), fontsize=9,
                 arrowprops=dict(arrowstyle="->", color="#E41A1C"), color="#E41A1C")
fig_legend(fig, axes)
finalize_figure(fig, f"{FIG}/Fig2_on_target_deletion", formats=["pdf"], dpi=300)
finalize_figure(fig, f"{TMP}/Fig2", formats=["png"], dpi=150); plt.close(fig)

# ── Fig 3：样本身份与污染 ────────────────────────────────────────────
M = tsv(f"{P}/analysis/stage0_identity_qc_20261123_v2/tables/identity_matrix.tsv")
mat = np.array([[float(r[c]) for c in NAMES] for r in M])
con = {r["sample"]: 100 * float(r["ref_frac_at_homalt"]) for r in tsv(f"{P}/analysis/stage0_identity_qc_20261123_v2/tables/contamination.tsv")}
fig, axes = plt.subplots(1, 2, figsize=(16.5, 5.8)); fig.subplots_adjust(wspace=0.42)
im = axes[0].imshow(mat, cmap="YlGnBu", vmin=0.995, vmax=1.0)
axes[0].set_xticks(x); axes[0].set_xticklabels(NAMES, rotation=45, ha="right")
axes[0].set_yticks(x); axes[0].set_yticklabels(NAMES)
axes[0].set_title("A  Pairwise genotype concordance", loc="left", pad=34)
axes[0].set_xlabel("26,674 common SNP sites", fontsize=10)
fig.colorbar(im, ax=axes[0], fraction=.046, label="Concordance")
bar_groups(axes[1], con, "REF reads at hom-ALT sites (%)", "B  Foreign-DNA contamination estimate",
           annot=lambda v: f"{v:.2f}")
fig_legend(fig, axes)
finalize_figure(fig, f"{FIG}/Fig3_sample_identity_and_contamination", formats=["pdf"], dpi=300)
finalize_figure(fig, f"{TMP}/Fig3", formats=["png"], dpi=150); plt.close(fig)

# ── Fig 4：chr12p 拷贝数（两法互证）────────────────────────────────
B = tsv(f"{P}/analysis/stage4c_cnv_20261123/baf_check/chr12p_baf_crosscheck.tsv")
arms = ["chr12p", "chr12q", "chr7(对照臂)"]
LAB = {"chr12p": "chr12p", "chr12q": "chr12q", "chr7(对照臂)": "chr7 (control arm)"}
fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))
for ai, arm in enumerate(arms):
    for g, c in GRP.items():
        v = [float(r["median_baf"]) for r in B if r["arm"] == arm and ROLE[r["sample"]] == g]
        axes[0].scatter([ai + (-.14 if g == "Control" else .14)] * len(v), v, s=85, color=c,
                        edgecolor="black", linewidth=.8, label=g if ai == 0 else None, zorder=3)
axes[0].axhline(0.5, color="#8C510A", linestyle=":", linewidth=1.4)
axes[0].text(2.45, .505, "0.50 = 2 copies", ha="right", fontsize=9, color="#8C510A")
axes[0].axhline(0.25, color="#8C510A", linestyle=":", linewidth=1.4)
axes[0].text(2.45, .255, "0.25 = 4 copies (3:1)", ha="right", fontsize=9, color="#8C510A")
axes[0].set_xticks(range(3)); axes[0].set_xticklabels([LAB[a] for a in arms])
axes[0].set_ylabel("Median B-allele fraction"); axes[0].set_ylim(.18, .55)
axes[0].set_title("A  Allele-ratio evidence", loc="left", pad=34); legend_outside(axes[0])
for ai, arm in enumerate(arms):
    for g, c in GRP.items():
        v = [100 * float(r["frac_below_0_40"]) for r in B if r["arm"] == arm and ROLE[r["sample"]] == g]
        axes[1].scatter([ai + (-.14 if g == "Control" else .14)] * len(v), v, s=85, color=c,
                        edgecolor="black", linewidth=.8, label=g if ai == 0 else None, zorder=3)
axes[1].set_xticks(range(3)); axes[1].set_xticklabels([LAB[a] for a in arms])
axes[1].set_ylabel("Heterozygous sites with BAF < 0.40 (%)")
axes[1].set_title("B  Same data, fraction of skewed sites", loc="left", pad=34); legend_outside(axes[1])
fig_legend(fig, axes)
finalize_figure(fig, f"{FIG}/Fig4_chr12p_copy_number", formats=["pdf"], dpi=300)
finalize_figure(fig, f"{TMP}/Fig4", formats=["png"], dpi=150); plt.close(fig)

# ── Fig 5：脱靶 ──────────────────────────────────────────────────────
OT = tsv(f"{P}/analysis/stage2b_offtarget_evidence_20261123_v2/tables/offtarget_evidence_n10.tsv")
off = [r for r in OT if r["is_ontarget"] != "1"]
cnt = collections.Counter(r["pattern"] for r in off)
order = ["none", "both", "treated_only", "control_only"]
LBL = {"none": "No evidence\nin any sample", "both": "Edited and\ncontrol clones",
       "treated_only": "Edited clones\nonly", "control_only": "Control clones\nonly"}
fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))
cols = [C_A, C_B, C_EDIT, C_CTRL]
bars = axes[0].bar(range(4), [cnt[k] for k in order], color=cols, edgecolor="black", linewidth=.8)
for b, k in zip(bars, order):
    axes[0].text(b.get_x() + b.get_width()/2, b.get_height(), str(cnt[k]), ha="center", va="bottom", fontsize=11)
axes[0].set_xticks(range(4)); axes[0].set_xticklabels([LBL[k] for k in order], fontsize=10)
axes[0].set_ylabel("Candidate off-target sites")
axes[0].set_title(f"A  Evidence pattern across {len(off):,} candidate sites", loc="left", pad=34)
pos = {s: sum(1 for r in off if r[f"{s}_positive"] == "1") for s in NAMES}
bar_groups(axes[1], pos, "Sites flagged", "B  Flagged sites per sample", annot=lambda v: str(v))
assert_no_palette_collision({"pattern": cols[:2], "group": [C_CTRL, C_EDIT]})
fig_legend(fig, axes)
finalize_figure(fig, f"{FIG}/Fig5_off_target_screen", formats=["pdf"], dpi=300)
finalize_figure(fig, f"{TMP}/Fig5", formats=["png"], dpi=150); plt.close(fig)

# ── Fig 6：载体拷贝数（归一后）──────────────────────────────────────
V = tsv(f"{P}/analysis/stage5_verify_20261126/tables/vector_copies_normalized.tsv")
fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.4))
for pi, (vec, ttl) in enumerate([("pLPC_mCherry_BP1_2", "A  53BP1-mCherry reporter vector"),
                                 ("pMSCV_Cas9_EGFP", "B  Cas9-EGFP vector (not used this round)")]):
    d = {r["sample"]: float(r["est_copies_normalized"]) for r in V
         if r["vector"] == vec and r["region"] == "WHOLE_VECTOR"}
    bar_groups(axes[pi], d, "Estimated copies per diploid genome", ttl, annot=lambda v: f"{v:.2f}")
axes[0].annotate("Sample_6: amplicon-like contamination", xy=(8, 7.09), xytext=(3.0, 5.6), fontsize=9,
                 arrowprops=dict(arrowstyle="->", color="#E41A1C"), color="#E41A1C")
fig_legend(fig, axes)
finalize_figure(fig, f"{FIG}/Fig6_vector_copy_number", formats=["pdf"], dpi=300)
finalize_figure(fig, f"{TMP}/Fig6", formats=["png"], dpi=150); plt.close(fig)

# ── Fig 7：私有变异负荷 ──────────────────────────────────────────────
PB = {r["sample"]: r for r in tsv(f"{P}/analysis/stage4d_private_burden_20261125/tables/private_burden_per_sample.tsv")}
PS = {r["sample"]: r for r in tsv(f"{P}/analysis/stage5_prep_20261125/tables/private_sv_per_sample.tsv")}
fig, axes = plt.subplots(1, 3, figsize=(19, 5.4))
for pi, (key, src, ylab, ttl) in enumerate([
        ("snv", PB, "Private SNVs", "A  Private single-nucleotide variants"),
        ("indel", PB, "Private indels", "B  Private insertions/deletions"),
        ("private_sv", PS, "Private structural variants", "C  Private structural variants")]):
    d = {s: int(src[s][key]) for s in NAMES}
    bar_groups(axes[pi], d, ylab, ttl, annot=lambda v: str(v))
fig_legend(fig, axes)
finalize_figure(fig, f"{FIG}/Fig7_private_variant_burden", formats=["pdf"], dpi=300)
finalize_figure(fig, f"{TMP}/Fig7", formats=["png"], dpi=150); plt.close(fig)

print("[done] PDF ->", FIG)
for f in sorted(os.listdir(FIG)):
    if f.endswith(".pdf"): print("   ", f, os.path.getsize(f"{FIG}/{f}"), "B")
print("[preview PNG, 仅我方查看，不进交付]", TMP)
