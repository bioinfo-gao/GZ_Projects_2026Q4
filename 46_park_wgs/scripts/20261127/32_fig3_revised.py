"""P46 Fig 3 修订版（2026-11-27），只重画 Fig 3；其余图沿用 scripts/32_figures.py 的 1126 产物。

修改（用户 2026-11-26 同意）：
  A：图内注明色标仅覆盖 0.995–1.000、全部成对一致率 ≥ min 值 ⇒ 同一个体；格子深浅差 < 0.1% 是噪声。
  B：标题改为 "Cross-individual DNA contamination (none detected)"；图内注明该方法测不到
     位点特异的扩增产物污染（Sample_6，见 §7.1）。
数据输入与 1126 版完全相同（只读），数字从源表现算并断言与报告正文一致。
"""
import matplotlib; matplotlib.use("Agg")
import sys, csv, os, textwrap
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/gao/.claude/skills/scientific-figure-making/assets")
sys.path.insert(0, "/home/gao/Utilities/color_and_format/python/color_figures_with_preconfigured_sets")
from pubfig import (apply_publication_style, FigureStyle, finalize_figure, assert_visible_on_white)
from palette_priority import select_palette

P = "/home/gao/projects_2026Q4/46_park_wgs"
FIG = f"{P}/analysis/stage5_figures_20261127"
TMP = f"{FIG}/tmp_png_preview"
os.makedirs(TMP, exist_ok=True)
out_pdf = f"{FIG}/Fig3_sample_identity_and_contamination.pdf"
if os.path.exists(out_pdf):
    sys.exit(f"FATAL 输出已存在，不覆盖：{out_pdf}")
apply_publication_style(FigureStyle(font_size=15, axes_linewidth=2))

st = list(csv.DictReader(open(f"{P}/scripts/0_sample_table.csv")))
NAMES = [r["safe_name"] for r in st]
ROLE = {r["safe_name"]: ("Control" if r["role"].startswith("control") else "Edited") for r in st}
pal = select_palette(4)
COL = pal["fill"] if isinstance(pal, dict) else list(pal)
C_CTRL, C_EDIT = COL[0], COL[1]          # 与 1126 版其余图同色
assert_visible_on_white([C_CTRL, C_EDIT], role="fill")
GRP = {"Control": C_CTRL, "Edited": C_EDIT}
def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))

M = tsv(f"{P}/analysis/stage0_identity_qc_20261123_v2/tables/identity_matrix.tsv")
mat = np.array([[float(r[c]) for c in NAMES] for r in M])
con = {r["sample"]: 100 * float(r["ref_frac_at_homalt"])
       for r in tsv(f"{P}/analysis/stage0_identity_qc_20261123_v2/tables/contamination.tsv")}
off = mat[~np.eye(len(NAMES), dtype=bool)]
lo, hi = off.min(), off.max()
# 与报告正文 "99.73–99.81%" 与 "0.12–0.26%" 对齐（同一事实两处出现 ⇒ 必须比对）
assert f"{100*lo:.2f}" == "99.73" and f"{100*hi:.2f}" == "99.81", (lo, hi)
assert f"{min(con.values()):.2f}" == "0.12" and f"{max(con.values()):.2f}" == "0.26", con
assert lo >= 0.995, "色标下限会截掉数据"

fig, axes = plt.subplots(1, 2, figsize=(16.5, 8.2))
fig.subplots_adjust(wspace=0.42, top=0.84, bottom=0.36)
x = np.arange(len(NAMES))
im = axes[0].imshow(mat, cmap="YlGnBu", vmin=0.995, vmax=1.0)
axes[0].set_xticks(x); axes[0].set_xticklabels(NAMES, rotation=45, ha="right")
axes[0].set_yticks(x); axes[0].set_yticklabels(NAMES)
axes[0].set_title("A  Pairwise genotype concordance", loc="left", pad=34)
fig.colorbar(im, ax=axes[0], fraction=.046, label="Concordance")

for g, c in GRP.items():
    idx = [i for i, s in enumerate(NAMES) if ROLE[s] == g]
    axes[1].bar(x[idx], [con[NAMES[i]] for i in idx], color=c, edgecolor="black", linewidth=0.8, label=g)
axes[1].set_xticks(x); axes[1].set_xticklabels(NAMES, rotation=45, ha="right")
axes[1].set_ylabel("REF reads at hom-ALT sites (%)")
axes[1].set_title("B  Cross-individual DNA contamination (none detected)", loc="left", pad=34)
for i, s in enumerate(NAMES):
    axes[1].text(i, con[s], f"{con[s]:.2f}", ha="center", va="bottom", fontsize=9)

H, L = axes[1].get_legend_handles_labels()
assert L == list(GRP), L                 # 图例与组别对齐
fig.legend(H, L, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))

noteA = (f"26,674 common SNP sites. Colour scale spans 0.995–1.000 only; every pair agrees at "
         f"{100*lo:.2f}–{100*hi:.2f}%,\nso all ten samples come from the same individual. "
         "Shade differences (<0.1%) are genotyping noise.")
noteB = (f"All samples {min(con.values()):.2f}–{max(con.values()):.2f}%: sequencing-error background, "
         "no DNA from another person.\nThis test cannot detect locus-specific PCR-product (amplicon) "
         "contamination,\nsuch as that found at the target locus in Sample_6 (Section 7.1).")
# finalize_figure 内部先 tight_layout 再存 ⇒ 会把注释位置作废（draft2 实测）。
# 故：先自己 tight_layout，再统一下移留底边，最后直接 savefig，重叠断言在最终布局上做。
fig.tight_layout(pad=2.0)
fig.subplots_adjust(bottom=0.36, top=0.84)
fig.canvas.draw(); R = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
notes = []
for ax, t in ((axes[0], noteA), (axes[1], noteB)):
    t = "\n".join(textwrap.wrap(t.replace("\n", " "), 62))
    tb = ax.get_tightbbox(R).transformed(inv)              # 含刻度标签的实际下沿
    notes.append((fig.text(tb.x0, tb.y0 - 0.025, t, ha="left", va="top", fontsize=10.5), ax))
fig.canvas.draw()
for txt, ax in notes:
    nb = txt.get_window_extent(R).transformed(inv); tb = ax.get_tightbbox(R).transformed(inv)
    assert nb.y1 < tb.y0 and nb.y0 > 0, ("注释与刻度标签重叠或出界", nb, tb)
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(f"{TMP}/Fig3.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("FIG3_1127_DONE", out_pdf)
