"""P46 Fig 5 修订版（2026-11-27），只重画 Fig 5；数据输入与 1126 版完全相同（只读）。

修改（用户 2026-11-27 指出）：1126 版用 fig_legend() 把 Control/Edited 图例放在【全图顶部】，
看起来像两个 panel 共用，实际只属于 Panel B（A 的四类已由 x 轴标签说明，不需要图例）。
⇒ 图例改挂在 Panel B 自己的绘图区右侧外侧；断言 A 无图例、B 图例项与组别一致。
"""
import matplotlib; matplotlib.use("Agg")
import sys, csv, os, collections
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/gao/.claude/skills/scientific-figure-making/assets")
sys.path.insert(0, "/home/gao/Utilities/color_and_format/python/color_figures_with_preconfigured_sets")
from pubfig import (apply_publication_style, FigureStyle, assert_no_palette_collision,
                    assert_visible_on_white, legend_outside)
from palette_priority import select_palette

P = "/home/gao/projects_2026Q4/46_park_wgs"
FIG = f"{P}/analysis/stage5_figures_20261127"
TMP = f"{FIG}/tmp_png_preview"
os.makedirs(TMP, exist_ok=True)
out_pdf = f"{FIG}/Fig5_off_target_screen.pdf"
if os.path.exists(out_pdf):
    sys.exit(f"FATAL 输出已存在，不覆盖：{out_pdf}")
apply_publication_style(FigureStyle(font_size=15, axes_linewidth=2))

st = list(csv.DictReader(open(f"{P}/scripts/0_sample_table.csv")))
NAMES = [r["safe_name"] for r in st]
ROLE = {r["safe_name"]: ("Control" if r["role"].startswith("control") else "Edited") for r in st}
pal = select_palette(4)
COL = pal["fill"] if isinstance(pal, dict) else list(pal)
C_CTRL, C_EDIT, C_A, C_B = COL[0], COL[1], COL[2], COL[3]   # 与 1126 版同色
assert_visible_on_white([C_CTRL, C_EDIT, C_A, C_B], role="fill")
GRP = {"Control": C_CTRL, "Edited": C_EDIT}
def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))

OT = tsv(f"{P}/analysis/stage2b_offtarget_evidence_20261123_v2/tables/offtarget_evidence_n10.tsv")
off = [r for r in OT if r["is_ontarget"] != "1"]
cnt = collections.Counter(r["pattern"] for r in off)
order = ["none", "both", "treated_only", "control_only"]
# 与报告 §7.2 表一致（同一事实两处出现 ⇒ 比对）
assert len(off) == 1584 and [cnt[k] for k in order] == [1556, 20, 7, 1], (len(off), cnt)
LBL = {"none": "No evidence\nin any sample", "both": "Edited and\ncontrol clones",
       "treated_only": "Edited clones\nonly", "control_only": "Control clones\nonly"}

fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.4))
cols = [C_A, C_B, C_EDIT, C_CTRL]
bars = axes[0].bar(range(4), [cnt[k] for k in order], color=cols, edgecolor="black", linewidth=.8)
for b, k in zip(bars, order):
    axes[0].text(b.get_x() + b.get_width()/2, b.get_height(), str(cnt[k]), ha="center", va="bottom", fontsize=11)
axes[0].set_xticks(range(4)); axes[0].set_xticklabels([LBL[k] for k in order], fontsize=10)
axes[0].set_ylabel("Candidate off-target sites")
axes[0].set_title(f"A  Evidence pattern across {len(off):,} candidate sites", loc="left", pad=20)

pos = {s: sum(1 for r in off if r[f"{s}_positive"] == "1") for s in NAMES}
x = np.arange(len(NAMES))
for g, c in GRP.items():
    idx = [i for i, s in enumerate(NAMES) if ROLE[s] == g]
    axes[1].bar(x[idx], [pos[NAMES[i]] for i in idx], color=c, edgecolor="black", linewidth=0.8, label=g)
axes[1].set_xticks(x); axes[1].set_xticklabels(NAMES, rotation=45, ha="right")
axes[1].set_ylabel("Sites flagged"); axes[1].set_title("B  Flagged sites per sample", loc="left", pad=20)
for i, s in enumerate(NAMES):
    axes[1].text(i, pos[s], str(pos[s]), ha="center", va="bottom", fontsize=9)
assert_no_palette_collision({"pattern": cols[:2], "group": [C_CTRL, C_EDIT]})

leg = legend_outside(axes[1], loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
assert axes[0].get_legend() is None and not fig.legends, "图例只能属于 Panel B"
assert [t.get_text() for t in leg.get_texts()] == list(GRP)

fig.tight_layout(pad=2.0)
fig.canvas.draw(); R = fig.canvas.get_renderer()
lb, ab = leg.get_window_extent(R), axes[1].get_window_extent(R)
assert lb.x0 > ab.x1, ("图例压到 Panel B 绘图区", lb, ab)
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(f"{TMP}/Fig5.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("FIG5_1127_DONE", out_pdf)
