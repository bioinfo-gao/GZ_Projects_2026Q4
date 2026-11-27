"""P46 Fig 6 修订版（2026-11-27）。1126 原图不改。

1126 版的问题（用户 2026-11-27 指出）：Panel B 标着「Cas9-EGFP vector (not used this round)」
却给出 0.02–2.71 拷贝。该数是 11,191 bp **全长平均**，而覆盖是零碎的：
拉高它的是 ① 与报告载体 99.1% 同源的质粒骨架、② EGFP CDS；
这个载体独有的 Cas9（4,101 bp = 全长 37%）全部 10 样本仅 0.003–0.235 拷贝。
⇒ 全长均值把「别处来的同源序列」记在了这个载体头上，读者会读成"载体还是在的"。

修订：三个 panel，每个只画**归属明确**的一类序列，且都用【单载体索引】测量（方法对称）：
  A 报告载体独有：mCherry CDS          （pLPC 单载体索引）
  B Cas9-EGFP 载体独有：Cas9 CDS        （pMSCV 单载体索引）
  C 两质粒共享的骨架 ori+AmpR            （pLPC 单载体索引；v2 双载体索引下恒为 0.00，无检出力）
C 的图注写明：10 个克隆骨架接头均为 0 条，而盒式序列 18–51 条 ⇒ 骨架未整合，是游离质粒 DNA。
"""
import matplotlib; matplotlib.use("Agg")
import sys, csv, os, textwrap
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/gao/.claude/skills/scientific-figure-making/assets")
sys.path.insert(0, "/home/gao/Utilities/color_and_format/python/color_figures_with_preconfigured_sets")
from pubfig import apply_publication_style, FigureStyle, assert_visible_on_white
from palette_priority import select_palette

P = "/home/gao/projects_2026Q4/46_park_wgs"
BB = f"{P}/analysis/stage3_backbone_v3_20261127/tables"
MS = f"{P}/analysis/stage3_pMSCVonly_v3_20261127/tables"
FIG = f"{P}/analysis/stage5_figures_20261127"
TMP = f"{FIG}/tmp_png_preview"
os.makedirs(TMP, exist_ok=True)
out_pdf = f"{FIG}/Fig6_vector_copy_number.pdf"
if os.path.exists(out_pdf):
    sys.exit(f"FATAL 输出已存在，不覆盖：{out_pdf}")
apply_publication_style(FigureStyle(font_size=15, axes_linewidth=2))

st = list(csv.DictReader(open(f"{P}/scripts/0_sample_table.csv")))
NAMES = [r["safe_name"] for r in st]
ROLE = {r["safe_name"]: ("Control" if r["role"].startswith("control") else "Edited") for r in st}
pal = select_palette(4)
COL = pal["fill"] if isinstance(pal, dict) else list(pal)
C_CTRL, C_EDIT = COL[0], COL[1]
assert_visible_on_white([C_CTRL, C_EDIT], role="fill")
GRP = {"Control": C_CTRL, "Edited": C_EDIT}
def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))

cop = tsv(f"{BB}/backbone_v3_copies.tsv")
pms = tsv(f"{MS}/pMSCVonly_copies.tsv")
jx = {r["sample"]: r for r in tsv(f"{BB}/backbone_v3_junction_crosscheck.tsv")}
def pick(rows, elem): return {r["sample"]: float(r["copies_per_diploid"]) for r in rows if r["element"] == elem}
mch, cas9, ampr = pick(cop, "mCherry"), pick(pms, "Cas9"), pick(cop, "AmpR")
for d in (mch, cas9, ampr):
    assert set(d) == set(NAMES) and len(d) == 10, sorted(d)

# 与互证表一致（同一事实两处出现 ⇒ 必须比对）
nb = sum(int(jx[s]["junctions_in_backbone"]) for s in NAMES)
nc = [int(jx[s]["junctions_in_cassette"]) for s in NAMES]
assert nb == 0 and min(nc) > 0, (nb, nc)

PANELS = [
    ("A  Reporter vector, mCherry CDS", mch, "%.2f",
     "Sequence unique to the 53BP1-mCherry reporter vector. The parental line carries this reporter; "
     "values near 1 are the inherited integration."),
    ("B  Cas9-EGFP vector, Cas9 CDS", cas9, "%.3f",
     "Sequence unique to the Cas9-EGFP vector (4,101 bp). Below 0.25 copies in every clone: "
     "this plasmid is not present, consistent with Cas9 supplied as RNP."),
    ("C  Plasmid backbone (ori + AmpR)", ampr, "%.2f",
     "Backbone sequence, 99.1%% identical between the two plasmids, so it cannot be assigned to either. "
     "Zero backbone junctions in all ten clones (cassette: %d-%d) means it is not integrated: "
     "it is residual free plasmid DNA in the preparation." % (min(nc), max(nc))),
]
fig, axes = plt.subplots(1, 3, figsize=(21.0, 6.4))
x = np.arange(len(NAMES))
for ax, (ttl, vals, fmt, note) in zip(axes, PANELS):
    for g, c in GRP.items():
        idx = [i for i, s in enumerate(NAMES) if ROLE[s] == g]
        ax.bar(x[idx], [vals[NAMES[i]] for i in idx], color=c, edgecolor="black", linewidth=0.8, label=g)
    ax.set_xticks(x); ax.set_xticklabels(NAMES, rotation=45, ha="right")
    ax.set_ylabel("Copies per diploid genome"); ax.set_title(ttl, loc="left", pad=20)
    top = max(vals.values())
    ax.set_ylim(0, top * 1.22 if top > 0 else 1)
    for i, s in enumerate(NAMES):
        ax.text(i, vals[s], fmt % vals[s], ha="center", va="bottom", fontsize=8.5)

H, L = axes[0].get_legend_handles_labels()
assert L == list(GRP), L
fig.legend(H, L, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(pad=2.0)
fig.subplots_adjust(top=0.84, bottom=0.34)
fig.canvas.draw(); R = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
notes = []
for ax, (_, _, _, note) in zip(axes, PANELS):
    tb = ax.get_tightbbox(R).transformed(inv)
    t = "\n".join(textwrap.wrap(note, 54))
    notes.append((fig.text(tb.x0, tb.y0 - 0.028, t, ha="left", va="top", fontsize=9.5), ax))
fig.canvas.draw()
for txt, ax in notes:
    b = txt.get_window_extent(R).transformed(inv); tb = ax.get_tightbbox(R).transformed(inv)
    assert b.y1 < tb.y0 and b.y0 > 0, ("注释与刻度标签重叠或出界", b, tb)
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(f"{TMP}/Fig6.png", dpi=150, bbox_inches="tight"); plt.close(fig)
print("FIG6_1127_DONE", out_pdf)
print("  mCherry", {k: round(v, 2) for k, v in mch.items()})
print("  Cas9   ", {k: round(v, 3) for k, v in cas9.items()})
print("  AmpR   ", {k: round(v, 2) for k, v in ampr.items()})
