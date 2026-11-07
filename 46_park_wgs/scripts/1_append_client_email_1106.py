#!/usr/bin/env python3
"""P46 —— 1106 去信：两个开跑前必须定下的问题（Gmail Drafts, APPEND only，不发送）。

问什么、为什么（详见 docs/P46_research_plan_and_review_1028.md §9）：
  Q1 本轮 Cas9 以何种形式进入细胞（RNP / 瞬时质粒 / 克隆自带整合的 Cas9 构建体），
     若是新构建体需要图谱 —— 决定 Stage 3 的 vectors.fa 组成，不给参考就查不到。
  Q2 样本分组确认：问卷用了两套编号（Sample 1–7 与 1–10），我方读作偏移 3 ⇒
     Ctr_1–3 对照 / Sample 1–4 单等位 / Sample 5–7 双等位。判错组别会贯穿整个组间比较。
  另附一句期望管理：十个克隆各自单细胞扩增，各带私有变异；三个对照的作用是估计这一本底。
     —— 先说清楚，避免结果阶段把任一私有变异读成「CRISPR 造成的」。

禁止项（与 hooks/checks/check_client_email.py 同一套判据）：
  长读长任何形式、认错/道歉词、数据传输（Dropbox/upload/fastq/还在等）、
  建库方式征询（PCR-free）、CJK 全角字符。
⚠ 本脚本只 APPEND 进 Drafts，由用户本人逐字过目后自行发出；任何情况下不执行发送。
"""
import json, imaplib, time, re
from email.message import EmailMessage

TO   = "Yeji.Park@childrens.harvard.edu"
CC   = ("Jing Shi <j.shi@athenomics.com>, Kuo-Chan Chung <kchung@athenomics.com>, "
        "Jianjun Sun <jianjun@athenomics.com>")
SUBJ = "Re: Inquiry About Whole Genome Sequencing Services [EXTERNAL]"

# 客户 2026-11-06 18:11 UTC 那封（带样本信息文件的）—— IMAP 实测取得
IRT  = "<SA1PR16MB6978A079A13EDD94617380AFAD762@SA1PR16MB6978.namprd16.prod.outlook.com>"
REFS = ("<SA1PR16MB697819243AE6B4CDE76753C3ADAC2@SA1PR16MB6978.namprd16.prod.outlook.com> "
        "<SA1PR16MB69787009E9FAF7F14F9FFF65AD762@SA1PR16MB6978.namprd16.prod.outlook.com> "
        + IRT)

BODY = """Hi Yeji,

Thank you for the detailed answers - they cover nearly everything we need. Two
points are worth settling before the analysis starts, because each one changes
what we look for.

First, how Cas9 reached the cells in this round. You wrote that the lines carry
53BP1-FFR-mCherry and that the deletion was induced with the dual caged guide
RNAs under light control, and that tells us how the cut was controlled. What it
does not tell us is the form Cas9 itself arrived in: as protein together with
the guides, from a plasmid transfected for that purpose, or from a Cas9
construct already integrated in the line. We ask because one part of the
analysis asks whether any exogenous sequence has joined the genome, and that
search can only find sequence we give it as a reference. From the previous
round we hold the maps for mCherry-BP1-2 pLPC-Puro and for pMSCVpuro-Cas9-EGFP,
and we will screen for both. If a different Cas9 construct was used this time,
please send its map. If Cas9 came in as protein and left nothing behind, then
there is nothing to add, and that is equally useful to know.

Second, a check on the sample groups. The information you sent uses two
numbering schemes - Sample 1 to 7 in some answers, 1 to 10 in others - and we
read them as the same ten samples with an offset of three. That gives:

  Ctr_1, Ctr_2, Ctr_3   no deletion seen by gel or Sanger sequencing
  Sample 1 to Sample 4  deletion on one of the two alleles
  Sample 5 to Sample 7  deletion on both alleles

Please confirm this, or correct whichever line does not match your records. The
grouping is what the genome-wide comparison rests on: the three controls set the
baseline the seven edited clones are measured against, so a sample sitting in
the other group would carry through every comparison downstream.

One point on how the results will read, worth saying before rather than after.
All ten lines are separately expanded single-cell clones, so each one carries
changes private to its own expansion history, quite apart from anything the
CRISPR system did. The three controls are what let us estimate that background,
and we will report the comparison at that level - whether the edited clones
carry more new changes than the controls do, and whether any of them fall at
sites the guides could reach. Single changes in a single clone cannot be traced
back to the editing on sequence alone, and we will not present them that way.

Best regards,

Zhen Gao, PhD
Principal Bioinformatics Scientist
Athenomics
"""

# ---- hard assertions -------------------------------------------------------
low = BODY.lower()
for w in ("fix", "mistake", "error", "wrong", "misleading", "sorry", "apolog",
          "you are right", "you're right", "ours to fix"):
    assert w not in low, f"BANNED (self-repudiation): {w}"
for w in ("nanopore", "pacbio", "long read", "long-read", "ont "):
    assert w not in low, f"BANNED (long-read platform): {w}"
for w in ("upload", "dropbox", "fastq", "data received", "still waiting"):
    assert w not in low, f"BANNED (data transfer): {w}"
assert "pcr-free" not in low, "BANNED (library prep is our decision)"
assert not re.search(r"[　-〿一-鿿＀-￯]", BODY), "CJK/fullwidth char"
assert "Zhen Gao, PhD" in BODY and "Athenomics" in BODY, "signature missing"
# 正向断言：两个问题都必须在，且分组表三行齐全
assert "Cas9" in BODY and "map" in BODY, "Q1 (Cas9 delivery / construct map) missing"
assert BODY.count("Sample 1") >= 1 and "both alleles" in BODY and "one of the two alleles" in BODY, \
    "Q2 (group confirmation table) incomplete"
# 不得把推断写成客户已确认的事实（断言前先把换行压成单空格，否则会被折行骗过）
FLAT = re.sub(r"\s+", " ", BODY)
assert "we read them as" in FLAT, "numbering alignment must be stated as our reading"
assert "offset of three" in FLAT, "must state the numbering offset we inferred"

msg = EmailMessage()
msg["To"], msg["Cc"], msg["Subject"] = TO, CC, SUBJ
msg["In-Reply-To"], msg["References"] = IRT, REFS
msg.set_content(BODY)

cfg = json.load(open("/home/gao/.claude.json"))
env = cfg["mcpServers"]["gmail"]["env"]
M = imaplib.IMAP4_SSL("imap.gmail.com")
M.login(env["GMAIL_ADDRESS"], env["GMAIL_APP_PASSWORD"])
M.select('"[Gmail]/Drafts"', readonly=True)
before = len(M.search(None, "TO", f'"{TO}"')[1][0].split())
M.append('"[Gmail]/Drafts"', r'(\Draft)', imaplib.Time2Internaldate(time.time()), bytes(msg))
M.select('"[Gmail]/Drafts"', readonly=True)
ids = M.search(None, "TO", f'"{TO}"')[1][0].split()
M.logout()
print(f"Drafts to client: before={before}  after={len(ids)}  uids={[i.decode() for i in ids]}")
assert len(ids) == before + 1, "APPEND did not land"
print("OK: draft landed in [Gmail]/Drafts (NOT sent)")
