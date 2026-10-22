#!/usr/bin/env python3
# transfer-answer-intended
"""Draft (NOT send) the reply to Xin's 2026-10-14 reduced-scope request, via IMAP APPEND to Gmail Drafts."""
import html, imaplib, json, re, time
from email.message import EmailMessage

TO = "Lin, Xin <Xin.Lin@childrens.harvard.edu>"
CC = ("Singh, Vedita <Vedita.Singh@childrens.harvard.edu>, Jing Shi <j.shi@athenomics.com>, "
      "Jianjun Sun <jianjun@athenomics.com>, Li, Julia <Julia.Li@childrens.harvard.edu>")
SUBJECT = "Re: Additional analysis at chr11:114,604,000-114,625,500 - feasibility, time estimate and budget"
BODY = """Dear Xin,

Thank you for narrowing the scope. I answer your three points in order.

1. Quote for Analysis 1 on the 8 samples
I read your list as S01 (DLD1_pop), S02 (DLD1_c1_ctrl), S04 (DLD1_c3_dox), S05 (DLD1_c4_dox), S06 (GM12878_pop), S07 (GM12878_c5), S08 (GM12878_c6) and S09 (GM12878_c7). This leaves out S03 (DLD1_c2_ctrl), S10 (GM12878_c8) and all Mutu I and DG75 samples; please tell me if that is not what you intended.

The earlier range of 12 to 22 hours for Analysis 1 was given before the number of samples and junctions was fixed; the figures below are for the defined scope. The work has two parts, because the first scales with the number of samples and the second with the number of junctions.

Read rescue for all candidate junctions (scales with the number of samples)
  - Shared by all samples (setting up the rescue and classification workflow for the locus, the repeat-content assessment, and write-up): 6 to 7 hours.
  - Per sample (re-extraction of all reads at the locus, junction-spanning read counts before and after rescue, a per-sample junction table): 0.5 to 1 hour each.
  | Sample | Candidate junctions in the target region (current analysis) | Incremental time |
  | S01 DLD1_pop | 21 | 0.5-1 h |
  | S02 DLD1_c1_ctrl | 14 | 0.5-1 h |
  | S04 DLD1_c3_dox | 18 | 0.5-1 h |
  | S05 DLD1_c4_dox | 12 | 0.5-1 h |
  | S06 GM12878_pop | 22 | 0.5-1 h |
  | S07 GM12878_c5 | 24 | 0.5-1 h |
  | S08 GM12878_c6 | 13 | 0.5-1 h |
  | S09 GM12878_c7 | 26 | 0.5-1 h |
  Candidate junctions are events with a junction end inside the target region that have at least one supporting read in that sample (150 in total).
Read rescue total for the 8 samples: 10 to 15 hours, $1,200 to $1,800.

Manual breakpoint curation (scales with the number of junctions)
Local assembly, microhomology or inserted sequence, and the call of "well supported" versus "ambiguous because of repeat content": 0.5 to 1 hour per junction. It applies only to junctions that keep their support after rescue, and that number is not known until the rescue has been run. In the current analysis most of the 150 candidate junctions rest on a single read, so many may not keep their support. I suggest that I do the read rescue first and then tell you how many junctions kept their support; you then choose how many to curate within your budget. For example, curating 10 junctions adds 5 to 10 hours ($600 to $1,200); read rescue plus manual curation of 10 junctions is 15 to 25 hours ($1,800 to $3,000).

2. Where to trim further
  - Each sample removed from the read rescue saves 0.5 to 1 hour ($60 to $120); the shared 6 to 7 hours do not change.
  - Each junction that is not curated saves 0.5 to 1 hour.
  - For example, Read rescue for the two +Dox clones and their no-dox comparators only (S04, S05, S01, S02) is 8 to 11 hours ($960 to $1,320).

3. Turnaround
Read rescue: about 2 weeks after approval; manual curation of 10 junctions adds about 1 week. If you tell me your deadline, I will say what fits before it.

4. BAM files now
Yes. The BAM files with indexes for all 18 samples (about 850 GB) will be shared in a Dropbox folder. A transfer of this size takes time to synchronise, so I will tell you when every file is in the folder, together with a checksum file (MD5) so that you can confirm that the files arrived intact.

I will not start the analysis until you confirm the scope and the approved amount.

Best regards,

Zhen Gao, PhD
Principal Bioinformatics Scientist
Athenomics
"""
BOLD = [
    "Read rescue total for the 8 samples: 10 to 15 hours, $1,200 to $1,800.",
    "read rescue plus manual curation of 10 junctions is 15 to 25 hours ($1,800 to $3,000).",
    "Read rescue for the two +Dox clones and their no-dox comparators only (S04, S05, S01, S02) is 8 to 11 hours ($960 to $1,320).",
    "I will not start the analysis until you confirm the scope and the approved amount.",
]
assert not re.search(r"[^\x00-\x7f]", BODY), "non-ASCII char in body"
assert sum([21, 14, 18, 12, 22, 24, 13, 26]) == 150

cfg = json.load(open("/home/gao/.claude.json"))["mcpServers"]["gmail"]["env"]
addr, pwd = cfg["GMAIL_ADDRESS"], cfg["GMAIL_APP_PASSWORD"]
msg = EmailMessage()
msg["From"], msg["To"], msg["Cc"], msg["Subject"] = addr, TO, CC, SUBJECT
msg.set_content(BODY)
H = html.escape(BODY)
for k in BOLD:
    ek = html.escape(k); assert H.count(ek) == 1, k[:40]
    H = H.replace(ek, "<b>" + ek + "</b>")
msg.add_alternative('<div style="white-space:pre-wrap;font-family:Arial,sans-serif;font-size:14px">' + H + "</div>", subtype="html")

M = imaplib.IMAP4_SSL("imap.gmail.com"); M.login(addr, pwd)
def n():
    typ, d = M.select('"[Gmail]/Drafts"'); return int(d[0])
before = n()
M.append('"[Gmail]/Drafts"', r"(\Draft)", imaplib.Time2Internaldate(time.time()), bytes(msg))
after = n()
print("Drafts before/after:", before, after)
assert after == before + 1
M.logout()
