# LatchBio — Question Creation — Zhen Gao

Submitted 2026-11-05. Paper: Gilbert et al., *N Engl J Med.* 2014;370:699–708 (RTOG 0825).

## What to read

**`LatchBio_Question_Creation_Zhen_Gao.md` is the submission.** It is self-contained: all five
required sections are in it, and the grading rules are stated in full in prose, so nothing
else in this archive is needed to read or reproduce the evaluation. The remaining files are
the supporting data and code the brief invites.

## Contents

| Path | What it is |
| :--- | :--- |
| `LatchBio_Question_Creation_Zhen_Gao.md` | **The submission.** Sections 1–5: paper choice, result critique, Problem 1, Problem 2, LLM outputs appendix |
| `scripts/grade.py` | Reference grader for both problems. Python 3 standard library only, no dependencies. `python3 -I scripts/grade.py p1 < response.txt` (or `p2`); exit 0 = pass |
| `appendix/prompt_problem1.txt` | Problem 1 prompt exactly as given to the model. Identical to the Task prompt block in Section 3 |
| `appendix/prompt_problem2.txt` | Problem 2 prompt exactly as given to the model. Identical to Section 4 |
| `appendix/p1_openai_run1.txt` | Complete response. ChatGPT, Problem 1, **original** prompt |
| `appendix/p1_openai_run2.txt` | Complete response. ChatGPT, Problem 1, **revised** prompt. This run identified a real over-claim in my option B; see Section 5 |
| `appendix/p1_gemini_run1.txt` | Complete response. Gemini, Problem 1, **revised** prompt. Structured block not retained; see Section 5 |
| `appendix/p2_gemini_run1.txt` | Complete response. Gemini, Problem 2 |
| `appendix/p2_openai_run1.txt` | Complete response. ChatGPT, Problem 2 |

All five are verbatim and unedited. Each begins with a header line giving model, date, problem,
run number and which prompt version produced it; everything after that line is the model's
output as received.

## Reproducing the reported grades

```bash
for f in appendix/p1_*.txt; do echo "== $f"; python3 -I scripts/grade.py p1 < "$f"; done
for f in appendix/p2_*.txt; do echo "== $f"; python3 -I scripts/grade.py p2 < "$f"; done
```

Four of the five runs retained their `<EVAL_ANSWER>` block; all four parse and score **PASS**.
The fifth (`p1_gemini_run1.txt`) is scientifically correct in its prose but its block was not
captured when the response was collected, so the grader reports it as unrecoverable and
Section 5 lists it as unscored rather than failed. **Every run reached the correct scientific
answer on every graded field, and each problem was answered correctly on both products.**

**Prompt modification.** In all runs the operator appended "in detailed English" to the prompt.
This is recorded in each response file.

**Note on the prompt files.** `appendix/prompt_problem1.txt` and `prompt_problem2.txt` are the
final versions. The Problem 1 option set went through three versions during testing; runs 1
used v1, runs 3–4 used v2, and the submitted file is v3, which incorporates the correction run
3 identified. Section 5 reports each revision and what prompted it.

## One constructed input, declared

Problem 1 supplies two tables. Table 1 is transcribed from the published trial report. **Table
2 is constructed by me for this exercise** — the split of progression events by the imaging
modality on which progression was first declared is not published. The prompt given to the
model states this explicitly, and Section 3's Inputs subsection explains the construction and
why it is the analysis I would request of the trialists (Section 2).

## Tools

Claude (Anthropic) was used to verify my recall of the trial's published figures against
public sources, to tighten wording, and to format this document. The paper choice, the
critique, both problems, the expected answers and every grading decision are my own.
