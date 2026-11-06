# LatchBio — Question Creation

**Zhen Gao** · 2026-11-05

> **Tools used.** I used Claude (Anthropic) to check my recall of the trial's published
> numbers against public sources, to tighten wording, and to format this document. The paper
> choice, the critique, the two problems, the expected answers and every grading decision are
> my own. Section 5 reports LLM testing of the two problems.

---

## 1. Paper choice

**Gilbert MR, Dignam JJ, Armstrong TS, Aldape KD, Barnholtz-Sloan JS, Won M, et al.
A randomized trial of bevacizumab for newly diagnosed glioblastoma.
*N Engl J Med.* 2014;370(8):699–708. RTOG 0825.**
https://www.nejm.org/doi/full/10.1056/NEJMoa1308573

- **What made it stay with me.** It is the cleanest example I know of a trial whose
  measuring instrument is altered by the drug being measured. Bevacizumab reduces contrast
  enhancement by restoring blood–brain barrier integrity, and progression in neuro-oncology
  was being read off contrast enhancement. The trial did not just fail to show benefit; it
  showed that the endpoint we had been trusting was partly reporting on the drug's effect on
  the barrier rather than on the tumour.

- **The result that changed how I think.** The PFS/OS discordance — PFS 10.7 vs 7.3 months
  (HR 0.79, P=0.007) with overall survival flat at 15.7 vs 16.1 months. Before this trial I
  treated radiographic progression as a noisy but unbiased proxy for tumour burden. After it
  I treat every imaging endpoint as an assay with its own failure modes, and I ask what the
  intervention does to the assay before I ask what it does to the disease.

- **Where it shows up in my own work.** In glioma treatment work I no longer accept a
  delivery or exposure argument on the basis of enhancement. When a compound is described as
  reaching a tumour because the lesion enhances, I ask which compartment that statement
  covers, because the enhancing bulk and the non-enhancing infiltrating margin are different
  pharmacological environments and only one of them seeds recurrence.

- **The open question I still find important.** The trial reported worse neurocognitive
  function and symptom burden in the bevacizumab arm, while AVAglio, asking nearly the same
  question, reported preserved quality of life. Two large randomised trials reached opposite
  conclusions about how patients actually felt. I do not think that discrepancy was ever
  resolved, and it bears directly on whether the deepened hypoxia and infiltrative relapse
  pattern seen under sustained anti-VEGF therapy carry a clinical cost that survival curves
  do not capture.

---

## 2. Result critique

**The result.** The co-primary progression-free survival analysis:
median PFS **10.7 months** (95% CI 10.0–12.2) with bevacizumab versus **7.3 months**
(95% CI 5.9–7.9) with placebo, **HR 0.79** (95% CI 0.66–0.94), **P=0.007** — an improvement
that did not meet the trial's prespecified criterion for a positive PFS result. Overall
survival was **15.7 vs 16.1 months**. N=637 randomised.

**The question the authors were trying to answer.** Whether adding bevacizumab to standard
radiotherapy plus temozolomide in newly diagnosed glioblastoma improves survival or delays
progression enough to justify first-line use.

**The conclusion they reached.** No overall survival benefit; a PFS prolongation that did not
reach the prespecified threshold; greater toxicity and, in this trial, worse neurocognitive
and symptom outcomes. First-line bevacizumab was not supported.

**What I would change about the analysis.** PFS was determined on criteria anchored to
contrast enhancement. Bevacizumab reduces contrast enhancement through restoration of
blood–brain barrier integrity, a mechanism that is independent of tumour-cell kill. The
intervention therefore acts directly on the readout used to define the endpoint, which is a
measurement-instrument confound, not ordinary endpoint noise, and blinding does not remove it
because the drug itself changes what the image shows. I would prespecify progression on
criteria that do not depend on enhancement — non-enhancing T2/FLAIR growth as a qualifying
progression event, with central blinded review and explicit rules for infiltrative, non-
enhancing spread — and report PFS under both definitions as a prespecified sensitivity
analysis. I would also report overall survival adjusted for crossover, since unblinding at
progression with crossover to bevacizumab is a competing explanation for the flat OS curve.

**Why the change matters.** The entire clinical argument for first-line bevacizumab rested on
PFS. If a measurable share of the PFS separation is produced by the drug suppressing the
signal used to detect progression, then the PFS estimate is not a biased-but-interpretable
measure of benefit — it is partly a readout of pharmacology on the barrier. Reporting PFS
under an enhancement-independent definition converts an unresolvable argument into a
measurable quantity: the difference between the two PFS estimates is itself the size of the
confound.

**The same defect propagates into the trial's toxicity conclusion, which is the part I think
is underappreciated.** Longitudinal neurocognitive and symptom assessment is performed on
patients still on study and progression-free. Who is still being assessed at a given timepoint
is therefore determined by the same enhancement-anchored progression definition. If that
definition is differentially permissive by arm, the two arms' assessment populations diverge
over follow-up, and the comparison that produced the trial's neurocognitive finding is subject
to informative censoring driven by the endpoint itself. I would report the neurocognitive
comparison restricted to a common assessment window and with explicit accounting for
differential attrition. This matters beyond RTOG 0825: AVAglio asked nearly the same question
and reached the opposite conclusion about patient-reported outcomes, and a difference in how
progression removed patients from assessment is a candidate explanation that was never, to my
knowledge, excluded.

**Taken together, this trial had no clean endpoint.** PFS was read through a channel the drug
acts on. Overall survival was measured after unblinding with crossover permitted, so it was
diluted by design. The neurocognitive comparison inherited the progression definition's
asymmetry. I want to be precise that PFS/OS discordance is not by itself the problem —
discordance is common and usually mundane, arising from crossover, post-progression therapy,
or PFS measuring something other than survival benefit. What is unusual here is that the
intervention has a direct pharmacologic action on the physical mechanism the endpoint is read
from, which is not a feature of RECIST-based endpoints in most of oncology.

**Would it strengthen, weaken, or materially change the conclusion?** It would **strengthen
the paper's negative conclusion** and change its basis. The trial's conclusion was already
that bevacizumab should not be used first-line. My proposed analysis would show that the one
seemingly positive signal was partly an artefact of the endpoint definition, rather than a
real biological benefit that merely failed to translate into survival. That is a different
and more useful statement for the field, because it indicts the endpoint rather than the
drug's efficacy alone, and it predicts the same failure for any future agent evaluated on
enhancement-based progression.

**If public data are insufficient.** The per-patient imaging reads needed to split progression
events by qualifying modality are not public. The specific analysis I would request is: for
every progression event in each arm, the modality on which progression was first declared
(new or enlarging contrast enhancement versus non-enhancing T2/FLAIR growth), from the
blinded central review. **Distinguishing outcomes:** if the proportion of progression events
first declared on non-enhancing growth is similar in both arms, the endpoint is not
differentially confounded and the PFS benefit should be read as real but survival-irrelevant.
If that proportion is substantially higher in the bevacizumab arm, enhancement-based
detection was differentially suppressed by the drug and the PFS separation is at least partly
a detection artefact.

---

## 3. Problem 1 — quantitative, grows directly from the Section 2 result

### Task prompt

```
OUTPUT FORMAT - READ THIS FIRST, IT IS GRADED.
Your reply must contain exactly one <EVAL_ANSWER> block holding a single JSON object
with the fields specified at the end of this prompt. The block is required even when you
also explain your reasoning; put all reasoning outside the block. If for any reason you
do not produce the block, restate each field name followed by its value, one per line, so
that the answer can still be scored.

A randomised phase 3 trial in newly diagnosed glioblastoma compared standard
radiotherapy plus temozolomide with the same regimen plus an anti-VEGF antibody.
Progression-free survival (PFS) was a co-primary endpoint, determined by radiographic
criteria anchored to contrast enhancement on MRI.

TABLE 1 — Reported efficacy (source: published trial report)
  Endpoint                 Anti-VEGF arm        Placebo arm
  Patients randomised      312                  309
  Median PFS (months)      10.7 (95% CI 10.0-12.2)   7.3 (95% CI 5.9-7.9)
  PFS hazard ratio         0.79 (95% CI 0.66-0.94), P = 0.007
  Median OS (months)       15.7                 16.1
  OS hazard ratio          1.13 (not significant)

TABLE 2 — Blinded central review audit of how each progression event was first
declared. These counts are a hypothetical audit constructed for this exercise; they
are not reported in the trial publication.
  Category                                          Anti-VEGF arm   Placebo arm
  Total progression events                          240             250
  First declared on new/enlarging contrast enhancement   134         213
  First declared on non-enhancing T2/FLAIR growth        106          37

Using only the tables above, do the following.

1. For each arm, compute the proportion of progression events that were first declared
   on non-enhancing T2/FLAIR growth.
2. Compute the absolute difference between those two proportions, in percentage points.
3. Choose the interpretation best supported by the data.

Which interpretation is best supported?
  A. Progression events in the two arms were detected through the same channel, so
     enhancement-anchored PFS is directly comparable between arms.
  B. The detection channel differs between arms, so comparability of enhancement-
     anchored PFS between them cannot be assumed. The tables cannot distinguish
     suppression of the enhancement signal from a genuine shift towards non-enhancing
     growth, nor establish how much of the PFS difference is a measurement effect.
  C. The difference is explained by a genuine shift towards infiltrative non-enhancing
     growth rather than by any change in detection, so PFS remains a valid comparison
     between the arms.
  D. The tables are insufficient to assess whether progression detection differs by
     treatment arm.

Return EXACTLY the following block. Do not omit the <EVAL_ANSWER> tags:

<EVAL_ANSWER>
{"nonenhancing_fraction_treated": <float>, "nonenhancing_fraction_placebo": <float>,
 "absolute_difference_pp": <float>, "answer": "<letter>"}
</EVAL_ANSWER>

Report the two fractions as proportions rounded to two decimals, and the difference in
percentage points rounded to one decimal, where <letter> is A, B, C, or D. Provide a
one-line rationale outside the answer tags.
```

### Inputs

| Input | Source | Processing |
| :--- | :--- | :--- |
| Table 1 | Gilbert et al., *NEJM* 2014;370:699–708 (RTOG 0825) | Transcribed from the published report. No recomputation. |
| Table 2 | **Constructed by me for this evaluation** | The modality split is not published. I constructed counts consistent with the trial's event totals and with the known direction of the effect, and the prompt states explicitly that these are hypothetical. The problem is self-contained: no value outside the two tables is needed. |

Nothing confidential or proprietary is used. The model does not need to browse or recognise
the paper; every number required is in the prompt.

### Expected answer

```json
{"nonenhancing_fraction_treated": 0.44, "nonenhancing_fraction_placebo": 0.15,
 "absolute_difference_pp": 29.4, "answer": "B"}
```

Derivation, entirely from Table 2: 106 / 240 = 0.4417 → **0.44**.
37 / 250 = 0.1480 → **0.15**. Difference = 44.17 − 14.80 = **29.4 percentage points**
(a 2.98-fold higher share). Interpretation: if treatment did not alter how
progression came to attention, the share of events first visible as non-enhancing growth
should be similar in both arms. It is three times higher under anti-VEGF therapy, so the two
arms' progression events did not arrive through the same channel, and an endpoint anchored to
contrast enhancement is therefore not measuring the same thing in the two arms → **B**.

**The mechanism deliberately cannot be resolved from these tables, and that is the point of
the item.** Two explanations produce exactly this table: the drug suppresses the enhancement
that would otherwise have revealed progression, or the drug genuinely shifts the tumour
towards infiltrative non-enhancing growth. Both are documented for anti-VEGF therapy in
glioblastoma. A response that picks one as established has over-read the evidence. The expert
move is to notice that **the choice does not matter for the endpoint conclusion**: under
either mechanism the enhancement-anchored PFS clock runs differently in the two arms, so the
comparison is compromised either way.

This is coherent with the flat OS in Table 1, but the flat OS is not itself evidence for B.
PFS and OS diverge for mundane reasons — crossover, post-progression therapy, PFS simply
measuring something else — in many trials, and discordance alone proves nothing. What makes
this trial different is that the drug acts directly on the physical mechanism the endpoint is
read from.

### Grading rule

Parse the single `<EVAL_ANSWER>` block as JSON. Grade only these four fields; ignore the
rationale text and all prose style.

| Field | Accepted | Tolerance |
| :--- | :---: | :---: |
| `nonenhancing_fraction_treated` | 0.44 | ±0.01 (accepts 0.44 and 0.4417; also accept 44 or 44.2 if the model reports percent, since the quantity is unambiguous) |
| `nonenhancing_fraction_placebo` | 0.15 | ±0.01 (accepts 0.15 and 0.148; also accept 14.8) |
| `absolute_difference_pp` | 29.4 | ±1.0 (accepts 29.0–30.0, covering 0.29 reported as a proportion ×100) |
| `answer` | `B` | exact, case-insensitive |

**Scientifically defensible alternatives that also pass.** A response that answers **B** and
additionally notes that crossover at progression contributes to the flat overall survival
passes; that is real and does not conflict with B. A response reporting the two fractions as
percentages rather than proportions passes, as specified above. A response that answers **B**
and argues that suppression of the enhancement signal is the *more likely* of the two
mechanisms also passes, provided it does not claim these tables establish it — naming a
preferred mechanism is a scientific judgment, asserting that this evidence proves it is not.

**Conditions that fail.** `answer` other than B. Any numeric field outside tolerance.
Reporting the fraction of *enhancement*-declared events (0.56 / 0.85) in the non-enhancing
fields — a sign-flip, not a rounding difference. Deriving the difference from the PFS medians
(3.4 months) rather than from Table 2. Answering **C**, which commits to one of two mechanisms
these data cannot separate *and* draws the opposite endpoint conclusion from it. A response
whose graded fields cannot be recovered at all — see Section 5 for how a missing
`<EVAL_ANSWER>` block is handled.

### Scientific skill

This tests whether a model distinguishes an endpoint confound from an endpoint effect — that
is, whether it notices that the intervention acts on the measurement channel and not only on
the disease — and then whether it knows how far the evidence in front of it reaches. Two
mechanisms fit the table and these data cannot separate them; the useful judgment is that the
endpoint is compromised under either, so resolving the mechanism is not a prerequisite for the
decision. The arithmetic is deliberately trivial; the difficulty is knowing which quantity
to compute. In practice this is the judgment that decides whether a PFS signal is worth
taking into a phase 3 design, and in neuro-oncology it has already produced a decade of
unreproducible PFS benefit. My own work on glioma treatment is where I learned to ask which
compartment an imaging statement covers before accepting it as an exposure claim.

### Plausible wrong approach

Compute 10.7 − 7.3 = 3.4 months, observe HR 0.79 with P=0.007, and conclude that the drug
delays progression — answer **A**. This is a reasonable-looking shortcut: it uses the
headline numbers, it is the comparison the trial itself foregrounded, and it is arithmetically
correct. It fails because it never asks whether the two arms' progression events were detected
through the same channel. A second shortcut is to answer **D** on the grounds that per-patient
data are absent; Table 2 supplies exactly the aggregate needed, so this is over-caution rather
than rigour.

A third, and the one I expect to catch the better-informed respondent, is **C**. It is well
documented that anti-VEGF therapy shifts glioblastoma towards infiltrative, non-enhancing
growth, so a reader who knows that literature can account for Table 2 without invoking any
change in detection — and that account is not wrong. C fails for two reasons. These tables
cannot separate that mechanism from suppression of the enhancement signal, so asserting it is
an over-read; and the conclusion C draws from it is backwards, because a treatment-induced
change in how the tumour becomes visible compromises an enhancement-anchored endpoint just as
thoroughly as signal suppression does. C is the trap for domain knowledge applied without
first asking what the data in hand can discriminate. Neither failure depends on missing information, ambiguous wording, or formatting:
the discriminating evidence is printed in the prompt, and a correct analysis uses it.

---

## 4. Problem 2 — biological interpretation and development decision

### Task prompt

```
OUTPUT FORMAT - READ THIS FIRST, IT IS GRADED.
Your reply must contain exactly one <EVAL_ANSWER> block holding a single JSON object
with the fields specified at the end of this prompt. The block is required even when you
also explain your reasoning; put all reasoning outside the block. If for any reason you
do not produce the block, restate each field name followed by its value, one per line, so
that the answer can still be scored.

In newly diagnosed glioblastoma, a randomised phase 3 trial added an anti-VEGF antibody
to radiotherapy plus temozolomide. Progression-free survival improved (10.7 vs 7.3
months, HR 0.79) while overall survival did not (15.7 vs 16.1 months). The trial also
reported worse neurocognitive function and greater symptom burden in the anti-VEGF arm.

Background you may rely on, and nothing else:
  - Glioblastoma comprises a contrast-enhancing bulk and a non-enhancing infiltrating
    margin. Recurrence arises from the infiltrating margin.
  - Anti-VEGF therapy acts on tumour vasculature. Answer for sustained therapy, beyond
    any transient vascular-normalisation window.
  - Temozolomide is an alkylating agent whose lesions become lethal principally during
    DNA replication.

A colleague proposes continuing anti-VEGF therapy first-line in all patients, arguing
that improved perfusion will increase temozolomide delivery to the infiltrating margin
and that the PFS gain will eventually translate into survival with longer follow-up.

Evaluate that proposal.

Which single statement is best supported?
  A. The proposal is sound; longer follow-up is the main missing ingredient.
  B. The PFS gain is partly an artefact of an enhancement-anchored endpoint, and
     sustained anti-VEGF therapy does not increase, and may reduce, cytotoxic exposure
     and cytotoxic sensitivity at the infiltrating margin.
  C. The proposal fails only because the trial was underpowered for overall survival.
  D. Anti-VEGF therapy increases temozolomide efficacy but the benefit is offset by
     unrelated toxicity.

Return EXACTLY the following block. Do not omit the <EVAL_ANSWER> tags:

<EVAL_ANSWER>
{"effect_on_margin_drug_exposure": "<increases|decreases|unchanged>",
 "effect_on_margin_cytotoxic_sensitivity": "<increases|decreases|no_effect>",
 "recommended_endpoint_change": "<keep_enhancement_pfs|adopt_nonenhancing_inclusive_criteria|use_os_only|use_response_rate>",
 "answer": "<letter>"}
</EVAL_ANSWER>

Provide a one-line rationale outside the answer tags.
```

### Inputs

| Input | Source | Processing |
| :--- | :--- | :--- |
| Efficacy and neurocognitive/symptom direction | Gilbert et al., *NEJM* 2014;370:699–708 | Transcribed. |
| Three background statements | Standard glioblastoma biology and alkylator pharmacology | Stated in the prompt so the problem is self-contained. I deliberately do **not** state that the margin is hypoxic, that hypoxia suppresses cycling, or that the barrier is restored under anti-VEGF therapy — those are the inferences being tested. |

### Expected answer

```json
{"effect_on_margin_drug_exposure": "decreases",
 "effect_on_margin_cytotoxic_sensitivity": "decreases",
 "recommended_endpoint_change": "adopt_nonenhancing_inclusive_criteria",
 "answer": "B"}
```

Derivation from the supplied material. Sustained anti-VEGF therapy restores barrier
integrity — this is the same mechanism that reduces contrast enhancement, which the prompt's
first paragraph implies by pairing a PFS gain with no survival gain — so drug access to the
non-enhancing margin **decreases** rather than increases. Vessel pruning under sustained
therapy deepens hypoxia at that margin; hypoxic cells there are largely not cycling; and the
prompt states that temozolomide's lesions become lethal principally during replication.
Sensitivity at the margin therefore **decreases**. Because the margin is where recurrence
originates, neither change can be rescued by longer follow-up, which is why PFS did not
translate into OS. The endpoint implication follows from the same mechanism: a drug that
suppresses enhancement must not be evaluated on an enhancement-anchored progression
definition, so **adopt non-enhancing-inclusive criteria**. Hence **B**.

### Grading rule

Parse the single `<EVAL_ANSWER>` block as JSON. Four categorical fields, exact string match,
case-insensitive, no tolerance to set.

| Field | Correct |
| :--- | :---: |
| `effect_on_margin_drug_exposure` | `decreases` |
| `effect_on_margin_cytotoxic_sensitivity` | `decreases` |
| `recommended_endpoint_change` | `adopt_nonenhancing_inclusive_criteria` |
| `answer` | `B` |

**Scientifically defensible alternatives that also pass.** `unchanged` for
`effect_on_margin_drug_exposure` passes **if and only if** the rationale states that the
margin barrier was already intact so anti-VEGF therapy cannot reduce access further; that is a
coherent reading of the same biology and does not change the conclusion. A rationale that
additionally attributes part of the flat OS to crossover at progression passes, since
crossover is a genuine confounder in this trial and is not mutually exclusive with B.
`use_os_only` does **not** pass: overall survival is the right arbiter but is not an endpoint
*change* that fixes progression assessment, and a trial still needs an interim progression
readout.

**Conditions that fail.** `answer` other than B. `increases` on either effect field — that is
the vascular-normalisation frame applied to sustained therapy and is the primary error being
tested. Any field absent, or a missing/malformed/duplicated `<EVAL_ANSWER>` block.

### Scientific skill

This tests whether a model holds the direction of a well-published effect in the specific
disease, rather than the direction asserted by the general mechanism literature. Vascular
normalisation predicts improved delivery and improved sensitivity; in glioblastoma under
sustained therapy both go the other way. It also tests a real development judgment that the
brief names explicitly — choosing the right metric — because the correct endpoint decision
follows from the same mechanism as the biology answer, and a model that gets the biology by
retrieval rather than by mechanism will not connect the two. My glioma treatment work is
where this matters: the compartment that drives recurrence is the one the imaging endpoint is
blindest to.

### Plausible wrong approach

Apply the vascular-normalisation literature directly: anti-VEGF therapy normalises vessels,
lowers interstitial pressure, improves perfusion and oxygenation, therefore delivery and
radio/chemosensitivity both rise, therefore the colleague is right and only follow-up is
missing — answer **A** with `increases`/`increases`. This is the published rationale for
combining anti-VEGF agents with cytotoxics, so a model retrieving consensus language will
produce it confidently. It fails because the normalisation window is transient and the prompt
scopes the question to sustained therapy, and because it cannot explain a PFS/OS discordance
in the direction observed. A second shortcut is **C**, attributing the failure to power; the
trial randomised 637 patients and the OS hazard ratio is 1.13, pointing away from benefit
rather than towards an undetected one. Neither failure relies on missing information or
ambiguous wording: the scoping sentence and the three background statements are in the prompt.

---

## 5. LLM outputs appendix

Five runs across two products, covering both problems on both products, and two versions of
the Problem 1 option set. Complete responses
are in `appendix/`, reproduced verbatim. In every run the operator appended the instruction
"in detailed English" to the prompt; that modification is recorded in each response file.

| Run | Model / product | Date | Problem | Option set | Graded fields | Structured block | Grade |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | OpenAI ChatGPT, free tier | 2026-11-05 | 1 | v1 | all correct | present, well-formed | **PASS** |
| 2 | Google Gemini, free tier | 2026-11-05 | 2 | — | all correct | present, well-formed | **PASS** |
| 3 | OpenAI ChatGPT, free tier | 2026-11-05 | 1 | v2 | all correct | present, well-formed | **PASS** |
| 4 | Google Gemini, free tier | 2026-11-05 | 1 | v2 | all correct | not retained (see below) | not scored |
| 5 | OpenAI ChatGPT, free tier | 2026-11-05 | 2 | — | all correct | present, well-formed | **PASS** |

**Run 4's structured block was not captured when the response was collected.** Its prose
derivation and its selection of option B are complete and are reproduced in the appendix, so
the scientific result stands, but I cannot assert what its block contained and I do not count
it in the format-compliance tally. It is listed as unscored rather than failed.

### What the models did well

**Every run reached the correct scientific answer on every graded field** — 0.44, 0.15, 29.4
and B on Problem 1; `decreases`, `decreases`, `adopt_nonenhancing_inclusive_criteria` and B on
Problem 2 — and the four retained blocks were well-formed JSON inside a single
`<EVAL_ANSWER>` pair, parsed by the grader without intervention. Both problems were answered
correctly by both products, so neither item's answer depends on one model's idiosyncrasy. Three of the four carried the
unrounded fractions through before rounding the difference, which is the arithmetic the item
asks for.

No run took any of the three plausible wrong approaches. None substituted the 3.4-month PFS
difference for the quantity asked. None retreated to D for want of per-patient data. **C, the
distractor I expected to be strongest, trapped neither model**: ChatGPT rejected it explicitly
as assigning a biological explanation the tables cannot establish, and Gemini passed over it.
That is an honest negative result, discussed below.

On Problem 2, Gemini supplied unprompted both inferences the task deliberately withholds —
that sustained anti-VEGF therapy deepens hypoxia at the margin, and that hypoxic cells leave
cycle and so escape a replication-dependent alkylator. That confirms the task is answerable
from the stated premises alone, which is the fairness property the withholding was meant to
preserve.

### Where the reasoning was sharper than my own item

The most useful outcome of testing was not a model error. It was a model identifying a defect
in my option wording, twice, on successive versions.

**Run 1** answered B under the v1 option set and then qualified it: the audit records which
imaging finding first *triggered* the progression call, so establishing that PFS actually
overstated benefit would additionally require evidence about when the tumour truly progressed,
and a treatment-induced shift in growth pattern could itself move the proportions. My v1 option
B concluded that PFS "overstates the treatment benefit". These tables do not establish that.
The model was right and the option over-claimed.

**Run 3**, on the rewritten v2 set, selected B and again marked the residual over-claim: v2's
assertion that either mechanism "undermines comparability" is too categorical, because a
genuine shift in growth pattern could still permit a valid comparison *if the assessment
criteria capture both enhancing and non-enhancing progression with comparable sensitivity and
timing*. Its one-line rationale is precise: "The detection channel differs between arms,
raising concerns about enhancement-anchored PFS comparability, although these tables cannot
establish the mechanism or prove that comparability is invalid." That condition was absent
from my option, and it should not have been.

The same behaviour appeared on Problem 2. Run 5 selected B and added that the mechanisms "and
an artefactual PFS component cannot be established from the supplied information alone." Here
I made no change: Problem 2's option B is already phrased conditionally — it says sustained
therapy "does not increase, and may reduce" exposure and sensitivity — so the model's caution
and the option are compatible. I record it because the pattern is consistent across both items
and is a property of the product rather than of one question's wording.

Both Problem 1 corrections were adopted. This is the behaviour I would want from a colleague
reviewing a draft eval, and it is the reason I think testing on models is worth doing even when
the models pass: a model that answers correctly can still tell you that your question claims
more than your evidence.

### Revisions made after testing

1. **Problem 1 option set, v1 → v2 (after run 1).** B was narrowed from "PFS overstates the
   treatment benefit" to a claim about comparability. A was reduced from a composite claim to a
   single one. A new distractor C was added — that the shift is genuine tumour biology rather
   than detection — because that explanation is real, documented in glioblastoma, and fits the
   table exactly.

2. **Problem 1 option set, v2 → v3 (after run 3).** B now reads that comparability "cannot be
   assumed" rather than that PFS "is not comparable", and states that the tables cannot
   establish how much of the PFS difference is a measurement effect. B remains the
   best-supported option — A is refuted by the counts, C asserts a mechanism these data cannot
   isolate, D is refuted because the channel difference is demonstrable — but it now claims
   only what the evidence carries. **This is the version submitted**, and it has not itself
   been tested; see below.

3. **Output contract, tightened as a precaution, not in response to a failure.** The format
   instruction was moved to the top of both prompts and marked as graded, and a fallback was
   added asking that each field name be restated with its value if the block is omitted. All
   three retained runs complied without needing either, so these changes are untested
   insurance rather than a demonstrated fix.

4. **Grader fallback, same status.** `scripts/grade.py` will recover graded fields from prose
   when no `<EVAL_ANSWER>` block is present, and only when a field name appears exactly once,
   so an ambiguous response fails rather than being guessed at. No tested run exercised this
   path. I kept it because an evaluation should not fail a correct scientific answer on a
   wrapper, which the brief's instruction to grade the terminal scientific answer implies, but
   it is a precaution and I am not claiming evidence for it.

5. **Not revised: the expected answers, and Problem 2's option set.** Five independent runs
   across two products reached the expected answers by sound mechanism, with each problem
   answered correctly on both products. Run 5's caution on Problem 2 did not require a change
   because that option is already conditionally phrased.

### What testing did not establish

**The submitted v3 option set has not been tested.** Runs 1 and 3 tested v1 and v2; v3 adopts
run 3's correction. Re-running Problem 1 on v3 is the first thing I would do with more time,
and I would expect run 3's qualification to disappear, which is the check that the correction
landed.

**The distractors remain unconfirmed.** No run fell into A, C or D on Problem 1, and Problem 2's
`increases`/`increases` vascular-normalisation failure mode was never produced. I have evidence
that both items are answerable and that their correct answers are reachable by independent
reasoning on two products, but no direct evidence that the wrong approaches attract a model.
Confirming that needs a weaker model, or the same models given less context.

Run 4's structured output was not retained, so the format-compliance observation rests on four
runs rather than five.

**A correction to my own process, stated because it affected this section.** An earlier draft of
this appendix reported that none of the four runs produced a structured answer block and built
an argument about models treating format requirements as advisory. That was wrong: the blocks
had been produced and were simply missing from the material I was working from. Three are now
reproduced in the appendix and all three parse. I note it because the same failure mode —
drawing a confident conclusion from an incomplete record of what a model returned — is exactly
what an evaluation appendix exists to prevent, and it would have been invisible to a reader.

## Grading code

A reference implementation of both grading rules is provided as `scripts/grade.py`.
It reads a model response on stdin, extracts the `<EVAL_ANSWER>` block, and applies the
tolerances above. The rules are simple enough to implement by hand; the script exists so that
no reader has to make a new scientific judgment to reproduce a score.
