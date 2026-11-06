# glioma_margin_rate_limit_v3

> ⚠ **本文是渲染产物，不是真身。** 内容真身在 `evals/glioma_margin_rate_limit_v3.json`；要改内容请改那个 JSON，然后重新渲染本文。
> 渲染时间：2026-11-05（渲染器：`scripts/2_render_eval.py`）

**Data:** https://github.com/bioinformatics-gao/CASE-BBB-prediction-Data

## Grader — MultipleChoice

| 字段 | 标准答案 | 容差 |
| :--- | :---: | :---: |
| `primary_determinant_of_margin_response` | `proliferative_state_of_the_target_cells` | (exact) |
| `agent_selection_by_bbb_label` | `lomustine_preferred_but_this_does_not_predict_response` | (exact) |
| `antivegf_effect_on_margin_exposure` | `decreases_via_barrier_restoration` | (exact) |
| `antivegf_effect_on_margin_cytotoxic_sensitivity` | `decreases` | (exact) |
| `relapse_pattern_after_antivegf` | `non_enhancing_and_distant` | (exact) |

## Task（交给模型的原文）

A patient has recurrent glioblastoma. Three cytotoxic agents are available. Their blood-brain barrier (BBB) permeability labels are recorded in the accompanying dataset (`1) Supplementary-data/TableS1-Drug-for-Train-and-Validation.xlsx`, column `BBB +/-`):

  * lomustine — dataset label `Y`; nitrosourea, alkylates and crosslinks DNA
  * methotrexate — dataset label `N`; antifolate, blocks thymidylate synthesis, S-phase specific
  * doxorubicin — dataset label `N`; topoisomerase II poison, converts enzyme-DNA complexes into double-strand breaks during replication

The following facts about this tumour hold:

  1. The contrast-enhancing tumour bulk has a disrupted BBB with increased vascular permeability.
  2. The non-enhancing infiltrating margin retains largely intact BBB properties.
  3. Recurrence arises from the infiltrating cells at the margin, not from the bulk.
  4. The margin contains a hypoxic niche.

Bevacizumab (anti-VEGF) is under consideration in addition to cytotoxic therapy. Answer for sustained anti-VEGF therapy, beyond any transient vascular-normalisation window.

Write your answer as JSON to `eval_answer.json` with exactly these fields:

  * `primary_determinant_of_margin_response` (string) — one of `bbb_permeability_of_the_agent`, `proliferative_state_of_the_target_cells`, `interstitial_transport_beyond_the_vessel_wall`, `agent_specific_dna_repair_capacity`, `immunosuppressive_microenvironment`. The factor that primarily determines whether these agents reduce the recurrence-driving cell population.
  * `agent_selection_by_bbb_label` (string) — one of `lomustine_preferred_and_this_predicts_response`, `lomustine_preferred_but_this_does_not_predict_response`, `no_label_based_preference_exists`, `methotrexate_preferred_and_this_predicts_response`. How the dataset `BBB +/-` labels bear on the choice among the three agents for this patient.
  * `antivegf_effect_on_margin_exposure` (string) — one of `increases_via_vascular_normalisation`, `increases_via_reduced_interstitial_pressure`, `decreases_via_barrier_restoration`, `unchanged_because_margin_barrier_already_intact`. The effect of sustained anti-VEGF therapy on cytotoxic drug exposure of the infiltrating margin.
  * `antivegf_effect_on_margin_cytotoxic_sensitivity` (string) — one of `increases`, `decreases`, `no_effect`. The effect of sustained anti-VEGF therapy on the intrinsic sensitivity of the margin cell population to these three agents, setting aside any change in the amount of drug delivered.
  * `relapse_pattern_after_antivegf` (string) — one of `enhancing_and_local`, `non_enhancing_and_distant`, `unchanged_from_baseline_pattern`. The characteristic radiographic pattern of relapse on sustained anti-VEGF therapy.

All fields are categorical. Report exactly one option per field.

## notes — platform

Reasoning task over drug-level tabular labels plus stated clinical facts. No computation, no randomness, no reference genome. Shares its dataset with bbb_label_implication_v1, which is its deliberate counterpart: that task requires the BBB labels to be taken seriously as data, this one requires recognising that they do not govern the clinical outcome.

## notes — grader

MultipleChoice across five closed option sets. No numeric field, therefore no tolerance to calibrate - the documented failure mode of numeric graders in this benchmark family is avoided by construction. Every field is independently discriminating; none is a restatement of another. Two fields (`antivegf_effect_on_margin_exposure` and `antivegf_effect_on_margin_cytotoxic_sensitivity`) have a correct answer that contradicts the dominant published frame, so they cannot be reached by retrieval of consensus language.

## notes — what_this_tests

Whether the agent can recognise that the most salient variable in the data is not the operative one. The dataset is a BBB permeability dataset and the prompt leads with BBB labels, so the pull is to treat BBB crossing as the decision variable. In glioblastoma it largely is not: the contrast-enhancing bulk has already disrupted its own barrier, so a drug's intact-BBB permeability does not determine exposure there, and at the margin, where the barrier is intact, exposure still does not determine kill because the target cells are not cycling. An agent that reasons from the labels alone produces a confident, data-grounded, wrong recommendation - it will prefer lomustine and expect that preference to translate into response.

The second half tests whether the agent holds the direction of the anti-VEGF effect in glioblastoma rather than the direction asserted by the vascular-normalisation literature. These two bodies of work disagree, and the benchmark-relevant question is which one the agent applies to this disease.

## notes — key_metrics

The decisive pairing is `primary_determinant_of_margin_response` against `agent_selection_by_bbb_label`. An agent can read the table correctly - lomustine is the only `Y` - and still be wrong, because the question is whether that reading predicts response. The correct answer concedes the label-based preference and denies its predictive value, which is the one combination that requires both halves of the picture. `no_label_based_preference_exists` is wrong on the data: the labels do differ and lomustine is the one that reaches an intact barrier.

The mechanism behind the denial is that all three agents are replication-contingent. Methotrexate is S-phase specific. Doxorubicin requires replication to convert topoisomerase II-DNA complexes into lethal double-strand breaks. The nitrosourea's crosslinks are likewise converted into lethal lesions principally during replication. A cell that is not cycling therefore does not transduce any of these drug-target interactions into death, however much drug arrives.

The task states the hypoxic niche but does not state that its cells are not cycling, and no field name contains the word arrest or quiescence. Both inference steps are left to the agent: that hypoxia suppresses cycling, and that all three agents are replication-contingent. Every premise needed is given.

## notes — challenges

The four traps are independent of one another.

First, the BBB decoy. The accompanying dataset is a BBB permeability dataset, every agent is introduced by its BBB label, and the first two stated facts are about barrier state. An agent that answers `bbb_permeability_of_the_agent` is not pattern-matching carelessly - it is reasoning from the most prominent feature of the data it was given. The option list for `primary_determinant_of_margin_response` offers no elimination path: all five options are agent-agnostic properties of the tumour or the delivery problem rather than properties of one drug, so the shared-property shortcut that solved an earlier version of this task does not apply. `interstitial_transport_beyond_the_vessel_wall` and `immunosuppressive_microenvironment` are both real and heavily cited; `agent_specific_dna_repair_capacity` is the MGMT answer in general form.

Second, the direction of the anti-VEGF effect on exposure. Three of the four options are the vascular-normalisation frame's own vocabulary, and `increases_via_reduced_interstitial_pressure` is a genuine, extensively published mechanism. In glioblastoma, sustained anti-VEGF therapy restores barrier integrity, which reduces drug access to the non-enhancing margin - the opposite direction. Reduced enhancement is the radiographic signature of that restoration, which is why it reads as improvement.

Third, the direction of the anti-VEGF effect on intrinsic sensitivity. The published rationale for combining anti-VEGF therapy with cytotoxics and radiation is that normalisation improves oxygenation and thereby sensitivity. Under sustained therapy the vessel pruning deepens margin hypoxia instead, which suppresses cycling further and so reduces sensitivity to three replication-contingent agents. The field is phrased to set aside the delivery change, so an agent cannot answer it by reusing its exposure answer; it has to carry the hypoxia-to-cycling-to-replication chain.

Fourth, `relapse_pattern_after_antivegf` is a coherence check. An agent that answered the two direction fields in the normalisation direction and then answers `non_enhancing_and_distant` here has stated a pattern its own mechanism cannot produce. Scoring the fields jointly therefore separates an agent that holds a mechanism from one that retrieves facts.

## notes — ground_truth_derivation

Each field traces to a documented statement rather than to a calculation. Dataset labels for the three agents were read directly from TableS1 (lomustine Y, methotrexate N, doxorubicin N; verified 2026-11-05). Barrier disruption in the contrast-enhancing bulk versus relative integrity at the non-enhancing infiltrating margin, and the margin as the origin of recurrence, are drawn from the glioblastoma drug-delivery literature. Replication dependence of alkylator cytotoxicity via the O6-methylguanine mismatch-repair futile cycle, and the requirement that lesions pass through a replication cycle to become lethal, are established mechanistically. Hypoxic maintenance of a non-cycling stem-like fraction and its association with therapeutic resistance and recurrence are likewise documented. Anti-VEGF restoration of barrier integrity with reduced enhancement, the characteristically non-enhancing and distant relapse pattern that follows, and vessel pruning with deepened tumour hypoxia under sustained therapy, are the reported glioblastoma observations. See docs/P48_research_plan_and_review_1105.md for the per-claim sourcing and for the claims deliberately not made.

## notes — known_gaps

Three gaps, none of which changes an answer, all of which bound the claims.

A transient vascular-normalisation window during which oxygenation and delivery do improve is documented. The task is therefore explicitly scoped to sustained therapy beyond that window, in both the preamble and every anti-VEGF field name's description. Without that scoping, `increases_via_reduced_interstitial_pressure` would be defensible and the item would be ambiguous rather than hard.

Deepened hypoxia under sustained anti-VEGF therapy rests more heavily on preclinical models than on direct human measurement at the infiltrating margin; the clinical inference is from the relapse phenotype rather than from oxygen electrodes in margin tissue. The answer `decreases` is the defensible direction, but it is an inference, not a measurement.

No retrieved source directly measured sparing of non-cycling glioma cells by these three specific agents. That step is an inference from the established replication requirement, which is why the reasoning is framed in terms of replication dependence rather than an asserted sparing fraction. One study in Xenopus extracts proposes a replication-independent component of alkylator action operating alongside the cell-cycle-dependent mode, which would attenuate but not reverse the conclusion.

Temozolomide is deliberately absent: it is the standard-of-care agent and the clearest example of replication-dependent cytotoxicity, but it is not present in the dataset, so including it would have meant asserting a BBB label the data does not carry.

## notes — version_history

v1 and v2 are archived under superseded_for_comparison/evals/. v1 stated in premise 4 that the margin cells were growth-arrested and named a field `arrest_effect_on_...`, which left the correct option recoverable by lexical match with no biology required. v2 removed that leak but remained solvable without glioma knowledge: four of its six fields were given in the prompt or were textbook, and its decisive field could be reached by elimination because the four distractors were all agent-specific while the correct option was the only agent-agnostic one - an elimination path that the instruction to choose the factor applying to all three agents actively signposted. v3 removes that instruction, makes every option in that field agent-agnostic, drops the two fields that were answerable by reading the prompt, merges the two redundant hypoxia fields, and replaces them with two direction-reversal fields whose correct answers contradict the dominant published frame.

