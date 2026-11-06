# glioma_margin_rate_limit_v1

> ⚠ **本文是渲染产物，不是真身。** 内容真身在 `evals/glioma_margin_rate_limit_v1.json`；要改内容请改那个 JSON，然后重新渲染本文。
> 渲染时间：2026-11-05（渲染器：`scripts/2_render_eval.py`）

**Data:** https://github.com/bioinformatics-gao/CASE-BBB-prediction-Data

## Grader — MultipleChoice

| 字段 | 标准答案 | 容差 |
| :--- | :---: | :---: |
| `agent_reaching_margin_unaided` | `lomustine` | (exact) |
| `enhancement_as_exposure_proxy` | `invalid` | (exact) |
| `delivery_alone_sufficient` | `no` | (exact) |
| `remaining_primary_limit` | `target_cell_quiescence` | (exact) |
| `arrest_effect_on_cytotoxic_sensitivity` | `decreases` | (exact) |
| `hypoxia_effect_on_radiotherapy` | `reduces_efficacy` | (exact) |

## Task（交给模型的原文）

A patient has recurrent glioblastoma. You are evaluating three cytotoxic agents whose blood-brain barrier (BBB) permeability labels are recorded in the accompanying dataset (`1) Supplementary-data/TableS1-Drug-for-Train-and-Validation.xlsx`, column `BBB +/-`):

  * lomustine — dataset label `Y`; nitrosourea, alkylates and crosslinks DNA
  * methotrexate — dataset label `N`; antifolate, blocks thymidylate synthesis, S-phase specific
  * doxorubicin — dataset label `N`; topoisomerase II poison, converts enzyme-DNA complexes into double-strand breaks during replication

The following clinical and histological facts about this tumour hold:

  1. The contrast-enhancing tumour core has a disrupted BBB with increased vascular permeability.
  2. The non-enhancing infiltrating margin retains largely intact BBB properties.
  3. Recurrence arises from the infiltrating cells at the margin, not from the core.
  4. The margin contains a hypoxic niche in which a substantial fraction of tumour cells are growth-arrested.

Answer the following. Then assume a hypothetical delivery technology that completely abolishes the BBB at the infiltrating margin, so that all three agents reach margin tissue at the same concentration they reach the core. Re-evaluate.

Write your answer as JSON to `eval_answer.json` with exactly these fields:

  * `agent_reaching_margin_unaided` (string) — one of `lomustine`, `methotrexate`, `doxorubicin`, `none`. Which of the three is expected to reach the infiltrating margin without the delivery technology, based on the dataset labels.
  * `enhancement_as_exposure_proxy` (string) — one of `valid`, `invalid`. Whether the extent of contrast enhancement on MRI is a valid proxy for effective drug exposure of the recurrence-driving cell population.
  * `delivery_alone_sufficient` (string) — one of `yes`, `no`. With the delivery technology in place and penetration no longer limiting, whether these agents are now expected to eradicate the recurrence-driving population.
  * `remaining_primary_limit` (string) — one of `target_cell_quiescence`, `mgmt_mediated_repair`, `drug_efflux_bcrp`, `insufficient_dose`, `mismatch_repair_deficiency`. The factor that remains primarily limiting for all three agents once penetration is solved. Choose the factor that applies to all three agents by virtue of a property they share.
  * `arrest_effect_on_cytotoxic_sensitivity` (string) — one of `increases`, `decreases`, `no_effect`. The effect of hypoxia-induced growth arrest on the sensitivity of those cells to these three agents.
  * `hypoxia_effect_on_radiotherapy` (string) — one of `increases_efficacy`, `reduces_efficacy`, `no_effect`. The effect of the same hypoxic state on the efficacy of concurrent radiotherapy.

All fields are categorical. Report exactly one option per field.

## notes — platform

Reasoning task over drug-level tabular labels plus stated clinical facts. No computation, no randomness, no reference genome. Shares its dataset with bbb_label_implication_v1.

## notes — grader

MultipleChoice across six closed option sets. No numeric field, therefore no tolerance to calibrate - the documented failure mode of numeric graders in this benchmark family is avoided by construction. Partial credit is meaningful here: `delivery_alone_sufficient` and `remaining_primary_limit` are the two fields that separate expert from non-expert reasoning, and can be weighted if the harness supports it.

## notes — what_this_tests

Whether the agent recognises that BBB penetration is not the rate-limiting step for curing this tumour, and that solving delivery does not solve the problem. Two independent barriers protect the recurrence-driving population, and an agent that holds only one of them produces a confidently wrong recommendation. Hold only the penetration barrier and you conclude that a delivery technology is curative. Hold only the cell-state barrier and you conclude that penetration is irrelevant. The correct answer requires both.

## notes — key_metrics

The decisive field is `remaining_primary_limit`. All three agents share one property: their cytotoxicity is contingent on DNA replication. Methotrexate is S-phase specific. Doxorubicin requires replication to convert topoisomerase II-DNA complexes into lethal double-strand breaks. The nitrosourea's crosslinks are likewise converted to lethal lesions principally during replication. A growth-arrested cell therefore does not transduce any of these drug-target interactions into death, however much drug arrives. The clinical corollary is the paradox worth seeing: the hypoxic stress that arrests these cells, and that would naively read as the tumour being under metabolic duress, is the same state that renders them refractory.

## notes — challenges

The strongest distractor is `mgmt_mediated_repair`. MGMT promoter methylation is the single most-cited resistance mechanism in glioblastoma, so an agent pattern-matching on 'glioblastoma drug resistance' retrieves it first. It is wrong here for two reasons: it is specific to O6-alkylguanine lesions and so does not explain failure of methotrexate or doxorubicin, and the question asks for the factor that applies to all three by virtue of a shared property. `mismatch_repair_deficiency` is a second trap of the same kind - it is a genuine temozolomide resistance mechanism and genuinely replication-linked, but again agent-specific rather than shared. `drug_efflux_bcrp` is plausible and real, but efflux inhibition is documented to improve penetration in the disrupted core rather than at the intact margin, so it does not explain the residual failure after delivery is solved.

`enhancement_as_exposure_proxy` carries an independent trap. Contrast enhancement maps the disrupted core, which is where the cells that do not drive recurrence live; the barrier in glioblastoma is a mosaic of states and radiographic leak is a poor proxy for effective exposure of the infiltrating population. The clinical demonstration is anti-VEGF therapy, which reduces enhancement by restoring barrier integrity - radiographically an improvement, pharmacologically a new obstacle - after which relapse is characteristically non-enhancing and distant from the original core.

`hypoxia_effect_on_radiotherapy` checks that the agent does not treat hypoxia as a one-sided liability for the tumour. The same niche reduces radiotherapy efficacy through the oxygen effect, since oxygen is required to fix radiation-induced DNA radicals into permanent damage. Reported glioma oxygen enhancement ratios are approximately 1.2 to 1.6, and roughly threefold dose is required in the absence of oxygen at the limit. So the arrested population is concurrently less sensitive to the cytotoxics and less sensitive to the radiation given alongside them.

## notes — ground_truth_derivation

Each field traces to a documented statement rather than to a calculation. Dataset labels for the three agents were read directly from TableS1 (lomustine Y, methotrexate N, doxorubicin N; verified 2026-11-05). Core-versus-margin barrier heterogeneity, the margin as the origin of recurrence, and the inadequacy of enhancement as an exposure proxy are drawn from the glioblastoma drug-delivery literature. Replication dependence of temozolomide-class alkylator cytotoxicity via the O6-methylguanine mismatch-repair futile cycle, and the requirement that lesions pass through a replication cycle to become lethal, are established mechanistically. Hypoxic maintenance of a quiescent stem-like fraction and its association with therapeutic resistance and recurrence, and the oxygen enhancement ratio, are likewise documented. See docs/P48_research_plan_and_review_1105.md for the per-claim sourcing and for the claims deliberately not made.

## notes — known_gaps

No retrieved source directly measured sparing of non-dividing glioma cells by these specific agents; that step is an inference from the established replication requirement, and the task is therefore phrased in terms of replication dependence rather than asserting a measured sparing fraction. One study in Xenopus extracts proposes a replication-independent component of alkylator action operating in addition to the cell-cycle-dependent mode, which would attenuate but not reverse the conclusion. Temozolomide itself is deliberately absent from the task: it is the standard-of-care agent and the clearest example of replication-dependent cytotoxicity, but it is not present in the dataset, so including it would have meant asserting a BBB label the data does not carry.

