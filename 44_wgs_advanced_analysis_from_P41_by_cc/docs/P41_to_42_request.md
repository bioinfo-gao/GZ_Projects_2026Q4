| **from:** | Lin, Xin <**Xin.Lin@childrens.harvard.edu**>                                                                         |
| --------------- | -------------------------------------------------------------------------------------------------------------------------- |
| **to:**   | 1. "Singh, Vedita" <Vedita.Singh@childrens.harvard.edu>,1. Zhen Gao <zorin.gao@gmail.com>                                  |
| **cc:**   | 1. Jianjun Sun <jianjun@athenomics.com>,1. "Li, Julia" <Julia.Li@childrens.harvard.edu>,1. Jing Shi <j.shi@athenomics.com> |
| **date:** | **Oct 13, 2026, 1:16 PM**                                                                                           |

Hi Zhen,
Thank you for the preliminary analysis. The report was thorough, and we took a few days to review and discuss it as a group.

We would like to request additional work focused on the events described in Section 7.3, specifically at the locus we are interested in - chr11:114,604,000-114,625,500. Two analyses:

1. Junction-level analysis at the locus, with read rescue

You mentioned that this locus is highly repetitive and that resolving junctions would require manually retrieving reads the standard pipeline discards. We would like to go ahead with that.

Our understanding is that this means re-extracting all reads at the locus — including MAPQ0 and multi-mapping reads, supplementary and secondary alignments, soft-clipped reads, and unmapped reads whose mates map to the region — followed by local reassembly to resolve breakpoints. Please correct us if you would approach it differently.

We are specifically interested in:

base-pair resolution breakpoints, including any microhomology or inserted sequence at the junction
junction-spanning read counts per event, before and after read rescue
your assessment of which junctions are well supported versus ambiguous because of repeat content
Before you begin, could you confirm this is feasible and give us a sense of what it would add beyond the current calls? Since both Manta and DELLY were run, could you also indicate which of the Section 7.3 events are called by both callers versus only one? Given the low read support on many of these events, caller concordance would help us judge which are worth following up.

2. Clone-specific event calling

Please compare the DLD1 +Dox clones (S04 DLD1_c3_dox, S05 DLD1_c4_dox) against the DLD1 parental population (S01_DLD1_pop) and the control clones (S02_DLD1_c1_ctrl, S03_DLD1_c2_ctrl) to identify events present only in the +Dox clones. Please also state the filtering criteria you apply — tier, supporting-read count, and allele-fraction thresholds — so we can interpret the output correctly.

To be clear, we understand the current run is already aligned to hg38. We are not requesting re-alignment, only these additional analyses on the existing alignments.

Deliverables

Once the analysis is complete, we would like the genome-wide files:

SV VCFs from both Manta and DELLY, unfiltered and filtered
SNV/indel VCFs
the annotated per-event summary table (TSV), regenerated for the new comparisons, with a column indicating which caller(s) support each event
Please let us know your preferred transfer method and the approximate total size.

Budget and timeline

We have approved up to $840 for the “Custom Data Analysis” service. Please treat this as a hard cap — if the scope above would exceed it, pause and let us know before proceeding so we can adjust priorities or seek additional approval.

Could you also give us an estimated turnaround time?

Please let us know if any part of this request needs clarification.

Best regards,
Xin
