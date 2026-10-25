"""Stage 12b (P45; from P44 12b_delivery_tables.py, request rows updated for the complete package) — delivery-level tables, written into the delivery folder after header sanitisation:
  tables/client_request_coverage.tsv  every question and deliverable in the client's 2026-10-13 request -> status ->
                                      report section -> delivered file(s) -> boundary   (numbers stay in the report)
  tables/external_input_status.tsv    inputs not available to this analysis and their effect on interpretation
  tables/data_inventory.tsv           every delivered file: relative path, folder category, bytes
Every file named in client_request_coverage.tsv must exist (exit 1 otherwise).
Usage: 12b_delivery_tables.py <delivery_dir>
"""
import glob
import os
import sys

D = sys.argv[1]
T = "tables/junction_table_target.tsv"
COV = [
    ("Analysis 1", "Base-pair resolution breakpoints, including any microhomology or inserted sequence at the junction",
     "answered", "7.4, 7.5", f"{T};figures/junction_views/;tables/key_junction_read_evidence.tsv;sequences/junction_and_reference_templates.fa",
     "Breakpoints are given as the interval within which the junction can slide over identical sequence; "
     "events with one supporting fragment are listed but not interpreted"),
    ("Analysis 1", "Junction-spanning read counts per event, before and after read rescue",
     "answered", "5, 7.2, 7.5", f"{T};tables/junction_by_sample_long.tsv;tables/rescue_source_summary.tsv",
     "Counts are independent DNA fragments; read-end counts are in junction_by_sample_long.tsv"),
    ("Analysis 1", "Which junctions are well supported versus ambiguous because of repeat content",
     "answered", "5, 7.1, 7.4, 7.5", f"{T};tables/key_junction_read_evidence.tsv;tables/structural_coherence_tests.tsv",
     "Evidence class and claim level use pre-specified criteria (Section 5); re-alignment of the same data is not an independent test"),
    ("Analysis 1", "Feasibility and what read rescue adds beyond the current calls",
     "answered", "3, 7.2, 7.4, 7.5", f"{T};tables/rescue_source_summary.tsv;tables/evidence_filter_impact.tsv", ""),
    ("Analysis 1", "Which Section 7.3 events are called by both Manta and DELLY versus only one",
     "answered", "7.5", f"{T};tables/caller_concordance_by_sample.tsv",
     "Two further callers (GRIDSS, SvABA) are reported alongside; agreement between callers run on the same data is not independent validation"),
    ("Analysis 2", "Events present only in the +Dox clones (S04, S05) and absent from S01-S03: at the locus",
     "answered", "5, 7.8", T, ""),
    ("Analysis 2", "Events present only in the +Dox clones (S04, S05) and absent from S01-S03: genome-wide",
     "answered", "5, 7.9", "genome_wide_DLD1_comparison/dld1_genome_wide_events.tsv;genome_wide_DLD1_comparison/dld1_dox_specific_events.tsv",
     "Genome-wide calls rest on DELLY re-genotyping and Manta per-sample calls; they were not read-rescued or curated like the locus"),
    ("Analysis 2", "Filtering criteria: tier, supporting-read count and allele-fraction thresholds",
     "answered", "5", f"{T};genome_wide_DLD1_comparison/dld1_genome_wide_events.tsv",
     "No allele-fraction minimum is applied; allele fraction is reported per sample"),
    ("Supporting analysis", "Copy-number step position and chr11 B-allele fraction",
     "added", "7.6, 7.7", "supporting_analyses/depth_boundary/;supporting_analyses/chr11_BAF/", ""),
    ("Supporting analysis", "Junction-spanning PCR primer designs for confirmation",
     "added", "7.11, 8.4", "tables/primer_design.tsv;tables/primer_binding_sites.tsv", "In-silico design only; no PCR was run"),
    ("Deliverable", "SV VCFs from both Manta and DELLY, unfiltered and filtered",
     "delivered", "9", "genome_wide_SV_VCF/", "As produced on 2026-10-09; local paths removed from headers"),
    ("Deliverable", "SNV/indel VCFs",
     "delivered", "7.10, 9", "genome_wide_SNV_indel_VCF/", "Joint genome-wide calls for S01-S10, chromosomes 1-22, X, Y and M"),
    ("Deliverable", "Annotated per-event summary table (TSV) with a column indicating which caller(s) support each event",
     "delivered", "7.5, 7.9, 9",
     f"{T};tables/junction_table_window_context.tsv;genome_wide_DLD1_comparison/dld1_genome_wide_events.tsv",
     "Locus tables cover the analysis window (target +- 50 kb, four callers); the genome-wide table covers DLD1 S01-S05 (DELLY and Manta)"),
]
EXT = [
    ("Junction-spanning PCR results", "not available", "The +Dox junctions rest on sequencing data only; "
     "candidate primers are in tables/primer_design.tsv"),
    ("Fluorescence in situ hybridisation images or scores", "not available",
     "The telomere-capped end in S04 and the fold-back in S05 are not cytogenetically confirmed"),
    ("Additional independent +Dox clones", "not available", "Recurrence of the homolog-specific loss rests on n = 2"),
    ("Independent assay for any junction", "not available",
     "No junction can reach the claim level 'confirmed'; the highest level reached is a limited candidate coherent with the copy-number step"),
]

missing = []
for row in COV:
    for p in row[4].split(";"):
        if p and not os.path.exists(f"{D}/{p}"):
            missing.append(p)
if missing:
    sys.exit(f"coverage table names missing files: {sorted(set(missing))}")
with open(f"{D}/tables/client_request_coverage.tsv", "w") as fh:
    fh.write("request_part\tclient_request\tstatus\treport_section\tdelivered_files\tboundary\n")
    for r in COV:
        fh.write("\t".join(r) + "\n")
with open(f"{D}/tables/external_input_status.tsv", "w") as fh:
    fh.write("input\tstatus\teffect_on_interpretation\n")
    for r in EXT:
        fh.write("\t".join(r) + "\n")
files = sorted(f for f in glob.glob(f"{D}/**/*", recursive=True) if os.path.isfile(f))
with open(f"{D}/tables/data_inventory.tsv", "w") as fh:
    fh.write("relative_path\tfolder\tbytes\n")
    for f in files:
        rel = os.path.relpath(f, D)
        fh.write(f"{rel}\t{os.path.dirname(rel) or '.'}\t{os.path.getsize(f)}\n")
    fh.write(f"tables/data_inventory.tsv\ttables\tself\n")
print(f"12b OK: coverage {len(COV)} rows, external {len(EXT)} rows, inventory {len(files) + 1} files")
