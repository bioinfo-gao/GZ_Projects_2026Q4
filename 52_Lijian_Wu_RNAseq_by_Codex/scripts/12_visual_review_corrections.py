#!/usr/bin/env python3
"""Archive first review package; remove defective optional static exports from client scope."""
from pathlib import Path
import shutil,json
P=Path(__file__).resolve().parents[1];D=P/'Data_Analysis_20261109';O=P/'custom_research_report_20261109';A=P/'superseded_for_comparison/delivery/v1_20261109_before_visual_review'
A.mkdir(parents=True,exist_ok=False)
# Archive original first package in full and the scripts that produced it.
shutil.move(str(O),str(A/O.name));(A/'scripts_snapshot').mkdir()
for name in ['6_generate_research_report.py','8_pre_delivery_audit.py','9_build_delivery.py','11_resume_delivery.sh']:
 shutil.copy2(P/'scripts'/name,A/'scripts_snapshot'/name)
report=D/'Crtc2_Line_Treatment_RNAseq_Report_1109.md';shutil.copy2(report,A/'source_report_before_corrections.md')
static=D/'QC/multiqc/star_salmon/multiqc_report_plots';shutil.move(str(static),str(A/'source_multiqc_report_plots'))
# Changes are language/output-scope only, no numbers/contrasts/models changed.
pairs=[('under untreated (un), TNF, PGE2 and adenosine conditions','under untreated (un), tumor necrosis factor (TNF), prostaglandin E2 (PGE2) and adenosine conditions'),('## 3. Key Findings\n\n','## 3. Key Findings\n\nDifferentially expressed genes (DEGs) are summarised below; log2FC denotes log2 fold change.\n\n'),('Every comparison supplies a full tested-gene table, significant subset, volcano plot, and top-50 DEG heatmap.','Every comparison supplies a full tested-gene table, significant subset, full-range volcano plot, and top-50 DEG heatmap. Two adenosine comparisons also include a volcano_zoom.pdf companion with |log2FC| ≤ 3; the caption states the number of points outside the view, and the full tables and original plots retain those genes.'),('Pipeline MultiQC HTML report and supporting data','Self-contained MultiQC HTML report and supporting data; use the HTML for interactive QC'),('Full and significant tables, volcano and top-50 heatmap for each of six comparisons','Full and significant tables, full-range volcano and top-50 heatmap for each of six comparisons; two adenosine comparisons also have a volcano_zoom.pdf companion')]
text=report.read_text()
for a,b in pairs:
 assert a in text,a;text=text.replace(a,b)
report.write_text(text)
script=P/'scripts/6_generate_research_report.py';t=script.read_text()
# Keep the visible generator aligned with the report corrections; report-only table intro is inserted explicitly.
t=t.replace('under untreated (un), TNF, PGE2 and adenosine conditions','under untreated (un), tumor necrosis factor (TNF), prostaglandin E2 (PGE2) and adenosine conditions')
t=t.replace("add('## 3. Key Findings');add()","add('## 3. Key Findings');add();add('Differentially expressed genes (DEGs) are summarised below; log2FC denotes log2 fold change.');add()")
for a,b in pairs[2:]:t=t.replace(a,b)
script.write_text(t)
audit=P/'scripts/8_pre_delivery_audit.py';t=audit.read_text().replace("len(pdfs)==174,'174 expected result PDFs'","len(pdfs)==106,'106 expected result PDFs (104 original analytical plots + two zoom companions)'").replace('all 174 PDFs decode successfully','all 106 PDFs decode successfully')
t=t.replace("pdfs=list(D.rglob('*.pdf'));", "check(not (D/'QC/multiqc/star_salmon/multiqc_report_plots').exists(),'client QC excludes defective optional static exports')\ncheck((D/'QC/multiqc/star_salmon/multiqc_report.html').stat().st_size>1000000,'self-contained MultiQC HTML present')\npdfs=list(D.rglob('*.pdf'));")
audit.write_text(t)
pack=P/'scripts/9_build_delivery.py';t=pack.read_text().replace("P/'logs/audit_source_20261109.log'","P/'logs/audit_source_final_20261109.log'").replace("P/'logs/delivery_sha256_20261109.json'","P/'logs/delivery_sha256_final_20261109.json'").replace("P/'logs/audit_delivery_20261109.log'","P/'logs/audit_delivery_final_20261109.log'")
pack.write_text(t)
(A/'README.md').write_text('''# P52 first internal package, superseded before final delivery audit\n\nThis is the original first-built package, not a reconstruction. No external sending or user final approval occurred. The automated content audit passed, but visual inspection found the optional 70 MultiQC static PDF exports had repeated/low-contrast content. All optional static plot exports are preserved here; the corrected client package retains the self-contained MultiQC HTML and supporting tables.\n\nThe report was also corrected to expand TNF, PGE2, DEG and log2FC at first prose/table use. Two new adenosine volcano zoom companions improve the view of the central distribution while preserving every original full-range volcano, table and DEG count. The producing scripts are in scripts_snapshot/. Compare the report here with the project-root current delivery report. The current valid package is the project-root custom_research_report_20261109/.\n''')
print('First package and original code archived:',A)
