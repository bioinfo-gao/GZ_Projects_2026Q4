#!/usr/bin/env python3
"""P52: write first report from existing P52 results; no refit, no other-project inputs."""
import csv,json,os,math
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[1]; D=P/'Data_Analysis_20261109'
def rows(path,sep=','):
 with Path(path).open() as f:return list(csv.DictReader(f,delimiter=sep))
facts=json.loads((P/'logs/report_facts_20261109.json').read_text()); fs=facts['filters']; pc=facts['pca'];qc=facts['qc']
st=rows(P/'scripts/0_sample_table.csv');de=rows(D/'DE_PCA_Results/DEG_summary_all_comparisons.csv');en=rows(D/'Enrichment/Enrichment_summary_all_comparisons.csv');metrics=rows(D/'Reads/library_metrics_per_sample.csv'); met={r['sample_id']:r for r in metrics}
assert len(st)==24 and len(de)==6 and len(en)==90 and all(r['status']=='ok' for r in en)
keys=[r['Comparison'] for r in de]; allgenes={k:rows(D/'DE_PCA_Results'/k/(k+'_ALL_genes.csv')) for k in keys};targets={k:next(r for r in allgenes[k] if r['gene_name']=='Crtc2') for k in keys}
hall={k:rows(D/'Enrichment'/k/'GSEA/GSEA_MSigDB_Hallmark.csv') for k in keys}
reportdate=os.environ.get('DELIVERY_DATE',datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d'))
filename='Crtc2_Line_Treatment_RNAseq_Report_'+reportdate[5:7]+reportdate[8:10]+'.md';report=D/filename
assert not report.exists(),f'Refuse overwrite: {report}'
n=lambda x:f'{int(float(x)):,}'
name=lambda k:k.replace('_vs_',' vs ')
line=[]
def add(x=''):line.append(x)
def table(headers,data):
 add('| '+' | '.join(headers)+' |');add('| :--- |'+' :---: |'*(len(headers)-1))
 for r in data:add('| '+' | '.join(str(x).replace('|','/') for x in r)+' |')
 add()
lfcs=[float(targets[k]['log2FoldChange']) for k in keys[:4]]
assert all(float(targets[k]['padj'])<=.05 and float(targets[k]['log2FoldChange'])<-.263 for k in keys[:4])
cr=rows(P/'scripts/rds_cache/crtc2_expression_internal.csv');means={}
for group in sorted({r['group'] for r in st}):
 ids={r['sample_id'] for r in st if r['group']==group};means[group]=sum(float(r['Crtc2_norm_count']) for r in cr if r['sample_id'] in ids)/len(ids)
ratios=[means['Crtc2-'+t]/means['control-'+t] for t in ['un','TNF','PGE2','adenosine']]
raw=[int(r['raw_read_pairs']) for r in st]
bad=[r for r in pc['nearest'] if not r['same_group']]
# Predeclared conclusions are guarded by directly read gene-set results, not old reports.
H=['TNFA_SIGNALING_VIA_NFKB','INFLAMMATORY_RESPONSE','INTERFERON_GAMMA_RESPONSE','E2F_TARGETS','G2M_CHECKPOINT']
def hrow(k,s):return next((r for r in hall[k] if r['ID']=='HALLMARK_'+s and float(r['p.adjust'])<.05),None)
assert all(hrow(k,'INFLAMMATORY_RESPONSE') and float(hrow(k,'INFLAMMATORY_RESPONSE')['NES'])<0 for k in keys[:4])
assert all(hrow(k,'E2F_TARGETS') and float(hrow(k,'E2F_TARGETS')['NES'])<0 for k in keys[4:])
claims=[
 ('Quality and grouping',f'All 24 libraries are retained. Unique alignment is {qc["unique"][0]:.1f}–{qc["unique"][1]:.1f}%; the mean between-group / within-group distance in the first two principal components is {pc["separation_ratio"]:.1f}. These checks support the expression analysis; they do not independently authenticate sample identity.','7.1','QC/multiqc/; DE_PCA_Results/PCA_coordinates_overview.csv','None; no identity assay was performed.'),
 ('Crtc2 transcript',f'Crtc2 messenger RNA (mRNA) is lower in the Crtc2 line in all four conditions (log2 fold changes {", ".join(f"{v:+.2f}" for v in lfcs)}). The transcript remains detectable; protein loss and the nature of the modification are unresolved.','7.2','DE_PCA_Results/<comparison>/<comparison>_ALL_genes.csv; row Crtc2','Protein and functional effects require independent confirmation.'),
 ('Untreated line difference',f'The untreated lines differ in {n(de[0]["Total"])} genes. Lower inflammatory/interferon signatures and higher cell-cycle signatures are associated with the Crtc2 line; these are expression associations rather than a demonstrated mechanism.','7.2 and 7.3','DE_PCA_Results/Crtc2-un_vs_control-un/; Enrichment/Crtc2-un_vs_control-un/GSEA/','An association with line identity; not proof of direct CRTC2 regulation.'),
 ('TNF condition',f'The largest line difference is under tumor necrosis factor (TNF): {n(de[1]["Total"])} genes. Inflammatory and interferon signatures are lower in the Crtc2 line under TNF. This comparison does not establish a difference in treatment responsiveness.','7.2 and 7.3','DE_PCA_Results/Crtc2-TNF_vs_control-TNF/; Enrichment/Crtc2-TNF_vs_control-TNF/GSEA/','A line-by-treatment interaction was not tested.'),
 ('PGE2 and adenosine',f'In control cells, prostaglandin E2 (PGE2) changes {n(de[4]["Total"])} genes and adenosine changes {n(de[5]["Total"])}. Both show inflammatory/interferon enrichment and lower cell-cycle signatures; fewer adenosine genes cross the chosen gene-level threshold.','7.2 and 7.3','DE_PCA_Results/control-PGE2_vs_control-un/; DE_PCA_Results/control-adenosine_vs_control-un/; Enrichment/','Pathway enrichment does not measure pathway activity or establish equivalent effect sizes.'),
 ('Practical next steps','Confirm CRTC2 protein abundance and selected transcript changes in independent cultures. A formal line-by-treatment interaction analysis would address whether the two lines respond differently to treatment.','8.4','Section 8.4','These are proposed follow-up tests, not completed validation.')]
add('# Bioinformatics Analysis Report: Control and Crtc2 Cell Lines under TNF, PGE2 and Adenosine');add()
for k,v in [('Report Date',reportdate),('Prepared by','Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics'),('Species','Mus musculus (mouse; GRCm39, GENCODE M35)'),('Tissue / Cell','Established mouse cell line (cell line name not specified by client)'),('Assay','Bulk RNA sequencing (RNA-seq), 24 samples'),('Analysis Platform','Linux HPC server')]:add(f'**{k}:** {v}  ')
add();add('## 1. Objectives');add()
add('Compare the client-labelled control and Crtc2 cell lines under untreated (un), tumor necrosis factor (TNF), prostaglandin E2 (PGE2) and adenosine conditions, with three biological replicates per condition. CRTC2 denotes CREB-regulated transcription coactivator 2; CREB is cyclic adenosine monophosphate response element-binding protein. The modification in the Crtc2 line was not specified. The six requested comparisons use the second named group as baseline.');add()
table(['Comparison','Numerator','Baseline'],[(name(r['Comparison']),r['Numerator'],r['Baseline']) for r in de])
add('## 2. Core Conclusions');add();add('This section explains the practical interpretation by experimental question. Section 3 summarises the measurements by topic. Section 8.2 adds the evidence files and inference boundaries for each conclusion; Section 8.3 describes limitations, and Section 8.4 identifies follow-up tests.');add()
for title,text,_,_,_ in claims:add(f'- **{title}.** {text}');add()
add('## 3. Key Findings');add();add('Differentially expressed genes (DEGs) are summarised below; log2FC denotes log2 fold change.');add()
table(['Comparison','Up','Down','Total DEGs','DEGs with absolute log2FC ≥ 1'],[(name(r['Comparison']),n(r['Up']),n(r['Down']),n(r['Total']),n(int(r['Up_lfc1'])+int(r['Down_lfc1']))) for r in de])
add('A differentially expressed gene (DEG) meets adjusted p-value (padj) ≤ 0.05 and absolute log2 fold change (log2FC) ≥ 0.263, using shrunken fold changes. The Benjamini–Hochberg (BH) procedure controls the false discovery rate (FDR) within each comparison. The 0.263 cutoff is approximately 1.2-fold; it was explicitly selected for this analysis. The stricter twofold criterion from the form is also summarised above. Positive log2FC means higher expression in the numerator.');add()
add(f'Quality control (QC) shows {min(raw)/1e6:.1f}–{max(raw)/1e6:.1f} million raw read pairs per sample, {qc["mapped"][0]:.1f}–{qc["mapped"][1]:.1f}% total alignment and {qc["rrna"][0]:.2f}–{qc["rrna"][1]:.2f}% ribosomal RNA (rRNA). Reverse library strandedness was inferred for all 24 samples.');add()
add('Gene-level and pathway-level significance answer different questions: a treatment can have relatively few threshold-crossing genes while showing coordinated changes across gene sets. The complete pathway summary is provided in Section 7.3.');add()
add('## 4. Sample Information');add()
table(['Sample ID','Client sample name','Cell line','Treatment','Replicate','Raw read pairs','Estimated gene counts (rounded sum)','Detected genes'],[(r['sample_id'],r['client_sample_name'],r['genotype'],r['treatment'],r['replicate'],n(r['raw_read_pairs']),n(met[r['sample_id']]['assigned_reads']),n(met[r['sample_id']]['genes_detected'])) for r in st])
add('Client sample names retain the form spelling PEG2; PGE2 is used for analysis labels, as written in the requested comparisons. FASTQ-to-client sample mapping was confirmed before analysis. The form specifies established mouse cells, snap-frozen. The specific cell line and Crtc2 modification were not provided. The assigned-count sum is based on rounded Salmon gene estimates and is not a direct raw-fragment count.');add()
add('## 5. Analysis Rationale and Decision Criteria');add()
table(['Step','Criterion','Reason and scope'],[
 ('Reference','GRCm39 primary assembly; GENCODE M35','Matches the stated species and the reference used for quantification.'),
 ('Gene-name filter',f'{n(fs["n_original"])} → {n(fs["n_after_regex"])} genes','Exclude ribosomal-protein names and specified mitochondrial/noncoding/predicted-name patterns. This is a name-based restriction, not a complete biotype classification; excluded loci remain in the supplied unfiltered matrices.'),
 ('Low-expression filter',f'{n(fs["n_after_regex"])} → {n(fs["n_final"])} genes in the model union','For each comparison, ≥10 rounded estimated counts in ≥3 of the 6 compared samples OR ≥8 of all 24 samples. This retains condition-restricted expression without requiring expression in most groups.'),
 ('Model','One DESeq2 ~ Group model, eight levels','Six predefined contrasts share the fitted model. Biological samples, not reads, are the replicate units; no samples were removed.'),
 ('Multiplicity','independentFiltering=FALSE; BH within the eligible set','Each comparison uses its own explicitly defined testable genes.'),
 ('Effect size','ashr shrinkage; absolute log2FC ≥ 0.263','Stabilise effect estimates and apply the confirmed effect-size threshold alongside padj.'),
 ('Over-representation analysis (ORA)','UP / DOWN / ALL; comparison-specific mapped background','Compare the selected gene set with genes that were testable and could be mapped to identifiers.'),
 ('Gene set enrichment analysis (GSEA)','Unshrunken DESeq2 Wald statistic; seed 42','Assess coordinated ranked expression shifts without restricting input to DEGs. Enrichment is an association, not an activity assay.')])
add('## 6. Methods');add()
v=facts['versions']
table(['Step','Software / reference','Parameters'],[
 ('Pipeline','nf-core/rnaseq 3.15.1','STAR alignment with Salmon quantification; paired-end reads; strandedness inferred automatically'),
 ('Alignment','STAR; GRCm39 / GENCODE M35','--twopassMode None --outFilterMultimapNmax 3'),
 ('Preprocessing and QC','Trim Galore!, FastQC, RSeQC, Qualimap, Picard, MultiQC','Pipeline-generated trimming, alignment and library QC summaries'),
 ('Differential expression','DESeq2 '+v[0]+'; ashr '+v[1],'Wald test; ~ Group; explicit per-comparison eligibility and BH correction; shrunken log2FC'),
 ('ORA / GSEA','clusterProfiler '+v[2]+'; fgsea '+v[3],'Gene Ontology (GO) biological process (BP), molecular function (MF), cellular component (CC); Kyoto Encyclopedia of Genes and Genomes (KEGG); Molecular Signatures Database (MSigDB) Hallmark'),
 ('Mouse identifier mapping','org.Mm.eg.db '+v[4],'Entrez identifier mapping; GSEA duplicate identifiers represented by the largest absolute Wald statistic')])
add('The supplied count matrix contains unrounded Salmon gene-level estimated counts. The expression filter and DESeq2 use integer-rounded estimates. Transcripts per million (TPM) are supplied separately. Gene-type and low-expression filters apply to statistical testing, not to the unfiltered delivered count/TPM matrices. ORA and GSEA tables contain terms passing their output thresholds; missing terms must not be read as demonstrated absence of a biological process.');add()
add('## 7. Results');add();add('### 7.1 Principal Component Analysis (PCA)');add()
add(f'The first principal component (PC1) explains {pc["percentVar"][0]}% of the variance and PC2 explains {pc["percentVar"][1]}%. PCA uses variance-stabilised expression of the top 500 variable model genes. The overview includes 24 points and 24 sample labels. Mean within-group distance is {pc["within_dist"]:.2f} and mean between-group distance is {pc["between_dist"]:.2f} in PC1–PC2 space (ratio {pc["separation_ratio"]:.1f}). Coordinates are in `DE_PCA_Results/PCA_coordinates_overview.csv`; the figure is `DE_PCA_Results/PCA_AllSamples_Overview.pdf`.');add()
coords=rows(D/'DE_PCA_Results/PCA_coordinates_overview.csv'); cent={}
for g in sorted({r['Group'] for r in coords}):
 z=[r for r in coords if r['Group']==g];cent[g]=[sum(float(r[k]) for r in z)/len(z) for k in ['PC1','PC2']]
table(['Group','Mean PC1','Mean PC2'],[(g,f'{xy[0]:.2f}',f'{xy[1]:.2f}') for g,xy in cent.items()])
add('TNF groups occupy the high-PC1 region, making treatment the dominant separation in this view. PC2 separates the two cell lines. These patterns support broad biological structure but do not substitute for a sample-identity assay.');add()
if bad:
 add('Nearest-neighbour checks found '+str(len(bad))+' sample outside its own group: '+ '; '.join(r['sample_id']+' → '+r['nearest'] for r in bad)+'. This involves the closely overlapping Crtc2 untreated and adenosine groups; the global two-component view does not resolve these conditions well. All other samples have a nearest neighbour in their own group. This isolated overlap is insufficient evidence of mislabelling, and all samples were retained.');add()
add('### 7.2 Differential Expression');add()
table(['Comparison','Genes tested','DEGs','DEGs / tested','Crtc2 log2FC','Crtc2 padj'],[(name(r['Comparison']),n(r['Genes_tested']),n(r['Total']),f'{100*int(r["Total"])/int(r["Genes_tested"]):.1f}%',f'{float(targets[r["Comparison"]]["log2FoldChange"]):+.2f}',f'{float(targets[r["Comparison"]]["padj"]):.3g}') for r in de])
add('The TNF line comparison has the largest DEG count. The two control-treatment comparisons show a much larger threshold-crossing response to PGE2 than to adenosine. These counts quantify the six specified comparisons; comparing DEG totals does not itself test an interaction.');add()
table(['Condition','Control mean Crtc2 normalised count','Crtc2-line mean','Crtc2 line / control'],[(t,f'{means["control-"+t]:.1f}',f'{means["Crtc2-"+t]:.1f}',f'{ratios[i]:.2f}') for i,t in enumerate(['un','TNF','PGE2','adenosine'])])
add('These target-gene group means use library-size normalisation estimated from the full count matrix, as provided in `Reads/Crtc2_expression_per_sample.csv`. The DE log2FC values above use the filtered DESeq2 model and shrinkage, so the two summaries need not be algebraically identical.');add()
add('Every comparison supplies a full tested-gene table, significant subset, full-range volcano plot, and top-50 DEG heatmap. Two adenosine comparisons also include a volcano_zoom.pdf companion with |log2FC| ≤ 3; the caption states the number of points outside the view, and the full tables and original plots retain those genes. Genes outside a comparison’s eligible set were not tested there, rather than demonstrated unchanged.');add()
add('### 7.3 Pathway Enrichment');add()
add('Hallmark normalised enrichment score (NES): positive means enrichment toward the numerator; negative toward the baseline. “n.s.” means the gene set did not appear among the significant outputs (padj < 0.05), not absence of the process.');add()
table(['Comparison']+[x.replace('_',' ') for x in H],[(name(k),*[f'{float(hrow(k,h)["NES"]):+.2f}' if hrow(k,h) else 'n.s.' for h in H]) for k in keys])
add('Inflammatory and interferon-related signatures are lower in the Crtc2 line in the four matched-condition comparisons. Higher cell-cycle signatures are seen in the untreated, PGE2 and adenosine line comparisons; the TNF condition does not show the same significant E2F/G2M pattern. In control cells, both PGE2 and adenosine show the opposite direction for these themes. These results motivate functional follow-up but do not prove direct CRTC2 regulation, pathway activity, or cell proliferation.');add()
E={(r['Comparison'],r['Analysis']):r for r in en}
cols=['GO_BP_UP','GO_BP_DOWN','KEGG_UP','KEGG_DOWN','GSEA_GO_BP','GSEA_KEGG','GSEA_Hallmark']
table(['Comparison']+cols,[(name(k),*[E[k,c]['N_significant'] for c in cols]) for k in keys])
add('Counts are significant terms at adjusted p-value < 0.05; overlapping terms are not independent biological discoveries. The Hallmark display above is a selected thematic view; all returned terms and the full 90-analysis summary are provided in `Enrichment/`. ORA UP / DOWN / ALL are separate tests, with the mapped testable genes as background.');add()
add('## 8. Conclusions, Limitations and Next Steps');add();add('### 8.1 Per-comparison interpretation');add()
interpret=['Baseline differences include inflammatory and cell-cycle expression themes.','Largest matched-condition line difference; lower inflammatory and interferon signatures in Crtc2.','Line difference under PGE2; it does not alone establish a difference in PGE2 responsiveness.','Line difference under adenosine; the Crtc2 untreated and adenosine groups overlap in PC1–PC2.','Broad control-cell PGE2-associated expression change, including inflammatory and cell-cycle themes.','Fewer gene-level DEGs than PGE2, with coordinated pathway shifts despite small gene-level effects.']
table(['Comparison','Meaning','Boundary'],[(name(r['Comparison']),interpret[i],'Expression association within this comparison; no formal interaction test.') for i,r in enumerate(de)])
add('### 8.2 Practical reading');add();add('The following evidence map gives the data location and inference boundary for each conclusion in Section 2.');add()
table(['Question','Measured evidence / delivered source','Interpretation and assumptions'],[(t,f'Section {sec}; {source}',boundary) for t,text,sec,source,boundary in claims])
add('### 8.3 Limitations');add()
for text in [
 'Three biological replicates per condition limit power for small or variable changes. A nonsignificant or untested gene is not shown to be unchanged.',
 'The specific cell line and the nature of the Crtc2 modification were not supplied. Reduced Crtc2 mRNA does not demonstrate loss of CRTC2 protein or function.',
 'Only six requested contrasts were tested. There is no formal line-by-treatment interaction test or delivered within-line TNF-response contrast. PCA overlap and differing DEG totals do not establish equivalent responses or loss of responsiveness.',
 'Name-based filtering restricts the statistical gene universe and is not a complete protein-coding biotype annotation. Low-expression and identifier-mapping filters also limit coverage. Filtered genes remain available in the full expression matrices.',
 'Pathway terms overlap and reflect gene-expression associations. Enrichment does not prove pathway activation, altered proliferation, or a direct molecular mechanism.',
 'Alignment and expression QC do not independently authenticate sample identity. The overlap of Crtc2 untreated and adenosine samples is reported explicitly rather than treated as proof of either mislabelling or equivalence.']:
 add('- '+text)
add();add('### 8.4 Suggested confirmation');add()
table(['Open question','Follow-up'],[
 ('Is CRTC2 protein reduced?','Western blot comparing control and Crtc2 cells under matched conditions.'),
 ('Do the transcript changes replicate?','Quantitative polymerase chain reaction (qPCR) for Crtc2 and selected DEGs in independent cultures.'),
 ('Does treatment responsiveness differ between lines?','A prespecified line-by-treatment interaction analysis on the existing samples; this would be an additional statistical analysis, not independent biological validation.'),
 ('Do pathway themes have functional consequences?','Independent proliferation measurements and cytokine-protein assays selected for the biological question.')])
add('## 9. Deliverable Files');add()
manifest=[('Reads/All_sample_gene_counts.tsv','Unfiltered Salmon gene-level estimated counts, 24 samples'),('Reads/All_sample_gene_tpm.tsv','Unfiltered gene-level TPM, 24 samples'),('Reads/library_metrics_per_sample.csv','Raw read pairs, estimated-count sums and detected genes'),('Reads/Crtc2_expression_per_sample.csv','Target-gene normalised counts for 24 samples'),('QC/multiqc/','Self-contained MultiQC HTML report and supporting data; use the HTML for interactive QC'),('DE_PCA_Results/PCA_AllSamples_Overview.pdf','24-sample PCA'),('DE_PCA_Results/PCA_coordinates_overview.csv','24 sample coordinates and grouping'),('DE_PCA_Results/SampleDistance_Heatmap.pdf','Sample distances with group annotation'),('DE_PCA_Results/normalized_counts_all_samples.csv','Filtered-model normalised counts'),('DE_PCA_Results/DEG_summary_all_comparisons.csv','Six-comparison gene counts at both effect-size thresholds'),('DE_PCA_Results/<comparison>/','Full and significant tables, full-range volcano and top-50 heatmap for each of six comparisons; two adenosine comparisons also have a volcano_zoom.pdf companion'),('Enrichment/<comparison>/{GO,KEGG,GSEA}/','Directional ORA and ranked GSEA tables and PDF plots'),('Enrichment/Enrichment_summary_all_comparisons.csv','90-analysis status and significant-term counts'),('Gene_filtering/','Full-gene and six per-comparison eligibility audits'),('Claim_evidence_map.csv','Measured facts, interpretation boundaries and source files for Section 2'),('mouse_Gene_annotation_20260628.xlsx','Gene annotation reference'),(filename,'Report and reading entry point')]
table(['Path','Content'],[('`'+a+'`',b) for a,b in manifest]);add('---');add();add('*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*');add()
text='\n'.join(line)
assert not any(x in text for x in ['/home/','DE_R45','P51','the same answers','restates the'])
# Output additions are new, exclusive creations; existing results are untouched.
for target,content in [(D/'Reads/Crtc2_expression_per_sample.csv',(P/'scripts/rds_cache/crtc2_expression_internal.csv').read_bytes())]:
 with target.open('xb') as f:f.write(content)
F=D/'Gene_filtering';F.mkdir(exist_ok=False)
for src in (P/'scripts/rds_cache').glob('filter_audit_*.csv'):
 if src.name=='filter_audit_added_vs_n-2.csv':continue
 with (F/src.name).open('xb') as f:f.write(src.read_bytes())
with (D/'Claim_evidence_map.csv').open('x') as f:
 w=csv.writer(f);w.writerow(['Claim','Measured_fact','Report_section','Data_source','Interpretation_boundary'])
 for t,txt,sec,source,boundary in claims:w.writerow([t,txt,sec,source,boundary])
with report.open('x') as f:f.write(text)
print('Report generated from P52 only:',report)
print('Six contrasts; 24 samples; 90 enrichment calls; no additional formal comparison or refit')
