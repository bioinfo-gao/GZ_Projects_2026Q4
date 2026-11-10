#!/usr/bin/env python3
"""P52 independent content/manifest checks against supplied matrices and sample SSOT."""
from pathlib import Path
import csv,json,math,re,sys,subprocess,collections
import numpy as np
P=Path(__file__).resolve().parents[1]; D=Path(sys.argv[1]) if len(sys.argv)>1 else P/'Data_Analysis_20261109'
def read(p,sep=','):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter=sep))
checks=[]
def check(ok,label):
 checks.append((bool(ok),label));print(('PASS ' if ok else 'FAIL ')+label,flush=True)
def number(s):
 try:return float(s)
 except (ValueError,TypeError):return float('nan')
st=read(P/'scripts/0_sample_table.csv'); ids=[r['sample_id'] for r in st];groups={r['sample_id']:r['group'] for r in st};de=read(D/'DE_PCA_Results/DEG_summary_all_comparisons.csv')
check(len(ids)==len(set(ids))==24,'24 distinct SSOT samples')
check(sorted(collections.Counter(groups.values()).values())==[3]*8,'eight groups with three samples each')
check(all(r['paired_on_disk']=='yes' and r['mapping_status']=='confirmed_by_user_2026-11-09' for r in st),'paired files and confirmed client mapping')
check(len(de)==6,'six requested contrasts')
raw=read(D/'Reads/All_sample_gene_counts.tsv','\t');tpm=read(D/'Reads/All_sample_gene_tpm.tsv','\t')
check(set(ids)<=set(raw[0]) and set(ids)<=set(tpm[0]),'both full expression matrices include all 24 samples')
check([r['gene_id'] for r in raw]==[r['gene_id'] for r in tpm],'count/TPM gene order agrees')
geneids=[r['gene_id'] for r in raw];genes=[r['gene_name'] for r in raw]
counts=np.rint(np.array([[float(r[s]) for s in ids] for r in raw]));check(np.isfinite(counts).all() and (counts>=0).all(),'finite nonnegative count matrix')
ribo=re.compile(r'^Rpl|^Rps|^Mrpl|^Mrps|^Rplp|^Rpsa',re.I);look=re.compile(r'^Rps6k|^Rps19bp1$|^Rps23rg1$',re.I);other=re.compile(r'^mt-|^Snord|^Snora|^Rnu|^Malat1|^Neat1|^Xist|^Mir[0-9]|^let-|^Linc|^Gm[0-9]|pseudogene|antisense',re.I)
regex=np.array([not(bool(ribo.search(g)) and not bool(look.search(g))) and not bool(other.search(g)) for g in genes]);globalpass=(counts>=10).sum(axis=1)>=8
full=read(D/'Gene_filtering/filter_audit_all_genes.csv');check(len(full)==len(geneids)==57132,'full filter audit covers 57,132 input genes')
check(all((r['passes_regex']=='TRUE')==bool(x) for r,x in zip(full,regex)),'gene-name inclusion independently reproduced')
union=set();alls={}
for summary in de:
 k=summary['Comparison'];num=summary['Numerator'];den=summary['Baseline'];sub=D/'DE_PCA_Results'/k
 a=read(sub/(k+'_ALL_genes.csv'));alls[k]=a;b=read(sub/(k+'_DEG_significant.csv'));cols=[i for i,s in enumerate(ids) if groups[s] in [num,den]]
 eligibility=regex & (((counts[:,cols]>=10).sum(axis=1)>=3)|globalpass);expected={g for g,v in zip(geneids,eligibility) if v};actual={r['gene_id'] for r in a};union|=actual
 check(len(cols)==6 and actual==expected and len(actual)==len(a),k+': tested IDs exactly equal independently computed eligible set')
 audit=read(D/'Gene_filtering'/('filter_audit_'+k+'.csv'));ad={r['gene_id'] for r in audit if r['eligible']=='TRUE'};check(ad==actual,k+': per-gene audit and tested IDs agree')
 pv=np.array([number(r['pvalue']) for r in a]);adj=np.array([number(r['padj']) for r in a]);lfc=np.array([number(r['log2FoldChange']) for r in a]);ok=np.isfinite(pv);ii=np.where(ok)[0];order=ii[np.argsort(pv[ok])];bh=np.full(len(a),np.nan)
 if len(order):bh[order]=np.minimum(1,np.minimum.accumulate((pv[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
 check(np.allclose(bh,adj,equal_nan=True,atol=1e-10,rtol=1e-8),k+': BH recalculated from exported p-values')
 sig=np.isfinite(adj)&(adj<=.05)&(abs(lfc)>=.263);sigids={r['gene_id'] for r,v in zip(a,sig) if v}
 check(sigids=={r['gene_id'] for r in b} and all((r['significant']=='TRUE')==bool(v) for r,v in zip(a,sig)),k+': threshold reproduces full and significant tables')
 check(int(summary['Up'])==int(sum(sig&(lfc>0))) and int(summary['Down'])==int(sum(sig&(lfc<0))) and int(summary['Total'])==int(sum(sig)) and int(summary['Genes_tested'])==len(a),k+': summary counts match tables')
 strict=np.isfinite(adj)&(adj<=.05)&(abs(lfc)>=1);check(int(summary['Up_lfc1'])==sum(strict&(lfc>0)) and int(summary['Down_lfc1'])==sum(strict&(lfc<0)),k+': twofold summary matches tables')
 check(all(math.isfinite(number(r['stat'])) for r in a),k+': all GSEA Wald statistics finite')
 check(all(all('normcount_'+s in a[0] for s in ids if groups[s] in [num,den]) for _ in [0]),k+': six normalised-count sample columns exist')
 # Verify top-50 heatmap sample labels in actual PDF text, not just existence.
 for kind in ['volcano','top50_heatmap']:
  f=sub/(k+'_'+kind+'.pdf');check(f.is_file() and f.stat().st_size>1000,k+': '+kind+' PDF exists')
  if kind=='top50_heatmap':
   txt=subprocess.check_output(['pdftotext',str(f),'-'],text=True);check(all(s in txt for s in ids if groups[s] in [num,den]),k+': six heatmap sample labels rendered in PDF')
check(len(union)==12327,'model gene universe is union of six tested sets')
check({r['gene_id'] for r in full if r['retained_for_model']=='TRUE'}==union,'full audit retained flags equal model universe')
coords=read(D/'DE_PCA_Results/PCA_coordinates_overview.csv');check(len(coords)==24 and {r['name'] for r in coords}==set(ids),'PCA coordinate count and sample set 24/24')
check(all(r['Group']==groups[r['name']] for r in coords),'PCA group labels agree with SSOT')
xy=np.array([[float(r['PC1']),float(r['PC2'])] for r in coords]);within=[];between=[]
for i in range(24):
 for j in range(i+1,24):
  (within if coords[i]['Group']==coords[j]['Group'] else between).append(float(np.linalg.norm(xy[i]-xy[j])))
facts=json.loads((P/'logs/report_facts_20261109.json').read_text());pc=facts['pca']
check(np.isclose(np.mean(within),pc['within_dist']) and np.isclose(np.mean(between),pc['between_dist']),'PCA distance summaries independently reproduced')
check(pc['plotted_points']==pc['plotted_labels']==24,'ggplot build: 24 finite points and 24 labels')
pdftext=subprocess.check_output(['pdftotext',str(D/'DE_PCA_Results/PCA_AllSamples_Overview.pdf'),'-'],text=True)
check(all(s in pdftext for s in ids),'actual PCA PDF renders all 24 sample labels')
en=read(D/'Enrichment/Enrichment_summary_all_comparisons.csv');check(len(en)==90 and all(r['status']=='ok' for r in en),'all 90 enrichment analyses recorded ok')
for r in en:
 k=r['Comparison'];an=r['Analysis'];kind=an.split('_')[0];suffix=an
 if an=='GSEA_Hallmark':suffix='GSEA_MSigDB_Hallmark'
 f=D/'Enrichment'/k/kind/(suffix+'.csv');t=read(f)
 check(len(t)==int(r['N_significant']),k+': '+an+' significant term count matches CSV')
 check(all(math.isfinite(number(x['p.adjust'])) and number(x['p.adjust'])<.05 for x in t),k+': '+an+' terms meet stated adjusted-p threshold')
reports=list(D.glob('*.md'));check(len(reports)==1,'exactly one report in package')
text=reports[0].read_text();date=re.search(r'\*\*Report Date:\*\* (\d{4}-\d{2}-\d{2})',text).group(1)
check(reports[0].name.endswith(date[5:7]+date[8:10]+'.md'),'report body and filename dates agree')
if D.name.startswith('custom_research_report_'):check(D.name=='custom_research_report_'+date.replace('-',''),'delivery folder and report dates agree')
for h in ['1. Objectives','2. Core Conclusions','3. Key Findings','4. Sample Information','5. Analysis Rationale and Decision Criteria','6. Methods','7. Results','8. Conclusions, Limitations and Next Steps','9. Deliverable Files']:check('## '+h in text,'report section '+h)
check(all('### '+x in text for x in ['7.1 Principal Component Analysis (PCA)','8.2 Practical reading','8.3 Limitations','8.4 Suggested confirmation']),'PCA interpretation and detailed conclusions/limitations present')
heads=set(re.findall(r'^#{2,3} (\d+(?:\.\d+)?)',text,re.M));refs=re.findall(r'Section (\d+(?:\.\d+)?)',text);check(all(x in heads for x in refs),'all section references resolve')
for r in de:
 row='| '+r['Comparison'].replace('_vs_',' vs ')+' | '+f'{int(r["Up"]):,}'+' | '+f'{int(r["Down"]):,}'+' | '+f'{int(r["Total"]):,}'+' |'
 check(row in text,r['Comparison']+': report DEG row matches delivered data')
check(all(s in text and f'{int(next(r["raw_read_pairs"] for r in st if r["sample_id"]==s)):,}' in text for s in ids),'report lists all samples and their read-pair values')
check('Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics' in text and 'cell line name not specified by client' in text,'signature and metadata boundary present')
check(not re.search(r'/home/|DE_R45|P51|long-read|Nanopore|PacBio|the same answers|restates the',text,re.I),'no local configuration, cross-project narrative or prohibited recommendations in report')
sep=[s for s in text.splitlines() if s.startswith('| :---')];check(all(re.fullmatch(r'\| :--- \|( :---: \|)+',s) for s in sep),'table separator format')
manifest=['Reads/All_sample_gene_counts.tsv','Reads/All_sample_gene_tpm.tsv','Reads/library_metrics_per_sample.csv','Reads/Crtc2_expression_per_sample.csv','QC/multiqc','DE_PCA_Results/PCA_AllSamples_Overview.pdf','DE_PCA_Results/PCA_coordinates_overview.csv','DE_PCA_Results/SampleDistance_Heatmap.pdf','DE_PCA_Results/normalized_counts_all_samples.csv','DE_PCA_Results/DEG_summary_all_comparisons.csv','Enrichment/Enrichment_summary_all_comparisons.csv','Gene_filtering','Claim_evidence_map.csv','mouse_Gene_annotation_20260628.xlsx']
for f in manifest:check((D/f).exists() and f in text,'fixed manifest and report index: '+f)
check(len(read(D/'Claim_evidence_map.csv'))==6,'six core conclusions have source and inference boundary')
check(not any(p.suffix in ['.rds','.RData','.log'] for p in D.rglob('*') if p.is_file()),'client package excludes internal caches/logs')
# Decode every PDF, including enrichment figures; catch damaged/empty pages.
pdfs=list(D.rglob('*.pdf'));check(len(pdfs)==174,'174 expected result PDFs')
for f in pdfs:
 proc=subprocess.run(['pdfinfo',str(f)],capture_output=True,text=True)
 if proc.returncode or not re.search(r'Pages:\s+[1-9]',proc.stdout):check(False,'PDF decodes: '+str(f.relative_to(D)))
check(all(subprocess.run(['pdfinfo',str(f)],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode==0 for f in pdfs),'all 174 PDFs decode successfully')
print(f'AUDIT SUMMARY: passed={sum(x[0] for x in checks)} failed={sum(not x[0] for x in checks)}',flush=True)
if not all(x[0] for x in checks):sys.exit(1)
