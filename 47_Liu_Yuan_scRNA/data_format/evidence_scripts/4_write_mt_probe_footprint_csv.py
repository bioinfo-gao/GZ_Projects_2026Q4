"""Emit the canonical machine-readable chrM probe-footprint table (SSOT).
Footprint = the 50 nt of chrM that the Flex probe pair hybridises to.
Ligation junction = between local offsets 25/26 (LHS and RHS half-probe, 25 nt each)."""
import csv
base='/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/'
ref='/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa'
out='/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/data_format/chrM_probe_footprints.csv'
seq=[]; on=False
with open(ref) as fh:
    for line in fh:
        if line.startswith('>'):
            if on: break
            on = (line[1:].split()[0]=='chrM'); continue
        if on: seq.append(line.strip())
chrM=''.join(seq).upper()
assert len(chrM)==16569, len(chrM)
comp=str.maketrans('ACGTN','TGCAN')
rows=[]
for line in open(base+'01_analysis/N1-ABE-N2-ABE-N1-CBE-N2-CBE/probe_set.csv'):
    if line.startswith('#') or line.startswith('gene_id'): continue
    p=line.rstrip('\n').split(',')
    gene_id,pseq,probe_id,included,region=p[0],p[1],p[2],p[3],p[4]
    gname=probe_id.split('|')[1]
    if not gname.startswith('MT-'): continue
    pos=chrM.find(pseq); strand='+'
    if pos<0:
        pos=chrM.find(pseq.translate(comp)[::-1]); strand='-'
    assert pos>=0, probe_id
    s1=pos+1
    rows.append(dict(gene_name=gname, gene_id=gene_id, probe_id=probe_id,
                     included_in_filtered_probe_set=included, probe_region=region,
                     chrM_start_1based=s1, chrM_end_1based=s1+49, strand=strand,
                     ligation_junction_between=f'{s1+24}/{s1+25}', probe_seq=pseq))
rows.sort(key=lambda r: r['chrM_start_1based'])
with open(out,'w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
cov=50*len(rows); inc=50*sum(r['included_in_filtered_probe_set']=='TRUE' for r in rows)
print(f'wrote {out}: {len(rows)} probes')
print(f'chrM covered by any probe : {cov} bp = {100*cov/16569:.2f}% of 16569')
print(f'chrM covered, included only: {inc} bp = {100*inc/16569:.2f}%')
