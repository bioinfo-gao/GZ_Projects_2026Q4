import gzip,sys
base='/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/'
ref='/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/GRCh38.primary_assembly.genome.fa'
# pull chrM
seq=[]; on=False
with open(ref) as fh:
    for line in fh:
        if line.startswith('>'):
            name=line[1:].split()[0]
            if on: break
            on = (name=='chrM')
            if on: print('found contig chrM', file=sys.stderr)
            continue
        if on: seq.append(line.strip())
chrM=''.join(seq).upper()
print('chrM length: %d'%len(chrM))
comp=str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(comp)[::-1]
rows=[]
for line in open(base+'01_analysis/N1-ABE-N2-ABE-N1-CBE-N2-CBE/probe_set.csv'):
    if line.startswith('#') or line.startswith('gene_id'): continue
    p=line.rstrip('\n').split(',')
    gid,pseq,pid,incl,region=p[0],p[1],p[2],p[3],p[4]
    gname=pid.split('|')[1]
    if not gname.startswith('MT-'): continue
    pos=chrM.find(pseq); strand='+'
    if pos<0:
        pos=chrM.find(rc(pseq)); strand='-'
    rows.append((gname,pid,incl,region,strand,pos))
print('%-8s %-10s %-9s %-6s %s'%('gene','included','strand','start','footprint chrM 1-based (ligation junction after 25nt)'))
for gname,pid,incl,region,strand,pos in sorted(rows,key=lambda r:r[5]):
    if pos<0:
        print('%-8s %-10s %-9s %-6s NOT FOUND in chrM'%(gname,incl,'?','-')); continue
    s1=pos+1; e1=pos+50
    jleft, jright = (s1+24, s1+25)
    print('%-8s %-10s %-9s %-6d chrM:%d-%d   junction between %d/%d'%(gname,incl,strand,s1,s1,e1,jleft,jright))
