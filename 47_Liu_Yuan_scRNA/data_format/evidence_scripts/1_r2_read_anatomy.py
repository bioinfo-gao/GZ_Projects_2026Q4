import gzip
base='/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/'
ps=base+'01_analysis/N1-ABE-N2-ABE-N1-CBE-N2-CBE/probe_set.csv'
probes={}
for line in open(ps):
    if line.startswith('#') or line.startswith('gene_id'): continue
    p=line.rstrip('\n').split(',')
    probes[p[1]]=p[2]
print('n probes', len(probes), 'probe len set', sorted({len(k) for k in probes}))
comp=str.maketrans('ACGTN','TGCAN')
rc={k.translate(comp)[::-1]:v for k,v in probes.items()}
f=base+'00_fastq/N1-ABE-N2-ABE-N1-CBE-N2-CBE_S1_L001_R2_001.fastq.gz'
n=0; hf=0; hr=0; ex=[]
with gzip.open(f,'rt') as fh:
    for i,line in enumerate(fh):
        if i%4!=1: continue
        s=line.strip(); n+=1
        a=s[:50]
        if a in probes:
            hf+=1
            if len(ex)<3: ex.append(('fwd',probes[a],s))
        elif a in rc:
            hr+=1
            if len(ex)<3: ex.append(('rc',rc[a],s))
        if n>=20000: break
print('R2 first50 exact match: fwd=%d rc=%d of n=%d (%.1f%%)'%(hf,hr,n,100*(hf+hr)/n))
for tag,pid,s in ex:
    print(tag,pid)
    print('   full90      :',s)
    print('   [0:50]probe :',s[:50])
    print('   [50:80]     :',s[50:80])
    print('   [80:90]     :',s[80:90])
