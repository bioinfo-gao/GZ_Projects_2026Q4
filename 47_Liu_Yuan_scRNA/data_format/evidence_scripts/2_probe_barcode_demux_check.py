import gzip, json, collections
base='/srv/external_data/gao_projects/47_Liu_Yuan_scRNA_raw/genewiz_30-1371703644_1124/30-1371703644/'
pools={'pool1_edited':'N1-ABE-N2-ABE-N1-CBE-N2-CBE',
       'pool2_control':'N1-Control-N2-Control-CBE-Control-ABE-Control'}
for pname,lib in pools.items():
    cpt=json.load(open(base+'01_analysis/'+lib+'/multiplexing_analysis/cells_per_tag.json'))
    # map probe-barcode-id -> the 10bp suffix carried in Cell Ranger cell barcodes
    tag2suf={}
    for tag,bcs in cpt.items():
        if not bcs: continue
        sufs={b.split('-')[0][16:] for b in bcs}
        assert len(sufs)==1, (tag,sufs)
        tag2suf[tag]=sufs.pop()
    # sample_id mapping from config.csv
    tag2samp={}
    insec=False
    for line in open(base+'01_analysis/'+lib+'/config.csv'):
        line=line.strip()
        if line=='[samples]': insec=True; continue
        if insec and line.startswith('['): break
        if insec and line and not line.startswith('sample_id'):
            p=line.split(','); tag2samp[p[1]]=p[0]
    print('=== %s  (%s)'%(pname,lib))
    for t in sorted(tag2suf): print('   %s -> %-12s probe-BC seq %s'%(t,tag2samp[t],tag2suf[t]))
    wanted={suf:t for t,suf in tag2suf.items()}
    cnt=collections.Counter(); n=0
    with gzip.open(base+'00_fastq/'+lib+'_S1_L001_R2_001.fastq.gz','rt') as fh:
        for i,line in enumerate(fh):
            if i%4!=1: continue
            s=line.strip(); n+=1
            tail=s[-20:]
            for suf,t in wanted.items():
                if suf in tail: cnt[t]+=1; break
            else: cnt['no_match']+=1
            if n>=200000: break
    print('   --- scan of %d R2 reads for those 4 probe barcodes:'%n)
    for t,c in sorted(cnt.items()):
        print('       %-10s %8d  %5.1f%%'%(t,c,100*c/n))
