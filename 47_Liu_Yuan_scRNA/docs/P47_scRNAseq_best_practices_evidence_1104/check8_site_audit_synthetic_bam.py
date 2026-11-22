import pysam, random
random.seed(2); L=2000; ref=''.join(random.choice('ACGT') for _ in range(L)); pos=1000
refb=ref[pos]; altb={'A':'G','G':'A','C':'T','T':'C'}[refb]
hdr={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chrT','LN':L}]}
R=[]
def rd(name, cb, ub, base, flag=0, mapq=255):
    st=pos-random.randint(10,80); s=list(ref[st:st+90]); s[pos-st]=base
    a=pysam.AlignedSegment(); a.query_name=name; a.query_sequence=''.join(s); a.flag=flag
    a.reference_id=0; a.reference_start=st; a.mapping_quality=mapq; a.cigar=((0,90),)
    a.query_qualities=pysam.qualitystring_to_array('I'*90)
    if cb: a.set_tag('CB',cb)
    if ub: a.set_tag('UB',ub)
    R.append(a)
A,B='AAAACCCCGGGGTTTT-1','TTTTGGGGCCCCAAAA-1'
for u in range(30): rd(f'A{u}',A,f'UA{u:04d}',altb)
rd('A0_dup1',A,'UA0000',altb); rd('A0_dup2',A,'UA0000',altb)          # 同一分子 3 条 read
rd('A1_conf1',A,'UA0001',refb); rd('A1_conf2',A,'UA0001',altb)        # 分子 UA0001: alt,ref,alt => 冲突, 多数=alt
for u in range(30): rd(f'B{u}',B,f'UB{u:04d}',refb)
rd('B0_sec',B,'UB0000',altb,flag=256)                                   # secondary
rd('B_lowq',B,'UB9999',altb,mapq=3)                                     # 低 MAPQ
rd('B_noUB',B,None,altb)                                                # 无 UB
rd('C_notcell','GGGGGGGGGGGGGGGG-1','UC0000',altb)                      # 不在 QC 细胞
# 共 30+2+2+30+4 = 68 ... 计数见下
R.sort(key=lambda r:r.reference_start)
with pysam.AlignmentFile('t2.bam','wb',header=hdr) as o:
    for r in R: o.write(r)
pysam.index('t2.bam')
open('cells.tsv','w').write(A+'\n'+B+'\n')
open('sites.tsv','w').write(f'chrT\t{pos+1}\t{refb}\t{altb}\tsiteT\n')
print('n_records', len(R))
