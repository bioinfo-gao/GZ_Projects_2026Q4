import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp, sys
print(sys.executable, sc.__version__)
rng=np.random.default_rng(0); n,g=1200,500
X=sp.csr_matrix(rng.negative_binomial(2,0.7,(n,g)).astype(np.int32))
a=ad.AnnData(X,obs=pd.DataFrame({'sample':np.repeat(['E1','E2','E3','C1','C2','C3'],200),'state':rng.choice(['s1','s2'],n)},index=[f'c{i}' for i in range(n)]),var=pd.DataFrame(index=[f'MT-{i}' if i<13 else f'G{i}' for i in range(g)]))
a.obs['group']=a.obs['sample'].str[0]; a.layers['counts']=a.X.copy()
a.var['mt']=a.var_names.str.startswith('MT-'); assert a.var['mt'].sum()==13
sc.pp.calculate_qc_metrics(a,qc_vars=['mt'],inplace=True,percent_top=None,log1p=True,layer='counts')
def mad_out(s,nm=5,side='both'):
    m=s.median(); d=(s-m).abs().median()*1.4826
    return {'both':(s<m-nm*d)|(s>m+nm*d),'low':s<m-nm*d,'high':s>m+nm*d}[side]
a.obs['qc_flag_low_lib']=a.obs.groupby('sample',observed=True)['log1p_total_counts'].transform(lambda s:mad_out(s,5,'low'))
pb=sc.get.aggregate(a,by=['sample','state'],func='sum',layer='counts')
print(pb.obs.head(3)); M=pd.DataFrame(pb.layers['sum'].toarray() if sp.issparse(pb.layers['sum']) else pb.layers['sum'],index=pb.obs_names,columns=pb.var_names)
ncell=a.obs.groupby(['sample','state'],observed=True).size().rename('n_cells').reset_index()
assert (M.values==M.values.astype(int)).all()
tot=np.asarray(a.layers['counts'].sum()); assert M.values.sum()==tot, (M.values.sum(),tot)
col=pb.obs.reset_index(drop=True).merge(ncell,on=['sample','state']); col.index=M.index
col['group']=col['sample'].str[0]; keep=col[(col.state=='s1')&(col.n_cells>=10)]
M.loc[keep.index].T.astype(int).to_csv('pb.csv'); keep.to_csv('col.csv')
print(M.shape,keep[['sample','n_cells','group']].to_string())
