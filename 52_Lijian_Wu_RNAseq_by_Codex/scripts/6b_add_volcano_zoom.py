#!/usr/bin/env python3
"""P52 display-only companions: retain all original volcano PDFs and tables."""
from pathlib import Path
import csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1];D=P/'Data_Analysis_20261109';review=P/'logs/figure_review_zoom_v2_20261109';review.mkdir(exist_ok=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.linewidth':1.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
audit=[]
for k in ['Crtc2-adenosine_vs_control-adenosine','control-adenosine_vs_control-un']:
 with (D/'DE_PCA_Results'/k/(k+'_ALL_genes.csv')).open() as f:rows=list(csv.DictReader(f))
 x=np.array([float(r['log2FoldChange']) for r in rows]);pv=np.array([float(r['pvalue']) if r['pvalue'] not in ['NA',''] else float('nan') for r in rows]);y=-np.log10(np.maximum(pv,np.finfo(float).tiny));valid=np.isfinite(x)&np.isfinite(y);visible=valid&(abs(x)<=3)
 fig,ax=plt.subplots(figsize=(8,6.5));fig.subplots_adjust(left=.13,right=.96,bottom=.23,top=.82)
 for group,color,label in [('NS','#B3B3B3','Not significant'),('DOWN','#3B5DA4','Lower in numerator'),('UP','#BC3535','Higher in numerator')]:
  m=visible&np.array([r['direction']==group for r in rows]);ax.scatter(x[m],y[m],s=8,c=color,alpha=.6,label=label,edgecolors='none')
 for cut in [-.263,.263]:ax.axvline(cut,color='#666666',ls='--',lw=1)
 ax.set(xlim=(-3,3),ylim=(0,float(y[visible].max())*1.08),xlabel='Shrunken log2 fold change',ylabel='−log10(p-value)')
 ax.grid(alpha=.15);ax.set_axisbelow(True)
 fig.suptitle(k.replace('_vs_',' vs '),fontsize=13,y=.96)
 fig.text(.13,.88,f'Zoom: |log2FC| ≤ 3; {int(visible.sum()):,} genes shown; {int((valid&~visible).sum()):,} outside x range; {int((~valid).sum()):,} unavailable p',fontsize=11)
 ax.legend(loc='upper center',bbox_to_anchor=(.5,-.15),ncol=3,frameon=False,fontsize=10)
 fig.text(.13,.035,'Significance: padj ≤ 0.05 and |log2FC| ≥ 0.263.\nFull-range volcano and all result rows are retained; this zoom does not change DEG counts.',fontsize=9)
 f=D/'DE_PCA_Results'/k/(k+'_volcano_zoom.pdf');assert not f.exists();fig.savefig(f,metadata={'CreationDate':None,'ModDate':None});fig.savefig(review/(k+'.png'),dpi=180);plt.close(fig)
 outside=[r for r,v in zip(rows,visible) if not v]
 with (review/(k+'_outside_view.csv')).open('x') as fh:
  w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(outside)
 audit.append({'comparison':k,'total_rows':len(rows),'visible_points':int(visible.sum()),'outside_x_range':int((valid&~visible).sum()),'unavailable_p':int((~valid).sum()),'outside_view':len(outside),'pdf':str(f.relative_to(D))})
 assert int(visible.sum())+len(outside)==len(rows)
with (review/'point_counts.json').open('x') as f:json.dump(audit,f,indent=2)
print('Added two display-only zoom PDFs; original plots/tables untouched',audit)
