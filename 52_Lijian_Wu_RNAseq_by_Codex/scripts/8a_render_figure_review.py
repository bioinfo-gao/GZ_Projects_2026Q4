#!/usr/bin/env python3
"""Render every P52 PDF to internal review directory; contact sheets are not deliverables."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from PIL import Image,ImageDraw,ImageFont
import subprocess,json
P=Path(__file__).resolve().parents[1];D=P/'Data_Analysis_20261109';O=P/'logs/figure_review_20261109';O.mkdir(exist_ok=False)
files=sorted(D.rglob('*.pdf'))
def render(pair):
 i,f=pair;out=O/f'{i:03d}'
 subprocess.run(['pdftoppm','-f','1','-singlefile','-scale-to','1100','-png',str(f),str(out)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
 return {'index':i,'source':str(f.relative_to(D)),'preview':str(out)+'.png'}
with ThreadPoolExecutor(max_workers=4) as pool:manifest=list(pool.map(render,enumerate(files,1)))
(O/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
# 12 figures/sheet, large enough to screen gross clipping; inspect source previews individually too.
font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',13)
for start in range(0,len(manifest),12):
 subset=manifest[start:start+12];canvas=Image.new('RGB',(1800,((len(subset)+2)//3)*540),'#dddddd');draw=ImageDraw.Draw(canvas)
 for j,r in enumerate(subset):
  img=Image.open(r['preview']).convert('RGB');img.thumbnail((594,480));x=(j%3)*600+(600-img.width)//2;y=(j//3)*540+48
  canvas.paste(img,(x,y));label=f'{r["index"]:03d} '+r['source'];label=label.replace('Enrichment/','').replace('DE_PCA_Results/','')
  for k in range(0,min(len(label),138),69):draw.text(((j%3)*600+4,(j//3)*540+2+(k//69)*16),label[k:k+69],font=font,fill='black')
 canvas.save(O/f'contact_{start//12+1:02d}.jpg',quality=90)
print(f'Rendered {len(manifest)} PDFs; {(len(manifest)+11)//12} contact sheets in {O}')
