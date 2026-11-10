import csv,gzip,concurrent.futures,subprocess,datetime,json
from pathlib import Path
P=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((P/'scripts/0_sample_table.csv').open()))
assert len(rows)==24
files=[r[k] for r in rows for k in ('fastq_r1','fastq_r2')]
assert len(set(files))==48
print('VERIFY_START',datetime.datetime.now().isoformat(),len(files),flush=True)
def check(f):
 subprocess.run(['gzip','-t',f],check=True)
 print('GZIP_OK',f,flush=True)
 return f
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool: list(pool.map(check,files))
with (P/'scripts/nf_core_samplesheet.csv').open('w') as h:
 w=csv.writer(h);w.writerow(['sample','fastq_1','fastq_2','strandedness'])
 for r in rows:w.writerow([r['sample_id'],r['fastq_r1'],r['fastq_r2'],'auto'])
(P/'logs/FASTQ_VERIFIED.marker').write_text(datetime.datetime.now().isoformat())
print('FASTQ_VERIFIED 48/48',flush=True)
