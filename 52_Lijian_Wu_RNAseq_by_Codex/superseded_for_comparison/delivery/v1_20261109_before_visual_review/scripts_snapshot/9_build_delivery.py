#!/usr/bin/env python3
"""P52 first delivery build: audited source -> new dated directory; never refresh in place."""
from pathlib import Path
import os,sys,shutil,subprocess,hashlib,json
from datetime import datetime
from zoneinfo import ZoneInfo
P=Path(__file__).resolve().parents[1];D=P/'Data_Analysis_20261109';date=sys.argv[1] if len(sys.argv)>1 else datetime.now(ZoneInfo('America/New_York')).strftime('%Y-%m-%d')
O=P/('custom_research_report_'+date.replace('-',''));assert not O.exists(),'Refuse existing delivery directory'
report=D/('Crtc2_Line_Treatment_RNAseq_Report_'+date[5:7]+date[8:10]+'.md');assert report.is_file()
# Source must have just passed the explicit content audit.
with (P/'logs/audit_source_20261109.log').open('x') as f:subprocess.run([sys.executable,str(P/'scripts/8_pre_delivery_audit.py'),str(D)],check=True,stdout=f,stderr=subprocess.STDOUT)
subprocess.run([sys.executable,'/home/gao/.claude/hooks/checks/check_report.py',str(report)],check=True)
shutil.copytree(D,O)
# Byte identity audit for EVERY file, not only report or a subset.
def hashes(root):return {str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in root.rglob('*') if f.is_file()}
a=hashes(D);b=hashes(O);assert a==b,'Source/delivery bytes differ'
with (P/'logs/delivery_sha256_20261109.json').open('x') as f:json.dump({'source':str(D),'delivery':str(O),'files':a},f,indent=2)
with (P/'logs/audit_delivery_20261109.log').open('x') as f:subprocess.run([sys.executable,str(P/'scripts/8_pre_delivery_audit.py'),str(O)],check=True,stdout=f,stderr=subprocess.STDOUT)
subprocess.run([sys.executable,'/home/gao/.claude/hooks/checks/check_report.py',str(O/report.name)],check=True)
print('DELIVERY_BUILT_AUDITED',O,'files',len(b),'bytes',sum(f.stat().st_size for f in O.rglob('*') if f.is_file()),flush=True)
# Human/agent visual inspection must be recorded before the active project becomes complete.
