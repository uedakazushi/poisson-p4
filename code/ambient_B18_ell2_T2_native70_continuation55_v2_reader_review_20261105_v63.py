#!/usr/bin/env python3
"""Metadata-only V2 repair peer; no finite/G8/image/rank arithmetic."""
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_V2_READER_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out):
 start=time.monotonic();out.mkdir(exist_ok=False,parents=True);pins={}
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 first='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_INPUT_GATE_20261105_v63_REVIEW/report.json';p=read(first);assert p['all_passed'] and pins[first]=='c7146d77a145b8e89c1e8568aa496799d1bc629408db744e38a54553bfa33786'
 oldpp='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63_REVIEW/report.json';oldpeer=read(oldpp);assert oldpeer['all_passed'] and p['input_sha256'][oldpp]==pins[oldpp]
 source='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63/report.json';g=read(source);assert g['all_passed'] and pins[source] in json.dumps(oldpeer)
 cg='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_INPUT_GATE_20261105_v63/report.json';cont=read(cg);assert cont['all_passed'] and p['input_sha256'][cg]==pins[cg]
 cp='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_V2_PREPARED_CUSTODY_20261105_v63/report.json';c=read(cp);assert c['all_passed'] and pins[cp]=='87ac74f5f2d4c81bab4fafdf31cbbfd674c0a2c1601e24fd7c80fb83c684c92b'
 oldc='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_PREPARED_CUSTODY_20261105_v63/report.json';oc=read(oldc);assert c['old_custody_sha256']==pins[oldc]=='0814d1eb45e09c34df1ee48a569376f9246ecabe0f73c4d90d4921cacee5ae7e'
 old=ROOT/'code/ambient_B18_ell2_T2_native70_continuation55_20261105_v63.py';new=ROOT/c['worker'];assert sha(old)==oc['worker_sha256'] and sha(new)==c['worker_sha256']=='d3989c58dc69d4628c75085aa578442571775e52abbc816893277da1a8a6d55e'
 a=old.read_text();b=new.read_text();fixes=[
 ("STEM='AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_20261105_v63'","STEM='AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_V2_20261105_v63'"),
 (" peer=read(peer_path);assert peer['all_passed'] and pins[gp] in json.dumps(peer)"," peer=read(peer_path);assert peer['all_passed']\n oldpp='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63_REVIEW/report.json';oldpeer=read(oldpp)\n assert oldpeer['all_passed'] and pins[gp] in json.dumps(oldpeer)\n oldpin=peer['input_sha256'][oldpp];assert (oldpin['sha256'] if isinstance(oldpin,dict) else oldpin)==pins[oldpp]"),
 ("count([g,peer,cr,ce,custody,pins,cg])","count([g,peer,oldpeer,cr,ce,custody,pins,cg])"),
 ("pins,g,peer,cr,ce,custody,cd,cg,integer_cols","pins,g,peer,oldpeer,cr,ce,custody,cd,cg,integer_cols")]
 for x,y in fixes:assert a.count(x)==1;a=a.replace(x,y)
 assert a==b
 fail='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_20261105_v63/BEFORE_ARITHMETIC_FAILURE.json';failure=read(fail);assert failure
 for path,h in c['worker_and_helper_sha256'].items():assert sha(ROOT/path)==h
 for path,h in c['literal_external_source_sha256'].items():assert sha(ROOT/path)==h
 h=subprocess.run([sys.executable,str(new),'--help'],cwd=ROOT,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=10);assert h.returncode==0 and not h.stderr and hashlib.sha256(h.stdout).hexdigest()==c['safe_help_stdout_sha256']
 live=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_LIVE_BOUND_20261105_v63/report.json');lev=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_LIVE_BOUND_20261105_v63/explicit-final-input-and-worker-overlap-count.json')
 metadata=count([g,p,oldpeer,live,lev,c,pins,cont])+100;assert metadata<900
 ev={'only_fixed_metadata_substitutions':len(fixes),'old_worker_and_custody_exact_preserved':True,'failure_before_all_finite_G8_55_rank_arithmetic':True,'explicit_old_source_peer_to_original_gate_chain':True,'continuation_gate_still_pinned':True,'safe_help_literal_PASS':True,'runtime_metadata_conservative_bound':metadata,'unchanged_metadata_cap':900,'mathematical_source_and_caps_unchanged':True,'new55_or_rank_values_queried':False};ep=out/'independent-V2-reader-chain-custody-and-unchanged-math.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'metadata_only_no_arithmetic_restart':True,'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'metadata':metadata,'seconds':report['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
