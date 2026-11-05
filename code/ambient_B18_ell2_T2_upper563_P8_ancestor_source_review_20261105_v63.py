#!/usr/bin/env python3
"""Literal bracket-stream C0 ancestor audit, with no coefficient evaluation."""
import json,resource,signal
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_UPPER563_P8_ANCESTOR_SOURCE_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));pins={}
 def read(p):pins[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 r=read(f'reports/{BASE}/report.json');assert r['all_passed'];p=f'reports/{BASE}/literal-two-P8-composite-program-and-exact-boundary-catalog.json';d=read(p);assert pins[p]['sha256']==r['evidence_sha256'][Path(p).name]
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 want=dict(d['selected_original_C0_arithmetic_nodes']);want.update({i:n for i,kind,n in d['exact_literal_cache_boundary_nodes']});cp='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';text=(ROOT/cp).read_text();pins[cp]={'path':cp,'sha256':sha(ROOT/cp)};pos=text.index('[',text.index('"DAG_scalar_nodes"'))+1;depth=1;quoted=False;escape=False;idx=-1;start=None;found={}
 for k in range(pos,len(text)):
  c=text[k]
  if quoted:
   if escape:escape=False
   elif c=='\\':escape=True
   elif c=='"':quoted=False
   continue
  if c=='"':quoted=True
  elif c=='[':
   depth+=1
   if depth==2:start=k;idx+=1
  elif c==']':
   if depth==2 and idx in want:found[idx]=json.loads(text[start:k+1])
   depth-=1
   if idx>=max(want) and depth==1:break
 assert found==want;del text
 fp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';pins[fp]={'path':fp,'sha256':sha(ROOT/fp)};known={i for i,v in fields(ROOT/fp,['first640'])['first640']['all79_atom_values_full10']};assert {20114,20138,40892,40916}<=known
 assert found[327603]==['a',20114,327602] and found[327602]==['n',40892] and found[327627]==['a',20138,327626] and found[327626]==['n',40916]
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'independent-eight-literal-nodes-and-two-difference-functions.json';ep.write_text(json.dumps({'literal_eight_nodes':sorted(found.items()),'two_exact_regular_functions':['327603=20114-40892','327627=20138-40916'],'all_four_boundaries_literal_cached79_members':True,'no_inverse_or_coefficient_or_availability_claim':True},separators=(',',':'))+'\n');rr={'id':STEM,'all_passed':True,'literal_source_only_no_function_values':True,'primary_sha256':pins[f'reports/{BASE}/report.json']['sha256'],'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)}};(out/'report.json').write_text(json.dumps(rr,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
