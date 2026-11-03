#!/usr/bin/env python3
"""Independently replay every J0 coefficient using the binomial formula."""
import argparse,hashlib,json,time
from fractions import Fraction as F
from math import comb
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def parse(line):
 z=line.split();n=int(z[0]);assert len(z)==1+4*n
 return [(F(z[1+4*k]),tuple(map(int,z[2+4*k:5+4*k])))for k in range(n)]
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',required=True,type=Path);ap.add_argument('--prepared',required=True,type=Path);ap.add_argument('--row-scalings',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();t=time.monotonic()
 scalings=[list(map(int,line.split()))for line in args.row_scalings.read_text().splitlines()];assert len(scalings)==36
 count=0;termcount=0;factor={}
 with args.source.open()as src,args.prepared.open()as dst:
  assert list(map(int,src.readline().split()))==[36,70,21];assert len(src.readline().split())==21;assert len(src.readline().split())==36
  assert list(map(int,dst.readline().split()))==[36,70,0]
  for i in range(36):
   row=[parse(src.readline())for _ in range(70)];M=max([0]+[d-a for entry in row for c,(a,e,d)in entry]);maxd=max([0]+[d for entry in row for c,(a,e,d)in entry])
   assert scalings[i][:3]==[i,M,maxd];terms=0
   for entry in row:
    actual={ex:c for c,ex in parse(dst.readline())};expected={}
    for c,(a,e,d)in entry:
     assert a+M-d>=0
     if d not in factor:factor[d]=[F(comb(d,j)*6**j*4**(d-j),5**d)for j in range(d+1)]
     for j,f in enumerate(factor[d]):
      ex=(a+M-d,e+j,0);v=expected.get(ex,F(0))+c*f
      if v:expected[ex]=v
      else:expected.pop(ex,None)
    assert actual==expected,(i,count%70);count+=1;terms+=len(expected)
   assert scalings[i][3]==terms;termcount+=terms;print('QQ_J0_ALL_ENTRY_ROW_PASS',i,flush=True)
 assert count==2520 and termcount==108530
 paths=[Path(__file__),args.source,args.prepared,args.row_scalings,ROOT/'code/s61_nonsquare_J0_pid_prepare_20261103_v63.cpp',ROOT/'work/nonsquare-v62/cubic27-v6/phase21-check.json',ROOT/'research/handoff/20261103_v62/v6_cubic27_coker/manifest.json']
 report={'id':'S61_NONSQUARE_J0_PID_SUBSTITUTION_20261103_V63','all_passed':True,'coefficient_field':'QQ',
  'shape':[36,70],'exact_entry_count':2520,'exact_term_count':termcount,'d_substitution':'(6e+4)/(5a)',
  'row_scaling_only_known_unit':'a^M','independent_method':'Full binomial expansion coefficient-by-coefficient, no point evaluation',
  'J_never_inverted':True,'safe_source_pivots':21,'source_known_units_only':['a','F','B','C'],
  'generic_finite_field_rank_used_as_proof':False,'seconds':time.monotonic()-t,
  'input_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}}
 args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
