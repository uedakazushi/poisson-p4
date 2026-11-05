#!/usr/bin/env python3
"""Formal20H3minor schedule: literal determinant signs are never applied twice."""
import argparse,hashlib,itertools,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_H3_SIGN_SCHEDULE_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 out.mkdir(parents=True,exist_ok=False);n=7;z={};one={(0,)*n:1}
 def mono(i,k=1):e=[0]*n;e[i]=k;return{tuple(e):1}
 def mul(a,b):
  ans={}
  for x,c in a.items():
   for y,d in b.items():
    e=tuple(i+j for i,j in zip(x,y));ans[e]=ans.get(e,0)+c*d
  return{x:c for x,c in ans.items()if c}
 def add(a,b,s=1):
  ans=a.copy()
  for x,c in b.items():ans[x]=ans.get(x,0)+s*c
  return{x:c for x,c in ans.items()if c}
 def det(m):
  ans={}
  for p in itertools.permutations(range(3)):
   t=one
   for i in range(3):t=mul(t,m[i][p[i]])
   s=-1 if sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))%2 else 1;ans=add(ans,t,s)
  return ans
 # Variables(g,h,h120,Delta,lambda505,lambda419,lambda428).
 g,h,hp,D,l1,l2,l3=[mono(i)for i in range(7)];ig=mono(0,-1);iD=mono(3,-1)
 H=[[mul(h,ig),z,z],[mul(hp,ig),z,z],[mul(l1,iD),mul(l2,iD),l3],[one,z,z],[z,one,z],[z,z,D]]
 expected={(0,2,4):mul(mul(h,ig),l3),(0,2,5):mul(mul(h,ig),l2),(0,4,5):mul(mul(h,ig),D),(1,2,4):mul(mul(hp,ig),l3),(1,2,5):mul(mul(hp,ig),l2),(1,4,5):mul(mul(hp,ig),D),(2,3,4):l3,(2,3,5):l2,(2,4,5):l1,(3,4,5):D}
 negative={(0,2,4),(1,2,4),(2,3,5)};actual={}
 for rr in itertools.combinations(range(6),3):
  v=det([H[i]for i in rr]);want=expected.get(rr,{})
  if rr in negative:want={e:-c for e,c in want.items()}
  assert v==want;actual[','.join(map(str,rr))]=[[list(e),c]for e,c in v.items()]
 schedule=[]
 for rr in expected:
  outer=-1 if sum(rr)%2 else 1
  if rr in[(0,2,4),(0,2,5),(1,2,4),(1,2,5)]:kind='actual_literal_H3_determinant';source=list(rr);mult=outer
  else:kind='frozen_cofactor';source={(0,4,5):'lambda415',(1,4,5):'lambda442',(2,3,4):'lambda428',(2,3,5):'lambda419',(2,4,5):'lambda505',(3,4,5):'Delta407'}[rr];mult=outer*(-1 if rr==(2,3,5)else 1)
  schedule.append(dict(selected_H3_lower6=list(rr),pure_remaining_lower6=[i for i in range(6)if i not in rr],coefficient_kind=kind,source=source,combined_multiplier=mult,outer_laplace_sign=outer))
 dump(out/'schedule.json',dict(H3_row_order=['A','B','D','C','E','F'],H3_column_order=[1,3,2],formal_variable_order=['g','h','h120','Delta','lambda505','lambda419','lambda428'],all20_formal_minors=actual,ten_scalar_terms=schedule,literal_minor_sign_is_intrinsic_NOT_extra_weight_sign=True))
 note=ROOT/f'research/notes/{ID}.md';gate=ROOT/'reports/S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63/report.json';report=dict(id=ID,all_passed=True,frozen=True,source_only=True,all20_Laurent_minor_coefficients_derived=True,ten_actual_signed_terms_explicit=True,original_guarded_gate_unchanged=True,clarification='ForfourcomputedliteralH3minors,theirminussignisalreadyinthedeterminant;onlyouterLaplacesignisapplied.MetadataH3_weight_signwasrelativeunsignedh*lambda/gheightrecipe.',input_gate_sha256=sha(gate),evidence_sha256={'schedule.json':sha(out/'schedule.json')},code_sha256=sha(Path(__file__)),note_sha256=sha(note));dump(out/'report.json',report);print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
