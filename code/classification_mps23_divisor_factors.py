#!/usr/bin/env python3
"""Exact generic divisor factors for every projective-dimension23 MPS family.

Uses whole accepted rational families, never a single specialization or just a graph.
"""
from pathlib import Path
import hashlib,json,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();xx=s.symbols('x0:5');data_path=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json';prior_path=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_PRIOR98_20261107.json';data=json.loads(data_path.read_text());prior=json.loads(prior_path.read_text());wanted={c['id']for c in prior['components']if c['id'].startswith('T')and c['projective_dimension']==23};assert len(wanted)==29
 cases=[];inputs=[data_path,prior_path,Path(__file__)]
 for record in data['diagrams']:
  cid=record['id']
  if cid not in wanted:continue
  path=ROOT/f'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/{cid}.json'
  if not path.exists():path=ROOT/f'reports/MPS_S61_RAW_CATALOGUE_20261103_v62/{cid}.json'
  inp=json.loads(path.read_text());inputs.append(path);assert inp.get('exact_projective_Jacobi_checked',inp.get('projective_Jacobi_identically_zero',False));pp=s.symbols('p0:'+str(record['period_dimension']));tt=s.symbols('t0:'+str(record['edge_count']));env={str(z):z for z in xx+pp+tt};P=s.zeros(5)
  for ij,f in inp['homogeneous_bivector'].items():
   i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
  pf=[]
  for missing in range(5):
   i,j,k,l=[a for a in range(5)if a!=missing];pf.append((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
  D=s.cancel(sum(xx[i]*pf[i]for i in range(5)));co,factors=s.factor_list(D,*xx);factor_data=[dict(factor=str(f),multiplicity=m,coordinate_degree=s.Poly(f,*xx).total_degree())for f,m in factors];assert sum(z['multiplicity']*z['coordinate_degree']for z in factor_data)==5;assert len(factor_data)>1;assert s.cancel(D-co*s.prod(f**m for f,m in factors))==0
  isolated=[i for i in range(5)if not any(i in e['edge']for e in record['edges'])]
  for i in isolated:assert s.cancel(D.subs(xx[i],0))==0
  case=dict(id=cid,input_path=str(path.relative_to(ROOT)),input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),all_general_parameters=list(map(str,pp+tt)),actual_quintic=str(D),factor_coefficient=str(co),factors=factor_data,irreducible_coordinate_factor_degrees=[s.Poly(f,*xx).total_degree()for f,m in factors for _ in range(m)],exact_product_identity=True,isolated_hyperplanes=isolated,generic_source_factorization=True)
  cases.append(case);print(cid,case['irreducible_coordinate_factor_degrees'],round(time.monotonic()-started,3),flush=True)
 assert {c['id']for c in cases}==wanted
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,scope='All29 accepted MPS dimension23 generic rational sources have reducible actual quintic. Exact whole-parameter factor product identities imply each full component closure lies in the closed reducible-quintic locus. No inference from SNC existence or a specialized point.',component_count=29,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),cases=cases,seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_MPS23_DIVISOR_FACTORS_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print('PASS',result['seconds'])
if __name__=='__main__':main()
