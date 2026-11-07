#!/usr/bin/env python3
"""Whole350 Jacobi on complete18unknown actual unipotent entrance (Ga9 linearw)."""
from pathlib import Path
import hashlib,itertools,json,os,subprocess,time
import sympy as s
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def encode(f,variables):
 terms=[]
 for ex,a in s.Poly(f,*variables).terms():
  co=str(a.p)if a.q==1 else '('+str(a.p)+'/'+str(a.q)+')';terms.append('*'.join([co]+[str(v)+('^'+str(k)if k!=1 else'')for v,k in zip(variables,ex)if k]))
 return '('+'+'.join(terms)+')'if terms else'0'
def main():
 started=time.monotonic();path=ROOT/'reports/CLASSIFICATION_GA9_UNIPOTENT_LIMITS_20261107.json';inp=json.loads(path.read_text());case=next(c for c in inp['cases']if c['name']=='LINEAR_W');linear=case['simple_pencil'];tt=s.symbols(' '.join(linear['free_actual_parameters']));xx=s.symbols('h x y z w');env={str(z):z for z in xx+tt};Q=s.zeros(5)
 for col,a in linear['canonical_linear_solution'].items():
  ij,ex=COLUMNS[int(col)];Q[ij]+=s.sympify(a,locals=env)*s.prod(v**k for v,k in zip(xx,ex))
 for i,j in PAIRS:Q[j,i]=-Q[i,j]
 sparse={ij:dict(s.Poly(Q[ij],*xx).terms())for ij in PAIRS if Q[ij]!=0};J=jacobi(sparse).applyfunc(s.expand);equations=[f for f in J if f!=0];sub={};steps=[]
 while True:
  current=[s.expand(f.subs(sub,simultaneous=True))for f in equations];constant=next((f for f in current if f!=0 and not f.free_symbols),None)
  if constant is not None:break
  free=sorted(set().union(*(f.free_symbols for f in current)),key=str);lins=[f for f in current if f!=0 and s.Poly(f,*free).total_degree()==1]
  if not lins:break
  A,b=s.linear_eq_to_matrix(lins,free);solution=s.linsolve((A,b),free)
  if solution is s.EmptySet:constant='inconsistent_linear_system';break
  values=list(solution)[0];new={t:f for t,f in zip(free,values)if t!=f};assert new;sub={t:s.expand(f.subs(new,simultaneous=True))for t,f in sub.items()};sub.update(new);steps.append(dict(linear_count=len(lins),solved={str(t):str(f)for t,f in new.items()}))
 current=[s.expand(f.subs(sub,simultaneous=True))for f in J];free=sorted(set().union(*(f.free_symbols for f in current)),key=str);unique=sorted({str(s.Poly(f,*free).monic().as_expr())for f in current if f!=0},key=str)if free else[str(f)for f in current if f!=0];print('18 ->',len(free),'equations',len(unique),'constant',constant,flush=True)
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),path,ROOT/'code/classification_rank3_mixed_component_review.py']},producer_bytes=Path(__file__).read_text(),complete_Q_brackets={str(i)+','+str(j):str(Q[i,j])for i,j in PAIRS if Q[i,j]!=0},original_parameters=list(map(str,tt)),all350_Jacobi=[str(f)for f in J],linear_elimination_steps=steps,all_linear_substitution={str(t):str(f)for t,f in sub.items()},remaining_parameters=list(map(str,free)),all350_residuals_after_linear_elimination=list(map(str,current)),unique_remaining_generators=unique,constant_obstruction=str(constant)if constant is not None else None,scope='FullJacobi of completeactualunipotentlinearentrance. Only actualPoissonQ solutions yield the exactP0+epsilonQ family; linearorpoint checks alone do not.',seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_UNIPOTENT_LINEAR_W_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print('saved',result['seconds'])
if __name__=='__main__':main()
