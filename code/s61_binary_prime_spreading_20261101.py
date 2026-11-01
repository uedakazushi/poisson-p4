#!/usr/bin/env python3
"""Check good reduction of the reviewed universal binary65 lower bound."""
import argparse,hashlib,time
import sys,json,sympy as s,numpy as np
from pathlib import Path
from sympy.polys.rings import ring
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def verify():
 start=time.time()
 p=32003;m=S61();root=ROOT;data=json.load((ROOT/'reports/S61_BINARY_QUADRATIC_PROJECTION_DATA_20261101.json').open());lift=json.load((ROOT/'reports/S61_BINARY_CUBIC_LIFT_DATA_20261101.json').open());cr=json.load((ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json').open());rr=json.load((ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json').open());qr=json.load((ROOT/'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json').open());names=data['variables'];R,*xx=ring(','.join(names),s.QQ);FF=R.to_field();sym=s.symbols(' '.join(names));env=dict(zip(names,sym));cache={}
 def audit(v):
  if v in cache:return cache[v]
  f=FF.from_expr(s.sympify(v,locals=env))
  for g in [f.numer,f.denom]:
   for c in g.values():assert int(c.denominator)%p
  assert any(int(c.numerator)%p for c in f.denom.values())
  cache[v]=f;return f
 for matrix in [data['projection'],data['tangent'],lift['annihilator'],cr['selected_Bianchi_multipliers']]:
  for row in matrix:
   for v in row:audit(v)
 encoded=list(m.reduced['equations'])+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
 for equation in encoded:
  for mon,num,den in equation:assert int(den)%p
 mods={}
 for key,v in [('slice_dominance',rr['slice_tangent_minor']),('leading_rank',rr['leading_minor_at_sample']),('full_rank',rr['full_minor_at_sample']),('quadratic_rank',qr['sample_minor']),('parameter_radical',cr['parameter_projection_minor']),('bianchi_rank',cr['Bianchi_minor'])]:
  q=s.Rational(v);assert q.q%p;mods[key]=int(q.p)*pow(int(q.q),-1,p)%p;assert mods[key]
 u=s.symbols('u1:62');one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];sample={u[i]:s.Rational(j+2) for j,i in enumerate(one)};sample.update({s.Symbol(k):s.Rational(v)for k,v in data['fixed_beta'].items()});sample[u[11]]=sample[u[12]]=0
 mod=lambda v:int(s.Rational(v).p)*pow(int(s.Rational(v).q),-1,p)%p
 def rank(M):
  M=np.array([[mod(v)for v in row]for row in s.Matrix(M).tolist()],dtype=np.int64);r=0
  for col in range(M.shape[1]):
   ii=np.flatnonzero(M[r:,col])
   if not len(ii):continue
   j=r+int(ii[0]);M[[r,j]]=M[[j,r]];M[r]=M[r]*pow(int(M[r,col]),-1,p)%p;ids=np.flatnonzero(M[:,col]);ids=ids[ids!=r];M[ids]=(M[ids]-M[ids,col,None]*M[r,None,:])%p;r+=1
   if r==len(M):break
  return r
 L=s.Matrix(lift['annihilator']).applyfunc(s.sympify).subs(sample);assert rank(L)==16;two,F=specialize(m,{i:sample[u[i]]for i in one},True);fi=data['chart_fixed_q'];cols=[i for i in range(45)if i!=fi];J=s.zeros(135,45)
 for row,f in enumerate(F):
  for mon,c in f.items():
   if len(mon)==2:
    a,b=mon
    if a==fi and b!=fi:J[row,cols.index(b)]+=c
    if b==fi and a!=fi:J[row,cols.index(a)]+=c
   elif len(mon)==1 and mon[0]==fi:J[row,44]+=c
 assert rank(J)==35;assert rank(J[data['jacobian_rows'],lift['support_columns']])==19
 for M in [data['projection'],data['tangent']]:
  MM=s.Matrix(M).applyfunc(s.sympify).subs(sample)
  for q in MM:assert s.Rational(q).q%p
 out=dict(prime=p,all_passed=True,audited_distinct_rational_functions=len(cache),equation_denominators_checked=135,witness_minors_mod_prime=mods,annihilator_rank=16,supported_jacobian_rank=19,full_jacobian_rank=35)
 paths=[Path(__file__)]+[ROOT/"reports"/(name+".json") for name in ["S61_BINARY_QUADRATIC_PROJECTION_DATA_20261101","S61_BINARY_CUBIC_LIFT_DATA_20261101","S61_BINARY_UNIFORM_CUBIC_20261101","S61_BINARY_UNIFORM_RANK_20261101","S61_BINARY_UNIFORM_QUADRATIC_20261101"]]+[ROOT/"code/s61_model.py",ROOT/"code/s61_flag_transverse_obstruction_20261101.py",ROOT/"code/s61_binary_uniform_cubic_20261101.py",ROOT/"reports/S61_UNIVERSAL_WEIGHT6_20261031.json"]
 out.update(id="S61_BINARY_PRIME_SPREADING_20261101",scope="The existing universal binary normal-length lower bound65 spreads to characteristic32003 on its stated constant-rank smooth locus.",proof_note="research/notes/S61_BINARY_PRIME_SPREADING_20261101.md",input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},seconds=time.time()-start)
 return out

if __name__=="__main__":
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--output",required=True,type=Path);args=parser.parse_args();r=verify()
 with args.output.open("x") as f:json.dump(r,f,indent=2);f.write("\n")
 print(json.dumps(dict(all_passed=r["all_passed"],prime=r["prime"],seconds=r["seconds"])))
