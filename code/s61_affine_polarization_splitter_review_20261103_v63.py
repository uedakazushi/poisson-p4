#!/usr/bin/env python3
"""Independent review of the global known-five polarization idempotent.

Universal four monomial templates audit every original coefficient. A separate
Hessian-array engine rebuilds all25 full134 midpoint matrices/determinants.
"""
import argparse,hashlib,json,subprocess,sys,time
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as s
from s61_affine_fiber_macaulay import load_fiber
ROOT=Path(__file__).resolve().parents[1];P=32003

def detmod(M):
 A=np.array(M,dtype=np.int64).copy()%P;n=A.shape[0];d=1
 for j in range(n):
  non=np.flatnonzero(A[j:,j])
  if not len(non):return 0
  k=j+int(non[0])
  if k!=j:A[[j,k]]=A[[k,j]];d=-d
  pivot=int(A[j,j]);d=d*pivot%P;A[j]=A[j]*pow(pivot,-1,P)%P
  for i in range(j+1,n):A[i,j:]=(A[i,j:]-A[i,j]*A[j,j:])%P
 return d%P

def cv(v):
 v=Fraction(v);assert v.denominator%P;return v.numerator*pow(v.denominator,-1,P)%P

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic();producer=ROOT/'code/s61_affine_polarization_splitter_20261103_v63.py';rp=ROOT/'reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63/report.json';r=json.loads(rp.read_text());assert r['all_passed']
 for f,h in r['input_sha256'].items():assert hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==h,f
 for f,h in r['evidence_sha256'].items():assert hashlib.sha256((rp.parent/f).read_bytes()).hexdigest()==h,f
 fresh=subprocess.run([sys.executable,str(producer),'--output',str(out/'producer-replay')],cwd=ROOT,check=True,capture_output=True,text=True);assert 'PASS globaldegree45 knownfive idempotent' in fresh.stdout
 source=ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json';data=json.loads(source.read_text());weights=data['parameter_weights'];assert weights.count(1)==16 and weights.count(2)==45
 # Exact symbolic template F(q)-F(p)=J((q+p)/2)(q-p), retaining any
 # beta-monomial as an untouched scalar coefficient. Every source term below
 # belongs to one of these FOUR templates, so this checks all its coefficients.
 q0,q1,p0,p1=s.symbols('q0 q1 p0 p1');templates=[1,q0,q0*q0,q0*q1]
 for f in templates:
  lhs=f-s.sympify(f).subs({q0:p0,q1:p1},simultaneous=True)
  rhs=sum(s.diff(f,v).subs({q0:(q0+p0)/2,q1:(q1+p1)/2},simultaneous=True)*(v-pv)for v,pv in [(q0,p0),(q1,p1)])
  assert s.expand(lhs-rhs)==0
 counts=[]
 for row in data['equations']:
  classes=Counter()
  for mon,n,d in row:
   qmon=[i for i in mon if weights[i]==2];beta=[i for i in mon if weights[i]==1];assert len(qmon)+len(beta)==len(mon) and len(qmon)<=2 and d%P
   classes['constant'if not qmon else 'linear'if len(qmon)==1 else 'square'if qmon[0]==qmon[1]else'mixed']+=1
  counts.append(dict(classes))
 assert len(counts)==134
 F,qids,center,*_=load_fiber();L=np.zeros((134,45),dtype=np.int64);H=np.zeros((134,45,45),dtype=np.int64)
 for i,f in enumerate(F):
  for mon,co in f.items():
   value=cv(co)
   if len(mon)==1:L[i,mon[0]]=(L[i,mon[0]]+value)%P
   elif len(mon)==2:
    u,v=mon;H[i,u,v]=(H[i,u,v]+value)%P;H[i,v,u]=(H[i,v,u]+value)%P
   else:assert co==0
 fivepath=ROOT/'reports/S61_KNOWN_FIVE_POINT_SECTION_20261101.json';five=json.loads(fivepath.read_text());vs=s.symbols('v0:3');slopes=np.array([[cv(s.sympify(f).coeff(v))for v in vs]for f in five['affine_plane_substitution']],dtype=np.int64)
 knownpath=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json';points=json.loads(knownpath.read_text())['known_five_mod_32003'];qpoints=[slopes@np.array(p[1:],dtype=np.int64)%P for p in points]
 separator=rp.parent/'known-five-separator.json';cert=json.loads(separator.read_text());rows=cert['selected_original_rows_zero_based'];assert len(set(rows))==45
 values=[]
 for p in qpoints:
  row=[]
  for q in qpoints:
   middle=(p+q)*pow(2,-1,P)%P;J=(L+np.einsum('ijk,k->ij',H,middle))%P;assert not np.any(J@((q-p)%P)%P);row.append(detmod(J[rows,:]))
  values.append(row)
 assert values==cert['midpoint25_determinants'];diag=[values[i][i]for i in range(5)];assert all(diag) and all(values[i][j]==0 for i in range(5)for j in range(5)if i!=j)
 hc=np.arange(1,46,dtype=np.int64);offset=np.array(list(map(cv,center)),dtype=np.int64);hvalues=[int(hc@((p+offset)%P))%P for p in qpoints];assert hvalues==cert['separating_h_values'] and len(set(hvalues))==5
 norm=1;vand=1
 for d in diag:norm=norm*d%P
 for i in range(5):
  for j in range(i+1,5):vand=vand*(hvalues[j]-hvalues[i])%P
 assert norm==cert['Jacobian_norm_mod_prime'] and vand==cert['h_Vandermonde_mod_prime'] and vand*vand%P==cert['h_discriminant_mod_prime']
 dependencies=[producer,rp,source,fivepath,knownpath,separator,ROOT/'research/notes/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63.md',Path(__file__)]
 result=dict(id='S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63_REVIEW',all_passed=True,producer_fresh_replay=True,universal_134_polarization_accepted_by_exact_monomial_templates=True,all16_beta_preserved_as_scalar_coefficients=True,original_term_classes=counts,independent_Hessian_134x45x45_engine=True,all25_full134_midpoint_kernel_identities=True,all25_45determinants=values,knownfive_derivative_norm=norm,separating_h_values=hvalues,separating_discriminant=vand*vand%P,generic_adjugate_trace_proof_accepted=True,known_quotient_source_audited_as_existing_input_not_rebuilt=True,accepted_result='On a beta-only nonempty open, GLOBAL actual affine A is E_known5 x A/(e), degree_q e<=45; no q/H/cyclic localization',unit_convention='E -> eA is unital for unit e; its inclusion into A is not unit-preserving',residual_empty_or_reduced_or_finite_claimed=False,input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in dependencies},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('GLOBAL POLARIZATION SPLITTER INDEPENDENT REVIEW PASS',result['seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
