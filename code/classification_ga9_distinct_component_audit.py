#!/usr/bin/env python3
"""Frozen-prior98 distinctness audit and exact stabilizer eigen-coordinate model."""
from pathlib import Path
import hashlib,json,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();names=['reports/CLASSIFICATION_GA9_DISTINCT_PRIOR98_20261107.json','reports/CLASSIFICATION_GA9_DISTINCT_CURVE_ORBIT_20261107.json','reports/CLASSIFICATION_MPS23_DIVISOR_FACTORS_20261107.json'];prior,orb,mps=[json.loads((ROOT/p).read_text())for p in names];assert prior['component_count']==98;assert orb['full_D_rank']==102 and orb['actual_group_Euler_scalar_rank']==48;mf={d['id']:d for d in mps['cases']};audit=[]
 for c in prior['components']:
  if c['projective_dimension']!=23:reason='Unequal full-component dimension';data=dict(prior_projective_dimension=c['projective_dimension'])
  elif c['id']in mf:reason='Whole accepted generic MPS source exact reducible quintic product; closed factorization locus contains whole component';data=dict(factor_report=names[2],degrees=mf[c['id']]['irreducible_coordinate_factor_degrees'])
  else:
   D=c['generic_actual_degeneracy_divisor'];degs=D['irreducible_factor_degrees'];mult=D['corresponding_multiplicities'];assert sum(a*b for a,b in zip(degs,mult))==5 and sum(mult)>1;reason='Accepted dense source has reducible quintic; closed factorization locus contains whole component';data=dict(degrees=degs,multiplicities=mult,proof_sources=c['proof_sources'])
  audit.append(dict(id=c['id'],reason=reason,**data))
 xx=s.symbols('h x y z w');uu=s.symbols('u0:5');env=dict(zip(map(str,xx),xx));A=s.Matrix(orb['full_projective_stabilizer_generators'][0]['matrix']);vectors=[];eigenvalues=[]
 for ev,m,vv in sorted(A.eigenvects(),key=lambda a:a[0]):
  assert m==len(vv)==1;v=vv[0];den=s.ilcm(1,*[a.q for a in v]);v*=den;g=s.igcd(*[abs(int(a))for a in v if a]);v/=g;vectors.append(v);eigenvalues.append(int(ev))
 B=s.Matrix.hstack(*vectors);sub=dict(zip(xx,B*s.Matrix(uu)));D=s.Poly(s.sympify(orb['actual_quintic_primitive'],locals=env).subs(sub,simultaneous=True),*uu);content,D=D.primitive();assert D.total_degree()==5
 P=s.zeros(5)
 for ij,f in orb['canonical_divergence_free_point'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
 P=(B.inv()*P*B.inv().T).subs(sub,simultaneous=True).applyfunc(s.expand);den=s.ilcm(1,*[a.q for f in P for a in s.Poly(f,*uu).coeffs()]);P*=den;g=s.igcd(*[abs(int(a))for f in P for a in s.Poly(f,*uu).coeffs()if a]);P/=g
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,scope='Distinctness of the smooth23 actualcurve component from every frozenprior98, CONDITIONAL only on separately checked absolute irreducibility. Exact eigencoordinates are optional readability evidence, not a new point.',input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in names+[str(Path(__file__).relative_to(ROOT))]},producer_bytes=Path(__file__).read_text(),candidate_dimension=23,candidate_generic_rank=4,all98_audit=audit,same23_MPS_count=len(mf),same23_nonMPS_count=sum(c['projective_dimension']==23 and c['id']not in mf for c in prior['components']),old_coordinates_from_eigen_coordinates=[[str(a)for a in rr]for rr in B.tolist()],eigenvalues=eigenvalues,canonical_eigen_coordinate_quintic=str(D.as_expr()),canonical_eigen_coordinate_bivector={str(i)+','+str(j):str(P[i,j])for i in range(5)for j in range(i+1,5)if P[i,j]!=0},eigen_quintic_original_scalar=str(content),bivector_original_scaling=str(den/g),seconds=time.monotonic()-started)
 out=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_COMPONENT_AUDIT_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k]for k in ['same23_MPS_count','same23_nonMPS_count','canonical_eigen_coordinate_quintic','canonical_eigen_coordinate_bivector','seconds']}))
if __name__=='__main__':main()
