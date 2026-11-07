#!/usr/bin/env python3
"""Actual cubicPoisson curve on a rational line of the COMPLETE Ga9 cut."""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,jacobi
ROOT=Path(__file__).resolve().parents[1]
def main():
 started=time.monotonic();frame_path=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_WEIGHT_BLOCKS_20261107.json';support_path=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_ANCHOR_SUPPORT_20261107.json';frame=json.loads(frame_path.read_text());packet=json.loads(support_path.read_text());xx=s.symbols('h x y z w');tt=s.symbols('t0:6');e=s.Symbol('epsilon');env={str(z):z for z in xx+tt};line={tt[0]:0,tt[1]:-3*e,tt[2]:e,tt[3]:0,tt[4]:0,tt[5]:0};P=s.zeros(5)
 for ij,f in frame['canonical_point'].items():
  i,j=map(int,ij.split(','));P[i,j]=s.sympify(f,locals=env);P[j,i]=-P[i,j]
 original=P.copy();v=s.zeros(150,1)
 for j,normal in enumerate(frame['true_normal_vectors']):
  for col,a in normal.items():v[int(col)]+=s.Rational(a)*line[tt[j]]
 for item in packet['complete_weight3_implicit_polynomials']:
  assert item['denominator_delta_power']>=1
  value=s.sympify(item['numerator'].replace('^','**'),locals=env).subs(line,simultaneous=True);v[item['column']]+=s.expand(value)
 for (ij,ex),co in zip(COLUMNS,v):P[ij]+=co*s.prod(z**a for z,a in zip(xx,ex))
 for i,j in PAIRS:P[j,i]=-P[i,j]
 P=P.applyfunc(s.expand);assert P.subs(e,0)==original
 for item in packet['complete_all350_remaining_equations_through_weight4']:assert s.expand(s.sympify(item['numerator'].replace('^','**'),locals=env).subs(line,simultaneous=True))==0
 sparse={ij:dict(s.Poly(P[ij],*xx).terms())for ij in PAIRS if P[ij]!=0};assert jacobi(sparse).applyfunc(s.expand)==s.zeros(350,1)
 weights=frame['weights'];c=frame['canonical_tensor_weight']
 for ij in PAIRS:
  for ex,co in s.Poly(P[ij],*xx).terms():
   relative=sum(a*w for a,w in zip(ex,weights))-sum(weights[i]for i in ij)-c
   assert s.expand(co-e**relative*co.subs(e,1))==0
 pf=[]
 for missing in range(5):
  i,j,k,l=[a for a in range(5)if a!=missing];pf.append((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k]))
 D=s.factor(sum(xx[i]*pf[i]for i in range(5)));point=P.subs(e,1);D1=s.factor(D.subs(e,1));factorization=s.factor_list(D1);print('D1',D1,flush=True)
 encode=lambda A:{str(i)+','+str(j):str(s.factor(A[i,j]))for i,j in PAIRS if A[i,j]!=0}
 recognition=[]
 for f,m in factorization[1]:
  if s.Poly(f,*xx).total_degree()==1:
   row=s.Matrix([[s.diff(f,x)for x in xx]]);pivot=next(i for i,a in enumerate(row)if a!=0);T=s.zeros(5);T[0,:]=row;order=[i for i in range(5)if i!=pivot]
   for j,i in enumerate(order,1):T[j,i]=1
   Q=(T*point*T.T).subs(dict(zip(xx,T.inv()*s.Matrix(xx))),simultaneous=True).applyfunc(s.expand);R=s.Matrix([0]+[s.cancel(Q[0,j]/xx[0])for j in range(1,5)]);assert all(not s.denom(a).free_symbols.intersection(xx)for a in R)
   Q=(Q-s.Matrix(xx)*R.T+R*s.Matrix(xx).T).applyfunc(s.expand);aff=Q.subs(xx[0],1);N=s.Matrix([s.expand(sum(s.diff(aff[i,j],xx[i])for i in range(1,5)))for j in range(1,5)]);A=N.jacobian(xx[1:]);cc=N.subs(dict.fromkeys(xx[1:],0))
   print('H',f,'N',N,'char',s.factor(A.charpoly().as_expr()),flush=True)
   recognition.append(dict(linear_factor=str(f),multiplicity=m,projective_matrix=[[str(a)for a in rr]for rr in T.tolist()],canonical_H_Casimir_brackets=encode(Q),modular=list(map(str,N)),linear_matrix=[[str(a)for a in rr]for rr in A.tolist()],constant=list(map(str,cc)),characteristic_polynomial=str(s.factor(A.charpoly().as_expr()))))
 inputs=[frame_path,support_path,Path(__file__),ROOT/'code/classification_rank3_mixed_component_review.py'];result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),
  complete_cut_rational_line=['0','-3*epsilon','epsilon','0','0','0'],all162_exact_cut_equations_vanish=True,all350_projective_Jacobi_identically_zero=True,original_source_anchor_at0=True,whole_actual_cubic_Poisson_curve=encode(P),nonzero_point_brackets=encode(point),actual_anticanonical_quintic=str(D),nonzero_point_actual_quintic=str(D1),nonzero_point_factorization=[dict(factor=str(f),multiplicity=m,degree=s.Poly(f,*xx).total_degree())for f,m in factorization[1]],
  all_nonzero_parameters_one_actual_PGL_scalar_orbit=True,equivariance='For g(t)=diag(t^0,t^1,t^0,t^1,t^1), t^(-1)*g(t)_*P_e=P_(e/t), so everye!=0 is equivalentto e=1.',simple_H_recognition=recognition,
  known_component_membership_claimed=False,scope='Exactactualcurve and cofactorrecognitiondata. Apply an acceptedwholeHquartic theorem with allhypotheses before claiming membership. The sourceorbit22 is dense in the entire5parameter sourceclosure by the separately accepted sourceupper22 and actualorbit22.',seconds=time.monotonic()-started)
 dest=ROOT/'reports/CLASSIFICATION_GA9_DISTINCT_CURVE_20261107.json'
 if dest.exists():raise FileExistsError(str(dest))
 dest.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(actual_curve=True,seconds=result['seconds'])))
if __name__=='__main__':main()
