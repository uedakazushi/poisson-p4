#!/usr/bin/env python3
"""One actual Cstar-filtered Ga9 curve along the unique weight1 normal.

The prior ambient grading bounds the ENTIRE curve by degree10. No claim
about the complete four-normal support is made by this geometric ansatz.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from classification_ga9_weighted_frame import build_frame
from classification_rank3_mixed_component_review import COLUMNS,PAIRS,ROWS,TRIPLES,QUADS
ROOT=Path(__file__).resolve().parents[1]
def main():
 frame=build_frame();started=frame['started'];K=frame['K'];D=frame['D'];xx=frame['xx'];P=frame['P'];rowindex={r:i for i,r in enumerate(ROWS)};colindex={c:i for i,c in enumerate(COLUMNS)}
 def entry(A,i,j):return A.get((i,j),{})if i<j else {ex:-co for ex,co in A.get((j,i),{}).items()}if i>j else{}
 def put(A,i,co):
  if co:
   ij,ex=COLUMNS[i];f=A.setdefault(ij,{});f[ex]=f.get(ex,K.zero)+co
   if not f[ex]:del f[ex]
 def bracket(A,B=None):
  terms=[(A,A)]if B is None else[(A,B),(B,A)];third={}
  for i,j,k in TRIPLES:
   f={}
   for AA,BB in terms:
    for u,v,w in[(i,j,k),(j,k,i),(k,i,j)]:
     for ell in range(5):
      for ex,co in entry(AA,u,ell).items():
       for ee,cc in entry(BB,v,w).items():
        if ee[ell]:
         power=tuple(ex[q]+ee[q]-(q==ell)for q in range(5));f[power]=f.get(power,K.zero)+co*cc*ee[ell]
   third[(i,j,k)]={ex:co for ex,co in f.items()if co}
  result=[K.zero for _ in range(350)]
  for quad in QUADS:
   for at,ell in enumerate(quad):
    for ex,co in third[quad[:at]+quad[at+1:]].items():
     power=tuple(ex[q]+(q==ell)for q in range(5));result[rowindex[(quad,power)]]+=(-1)**at*co
  return result
 dcols=[[(i,K.from_sympy(D[i,j]))for i in range(350)if D[i,j]!=0]for j in range(150)]
 def linear(A):
  result=[K.zero for _ in range(350)]
  for ij,f in A.items():
   for ex,co in f.items():
    for row,a in dcols[colindex[(ij,ex)]]:result[row]+=a*co
  return result
 W={};obstruction=None
 for degree in range(1,17):
  nonlin=[K.zero for _ in range(350)]
  for i in range(1,degree//2+1):
   j=degree-i
   if i in W and j in W:
    terms=bracket(W[i])if i==j else bracket(W[i],W[j]);nonlin=[a+b for a,b in zip(nonlin,terms)]
  Wd={}
  if degree<=10:
   if degree==1:
    for i,a in enumerate(frame['normals'][0]):put(Wd,i,K.from_sympy(a))
   if degree in frame['solvers']:
    cs,rs,inverse=frame['solvers'][degree]
    for j,col in enumerate(cs):put(Wd,col,-sum((nonlin[row]*inverse[j][k]for k,row in enumerate(rs)),K.zero))
   W[degree]=Wd;ll=linear(Wd);current=[a+b for a,b in zip(nonlin,ll)]
  else:current=nonlin
  assert all(not current[i]for i in frame['selected_rows'])
  nz=[i for i,f in enumerate(current)if f]
  print(json.dumps(dict(weight=degree,nonzero_curve_columns=sum(map(len,Wd.values())),nonzero_remaining_equations=len(nz),seconds=time.monotonic()-started)),flush=True)
  if nz:
   obstruction=dict(weight=degree,row=nz[0],coefficient=str(current[nz[0]]));break
 e=s.Symbol('epsilon');family=P.copy()
 for degree,WW in W.items():
  for ij,f in WW.items():
   for ex,co in f.items():family[ij]+=K.to_sympy(co)*e**degree*s.prod(z**power for z,power in zip(xx,ex))
 for i,j in PAIRS:family[j,i]=-family[i,j]
 encode=lambda A:{str(i)+','+str(j):str(s.factor(A[i,j]))for i,j in PAIRS if A[i,j]!=0}
 inputs=[Path(__file__),ROOT/'code/classification_ga9_weighted_frame.py',frame['prior_path']]
 result=dict(task_id='E_reduced_nonnormal',all_passed=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_bytes=Path(__file__).read_text(),prior_finite_degree_bound10=True,actual_lowest_weight1_normal_axis=True,
  complete_curve_satisfies_all350_Jacobi=obstruction is None,first_exact_obstruction=obstruction,
  solved_curve_coefficients_by_weight={str(d):{str(colindex[(ij,ex)]):str(co)for ij,f in WW.items()for ex,co in f.items()}for d,WW in W.items()},
  scope='One fixed weight1 filteredcurve, with allothertrue normalparameters zero in the prescribed genuinecut. A positive result is an exactfinite polynomial Poissonfamily; a negative result excludes this axis only. No fullsupport or knownmembership follows automatically.',seconds=time.monotonic()-started)
 if obstruction is None:
  result['constructed_polynomial']=encode(family)
  pf=[]
  for missing in range(5):
   i,j,k,l=[i for i in range(5)if i!=missing];pf.append((-1)**missing*(family[i,j]*family[k,l]-family[i,k]*family[j,l]+family[i,l]*family[j,k]))
  result['actual_anticanonical_quintic']=str(s.factor(sum(xx[i]*pf[i]for i in range(5))))
 out=ROOT/'reports/CLASSIFICATION_GA9_CONSTANT_FILTERED_CURVE_20261107.json'
 if out.exists():raise FileExistsError(str(out))
 out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(actual_Poisson_curve=obstruction is None,obstruction=obstruction,seconds=result['seconds'])))
if __name__=='__main__':main()
