#!/usr/bin/env python3
"""The one-parameter second-order vertex flag for the short 2/3 source.

Complete parameter treatment is via exact rational-function elimination and
all roots of the residual gcd or pivot denominator, never a point sweep.
"""
from pathlib import Path
import hashlib,itertools,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from classification_nil31_short_flag_confluence import PAIRS,wedge,canonical
ROOT=Path(__file__).resolve().parents[1]
def main():
 begun=time.monotonic();xx=s.symbols('h x y z w');h,x,y,z,w=xx;e,beta=s.symbols('epsilon beta');ss=s.symbols('hs X Y Z W');hs,X,Y,Z,W=ss
 D1=s.Matrix([0,X,0,0,W]);D2=s.Matrix([0,0,Y,Z,0]);N=s.Matrix([0,0,0,Y,-X]);bases=[wedge(D1,D2),(X*Z+Y*W)*wedge(s.eye(5)[:,3],s.eye(5)[:,4]),wedge(N,D2),wedge(D1,N)]
 forms=s.Matrix([h+e*y,h,z+beta*e*y,y+e*x,-w]);L=forms.jacobian(s.Matrix(xx));Li=L.inv();assert s.factor(L.det()) in [e*e,-e*e]
 pulled=[canonical((Li*P.subs(dict(zip(ss,forms)),simultaneous=True)*Li.T).applyfunc(s.expand),xx)for P in bases]
 target=s.zeros(5);br=[h*z-h*y+2*z*w-y*y/3,h*z-2*y*z/3,h*h+y*w/3,-2*z*z/3,(h*y+z*w)/3,2*h*z/3]
 for (i,j),a in zip(itertools.combinations(range(1,5),2),br):target[i,j]=a;target[j,i]=-a
 target=canonical(target,xx);cs=s.symbols('c0:20');periods=[sum(cs[5*i+k]*e**k for k in range(5))for i in range(4)];candidate=sum((a*P for a,P in zip(periods,pulled)),s.zeros(5));rows=[];rhs=[];labels=[]
 for i,j in PAIRS:
  for (k,*ex),coef in s.Poly(s.expand(e**4*(candidate-target)[i,j]),e,*xx).terms():
   if k<=4:
    M,v=s.linear_eq_to_matrix([coef],cs);rows.append(list(M));rhs.append(v[0]);labels.append(dict(pair=[i,j],epsilon_power=k-4,exponents=ex))
 A=s.Matrix(rows);b=s.Matrix(rhs);field=s.QQ.frac_field(beta);dm=lambda M:DomainMatrix.from_Matrix(M).convert_to(field)
 piv=list(dm(A).rref()[1]);rr=list(dm(A[:,piv].T).rref()[1]);block=A[rr,piv];det=s.factor(dm(block).det().as_expr());inv=dm(block).inv().to_Matrix();res=(b-A[:,piv]*inv*b[rr,:]).applyfunc(s.cancel);nonzero=[a for a in res if a!=0];gcd=s.factor(s.gcd_list([s.fraction(a)[0]for a in nonzero]))if nonzero else s.Integer(0)
 # Every beta outside zero(det) is possible only if it zeros all numerators.
 exceptional=s.factor(det*s.prod(s.denom(a)for a in inv));exceptional=s.fraction(exceptional)[0];rootpoly=s.Poly(s.prod(f for f,m in s.factor_list(exceptional* (gcd if gcd!=0 else 1),beta)[1]),beta);roots=s.polys.polytools.ground_roots(rootpoly)
 unresolved=s.div(rootpoly,s.prod(beta-a for a in roots),beta)[0];remaining=[f for f,m in s.factor_list(unresolved,beta)[1]if s.Poly(f,beta).degree()>0]
 candidates=[];positive=[]
 for value in sorted(roots,key=str):
  Av=A.subs(beta,value);bv=b.subs(beta,value);dv=lambda M:DomainMatrix.from_Matrix(M).convert_to(s.QQ);pc=list(dv(Av).rref()[1]);pr=list(dv(Av[:,pc].T).rref()[1]);bi=dv(Av[pr,pc]).inv().to_Matrix();rv=(bv-Av[:,pc]*bi*bv[pr,:]).applyfunc(s.cancel);bad=next((i for i,v in enumerate(rv)if v!=0),None);record=dict(beta=str(value),rank=len(pc),rank_augmented=len(pc)+(bad is not None),consistent=bad is None)
  if bad is not None:
   lam=s.zeros(1,Av.rows);lam[bad]=1
   for j,k in enumerate(pr):lam[k]-=(Av[bad,pc]*bi)[j]
   assert lam*Av==s.zeros(1,20) and(lam*bv)[0]!=0;record.update(left_functional={str(j):str(a)for j,a in enumerate(lam)if a},left_rhs=str((lam*bv)[0]))
  else:
   solution=s.zeros(20,1)
   for j,k in enumerate(pc):solution[k]=(bi*bv[pr,:])[j]
   bind=dict(zip(cs,solution));fam=candidate.subs(beta,value).subs(bind).applyfunc(s.cancel);assert all(s.denom(a).subs(e,0)!=0 for a in fam);assert(fam.subs(e,0)-target).applyfunc(s.cancel)==s.zeros(5)
   triples={inds:s.expand(sum(fam[i,l]*s.diff(fam[j,k],xx[l])+fam[j,l]*s.diff(fam[k,i],xx[l])+fam[k,l]*s.diff(fam[i,j],xx[l])for l in range(5)))for inds in itertools.combinations(range(5),3)for i,j,k in [inds]};assert all(s.expand(sum((-1)**k*xx[inds[k]]*triples[inds[:k]+inds[k+1:]]for k in range(4)))==0 for inds in itertools.combinations(range(5),4))
   record.update(actual_source_periods=[str(s.expand(p.subs(bind)))for p in periods],actual_family={str(i)+','+str(j):str(fam[i,j])for i,j in PAIRS if fam[i,j]!=0},all350_Jacobi=True,endpoint_exact=True);positive.append(str(value))
  candidates.append(record)
 paths=[Path(__file__),ROOT/'code/classification_nil31_short_flag_confluence.py',ROOT/'reports/CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_CANONICAL_20261107.json',ROOT/'research/classification/THREE_HYPERPLANES_QUADRIC_INDEPENDENT.md',ROOT/'reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_REVIEW_20261107.json']
 out=dict(task_id='E_multiple_hyperplane',all_passed=True,status='ACTUAL_CONFLUENCE_FOUND'if positive else('FIXED_CHAIN_FLAG_CLASS_EXCLUDED'if nonzero and not remaining else'PARTIAL_PARAMETER_ELIMINATION_ONLY'),source_forms=list(map(str,forms)),source_inverse=[[str(a)for a in r]for r in Li.tolist()],source_frame_determinant=str(s.factor(L.det())),geometric_diagnosis='rank4 quadric vertex moves as e_x-epsilon*e_y+beta*epsilon²*e_z, incorporating the length3 Jordan chain absent in the first flag.',full_linear_matrix=[[str(a)for a in r]for r in A.tolist()],full_rhs=list(map(str,b)),row_labels=labels,matrix_shape=list(A.shape),generic_rank=len(piv),generic_rank_augmented=len(piv)+(bool(nonzero)),minor_rows=rr,minor_columns=piv,minor_determinant=str(det),generic_nonzero_residuals=list(map(str,nonzero)),generic_residual_numerator_gcd=str(gcd),all_candidate_parameter_polynomial=str(rootpoly.as_expr()),rational_candidates=candidates,remaining_algebraic_parameter_factors=list(map(str,remaining)),positive_parameters=positive,whole_source_membership='Candidate actual curve requires independent review'if positive else'UNKNOWN: only this single geometric flag class was examined. No arbitrary arc or component exclusion.',bound='Only four genuine source coefficients, polynomial epsilon degree≤4. Inverse pole2/wedge pole4. No stronger completeness is used for a positive curve.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun)
 path=ROOT/'reports/CLASSIFICATION_NIL31_SHORT_CHAIN_FLAG_20261107.json';assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ['status','matrix_shape','generic_rank','generic_residual_numerator_gcd','positive_parameters','seconds']}))
if __name__=='__main__':main()
