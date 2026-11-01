#!/usr/bin/env python3
"""Exact b8=0 fourth initial space on the T7 flag core; supplied directory only."""
import json,pickle,time,sys
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
ROOT=Path(__file__).resolve().parents[1];WD=Path(sys.argv[1]);sys.path.insert(0,str(ROOT/'code'))
from s61_initial_quartic_from_jets_20261101 import initial_quartic
start=time.time();cache=pickle.load(open(WD/'flag-rankone-socle-flint-raw4-cache.pkl','rb'));old=list(map(s.Symbol,cache['variables']));removed=old.index(s.Symbol('b8'));variables=old[:removed]+old[removed+1:];K=s.QQ.frac_field(*variables);ud=s.Rational(cache['ud']);td=s.Rational(cache['td'])
def matrix(d,key='A'):
 f=1/(ud**(d-1)*td**d)
 return s.Matrix([[sum(s.Rational(c)*f*s.prod(x**k for i,(x,k)in enumerate(zip(old,ex))if i!=removed)for ex,c in p if ex[removed]==0)for p in row]for row in cache[key+str(d)]])
r=initial_quartic(*[matrix(d)for d in [2,3,4]],9,K,progress=lambda x:print(x,time.time()-start,flush=True));L,E=r['quartic_lower'],r['quartic_extra'];rr,piv=L.rref();rank=len(piv);np=[i for i in range(495)if i not in piv]
def project(A):return A.extract(range(A.shape[0]),np)-A.extract(range(A.shape[0]),piv)*rr.extract(range(rank),np)
N=project(E);nr,pc=N.rref();rows=list(N.transpose().rref()[1]);print('quartic ranks',rank,len(pc),time.time()-start,flush=True)
mi=s.factor(K.to_sympy(N.extract(rows,pc).det()));Lform=-200*s.Symbol('b10')-91*s.Symbol('b11')+200*s.Symbol('b9');assert s.cancel(mi/Lform**2).is_Rational
print('extra minor constant times L squared',flush=True)
B3=r['degree3_leading'].extract(r['cubic_basis_source_rows'],range(165));tail=r['degree3_tail'].extract(r['cubic_basis_source_rows'],range(495));p3=r['cubic_pivots'];t3=matrix(3,'W')[-1,:].to_DM().convert_to(K);t4=matrix(4,'W')[-1,:].to_DM().convert_to(K);ct=t3.extract([0],p3)*B3.extract(range(B3.shape[0]),p3).inv();assert ct*B3==t3;nt=project(t4-ct*tail);nc=[i for i in range(nt.shape[1])if i not in pc];tc=nt.extract([0],nc)-nt.extract([0],pc)*nr.extract(range(len(pc)),nc)
witness=next(((j,str(K.to_sympy(c)))for(i,j),c in tc.to_dok().items()if c and K.to_sympy(c).is_Rational),None)
den=set()
for D in [rr,nr,ct,tc,r['quadratic_row_reduction'],r['cubic_row_reduction']]:
 for v in D.to_dok().values():
  denominator=K.to_sympy(v).as_numer_denom()[1];degree=s.degree(denominator,s.Symbol('b9'));constant=s.cancel(denominator/Lform**degree);assert constant.is_Rational and constant;den.add(int(degree))
lower_rows=list(L.transpose().rref()[1]);lower_minor=K.to_sympy(L.extract(lower_rows,piv).det());cubic_minor=K.to_sympy(B3.extract(range(B3.shape[0]),p3).det());assert lower_minor.is_Rational and lower_minor and cubic_minor.is_Rational and cubic_minor
assert witness is not None and r['quadratic_rank']==9 and r['cubic_initial_rank']==66 and rank==275 and len(pc)==14

out=dict(variables=list(map(str,variables)),I2=r['quadratic_rank'],I3=r['cubic_initial_rank'],I4=rank+len(pc),lower4=rank,extra4=len(pc),HF4=495-rank-len(pc),extra_minor=str(mi),denominator_L_powers=sorted(den),lower_minor=str(lower_minor),cubic_minor=str(cubic_minor),all_L_nonzero_points=True,t4_nonzero=not tc.is_zero_matrix,t4_constant_witness=witness,seconds=time.time()-start)
(WD/'flag-rankone-socle-zero8-four.json').write_text(json.dumps(out,indent=2)+'\n');print('ZERO8 FOURTH PASS',out['seconds'],flush=True)
