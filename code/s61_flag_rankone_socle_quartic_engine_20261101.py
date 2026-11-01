#!/usr/bin/env python3
"""True fourth initial space on the seven-parameter rank-one socle core."""
import sys,pickle,time,json
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_initial_quartic_from_jets_20261101 import initial_quartic
WD=Path(sys.argv[1]);start=time.time();cache=pickle.load(open(WD/'flag-rankone-socle-flint-raw4-cache.pkl','rb'))
variables=list(map(s.Symbol,cache['variables']));ud=s.Rational(cache['ud']);td=s.Rational(cache['td']);K=s.QQ.frac_field(*variables)
def matrix(d,key='A'):
    factor=1/(ud**(d-1)*td**d)
    return s.Matrix([[sum(s.Rational(c)*factor*s.prod(x**k for x,k in zip(variables,ex)) for ex,c in p) for p in row] for row in cache[key+str(d)]])
A2,A3,A4=map(matrix,[2,3,4]);result=initial_quartic(A2,A3,A4,9,K,progress=lambda x:print(x,time.time()-start,flush=True))
assert result['quadratic_rank']==9 and result['cubic_initial_rank']==66
L,E=result['quartic_lower'],result['quartic_extra'];rr,piv=L.rref();rank=len(piv);assert rank==275
nonp=[i for i in range(495) if i not in piv]
def project(X):return X.extract(range(X.shape[0]),nonp)-X.extract(range(X.shape[0]),piv)*rr.extract(range(rank),nonp)
N=project(E);nr,pc=N.rref();assert len(pc)==14;extra_rows=list(N.transpose().rref()[1]);assert len(extra_rows)==14
print('quartic lower275 plus14 PASS',time.time()-start,flush=True)
minor=s.factor(K.to_sympy(N.extract(extra_rows,pc).det()))
assert s.cancel(minor/s.Symbol('b8')**25).is_Rational
lower_rows=list(L.transpose().rref()[1]);lower_minor=s.factor(K.to_sympy(L.extract(lower_rows,piv).det()))
B3=result['degree3_leading'].extract(result['cubic_basis_source_rows'],range(165));tail=result['degree3_tail'].extract(result['cubic_basis_source_rows'],range(495));p3=result['cubic_pivots']
t3=matrix(3,'W')[-1,:].to_DM().convert_to(K);t4=matrix(4,'W')[-1,:].to_DM().convert_to(K)
ct=t3.extract([0],p3)*B3.extract(range(B3.shape[0]),p3).inv();assert ct*B3==t3
nt=project(t4-ct*tail);nc=[i for i in range(nt.shape[1]) if i not in pc];tc=nt.extract([0],nc)-nt.extract([0],pc)*nr.extract(range(14),nc)
assert not tc.is_zero_matrix
witness=next(((j,str(K.to_sympy(c)))for (i,j),c in tc.to_dok().items()if K.to_sympy(c).is_Rational and c),None)
enc=lambda M:[[str(v)for v in row]for row in M.to_Matrix().tolist()]
denominators=set()
for D in [rr,nr,ct,tc,result['quadratic_row_reduction'],result['cubic_row_reduction']]:
    for c in D.to_dok().values():
        terms=list(c.denom.terms());assert len(terms)==1
        ex,_=terms[0];assert all(k==0 for j,k in enumerate(ex)if variables[j]!=s.Symbol('b8'))
        denominators.add(ex[variables.index(s.Symbol('b8'))])
for name,value in [('lower',lower_minor),('cubic',K.to_sympy(B3.extract(range(66),p3).det()))]:
    poly=s.Poly(value,*variables);assert len(poly.terms())==1
    ex,co=poly.terms()[0];assert co and all(k==0 for j,k in enumerate(ex)if variables[j]!=s.Symbol('b8'))
    print(name,'minor b8 exponent',ex[variables.index(s.Symbol('b8'))],flush=True)
    if name=='cubic':cubic_minor=str(value)
assert witness is not None
report=dict(id='S61_FLAG_RANKONE_SOCLE_QUARTIC_20261101',all_passed=True,field=str(K),parameters=list(map(str,variables)),core='r=-91/400',quadratic_rank=9,cubic_initial_rank=66,quartic_lower_rank=275,quartic_additional_rank=14,quartic_initial_rank=289,hilbert_function=[1,9,36,99,206],quadratic_basis_rows=result['quadratic_basis_rows'],cubic_basis_source_rows=result['cubic_basis_source_rows'],quartic_extra_source_rows=extra_rows,quartic_extra_pivot_columns=list(pc),quartic_extra_minor=str(minor),quartic_extra_rref=enc(nr.extract(range(14),range(nr.shape[1]))),quartic_quotient_columns=nonp,t3_in_true_cubic_ideal=True,corrected_t4_nonzero=True,constant_t4_witness=witness,quartic_t_class=[str(v) for v in tc.to_Matrix()],row_reduction_denominator_b8_exponents=sorted(denominators),lower_minor=str(lower_minor),cubic_minor=cubic_minor,all_b8_nonzero_points=True,seconds=time.time()-start)
(WD/'flag-rankone-socle-quartic.json').write_text(json.dumps(report,indent=2)+'\n');print('TRUE QUARTIC AND ORD T4 PASS',report['seconds'],flush=True)
