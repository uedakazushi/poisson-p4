#!/usr/bin/env python3
"""Complete residual-Jacobi support for q4 semisimple2 cofactor sources."""
import argparse,hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
PREFIX='reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_SEMISIMPLE2_'
def decode(r):
    P=s.zeros(5)
    for key,value in r['brackets'].items():
        i,j=map(int,key);P[i,j]=s.sympify(value);P[j,i]=-P[i,j]
    return P
def jacs(P,x):
    return [s.expand(sum(P[i,l]*s.diff(P[j,k],x[l])+P[j,l]*s.diff(P[k,i],x[l])+P[k,l]*s.diff(P[i,j],x[l])for l in range(5)))for i,j,k in itertools.combinations(range(5),3)]
def divisor(P,x,check=True):
    pf=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5)if a!=missing]
        pf.append(s.expand((-1)**missing*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])))
    D=s.factor(sum(x[i]*pf[i]for i in range(5))/5)
    if check:assert all(s.expand(s.diff(D,x[i])-pf[i])==0 for i in range(5))
    return D
def rank(A):return len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1])
def coefficient_span(polys,variables):
    mons=sorted(set().union(*(set(s.Poly(p,*variables).monoms())for p in polys)))
    return s.Matrix([[s.Poly(p,*variables).coeff_monomial(m)for m in mons]for p in polys])
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(70,71));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();xx=s.symbols('x0:5');q=xx[0]*xx[1]+xx[2]*xx[3];n=s.symbols('n0:11');a,b,c,d,e,f,g,h,i,j,k=n
    paths=[PREFIX+kind+'_20261107.json'for kind in ['ZERO','ISOTROPIC','NONISOTROPIC']]
    zero,iso,noniso=[json.loads((ROOT/p).read_text())for p in paths]
    original=[s.sympify(v)for v in zero['residual_jacobi_coefficients']]
    conditions=[a*e+b*h,a*h+c*f,b*e+c*d,c*e,c*h,a*f,b*d,b*f-a*(k+2),a*d-b*(k+2),c*(k+2),c*i-a*(g-3),c*j+b*(k+g+4)]
    together=coefficient_span(original+conditions,n)
    assert rank(together)==rank(together[:len(original),:])==rank(together[len(original):,:])
    Pzero=decode(zero);zero_records=[]
    branches=[dict(name='c_nonzero',sub={d:0,e:0,f:0,h:0,k:-2,i:a*(g-3)/c,j:-b*(g+2)/c},expected=xx[0]*xx[1]*(a*xx[2]+b*xx[3]+c*xx[4])),
              dict(name='a_nonzero_c_b_zero',sub={b:0,c:0,d:0,e:0,f:0,h:0,k:-2,g:3},expected=a*xx[0]*xx[1]*xx[2]),
              dict(name='b_nonzero_c_a_zero',sub={a:0,c:0,d:0,e:0,f:0,h:0,k:-2,g:-2},expected=b*xx[0]*xx[1]*xx[3]),
              dict(name='abc_zero_rank_at_most2',sub={a:0,b:0,c:0},expected=0)]
    for branch in branches:
        P=Pzero.subs(branch['sub'],simultaneous=True).applyfunc(s.cancel)
        assert all(s.cancel(v)==0 for v in jacs(P,xx))
        D=divisor(P,xx);assert s.cancel(D-q*branch['expected'])==0
        zero_records.append(dict(case=branch['name'],substitution={str(k):str(v)for k,v in branch['sub'].items()},actual_divisor=str(D),cubic=str(branch['expected'])))
    # c=0,a*b!=0 forces g=3 and g=-2, so this missing stratum is empty.
    iso_j=[s.sympify(v)for v in iso['residual_jacobi_coefficients']]
    span=coefficient_span(iso_j+[b],n[:6]);assert rank(span)==rank(span[:-1,:])
    P_iso=decode(iso).subs(b,0);potential=divisor(P_iso,xx,check=False)
    assert s.cancel(potential/(q*xx[2])).is_polynomial(*xx,*n[:6])
    noniso_records=[];Pgen=decode(noniso)
    noniso_j=[s.sympify(v)for v in noniso['residual_jacobi_coefficients']]
    U=c+2*f+1;V=c-f
    factors=[(3*a-4)*U,(3*b+1)*V,(a-2*b-2)*U,(2*a-b-3)*V]
    assert all(any(s.expand(v-p)==0 or s.expand(v+p)==0 for v in noniso_j)for p in factors)
    for substitution,required in [({a:s.Rational(4,3),b:s.Rational(-1,3)},[f-2*c-2,2*c+3*f]),
                                 ({c:s.Rational(-1,3),f:s.Rational(-1,3)},[a-3,b+2])]:
        values=[s.expand(v.subs(substitution))for v in noniso_j]
        matrix=coefficient_span(values+required,n[:6])
        assert rank(matrix)==rank(matrix[:len(values),:])==rank(matrix[len(values):,:])==2
    for branch_id,values in enumerate([(s.Rational(4,3),s.Rational(-1,3),s.Rational(-3,4),s.Rational(1,2)),(3,-2,s.Rational(-1,3),s.Rational(-1,3))],1):
        sub=dict(zip([a,b,c,f],values));P=Pgen.subs(sub)
        assert all(v==0 for v in jacs(P,xx))
        D=divisor(P,xx);cubic=s.factor(s.cancel(D/q))
        quotient=s.cancel(cubic/(xx[2]-xx[3]))
        assert quotient.is_polynomial(*xx)and s.Poly(quotient,*xx).total_degree()==2
        noniso_records.append(dict(branch=branch_id,parameter_substitution={str(k):str(v)for k,v in sub.items()},
            remaining_all_parameters=['n3','n4'],full_brackets={''.join(map(str,ij)):str(P[ij])for ij in itertools.combinations(range(5),2)},
            all_ten_universal_Jacobi_zero=True,all_signed_Pfaffians_equal_divisor_gradient=True,
            actual_divisor=str(D),cubic_factorization=str(cubic),common_linear_factor=str(xx[2]-xx[3]),
            cubic_quotient=str(quotient),irreducible_cubic_locus_empty=True))
        print(json.dumps(dict(branch=branch_id,cubic=str(cubic),irreducible_cubic_locus_empty=True)),flush=True)
    inputs=paths+['code/classification_degenerate_quadric_semisimple2_support.py']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,
        zero_complete_reduced_conditions=list(map(str,conditions)),zero_conditions_equal_full_Jacobi_constant_span=True,
        zero_all_support_branches=zero_records,zero_c0_ab_nonzero_empty_reason='The equations imply g=3 and g=-2.',zero_irreducible_cubic_locus_empty=True,
        isotropic_n1_in_exact_Jacobi_constant_span=True,isotropic_potential_after_n1_zero=str(potential),isotropic_irreducible_cubic_locus_empty=True,
        nonisotropic_complete_branches=noniso_records,
        nonisotropic_completeness_reason='Set U=n2+2*n5+1,V=n2-n5. If either nonzero, factored Jacobi forces n0=4/3,n1=-1/3 and the remaining linear equations force n2=-3/4,n5=1/2. If U=V=0, n2=n5=-1/3 and the remaining equations force n0=3,n1=-2.',
        nonisotropic_two_case_linear_reductions_exactly_verified=True,
        all_three_fixed_cofactor_types_irreducible_cubic_locus_empty=True,
        ambient_rank_tests_unnecessary_and_not_run=True,global_exhaustion_proved=False,
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
