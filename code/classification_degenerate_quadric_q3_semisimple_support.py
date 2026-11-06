#!/usr/bin/env python3
"""Weakest support certificates: all six nonregular q3 semisimple types are empty."""
import argparse,hashlib,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
PREFIX='reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q3_'
CASES=['ASEMISIMPLE_CNIL2_BZERO','ASEMISIMPLE_CZERO_BZERO','CSEMISIMPLE_ANIL3',
       'ASEMISIMPLE_CNIL2_BNONZERO','CSEMISIMPLE_AZERO','ASEMISIMPLE_CZERO_BNONZERO']
def rank(M):return len(DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref()[1])
def span_contains(original,extra,n):
    polys=original+extra;mons=sorted(set().union(*(set(s.Poly(p,*n).monoms())for p in polys)))
    matrix=s.Matrix([[s.Poly(p,*n).coeff_monomial(m)for m in mons]for p in polys])
    assert rank(matrix)==rank(matrix[:len(original),:])
def main(output):
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    begun=time.monotonic();x=s.symbols('x0:5');n=s.symbols('n0:12');records=[];paths=[]
    for case in CASES:
        path=PREFIX+case+'_20261107.json';paths.append(path);source=json.loads((ROOT/path).read_text())
        J=[s.sympify(v)for v in source['residual_jacobi_coefficients']]
        Q=s.sympify(source['q']);F=s.cancel(s.sympify(source['euler_contracted_pfaffian_not_yet_actual_divisor'])/Q)
        record=dict(case=case,complete_source_dimension=source['source_dimension'],full_residual_Jacobi_retained=True)
        if case=='CSEMISIMPLE_ANIL3':
            quotient=s.factor(F/x[0]);assert quotient.is_polynomial(*x,*n)
            record.update(argument='The cubic has an x0 factor before any further Jacobi restriction.',cubic_factorization=str(s.factor(F)),linear_factor='x0',quotient=str(quotient))
        elif case=='ASEMISIMPLE_CNIL2_BNONZERO':
            quotient=s.factor(F/x[2]);assert quotient.is_polynomial(*x,*n)
            record.update(argument='The cubic has an x2 factor before any further Jacobi restriction.',cubic_factorization=str(s.factor(F)),linear_factor='x2',quotient=str(quotient))
        elif case=='ASEMISIMPLE_CNIL2_BZERO':
            span_contains(J,[n[3]],n)
            reduced=[s.expand(p.subs(n[3],0))for p in J];span_contains(reduced,[n[2]],n)
            reducedF=s.factor(F.subs({n[2]:0,n[3]:0}));assert s.cancel(reducedF/x[2]).is_polynomial(*x,*n)
            record.update(argument='Full Jacobi forces n3=0, then n2=0, after which the cubic has an x2 factor.',forced_zero_parameters=['n3','n2'],sequential_constant_span_certified=True,cubic_after_forced_conditions=str(reducedF),linear_factor='x2')
        elif case=='ASEMISIMPLE_CZERO_BNONZERO':
            span_contains(J,[(n[2]-2)*(n[6]+2)],n)
            first=s.factor(F.subs(n[2],2));second=s.factor(F.subs(n[6],-2))
            expected1=-(n[6]+2)*x[2]*(2*n[0]*x[0]*x[1]+2*n[1]*x[2]**2+5*x[2]*x[4])/10
            expected2=-(n[2]-2)*x[0]*x[1]*(2*n[5]*x[2]-5*x[4])/10
            assert s.expand(first-expected1)==s.expand(second-expected2)==0
            record.update(argument='The forced product (n2-2)(n6+2)=0 gives two exhaustive branches, both with a proper cubic factor.',forced_product='(n2-2)*(n6+2)',constant_span_certified=True,branches=[{'substitution':{'n2':'2'},'cubic':str(first)},{'substitution':{'n6':'-2'},'cubic':str(second)}])
        elif case=='ASEMISIMPLE_CZERO_BZERO':
            p,q,r,t=n[2],n[3],n[4],n[5];trace=n[7]+n[11]
            span_contains(J,[p*r,q*t,p*t+q*r],n)
            without_rs=[s.expand(v.subs({r:0,t:0}))for v in J]
            span_contains(without_rs,[(2*trace+3)*p,(2*trace+3)*q],n)
            first=s.factor(F.subs({p:0,q:0}))
            expected=-(n[0]*x[0]*x[1]+n[1]*x[2]**2)*((4*trace+6)*x[2]-4*r*x[3]-4*t*x[4])/10
            assert s.expand(first-expected)==0
            second=s.factor(F.subs({r:0,t:0,n[11]:-s.Rational(3,2)-n[7]},simultaneous=True))
            assert s.cancel(second/(x[0]*x[1])).is_polynomial(*x,*n)
            record.update(argument='The symmetric outer-product equations imply (p,q)=0 or (r,s)=0. In the latter branch, (p,q) nonzero forces trace=-3/2. These two cases cover every point and force a reducible or zero cubic.',
                first_conditions=list(map(str,[p*r,q*t,p*t+q*r])),second_conditions=list(map(str,[(2*trace+3)*p,(2*trace+3)*q])),sequential_constant_spans_certified=True,
                branches=[{'substitution':{'n2':'0','n3':'0'},'cubic':str(first)},{'substitution':{'n4':'0','n5':'0','n11':'-3/2-n7'},'open_condition':'(n2,n3) != (0,0)','cubic':str(second),'proper_factor':'x0*x1'}])
        else:
            K=(n[1]+n[3])*x[0]-(n[2]+n[4])*x[1]+(n[0]+n[11])*x[2]
            R=n[5]*x[0]**2+n[6]*x[0]*x[1]+n[7]*x[0]*x[2]+n[8]*x[1]**2+n[9]*x[1]*x[2]+n[10]*x[2]**2
            j1=(2*n[0]*n[2]-n[0]-2*n[11]*n[4]-n[11])/2
            j2=(2*n[1]*n[4]+n[1]-2*n[2]*n[3]+n[3])/2
            j3=-n[0]*n[1]+n[11]*n[3]
            span_contains(J,[j1,j2,j3],n)
            target=-K*(R/5+x[3]*x[4]/2)
            difference=x[3]*x[4]*(j2*x[0]-j3*x[1]+j1*x[2])/5
            assert s.expand(F-target-difference)==0
            record.update(argument='Three exact Jacobi coefficients make the entire cubic equal a linear form times a quadratic. If the linear form vanishes identically, the cubic is zero.',required_Jacobi_conditions=list(map(str,[j1,j2,j3])),constant_span_certified=True,
                linear_factor=str(K),quadratic_factor=str(-R/5-x[3]*x[4]/2),exact_remainder_identity=str(difference),cubic_mod_Jacobi=str(target))
        record['irreducible_cubic_locus_empty']=True;records.append(record)
        print(json.dumps({'case':case,'irreducible_cubic_locus_empty':True}),flush=True)
    inputs=paths+['reports/CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json','code/classification_degenerate_quadric_q3_semisimple_support.py']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,cases=records,
        all_six_nonregular_q3_semisimple_irreducible_cubic_loci_empty=True,
        source_completeness_independent_review='reports/CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json',
        method='Only necessary consequences of every original Jacobi source and exact cubic factorizations; no coefficient discarded and no genericity or finite-field inference.',
        global_exhaustion_proved=False,input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in inputs},producer_source=Path(__file__).read_text(),seconds=time.monotonic()-begun)
    with output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
