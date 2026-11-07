#!/usr/bin/env python3
"""Independent full-source nil31-short dimension bound; no producer/CAS import."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

from classification_quadratic_cofactor_review import Ring,verify_source
from classification_minimal_quartic_review import LaurentRing
from classification_unimodular_logarithmic_review import RF,U,ug,ud,rank,det

ROOT=Path(__file__).resolve().parents[1]
LOWER=[3,4,8,12,13,14]
MAG=[0,1,2,5,6,7]
CORE=[10,11,9,15,16,17]

def rfparse(text,parameter):
    def read(n):
        if isinstance(n,ast.Constant):return RF(F(n.value))
        if isinstance(n,ast.Name):
            assert n.id==parameter,n.id
            return U
        if isinstance(n,ast.UnaryOp):return -read(n.operand) if isinstance(n.op,ast.USub) else read(n.operand)
        a,b=read(n.left),read(n.right)
        if isinstance(n.op,ast.Add):return a+b
        if isinstance(n.op,ast.Sub):return a-b
        if isinstance(n.op,ast.Mult):return a*b
        if isinstance(n.op,ast.Div):return a/b
        if isinstance(n.op,ast.Pow):
            assert len(b.p)<=1 and b.q==(F(1),)
            k=int(b.p[0]) if b.p else 0
            out=RF(1)
            for _ in range(abs(k)):out=out*a
            return out if k>=0 else 1/out
        raise ValueError(ast.dump(n))
    return read(ast.parse(text,mode='eval').body)

def core_completeness():
    R=LaurentRing(['a','b','t','j','k','c'],[],invert=['c'])
    raw=['(a+b+1/6)*t+c*k','-j*t-t*t+(b+1/3)*k','c*(a-b+1/3)',
         'c*j-a*b/2+a/6-c*t/2-b*b/2+b/6',
         'j/9-2*(a+b)*t/3+t/9+c*k','j/12-(a+b)*t/4+t/8+c*k']
    def subst(p,values):
        result={}
        for ex,co in p.items():
            term=R.constant(co)
            for name,power in zip(R.names,ex):
                for _ in range(power):term=R.mul(term,values.get(name,R.variable(name)))
            result=R.add(result,term)
        return result
    # Scaled c!=0 chart: original t=tau/c,j=J/c,k=K/c^2.
    values={'a':R.parse('b-1/3'),'j':R.parse('(t/2+b*b-b/2+1/18)/c'),
            'k':R.parse('-(2*b-1/6)*t/c**2'),'t':R.parse('t/c')}
    actual=[subst(R.parse(q),values) for q in raw]
    L='(b-1/3)*(b-1/6)-30*(b-1/6)*t'
    expected=['0','-3*t*(t+2*b*b)/(2*c**2)','0','0','('+L+')/(9*c)','('+L+')/(12*c)']
    assert actual==list(map(R.parse,expected))
    assert R.parse('(b-1/3)*(b-1/6)+60*(b-1/6)*b*b')==R.parse('60*(b-1/6)*(b+1/12)*(b-1/15)')
    # The c=0 chart has these four exact generators, without a localization.
    zero=[R.parse(q) for q in ['(a+b+1/6)*t','-j*t-t*t+(b+1/3)*k',
            '-(a+b)*(b-1/3)/2','j/9-2*(a+b)*t/3+t/9','j/12-(a+b)*t/4+t/8']]
    assert R.add(R.scale(zero[3],9),R.scale(zero[0],6))==R.parse('j+2*t')
    return dict(nonzero_chart_two_remaining_equations=[L,'t*(t+2*b*b)'],
                cubic_factorization='60*(b-1/6)*(b+1/12)*(b-1/15)',
                no_division_by_b_minus_one_sixth=True,
                zero_chart_reduced_by_j_equals_minus_two_t=True)

def branch_data():
    out=[]
    for b,tau in [(F(1,3),0),(F(1,6),0),(F(1,6),F(-1,18)),(F(-1,12),F(-1,72)),(F(1,15),F(-2,225))]:
        a=b-F(1,3);j=(RF(tau/2)+b*b-b/2+F(1,18))/U;k=-(2*b-F(1,6))*tau/(U*U)
        out.append(('n17',True,{17:U,16:RF(b),15:tau/U,11:RF(a),10:j,9:k}))
    out += [('n15',True,{17:RF(),15:U,10:-2*U,11:RF(F(-1,2)),16:RF(F(1,3)),9:RF(F(-3,2))*U*U}),
            ('n9',True,{17:RF(),15:RF(),10:RF(),11:RF(F(1,3)),16:RF(F(-1,3)),9:U}),
            ('n11',False,{17:RF(),15:RF(),10:RF(),11:U,16:RF(F(1,3)),9:RF()}),
            ('n16',False,{17:RF(),15:RF(),10:RF(),11:-U,16:U,9:RF()})]
    return out

def substitute_core(poly,values):
    out={}
    for ex,co in poly.items():
        v=RF(co);rest=[]
        for i,power in enumerate(ex[:18]):
            if i in values:
                for _ in range(power):v=v*values[i]
            else:rest.append(power)
        assert not any(ex[18:])
        key=tuple(rest);out[key]=out.get(key,RF())+v
    return {k:v for k,v in out.items() if v}

def extract_matrices(equations,values):
    remaining=[i for i in range(18) if i not in values]
    assert set(remaining)==set(LOWER+MAG)
    lower=[];magnetic=[]
    for p in equations:
        q=substitute_core(p,values)
        if not q:continue
        ismag=any(any(ex[remaining.index(i)] for i in MAG) for ex in q)
        variables=MAG if ismag else LOWER
        row=[RF() for _ in variables]
        for ex,co in q.items():
            if ismag:
                degree=sum(ex[remaining.index(i)] for i in MAG)
                assert degree<=1
                if not degree:continue
                assert sum(ex)==1
            else:assert sum(ex)<=1
            if sum(ex):
                index=remaining[ex.index(1)]
                assert index in variables
                row[variables.index(index)]+=co
        (magnetic if ismag else lower).append(row)
    return lower,magnetic

def encode_matrix(M):
    return [[dict(numerator=list(map(str,v.p)),denominator=list(map(str,v.q))) for v in row]for row in M]

def rank_certificate(M,reported,parameter,excluded_zero):
    r=rank(M);assert r==reported['generic_rank']
    G=()
    for rr in it.combinations(range(len(M)),r):
        for cc in it.combinations(range(len(M[0])),r):
            value=det([[M[i][j]for j in cc]for i in rr])
            if value:G=ug(G,value.p)
    target=rfparse(reported['maximum_minor_gcd'],parameter)
    assert target.q==(F(1),)
    assert G==tuple(v/target.p[-1] for v in target.p)
    roots=[F(row['value'])for row in reported['exceptional']]
    rem=G
    for value in roots:
        while len(rem)>1:
            quotient,remainder=ud(rem,(-value,F(1)))
            if remainder:break
            rem=quotient
    assert len(rem)<=1,'unlisted algebraic rank-drop roots'
    def evaluate(v,x):
        a=sum((co*x**i for i,co in enumerate(v.p)),F())
        b=sum((co*x**i for i,co in enumerate(v.q)),F())
        assert b
        return RF(a/b)
    specials={}
    for row in reported['exceptional']:
        value=F(row['value'])
        if excluded_zero and not value:continue
        actual=rank([[evaluate(v,value)for v in rr]for rr in M])
        assert actual==row['rank']
        specials[value]=actual
    return r,G,specials

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    source_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_SUPPORT_20261107.json'
    data=json.loads(path.read_text());source=json.loads(source_path.read_text())['full_nilpotent_sources'][3]
    assert data['field']==source['field_linear_matrix'] and data['translation']==source['field_constant']
    full_source=verify_source(source,4);assert len(full_source['free_indices'])==18
    R=Ring(18)
    monic=lambda p:R.scale(p,1/max(p.items())[1])
    equations=list(map(R.parse,data['all_reduced_jacobi_equations']))
    assert {tuple(sorted(monic(p).items()))for p in equations}=={tuple(sorted(monic(R.parse(p)).items()))for p in source['distinct_nonzero_jacobi_coefficients']}
    core=[p for p in equations if all(not any(ex[i]for i in range(18)if i not in CORE)for ex in p)]
    assert len(core)==6
    rawcore=['(n11+n16+1/6)*n15+n17*n9','-n10*n15-n15*n15+(n16+1/3)*n9',
             'n17*(n11-n16+1/3)',
             'n17*n10-n11*n16/2+n11/6-n17*n15/2-n16*n16/2+n16/6',
             'n10/9-2*(n11+n16)*n15/3+n15/9+n17*n9',
             'n10/12-(n11+n16)*n15/4+n15/8+n17*n9']
    assert {tuple(sorted(monic(p).items()))for p in core}=={tuple(sorted(monic(R.parse(p)).items()))for p in rawcore}
    complete=core_completeness();records=[]
    assert len(data['whole_source_dimension_bounds'])==9
    for (parameter,nonzero,values),claim in zip(branch_data(),data['whole_source_dimension_bounds']):
        assert parameter==claim['parameter'] and nonzero==claim['parameter_nonzero']
        assert all(not substitute_core(p,values)for p in core)
        for name,value in claim['substitution'].items():assert values[int(name[1:])]==rfparse(value,parameter)
        lower,magnetic=extract_matrices(equations,values)
        for M,reported in [(lower,claim['lower_rank']),(magnetic,claim['magnetic_rank'])]:
            assert M==[[rfparse(v,parameter)for v in row]for row in reported['matrix']]
        lr,lg,ls=rank_certificate(lower,claim['lower_rank'],parameter,nonzero)
        mr,mg,ms=rank_certificate(magnetic,claim['magnetic_rank'],parameter,nonzero)
        generic=13-lr-mr
        specials=[]
        def at_rank(M,x):
            evaluate=lambda v:RF(sum((co*x**i for i,co in enumerate(v.p)),F())/sum((co*x**i for i,co in enumerate(v.q)),F()))
            return rank([[evaluate(v)for v in row]for row in M])
        for value in sorted(set(ls)|set(ms)):
            a,b=at_rank(lower,value),at_rank(magnetic,value)
            specials.append(dict(parameter=str(value),lower_rank=a,magnetic_rank=b,bound=12-a-b))
        upper=max([generic]+[row['bound']for row in specials])
        assert upper==claim['fixed_source_dimension_upper_bound']
        records.append(dict(parameter=parameter,nonzero=nonzero,lower_matrix=encode_matrix(lower),
            magnetic_matrix=encode_matrix(magnetic),generic_ranks=[lr,mr],gcds=[list(map(str,lg)),list(map(str,mg))],
            all_special_ranks=specials,source_dimension_upper_bound=upper))
    # All20 infinitesimal affine transformations, in all20 affine-field coords.
    A=[[F(v)for v in row]for row in data['field']];c=list(map(F,data['translation']))
    columns=[]
    for i,j in it.product(range(4),repeat=2):
        B=[[F(a==i and b==j)for b in range(4)]for a in range(4)]
        columns.append([sum(B[a][k]*A[k][b]-A[a][k]*B[k][b]for k in range(4))for a,b in it.product(range(4),repeat=2)]+[sum(B[a][k]*c[k]for k in range(4))for a in range(4)])
    for j in range(4):columns.append([F()]*16+[-A[a][j]for a in range(4)])
    affine_rank=rank([[RF(v)for v in row]for row in zip(*columns)])
    assert affine_rank==13
    maximum=max(row['source_dimension_upper_bound']for row in records);assert maximum==7
    result=dict(task_id='E_reduced_nonnormal',no_producer_import=True,no_CAS_import=True,all_passed=True,
        source=full_source,complete_core_reduction=complete,all_nine_branch_bounds=records,
        affine_cofactor_orbit_rank=affine_rank,full20by20_affine_orbit_matrix=[[str(v)for v in row]for row in zip(*columns)],
        fixed_source_dimension_upper_bound=maximum,projective_saturated_dimension_upper_bound=maximum+affine_rank+4-1,
        E_membership_proved=False,G_generic_type_excluded=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [path,source_path]},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','fixed_source_dimension_upper_bound','projective_saturated_dimension_upper_bound','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
