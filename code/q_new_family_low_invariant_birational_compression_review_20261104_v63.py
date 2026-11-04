#!/usr/bin/env python3
"""Independent reverse sparse substitution and Laurent coordinate-isomorphism audit."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
import math
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW'
def read(p):
    b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
    d=read(p);assert d['all_passed']
    for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in d.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return d
def mul(A,B):
    C={}
    for x,c in A.items():
        for y,d in B.items():
            z=tuple(i+j for i,j in zip(x,y));C[z]=C.get(z,Q(0))+c*d
    return {x:c for x,c in C.items() if c}
def scale(A,c):return {x:v*c for x,v in A.items() if v*c}
def shift(A,e):return {tuple(x+y for x,y in zip(a,e)):v for a,v in A.items()}
def add(*As):
    d={}
    for A in As:
        for x,c in A.items():d[x]=d.get(x,Q(0))+c
    return {x:c for x,c in d.items() if c}
def enc(A):return [[str(c.numerator),str(c.denominator),*e] for e,c in sorted(A.items(),reverse=True)]

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63'
    old=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63'
    ap=base/'report.json';author=accepted(ap)
    rev=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63_REVIEW/report.json';accepted(rev)
    fieldrev=ROOT/'reports/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW/report.json';fields=accepted(fieldrev)
    assert fields['known_full_rawjet_order2_equivalence_and_rational_section_fixed_field']
    recipepath=base/'coordinate-recipe.json';recipe=read(recipepath)
    assert not recipe['actual_NEW_target_field_equality_claimed']
    # Monomial maps between old(a,Y,k,B) and new(xi,r,s,B).
    forward=[(Q(40),(-1,0,1,-1)),(Q(1,1600),(2,1,0,0)),(Q(48),(1,1,0,0)),(Q(1),(0,0,0,1))]
    inverse=[(Q(76800),(0,1,-1,0)),(Q(1,3686400),(0,-1,2,0)),(Q(1920),(1,1,-1,1)),(Q(1),(0,0,0,1))]
    def compose(m,images):
        c,ex=m;v=c;power=[0]*4
        for n,(k,e) in zip(ex,images):
            v*=k**n
            for i in range(4):power[i]+=n*e[i]
        return v,tuple(power)
    identity=[(Q(1),tuple(int(i==j) for i in range(4))) for j in range(4)]
    assert [compose(m,forward) for m in inverse]==identity
    assert [compose(m,inverse) for m in forward]==identity
    assert Q(3072)*Q(1,3072)==1 and Q(3072)*Q(1,3072)-1==0
    H={(0,1,0,0):Q(147456),(0,0,2,0):Q(-9437184),(0,0,1,0):Q(6144),(0,0,0,0):Q(-1)}
    h={(0,2,0,0):Q(1),(1,0,2,0):Q(-25)}
    def reverse(A):
        out={}
        for (i,j,k,l),c in A.items():
            c*=Q(1,1600)**i*Q(48)**j
            for w in range(k+1):
                ex=(2*i+j,i+j,w,l);v=c*math.comb(k,w)*Q(3072)**w*(-1)**(k-w)
                out[ex]=out.get(ex,Q(0))+v
        return {e:c for e,c in out.items() if c}
    assert reverse(h)==scale(shift(H,(2,1,0,0)),Q(1,64))
    paths=[ap,rev,fieldrev,recipepath];maps={};records=[]
    names=['F4','p0','p1','p2','p3','p4','N12'];expected=[5,6,29,33,31,30,43]
    for i,n in enumerate(names):
        np=base/(n+'.json.gz');op=old/(n+'.json.gz');paths.extend([np,op]);nd=read(np);od=read(op)
        M={tuple(map(int,t[1:])):Q(t[0])/Q(nd['denominator_integer']) for t in nd['terms']}
        assert len(M)==expected[i] and len(M)==len(nd['terms'])
        cp,hp,bp=od['source_c_power'],od['denominator_H_power'],od['denominator_B_power']
        assert nd['xi_power']==cp and nd['denominator_h_power']==hp
        oldq={tuple(map(int,t[2:])):Q(int(t[0]),int(t[1])) for t in od['terms']}
        lhs=scale(shift(reverse(M),(0,0,0,bp)),Q(40)**cp*Q(64)**hp*Q(od['denominator_integer']))
        rhs=shift(oldq,(cp+2*hp,cp//2+hp,0,cp));assert lhs==rhs
        maps[n]=M
        deg=[max(e[j] for e in M) for j in range(4)]
        assert deg==author['results'][i]['multidegree']
        assert author['results'][i]['terms']==len(M)
        records.append(dict(name=n,terms=len(M),h_denominator=hp,xi_power=cp,multidegree=deg,
                            inverse_substitution_all_coefficients_equal=True))
    F={(2,0,2,0):Q(100),(1,0,2,0):Q(-100),(1,1,0,0):Q(-4),(0,2,0,0):Q(1),(1,0,0,0):Q(4)}
    assert maps['F4']==scale(F,Q(2,405))
    minus={(1,0,0,0):Q(1),(0,0,0,0):Q(-1)}
    linear={(0,1,0,0):Q(-1),(0,0,1,0):Q(5),(0,0,0,0):Q(1)}
    D=scale(shift(mul(minus,linear),(4,0,2,0)),Q(2359296000))
    assert maps['p0']==scale(shift(D,(0,0,0,1)),Q(4,405**2))
    assert add(maps['p0'],scale(maps['p4'],Q(4)))==scale(mul(h,mul(maps['F4'],maps['F4'])),Q(1,400))
    assert D and all(e[3]==0 for e in D)
    cert=dict(old_monomial_coordinates=['a','Y','k=c/X','B'],new_monomial_coordinates=['xi','r','s','B'],
              forward=[dict(constant=str(c),exponents=e) for c,e in forward],inverse=[dict(constant=str(c),exponents=e) for c,e in inverse],
              both_compositions_exact_identity=True,all7_reverse_substitution_checks=records,F=enc(F),h=enc(h),D=enc(D),
              B_recovery_identity='Q0*h*F^2=D*B',normalized_p0_identity='Q0=1/400-4*p4/F4^2',
              domain='original source open intersect r*s*B*h*F*D!=0',target_field_index_not_determined=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(cert,indent=2)+'\n')
    paths += [Path(__file__),ROOT/f'research/notes/{STEM}.md',ROOT/'code/q_new_family_low_invariant_birational_compression_20261104_v63.py',
              ROOT/'research/notes/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63.md']
    result=dict(id=STEM,all_passed=True,
        both_rational_coordinate_compositions_exact_and_generate_entire_known_tau_fixed_source_field=True,
        independent_reverse_new_to_old_sparse_substitution_matches_all7_polynomials=True,
        all_h_denominators_and_xi_scalars_preserved=True,termcounts=expected,
        independent_h_identity_F4_five_terms_and_p0_factorization=True,
        exact_Q0_p0_p4_relation_and_B_recovery_on_displayed_dense_open=True,
        actual_NEW_target_field_equality_birationality_or_equation_claimed=False,
        no_new_point_equation_search_degree_job_GB_native_or_author_rerun=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','termcounts','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
