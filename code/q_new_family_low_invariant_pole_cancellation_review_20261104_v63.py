#!/usr/bin/env python3
"""Independent Fraction sparse remultiplication and H-remainder audit; no SymPy."""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63_REVIEW'
def read(p):
    b=p.read_bytes();return json.loads(gzip.decompress(b) if p.suffix=='.gz' else b)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
    r=read(p);assert r['all_passed']
    for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r
def mul(A,B):
    C={}
    for x,c in A.items():
        for y,d in B.items():
            z=tuple(a+b for a,b in zip(x,y));C[z]=C.get(z,Q(0))+c*d
    return {x:c for x,c in C.items() if c}
def enc(A):return [[str(c.numerator),str(c.denominator),*x] for x,c in sorted(A.items(),reverse=True)]

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63'
    old=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'
    ap=base/'report.json';author=accepted(ap)
    rp=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW/report.json';accepted(rp)
    H={(0,1,0,0):Q(147456),(0,0,2,0):Q(-9437184),(0,0,1,0):Q(6144),(0,0,0,0):Q(-1)}
    hp=[{(0,0,0,0):Q(1)}]
    for k in range(10):hp.append(mul(hp[-1],H))
    # Y=(3072W-1)^2/147456, univariate sparse W-powers.
    Y=[{(0,):Q(1)}];linearY={(2,):Q(9437184,147456),(1,):Q(-6144,147456),(0,):Q(1,147456)}
    for k in range(5):Y.append(mul(Y[-1],linearY))
    records=[];polys={};paths=[ap,rp];names=['F4','p0','p1','p2','p3','p4','N12']
    expect_k=[4,7,7,7,7,7,10];expect_n=[8,14,74,81,77,74,110];expect_B=[0,0,1,1,1,1,1]
    for i,n in enumerate(names):
        ip=old/(n+'.json.gz');qp=base/(n+'.json.gz');paths.extend([ip,qp]);original=read(ip);reduced=read(qp)
        P={tuple(map(int,v[1:])):Q(v[0]) for v in original['terms']}
        q={tuple(map(int,v[2:])):Q(int(v[0]),int(v[1])) for v in reduced['terms']}
        assert len(q)==len(reduced['terms'])==expect_n[i]
        assert all(c for c in q.values())
        k=reduced['cancelled_H_power'];assert k==expect_k[i] and mul(q,hp[k])==P
        assert reduced['denominator_H_power']+k==original['denominator_D_power']
        for field in ['denominator_B_power','denominator_integer','source_c_power']:assert reduced[field]==original[field]
        rem={}
        for (a,y,w,b),c in q.items():
            for (j,),v in Y[y].items():
                ex=(a,w+j,b);rem[ex]=rem.get(ex,Q(0))+c*v
        rem={x:c for x,c in rem.items() if c};assert rem
        mindeg=min(x[3] for x in q);maxdeg=max(x[3] for x in q)
        assert mindeg==0 and maxdeg==expect_B[i]
        assert reduced['source_c_power']==[2,4,4,4,4,4,6][i]
        assert author['results'][i]['name']==n and author['results'][i]['cancelled_H_power']==k
        assert author['results'][i]['reduced_terms']==len(q)
        records.append(dict(name=n,maximal_H_power=k,remaining_H_denominator=reduced['denominator_H_power'],
            original_terms=len(P),reduced_terms=len(q),B_degree=maxdeg,min_B_degree=mindeg,
            independent_H_remainder_variables=['a','W','B'],independent_nonzero_H_remainder=enc(rem)))
        polys[n]=(q,reduced)
    qF,df=polys['F4'];assert df['denominator_H_power']==0 and df['denominator_integer']=='3240' and df['denominator_B_power']==2
    dF=Q(df['denominator_integer']);ratios=[]
    for n in names[1:6]:
        q,meta=polys[n];shift=4-meta['denominator_B_power'];assert shift in [0,1]
        A={(a,y,w,b+shift):c for (a,y,w,b),c in q.items()}
        assert max(x[3] for x in A)==1
        ratios.append(dict(name=n,constant=str(dF**2/Q(meta['denominator_integer'])),
            denominator='H*qF^2',B_degree=1,numerator_A=enc(A)))
    qN,dn=polys['N12'];assert dn['denominator_H_power']==2 and dn['denominator_B_power']==5
    N={(a,y,w,b+1):c for (a,y,w,b),c in qN.items()};assert max(x[3] for x in N)==2
    certificate=dict(variables=['a=u/X','Y=X^2','W','B'],H=enc(H),qF=enc(qF),all7_exact_sparse_checks=records,
        normalized_p_ratios=ratios,normalized_N12_ratio=dict(constant=str(dF**3/Q(dn['denominator_integer'])),
            denominator='H^2*qF^3',factorization='B*(qN0+B*qN1)',B_degree=2,numerator=enc(N)),
        normalization_domain='original source domain intersect qF!=0; no whole H=0 chart claim')
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(certificate,indent=2)+'\n')
    paths += [Path(__file__),ROOT/f'research/notes/{STEM}.md',ROOT/'research/notes/Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63.md',
              ROOT/'code/q_new_family_low_invariant_pole_cancellation_20261104_v63.py']
    report=dict(id=STEM,all_passed=True,
        independent_Python_Fraction_sparse_engine_no_SymPy_producer_or_division=True,
        all7_entire_original_integer_polynomials_exactly_remultiplied=True,
        all7_maximal_H_valuations_by_independent_nonzero_Y_substitution_remainders=True,
        cancellation_H_powers=expect_k,reduced_termcounts=expect_n,reduced_B_degrees=expect_B,
        all_min_B_exponents_zero=True,normalized_all5_p_over_F4_squared_B_linear=True,
        normalized_N12_over_F4_cubed_B_times_B_linear=True,
        constants_and_c_B_H_denominator_metadata_bound_to_reviewed_source=True,
        no_new_point_relation_equation_degree_bound_or_target_field_index_claim=True,
        no_source26_higher_degree_probe_GB_native_or_author_rerun=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','cancellation_H_powers','reduced_termcounts','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
