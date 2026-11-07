#!/usr/bin/env python3
"""Independent same-anchor fullD/actual joint-family smoothness review.

Only existing independent Fraction and sparse helpers are imported. No producer,
CAS, local ideal, new point or extra jet is used.
"""
from fractions import Fraction as F
from itertools import product
import hashlib,json,resource,signal,time
from pathlib import Path
import classification_quadratic_a_symplectic_review as S
import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_unimodular_logarithmic_review import RF
ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_QUADRATIC_A_SMOOTH_ANCHOR_20261107.json'
FRAME='reports/CLASSIFICATION_NIL22_QUADRATIC_A_WEIGHT_FRAME_20261107.json'
OUT='reports/CLASSIFICATION_QUADRATIC_A_SMOOTH_REVIEW_20261107.json'
JOB='reports/CLASSIFICATION_QUADRATIC_A_SMOOTH_REVIEW_JOB_20261107.json'


def at7(v):
    v=RF(v)
    def p(poly):
        out=F(0)
        for a in reversed(poly):out=out*7+a
        return out
    return p(v.p)/p(v.q)
def rf_der(v):
    v=RF(v);dp=tuple((i+1)*v.p[i+1] for i in range(max(0,len(v.p)-1)))
    dq=tuple((i+1)*v.q[i+1] for i in range(max(0,len(v.q)-1)))
    return RF(dp)/RF(v.q)-RF(v.p)*RF(dq)/(RF(v.q)*RF(v.q))
def numeric(poly):
    out={}
    for ex,c in poly.items():
        value=at7(c)
        for k,n in enumerate(ex[5:],start=5):value*=F(k==6)**n
        e=ex[:5]+(0,0,0);out[e]=out.get(e,F(0))+value
    return {e:c for e,c in out.items() if c}
def point(P):return q.antisymmetric({ij:numeric(v) for ij,v in P.items() if ij in q.PAIRS and v})
def pplus(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{}) for P in Ps)) for ij in q.PAIRS})
def vec(P):return q.vec(P)
def mm(A,B):return [[sum((a*b for a,b in zip(row,col) if a and b),F(0)) for col in zip(*B)] for row in A]


def main():
    started=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    dest=ROOT/OUT
    if dest.exists():raise FileExistsError(dest)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    data=json.loads((ROOT/INPUT).read_text());old=json.loads((ROOT/FRAME).read_text())['frame']
    original,_=S.original_source();whole=S.homogenize({ij:S.scale(v,2) for ij,v in original.items() if ij in S.PAIRS})
    source=json.loads((ROOT/S.ORIGINAL).read_text())['cases'][7]
    aliases={f'x{k}':S.a(1,k) for k in range(1,5)}
    saved=S.decode(source['universal_brackets'],aliases)
    assert whole==saved
    P=point(whole)
    legacy=q.decode(old['canonical_point']);assert vec(P)==vec(legacy)
    assert not any(q.jacobi(P).values())
    # True7 source derivatives at the SAME anchor, using I=m0+m4 chain rule.
    dparams=[]
    for index in (6,7,8,9,10,11):
        derivative={ij:S.der(v,index) for ij,v in whole.items() if ij in q.PAIRS}
        if index==10:
            derivative={ij:S.add(derivative[ij],S.der(whole.get(ij,{}),6)) for ij in derivative}
        dparams.append(point(derivative))
    dparams.append(point({ij:{ex:rf_der(c) for ex,c in v.items() if rf_der(c)} for ij,v in whole.items() if ij in q.PAIRS}))
    # Undo C=(q1,x,y,y-q2) on the already independently verified literal curve.
    arc=json.loads((ROOT/S.ARC).read_text())['curves'][1]
    family=S.decode(arc['actual_regular_affine_family'])
    Cback=[S.parse(t) for t in ['u1','v1-v2','u2','v1']]
    Cforward=[S.parse(t) for t in ['u1','v1','v2','v2-u2']]
    moved=S.push(family,Cback,Cforward,S.a(2))
    H=S.homogenize(moved)
    assert vec(point({ij:S.subst(v,{5:{}}) for ij,v in H.items() if ij in q.PAIRS}))==vec(P)
    R=point({ij:S.der(v,5) for ij,v in H.items() if ij in q.PAIRS})
    expected=q.decode(data['new_arc_direction']);assert vec(R)==vec(expected)
    DR=q.jvec({quad:q.add(q.jacobi(pplus(P,R))[quad],q.scale(q.jacobi(R)[quad],-1)) for quad in q.QUADS})
    assert not any(DR)
    # Full ambient350x150D independently reconstructed from elementary raw lifts.
    columns=[]
    for ij,ex in q.COLS:
        B=q.antisymmetric({ij:{ex+(0,0,0):F(1)}})
        columns.append(q.jvec(q.jacobi(pplus(P,B))))
    D=q.transpose(columns)
    weights=old['weights'];cw=[sum(a*b for a,b in zip(ex,weights))-sum(weights[k] for k in ij) for ij,ex in q.COLS]
    rw=[sum(a*b for a,b in zip(ex,weights))-sum(weights[k] for k in quad) for quad,ex in q.ROWS]
    assert all(not v or cw[j]==rw[i] for i,row in enumerate(D) for j,v in enumerate(row))
    ranks=[]
    for block in old['complete_blocks']:
        rr,cc=block['implicit_rows'],block['implicit_columns'];minor=[[D[i][j] for j in cc] for i in rr]
        assert len(rr)==len(cc)==block['D_rank']
        value=q.determinant(minor);assert value and value==F(block['minor_determinant'])
        inverse=[[F(v) for v in row] for row in block['implicit_inverse']]
        eye=[[F(i==j) for j in range(len(rr))] for i in range(len(rr))]
        assert mm(minor,inverse)==eye and mm(inverse,minor)==eye
        ranks.append(block['D_rank'])
    assert sum(ranks)==101
    # Actual group/gauge/scalar and all7 source directions, in exact saved order.
    X=[q.atom(1,i) for i in range(5)];actual=[]
    for i,j in product(range(5),repeat=2):
        field=[{} for _ in range(5)];field[i]=X[j]
        E=q.antisymmetric({(a,b):q.add(q.mul(X[a],field[b]),q.scale(q.mul(X[b],field[a]),-1)) for a,b in q.PAIRS})
        orbit=q.antisymmetric({(a,b):q.add(P.get((j,b),{}) if a==i else {},P.get((a,j),{}) if b==i else {},q.scale(q.mul(X[j],q.der(P.get((a,b),{}),i)),-1)) for a,b in q.PAIRS})
        actual.extend([vec(E),vec(orbit)])
    actual.append(vec(P));actual.extend(vec(v) for v in dparams)
    oldF=q.transpose(actual);newF=q.transpose(actual+[vec(R)])
    assert newF==[[F(v) for v in row] for row in data['complete_actual_matrix']]
    assert not any(any(row) for row in mm(D,newF))
    rr,cc=data['rank49_minor_rows'],data['rank49_minor_columns']
    assert len(rr)==len(cc)==49
    value=q.determinant([[newF[i][j] for j in cc] for i in rr]);assert value==F(data['rank49_minor_determinant']) and value
    # Lower rankD101 and lower actualrank49 now force exact ranks101/49.
    assert q.rank(oldF)[0]==48
    assert q.rank([row[:51] for row in newF])[0]==47
    assert q.rank([[row[j] for j in range(0,50,2)] for row in newF])[0]==24
    normals=old['true_normal_vectors'];assert len(normals)==1
    Nvec=[F(normals[0].get(str(i),'0')) for i in range(150)]
    assert not any(row[0] for row in mm(D,[[v] for v in Nvec]))
    assert q.rank([row+[v] for row,v in zip(oldF,Nvec)])[0]==49
    assert len(set(old['implicit_columns']))==101
    # Independent D-pivot columns plus kerD vectors F48+N give the full150
    # direct sum, proving the old102-dimensional cut transverse to F48.
    assert 150-101==49 and 49-24-1==24
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,
        same_anchor='lambda7,m=(1,0,0,0,0,0),universal_integer_scale2',
        all7_source_derivatives_chain_rule_verified=True,all150_point_and_curve_direction=True,
        full350_DR=True,full350_by150_D_reconstructed=True,
        independent_weight_blocks_nonzero_minors=ranks,block_inverse_both_products=True,
        full150_by59_actual_matrix=True,all_D_times_actual_columns_zero=True,
        actual49_nonzero_minor=str(value),rank_D=101,rank_oldF=48,rank_new_joint_actual=49,
        raw_embedding_dimension=49,actual_single_family_dimension_lower=49,
        Euler_rank=24,projective_dimension=24,smooth_unique_full_component=True,
        known_component='One of accepted MPS40; T number not identified',
        joint_family_proof='Apply inverse full6mag normalizer and inverse scalar to the all-lambda actual epsiloncurve. Near this same anchor alpha,beta,I units, the one jointfamily is regular; epsilon0 equals whole7source, yielding all7 tangent columns plusR simultaneously. Add actualPGL/scalar/Euler to obtain all59 columns.',
        old_true48cut_support_dimension=1,
        old_cut_full150_directsum_rechecked=True,
        old_cut_scope='Using accepted true48cut identification, its differentials have rank48 on the nowregular49dim raw germ; localcut dimension1. No exact oldideal thickness/model parametrization claimed.',
        no_producer_or_CAS_imports=True,no_new_point_or_local_ideal_run=True,
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)
if __name__=='__main__':main()
