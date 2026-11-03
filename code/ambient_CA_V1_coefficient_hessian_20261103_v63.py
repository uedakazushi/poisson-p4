#!/usr/bin/env python3
"""Numeric ordered-exponential Lie Hessian and univariate QQ source jets.

No multivariable physical differentiation is used. The source order is
20 ordered GL exponentials followed by9 linear fixed-fiber coordinates.
"""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))
INDEX={(pair,mon):i for i,(pair,mon) in enumerate(itertools.product(PAIRS,MONS))}


def operator(a,b):
    entries={}
    def put(rowpair,mon,col,value):
        i,j=rowpair
        if i==j:return
        if i>j:i,j,value=j,i,-value
        row=INDEX[((i,j),tuple(sorted(mon)))];key=(row,col)
        entries[key]=entries.get(key,0)+value
    for col,((i,j),mon) in enumerate(itertools.product(PAIRS,MONS)):
        if i==b:put((a,j),mon,col,1)
        if j==b:put((i,a),mon,col,1)
        if a in mon:
            other=list(mon);other.remove(a);other.append(b)
            put((i,j),other,col,-mon.count(a))
    return [(row,col,value) for (row,col),value in entries.items() if value]


def apply(L,vector):
    result=[Q(0)]*150
    for row,col,value in L:
        if vector[col]:result[row]+=value*vector[col]
    return result


def scaleadd(target,vector,scalar):
    if scalar:
        for i,value in enumerate(vector):
            if value:target[i]+=scalar*value


def polynomial_product(a,b,N):
    result=[Q(0)]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b[:N+1-i]):
                if y:result[i+j]+=x*y
    return result


def evaluator(point150,basis9_columns,group_GL_indices,degree):
    """Return q29-series ->150 coefficient series for orderedexp source.

    The group is exp(q0*E0)...exp(q19*E19), so coefficient-Lie
    operators apply in reverse order. Diagonal exponentials are retained.
    All input source coordinates have zero constant term.
    """
    N=degree;point=list(map(Q,point150));basis=[list(map(Q,v)) for v in basis9_columns]
    ops=[operator(*divmod(i,5)) for i in group_GL_indices]
    def run(q):
        assert len(q)==29 and all(Q(row[0])==0 for row in q)
        q=[[Q(v) for v in row[:N+1]]+[Q(0)]*max(0,N+1-len(row)) for row in q]
        current=[[Q(0)]*150 for _ in range(N+1)];current[0]=point.copy()
        for j,vector in enumerate(basis):
            for d in range(1,N+1):scaleadd(current[d],vector,q[20+j][d])
        for j in range(19,-1,-1):
            source=q[j];minimum=next((i for i,v in enumerate(source) if v),None)
            if minimum is None:continue
            result=[row.copy() for row in current]
            differentiated=[row.copy() for row in current];power=[Q(1)]+[Q(0)]*N
            for k in range(1,N//minimum+1):
                limit=N-k*minimum
                differentiated=[apply(ops[j],row) for row in differentiated[:limit+1]]
                if not any(any(row) for row in differentiated):break
                power=[v/k for v in polynomial_product(power,source,N)]
                for a in range(k*minimum,N+1):
                    if power[a]:
                        for b in range(min(limit,N-a)+1):scaleadd(result[a+b],differentiated[b],power[a])
            current=result
        return [[current[d][i] for d in range(N+1)] for i in range(150)]
    return run


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def main(checkpoint,output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    with gzip.open(checkpoint,'rt') as f:data=json.load(f)
    assert data['all_passed'] and data['checkpoint_only_no_cusp_order_claim']
    point=list(map(Q,data['point150']));C=[[Q(v) for v in row] for row in data['fixedCA150_by20']]
    indices=data['group20_GL_indices'];extra=data['extra9_fiber_indices']
    F=[[row[j] for row in C] for j in extra];ops=[operator(*divmod(i,5)) for i in indices]
    LP=[apply(L,point) for L in ops]
    J=LP+F
    assert [[J[j][i] for j in range(29)] for i in range(150)]==[[Q(v) for v in row] for row in data['full150_by29_target_Jacobian']]
    kappa=list(map(Q,data['kernel29']));linear=[Q(0)]*150
    for c,column in zip(kappa,J):scaleadd(linear,column,c)
    assert not any(linear)
    LF=[[apply(L,vector) for vector in F] for L in ops]
    LL={(i,j):apply(ops[i],LP[j]) for i in range(20) for j in range(i,20)}
    H=[]
    for j in range(20):
        column=[Q(0)]*150
        for i in range(20):scaleadd(column,LL[min(i,j),max(i,j)],kappa[i])
        for f in range(9):scaleadd(column,LF[j][f],kappa[20+f])
        H.append(column)
    for f in range(9):
        column=[Q(0)]*150
        for i in range(20):scaleadd(column,LF[i][f],kappa[i])
        H.append(column)
    Hkk=[Q(0)]*150
    for c,column in zip(kappa,H):scaleadd(Hkk,column,c)
    run=evaluator(point,F,indices,2)
    base=run([[Q(0),c,Q(0)] for c in kappa]);assert all(row[1]==0 for row in base)
    assert [row[2]*2 for row in base]==Hkk
    # Full29 polarization checks every mixed Hessian column, rather than
    # checking only its kernel-kernel contraction.
    for j in range(29):
        unit=[Q(int(i==j)) for i in range(29)]
        alone=run([[Q(0),c,Q(0)] for c in unit])
        mixed=run([[Q(0),a+b,Q(0)] for a,b in zip(kappa,unit)])
        assert [mixed[i][2]-base[i][2]-alone[i][2] for i in range(150)]==H[j]
    # All diagonal operators throughdegree7 match their closed weight formula.
    run7=evaluator(point,F,indices,7);diagonal_checks=[]
    for j,index in enumerate(indices):
        a,b=divmod(index,5)
        if a!=b:continue
        series=[Q(0),Q(1),Q(1,3)]+[Q(0)]*5
        q=[[Q(0)]*8 for _ in range(29)];q[j]=series
        got=run7(q)
        for row,((i,k),mon) in enumerate(itertools.product(PAIRS,MONS)):
            weight=int(i==a)+int(k==a)-mon.count(a)
            powers=[Q(1)]+[Q(0)]*7;want=powers.copy()
            for degree in range(1,8):
                powers=[v*weight/degree for v in polynomial_product(powers,series,7)]
                want=[x+y for x,y in zip(want,powers)]
            assert got[row]==[point[row]*v for v in want]
        diagonal_checks.append(index)
    payload=dict(id='AMBIENT_CA_V1_COEFFICIENT_HESSIAN_20261103_v63',all_passed=True,
        source_orderedexp=True,source_order='group20_GL_indices then extra9 fixedfiber indices',
        group20_GL_indices=indices,extra9_fiber_indices=extra,
        point150=list(map(str,point)),J_original150_by29=[[str(J[j][i]) for j in range(29)] for i in range(150)],
        kernel29=list(map(str,kappa)),complement28_columns=data['complement28_columns'],
        H_kappa_original150_by29=[[str(H[j][i]) for j in range(29)] for i in range(150)],
        kernel_kernel_quadratic150=[str(c/2) for c in Hkk],
        ordered_group_group_second_derivative='L_min(i,j) L_max(i,j) P, including diagonal L_i^2 P',
        group_fiber_second_derivative='L_i F_j',all29_numeric_polarization_checks=True,
        diagonal_degree7_closed_weight_checks=diagonal_checks,
        no_symbolic_multivariate_Hessian=True,no57block_or_cusp_order_claim=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),checkpoint]},seconds=time.monotonic()-started)
    target=output/'hessian.json.gz'
    with target.open('wb') as f:f.write(gzip.compress((json.dumps(payload,separators=(',',':'))+'\n').encode(),mtime=0))
    report={k:v for k,v in payload.items() if k not in ['point150','J_original150_by29','kernel29','H_kappa_original150_by29','kernel_kernel_quadratic150']}
    report['evidence_sha256']={'hessian.json.gz':sha(target)}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','all29_numeric_polarization_checks','no57block_or_cusp_order_claim','seconds']},indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();main(args.checkpoint.resolve(),args.output.resolve())
