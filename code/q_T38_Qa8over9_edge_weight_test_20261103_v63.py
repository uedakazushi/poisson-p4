#!/usr/bin/env python3
"""Exact bounded obstruction to aligned level-minus-one T38 edge limits."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from q_t01_membership_20261103_v63 import identity,multiply,plus,scale,wedge_coefficients,PAIRS,MONS

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_T38_QA8OVER9_EDGE_WEIGHT_TEST_20261103_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/Q_T38_QA8OVER9_COMPACT_FRAME_20261103_v63'
    receipt=json.loads((folder/'report.json').read_text());assert receipt['all_passed']
    for name,want in receipt['input_sha256'].items():assert sha(ROOT/name)==want
    assert sha(folder/'certificate.json.gz')==receipt['evidence_sha256']['certificate.json.gz']
    data=json.load(gzip.open(folder/'certificate.json.gz','rt'))
    K=s.QQ.algebraic_field(s.sqrt(25-8*s.sqrt(5)),alias='theta');theta=K.unit
    def dec(v):return sum((K.convert(s.Rational(c))*theta**i for i,c in enumerate(v)),K.zero)
    def enc(v):
        c=list(reversed(v.to_list()));return list(map(str,c))+['0']*(4-len(c))
    S=[[dec(v)for v in row]for row in data['eigenframe_coefficient_S0']]
    I=identity(K.zero,K.one);U=plus(S,scale(I,-K.one));inv=I;power=I
    for k in range(1,5):power=multiply(power,U,K.zero);inv=plus(inv,scale(power,K.convert((-1)**k)))
    assert multiply(S,inv,K.zero)==I
    basis=list(itertools.product(PAIRS,MONS));position={v:i for i,v in enumerate(basis)}
    def push(a,b,p,q):
        answer=[K.zero]*150
        for i,j in PAIRS:
            vector=S[i][a]*S[j][b]-S[i][b]*S[j][a]
            if vector:
                for u in range(5):
                    for v in range(5):
                        coefficient=vector*inv[p][u]*inv[q][v]
                        if coefficient:answer[position[(i,j),tuple(sorted((u,v)))]]+=coefficient
        div=[[K.zero]*5 for _ in range(5)]
        for ((i,j),(p,q)),value in zip(basis,answer):
            for a,b,c in [(i,j,value),(j,i,-value)]:
                if a==p:div[b][q]+=c
                if a==q:div[b][p]+=c
        assert div==[[K.zero]*5 for _ in range(5)]
        return answer
    records=[(0,1,2,4),(0,2,1,3),(1,3,4,4),(2,4,3,3)]
    E=[push(*v)for v in records]
    weights=[i+j-p-q for ((i,j),(p,q))in basis]
    source_weights=[a+b-p-q for a,b,p,q in records];assert source_weights==[-5,-2,-4,0]
    for W,poly in zip(weights,data['canonical150_arc']):assert all(power==5-W for power,coefficient in poly)
    for pair,record in zip([(0,3),(1,2),(2,3)],[(0,1,3,3),(0,2,4,4),(1,2,3,4)]):
        a,b,p,q=record;assert a+b-p-q==sum(source_weights[j]for j in pair)
    pivots={};selected=[]
    for index,W in enumerate(weights):
        if W<0:continue
        row=[E[k][index]for k in range(4)]
        for col,pivot in sorted(pivots.items()):
            coefficient=row[col]
            if coefficient:row=[a-coefficient*b for a,b in zip(row,pivot)]
        col=next((j for j,c in enumerate(row)if c),None)
        if col is not None:
            value=row[col];pivots[col]=[c/value for c in row];selected.append(index)
    assert len(pivots)==len(selected)==4
    M=[[E[j][i]for j in range(4)]for i in selected]
    determinant=K.zero
    for perm in itertools.permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j]for i in range(4)for j in range(i+1,4))
        determinant+=sign*s.prod(M[i][perm[i]]for i in range(4))
    assert determinant
    coeff=[enc(determinant*theta**i)for i in range(4)]
    norm=s.Matrix([[s.Rational(coeff[j][i])for j in range(4)]for i in range(4)]).det();assert norm
    certificate=dict(source_edge_basis_pair_monomial=records,source_weights=source_weights,
        target_weights=weights,canonical150_by4_edge_pushforward=[[enc(E[j][i])for j in range(4)]for i in range(150)],
        selected4_rows=selected,selected4_target_weights=[weights[i]for i in selected],
        exact4minor_theta_basis=enc(determinant),exact4minor_inverse=enc(K.one/determinant),exact4minor_norm=str(norm),
        Laurent_pushforward_rule='entry(row,j)=E(row,j)*t^(source_weight_j-target_weight_row)',
        aligned_source_edges='tau_j=c_j*t^(-1-source_weight_j)',
        distinct_levels_in_full_T38_source=dict(toric=5,linear=-1,cross=-7))
    (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',folder/'report.json',folder/'certificate.json.gz',
           ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json']
    result=dict(id=STEM,all_passed=True,exact_quartic_field_rank_of_all_W_nonnegative_edge_rows=4,
        all4_edge_tensors_canonical_divergence_zero=True,
        all4_edge_highest_weights=[max(W for W,c in zip(weights,E[j])if c)for j in range(4)],
        exact_nonzero4minor_and_inverse_and_QQ_norm_saved=True,
        aligned_level_minus1_edge_ansatz_regular_only_for_zero_amplitudes=True,
        toric_linear_cross_levels_distinct_and_cannot_cancel_poles_in_fixed_coefficient=True,
        scope='fixed frozen quartic frame and constant aligned edge amplitudes only',
        T38_closure_or_P2_actual_branch_excluded=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('T38 ALIGNED EDGE TEST RANK4 / NO MEMBERSHIP EXCLUSION PASS',result['seconds'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
