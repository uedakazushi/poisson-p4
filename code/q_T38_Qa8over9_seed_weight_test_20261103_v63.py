#!/usr/bin/env python3
"""Exact seven-column certificate: the fixed-frame Qa8/9 seed is uniquely zero-edge."""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from q_t01_membership_20261103_v63 import identity,multiply,plus,scale,PAIRS,MONS
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_T38_QA8OVER9_SEED_WEIGHT_TEST_20261103_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    frame=ROOT/'reports/Q_T38_QA8OVER9_COMPACT_FRAME_20261103_v63'
    edge=ROOT/'reports/Q_T38_QA8OVER9_EDGE_WEIGHT_TEST_20261103_v63'
    for folder in [frame,edge]:
        r=json.loads((folder/'report.json').read_text());assert r['all_passed']
        for name,want in r['input_sha256'].items():assert sha(ROOT/name)==want
        for name,want in r['evidence_sha256'].items():assert sha(folder/name)==want
    D=json.load(gzip.open(frame/'certificate.json.gz','rt'));E=json.loads((edge/'certificate.json').read_text())
    K=s.QQ.algebraic_field(s.sqrt(25-8*s.sqrt(5)),alias='theta');theta=K.unit
    def dec(v):return sum((K.convert(s.Rational(c))*theta**i for i,c in enumerate(v)),K.zero)
    def enc(v):
        c=list(reversed(v.to_list()));return list(map(str,c))+['0']*(4-len(c))
    S=[[dec(v)for v in row]for row in D['eigenframe_coefficient_S0']]
    I=identity(K.zero,K.one);U=plus(S,scale(I,-K.one));inv=I;power=I
    for k in range(1,5):power=multiply(power,U,K.zero);inv=plus(inv,scale(power,K.convert((-1)**k)))
    basis=list(itertools.product(PAIRS,MONS));position={v:i for i,v in enumerate(basis)}
    def push(a,b,p,q):
        answer=[K.zero]*150
        for i,j in PAIRS:
            coefficient=S[i][a]*S[j][b]-S[i][b]*S[j][a]
            for u in range(5):
                for v in range(5):
                    value=coefficient*inv[p][u]*inv[q][v]
                    if value:answer[position[(i,j),tuple(sorted((u,v)))]]+=value
        return answer
    cross=[push(*v)for v in [(0,1,3,3),(0,2,4,4),(1,2,3,4)]]
    scale0=dec(D['mu'])*dec(D['gamma']['1,2']);assert scale0
    rows=[i for i,W in enumerate(E['target_weights'])if W>=5];assert rows==[90,105,120,121,135,136,137,140]
    A=[[dec(v)for v in E['canonical150_by4_edge_pushforward'][i]]+
       [factor*cross[j][i]/scale0 for j,factor in enumerate([-2,2,4])]for i in rows]
    # The first seven actual coefficient equations already span all seven monomials.
    augmented=[row[:]+[K.one if i==j else K.zero for j in range(7)]for i,row in enumerate(A[:7])]
    determinant=K.one
    for col in range(7):
        pivot=next(i for i in range(col,7)if augmented[i][col])
        if pivot!=col:augmented[pivot],augmented[col]=augmented[col],augmented[pivot];determinant=-determinant
        value=augmented[col][col];determinant*=value;augmented[col]=[v/value for v in augmented[col]]
        for i in range(7):
            if i!=col:
                value=augmented[i][col]
                if value:augmented[i]=[a-value*b for a,b in zip(augmented[i],augmented[col])]
    inverse=[row[7:]for row in augmented];assert determinant
    for i in range(7):
        for j in range(7):assert sum((inverse[i][k]*A[k][j]for k in range(7)),K.zero)==(K.one if i==j else K.zero)
    columns=[enc(determinant*theta**i)for i in range(4)]
    norm=s.Matrix([[s.Rational(columns[j][i])for j in range(4)]for i in range(4)]).det();assert norm
    cert=dict(actual_raw_rows=rows,monomials=['c0','c1','c2','c3','c0*c3','c1*c2','c2*c3'],
        exact8_by7_constant_coefficient_matrix=[[enc(v)for v in row]for row in A],
        first7_equations_inverse=[[enc(v)for v in row]for row in inverse],
        selected7_raw_rows=rows[:7],exact7minor=enc(determinant),exact7minor_inverse=enc(K.one/determinant),exact7minor_norm=str(norm),
        inverse_direction='monomials7 = inverse7 * first7_actual_coefficient_equations',
        normalized_seed_edges='tau_j=c_j*t^(5-source_weight_j)',
        full_coefficient_equation_ideal_equals='(c0,c1,c2,c3)')
    (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',frame/'report.json',frame/'certificate.json.gz',edge/'report.json',edge/'certificate.json',
           ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json']
    result=dict(id=STEM,all_passed=True,all8_full_source_W_ge5_coefficient_equations_retained=True,
        exact_constant7_by7_unit_and_inverse=True,full_quadratic_equation_ideal_exactly_maximal_c0_c1_c2_c3=True,
        only_zero_edge_constant_seed_gives_exact_Qa8over9_limit_in_this_frame=True,
        no_Groebner_basis_needed_for_frozen_proof=True,
        scope='fixed frozen quartic frame and normalized same-level constant seed only',
        general_T38_membership_or_P2_branch_excluded=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('T38 EXACT FIXED-FRAME SEED IDEAL MAXIMAL / CONSTANT UNIT7 PASS',result['seconds'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
