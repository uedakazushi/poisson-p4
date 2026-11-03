#!/usr/bin/env python3
"""Independent full period/edge torus stabilizers and closed action-rank audit."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from s61_model import X
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replay',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    directory=ROOT/'reports/Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63'
    d=read(directory/'report.json');fresh=read(args.replay/'report.json')
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
    assert (directory/'certificate.json').read_bytes()==(args.replay/'certificate.json').read_bytes()
    for filename,digest in d['input_sha256'].items():assert sha(ROOT/filename)==digest,filename
    cert=read(directory/'certificate.json')
    catalog=read(ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json')
    pairs=list(itertools.combinations(range(1,5),2))
    cases=[];sourcepaths=[]
    for name in ['T01','T02','T03']:
        record=next(r for r in catalog['records']if r['component_id']==name)
        assert record['full_s61_component_dominance']
        diagram=record['source_diagram'];period=s.symbols(f'p0:{diagram["period_dimension"]}')
        six=[sum(p*s.Rational(row[k])for p,row in zip(period,diagram['period_basis']))for k in range(6)]
        dual=[six[5],-six[4],six[3],six[2],-six[1],six[0]]
        P=s.zeros(5)
        for (i,j),coefficient in zip(pairs,dual):P[i,j]=coefficient*X[i]*X[j];P[j,i]=-P[i,j]
        characters=[]
        for edge in diagram['edges']:
            i,j=edge['edge'];powers=[0]*5
            for k,power in edge['orders'].items():powers[int(k)]=int(power)
            assert sum(powers)==2
            character=[int(k==i)+int(k==j)-powers[k]for k in range(5)]
            assert sum(character)==0
            characters.append(character)
            term=s.Symbol('t0')*s.prod(z**power for z,power in zip(X,powers))
            P[i,j]+=term;P[j,i]-=term
        assert characters==cert[name]['diagonal_characters']
        if name=='T01':
            weights=[s.eye(5)[:,i]for i in range(1,5)]
        else:
            assert len(characters)==1
            C=s.Matrix(characters)
            weights=[s.Matrix([0,*list(v)])for v in C[:,1:5].nullspace()]
            source=ROOT/f'reports/MPS_S61_RAW_CATALOGUE_20261103_v62/{name}.json'
            sourcepaths.append(source);saved=read(source)
            assert saved['projective_Jacobi_identically_zero']
            local={str(z):z for z in [*X,*period,s.Symbol('t0')]}
            expected=s.zeros(5)
            for pair,text in saved['homogeneous_bivector'].items():
                i,j=map(int,pair.split(','));expected[i,j]=s.sympify(text,locals=local);expected[j,i]=-expected[i,j]
            assert (P-expected).applyfunc(s.expand)==s.zeros(5)
        assert len(weights)==(4 if name=='T01'else 3)
        W=s.Matrix.hstack(*weights)
        assert W.rank()==len(weights) and W[0,:]==s.zeros(1,len(weights))
        for weight in weights:
            A=s.diag(*weight)
            value=A*P+P*A-sum((weight[i]*X[i]*P.diff(X[i])for i in range(5)),s.zeros(5))
            assert value.applyfunc(s.expand)==s.zeros(5)
        upper=24-len(weights)
        assert upper==cert[name]['whole_component_projective_orbit_dimension_at_most']<22
        cases.append(dict(component=name,full_period_and_all_edge_variables_retained=True,
                          diagonal_characters=characters,independent_torus_dimension=len(weights),
                          whole_closed_component_projective_orbit_upper=upper))
    orbit=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json')
    assert orbit['all_passed'] and orbit['projective_PGL_orbit_dimension']==22
    paths=[Path(__file__),directory/'report.json',directory/'certificate.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json',
           ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json',
           ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
           ROOT/'code/q_smith_global_model_review_20261103_v63.py',ROOT/'code/s61_model.py',
           ROOT/'research/notes/Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63_REVIEW.md',*sourcepaths]
    result=dict(id='Q_POSITIVE_WEIGHT6_TORUS_EXCLUSION_20261103_v63_REVIEW',all_passed=True,
        fresh_producer_and_complete_character_certificate_match=True,
        independent_general_period_Pfaffian_dual_toric_tensor_and_actual_single_edge_sources=True,
        all_parameters_retained_and_all_diagonal_kernel_actions_zero=True,
        kernel_tori_descend_to_PGL_because_character_sum_zero=True,cases=cases,
        action_rank_upper_bound_extends_from_dense_PGL_source_to_whole_closed_component=True,
        positive_ray_actual_projective_orbit22_excludes_entire_T01_T02_T03=True,
        no_sample_rank_as_generic_upper_bound_or_other_MPS_exclusion_claimed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as f:
        json.dump(result,f,indent=2);f.write('\n')
    print('PASS full period/edge torus action and closed-rank exclusion of entire T01/T02/T03')


if __name__=='__main__':main()
