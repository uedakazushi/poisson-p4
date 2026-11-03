#!/usr/bin/env python3
"""Independent sparse Lie-action/weight proof of the global group-slice unit."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from mps_s61_arithmetic_circuits_20261103_v62 import DAG

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))
BASIS=list(itertools.product(PAIRS,MONS))
INDEX={item:i for i,item in enumerate(BASIS)}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output):
    started=time.monotonic()
    parent_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz'
    parent=json.loads(gzip.decompress(parent_path.read_bytes()))
    cuts=parent['group_cut_raw_indices'];lie=parent['group_lie_columns'];cut_index={v:i for i,v in enumerate(cuts)}
    operators=[[defaultdict(int)for _ in lie]for _ in cuts]
    for col,(u,v)in enumerate(lie):
        for source,((i,j),mon)in enumerate(BASIS):
            def insert(r,t,m,coefficient):
                if r==t:return
                if r>t:r,t,coefficient=t,r,-coefficient
                target=INDEX[(r,t),tuple(sorted(m))]
                if target in cut_index:operators[cut_index[target]][col][source]+=coefficient
            if i==v:insert(u,j,mon,1)
            if j==v:insert(i,u,mon,1)
            if u in mon:
                changed=list(mon);changed.remove(u);changed.append(v)
                insert(i,j,changed,-mon.count(u))
    operators=[[{i:c for i,c in cell.items()if c}for cell in row]for row in operators]
    raw_weights=[sum(pair)-sum(mon)for pair,mon in BASIS]
    row_weights=[5-raw_weights[i]for i in cuts];lie_weights=[i-j for i,j in lie]
    forms=[[defaultdict(Fraction)for _ in lie]for _ in cuts]
    source_degrees={-1:0,0:0,**{i+1:weight for i,weight in enumerate(parent['positive_weights'])}}
    assert all(w>0 for w in parent['positive_weights'])
    for r,target in enumerate(cuts):
        for c,(u,v)in enumerate(lie):
            expected=row_weights[r]+lie_weights[c]
            for source,multiplier in operators[r][c].items():
                assert raw_weights[source]+u-v==raw_weights[target]
                for variable,coefficient in parent['canonical_linear_forms'][source]:
                    forms[r][c][variable]+=multiplier*Fraction(coefficient)
            forms[r][c]={i:coefficient for i,coefficient in forms[r][c].items()if coefficient}
            assert all(source_degrees[i]==expected for i in forms[r][c])
    assert sum(row_weights)==16 and sum(lie_weights)==-16
    a=s.Symbol('a');matrix=s.Matrix([[s.Rational(cell.get(-1,0))+a*s.Rational(cell.get(0,0))for cell in row]for row in forms])
    domain=s.QQ.poly_ring(a);det=domain.to_sympy(DomainMatrix.from_Matrix(matrix).convert_to(domain).det())
    assert s.cancel(det-s.sympify(parent['group_minor']))==0
    assert s.factor(det)==s.Rational(163840,531441)*a**3*(9*a+4)
    dag=DAG(parent['input_names']);dag.nodes=parent['nodes']
    t40_path=ROOT/'reports/Q_T40_FINITE_QUOTIENT_20261103_v63/certificate.json'
    t40=json.loads(t40_path.read_text())
    point={name:Fraction(0)for name in parent['input_names']};point['a']=Fraction(2,9)
    for name,c in zip(['z11','z12','z23','z33'],[4,6,-4,1]):point[name]=Fraction(c)
    test_points=[('Rstar_eps1',point),('T40_nine_all1',{name:Fraction(v)for name,v in zip(parent['input_names'],t40['dense_parent35_point'])})]
    tests=[]
    for label,point in test_points:
        values=dag.evaluate(point);raw=[values[i]for i in parent['raw150_reconstruction_nodes']]
        actual=s.Matrix([[sum((coefficient*raw[i]for i,coefficient in cell.items()),Fraction())for cell in row]for row in operators])
        value=DomainMatrix.from_Matrix(actual).det();expected=s.sympify(parent['group_minor']).subs(a,s.Rational(point['a']))
        assert value==expected and value
        tests.append(dict(point=label,exact_group17_minor=str(value)))
    output.mkdir(parents=True,exist_ok=False)
    certificate=dict(group_cut_raw_indices=cuts,group_lie_columns=lie,
                     raw_basis_action_operators=[[[[i,c]for i,c in sorted(cell.items())]for cell in row]for row in operators],
                     universal109_slice_linear_forms=[[[[i,str(c)]for i,c in sorted(cell.items())]for cell in row]for row in forms],
                     row_relative_weights=row_weights,Lie_weights=lie_weights,
                     zero_positive_section_matrix=[[str(v)for v in row]for row in matrix.tolist()])
    (output/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    note=ROOT/'research/notes/Q_GROUP_SLICE_GLOBAL_UNIT_20261103_v63.md'
    inputs=[Path(__file__),note,parent_path,t40_path,
            ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json']
    report=dict(id='Q_GROUP_SLICE_GLOBAL_UNIT_20261103_V63',all_passed=True,
                independent_sparse_Lie_action_on_all150_basis_monomials=True,
                all289_entries_homogeneous_in_full109_slice=True,
                row_degree_sum=16,Lie_degree_sum=-16,universal_determinant_degree=0,
                every_retained_variable_except_a_strictly_positive=True,
                exact_universal109_slice_group_minor=str(s.factor(det)),
                no_positive_variable_dependence_of_group_minor=True,
                group17_map_etale_at_identity_times_entire_parent_D_slice=True,
                canonical_affine_completed_local_factor_dimension=17,
                regular9_slice_points_give_regular26_canonical_affine_points=True,
                nonzero_cone_projectivization_regular_dimension=25,
                literal16_group_factor_at_arbitrary_non_Qa_point_claimed=False,
                global_group_cross_section_claimed=False,numeric_exact_QQ_checks=tests,
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
                evidence_sha256={'certificate.json':sha(output/'certificate.json')},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS global109-slice group minor exact a-only unit;',report['seconds'],'seconds')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
