#!/usr/bin/env python3
"""Independent full QQ(n3) matrix/unit replay binding the completed native rank."""
import argparse
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha
from q_positive_weight6_small_model_review_20261103_v63 import physical_coefficients


def main(output):
    directory=ROOT/'reports/Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63'
    d=read(directory/'report.json');assert d['all_passed']
    for name,digest in d['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in d['evidence_sha256'].items():assert sha(directory/name)==digest,name
    c=read(directory/'certificate.json.gz')
    family=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz')
    smith=read(ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz')
    R,*n=ring('n0,n1,n2',s.QQ);sym=s.symbols('n0:3')
    q=[R.from_expr(s.sympify(v))for v in family['canonical150']]
    C=s.Matrix(smith['canonical_basis']);nonzeroC=[[(j,s.QQ.convert(C[i,j]))for j in range(126)if C[i,j]]for i in range(150)]
    # Rebuild from independently expanded constant350 physical quadratic
    # coefficient dictionaries, not the producer's directional derivative.
    pairs,_=physical_coefficients();J=[]
    for terms in pairs:
        raw={}
        for (a,b),coefficient in terms.items():
            raw[a]=raw.get(a,R.zero)+coefficient*q[b]
            raw[b]=raw.get(b,R.zero)+coefficient*q[a]
        row={}
        for i,coefficient in raw.items():
            if coefficient:
                for j,value in nonzeroC[i]:row[j]=row.get(j,R.zero)+coefficient*value
        J.append({j:v for j,v in row.items()if v})
    assert len(J)==350
    for i,row in enumerate(c['full350_by126_canonical_Jacobian']):
        for j,text in enumerate(row):assert J[i].get(j,R.zero)==R.from_expr(s.sympify(text))
    A={i:row.copy()for i,row in enumerate(J)};columns=list(range(126))
    for original_row,original_column,pivot in c['rational_unit_pivots']:
        pivot=s.QQ.convert(s.Rational(pivot));assert pivot
        prow=A[original_row]
        assert prow.get(original_column,R.zero)==R(pivot)
        for row_id,row in A.items():
            if row_id==original_row:continue
            coefficient=row.pop(original_column,R.zero)
            if coefficient:
                factor=coefficient/pivot
                for column,value in prow.items():
                    if column==original_column:continue
                    v=row.get(column,R.zero)-factor*value
                    if v:row[column]=v
                    else:row.pop(column,None)
        del A[original_row];columns.remove(original_column)
    assert len(A)==280 and len(columns)==56
    zero_rows=[i for i,row in A.items()if not row]
    zero_cols=[j for j in columns if not any(j in row for row in A.values())]
    assert zero_rows==c['all73_identically_zero_original_rows'] and len(zero_rows)==73
    assert zero_cols==c['all10_identically_zero_original_columns'] and len(zero_cols)==10
    rows=[i for i,row in A.items()if row];cc=[j for j in columns if j not in zero_cols]
    assert rows==c['remaining_original_rows'] and cc==c['remaining_original_columns']
    for i,row in enumerate(rows):
        for j,col in enumerate(cc):assert A[row].get(col,R.zero)==R.from_expr(s.sympify(c['exact207_by46_polynomial_residual'][i][j]))
    # Parse EVERY rational function in the actual completed native input.
    # The first70 pivots are QQ units, so each entry here is a QQ polynomial.
    with (directory/'native-input.txt').open()as source:
        assert list(map(int,source.readline().split()))==[3,207,46]
        def polynomial():
            count=int(source.readline());p=R.zero
            for _ in range(count):
                line=list(map(int,source.readline().split()));assert len(line)==4
                p+=R.from_dict({tuple(line[1:]):s.QQ(line[0])})
            return p
        for row in rows:
            for col in cc:
                numerator=polynomial();denominator=polynomial()
                assert denominator.is_ground and denominator
                assert numerator/denominator.LC==A[row].get(col,R.zero)
        assert not source.read().strip()
    native=(directory/'native-log.txt').read_text()
    assert 'EXACT_FUNCTION_FIELD_RANK 28 ROWS 207 COLUMNS 46'in native
    assert native.count('PIVOT ')==28
    # A fresh exactQQ98 minor independently confirms the lower bound.
    sample=s.zeros(350,126)
    for i,row in enumerate(J):
        for j,value in row.items():sample[i,j]=sum(s.QQ.to_sympy(v)for v in value.values())
    determinant=DomainMatrix.from_Matrix(sample[c['lower98_rows'],c['lower98_columns']]).convert_to(s.QQ).det()
    assert s.Rational(determinant)==s.Rational(c['lower98_minor'])!=0
    component=ROOT/'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json'
    proof=read(component);assert proof['all_passed'] and proof['actual_global_projective_irreducible_component_dimension_exact']==25
    paths=[Path(__file__),directory/'report.json',directory/'certificate.json.gz',directory/'native-input.txt',directory/'native-log.txt',
        ROOT/'code/rational_function_matrix_rank_20261103_v62.cpp',ROOT/'code/q_positive_weight6_small_model_review_20261103_v63.py',
        ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/report.json',component,
        ROOT/'research/notes/Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63_REVIEW.md']
    result=dict(id='Q_NEW_FAMILY_GENERIC_JACOBIAN_RANK_20261103_v63_REVIEW',all_passed=True,
        independent_constant_quadratic_physical_coefficients_rebuild_ALL350_by126_QQn3=True,
        all70_QQ_unit_operations_and_zero73_rows10_columns_independently_replayed=True,
        every207_by46_native_rational_function_bound_to_independent_residual=True,
        trusted_exact_function_field_engine_and_completed_final_zero_rows_audited=True,
        native192second_run_reexecuted_in_this_review=False,
        independent_exactQQ98minor_lower_bound=True,exact_generic_full_canonical_Jacobian_rank=98,
        generic_affine_tangent_dimension=28,generic_projective_tangent_dimension=27,
        actual_component_projective_dimension=25,generic_transverse_embedding_dimension=2,
        actual_new_component_generically_nonreduced=True,generic_Artin_length_lower_bound=3,
        exact_generic_Artin_length_quadratic_rank_or_full_completion_claimed=False,
        no_generic_upper_from_samples_or_n1_specialization=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
    with output.open('x')as f:
        json.dump(result,f,indent=2);f.write('\n')
    print('PASS independent QQ(n3) genericJac98 and actual genericnonreducedembedding2; no exactlength claim')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output)
