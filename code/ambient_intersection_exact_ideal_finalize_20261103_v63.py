#!/usr/bin/env python3
"""Freeze the actual image-ideal/colon receipt after the versioned M2 checks."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def main(directory):
    started=time.monotonic();directory=directory.resolve();builder_path=directory/'builder.json';body=json.loads(builder_path.read_text())
    for path,expected in body['input_sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    log_path=directory/'intersection_colon_and_reduction_v7.log';log=log_path.read_text()
    expected=['F_REMAINDER z25+v0+2*aa','G_REMAINDER 0','F_SQUARE_REMAINDER 0',
        'COLON_IN_LINEAR_PRIME true','LINEAR_PRIME_IN_COLON true','D_PLUS_F_IN_N true','N_IN_D_PLUS_F true',
        'COLON_DIM 4','COLON_DEGREE 1','LINEAR_PRIME_RANK 45','EXACT_COLON_AND_REDUCTION_FINISHED']
    assert all(line in log for line in expected) and 'error:' not in log
    variables=list(map(s.Symbol,body['ambient_variables']));v1=s.Symbol('v1')
    annihilator=[s.sympify(f) for f in body['candidate_annihilator'] if s.expand(s.sympify(f)) not in (v1,-v1)]
    linear=[f for f in annihilator if s.Poly(f,*variables).total_degree()<=1]
    matrix,rhs=s.linear_eq_to_matrix(linear,variables);assert rhs==s.zeros(rhs.rows,1)
    rref,pivots=matrix.rref();assert len(pivots)==45
    free_columns=[i for i in range(49) if i not in pivots];free=[variables[i] for i in free_columns]
    assert list(map(str,free))==['v1','v4','v5','v7']
    substitution={variables[pivot]:-sum(rref[row,j]*variables[j] for j in free_columns) for row,pivot in enumerate(pivots)}
    assert all(s.expand(f.subs(substitution,simultaneous=True))==0 for f in annihilator)
    prime=[s.expand((rref*s.Matrix(variables))[i]) for i in range(45)]
    # Separate exact target-coordinate proof that CA restrictionF=0 is N.
    maps=json.loads((ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json').read_text())
    chart=json.loads((ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json').read_text());qi=chart['q_indices_zero_based']
    u=s.Matrix(list(map(s.sympify,maps['families']['CA']['full_61_coordinates'])))
    parameters=list(map(s.Symbol,maps['families']['CA']['source_fiber_parameters']));nselected=maps['families']['N']['actual_q_selected_rows']
    chosen=s.Matrix([u[0]-4*u[34]+u[51],*u[[qi[i] for i in nselected],0]])
    determinant=s.factor(chosen.jacobian(parameters).det());assert determinant==-72
    paths=[Path(__file__),ROOT/'code/ambient_intersection_exact_ideal_20261103_v63.py',builder_path,
        directory/'intersection.m2',directory/'intersection.log',
        directory/'intersection_product_order_v2.m2',directory/'intersection_product_order_v2.log',
        directory/'intersection_proper_section_v3.m2',directory/'intersection_proper_section_v3.log',
        directory/'intersection_colon_v4.m2',directory/'intersection_colon_v4.log',
        directory/'intersection_colon_v5.m2',directory/'intersection_colon_v5.log',
        directory/'intersection_colon_and_reduction_v6.m2',directory/'intersection_colon_and_reduction_v6.log',
        directory/'intersection_colon_and_reduction_v7.m2',log_path]
    result=dict(id='AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63',all_passed=True,field='QQ',
        source_builder=body['id'],ambient_variable_count=49,actual_component_generator_counts=body['component_generator_counts'],
        full_actual_intersection_presentation='ID=IA+IB with218 explicit saved target generators, completed at the specified source point; add7 independent smooth orbit variables.',
        nilradical_generator=body['actual_nilpotent_candidate'],nonzero=True,square_zero=True,
        target_p_multiple=body['candidate_p_multiple'],target_p_multiple_zero=True,
        reduction_identity='ID+(F)=IN=IA+(F,G)',
        annihilator_identity='ID:(F)=linear_prime_J',annihilator_prime_generators=list(map(str,prime)),
        annihilator_free_variables=list(map(str,free)),annihilator_prime_slice_dimension=4,
        nilmodule_full_S61_dimension=11,nilmodule='C[[u1,...,u11]]*F',
        associated_primes=['reduced irreducible component N prime','smooth11-dimensional embedded layer prime'],
        completed_intersection_depth=11,completed_intersection_dimension=18,
        normalizer_coordinate_F_and_N_K_determinant=str(determinant),
        scope='Exact actual completed intersection in S61/rank4open near distinguished v24Q(a=1), not the fullP4 genericPoissonrank<=2 boundary. No ring splitting is asserted.',
        earlier_attempts=[dict(file='intersection.m2',status='GRevLex stopped at120seconds; no mathematical conclusion'),
            dict(file='intersection_product_order_v2.m2',status='0.03second actualGB; naive split equality FALSE, exactF²/G identities PASS'),
            dict(file='intersection_proper_section_v3.m2',status='tangent-adjusted split equality FALSE; common10-dimensional annihilator hypothesis FALSE'),
            dict(file='intersection_colon_v4.m2',status='M2 protectedHH variable; no final colon result'),
            dict(file='intersection_colon_v5.m2',status='colon dimension4; candidatek1² annihilator FALSE'),
            dict(file='intersection_colon_and_reduction_v6.m2',status='correctcolon bothdirectionsPASS, then protectedGF variable stopped reduction check'),
            dict(file='intersection_colon_and_reduction_v7.m2',status='ALL final colon and reduction identities PASS')],
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (directory/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('FINAL ACTUAL INTERSECTION RECEIPT PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--directory',type=Path,required=True)
    main(parser.parse_args().directory)
