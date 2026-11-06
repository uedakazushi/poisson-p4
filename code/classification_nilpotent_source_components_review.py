#!/usr/bin/env python3
"""Independent QQ ambient ranks for the four nilpotent modular sources.

Uses the separately authored sparse coefficient-loop engine from the mixed
quadric review; imports no arithmetic from either nilpotent producer.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

import sympy as S

import classification_rank3_mixed_component_review as independent

ROOT=Path(__file__).resolve().parents[1]


def sparse(matrix, coordinates):
    return {ij:dict(S.Poly(matrix[ij], *coordinates).terms())
            for ij in independent.PAIRS if matrix[ij] != 0}


def main(output):
    began=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS, (1536*1024**2, 1536*1024**2))
    signal.alarm(90)
    source_path=ROOT/'reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT_20261107.json'
    producer_path=ROOT/'reports/CLASSIFICATION_NILPOTENT_SOURCE_COMPONENTS_20261107.json'
    source=json.loads(source_path.read_text())
    producer=json.loads(producer_path.read_text())
    h, x1, x2, x3, x4=coords=S.symbols('h x1 x2 x3 x4')
    nn=S.symbols('n0:9')
    env={str(v):v for v in (*coords, *nn)}
    base=[S.sympify(f, locals=env) for f in source['charts'][0]['all_six_brackets']]
    results={}
    for name in ('s_1_3', 's_1_7', 's_1_5', 's_1_15'):
        branch=source['modular_zero_chart']['five_full_sources'][name]
        substitution={env[k]:S.sympify(v, locals=env)
                      for k,v in branch['all_point_substitution'].items()}
        parameters=[env[v] for v in branch['free_parameters']]
        point={v:int(str(v)==branch['rank4_iff_nonzero']) for v in parameters}
        family=S.zeros(5)
        # Homogenize by h^2 f(x/h), independently of the producer's term loop.
        for ij, entry in zip(itertools.combinations(range(1,5),2),base):
            affine=S.expand(entry.subs(substitution, simultaneous=True))
            hom=S.expand(h*h*affine.subs({v:v/h for v in coords[1:]}, simultaneous=True))
            family[ij]=hom
            family[ij[::-1]]=-hom
        actual=family.subs(point).applyfunc(S.expand)
        scale=S.ilcm(*[c.q for ij in independent.PAIRS
                       for _,c in S.Poly(actual[ij], *coords).terms()])
        family*=scale
        actual*=scale
        assert str(scale)==producer['results'][name]['integer_scale']
        for ij in independent.PAIRS:
            saved=producer['results'][name]['full_witness_brackets'].get(f'{ij[0]},{ij[1]}','0')
            assert S.expand(actual[ij]-S.sympify(saved, locals=env))==0
        universal_raw=independent.jacobi(sparse(family,coords),raw=True)
        assert all(S.expand(c)==0 for polynomial in universal_raw.values() for c in polynomial.values())
        raw=sparse(actual,coords)
        assert all(not q for q in independent.jacobi(raw,raw=True).values())
        assert independent.jacobi(raw)==S.zeros(350,1)

        columns=[]
        for ij, exponent in independent.COLUMNS:
            perturbed={pair:dict(coeffs) for pair,coeffs in raw.items()}
            perturbed.setdefault(ij,{})[exponent]=perturbed.get(ij,{}).get(exponent,0)+1
            columns.append(independent.jacobi(perturbed))
        derivative=S.Matrix.hstack(*columns)
        vertical=[]
        orbit=[]
        ee=S.Matrix(coords)
        for row,col in itertools.product(range(5),repeat=2):
            linear=S.zeros(5)
            linear[row,col]=1
            field=linear*ee
            vertical.append(independent.vector(sparse(ee*field.T-field*ee.T,coords)))
            tangent=linear*actual+actual*linear.T-coords[col]*actual.diff(coords[row])
            orbit.append(independent.vector(sparse(tangent,coords)))
        vm=S.Matrix.hstack(*vertical)
        group=S.Matrix.hstack(*vertical,*orbit,independent.vector(raw))
        directions=S.Matrix.hstack(*[independent.vector(sparse(family.diff(v).subs(point),coords))
                                    for v in parameters])
        full=group.row_join(directions)
        assert derivative*full==S.zeros(350,full.cols)
        old=producer['results'][name]
        certificates={label:independent.minor_certificate(mat,rank)
                      for label,mat,rank in
                      [('D',derivative,103),('vertical',vm,24),
                       ('group_vertical_scale',group,46 if name=='s_1_7' else 47),
                       ('actual_source_group_vertical_scale',full,47)]}
        assert old['full_derivative']['rank']==103
        assert old['source_orbit_vertical_scale']['rank']==47
        assert old['orbit_vertical_scale']['rank']==certificates['group_vertical_scale']['rank_QQ']
        pf=S.factor(actual[1,2]*actual[3,4]-actual[1,3]*actual[2,4]+actual[1,4]*actual[2,3])
        assert pf!=0
        results[name]=dict(all_passed=True, universal_full_raw_jacobi_zero=True,
            full350_projective_jacobi_zero=True, all55_actual_columns_annihilated=True,
            source_parameter_count=len(parameters), source_point={str(k):v for k,v in point.items()},
            integer_scale=str(scale), certificates=certificates,
            source_derivatives_sparse=[{str(i):str(directions[i,j]) for i in range(150) if directions[i,j]}
                                       for j in range(directions.cols)],
            matrix_sha256={label:hashlib.sha256(json.dumps(mat.tolist(),default=str).encode()).hexdigest()
                           for label,mat in [('D',derivative),('group',group),('source',directions)]},
            full_kernel_dimension=47, projective_tangent_dimension=22,
            projective_actual_family_rank=22, witness_smooth=True,
            divisor=str(h*pf), full_witness_brackets={f'{i},{j}':str(actual[i,j])
                                                    for i,j in independent.PAIRS if actual[i,j]!=0})
        print(json.dumps(dict(source=name, D=103,actual_source=47,
                              group=certificates['group_vertical_scale']['rank_QQ'])),flush=True)
    paths=[Path(__file__).resolve(),Path(independent.__file__).resolve(),source_path,producer_path]
    report=dict(all_passed=True,task_id='E_multiple_hyperplane',
        producer_arithmetic_imported=False, independent_sparse_engine_reused=True,
        all150_ambient_variables=True, all350_projective_equations=True,
        exact_field='QQ', results=results,
        independent_source_surjectivity_or_distinctness_claimed=False,
        smooth_full_component_dimension_for_each_source=22,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.monotonic()-began)
    with output.open('x') as handle:
        json.dump(report,handle,indent=2)
        handle.write('\n')
    signal.alarm(0)
    print(json.dumps(dict(all_passed=True,sources=4,seconds=report['seconds'])),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
