#!/usr/bin/env python3
"""Rebuild all quadratic and cubic initial classes on the eight-variable flag chart.

The claims concern its function field. No claim about all core divisors,
fourth-order obstructions, or t-saturation is made by this certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports/S61_FLAG_BOUNDARY_EIGHT_CUBIC_20261101.json')
    args=parser.parse_args(); start=time.time()
    engines=[f'code/s61_flag_boundary_eight_{stage}_engine_20261101.py' for stage in ['frame','inverse','quadratic','cubic']]
    with tempfile.TemporaryDirectory(prefix='flag-eight-') as directory:
        wd=Path(directory)
        for engine in engines:
            subprocess.run([sys.executable,str(ROOT/engine),str(wd)],cwd=ROOT,check=True)
        frame=json.loads((wd/'flag-b1-normalized-frame.json').read_text())
        inverse=json.loads((wd/'flag-b1-triangular-inverse.json').read_text())
        quadratic=json.loads((wd/'flag-b1-universal-quadratic.json').read_text())
        cubic=json.loads((wd/'flag-b1-universal-cubic.json').read_text())
    assert len(frame['pivot_columns'])==36
    assert len(quadratic['pivot_columns'])==9
    assert cubic['rank']==3 and cubic['rref_symbols']==['b2','b8']
    assert cubic['pivot_columns']==[96,97,98]
    sparse=[]
    for row in cubic['cubic_rref']:
        sparse.append([{'monomial':cubic['monomials'][j],'coefficient':v} for j,v in enumerate(row) if v!='0'])
    dependencies=engines+[
        'code/s61_flag_boundary_eight_cubic_20261101.py',
        'code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py',
        'code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py',
        'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',
        'reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json']
    report=dict(all_passed=True,seconds=time.time()-start,
        scope='Function field QQ(b2,b4,b7,b8,b9,b10,b11,b13), b0=b12=1 and b1=b3=b5=b6=b14=0.',
        raw_equations=135,linear_rank=36,tangent_variables=9,
        quadratic_rank=9,cubic_new_initial_rank=3,
        hilbert_function_through_degree_three=[1,9,36,97],
        active_variables=[4,5,7,8],radical_variables=[0,1,2,3,6],
        quadratic_quotient_degree_three_dimension=100,
        quadratic_rref=quadratic['quadratic_rref'],quadratic_monomials=quadratic['monomials'],
        quadratic_pivot_rows=quadratic['pivot_rows'],quadratic_pivot_columns=quadratic['pivot_columns'],
        cubic_rref_sparse=sparse,
        cubic_rref_symbols=cubic['rref_symbols'],
        cubic_pivot_rows=cubic['pivot_rows'],
        denominators=dict(linear_minor=frame['minor'],tangent=frame['tangent_denominator'],
            right_inverse=quadratic['right_inverse_denominator'],quadratic_rref=cubic['quadratic_denominator'],
            quadratic_row_relations=cubic['row_relation_denominator']),
        inverse_dependency_order=inverse['dependency_order'],
        exact_checks=[
            'Acyclic triangular inverse multiplied by the original 36-square block equals identity.',
            'All135 quadratic rows equal their selected nine-row span over the full parameter field.',
            'The active quadratic multiplication form has rank at least two at b2=1.',
            'Every constant quadratic-row relation is applied to the projected cubic lift.',
            'All135 resulting cubic rows equal the selected three-row span over the full parameter field.'],
        limitations=[
            'Common row factors and core denominators are inverted; no all-point quadratic/cubic rank assertion follows.',
            'Reducing the second implicit lift modulo quadratics is valid for these cubic classes only.',
            'Unmodified fourth-order jets cannot use the reduced second lift.',
            'No formal ring, full initial ideal, or saturation removal is established.'],
        input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in dependencies})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('PASS eight-variable genuine cubic initial classes',report['seconds'],flush=True)

if __name__=='__main__':main()
