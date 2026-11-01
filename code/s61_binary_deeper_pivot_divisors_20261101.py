#!/usr/bin/env python3
"""Two deeper binary pivot divisors have a free65 open and are removed generically."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_binary_deeper_boundary_removal_20261101 import verify_point
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time()
    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    rank_path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    pivot_path=ROOT/'reports/S61_BINARY_DEEPER_RANK_SIX_EXCLUSION_20261101.json'
    family,rank,pivots=[json.loads(p.read_text()) for p in [family_path,rank_path,pivot_path]]
    assert all(r['all_passed'] for r in [family,rank,pivots])
    u=s.symbols('u1:62');uu=list(map(s.Symbol,family['sample']))
    a,b,c,d,e,f=[u[i] for i in [3,13,34,35,48,49]]
    K=5*c+f;Q=600*b*b-255*b*e-63*e*e+100*a*f
    R=d*(100*b*b-505*b*e+187*e*e)-40*f*f*(2*b-e)
    action=s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    stab=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    results=[]
    for saved in pivots['nonzero_minor_on_each_remaining_divisor']:
        label=saved['divisor']
        if label not in ['5b-11e','R']:continue
        sample=dict(family['sample']);sample.update(saved['sample'])
        point={s.Symbol(k):s.Rational(v) for k,v in sample.items()}
        factor=5*b-11*e if label=='5b-11e' else R
        cons=s.Matrix([u[11],u[12],K,Q,factor])
        assert cons.subs(point)==s.zeros(5,1) and point[f]!=0
        assert cons.jacobian(uu).subs(point).rank()==5
        gb=s.groebner(list(cons),*uu,domain=s.QQ)
        assert all(gb.reduce(s.expand(v))[1]==0 for v in cons.jacobian(uu)*stab)
        normal=(cons.jacobian(uu)*action).subs(point)
        columns=list(normal.to_DM().rref()[1]);assert len(columns)==5
        determinant=normal[:,columns].det();assert determinant!=0
        data=verify_point(sample)
        data.update(label=label,source_equations=list(map(str,cons)),source_open='u50 != 0',
            source_dimension=11,source_geometrically_irreducible=True,
            beta_normal_matrix=[[str(v) for v in row] for row in normal.tolist()],
            beta_minor_columns=columns,beta_dominance_minor=str(determinant),
            projective_stabilizer_preserves_source=True,joint_image_dimension=16,
            relative_binary_dimension=14,entire_general_beta_stratum_removed=True)
        results.append(data)
    assert len(results)==2
    paths=[Path(__file__),family_path,rank_path,pivot_path,
        ROOT/'code/s61_binary_deeper_boundary_removal_20261101.py',
        ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',
        ROOT/'research/notes/S61_BINARY_DEEPER_PIVOT_DIVISORS_20261101.md']
    return dict(id='S61_BINARY_DEEPER_PIVOT_DIVISORS_20261101',all_passed=True,prime=32003,cases=results,
        proof_note='research/notes/S61_BINARY_DEEPER_PIVOT_DIVISORS_20261101.md',
        limitations=['The two entire16D strata are removed only over sufficiently general beta.',
            'The factor5b-3e is separate: a sample has transverse length67 and is not claimed free65.',
            'No exhaustion of the remaining17D deeper stratum or the full affine fiber.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],cases=[{k:r[k] for k in ['label','beta_dominance_minor','relation_minor_mod_prime']} for r in report['cases']])))
