#!/usr/bin/env python3
"""Check integral coordinate maps and the actual d=0 pencil flag endpoint.

This is an arithmetic bridge certificate, not a generic-fiber exhaustion test.
Uses the published Cloud tool prefix through PATH and never edits old evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

import sympy as s
from s61_model import S61, vertical_rref
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_marked_sigma_core_axis_engine_20261102 import build_axis
from toric_closed_jet_coordinates_20261102 import coords
from weight1_structure import sp_basis, COORDS, J

ROOT = Path(__file__).resolve().parents[1]
P = 32003


def rational_integral(value):
    value = s.Rational(value)
    assert int(value.q) % P
    return value


def unit(value):
    value = rational_integral(value)
    assert int(value.p) % P
    return value


def polynomial_integral(expression):
    expression = s.cancel(expression)
    numerator, denominator = s.fraction(expression)
    symbols = sorted(expression.free_symbols, key=str)
    for part in [numerator, denominator]:
        coefficients = s.Poly(part, *symbols).coeffs() if symbols else [part]
        for value in coefficients:
            rational_integral(value)
    if not symbols:
        unit(denominator)


def verify(output):
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    model = S61()
    inputs = [Path(__file__), ROOT/'code/s61_model.py',
              ROOT/'code/s61_marked_sigma_core_axis_engine_20261102.py',
              ROOT/'code/s61_flag_transverse_obstruction_20261101.py']
    def load(name):
        path = ROOT/'reports'/f'{name}.json'
        inputs.append(path)
        data = json.loads(path.read_text())
        assert data['all_passed']
        for relative,expected in data.get('input_sha256',{}).items():
            dependency = ROOT/relative
            assert hashlib.sha256(dependency.read_bytes()).hexdigest() == expected, relative
        return data
    coeff_count = 0
    for row in model.reduced['equations']:
        for _, numerator, denominator in row:
            rational_integral(s.Rational(numerator, denominator))
            coeff_count += 1
    for row in model.data['eliminated_raw_variables']:
        for _, numerator, denominator in row['polynomial']:
            rational_integral(s.Rational(numerator, denominator))
            coeff_count += 1
    g6 = load('S61_UNIVERSAL_WEIGHT6_20261031')['universal_new_equation']
    for _, numerator, denominator in g6:
        rational_integral(s.Rational(numerator, denominator))
        coeff_count += 1
    V, vpiv = vertical_rref()
    for value in V:
        rational_integral(value)
    perspective = model.perspective_matrix()
    unit(perspective.det())
    for value in perspective.inv():
        rational_integral(value)
    coord = coords()
    for key in ['E', 'Pi']:
        for value in coord[key]:
            rational_integral(value)

    square = load('S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101')
    K = s.Matrix(square['kernel']).applyfunc(s.sympify)
    for value in K:
        polynomial_integral(value)
    u = s.symbols('u1:62')
    for value in K:
        den = s.denom(s.cancel(value))
        assert den.free_symbols <= {u[4]}
        if den.free_symbols:
            terms = s.Poly(den, u[4]).terms()
            assert len(terms) == 1
            unit(terms[0][1])
        else:
            unit(den)
    unit(s.Rational(80, 21)); unit(-4); unit(-1)

    true = load('S61_FLAG_TRUE_CORE_WEIGHTED_REMOVAL_20261102')['models'][0]
    frame = true['frame']
    unit(s.Rational(frame['td'])); unit(s.Rational(frame['ud']))
    T = s.Matrix(frame['T']).applyfunc(s.Rational)
    U = s.Matrix(frame['U']).applyfunc(s.Rational)
    for value in list(T)+list(U):
        rational_integral(value)
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    beta = dict.fromkeys(one,s.S.Zero)
    beta.update({4:1,48:1,12:-s.Rational(4,25)})
    _, forms = specialize(model,beta,True)
    chart = [i for i in range(45) if i != 25]
    jacobian = s.Matrix([[-f.get(tuple(sorted((25,j))),0) for j in chart]
                         +[-f.get((25,),0)] for f in forms])
    assert jacobian*T == s.zeros(135,8)
    assert jacobian[frame['selected'],:]*U == s.eye(35)
    witnesses = []
    for record in true['rank_witnesses']:
        unit(record['minor_mod'])
        assert record['rank'] == record['columns']-5
        witnesses.append({key:record[key] for key in ['degree','rank','columns','minor_mod']})
    unit(s.Rational(54,25)); unit(2)

    flag = load('S61_FLAG_UNIFORM_OPEN_20261101')
    H = s.Matrix(flag['flag']['beta_kernel']).applyfunc(s.Rational)
    db = s.Matrix(flag['frame_change']['beta_derivative']).applyfunc(s.sympify)
    stabilizer = s.Matrix(flag['frame_change']['stabilizer']).applyfunc(s.Rational)
    b = s.symbols('b0:15')
    point = s.zeros(15,1)
    point[0] = point[12] = point[5] = 1
    point[3],point[6],point[14] = s.Rational(4,5),3,2
    invariant = (point[0]*point[2]-point[3]*point[12]+point[3]**2)/point[12]**2
    assert invariant == -s.Rational(4,25)
    spatial, parameters, primitive = primitive_beta(model)
    vector = lambda f:s.Matrix([f[i,j].diff(spatial[k]) for (i,j),k in COORDS])
    embedding = vector(primitive).jacobian(parameters)
    basis_rows = list(embedding.T.to_DM().rref()[1])
    inverse_embedding = embedding[basis_rows,:].inv()
    group = s.eye(4)
    stages = []
    for column, choose in [(4,lambda v:-v[3]/v[0]),
                           (2,lambda v:(v[14]-v[5])/(3*v[0])),
                           (3,lambda v:v[5]/v[12]),
                           (0,lambda v:-v[6]/v[0])]:
        parameter = choose(point)
        rational_integral(parameter)
        A = sum((stabilizer[k,column]*sp_basis()[k] for k in range(10)),s.zeros(4))
        assert A*A == s.zeros(4)
        g = s.eye(4)+parameter*A
        assert g.T*J*g == J
        D = db[:,column].jacobian(b)
        assert D**4 == s.zeros(15)
        move = sum((parameter**k/s.factorial(k)*D**k for k in range(4)),s.zeros(15))
        moved = (move*point).applyfunc(s.cancel)
        tensor = primitive.subs(dict(zip(parameters,H*point)))
        pulled = g.T*tensor.subs(dict(zip(spatial,g*s.Matrix(spatial))),simultaneous=True)*g
        actual = (inverse_embedding*vector(pulled)[basis_rows,:]).applyfunc(s.cancel)
        assert actual == H*moved
        point = moved; group = group*g
        stages.append(dict(column=column,parameter=str(parameter)))
    assert point[0] == point[12] == 1 and point[2] == -s.Rational(4,25)
    assert all(point[i] == 0 for i in [1,3,5,6,14])
    for value in list(group)+list(group.inv()):
        rational_integral(value)

    endpoint_core = s.zeros(15,1)
    endpoint_core[3],endpoint_core[12],endpoint_core[5],endpoint_core[6],endpoint_core[14] = s.Rational(4,5),1,1,3,2
    old = load('S61_MARKED_SIGMA_RESIDUAL_FINITE_20261102')
    names, weights, polynomials, axis_frame = build_axis(endpoint_core,[9,19,27],old['frame']['q_weights'])
    for polynomial in polynomials:
        for value in polynomial.values():
            rational_integral(value)
    endpoint_beta = {j:value for j,value in zip(one,[0,*endpoint_core])}
    _, endpoint_forms = specialize(model,endpoint_beta,True)
    endpoint_jacobian = s.Matrix([[-f.get(tuple(sorted((25,j))),0) for j in chart]
                                  +[-f.get((25,),0)] for f in endpoint_forms])
    endpoint_minor = unit(endpoint_jacobian[axis_frame['selected_rows'],
                                           axis_frame['selected_columns']].det(method='domain-ge'))
    orbit = s.Matrix(flag['flag']['q_orbit_differential']).applyfunc(s.Rational)
    # q32 is a nonpositive base coordinate, not a fourth group cut.
    # Restore the three q-cuts together with the H15 normal u4.
    orbit_rows = [chart.index(j) for j in [9,19,27]]
    normal = s.Matrix([[0,0,0,0,-2*endpoint_core[5],0,
                        -endpoint_core[3]+2*endpoint_core[12],0,0,0]])
    orbit_block = orbit[orbit_rows,:].col_join(normal)
    orbit_columns = list(orbit_block.to_DM().rref()[1])
    assert len(orbit_columns) == 4
    orbit_minor = unit(orbit_block[:,orbit_columns].det())
    def program(prime):
        def expression(polynomial):
            terms = []
            for monomial, coefficient in polynomial.items():
                value = s.Rational(coefficient)
                value = int(value.p)*pow(int(value.q),-1,prime)%prime if prime else f'({value})'
                if value:
                    terms.append(str(value)+''.join('*'+name+('^'+str(power) if power>1 else '')
                                                  for name,power in zip(names,monomial) if power))
            return '+'.join(terms) or '0'
        field = f'(ZZ/{prime})' if prime else 'QQ'
        return ('R='+field+'['+','.join(names)+',Degrees=>{'+','.join(map(str,weights))+'}];\n'
                +'I=ideal('+','.join(expression(p) for p in polynomials if p)+');\n'
                +'P=ideal apply(0..12,i->R_i); A=saturate(I,ideal(rho)); U=I+ideal(rho);\n'
                +'assert(dim(R/I)==1 and dim(R/U)==0 and numColumns basis(R/U)==76);\n'
                +'assert(radical A==P); assert(numColumns basis(R/(A+ideal(rho)))==31);\n'
                +'assert(v12^2%I==0);print toString hilbertSeries(R/I,Reduce=>true);\n'
                +'print "D0 FLAG CENTRAL PASS"; exit 0;\n')
    m2 = shutil.which('M2')
    assert m2
    central = []
    for prime in [P,0]:
        path = output/f'endpoint-{prime}.m2'; path.write_text(program(prime))
        log = output/f'endpoint-{prime}.log'
        with log.open('w') as stream:
            subprocess.run([m2,'--script',str(path.resolve())], stdout=stream,stderr=subprocess.STDOUT,
                           check=True,timeout=120)
        assert 'D0 FLAG CENTRAL PASS' in log.read_text()
        central.append(dict(coefficient_field='QQ' if not prime else f'GF({prime})',
                            output=log.read_text(),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    report = dict(id='S61_BINARY_D_ZERO_PRIME_BRIDGE_20261103',all_passed=True,prime=P,
                  full_reconstruction_and_135_coefficients_checked=coeff_count,
                  vertical_rref_p_integral=True,perspective_determinant=str(perspective.det()),
                  E_and_Pi_p_integral=True,squarezero_kernel_denominators_only_p_units_and_u5=True,
                  true_core_frame_p_integral=True,true_core_JT_and_JU_rebuilt=True,
                  existing_true_core_minor_witnesses=witnesses,
                  true_core_minor_replay_scope='Stored nonzero minors and input hashes checked; native rank matrices not rerun.',
                  actual_endpoint_invariant=str(invariant),endpoint_Borel_steps=stages,
                  endpoint_Borel_matrix=[[str(x) for x in row] for row in group.tolist()],
                  endpoint_normal_form=list(map(str,point)),endpoint_axis_frame=axis_frame,
                  endpoint_Jacobian_minor=str(endpoint_minor),endpoint_orbit_minor=str(orbit_minor),
                  endpoint_orbit_columns=orbit_columns,
                  central_cases=central,
                  limitations=['This script checks arithmetic prerequisites and central endpoint equations.',
                               'Ambient nilpotence requires flat square-zero quotient and saturated-support proof.',
                               'No generic five-point exhaustion is claimed.'],
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in dict.fromkeys(inputs)},seconds=time.monotonic()-start)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=verify(args.output)
    print('PASS',result['id'],result['seconds'])
