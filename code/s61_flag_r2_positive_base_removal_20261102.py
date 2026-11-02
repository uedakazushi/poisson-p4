#!/usr/bin/env python3
"""Prove removal of r=2 using a positive base variable and three transverse cuts."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import sympy as sp
from sympy.polys.rings import ring
from s61_model import BASIS, S61, X
from s61_flag_transverse_obstruction_20261101 import specialize
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
STEM = 'S61_FLAG_R2_POSITIVE_BASE_REMOVAL_20261102'
SPECIAL = 'S61_FLAG_REMAINING_CORE_WEIGHTED_POINTS_20261102'
ORDINARY = 'S61_FLAG_ORDINARY_WEIGHTED_FINITE_KERNEL_20261102'
SQUARE = 'S61_FLAG_BOUNDARY_SQUAREZERO_QUOTIENT_20261101'
CUT = [19,26,27]
WEIGHTS = [1,1,2,2,1,2,2]


def read(stem):
    data = json.loads((ROOT/'reports'/f'{stem}.json').read_text())
    assert data['all_passed']
    return data


def positive_rows(model, ordinary):
    rho = sp.Symbol('p')
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    base = dict.fromkeys(one,sp.S.Zero)
    base.update({4:1, 48:1, 12:2, 11:rho})
    _, equations = specialize(model,base,True)
    frame = read(SPECIAL)['points'][0]['frame']
    columns = frame['chart']
    fixed = frame['fixed']
    selected = frame['selected']
    assert fixed == 25
    tangent = sp.Matrix(frame['T']).applyfunc(sp.Rational)
    inverse = sp.Matrix(frame['U']).applyfunc(sp.Rational)
    jac = sp.Matrix([[-f.get(tuple(sorted((fixed,j))),0) for j in columns]
                     +[-f.get((fixed,),0)] for f in equations]).subs(rho,0)
    assert jac*tangent == sp.zeros(135,7)
    assert jac[selected,:]*inverse == sp.eye(36)
    assert tangent[columns.index(27),:] == sp.eye(7)[0,:]
    for i in CUT:
        assert not any(inverse[columns.index(i),:])
    ring_, *variables = ring('p,x1,x2,x3,x4,x5,x6',sp.QQ)
    coordinate = [sum((ring_(a)*v for a,v in zip(tangent.row(i)[1:],variables[1:]) if a),
                      ring_.zero) for i in range(45)]
    assert all(not coordinate[columns.index(i)] for i in CUT)
    rows = []
    for equation in equations:
        terms = []
        for mon,coefficient in equation.items():
            polynomial = sum((ring_(c)*variables[0]**k for (k,),c in sp.Poly(coefficient,rho).terms()),
                             ring_.zero)
            indices = tuple(columns.index(i) for i in mon if i != fixed)+(44,)*(2-len(mon))
            terms.append((polynomial*(-1)**mon.count(fixed),indices))
        rows.append(terms)

    def evaluate(terms):
        value = ring_.zero
        for coefficient,indices in terms:
            for index in indices:
                coefficient *= coordinate[index]
            value += coefficient
        return value

    row_weights = ordinary['row_weights']
    raw_weights = ordinary['raw_coordinate_weights']
    for weight in range(1,max(raw_weights)+1):
        residuals = [(j,evaluate(rows[i])) for j,i in enumerate(selected) if row_weights[i] == weight]
        for i,grade in enumerate(raw_weights):
            if grade == weight:
                for j,residual in residuals:
                    if inverse[i,j]:
                        coordinate[i] -= ring_(inverse[i,j])*residual
    equations = [evaluate(terms) for terms in rows]
    assert all(not equations[i] for i in selected)
    assert all(not coordinate[columns.index(i)] for i in CUT)
    for grade,polynomial in list(zip(row_weights,equations))+list(zip(raw_weights,coordinate)):
        assert all(sum(a*b for a,b in zip(mon,WEIGHTS)) == grade for mon in polynomial)
    return equations, coordinate[-1]


def m2_polynomial(polynomial):
    terms = []
    for monomial,coefficient in polynomial.items():
        factors = [f'({coefficient})']+[f'x{i}^{power}' for i,power in enumerate(monomial) if power]
        terms.append('*'.join(factors))
    return '+'.join(terms) or '0'


def check_central_ring(rows,homogenizer):
    executable = os.getenv('M2_BINARY') or shutil.which('M2-binary') or shutil.which('M2')
    if not executable:
        raise RuntimeError('Macaulay2 is required to replay the exact QQ ideal computations.')
    environment = os.environ.copy()
    runtime = ROOT/'work/m2-runtime/lib'
    if runtime.is_dir():
        environment['DYLD_LIBRARY_PATH'] = ':'.join([str(runtime),'/opt/homebrew/lib',
                                                     '/opt/homebrew/opt/givaro/lib',
                                                     environment.get('DYLD_LIBRARY_PATH','')])
    program = '\n'.join([
        'R=QQ[x0,x1,x2,x3,x4,x5,x6,Degrees=>{1,1,2,2,1,2,2}];',
        'I=ideal('+','.join(map(m2_polynomial,rows))+'); tt='+m2_polynomial(homogenizer)+';',
        'G=gb I; D=saturate(I,ideal(x0));',
        'print toString hilbertSeries(D/I,Reduce=>true);',
        'print toString hilbertSeries(R/D,Reduce=>true);',
        'print toString hilbertSeries(R/(I+ideal(x0)),Reduce=>true);',
        'assert(tt%G!=0);assert((tt^2)%G==0);assert(tt%(gb D)==0);',
        'A=QQ[y1,y2,y3,y4,y5,y6,MonomialOrder=>{Weights=>{1,2,2,1,2,2},GRevLex}];',
        'atOne=map(A,R,{1_A,y1,y2,y3,y4,y5,y6}); S=atOne I;',
        'N=ideal(y1,y4,y5,y6);',
        'S5=N^2+ideal(1872*y3-100*y5+1808*y6-285*y1-420*y4,1170*y2-530*y5+3592*y6-399*y1-588*y4);',
        'assert(S==S5);assert(degree S==5);assert((atOne tt)%(gb S)==0);',
        'print "POSITIVE CENTRAL RING PASS";'])
    with tempfile.TemporaryDirectory(prefix='s61-r2-positive-') as temp:
        script = Path(temp)/'verify.m2'
        script.write_text(program)
        completed = subprocess.run([executable,'--script',str(script)],cwd=temp,env=environment,
                                   capture_output=True,text=True,timeout=80)
    if completed.returncode or 'POSITIVE CENTRAL RING PASS' not in completed.stdout:
        raise RuntimeError(completed.stdout+completed.stderr)
    z = sp.Symbol('T')
    expected = [5*z**2+17*z**3+35*z**4+44*z**5+10*z**6,
                (1+2*z+2*z**2)/(1-z), 1+2*z+7*z**2+12*z**3+18*z**4+9*z**5]
    series = [line for line in completed.stdout.splitlines() if 'T' in line and '/' in line]
    assert len(series) == 3
    for expression,target in zip(series,expected):
        assert sp.cancel(sp.sympify(expression.replace('^','**'),locals={'T':z})-target) == 0
    return hashlib.sha256(program.encode()).hexdigest()


def geometry(model):
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    weights = {}
    torus = [0,1,-1,0,0]
    for i in one:
        pair,mon = BASIS[model.certificate['raw_free_indices'][i]]
        exponents = sp.Poly(mon,*X).monoms()[0]
        weights[i+1] = sum(torus[j] for j in pair)-sum(a*b for a,b in zip(torus,exponents))-1
    assert {i:w for i,w in weights.items() if w>0} == {4:1,12:1}
    data = read(SQUARE)
    names = sp.symbols('u1:62')
    rho = sp.Symbol('p')
    values = dict.fromkeys(names,0)
    values.update({names[4]:1,names[12]:2,names[48]:1,names[11]:rho})
    cut = sp.Matrix(data['actual_C_tangent']).applyfunc(sp.sympify).subs(values)[CUT,:]
    assert cut.det() == 8

    y,parameters,beta = primitive_beta(model)
    entries = sp.Matrix([beta[i,j].diff(y[k]) for (i,j),k in COORDS])
    basis = entries.jacobian(parameters)
    independent = list(basis.T.to_DM().rref()[1])
    beta_motion = sp.Matrix.hstack(*[basis[independent,:].inv()*infinitesimal_action(entries,g)[independent,:]
                                    for g in sp_basis()]).subs(values)
    change = coords()['Pi']
    y = sp.Matrix(y)
    tensor = sp.zeros(4)
    tensor[1,3] = -y[0]**2
    tensor -= tensor.T
    motion = []
    for generator in sp_basis():
        # Use pullback for both beta and the contravariant tensor.
        derivative = -generator*tensor-tensor*generator.T+sum(
            ((generator*y)[i]*tensor.diff(y[i]) for i in range(4)),sp.zeros(4))
        physical = sp.Matrix([sp.Poly(derivative[i,j],*y).coeff_monomial(y[a]*y[b])
                              for i,j in pairs for a,b in m2])
        motion.append((change*physical)[CUT,:])
    full = beta_motion[0,:].col_join(sp.Matrix.hstack(*motion))
    group_minor = full[:,[4,5,6,9]].det()
    assert group_minor == 8
    return weights,cut,group_minor


def verify():
    start = time.monotonic()
    model = S61()
    rows,t = positive_rows(model,read(ORDINARY))
    certificate = check_central_ring(rows,t)
    weights,cut,minor = geometry(model)
    paths = [Path(__file__),ROOT/'research/notes'/f'{STEM}.md',
             ROOT/'code/toric_closed_jet_coordinates_20261102.py',
             ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
             ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
             ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    for stem in [SQUARE,SPECIAL,ORDINARY,'S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102',
                 'S61_UNIVERSAL_WEIGHT6_20261031','S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101']:
        assert read(stem)['all_passed']
        paths.append(ROOT/'reports'/f'{stem}.json')
    return dict(id=STEM,all_passed=True,core=2,positive_base_coordinate='u12',
                fixed_zero_coordinate='u4',cut_q_zero_based=CUT,central_graded_weights=WEIGHTS,
                raw_rows=135,exact_residual_terms=sum(map(len,rows)),central_QQ_program_sha256=certificate,
                central_fiber_length_at_rho0=49,central_fiber_length_at_nonzero_rho=5,
                central_torsion_weight_dimensions=[0,0,5,17,35,44,10],central_torsion_length=111,
                relative_base_weights=weights,actual_section_cut_determinant=str(cut.det()),
                smooth_source_minor=str(minor),whole_generic_beta_r2_core_t_squared_zero=True,
                full_affine_exhaustion=False,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                                         for p in paths},seconds=time.monotonic()-start)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    result = verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print('PASS',result['central_torsion_length'],result['smooth_source_minor'],result['seconds'])
