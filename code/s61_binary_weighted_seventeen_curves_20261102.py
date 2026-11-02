#!/usr/bin/env python3
"""Seventeen explicit deeper-core curves: raw weighted minors and same-C lower bound.

The supplied sparse polynomials define the curves; no completeness of their
list is an input or conclusion. Large implicit data are rebuilt and hashed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61, BASIS, X
from s61_flag_transverse_obstruction_20261101 import specialize
import s61_binary_deeper_core_artifacts_20261101 as geometry
ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'research/notes/S61_BINARY_WEIGHTED_SEVENTEEN_CURVES_20261102.md'
DATA=ROOT/'code/data/s61_binary_weighted_seventeen_curves_20261102.json'
NW=[1,1,1,1,1,2,1,2,1]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def prepare(wd):
    model=S61()
    chart=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    quad_path=ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
    data=json.loads(chart.read_text())['linear']
    R,*_=ring(','.join(data['variables']),s.QQ);field=R.to_field()
    fractions=lambda rows:[[field.from_expr(s.sympify(v)) for v in row] for row in rows]
    full_T=fractions(data['T']);full_U=fractions(data['U'])
    den_U=R.one
    for row in full_U:
        for v in row:den_U=den_U.lcm(v.denom)
    den_T=[]
    for j in data['normal']:
        den_j=R.one
        for row in full_T:den_j=den_j.lcm(row[j].denom)
        den_T.append(den_j)
    Q=dict(Tn=[[row[j].numer*den_T[k].exquo(row[j].denom) for k,j in enumerate(data['normal'])] for row in full_T],
           Un=[[v.numer*den_U.exquo(v.denom) for v in row] for row in full_U],
           Uden=den_U,Tden=den_T)
    S,e,d=ring('e,d',s.QQ);z=S.zero
    core=lambda p:S.from_dict({m[:2]:c for m,c in p.items() if not any(m[2:])})
    pol=lambda v:core(R.from_expr(s.sympify(v)))
    T=[[core(v) for v in row] for row in Q['Tn']]
    U=[[core(v) for v in row] for row in Q['Un']]
    J=[[pol(v) for v in row] for row in data['J']]
    den=core(Q['Uden']);dens=[core(v) for v in Q['Tden']]
    known=[S.from_expr(s.sympify(v)) for v in json.loads(quad_path.read_text())['irreducible_minor_factors']]
    for value in [den,*dens]:
        assert value
        for factor in known:
            while value.rem(factor)==0:value=value.exquo(factor)
        assert max(map(sum,value))==0
    sel=data['pivot_rows']
    for i,row in enumerate(J):
        for k in range(9):assert sum((v*T[j][k] for j,v in enumerate(row) if v),z)==0
    for a,i in enumerate(sel):
        for b in range(33):assert sum((v*U[j][b] for j,v in enumerate(J[i]) if v),z)==(den if a==b else z)
    chart_basis=s.Matrix(data['chart_basis']).applyfunc(s.sympify)
    assert not chart_basis.free_symbols
    chart_inverse=chart_basis.inv()
    inverse=[[S.domain.convert(v) for v in row] for row in chart_inverse.tolist()]
    Tc=[[sum((a*T[j][k] for j,a in enumerate(row) if a),z) for k in range(9)] for row in inverse]
    Uc=[[sum((a*U[j][k] for j,a in enumerate(row) if a),z) for k in range(33)] for row in inverse]
    for row in data['free']:assert all(v==0 for v in Uc[row])
    for k,idx in enumerate(data['normal']):
        assert Tc[data['free'][idx]]==[dens[k] if j==k else z for j in range(9)]
    cuts=[i for i in range(12) if i not in data['normal']]
    assert all(not any(Tc[data['free'][i]]) for i in cuts)
    tw=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    torus=[0,1,-1,-1,1];raw=[]
    for idx in model.certificate['raw_free_indices']:
        pair,mon=BASIS[idx]
        raw.append(sum(torus[i] for i in pair)-sum(a*b for a,b in zip(torus,s.Poly(mon,*X).monoms()[0])))
    assert all((raw[tw[i]]+4)%2==0 for i in data['cols'])
    cw=[(raw[tw[i]]+4)//2 for i in data['cols']]+[3]
    assert min(cw)>=0 and max(cw)==4
    beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()}
    _,F=specialize(model,beta,True);H=[];rw=[]
    for i,equation in enumerate(F):
        terms=[];weights={cw[j] for j,v in enumerate(J[i]) if v}
        for mon,c in equation.items():
            if data['fixed'] in mon:continue
            ij=[data['cols'].index(j) for j in mon]+[44]*(2-len(mon));c=pol(c)
            if c:terms.append((ij,c));weights.add(sum(cw[j] for j in ij))
        assert len(weights)==1
        rw.append(next(iter(weights)));H.append(terms)
    assert all(not T[i][j] or cw[i]==NW[j] for i in range(45) for j in range(9))
    assert all(not U[i][j] or cw[i]==rw[sel[j]] for i in range(45) for j in range(33))
    assert max(rw[i] for i in sel)<=4
    inp=wd/'input';rawout=wd/'raw';rowsfile=wd/'rows';meta=wd/'meta';weightsfile=wd/'weights';coeff=wd/'coeff'
    def write_poly(stream,p,normals=0):
        stream.write(str(len(p))+'\n')
        for m,c in p.items():stream.write(str(c)+' '+' '.join(map(str,m+(0,)*normals))+'\n')
    with inp.open('w') as stream:
        stream.write('135 2 9 4 33\n'+' '.join(map(str,sel))+'\n');write_poly(stream,den,9)
        for M in [J,T,U]:
            for row in M:
                for p in row:write_poly(stream,p,9)
        for terms in H:
            stream.write(str(len(terms))+'\n')
            for (i,j),p in terms:stream.write(f'{i} {j}\n');write_poly(stream,p,9)
    weightsfile.write_text(' '.join(map(str,cw+rw+NW))+'\n')
    with coeff.open('w') as stream:
        write_poly(stream,den);stream.write(str(len(known))+'\n')
        for p in known:write_poly(stream,p)
    return inp,rawout,rowsfile,meta,weightsfile,coeff,known


def same_component(cases,linear):
    model=S61(); cache={};result=[]
    variables=s.symbols('e d p q r h k')
    beta={int(name[1:])-1:s.sympify(value) for name,value in linear['beta'].items()}
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    bv=s.Matrix([beta[i] for i in one])
    action=s.Matrix(json.loads((ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json').read_text())['primitive_action']).applyfunc(s.sympify)
    action=action.subs({s.Symbol('u'+str(i+1)):beta[i] for i in one},simultaneous=True)
    source=action.row_join(bv.jacobian(variables))
    for case in cases:
        prime=case['prime']; geometry.PRIME=prime
        if prime not in cache:
            # Q1 only needs the original Jacobian and Hessian. Do not build
            # the large rational T/U inverses used by the older quartic test.
            nuisance=dict(zip(variables[2:],geometry.NUISANCE))
            def encode(expression):
                numerator,denominator=s.cancel(s.sympify(expression).subs(nuisance)).as_numer_denom()
                return [[(mon,geometry.mod(c)) for mon,c in s.Poly(poly,*variables[:2]).terms()]
                        for poly in (numerator,denominator)]
            raw={'T':[],'U':[],'J':[[encode(v) for v in row] for row in linear['J']],'H':[]}
            _,equations=specialize(model,{i:v.subs(nuisance) for i,v in beta.items()},True)
            for i,equation in enumerate(equations):
                for mon,c in equation.items():
                    if linear['fixed'] in mon:continue
                    a,b=[linear['cols'].index(v) for v in mon]+[44]*(2-len(mon))
                    raw['H'].append((i,a,b,encode(c if a==b else c/2)))
            cache[prime]=raw
        core=case['point']; _,_,J,H=geometry.evaluate(cache[prime],core)
        ordinary=geometry.ordinary_same_component(linear,J,H)
        values=dict(zip(variables,core+geometry.NUISANCE))
        B=np.array([[geometry.mod(v) for v in row] for row in source.subs(values).tolist()],dtype=np.int64)
        gradient=case['gradient']
        domain=np.column_stack([B[:,:10],B[:,12:],(gradient[1]*B[:,10]-gradient[0]*B[:,11])%prime])
        dominance=geometry.rank_witness(domain,16)
        cj=np.array([[geometry.mod(s.sympify(v).subs(values)) for v in row]for row in linear['curve_jacobian']],dtype=np.int64)
        vertex=geometry.rank_witness(cj,2)
        family_minor=geometry.mod(s.sympify(linear['family_parameter_minor']).subs(values));assert family_minor
        slice_denominator=geometry.mod(4-33*values[variables[0]]);assert slice_denominator
        result.append(dict(id=case['id'],prime=prime,core=core,nuisance=geometry.NUISANCE,
            ordinary_same_P2_component=ordinary,source_beta=dominance,
            source_beta_matrix=domain.tolist(),vertex=vertex,family_parameter_minor=family_minor,
            seven_slice_denominator=slice_denominator))
        print('geometry',case['id'],'PASS',flush=True)
    return result


def verify():
    started=time.time();definition=json.loads(DATA.read_text());curves=definition['curves'];assert len(curves)==17
    finite=ROOT/'reports/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.json'
    frozen=json.loads(finite.read_text());assert frozen['all_passed']
    native=['code/parametric_implicit_jets_20261101.cpp',
            'code/s61_binary_deeper_weighted_extract_20261102.cpp',
            'code/s61_binary_weighted_seventeen_curve_minors_20261102.cpp']
    logs=[]
    with tempfile.TemporaryDirectory(prefix='s61-seventeen-curves-') as temp:
        wd=Path(temp);inp,rawout,rowsfile,meta,weightsfile,coeff,known=prepare(wd)
        factorsfile=wd/'factors';pointsfile=wd/'points';casefile=wd/'cases.json'
        def write_sparse(stream,sparse):
            stream.write(str(len(sparse))+'\n')
            for m,c in sparse:stream.write(str(c)+' '+' '.join(map(str,m))+'\n')
        with factorsfile.open('w') as stream:
            stream.write('17\n')
            for curve in curves:write_sparse(stream,curve['sparse'])
            stream.write('13\n')
            for p in known:write_sparse(stream,list(p.items()))
        with pointsfile.open('w') as stream:
            stream.write('17\n')
            for curve in curves:
                p=curve['candidate_point'];assert curve['id']==p['curve']
                stream.write(' '.join(map(str,[curve['id'],p['prime'],*p['point']]))+'\n')
        env=dict(os.environ);env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
        binaries=[]
        for k,path in enumerate(native):
            binary=wd/f'engine{k}';binaries.append(binary)
            subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),
                '-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(ROOT/path),
                str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
        for command in [[binaries[0],inp,rawout],[binaries[1],rawout,coeff,rowsfile,meta,weightsfile],
                        [binaries[2],rowsfile,pointsfile,factorsfile,casefile]]:
            run=subprocess.run(list(map(str,command)),check=True,capture_output=True,text=True,env=env);logs.append(run.stderr)
        generated={k:sha(p) for k,p in [('raw_input',inp),('raw_output',rawout),('weighted_rows',rowsfile)]}
        assert generated==frozen['generated_sha256']
        exact=json.loads(casefile.read_text());assert exact['all_passed']
    linear=json.loads((ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json').read_text())['linear']
    geometry_cases=same_component(exact['cases'],linear)
    paths=[Path(__file__),NOTE,DATA,finite,*[ROOT/p for p in native],
        ROOT/'code/s61_binary_deeper_quartic_field_20261101.hpp',
        ROOT/'code/s61_binary_deeper_core_artifacts_20261101.py',
        ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py',
        ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json',
        ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json',
        ROOT/'reports/S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101.json',
        ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
        ROOT/'research/notes/S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101.md',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_BINARY_WEIGHTED_SEVENTEEN_CURVES_20261102',all_passed=True,
        curve_definitions=str(DATA.relative_to(ROOT)),QQ_irreducible_distinct_and_coprime_to_D=True,
        weighted_cases=exact['cases'],same_component_geometry=geometry_cases,
        generated_sha256=generated,native_logs=logs,
        seventeen_explicit_curve_strata_removed_over_sufficiently_general_beta=True,
        proof_note=str(NOTE.relative_to(ROOT)),
        limitations=['No completeness of the seventeen-curve list is asserted.',
            'No all-core weighted-rank cover or affine exhaustion is asserted.',
            'The weighted upper65 witnesses have nuisance zero; the same-C lower65 and dominance witnesses use separate nuisance values.',
            'No ambient nilpotence is asserted at the particular nuisance-zero Hensel lift.'],
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},seconds=time.time()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    assert not args.output.exists()
    report=verify()
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
