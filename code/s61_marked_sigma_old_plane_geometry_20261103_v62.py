#!/usr/bin/env python3
"""Exact orbit/incidence beta-image ranks on the h^2+2 heavy plane."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce
from toric_closed_jet_coordinates_20261102 import coords,pairs,m2
from s61_model import X
from weight1_structure import sp_basis

ROOT=Path(__file__).resolve().parents[1]


def dm(values,domain):
    return DomainMatrix.from_rep(DDM(values,(len(values),len(values[0])),domain))


def certificate(lifter,K,h,u,v,name):
    zero=K.zero
    normals=[zero]*13;normals[9]=h*u;normals[10]=u;normals[12]=v
    point=lifter.lift(h,normals,domain=K)
    assert point['all135_zero'] and not point['t']
    q=point['q'];one,two=lifter.one,lifter.two
    linear=[[zero]*16 for _ in range(135)]
    quadratic=[[zero]*16 for _ in range(16)]
    constant=[zero]*135
    for row,form in enumerate(lifter.original):
        for mon,numerator,denominator in form:
            if sum(i in two for i in mon)!=2:continue
            value=coerce(K,s.Rational(numerator,denominator))
            beta=[]
            for index in mon:
                if index in two:value*=q[two.index(index)]
                else:beta.append(one.index(index))
            if len(beta)==0:constant[row]+=value
            elif len(beta)==1:linear[row][beta[0]]+=value
            elif len(beta)==2:quadratic[beta[0]][beta[1]]+=value
            else:raise AssertionError(beta)
    assert not any(constant)
    M=dm(linear,K);H=M.nullspace().transpose();hdim=H.shape[1]
    Q=dm(quadratic,K)
    isotropic=H.transpose()*(Q+Q.transpose())*H
    full_linear_incidence=isotropic.is_zero_matrix
    Hsym=H.to_Matrix();qexpr=[K.to_sympy(value) for value in q]
    physical=coords()['E']*s.Matrix(qexpr);y=s.Matrix(X[1:]);tensor=s.zeros(4)
    for z,(i,j) in enumerate(pairs):
        tensor[i,j]=sum(physical[z*10+t]*y[a]*y[b] for t,(a,b) in enumerate(m2))
    tensor-=tensor.T
    orbit=[]
    for g in sp_basis():
        derivative=g*tensor+tensor*g.T-sum(((g*y)[i]*tensor.diff(y[i]) for i in range(4)),s.zeros(4))
        coefficients=s.Matrix([s.Poly(derivative[i,j],*y).coeff_monomial(y[a]*y[b])
                               for i,j in pairs for a,b in m2])
        raw=coords()['Pi']*coefficients
        column=[K.from_sympy(value) for value in raw]
        # Projectivize in q25=-1 chart.
        orbit.append([value+q[i]*column[25] for i,value in enumerate(column)])
    O=dm([list(column) for column in zip(*orbit)],K)
    orbit_dimension=O.rank()
    print('GEOMETRY',name,'H',hdim,'ORBIT',orbit_dimension,flush=True)
    action_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    action=s.Matrix(json.loads(action_path.read_text())['primitive_action']).applyfunc(s.sympify)
    symbols=[s.Symbol(f'u{i+1}') for i in one]
    representation=[]
    for j in range(10):
        R=action[:,j].jacobian(symbols)
        assert (R*s.Matrix(symbols)-action[:,j]).applyfunc(s.expand)==s.zeros(16,1)
        representation.append(R)
    z=s.symbols(f'z0:{hdim}');F=K.frac_field(*z)
    Hf=H.convert_to(F)
    beta=Hf*dm([[F.gens[i]] for i in range(hdim)],F)
    beta_list=[row[0] for row in beta.to_list()]
    actions=[]
    for R in representation:
        actions.append([sum((coerce(F,c)*value for c,value in zip(row,beta_list) if c and value),F.zero)
                        for row in R.tolist()])
    A=dm([list(row) for row in zip(*actions)],F)
    annihilator=H.transpose().nullspace().convert_to(F)
    normal=annihilator*A
    R=K.poly_ring(*z)
    normal_poly=normal.convert_to(R)
    _,_,pivot_columns=normal_poly.rref_den(method='FF')
    normal_rank=len(pivot_columns)
    print('BETA_RANK',name,hdim+normal_rank,flush=True)
    if normal_rank:
        for seed in range(1,12):
            values=[K.convert(seed+i) for i in range(hdim)]
            def evaluate(value):
                return sum((c*s.prod(values[i]**power for i,power in enumerate(mon) if power)
                            for mon,c in value.items()),K.zero)
            numeric=dm([[evaluate(value) for value in row] for row in normal_poly.to_list()],K)
            _,witness_columns=numeric.rref()
            if len(witness_columns)==normal_rank:break
        else:raise AssertionError('No full-rank rational witness')
        selected_columns=numeric.extract(range(normal.shape[0]),witness_columns)
        _,selected_rows=selected_columns.transpose().rref()
        minor=normal_poly.extract(selected_rows,witness_columns).det()
        assert minor
        minor_expression=str(lifter.expression(F.convert(minor,R),F,s.Symbol('h')))
    else:
        selected_rows=[];witness_columns=[];minor_expression='1'
    left=normal_poly.transpose().nullspace()
    covectors=left*annihilator.convert_to(R)
    assert (covectors*H.convert_to(R)).is_zero_matrix
    assert (covectors*A.convert_to(R)).is_zero_matrix
    return dict(case=name,u=str(lifter.expression(u,K,s.Symbol('h'))),
                v=str(lifter.expression(v,K,s.Symbol('h'))),
                highest_linear_rank=M.rank(),H_dimension=hdim,
                full135_vanishes_on_entire_linear_H=full_linear_incidence,
                orbit_dimension=orbit_dimension,total_incidence_dimension=hdim+orbit_dimension,
                generic_beta_differential_rank=hdim+normal_rank,
                normal_rank=normal_rank,normal_minor=minor_expression,
                normal_minor_rows=list(selected_rows),normal_minor_columns=list(witness_columns),
                annihilating_covectors=[[str(lifter.expression(F.convert(value,R),F,s.Symbol('h')))
                                         for value in row] for row in covectors.to_list()],
                H=[[str(lifter.expression(value,K,s.Symbol('h'))) for value in row] for row in H.to_list()],
                q=[str(lifter.expression(value,K,s.Symbol('h'))) for value in q])


def verify(model_folder,output):
    start=time.monotonic();lifter=load_model(model_folder)
    K=s.QQ.algebraic_field(s.sqrt(-2),alias='h');h=K.unit
    cases=[]
    for u,v,name in [(K.one,K.zero,'distinct-roots-u'),(K.zero,K.one,'distinct-roots-v'),
                     (h,coerce(K,s.Rational(2,3))*h,'double-root'),(K.zero,K.zero,'origin-flag')]:
        cases.append(certificate(lifter,K,h,u,v,name))
        output.with_suffix('.partial.json').write_text(json.dumps(cases,indent=2)+'\n')
    report=dict(id='S61_MARKED_SIGMA_OLD_PLANE_GEOMETRY_20261103_v62',all_passed=True,
                h_square=-2,cubic='-Z^3/3+h*u*Z*W^2/4+v*W^3/4',
                discriminant_of_12G='-432*(2*h*u^3+9*v^2)',cases=cases,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',
                                        ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json']},
                seconds=time.monotonic()-start)
    with output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-folder',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=verify(args.model_folder,args.output)
    print('PASS',[(c['case'],c['H_dimension'],c['orbit_dimension'],c['generic_beta_differential_rank'])
                  for c in report['cases']],report['seconds'])
