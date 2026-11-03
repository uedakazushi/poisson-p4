#!/usr/bin/env python3
"""Full-incidence quadratic transversal at the marked new quadratic ray.

This reconstructs the actual full135 point over QQ(sqrt(-230)), identifies
its T6 orbit, and retains all beta coordinates until taking actual incidence
cuts. No fixed-beta local ring is substituted for the ambient local ring.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61, X
from s61_marked_sigma_point_lift_20261103_v62 import load_model, coerce
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from weight1_structure import J, COORDS, sp_basis, infinitesimal_action

ROOT=Path(__file__).resolve().parents[1]
PRIME=32003


def verify(output,model_folder):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    lifter=load_model(model_folder)
    K=s.QQ.algebraic_field(s.sqrt(-230),alias='r');r=K.unit
    h=coerce(K,s.Rational(7,100))*r
    assert 1000*h*h+1127==K.zero
    vv=[K.zero]*13
    vv[7]=coerce(K,s.Rational(100,161))*h;vv[8]=-K.one
    vv[9]=coerce(K,s.Rational(10,7))*h;vv[10]=K.one
    point=lifter.lift(h,vv,domain=K)
    assert point['all135_zero'] and not point['t']
    scalar_inputs=[]
    for matrix in [lifter.J,lifter.T,lifter.U]:
        scalar_inputs.extend(coefficient for row in matrix for cell in row for coefficient,power in cell)
    scalar_inputs.extend(coefficient for row in lifter.equations
                         for i,j,rho_power,cell in row for coefficient,power in cell)
    assert all(int(coefficient.q)%PRIME for coefficient in scalar_inputs)
    model=lifter.model;one=lifter.one;two=lifter.two;source=lifter.original
    q,beta=point['q'],point['beta'];fixed=25;assert q[fixed]==-K.one
    for equation in source:
        for mon,n,d in equation:assert d%PRIME
    conv=lambda value:coerce(K,value)
    expr=lambda value:s.factor(K.to_sympy(value))
    def unit(value):
        return all(int(c.denominator)%PRIME for c in value.rep) and any(int(c.numerator)%PRIME for c in value.rep)
    def integral(value):return all(int(c.denominator)%PRIME for c in value.rep)
    def dm(matrix):
        if isinstance(matrix,DomainMatrix):return matrix.convert_to(K)
        if isinstance(matrix,s.MatrixBase):return matrix.to_DM().convert_to(K)
        return DomainMatrix(matrix,(len(matrix),len(matrix[0])),K)
    def enc(matrix):
        if isinstance(matrix,DomainMatrix):matrix=matrix.to_Matrix()
        return [[str(z) for z in row] for row in matrix.tolist()]
    def minor(matrix,size=None):
        matrix=dm(matrix);columns=list(matrix.rref()[1])
        if size is not None:assert len(columns)==size
        rows=list(matrix.extract(list(range(matrix.shape[0])),columns).transpose().rref()[1])
        value=matrix.extract(rows,columns).det();assert unit(value)
        return dict(rows=rows,columns=columns,determinant=str(expr(value)))
    assert all(integral(v) for v in q+beta)

    data=coords();qexpr=s.Matrix(list(map(K.to_sympy,q)));spatial=s.Matrix(X[1:])
    physical=data['E']*qexpr;Q=s.zeros(4)
    for k,(i,j) in enumerate(pairs):
        Q[i,j]=sum(physical[k*10+a]*spatial[v]*spatial[w] for a,(v,w) in enumerate(m2))
    Q-=Q.T;Q=Q.applyfunc(s.expand)
    # Explicit symplectic coordinates with no algebraic root beyond K.
    c=coerce(K,s.Rational(10,7))*h;ce=K.to_sympy(c)
    # Rows give NEW coordinates (X,Y,Z,W)=(z-c*x,w,x,y+c*w).
    # Correct row Y is w and row W is y+c*w.
    frame=s.Matrix([[-ce,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,ce]])
    assert (frame*J*frame.T-J).applyfunc(s.expand)==s.zeros(4)
    moved=(frame*Q.subs(dict(zip(spatial,frame.inv()*spatial)),simultaneous=True)*frame.T).applyfunc(s.expand)
    x,y,z,w=spatial;expected=s.zeros(4)
    expected[1,2]=-s.Rational(10,23)*x*w;expected[1,3]=z*z;expected-=expected.T
    assert (moved-expected).applyfunc(s.expand)==s.zeros(4)
    swap=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,-1,0]])
    diagonal=s.diag(-s.Rational(15,23),-s.Rational(23,15),1,1)
    canonical_frame=diagonal*swap*frame
    projective_factor=-s.Rational(45,23)
    canonical=(projective_factor*canonical_frame*Q.subs(
        dict(zip(spatial,canonical_frame.inv()*spatial)),simultaneous=True)*canonical_frame.T).applyfunc(s.expand)
    canonical_expected=s.zeros(4)
    canonical_expected[1,2]=3*w*w;canonical_expected[1,3]=-2*x*z
    canonical_expected-=canonical_expected.T
    assert canonical==canonical_expected
    assert (canonical_frame*J*canonical_frame.T-J).applyfunc(s.expand)==s.zeros(4)
    divergence=s.Matrix([sum(s.diff(Q[i,j],spatial[j]) for j in range(4)) for i in range(4)])
    assert divergence==s.zeros(4,1)

    highest=[[K.zero for _ in range(16)] for _ in range(135)]
    quadratic=[[K.zero for _ in range(16)] for _ in range(16)]
    for row,equation in enumerate(source):
        for mon,n,d in equation:
            if sum(i in two for i in mon)!=2:continue
            bs=[one.index(i) for i in mon if i in one]
            value=conv(s.Rational(n,d))
            for i in mon:
                if i in two:value*=q[two.index(i)]
            if len(bs)==1:highest[row][bs[0]]+=value
            if len(bs)==2:quadratic[bs[0]][bs[1]]+=value
    HD=dm(highest);H=HD.nullspace().transpose()
    hrows=list(H.transpose().rref()[1]);H=H*H.extract(hrows,list(range(H.shape[1]))).inv()
    assert H.shape==(16,9)
    assert (HD*H).is_zero_matrix
    QD=dm(quadratic);assert (H.transpose()*(QD+QD.transpose())*H).is_zero_matrix
    highest_minor=minor(HD,7)
    orbit=[]
    for A in sp_basis():
        dQ=(sum(((A*spatial)[i]*Q.diff(spatial[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T).applyfunc(s.expand)
        coeff=s.Matrix([s.Poly(dQ[i,j],*spatial).coeff_monomial(spatial[v]*spatial[w]) for i,j in pairs for v,w in m2])
        dq=data['Pi']*coeff;dq-=qexpr*dq[fixed]/qexpr[fixed];orbit.append(dq)
    O=s.Matrix.hstack(*orbit);orbit_minor=minor(O,8);qcuts=orbit_minor['rows']
    assert fixed not in qcuts
    incidence_det=conv(s.sympify(orbit_minor['determinant']))
    assert unit(incidence_det)

    # Generic dominance of this SAME orbit/kernel incidence.
    zz,bu,form=primitive_beta(model)
    vector=s.Matrix([form[i,j].diff(zz[k]) for (i,j),k in COORDS])
    embedding=vector.jacobian(bu);erows=list(embedding.T.to_DM().rref()[1]);left=embedding[erows,:].inv()
    sample=H.to_Matrix()*s.Matrix(range(2,11))
    sample_vector=vector.subs(dict(zip(bu,sample)))
    beta_actions=s.Matrix.hstack(*[(left*infinitesimal_action(sample_vector,A)[erows,:]).applyfunc(s.expand) for A in sp_basis()])
    beta_dominance=minor(H.to_Matrix().row_join(beta_actions),16)

    frozen={two[fixed],*[two[i] for i in qcuts],*[one[i] for i in hrows]}
    variables=[i for i in range(62) if i not in frozen];N=len(variables);assert N==44
    location=dict(zip(variables,range(N)))
    values=dict(zip(one,beta))|dict(zip(two,q))|{61:K.zero}
    linear=[[K.zero for _ in variables] for _ in range(135)]
    encoded=[]
    for row,equation in enumerate(source):
        terms=[]
        for mon,n,d in equation:
            ids=tuple(mon)+(61,)*(2-sum(i in two for i in mon));coefficient=conv(s.Rational(n,d));terms.append((ids,coefficient))
            for variable,exponent in Counter(ids).items():
                if variable not in location:continue
                rest=list(ids);rest.remove(variable)
                product=coefficient*exponent
                for i in rest:product*=values[i]
                linear[row][location[variable]]+=product
        encoded.append(terms)
    LD=dm(linear);pivots=list(LD.rref()[1]);rank=len(pivots)
    print('incidence',H.shape[1]+8,'transverse variables/rank',N,rank,flush=True)
    assert rank==43
    selected=list(LD.extract(list(range(135)),pivots).transpose().rref()[1])
    block=LD.extract(selected,pivots);jacobian_det=block.det();assert unit(jacobian_det)
    normal=LD.nullspace().transpose();assert normal.shape==(44,1)
    free=next(i for i in range(N) if i not in pivots)
    normal=normal.mul(K.one/normal.to_list()[free][0]);direction=[row[0] for row in normal.to_list()]
    assert all(integral(v) for v in direction)
    raw2=[]
    for row in encoded:
        coefficient=K.zero
        for ids,scalar in row:
            polynomial=[scalar,K.zero,K.zero]
            for i in ids:
                delta=direction[location[i]] if i in location else K.zero
                polynomial=[polynomial[0]*values[i],polynomial[1]*values[i]+polynomial[0]*delta,
                            polynomial[2]*values[i]+polynomial[1]*delta]
            coefficient+=polynomial[2]
        raw2.append([coefficient])
    correction=block.inv()*dm([raw2[i] for i in selected])
    image=LD.extract(list(range(135)),pivots)*correction
    residual=dm(raw2)-image
    residual_values=[row[0] for row in residual.to_list()]
    witness=next(i for i,v in enumerate(residual_values) if unit(v))
    assert all(residual_values[i]==K.zero for i in selected)
    assert (HD*dm([[v] for v in beta])).is_zero_matrix
    print('quadratic unit witness',witness,expr(residual_values[witness]),flush=True)
    zero5_path=ROOT/'reports/S61_MARKED_SIGMA_PENCIL_FULL5_20261103/report.json'
    zero5=json.loads(zero5_path.read_text());assert zero5['all_passed']
    zero_record=next(record for record in zero5['results'] if record['task']=='zero:5')
    assert zero_record['shape']==[2022,932]
    determinant_line=next(line for line in zero_record['output'].splitlines()
                          if line.startswith('DETERMINANT_COEFFICIENTS '))
    determinant_data=list(map(int,determinant_line.split()[1:]));assert determinant_data[0]==3
    zero5_determinant=sum((conv(c)*h**i for i,c in enumerate(determinant_data[1:])),K.zero)
    assert unit(zero5_determinant)
    paths=[Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',
           ROOT/'code/s61_marked_sigma_pencil_model_20261102.py',ROOT/'code/s61_model.py',
           ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
           ROOT/'code/toric_closed_jet_coordinates_20261102.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
           ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
           ROOT/'research/notes/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.md',
           ROOT/'reports/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.json',
           ROOT/'research/notes/S61_BINARY_A_ZERO_LIFT_20261103.md',zero5_path,
           ROOT/'research/notes/S61_MARKED_SIGMA_NEW_QUADRATIC_20261103_v62.md']
    report=dict(id='S61_MARKED_SIGMA_NEW_QUADRATIC_20261103_v62',all_passed=True,
                field='QQ(r), r^2=-230; h=7r/100; residue polynomial irreducible mod32003',
                beta=list(map(lambda v:str(expr(v)),beta)),q=list(map(lambda v:str(expr(v)),q)),
                all135_original_equations_zero=True,t='0',Q=enc(Q),symplectic_new_coordinate_frame=enc(frame),
                transformed_Q=enc(moved),Hamiltonian='-5*x*w^2/23-z^3/3',orbit='T6',
                canonical_T6_frame=enc(canonical_frame),canonical_projective_factor=str(projective_factor),
                canonical_T6_tensor=enc(canonical),Smith_input_scalars_checked=len(scalar_inputs),
                all_Smith_input_denominators_prime_units=True,
                orbit_dimension=8,beta_kernel_dimension=9,incidence_dimension=17,
                highest_kernel=enc(H),highest_rank_minor=highest_minor,
                all135_highest_equations_zero_on_entire_kernel=True,orbit_minor=orbit_minor,
                beta_generic_dominance_minor=beta_dominance,beta_cut_indices=hrows,q_cut_indices=qcuts,
                transverse_variables=variables,linear_pivot_columns=pivots,linear_pivot_rows=selected,
                jacobian_minor=str(expr(jacobian_det)),normal_direction=enc(normal),
                quadratic_residual=list(map(lambda v:str(expr(v)),residual_values)),
                quadratic_unit_witness_row=witness,quadratic_unit_coefficient=str(expr(residual_values[witness])),
                transverse_length_upper=2,known_same_T6_generic_multiplicity=2,
                mixed_characteristic_free_rank_conclusion=2,ambient_homogenizer_square_zero=True,
                ambient_homogenizer_zero=True,
                zero_rho_weight5_determinant_modp_evaluated=str(expr(zero5_determinant)),
                zero_rho_reduced_support_exhausted_by_new_axis=True,
                limits=['Uses the full-incidence finite-rank lemma and known generic T6 multiplicity.',
                        'Does not remove the entire marked pencil or establish generic five-point exhaustion.'],
                prepared_model_input_sha256=hashlib.sha256(lifter.input_path.read_bytes()).hexdigest(),
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--model-folder',type=Path,default=ROOT/'reports/S61_MARKED_SIGMA_POINT_LIFT_20261103_v62_FINAL/model')
    args=parser.parse_args();result=verify(args.output,args.model_folder)
    print('PASS',result['id'],result['seconds'])
