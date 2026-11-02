#!/usr/bin/env python3
"""Compute first-order contact at binary/flag infinity and the intrinsic degree-ten pencil."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as sp

from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import J, COORDS

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_AFFINE_SURFACE_DEGREE_OBSTRUCTION_20261102'


def homogeneous_column(data, coordinate, scale, forced_zero):
    one=[i for i,w in enumerate(data['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(data['parameter_weights']) if w==2]
    b={j:sp.Integer(0) if j+1 in forced_zero else sp.Symbol('u'+str(j+1)) for j in one}
    leading=sp.zeros(70,45)
    time=sp.zeros(70,1)
    for i,equation in enumerate(data['equations'][:70]):
        for mon,num,den in equation:
            assert sum(data['parameter_weights'][j] for j in mon)==4
            quadratic=[j for j in mon if j in two]
            coeff=sp.Rational(num,den)
            if len(quadratic)==2:
                for h,j in enumerate(quadratic):
                    if quadratic[1-h]==coordinate:
                        leading[i,two.index(j)]+=scale*coeff
            if quadratic==[coordinate]:
                time[i]+=scale*coeff*sp.prod(b[j] for j in mon if j in b)
    pc=list(leading.to_DM().rref()[1])
    pr=list(leading[:,pc].T.to_DM().rref()[1])
    residual=time-leading[:,pc]*leading.extract(pr,pc).inv()*time[pr,:]
    return len(pc),list(dict.fromkeys(sp.factor(f) for f in residual if f!=0))


def check_radical(polynomials, linear_generators, variables):
    def represent(polys,target):
        rows=[sp.Poly(f,*variables) for f in list(polys)+[target]]
        monomials=sorted(set().union(*(set(row.monoms()) for row in rows)))
        matrix=sp.Matrix([[f.coeff_monomial(m) for m in monomials]for f in rows])
        reduced=matrix[:-1,:].to_DM().rref()[0].to_Matrix()
        assert matrix.rank()==reduced.rank()
    substitutions={}
    for linear in linear_generators:
        value=linear.subs(substitutions)
        variable=next(x for x in variables if value.diff(x)!=0)
        coefficient=value.diff(variable)
        assert coefficient.is_number and coefficient!=0
        substitutions[variable]=sp.expand(variable-value/coefficient)
    assert all(sp.expand(f.subs(substitutions))==0 for f in polynomials)
    # Squares in the previous zero locus give the converse, with no generic division.
    previous={}
    for linear in linear_generators:
        represent([f.subs(previous) for f in polynomials],linear.subs(previous)**2)
        value=linear.subs(previous)
        variable=next(x for x in variables if value.diff(x)!=0)
        previous[variable]=sp.expand(variable-value/value.diff(variable))


def verify():
    model=S61()
    y,u,beta=primitive_beta(model)
    by_number={int(str(x)[1:]):x for x in u}
    us=lambda i:by_number[i]
    binary_rank,binary=homogeneous_column(model.reduced,29,1,{12,13})
    flag_rank,flag=homogeneous_column(model.reduced,37,-1,{4})
    assert (binary_rank,flag_rank)==(16,15)
    K,L=5*us(35)+us(50),5*us(14)-us(49)
    check_radical(binary,[K,L],u)
    check_radical(flag,[us(5),5*us(14)-4*us(49)],u)
    # Transform the marked factor x from the flag plane <e1,e3> to the binary plane <e1,e2>.
    y=sp.Matrix(y)
    A=sp.eye(4)
    A[2,2]=A[3,3]=0
    A[2,3],A[3,2]=-1,1
    assert A.T*J*A==J
    tensor=lambda M:sp.Matrix([M[i,j].diff(y[k]) for (i,j),k in COORDS])
    basis=tensor(beta).jacobian(u)
    rows=list(basis.T.to_DM().rref()[1])
    transformed=A.T*beta.subs(dict(zip(y,A*y,strict=True)),simultaneous=True)*A
    new=basis[rows,:].inv()*tensor(transformed)[rows,:]
    marked=L.subs(dict(zip(u,new,strict=True)),simultaneous=True)
    assert sp.expand(marked+5*us(14)-4*us(49))==0
    # Reconstruct the vertex quadratic on e1-perp/e1 from all original primitive coefficients.
    v=sp.eye(4)[:,1]
    z,w=sp.symbols('z w')
    vector=sp.Matrix([0,0,z,w])
    qv=sp.expand((vector.T*beta.subs(dict(zip(y,vector,strict=True)),simultaneous=True)*v)[0])
    assert sp.expand(qv-us(13)*z*z-(2*us(14)-us(49))*z*w+us(5)*w*w)==0
    q=sp.Poly(qv,z,w)
    coeff=[q.coeff_monomial(m) for m in [z*z,z*w,w*w]]
    eta=us(49)
    six=sp.Matrix([us(4),us(12),*coeff,eta]).jacobian(u)
    assert six.rank()==6
    disc=sp.expand(coeff[1]**2-4*coeff[0]*coeff[2])
    r=(us(5)*us(13)-us(14)*eta+us(14)**2)/eta**2
    assert sp.cancel(disc/eta**2-1-4*r)==0
    c1,c2,c3,extra,root=sp.symbols('c1 c2 c3 eta root')
    hessian=sp.hessian(c2*c2-4*c1*c3-root*extra**2,[c1,c2,c3,extra])
    assert hessian[:3,:3].det()!=0
    sources=[
        'S61_GENERIC_BOUNDARY_INDEX_20261102','S61_AFFINE_SURFACE_BOUND_20261102',
        'S61_AFFINE_BOUNDARY_DIMENSION_20261102','S61_BINARY_DEEPER_REDUCTION_20261101',
        'S61_NONBINARY_FORMAL_RING_20261101','S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101',
        'S61_R0_MIXED_CURVE_EXHAUSTION_20261101','S61_R0_T5_ZERO_L_EXHAUSTION_20261101',
        'S61_R0_MIXED_SURFACE_T7_EXHAUSTION_20261101','S61_MINIMAL_VERTEX_EXHAUSTION_20261101',
        'S61_MINIMAL_MAIN_REMOVAL_20261101','S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101',
        'S61_MINIMAL_RANKONE_A1_REMOVAL_20261101','S61_MINIMAL_M0_BASELINE_REMOVAL_20261101',
        'S61_MINIMAL_M1_BASELINE_REMOVAL_20261101','S61_MINIMAL_MIXED_TAILS_REMOVAL_20261101',
        'S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101','S61_MINIMAL_PURE_CUBIC_TAILS_REMOVAL_20261101',
        'S61_MINIMAL_SEMISIMPLE_TAILS_REMOVAL_20261101','S61_R0_GENERIC_EXHAUSTION_20261101',
        'S61_FLAG_ORDINARY_CORE_SUPPORT_20261102','S61_FLAG_ZERO12_RELATIVE_WEIGHT_REMOVAL_20261102',
        'NULL_CORRELATION_DUALITY_20261101']
    files=[Path(__file__),ROOT/'research/notes'/(STEM+'.md')]
    for name in sources:
        path=ROOT/'reports'/(name+'.json')
        assert json.loads(path.read_text())['all_passed']
        files.extend([path,ROOT/'research/notes'/(name+'.md')])
    files.extend(ROOT/'code'/n for n in ['s61_model.py','s61_infinity_flags.py','weight1_structure.py',
                                        'primitive_jet_elliptic_curve_20261101.py','elliptic_weight1_dominance.py'])
    files.extend(ROOT/'legacy/v02/data'/n for n in ['symplectic_slice_certificate.json',
                'symplectic_slice_equations.json','symplectic_slice_reduced_equations.json'])
    return dict(id=STEM,all_passed=True,field='QQ identities, generic base K=C(W16)',
                q_linearization_ranks=[binary_rank,flag_rank],binary_residuals=list(map(str,binary)),
                flag_residuals=list(map(str,flag)),binary_radical=[str(K),str(L)],
                flag_radical=[str(us(5)),str(5*us(14)-4*us(49))],marked_flag_factor=str(marked),
                simple_contact_boundary_curve_only_possible='marked square section Sigma',
                marked_square_section_degree=15,positive_dimensional_surface_degree_lower=10,
                surface_degree_divisor=5,degree_ten_surface_boundary='twice a degree-five closed curve consisting geometrically of five conjugate vertical binary lines over C_beta',
                vertex_quadratic=str(qv),vertex_discriminant=str(disc),vertex_six_evaluation_rank=6,
                finite_core_vertex_closed_point_degree=10,ordinary_pencil_ratio='disc(q_v)/eta^2=1+4r',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print('PASS',report['q_linearization_ranks'],len(report['input_sha256']))
