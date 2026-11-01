#!/usr/bin/env python3
"""Reconstruct the q(5,2) lift at the saved toric first jet by linear algebra.

The 150 linear equations use the known-component constant kernel condition.
They are not asserted to hold on arbitrary S61 components.
"""
import argparse,hashlib,json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import S61,X
from elliptic_weight1_dominance import projective_affine
from weight1_structure import J,minor_witness
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'reports/S61_KNOWN_FIBER_AUDIT_20261031.json'
PAIRS=list(combinations(range(4),2))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    model=S61();record=next(p for p in json.loads(SOURCE.read_text())['points'] if p['component_id']=='T01')
    point=list(map(s.Rational,record['S61_coordinates_u1_through_u61']))
    weights=model.reduced['parameter_weights'];two=[i for i,w in enumerate(weights) if w==2]
    q=s.symbols('q0:45');values=list(point)
    for i,v in zip(two,q,strict=True):values[i]=v
    y=s.symbols('y0:4');p=projective_affine(model.reconstruct(values),X,y);zero=dict.fromkeys(y,0)
    linear=sum((p.diff(yi).subs(zero)*yi for yi in y),s.zeros(4));beta=(-J*linear*J).applyfunc(s.expand)
    # This beta is the full closed first jet in the saved gauge. The graph
    # reconstruction works in dimension20 before the primitive gauge section.
    n=s.Matrix([[0,p[2,3],-p[1,3],p[1,2]],[-p[2,3],0,p[0,3],-p[0,2]],
                [p[1,3],-p[0,3],0,p[0,1]],[-p[1,2],p[0,2],-p[0,1],0]])
    r=s.symbols('r0:5');z=s.Matrix([[0,r[0],r[1],r[2]],[-r[0],0,r[3],r[4]],
                                  [-r[1],-r[3],0,-r[0]],[-r[2],-r[4],r[0],0]])
    ell=s.Matrix([sum(s.diff(beta[i,j],y[k])*z[i,j] for i,j in PAIRS) for k in range(4)])
    w=-J*ell;yy=s.Matrix(y);kernel=z+yy*w.T-w*yy.T
    residual=s.expand(sum(n[i,j]*kernel[i,j] for i,j in PAIRS))
    equations=s.Poly(residual,*y,*r).coeffs();a,b=s.linear_eq_to_matrix(equations,q)
    assert a.shape==(150,45)
    reduced,pivots=a.row_join(b).to_DM().rref();reduced=reduced.to_Matrix()
    assert len(pivots)==45 and all(i<45 for i in pivots)
    solution=reduced[:45,-1];assert a*solution==b
    lift=list(point)
    for i,v in zip(two,solution,strict=True):lift[i]=v
    assert lift!=point
    assert all(model.evaluate(eq,lift)==0 for eq in model.reduced['equations'])
    # Full fixed-jet tangent rank, independent of the extra kernel equations.
    index={old:new for new,old in enumerate(two)};jacobian=s.zeros(134,45)
    for row,equation in enumerate(model.reduced['equations']):
        for mon,num,den in equation:
            for position,old in enumerate(mon):
                if old in index:
                    jacobian[row,index[old]]+=s.Rational(num,den)*s.prod(lift[v] for j,v in enumerate(mon) if j!=position)
    witness=minor_witness(jacobian);assert witness['rank']==45
    gradients=ROOT/'reports/PRIMITIVE_CURVE_INVARIANT_GRADIENTS_20261101.json'
    invariant_checks={}
    if gradients.exists():
        data=json.loads(gradients.read_text());vector=solution.col_join(s.ones(1,1))
        for name in ('c4','c6','j'):
            value=s.Matrix(data['candidate_'+name+'_transport_affine_rows'])*vector
            assert value==s.zeros(4,1);invariant_checks[name]=True
    report={'id':'Q52_FIXED_JET_LIFT_20261101','all_passed':True,'field':'QQ',
        'fixed_first_jet':'T01 log parameters(2,3,5,7,11,13)',
        'S61_coordinates_u1_through_u61':list(map(str,lift)),'q_to_u_indices_zero_based':two,
        'linear_kernel_condition_matrix':[[str(v) for v in row] for row in a.tolist()],
        'linear_kernel_condition_rhs':list(map(str,b)),'linear_rank':45,
        'all134_equations_checked':True,'full_fixed_jet_jacobian_witness':witness,
        'distinct_from_saved_toric_point':True,'known_invariant_transport_checks':invariant_checks,
        'scope':'One additional explicit isolated reduced lift in the same fixed fiber. Kernel equations characterize the known q52 locus, not the full S61 fiber.',
        'source_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in (Path(__file__),SOURCE,ROOT/'code/s61_model.py',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'kernel_linear_rank':45,'full_fiber_jacobian_rank':45,'output':str(args.output)}))

if __name__=='__main__':main()
