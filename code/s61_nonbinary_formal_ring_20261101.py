#!/usr/bin/env python3
"""Constant generic leading ranks that turn the local bound into a formal ring.

Combined with the smooth fixed-beta curve and local tangent-cone certificate,
these ranks force the homogenizer to vanish and determine square-zero rings.
No assertion is made about all other infinity points or affine-fiber finiteness.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize,EXTRA
from weight1_structure import J,minor_witness
ROOT=Path(__file__).resolve().parents[1]


def calculate(model,extra):
    u=s.symbols('u1:62')
    beta={i:u[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    beta.update({3:0,11:0,4:0,12:0,49:-5*u[34],50:5*u[25]})
    two,polynomials=specialize(model,beta,extra);fixed=two.index(40)
    derivative=s.zeros(len(polynomials),45)
    for row,f in enumerate(polynomials):
        assert f.get((fixed,fixed),0)==0
        for j in range(45):
            derivative[row,j]=-f.get(tuple(sorted((fixed,j))),0)*(2 if j==fixed else 1)
    constant=derivative[:70,:];kernel=s.Matrix.hstack(*constant.nullspace())
    assert constant.rank()==24 and constant*kernel==s.zeros(70,21)
    residual=derivative[70:,:]*kernel
    reduced,pivots=residual.to_DM().rref();reduced=reduced.to_Matrix()
    free=[j for j in range(21) if j not in pivots]
    relations=s.zeros(21,len(free))
    for k,j in enumerate(free):
        relations[j,k]=1
        for row,col in enumerate(pivots):relations[col,k]=-reduced[row,j]
    assert all(s.cancel(f)==0 for f in residual*relations)
    assert relations[free,:]==s.eye(len(free))
    variables=sorted(residual.free_symbols,key=str)
    sample={x:s.Integer((-1)**i*(i*i+3*i+5)) for i,x in enumerate(variables)}
    witness=minor_witness(residual.subs(sample))
    expected=18 if extra else 17
    assert witness['rank']==len(pivots)==expected
    def rows(matrix):return [[str(x) for x in row] for row in matrix.tolist()]
    return dict(equation_count=len(polynomials),canonical_Q='d_y wedge (x*z*d_z-x*w*d_w)',
        beta_linear_subspace={str(u[i]):str(value) for i,value in beta.items()},
        constant_jacobian_rank=24,constant_kernel=rows(kernel),residual_matrix=rows(residual),
        independent_rational_kernel=rows(relations),residual_sample={str(k):str(v) for k,v in sample.items()},
        residual_minor=witness,generic_full_leading_rank=24+expected,
        rank_upper_bound_holds_everywhere='All larger minors are polynomial and vanish on the dense fraction-field open.',
        fixed_beta_projective_curve_generic_normal_embedding_dimension=3-int(extra))


def verify():
    model=S61();primary=calculate(model,False);augmented=calculate(model,True)
    a,b=s.symbols('a b');x,y,z,w=s.symbols('x y z w');v=s.Matrix([x,y,z,w])
    change=s.Matrix([[1,0,0,0],[0,1,a,-b],[b,0,1,0],[a,0,0,1]])
    assert change.T*J*change==J and change.det()==1
    vector=s.Matrix([0,0,b*x*x+x*z,-a*x*x-x*w]);factor=s.Matrix([0,1,0,0])
    old=factor*vector.T-vector*factor.T
    inverse=change.inv()*v
    pushed=(change*old*change.T).subs(dict(zip(v,inverse)),simultaneous=True).applyfunc(s.expand)
    target=factor*s.Matrix([0,0,x*z,-x*w]).T-s.Matrix([0,0,x*z,-x*w])*factor.T
    assert pushed==target
    local_path=ROOT/'reports/S61_NONBINARY_LOCAL_DIMENSION_20261101.json'
    local=json.loads(local_path.read_text());assert local['all_passed'] and local['local_dimension']==1
    assert local['primary134']['common_hessian_rank']==2
    assert local['augmented135']['common_hessian_rank']==1
    inputs=[Path(__file__),ROOT/'code/s61_flag_transverse_obstruction_20261101.py',EXTRA,local_path,
        ROOT/'reports/S61_R0_INFINITY_20261101.json',ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md']
    return dict(id='S61_NONBINARY_FORMAL_RING_20261101',all_passed=True,field='QQ',
        primary134=primary,augmented135=augmented,
        symplectic_translation=[[str(t) for t in row] for row in change.tolist()],
        G_translation='G=x*(z+b*x)*(w+a*x) modulo x^3; the displayed change sends Q to the canonical Q.',
        completed_local_ring_134='C[[u,e1,e2]]/(e1^2,e1*e2,e2^2)',
        completed_local_ring_135='C[[u,e]]/(e^2)',homogenizer_zero_in_both_completions=True,
        proof_note='research/notes/S61_NONBINARY_FORMAL_RING_20261101.md',
        proof_method='Tangent-cone upper length plus generic nonreduced length lower bound gives a finite free C[[u]]-module. The ideal(t) has rank0 in a torsion-free module and vanishes. The curve ideal is free of rank2 or1 and has square zero at the generic point, hence everywhere.',
        limitations=['At the displayed nonbinary sample and the open family where the same local rank conditions hold.','Not the full infinity scheme, all R0 branches, or full S61 generic-fiber exhaustion.'],
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs})


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=verify();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,leading_ranks=[41,42],homogenizer_zero=True)))
if __name__=='__main__':main()
