#!/usr/bin/env python3
"""Exact projective Lie actions and global magnetic quotients for stable2.

Only constant changes, full six-bracket identities, and exact matrix algebra
are used; no component recognition is inferred from a divisor coincidence.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from classification_nil22_image_source_review import build_P,coords,mag,x,y,q1,q2

ROOT=Path(__file__).resolve().parents[1]


def encode(P):
    return {f'{i},{j}':str(s.factor(P[i,j]))for i,j in itertools.combinations(range(P.rows),2)if P[i,j]!=0}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    specs=[('stable_line_conic',s.Matrix([[3,-15],[-15,3]]),s.Matrix([s.Rational(3,2),-s.Rational(3,2)]),(4,5)),
        ('stable_three_lines',s.Matrix([[-6,3],[3,-6]]),s.zeros(2,1),(3,4,5))]
    ell=s.symbols('l0:6');L=s.Matrix([ell[0]+ell[1]*x+ell[2]*y,ell[3]+ell[4]*x+ell[5]*y]);mons=[1,x,y,x*x,x*y,y*y];records=[]
    for name,S,V0,residual_indices in specs:
        P,M,D=build_P(x*y*(x-y),S,V0,-s.eye(2)/3)
        images=s.Matrix([q1+L[0],q2+L[1],x,y]);J=images.jacobian(s.Matrix(coords))
        shifted=(J*P*J.T).subs({q1:q1-L[0],q2:q2-L[1]},simultaneous=True).applyfunc(s.expand)
        difference=s.Poly(shifted[0,1]-P[0,1],x,y);action=s.Matrix([difference.coeff_monomial(m)for m in mons]).jacobian(ell)
        coker=action.T.nullspace();assert len(coker)==len(residual_indices)
        Q=s.Matrix.hstack(*(s.eye(6)[:,i]for i in residual_indices));C=s.Matrix.vstack(*(v.T for v in coker))
        assert C*Q==s.eye(len(coker));invariants=C*s.Matrix(mag);target_coefficients=Q*invariants
        solutions=list(s.linsolve((action,target_coefficients-s.Matrix(mag)),ell));assert len(solutions)==1
        solution=solutions[0];free=set().union(*(value.free_symbols for value in solution))&set(ell)
        solution=s.Matrix([v.subs(dict.fromkeys(free,0))for v in solution]);translation=L.subs(dict(zip(ell,solution)),simultaneous=True)
        actual=shifted.subs(dict(zip(ell,solution)),simultaneous=True).applyfunc(s.expand)
        expected=P.subs(dict(zip(mag,target_coefficients)),simultaneous=True)
        assert (actual-expected).applyfunc(s.expand)==s.zeros(4)
        A=[]
        for i in range(2):
            # A projective field at h=1 is (A*p)_spatial-p*(A*p)_h.
            matrix=s.zeros(3);matrix[0,i+1]=s.Rational(1,3)
            for j in range(2):
                linear=s.expand(M[i,j]+coords[i+2]*coords[j+2]/3)
                polynomial=s.Poly(linear,x,y);assert polynomial.total_degree()<=1
                matrix[j+1,0]=polynomial.coeff_monomial(1);matrix[j+1,1]=polynomial.coeff_monomial(x);matrix[j+1,2]=polynomial.coeff_monomial(y)
            matrix-=s.trace(matrix)*s.eye(3)/3
            v=matrix*s.Matrix([1,x,y]);assert all(s.expand(v[j+1]-coords[j+2]*v[0]-M[i,j])==0 for j in range(2))
            A.append(matrix)
        commutator=A[0]*A[1]-A[1]*A[0]
        c1,c2=s.symbols('c1 c2');relation=list(s.linsolve(list(commutator-c1*A[0]-c2*A[1]),(c1,c2)))[0]
        assert all(not v.free_symbols for v in relation)
        lie=dict(traceless_projective_matrices=[v.tolist()for v in A],usual_matrix_commutator=commutator.tolist(),commutator_coefficients=list(map(str,relation)))
        if name=='stable_three_lines':
            R=s.Matrix([[3,1,0],[-3,0,1],[-3,1,-1]])
            diagonal=[(R*v*R.inv()).applyfunc(s.factor)for v in A]
            assert all(v.is_diagonal()for v in diagonal)
            assert s.Matrix([[v[1,1]-v[0,0],v[2,2]-v[0,0]]for v in diagonal]).det()!=0
            lie.update(type='Cartan: two commuting diagonalizable independent projective fields',simultaneous_base_coordinates=['3h+x','-3h+y','-3h+x-y'],base_coordinate_matrix=R.tolist(),simultaneous_diagonal_matrices=[v.tolist()for v in diagonal])
        else:
            assert commutator!=s.zeros(3) and commutator**3==s.zeros(3) and (commutator**2).rank()==1
            semisimple=A[0]/relation[1]
            assert (semisimple*commutator-commutator*semisimple-commutator).applyfunc(s.expand)==s.zeros(3)
            assert sorted(semisimple.eigenvals())==[-1,0,1]
            lie.update(type='Two-dimensional affine Lie algebra acting as the restriction of the three-dimensional irreducible sl2 representation to a Borel subalgebra',regular_nilpotent_commutator_index=3,semisimple_matrix=semisimple.tolist(),semisimple_eigenvalues=[-1,0,1],bracket_D_K_equals_K=True)
        records.append(dict(name=name,all_original_magnetic_parameters=list(map(str,mag)),magnetic_action_rank=action.rank(),invariants=list(map(str,invariants)),
            exact_parameter_polynomial_translation=list(map(str,translation)),no_parameter_localization=True,
            canonical_magnetic_polynomial=str(s.expand(sum(a*b for a,b in zip(target_coefficients,mons)))),all_six_bracket_identity=True,
            entire_residual_source_brackets=encode(expected),projective_Lie_action=lie,
            scope='An exact whole-source normal form and projective Lie-action reduction; no ambient component or MPS membership is asserted.'))
    paths=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json',ROOT/'reports/CLASSIFICATION_NIL22_STABLE_INFINITY_20261107.json']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,cases=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6))
    with args.output.open('x')as handle:json.dump(result,handle,indent=2,default=str);handle.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))


if __name__=='__main__':main()
