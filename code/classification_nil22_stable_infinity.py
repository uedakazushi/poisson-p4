#!/usr/bin/env python3
"""Whole-source changes of infinity for both stable nil22 image sources.

All six magnetic variables, all 150 homogeneous coefficients, and every
simple linear factor are retained. No parameter samples or deformations.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import resource
import time

import sympy as s
from classification_nil22_image_source_review import build_P,coords,p,F0,mag,x,y,q1,q2

ROOT=Path(__file__).resolve().parents[1]


def jacobi(P,variables):
    return [s.expand(sum(P[i,l]*s.diff(P[j,k],variables[l])+P[j,l]*s.diff(P[k,i],variables[l])+P[k,l]*s.diff(P[i,j],variables[l])for l in range(len(variables))))for i,j,k in itertools.combinations(range(len(variables)),3)]


def encode(P):
    return {f'{i},{j}':str(s.factor(P[i,j]))for i,j in itertools.combinations(range(P.rows),2)if P[i,j]!=0}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    h,H,Q1,Q2,X,Y=s.symbols('h H Q1 Q2 X Y');old=s.Matrix([h,q1,q2,x,y]);new=s.Matrix([H,Q1,Q2,X,Y])
    specs=[('stable_line_conic',s.Matrix([[3,-15],[-15,3]]),s.Matrix([s.Rational(3,2),-s.Rational(3,2)]),[(x-y-12*h,y)]),
        ('stable_three_lines',s.Matrix([[-6,3],[3,-6]]),s.zeros(2,1),[(x+3*h,y),(y-3*h,x),(x-y-3*h,y)])]
    records=[]
    for name,S,V0,factors in specs:
        P,M,D=build_P(x*y*(x-y),S,V0,-s.eye(2)/3)
        P5=s.zeros(5)
        for i,j in itertools.combinations(range(4),2):
            homogeneous=sum(co*h**(2-sum(e))*s.prod(z**k for z,k in zip(coords,e))for e,co in s.Poly(P[i,j],*coords).terms())
            P5[i+1,j+1]=homogeneous;P5[j+1,i+1]=-homogeneous
        Dhom=s.expand(sum(co*h**(3-sum(e))*x**e[0]*y**e[1]for e,co in s.Poly(D,x,y).terms()))
        actual_D=s.factor(-h*h*Dhom)
        # Whole magnetic action q'=q+L(p), solved without deleting a coefficient.
        ell=s.symbols('l0:6');L=s.Matrix([ell[0]+ell[1]*x+ell[2]*y,ell[3]+ell[4]*x+ell[5]*y])
        moved=s.Matrix([q1+L[0],q2+L[1],x,y]);J=moved.jacobian(s.Matrix(coords))
        shifted=(J*P*J.T).subs({q1:q1-L[0],q2:q2-L[1]},simultaneous=True).applyfunc(s.expand)
        change=s.Poly(shifted[0,1]-P[0,1],x,y)
        mons=[1,x,y,x*x,x*y,y*y]
        action=s.Matrix([change.coeff_monomial(mon)for mon in mons]).jacobian(ell)
        assert all(s.expand(shifted[i,j]-P[i,j])==0 for i,j in itertools.combinations(range(4),2)if(i,j)!=(0,1))
        source=dict(name=name,all_six_magnetic_parameters=list(map(str,mag)),actual_anticanonical_divisor=str(actual_D),base_D=str(s.factor(Dhom)),
            full_magnetic_action_matrix=[list(map(str,row))for row in action.tolist()],magnetic_action_rank=action.rank(),
            full_magnetic_cokernel=[list(map(str,v))for v in action.T.nullspace()],simple_infinity=[])
        for linear,other in factors:
            image=s.Matrix([linear,q1,q2,h,other]);G=image.jacobian(old);assert G.det()!=0
            raw=(G*P5*G.T).subs(dict(zip(old,G.inv()*new)),simultaneous=True).applyfunc(s.expand)
            R=s.Matrix([0,*[s.cancel(raw[0,j]/H)for j in range(1,5)]])
            assert all(not s.denom(value).has(H)for value in R)
            assert all(not value.has(*mag)for value in R)
            transformed=(raw-(new*R.T-R*new.T)).applyfunc(s.expand)
            assert all(transformed[0,j]==0 for j in range(5))
            # Explicit inverse reconstruction proves all 150 quotient coefficients.
            recovered=(G.inv()*(transformed+new*R.T-R*new.T)*G.inv().T).subs(dict(zip(new,G*old)),simultaneous=True).applyfunc(s.expand)
            assert recovered==P5
            affine=transformed[1:,1:].subs(H,1)
            variables=(Q1,Q2,X,Y)
            assert not any(jacobi(affine,variables))
            Z=s.Matrix([s.expand(sum(s.diff(affine[i,j],variables[i])for i in range(4)))for j in range(4)])
            A=Z.jacobian(s.Matrix(variables));c=Z.subs(dict.fromkeys(variables,0))
            assert not any(value.has(*mag,*variables)for value in A)
            assert not any(value.has(*mag,*variables)for value in c)
            assert (Z-A*s.Matrix(variables)-c).applyfunc(s.expand)==s.zeros(4,1)
            polynomial=s.factor(A.charpoly(s.Symbol('lambda')).as_expr())
            zero_set=s.linsolve(list(Z),variables)
            zero=None;value_pf=None;C=None
            if zero_set is not s.EmptySet and bool(zero_set):
                representative=list(zero_set)[0]
                free=set().union(*(z.free_symbols for z in representative))
                zero=s.Matrix([z.subs(dict.fromkeys(free,0))for z in representative])
                C=affine.subs(dict(zip(variables,zero)),simultaneous=True).applyfunc(s.expand)
                value_pf=s.factor(C[0,1]*C[2,3]-C[0,2]*C[1,3]+C[0,3]*C[1,2])
            square_scalar=s.factor(s.trace(A*A)/4)
            opposite_pair=square_scalar!=0 and (A*A-square_scalar*s.eye(4)).applyfunc(s.expand)==s.zeros(4)
            record=dict(simple_factor=str(linear),other_base_coordinate=str(other),actual_PGL_matrix=G.tolist(),determinant=str(G.det()),
                whole_Euler_correction=list(map(str,R)),all150_inverse_reconstruction=True,whole_affine_six_brackets=encode(affine),
                full_source_strict_Jacobi_zero=True,new_modular_field=list(map(str,Z)),new_modular_matrix=A.tolist(),new_modular_translation=list(map(str,c)),
                modular_characteristic_polynomial=str(polynomial),modular_rank=A.rank(),modular_eigenspace_dimensions={str(k):4-(A-k*s.eye(4)).rank()for k in A.eigenvals()},
                opposite_eigenvalue_pairs=bool(opposite_pair),square_scalar=str(square_scalar),
                modular_zero=list(map(str,zero))if zero is not None else None,full_constant_value=encode(C)if C is not None else None,
                constant_Pf_at_modular_zero=str(value_pf)if value_pf is not None else None,
                accepted_whole_SS20_rank2_hypotheses=bool(opposite_pair and value_pf is not None and value_pf!=0),
                parameters_localized=[])
            source['simple_infinity'].append(record)
        records.append(source)
    paths=[Path(__file__),ROOT/'code/classification_nil22_image_source_review.py',ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_SOURCE_REVIEW_V2_20261107.json',ROOT/'research/classification/HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md']
    result=dict(task_id='E_multiple_hyperplane',all_passed=True,whole_source_changes_of_infinity=records,
        source_parameters_dropped=False,no_parameter_localization=True,
        entire_sources_MPS_covered=all(any(v['accepted_whole_SS20_rank2_hypotheses']for v in r['simple_infinity'])for r in records),
        scope='Only apply a previously proved whole-source theorem when the recorded exact modular and constant-value hypotheses hold. A negative flag is not a nonmembership conclusion.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths},producer_source_bytes=Path(__file__).read_text(),
        seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2,default=str);handle.write('\n')
    print(json.dumps(dict(all_passed=True,entire_sources_MPS_covered=result['entire_sources_MPS_covered'],seconds=result['seconds'])))


if __name__=='__main__':main()
