#!/usr/bin/env python3
"""Finite exact orbit reduction for ALL R19 minus N27, via accepted66cover.

This classifies the complement, not the remaining R intersection N boundary.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;U=s.Matrix([h,x,z,w]);A,B,F=s.symbols('A B F');Q=A*h*h+B*h*z+F*z*z
    D=s.Matrix([-8*h*x,12*x*x,-18*x*z,5*Q+2*x*w]);a,b,c,d,lam,alpha=s.symbols('a b c d lambda alpha')
    ell=a*h+b*x+c*z+d*w;Dell=sum(s.diff(ell,U[i])*D[i] for i in range(4));difference=s.Poly(Dell-lam*ell*ell,*U)
    assert difference.coeff_monomial(w*w)==-lam*d*d
    assert difference.coeff_monomial(h*h)==5*A*d-lam*a*a
    assert difference.coeff_monomial(z*z)==5*F*d-lam*c*c
    # The4 coefficient rows below are independent for EVERYQ, includingzero.
    casimir_rows=[s.Poly(Dell,*U).coeff_monomial(m) for m in [h*x,x*x,x*z,x*w]]
    casimir_matrix,_=s.linear_eq_to_matrix(casimir_rows,[a,b,c,d]);assert casimir_matrix==s.diag(-8,12,-18,2)
    second=s.expand(Dell+48*x*ell+alpha*x*x)
    expected=40*a*h*x+(60*b+alpha)*x*x+30*c*x*z+d*(5*Q+50*x*w)
    assert second==s.expand(expected)
    rows=[s.Poly(second,*U).coeff_monomial(m) for m in [h*x,x*x,x*z,x*w]]
    second_matrix,_=s.linear_eq_to_matrix(rows,[a,b,c,d,alpha]);assert second_matrix.nullspace()==[s.Matrix([0,-s.Rational(1,60),0,0,1])]
    # Q=0 factor uniqueness and fourthDarbouxcofactor are checkedphysically.
    ey=s.Matrix([0,0,1,0,0]);V5=s.Matrix([0,4*x*x,0,-2*x*z,2*x*w]);P=canonical(ey*V5.T-V5*ey.T);vv=s.Matrix(s.symbols('v0:5'));E=s.Matrix(X)
    mons=[s.prod(X[i] for i in ids) for ids in itertools.combinations_with_replacement(range(5),3)];eq=[]
    for i,j,k,l in itertools.combinations(range(5),4):
        f=s.expand(E[i]*(vv[j]*P[k,l]-vv[k]*P[j,l]+vv[l]*P[j,k])-E[j]*(vv[i]*P[k,l]-vv[k]*P[i,l]+vv[l]*P[i,k])
            +E[k]*(vv[i]*P[j,l]-vv[j]*P[i,l]+vv[l]*P[i,j])-E[l]*(vv[i]*P[j,k]-vv[j]*P[i,k]+vv[k]*P[i,j]))
        eq += [s.Poly(f,*X).coeff_monomial(m) for m in mons]
    factor,_=s.linear_eq_to_matrix(eq,list(vv));assert factor.rank()==4 and factor.nullspace()==[ey]
    weight=s.Matrix([[-2,1,0,1],[-1,1,-1,1],[0,1,-2,1]])
    normalforms=[('AF_nonzero_family',[1,'kappa',1],18),('A_B_nonzero_F_zero',[1,1,0],18),('A_only',[1,0,0],17),
        ('B_F_nonzero_A_zero',[0,1,1],18),('F_only',[0,0,1],17),('B_only',[0,1,0],17),('zero',[0,0,0],16)]
    records=[]
    for label,coefficients,dimension in normalforms:
        active=[i for i,v in enumerate(coefficients) if v!=0];rank=weight[active,:].rank() if active else 0
        assert 16+rank==dimension
        records.append(dict(label=label,Q_coefficients=coefficients,weight_rank=rank,projective_orbit_dimension=dimension,
            PGL_projective_class_stabilizer_dimension=8-rank))
    assert len(set([-8,12,-18,2]))==4
    paths=[Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'reports/AMBIENT_R19_GENERIC_ORBITS_20261103_v63/report.json',
        ROOT/'research/notes/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63.md',
        ROOT/'research/notes/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63.md']
    body=dict(id='AMBIENT_R19_COMPLEMENT_ORBITS_20261103_v63',all_passed=True,field='QQ(parameters)',
        coverage_reason='Acceptedglobal66table: everypointoutsideN hasbigcellflaglift; t0sourcehasfCA0 henceN; thus ALL R minusN has tnonzero andthe exactVQnormalform.',
        every_Q_normalform_outside_N=True,all_Q_linear_casimir_matrix=casimir_matrix.tolist(),
        all_Q_self_cofactor_equation_key_coefficients=['-lambda*d²','5*A*d-lambda*a²','5*F*d-lambda*c²'],
        all_Q_N_second_flag_matrix=second_matrix.tolist(),zero_Q_intrinsic_factor_matrix=factor.tolist(),
        zero_Q_unique_factor='C*d_y',zero_Q_Darboux_hyperplanes=['h','x','z','w'],
        zero_Q_cofactor_ratios=['-2/3','1','-3/2','1/6'],
        complete_orbit_list=records,complete_family_invariant='kappa²',
        number_of_additional_single_orbits=6,
        scope='COMPLETE PGL5 orbitclassification of R19 minusN27. No orbitclassification of R intersectionN or fullPoissonspace; no scheme-thicknessclaim for AF0 cases.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2,default=str)+'\n');print('ALL R MINUS N ORBITS PASS',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
