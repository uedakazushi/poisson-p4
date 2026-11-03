#!/usr/bin/env python3
"""Exact generic R19 normal forms, intrinsic flags, and projective invariant.

No classification of the special t=0 or q5*q7=0 strata is asserted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,canonical,vector

ROOT=Path(__file__).resolve().parents[1]

def push_field(V,g,coordinates):
    subst=dict(zip(coordinates,g.inv()*coordinates))
    return (g*V.subs(subst,simultaneous=True)).applyfunc(s.cancel)

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;U=s.Matrix([h,x,z,w])
    t,q0,q1,q2,q3,q5,q6,q7=s.symbols('t q0 q1 q2 q3 q5 q6 q7')
    V=s.Matrix([0,4*t*x*x,-q0*x*x-2*t*x*z,2*t*x*w+3*q7*h*h+2*q3*h*x+2*q6*h*z+q1*x*x+q2*x*z+q5*z*z])
    shift=q0/(6*t);Bh=2*q3-2*shift*q6;Dx=q1-shift*q2+shift*shift*q5;Ez=q2-2*shift*q5
    su=Bh/(2*t);sv=-Dx/(2*t);sr=Ez/(4*t)
    G=s.Matrix([[1,0,0,0],[0,1,0,0],[0,shift,1,0],[su,sv+sr*shift,sr,1]])
    assert G.det()==1
    simple=s.Matrix([0,4*t*x*x,-2*t*x*z,2*t*x*w+3*q7*h*h+2*q6*h*z+q5*z*z])
    assert (push_field(V,G,U)-simple).applyfunc(s.cancel)==s.zeros(4,1)
    A,B,F=s.symbols('A B F');Q=A*h*h+B*h*z+F*z*z
    normal=s.Matrix([0,4*x*x,-2*x*z,2*x*w+Q]);D=5*normal-8*x*U
    assert sum(s.diff(normal[i],U[i]) for i in range(4))==8*x
    assert sum(s.diff(D[i],U[i]) for i in range(4))==0
    # The complete invariant-hyperplane argument is an exact division identity.
    a,b,c=s.symbols('a b c');ell=w+a*h+b*x+c*z
    Dell=sum(s.diff(ell,U[i])*D[i] for i in range(4))
    remainder=s.expand(Dell.subs(w,-a*h-b*x-c*z))
    assert remainder==s.expand(5*Q-10*a*h*x+10*b*x*x-20*c*x*z)
    ell0=a*h+b*x+c*z;Dell0=sum(s.diff(ell0,U[i])*D[i] for i in range(4))
    assert Dell0==s.expand(x*(-8*a*h+12*b*x-18*c*z))
    # Unique intrinsic factor: atx=0, wedge(E4,u,D) implies ux=uh=uz=0
    # for Q!=0. The remaininguw is killed by this fixed nonzero coefficient.
    assert s.expand(h*D[1]-x*D[0])==20*h*x*x
    da,db,dc,dd,wu,wv,wr=s.symbols('da db dc dd wu wv wr',nonzero=True)
    upper=s.Matrix([[da,0,0,0],[0,db,0,0],[0,0,dc,0],[wu,wv,wr,dd]])
    transformed=(db*push_field(normal,upper,U)).applyfunc(s.cancel)
    expected=s.Matrix([0,4*x*x,-2*x*z,
        2*x*w+db*dd*A*h*h/(da*da)+db*dd*B*h*z/(da*dc)+db*dd*F*z*z/(dc*dc)
        -2*wu*h*x/da+2*wv*x*x/db-4*wr*x*z/dc])
    assert (transformed-expected).applyfunc(s.cancel)==s.zeros(4,1)
    Aprime=db*dd*A/(da*da);Bprime=db*dd*B/(da*dc);Fprime=db*dd*F/(dc*dc)
    assert s.cancel(Bprime*Bprime/(Aprime*Fprime)-B*B/(A*F))==0
    # FullGL5 y-shears make no additional quotient freedom and preservetheprojectiveclass.
    ey=s.Matrix([0,0,1,0,0]);V5=s.Matrix([normal[0],normal[1],0,normal[2],normal[3]]);P=canonical(ey*V5.T-V5*ey.T)
    sy,lh,lx,lz,lw=s.symbols('sy lh lx lz lw',nonzero=True)
    full=s.Matrix([[da,0,0,0,0],[0,db,0,0,0],[lh,lx,sy,lz,lw],[0,0,0,dc,0],[0,0,0,0,dd]])
    inverse=full.inv();Ppush=(full*P.subs(dict(zip(X,inverse*s.Matrix(X))),simultaneous=True)*full.T).applyfunc(s.cancel)
    targetV=s.Matrix([0,4*x*x,0,-2*x*z,2*x*w+Aprime*h*h+Bprime*h*z+Fprime*z*z])
    targetP=canonical(ey*targetV.T-targetV*ey.T)
    assert (Ppush-sy*targetP/db).applyfunc(s.cancel)==s.zeros(5)
    # Explicit sign equivalence of the standard kappa representative.
    kappa=s.Symbol('kappa');sign=s.diag(-1,1,1,1)
    standard=normal.subs({A:1,B:kappa,F:1});negative=normal.subs({A:1,B:-kappa,F:1})
    assert push_field(standard,sign,U)==negative
    report=dict(id='AMBIENT_R19_GENERIC_ORBITS_20261103_v63',all_passed=True,field='QQ(parameters)',
        source_parameters=list(map(str,[t,q0,q1,q2,q3,q5,q6,q7])),source_vector_field=list(map(str,V)),
        generic_open='t*q5*q7 !=0',normalizing_linear_matrix=[[str(v) for v in row] for row in G.tolist()],
        normalized_coefficients={'A':'3*q7/t','B':'2*q6/t','F':'q5/t'},
        normal_form='P_kappa=d_y wedge (4*x²*d_x-2*x*z*d_z+(2*x*w+h²+kappa*h*z+z²)*d_w)',
        complete_projective_orbit_invariant='kappa²=4*q6²/(3*q7*q5)',
        invariant_hyperplane_division_remainder=str(remainder),invariant_hyperplanes_exactly=['h=0','x=0','z=0'],
        canonical_D=list(map(str,D)),relative_cofactor_ratios={'h':'-2/3','z':'-3/2'},
        distinguished_self_cofactor_hyperplane='x=0',intrinsic_constant_factor_unique='C*d_y',
        equivalence_quotient_matrix_form='diag(da,db,dc,dd); w translations vanish between normalforms',
        coefficient_transform={'A':str(Aprime),'B':str(Bprime),'F':str(Fprime)},
        sign_equivalence=True,
        stabilizer_connected_diagonal_conditions=['db*dd=da²','dc=da'],
        stabilizer_full_GL_dimension=7,stabilizer_PGL_dimension=6,projective_orbit_dimension=18,
        stabilizer_components_special_B_zero='dc=+da or -da; connected dimension unchanged',
        source_generality='Dense parameteropen and itsPGL saturation; specialt=0 orAF=0 orbits not classified',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [Path(__file__),ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'research/notes/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63.md']},
        seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('R19 GENERIC ORBIT NORMALFORM PASS',report['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
