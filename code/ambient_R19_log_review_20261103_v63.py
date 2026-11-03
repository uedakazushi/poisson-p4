#!/usr/bin/env python3
"""Independent actual legacy log tensor, rational Poisson map and leaf audit."""
import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import IDS, canonical

ROOT = Path(__file__).resolve().parents[1]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text())

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replay',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    report=ROOT/'reports/AMBIENT_R19_LOG_AND_BIRATIONAL_20261103_v63/report.json'
    d=load(report);fresh=load(args.replay)
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
    paths=[Path(__file__),report,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
    for file,hsh in d['input_sha256'].items():
        p=ROOT/file;assert sha(p)==hsh,file;paths.append(p)
    h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);A,B,F=s.symbols('A B F')
    Q=A*h*h+B*h*z+F*z*z;V=s.Matrix([0,4*x*x,0,-2*x*z,2*x*w+Q]);ey=s.Matrix([0,0,1,0,0]);P=ey*V.T-V*ey.T
    K=120*x*w+20*A*h*h+15*B*h*z+12*F*z*z
    deriv=lambda f:s.cancel(sum(V[i]*s.diff(f,X[i])for i in range(5)))
    assert s.expand(deriv(K)-6*x*K)==0
    phi=x*z*z/h**3;psi=z*K/(120*x*h*h)
    assert deriv(phi)==deriv(psi)==0 and s.diff(phi,y)==s.diff(psi,y)==0
    # Call the original 2026-11-01 family constructor, not the new producer's
    # reconstructed exterior-form implementation.  No component-rank rerun.
    legacy=ROOT/'legacy/v14/code/rank2_lower_bounds'
    sys.path.insert(0,str(legacy))
    spec=importlib.util.spec_from_file_location('legacy_log_actual_builder',legacy/'verify_rank2_components.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    divisors=[h,x,z,K]
    polys=[dict(s.Poly(f,*X).terms())for f in divisors]
    old=module.log_bivector(polys,[s.Rational(5,4),s.Rational(1,4),s.Rational(3,4)])
    oldvector=s.zeros(150,1)
    rawvector=s.zeros(150,1)
    for (ij,m),index in IDS.items():
        powers=tuple(m.count(k)for k in range(5))
        oldvector[index]=old[ij].get(powers,0)
        rawvector[index]=s.Poly(P[ij[0],ij[1]],*X).coeff_monomial(s.prod(X[k]for k in m))
    assert(oldvector+300*canonical(rawvector)).applyfunc(s.expand)==s.zeros(150,1)
    ra=s.Matrix([-3,1,2,0]);rb=s.Matrix([-2,-1,1,1]);deg=s.Matrix([1,1,1,2])
    assert ra.dot(deg)==rb.dot(deg)==0
    assert [int((ra*rb.T-rb*ra.T)[i,j])for i,j in itertools.combinations(range(4),2)]==[5,1,-3,3,1,2]
    # Exact rational inverse, full five-dimensional pushforward, Euler descent.
    change=s.Matrix([h,x,y,z,K/(120*x)])
    inverse={w:w-A*h*h/(6*x)-B*h*z/(8*x)-F*z*z/(10*x)}
    assert(change.subs(inverse,simultaneous=True)-X).applyfunc(s.cancel)==s.zeros(5,1)
    J=change.jacobian(X)
    assert s.cancel(J.det()-1)==0 and(J*X-change).applyfunc(s.cancel)==s.zeros(5,1)
    V0=s.Matrix([0,4*x*x,0,-2*x*z,2*x*w])
    assert((J*P*J.T).subs(inverse,simultaneous=True)-(ey*V0.T-V0*ey.T)).applyfunc(s.cancel)==s.zeros(5)
    r,t,u,Y=s.symbols('r t u Y')
    parametrization={h:1,x:r/u**2,y:Y,z:u,w:t/u-A*u*u/(6*r)-B*u**3/(8*r)-F*u**4/(10*r)}
    assert s.cancel(phi.subs(parametrization,simultaneous=True)-r)==0
    assert s.cancel(psi.subs(parametrization,simultaneous=True)-t)==0
    recovered={r:phi.subs(h,1),t:psi.subs(h,1),u:z,Y:y}
    for coordinate in [x,y,z,w]:
        assert s.cancel(parametrization[coordinate].subs(recovered,simultaneous=True)-coordinate)==0
    # The general projective leaf is the cone over this basepoint-free sextic.
    U,T=s.symbols('U T')
    curve=[U*U*T**4,r*T**6,U**3*T**3,t*U*T**5-A*U**4*T**2/(6*r)-B*U**5*T/(8*r)-F*U**6/(10*r)]
    assert all(s.Poly(f,U,T).total_degree()==6 for f in curve)
    assert s.cancel(curve[2]/curve[0]-U/T)==0
    assert curve[1].subs(U,0)==r*T**6 and curve[3].subs(T,0)==-F*U**6/(10*r)
    # Every common zero has T=0 from the second coordinate; the fourth then
    # excludes U!=0 when r*F!=0.  The inverse U/T proves degree-one onto image.
    oldcert=d['existing_Log1112_component_certificate']
    assert oldcert['affine_jacobi_exact_zero'] and oldcert['J_vertical_exact_zero'] and oldcert['J_parameter_derivative_exact_zero']
    assert oldcert['family_affine_differential_rank']==oldcert['ambient_affine_tangent_dimension']==29
    assert oldcert['projective_component_dimension']==28
    paths += [legacy/'poisson_p4_components_check.py',legacy/'lpb_family.py',
        ROOT/'reports/AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_R19_GEOMETRY_20261103_v63_REVIEW.json',
        ROOT/'research/notes/AMBIENT_R19_LOG_AND_BIRATIONAL_20261103_v63.md',
        ROOT/'research/notes/AMBIENT_R19_LOG_AND_BIRATIONAL_20261103_v63_REVIEW.md']
    result=dict(id='AMBIENT_R19_LOG_AND_BIRATIONAL_20261103_v63_REVIEW',all_passed=True,
        fresh_QQ_producer_matches=True,original_legacy_log_builder_all150_equals_minus300_canonical_P=True,
        old_full_Log1112_component_certificate_reused_not_rerun=True,
        entire_closed_R_contained_in_existing_Log1112_28=True,
        independent_actual_rational_Poisson_pushforward_Euler_and_inverse=True,
        all_R_minus_N_birationally_Poisson_equivalent_to_Q_zero=True,
        independent_two_Casimirs_and_full_function_field_inverse=True,
        generic_leaves_rational_for_all_R_minus_N=True,
        generic_rF_nonzero_leaf_closure='degree6 rational cone, vertex C*d_y',
        birational_or_leaf_claim_for_entire_R_intersection_N=False,new_component_count=0,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print('PASS original legacy Log1112 tensor, Cremona and rational leaves')

if __name__=='__main__':main()
