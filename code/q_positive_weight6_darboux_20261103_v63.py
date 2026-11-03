#!/usr/bin/env python3
"""Whole positive ray: polynomial Darboux map and (2,3) Poisson Cremona map."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from s61_model import X as HOM,bivector_from_coefficients
from q51_universal_section_20261101 import affine

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def strings(M):return [[str(v)for v in row]for row in M.tolist()]
def zero(M):return M.applyfunc(s.expand)==s.zeros(M.rows,M.cols)

def build(out):
    started=time.monotonic();eta=s.Symbol('eta')
    source=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz'
    report_path=source.parent/'report.json';old=json.loads(report_path.read_text());assert old['all_passed']and sha(source)==old['certificate_sha256']
    saved=json.loads(gzip.decompress(source.read_bytes()));P=bivector_from_coefficients(list(map(s.sympify,saved['raw150'])))
    x,y,z,w=s.symbols('x y z w');V=s.Matrix([x,y,z,w]);Paff=affine(P,V)
    phi=s.Matrix([x,y-s.Rational(3,4)*x*x+48*eta*w*w,z-s.Rational(3,4)*x*y,w])
    target=s.zeros(4);target[0,3]=s.Rational(8,9);target[1,2]=-s.Rational(2,9);target=target-target.T
    Jac=phi.jacobian(V);assert Jac.det()==1
    assert zero(Jac*Paff*Jac.T-target)
    XX,YY,ZZ,WW=s.symbols('X Y Z W');U=s.Matrix([XX,YY,ZZ,WW])
    inverse=s.Matrix([XX,YY+s.Rational(3,4)*XX*XX-48*eta*WW*WW,
       ZZ+s.Rational(3,4)*XX*YY+s.Rational(9,16)*XX**3-36*eta*XX*WW*WW,WW])
    assert zero(inverse.subs(dict(zip(U,phi)),simultaneous=True)-V)
    assert zero(phi.subs(dict(zip(V,inverse)),simultaneous=True)-U)
    assert inverse.jacobian(U).det()==1
    assert zero(inverse.jacobian(U)*target*inverse.jacobian(U).T-Paff.subs(dict(zip(V,inverse)),simultaneous=True))
    Hamil=s.Rational(9,2)*y-s.Rational(27,8)*x*x+216*eta*w*w
    assert zero(Paff*s.Matrix([s.diff(Hamil,v)for v in V])-s.Matrix([0,0,1,0]))

    h,xx,yy,zz,ww=HOM
    F=s.Matrix([h*h,h*xx,h*yy-s.Rational(3,4)*xx*xx+48*eta*ww*ww,h*zz-s.Rational(3,4)*xx*yy,h*ww])
    HH,TX,TY,TZ,TW=s.symbols('H Xh Yh Zh Wh');T=s.Matrix([HH,TX,TY,TZ,TW])
    G=s.Matrix([HH**3,TX*HH**2,TY*HH**2+s.Rational(3,4)*TX*TX*HH-48*eta*TW*TW*HH,
       TZ*HH**2+s.Rational(3,4)*TX*TY*HH+s.Rational(9,16)*TX**3-36*eta*TX*TW*TW,TW*HH**2])
    assert zero(G.subs(dict(zip(T,F)),simultaneous=True)-h**5*s.Matrix(HOM))
    assert zero(F.subs(dict(zip(HOM,G)),simultaneous=True)-HH**5*T)
    assert [s.Poly(v,*HOM).total_degree()for v in F]==[2]*5
    assert [s.Poly(v,*T).total_degree()for v in G]==[3]*5
    assert s.gcd_list(list(F))==1 and s.gcd_list(list(G))==1
    # This gcd statement persists for EVERY eta specialization: h is the
    # only candidate common factor and fixed x^2 / X^3 terms exclude it.
    assert s.Poly(F[2],h).coeff_monomial(1).coeff(xx,2)==-s.Rational(3,4)
    assert s.Poly(G[3],HH).coeff_monomial(1).coeff(TX,3)==s.Rational(9,16)
    constant=s.zeros(5);constant[1,4]=s.Rational(8,9)*HH**2;constant[2,3]=-s.Rational(2,9)*HH**2;constant=constant-constant.T
    assert zero(F.jacobian(HOM)*P*F.jacobian(HOM).T-constant.subs(dict(zip(T,F)),simultaneous=True))
    assert zero(G.jacobian(T)*constant*G.jacobian(T).T-P.subs(dict(zip(HOM,G)),simultaneous=True))
    out.mkdir(parents=True,exist_ok=False)
    certificate=dict(affine_variables=list(map(str,V)),affine_Darboux_map=list(map(str,phi)),affine_inverse_variables=list(map(str,U)),affine_inverse=list(map(str,inverse)),
       full_affine_source_P=str(Paff),full_affine_target_matrix=str(target),Jacobian= strings(Jac),Jacobian_determinant='1',translation_Hamiltonian=str(Hamil),
       homogeneous_source_variables=list(map(str,HOM)),homogeneous_forward=list(map(str,F)),homogeneous_target_variables=list(map(str,T)),homogeneous_inverse=list(map(str,G)),
       homogeneous_target_constant_bivector=str(constant),forward_inverse_composition_scalar='x0^5',inverse_forward_composition_scalar='H^5',
       degree_two_gcd_one_uniform_witness='F0=h^2 and F1=h*x force any common factor to divide h; coefficient x^2 of F2 modh=-3/4',
       degree_three_gcd_one_uniform_witness='G0=H^3 and G1=Xh*H^2 force any common factor to divide H; coefficient Xh^3 of G3 modH=9/16')
    (out/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'research/notes/Q_POSITIVE_WEIGHT6_DARBOUX_20261103_v63.md',source,report_path,ROOT/'code/s61_model.py',ROOT/'code/q51_universal_section_20261101.py']
    result=dict(id='Q_POSITIVE_WEIGHT6_DARBOUX_20261103_v63',all_passed=True,field='QQ[eta]',parameter_scope='Every complex eta, including0',
       actual_positive_ray_full150_input=True,affine_polynomial_Darboux_map=True,affine_forward_inverse_degrees=[2,3],affine_forward_and_inverse_Jacobian_determinants='1',
       all_six_Darboux_brackets_exact=True,all_forward_and_inverse_polynomial_compositions_exact=True,unique_translation_has_quadratic_Hamiltonian=str(Hamil),
       homogeneous_birational_bidegree=[2,3],homogeneous_both_compositions_scalar_fifth_power_identity=True,homogeneous_both_coordinate_gcds_one_at_every_eta=True,
       both_full_homogeneous_Poisson_pushforward_identities_exact=True,
       PGL_equivalence_or_known_component_membership_inferred=False,new_component_or_exhaustion_claimed=False,
       input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs},certificate_sha256=sha(out/'certificate.json'),seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('POSITIVE WHOLE RAY POLYNOMIAL DARBOUX + (2,3) POISSON CREMONA PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args().output)
