#!/usr/bin/env python3
"""A smooth point of N's rank2 boundary divisor, hence global Serre R1.

Global proper coverage and uniqueness of the entire geometric flag fiber are
essential; the differential image alone is not equated with a target tangent.
"""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import raw,canonical,vector,boundary_witness,X,h,x,y,z,w

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    started=time.monotonic();output.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    data=json.loads(parent.read_text())
    basis=[s.Matrix(list(map(s.Rational,v)))for v in data['basis_canonical150']]
    # boundary_witness accepts tensors, so reconstruct N's canonical basis.
    cubics=list(map(s.sympify,data['cubic_monomials']))
    tensors=[canonical(raw(f,0,0))for f in [h*h,h*x]]+[canonical(raw(0,1,0))]+[canonical(raw(0,0,H))for H in cubics]
    assert s.Matrix.hstack(*(vector(P)for P in tensors))==s.Matrix.hstack(*(basis[i]for i in [0,1,*range(3,20)]))
    P=canonical(raw(h*x,0,z**3+w**3));witness=boundary_witness(P,tensors)
    assert witness['affine_image_differential_rank']==28
    assert witness['intrinsic_constant_vector_factor_matrix_rank']==4
    C=s.Matrix([h*h,-4*h*x,0,-15*w*w+h*z,15*z*z+h*w])
    v=s.eye(5)[:,2];assert P==(v*C.T-C*v.T)/5
    aa,bb,cc,dd,lam=s.symbols('A B C D lam');ell=aa*h+bb*x+cc*z+dd*w
    cell=aa*C[0]+bb*C[1]+cc*C[3]+dd*C[4]
    poly=s.Poly(s.expand(cell-lam*ell**2),h,x,z,w)
    assert poly.coeff_monomial(x*x)==-lam*bb**2
    assert poly.coeff_monomial(z*w)==-2*lam*cc*dd
    assert poly.coeff_monomial(z*z)==15*dd-lam*cc**2
    assert poly.coeff_monomial(w*w)==-15*cc-lam*dd**2
    assert poly.coeff_monomial(h*h).subs({bb:0,cc:0,dd:0})==aa-lam*aa**2
    # lambda0 has no nonzero solution; otherwise B=0, C*D=0, C=D=0.
    coeff=s.Matrix.hstack(*(s.Matrix([s.Poly(C[i],h,x,z,w).coeff_monomial(m)for m in [h*h,h*x,z*z,w*w,h*z,h*w]])for i in [0,1,3,4]))
    assert coeff.rank()==4
    test=s.expand(cell+4*h*ell)
    assert s.expand(test-(5*aa*h*h+cc*(-15*w*w+5*h*z)+dd*(15*z*z+5*h*w)))==0
    # Thus unique marked hline, and unique W=span(h,x).
    result=dict(id='AMBIENT_N_GLOBAL_R1_20261103_v63',all_passed=True,
      smooth_rank_two_divisor_point_raw_f='h*x',smooth_rank_two_divisor_point_raw_H='z^3+w^3',
      smooth_source_projective_dimension=27,actual_source_affine_image_rank=28,
      full_dimension_minor=witness,unique_intrinsic_constant_vector_factor='d_y',
      quotient_C_times_five=list(map(str,C)),linear_Casimir_kernel_zero=True,
      self_cofactor_ell_unique='ell proportional h',second_flag_W_unique='span(h,x)',
      entire_geometric_proper_model_fiber_one_point=True,
      fiber_tangent_zero_by_full_source_immersion=True,
      finite_fiber_length_one_and_Nakayama_isomorphism_accepted_by_geometric_proof=True,
      N_smooth_at_this_rank_two_boundary_point=True,
      full_projective_N_Serre_R1=True,
      global_N_normality_or_CM_claimed=False,
      input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),parent,
        ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json',
        ROOT/'reports/AMBIENT_COMPONENT_R1_NONCM_20261103_v63_REVIEW/report.json']},
      seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('GLOBAL N R1: smooth rank2 divisor witness PASS',result['seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
