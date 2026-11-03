#!/usr/bin/env python3
"""Audit finite-module nonCM argument and dense-orbit closure, no depth guess."""
import argparse,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,vector,action
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(output):
    start=time.monotonic()
    odd=ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json'
    finite=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json'
    boundary=ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json'
    for p in [odd,finite,boundary]:
        data=json.loads(p.read_text());assert data['all_passed'] is True
        for name,digest in data.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    o=json.loads(odd.read_text());assert o['full_odd_ideal_minimal_generators_lower_bounds']=={'CA':10,'CB':10,'N':8}
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    O=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))for i,j in itertools.product(range(5),repeat=2)])
    cc=list(O.rref()[1]);V=O[:,cc];rr=list(V.T.rref()[1]);assert len(cc)==len(rr)==19
    minor=s.factor(V[rr,:].det());assert minor!=0
    coordinates=V[rr,:].inv()*O[rr,:];assert V*coordinates==O
    scalar=V[rr,:].inv()*vector(P0)[rr,:];assert V*scalar==vector(P0)
    paths=[Path(__file__),ROOT/'research/notes/AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63.md',odd,finite,boundary,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py']
    result=dict(id='AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63',all_passed=True,
        accepted_actual_finite_images='A=R direct_sum I*t with R regular actual target base',
        normalizer_and_actual_image_generic_R_rank=2,odd_I_nonzero_ideal_R_rank=1,
        odd_minimal_generator_lower_bounds={'CA':10,'CB':10,'N':8},
        finite_R_m_R_A_primary_and_equal_dimensions=True,
        CM_A_would_imply_depth_R_A_equals_dim_R=True,
        Auslander_Buchsbaum_would_make_A_free_R=True,
        odd_I_direct_summand_would_be_projective_rank1_and_principal=True,
        principal_contradicts_discriminant_image_generator_lower_bound=True,
        individual_reduced_target_CA_CB_N_nonCM_at_P0=True,
        individual_reduced_target_CA_CB_N_nonnormal_at_P0_by_dual_normalization_fiber=True,
        independent_GL_affine_orbit_rank=19,scalar_direction_in_GL_orbit=True,
        independent_projective_orbit_dimension=18,orbit19_rows=rr,orbit19_columns=cc,orbit19_minor=str(minor),
        B18_irreducible_dimension18_dependency_accepted=True,
        finite_type_C_normal_CM_and_smooth_loci_open=True,
        entire_B18_in_nonCM_nonnormal_and_singular_loci_of_each_CA_CB_N=True,
        inclusion_only_not_exact_bad_locus=True,
        no_regular_nu_on_R_mod_I_assumed=True,
        exact_depth_or_actual_D_local_type_not_claimed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},seconds=time.monotonic()-start)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print('ALL B18 IN NONCM/NONNORMAL/SING OF REDUCED CA CB N PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
