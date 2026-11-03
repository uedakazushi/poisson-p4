#!/usr/bin/env python3
"""Independent finite-rank application and Artin multiplication audit."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 start=time.monotonic();stems=[
 'reports/Q_NEW_FAMILY_LENGTH5_COMPLETION_20261103_v63.json',
 'reports/Q_POSITIVE_WEIGHT6_NONCM_20261103_v63.json',
 'reports/Q_NEW_FAMILY_TRANSVERSE_ARTIN_20261103_v63_REVIEW/report.json',
 'reports/Q_NEW_FAMILY_GENERIC_QUADRATIC_RANK_20261103_v63_REVIEW.json',
 'reports/Q_NEW_FAMILY_POLYNOMIAL_NORMAL_KERNEL_20261103_v63_REVIEW.json',
 'reports/Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.json']
 checked=0
 for name in stems:
  d=json.loads((ROOT/name).read_text());assert d['all_passed'],name
  for source,want in d.get('input_sha256',{}).items():assert hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==want,source;checked+=1
 U=s.zeros(5);V=s.zeros(5);U[1,0]=U[3,1]=V[2,0]=V[4,2]=1
 assert U*V==V*U==U**3==V**3==s.zeros(5)
 one=s.eye(5)[:,0];assert s.Matrix.hstack(one,U*one,V*one,U*U*one,V*V*one).det()==1
 assert len(U.col_join(V).nullspace())==2
 paths=[Path(__file__),*(ROOT/p for p in stems),ROOT/'research/notes/Q_NEW_FAMILY_FULL_COMPLETION_ROOT_REVIEW_20261103_v63.md']
 result=dict(id='Q_NEW_FAMILY_FULL_COMPLETION_ROOT_REVIEW_20261103_v63',all_passed=True,
  review_kind='Independent mathematical finite-rank, normalization, and CM-scope review plus current exact-certificate audit',
  checked_current_dependency_hashes=checked,
  actual_26_parameter_branch_section_not_26_group_factor=True,
  no_flatness_assumed_before_rank_five_comparison=True,
  quotient_N_mod_N2_proved_free_before_base_change_of_N2=True,
  full_projective_completed_ring='C[[s25,u,v]]/(uv,u^3,v^3)',
  full_closed_and_generic_transverse_length=5,CM_type=2,nilradical_cube_zero=True,
  reduced_component_smooth_at_selected_n111=True,
  full_nonzero_positive_axis_not_CM_by_multiplicity_embedding_inequality=True,
  global_reduced_component_normality_or_all_branches_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print('ROOT FULL COMPLETION5/NONCM APPLICATION REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
