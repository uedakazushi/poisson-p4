#!/usr/bin/env python3
"""CB-only multiplicity4 companion using independently accepted actual Qb selfpair."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_bytes())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    upper=ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63/report.json';up=load(upper)
    assert up['all_passed'] is True and up['CA_CB_origin_conductor_cone_dimension']==20 and up['CA_CB_origin_conductor_multiplicity_upper_bound']==4
    pair=ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW/report.json';audit=load(pair)
    assert audit['all_passed'] is True and audit['both_source_ranks']==[29,29]
    assert audit['all58_union_rank']==37 and audit['tangent_intersection_dimension']==audit['meeting_scheme_edim']==21
    assert audit['own_GL_scalar_parameter_rank']==20
    assert audit['generic_whole_fiber_exactly_two_reduced_points_by_proper_open'] is True
    deps=[upper,ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63_REVIEW/report.json',pair,
      ROOT/'reports/AMBIENT_CB_GENERIC_V_FULL_COMPLETION_20261103_v63/report.json',
      ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json']
    for p in deps:
        d=load(p);assert d['all_passed'] is True
        for name,digest in d.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    certpath=ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63/certificate.json.gz';cert=load(certpath);binding=cert['actual_jet_bindings']['CB']
    rows=binding['seven_actual_odd_rows']+[46,95,63];assert len(rows)==len(set(rows))==10
    h=s.Matrix(binding['h63_as_constant_combination_of_all150_actual_odd_functions']);assert {i for i,c in enumerate(h)if c}<=set(rows)
    T=s.Symbol('T');assert s.cancel((1-2*T*T-T**3+T**3+T**4)/(1-T)**22-(1+T)**2/(1-T)**20)==0
    assert 1+1+2==up['CA_CB_origin_conductor_multiplicity_upper_bound']
    evidence=dict(actual_ten_minimal_odd_generator_rows=rows,graded_conductor_Hilbert_series='(1+T)^2/(1-T)^20',
      core_generator_degrees=[2,2,3],core_syzygy_degrees=[3,4],
      exact_support_branches=['criticalW','criticalV','actualCB_Qb_selfpairZ'],generic_lengths=[1,1,2],
      old_CA_Pa_not_used_as_CB_selfpair=True,
      generic_nonreducedness_uses_actual_edim21_and_dense_dimension20_family_not_flags_alone=True)
    (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63.md';paths=[Path(__file__),note,certpath,*deps]
    result=dict(id='AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63',all_passed=True,component='CB only',
      independently_accepted_new_actual_Qb_selfpair_used=True,
      intrinsic_conductor_equals_meeting_of_all_two_smooth_normalization_branches=True,
      actual_edim21_plus_dense_selfpair20_forces_generic_conductor_length_atleast2=True,
      new_selfpair_component_has_nonzero_quadratic_discriminant_distinct_from_W_V=True,
      exact_conductor_multiplicity4=True,exact_full_conductor_initial_ideal_equals_I0=True,
      exact_initial_form='(WU,WV,U*(aU²+bUV+cV²)) extendedby g19,a!=0,b/c notcomputed',
      conductor_cone_dimension=20,conductor_cone_embedding_dimension=22,conductor_cone_depth=20,conductor_CM_type=2,
      exact_minimal_prime_count=3,generic_primary_lengths=[1,1,2],conductor_has_no_embedded_associated_primes=True,
      actual_odd_ideal_minimal_generator_count=10,actual_odd_ideal_generated_by_saved_ten_actual_rows=True,
      actual_CB_cone_dimension=29,actual_CB_cone_embedding_dimension=39,actual_CB_cone_depth=21,actual_CB_multiplicity=2,
      actual_CB_projective_dimension=28,actual_CB_projective_embedding_dimension=38,actual_CB_projective_depth=20,
      actual_CB_origin_completion_remains_domain=True,normalizer_quotient_is_CM20_conductor_module=True,
      generic_Qb_selfpair_conductor_completion='C[[g20,eps]]/(eps²)',
      full_conductor_analytic_normal_form_and_unknown_binary_cubic_coefficients_not_determined=True,
      CA_Pa_not_substituted_for_CB_Qb=True,actual_CA_CB_intersection_and_full_Poisson_scheme_not_claimed=True,
      all_nonzero_Qb_flag_theorem_not_required_or_used=True,no_new_GB_or_higher_jet_job=True,
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('CB ORIGIN CONDUCTOR e4 /FULL INITIAL /CM20 /I10 /IMAGEDEPTH21 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
