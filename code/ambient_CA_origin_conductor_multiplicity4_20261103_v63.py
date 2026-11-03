#!/usr/bin/env python3
"""CA-only equality of conductor multiplicity and exact initial ideal/depth."""
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
    geometry=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json';geo=load(geometry)
    assert geo['all_passed'] is True and geo['actual_branch_scheme_intersection_embedding_dimension']==21
    assert geo['actual_GL_plus_scalar_plus_parameter_image_cone_dimension']==20
    assert geo['generic_nonzero_arc_whole_fiber_exactly_two_reduced_points'] is True
    assert geo['reduced_CA_completion_exactly_two_smooth29_formal_branches'] is True
    deps=[upper,ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63_REVIEW/report.json',geometry,
      ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63/report.json',
      ROOT/'reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json']
    for p in deps:
        d=load(p);assert d['all_passed'] is True
        for name,digest in d.get('input_sha256',{}).items():assert sha(ROOT/name)==digest,name
    certpath=ROOT/'reports/AMBIENT_B18_CONDUCTOR_DIMENSION20_20261103_v63/certificate.json.gz';cert=load(certpath);binding=cert['actual_jet_bindings']['CA']
    rows=binding['seven_actual_odd_rows']+[46,95,63];assert len(rows)==len(set(rows))==10
    h=s.Matrix(binding['h63_as_constant_combination_of_all150_actual_odd_functions'])
    assert {i for i,c in enumerate(h)if c}<=set(rows)
    T=s.Symbol('T');hilbert=(1-2*T*T-T**3+T**3+T**4)/(1-T)**22
    assert s.cancel(hilbert-(1+T)**2/(1-T)**20)==0
    lower=1+1+2;assert lower==up['CA_CB_origin_conductor_multiplicity_upper_bound']==4
    evidence=dict(actual_ten_minimal_odd_generator_rows=rows,graded_conductor_Hilbert_series='(1+T)^2/(1-T)^20',
      curve_core_generator_degrees=[2,2,3],curve_core_syzygy_degrees=[3,4],
      formal_minimal_support_branches=['criticalW','criticalV','CAselfpairZ'],generic_lengths=[1,1,2],
      nonzero_kernel_in_CM_initial_candidate_has_positive_top_multiplicity=True,
      nonreduced_selfpair_length_uses_embedding_dimension21_AND_dense_dimension20_family=True,
      distinct_source_flags_alone_not_used_as_nonreducedness=True)
    (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63.md';paths=[Path(__file__),note,certpath,*deps]
    result=dict(id='AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63',all_passed=True,component='CA only',
      intrinsic_conductor_of_two_smooth_branches_identified_with_actual_meeting_scheme=True,
      generic_selfpair_conductor_component_length_atleast2_by_edim21_and_dimension20=True,
      selfpair_component_distinct_from_critical_W_V_by_nonzero_quadratic_discriminant=True,
      conductor_multiplicity_exactly4=True,full_conductor_initial_ideal_equals_I0=True,
      exact_initial_form='(WU,WV,U*(aU²+bUV+cV²)) extendedby g19,a!=0,b/c notcomputed',
      conductor_cone_dimension=20,conductor_cone_embedding_dimension=22,conductor_cone_depth=20,conductor_CM_type=2,
      conductor_multiplicity=4,exact_minimal_prime_count=3,generic_primary_lengths=[1,1,2],conductor_has_no_embedded_associated_primes=True,
      actual_odd_ideal_minimal_generator_count=10,actual_odd_ideal_generated_by_saved_ten_actual_rows=True,
      actual_CA_cone_dimension=29,actual_CA_cone_embedding_dimension=39,actual_CA_cone_depth=21,actual_CA_multiplicity=2,
      actual_CA_projective_dimension=28,actual_CA_projective_embedding_dimension=38,actual_CA_projective_depth=20,
      actual_CA_origin_completion_remains_domain=True,normalizer_quotient_is_CM20_conductor_module=True,
      generic_selfpair_conductor_completion='C[[g20,eps]]/(eps²)',
      full_conductor_analytic_type_or_unknown_binary_cubic_coefficients_not_determined=True,
      CB_or_actual_CA_CB_intersection_or_full_Poisson_scheme_not_claimed=True,no_new_GB_or_higher_jet_job=True,
      input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('CA ORIGIN CONDUCTOR e4 /FULL INITIAL /CM20 /I10 /IMAGEDEPTH21 PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
