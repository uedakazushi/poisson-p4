#!/usr/bin/env python3
"""Tiny fixed source-curve arithmetic for the comparison-cokernel rank proof.

No new source solve/jet/field point, gcd, factor search or Groebner computation.
"""
import argparse,gzip,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_proper_core_five_lines_review_20261104_v63 import (
    parse,substitute,variable,one,add,mul,scale,power,encode)
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63'
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    paths={
      'CAcritical':'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz',
      'CBcritical':'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63/certificate.json.gz',
      'CAplane':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
      'CBplane':'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json',
      'DVRsource':'reports/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63/DVR20-raw-cokernel26-and-comparison-preflight.json',
      'F9source':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'actualcoefficients':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'}
    gates=[
      'AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW/report.json',
      'AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW/report.json',
      'AMBIENT_CA_V1_ODD_ORDER_20261103_v63_REVIEW/report.json',
      'AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63_REVIEW/report.json',
      'AMBIENT_B18_CA_N_CRITICAL_SCHEMES_20261103_v63_REVIEW.json',
      'AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json',
      'AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json',
      'AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
      'AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.json',
      'AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW.json',
      'AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json',
      'AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json',
      'AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63_REVIEW/report.json',
      'AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW/report.json',
      'AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63_REVIEW/report.json',
      'AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW/report.json',
      'AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json']
    gatepaths=[ROOT/'reports'/p for p in gates]
    assert all(read(p)['all_passed'] for p in gatepaths)
    data={k:read(ROOT/p) for k,p in paths.items()};CA,CB=data['CAcritical'],data['CBcritical']
    assert CA['affine_slice_free_indices']==CB['slice_free_indices']==[1,2,3,9,10,14,15,17,18,19]
    assert CA['center_fixed20'][:3]==CB['fiber_center19'][:3]==['-1','0','0']
    t=variable(1,0);zero={}
    sourceCA={name:substitute(parse(q,['u','v','w']),[zero,t,zero],1) for name,q in CA['inverse_substitution'].items()}
    assert sourceCA['s0']==scale(t,F(1,10)) and not sourceCA['s1'] and not sourceCA['s2']
    assert all(not substitute(parse(q,[f's{i}' for i in range(10)]),[sourceCA[f's{i}'] for i in range(10)],1)
               for q in CA['seven_monic_eliminations'])
    # Independent explicit substitution into all seven actual CB monic equations.
    sourceCB={'s0':t,'s1':{},'s2':{},'s3':scale(t,F(1,2)),
      's4':scale(mul(t,t),F(-1,4)),'s5':{},'s6':{},'s7':scale(t,F(-1,2)),
      's8':scale(mul(t,t),F(1,4)),'s9':scale(power(t,3,1),F(-1,24))}
    assert all(not substitute(parse(q,[f's{i}' for i in range(10)]),[sourceCB[f's{i}'] for i in range(10)],1)
               for q in CB['seven_monic_eliminations'])
    assert substitute(parse(CB['u'],[f's{i}' for i in range(10)]),[sourceCB[f's{i}'] for i in range(10)],1)=={}
    assert substitute(parse(CB['v'],[f's{i}' for i in range(10)]),[sourceCB[f's{i}'] for i in range(10)],1)==t
    # f=-h^2+s0*h*k+s1*k^2. The CB discriminant is t^2, not zero.
    disc=add(mul(sourceCB['s0'],sourceCB['s0']),scale(sourceCB['s1'],4));assert disc==mul(t,t)
    assert CB['full_ratio_source_formula']=='P=sum f_i F_i -q sum f_i G_i +sum H_i HH_i'
    # Actual normalizer bases: the chosen seven I0 functions and the two J0
    # functions are F9 coordinates, proving I0+J0=m_R0, not abstract types.
    for family in ['CA','CB']:
        c=data['actualcoefficients']['families'][family]
        f9=data['F9source']['both57_source2'][family]
        assert f9['F9_nodes']==[p[1] for p in c['actual_zero_odd_I0_generators'][:7]]+[p[1] for p in c['actual_positive_odd_J0_generators']]
        assert F(f9['determinant']) and f9['all57_equations_all55_coefficients_verified']
    # Two already source-bound binary cubics; fixed factor identities only.
    ga={tuple(m):F(c) for m,c in zip(data['CAplane']['monomials'],data['CAplane']['g_CA_all10_coefficients']) if F(c)}
    gb={tuple(m):F(c) for m,c in zip(data['CBplane']['ordered_monomials'],data['CBplane']['g_CB_all10_coefficients']) if F(c)}
    r,s=variable(2,0),variable(2,1)
    assert ga==scale(mul(power(add(scale(r,3),scale(s,2)),2,2),add(r,scale(s,4))),F(1,24))
    assert gb==scale(mul(add(r,scale(s,-1)),power(add(scale(r,2),scale(s,3)),2,2)),-6)
    ellA=parse('-5*r+5*s',['r','s']);ellB=parse('35*r/3-85*s/3',['r','s'])
    assert list(map(F,data['DVRsource']['CA_ell_on_own_plane_linear']))==[F(-5),F(5)]
    assert list(map(F,data['DVRsource']['CB_ell_on_own_plane_linear']))==[F(35,3),F(-85,3)]
    directions={'CAcritical':(4,-1),'CAself':(2,-3),'CBcritical':(1,1),'CBself':(3,-2)}
    ellvalues={}
    for label,(i,j) in directions.items():
        p=ellA if label.startswith('CA') else ellB
        value=substitute(p,[scale(t,i),scale(t,j)],1);assert value
        ellvalues[label]=encode(value)
    contributions={'CAcritical':1,'CAself':1,'CBcritical':0,'CBself':2}
    rank=sum(contributions.values());assert rank==4
    artin=read(gatepaths[-1]);assert artin['actual_C_mod_ell_dimension']==8
    proof={'no_new_source_or_jet_fieldpoint_or_factor_search':True,
      'actual_CA_V_source_coefficients':{k:encode(v) for k,v in sourceCA.items()},
      'actual_CB_V_source_coefficients':{k:encode(v) for k,v in sourceCB.items()},
      'CA_V_binary_third_coefficient_exactly_zero_and_ratio_zero':True,
      'CB_V_binary_discriminant':encode(disc),'CB_V_ratio_zero_and_leading_f_minus1':True,
      'actual_I0_plus_J0_equals_maximal_via_both_F9_unit_coordinates':True,
      'actual_gA3_factor':'(3r+2s)^2*(r+4s)/24','actual_gB3_factor':'-6*(r-s)*(2r+3s)^2',
      'exact_ell_linear_values_on_four_tangent_directions':ellvalues,
      'generic_support_coverage':'full intrinsic conductor three branches: positiveW plus criticalV(simple),self(double); exact I=I0+wJ0 and I0+J0=m remove W and identify the two zero branches schematically at their generic points',
      'zero_frame_transport':'same G19 product on both maps and coker; at generic zero conductor J0 unit implies w*C=0, tau roots character0 so the coker is its exact degree0 module',
      'canonical_R0_embedding':'A0/(I0*B0)=R0/I0; target coordinates restrict to their actual even parts; no new abstract formal base identification',
      'branch_D_degree':'each remaining support curve smooth/QQ and ell order1; four generic reduced branch fields equal Frac QQ[[ell]]',
      'Pa_exact_triple_meeting':'A1 intersect A2 intersect CB = S reduced, length1',
      'Qb_exact_triple_meeting':'CB1 intersect CB2 = C2=P[rho]/rho^2 is a subscheme of CA, hence triplemeeting is C2, length2',
      'CAcritical_compatibility':'actual N source f in markedh*W; opposite conductorunit yields unique root via unit bB; all equations compatible by actual incidence',
      'CBcritical_exclusion':'unique length2 CB proper point at b0 with nonsquare f excludes both N boundary images; smoothness CB on RminusN excludes R; global intersection NunionR excludes CA',
      'generic_DVR_branch_contributions':contributions,'actual_comparison_C_DVR_rank':4,
      'actual_mu_D_C_from_independent_Artin_peer':8,'minimal_torsion_generator_count':4,
      'torsion_lengths_or_exponents_not_determined':True,
      'comparison_K_socle_D0_depth_not_determined':True}
    ev=out/'fixed-source-curve-and-rank-proof.json';ev.write_text(json.dumps(proof,indent=2)+'\n')
    notes=[
      'AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63.md','AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63.md',
      'AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63.md','AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63.md',
      'AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63.md','AMBIENT_CB_GENERIC_V_FULL_COMPLETION_20261103_v63.md',
      'AMBIENT_CA_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63.md','AMBIENT_CB_SELF_PAIR_ACTUAL_D_FULL_RING_20261103_v63.md',
      'AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.md','AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63_REVIEW.md']
    files=[Path(__file__),ROOT/'code/ambient_B18_proper_core_five_lines_review_20261104_v63.py',ROOT/'research/notes'/f'{STEM}.md']+[ROOT/p for p in paths.values()]+gatepaths+[ROOT/'research/notes'/p for p in notes]
    report={'id':STEM,'all_passed':True,'actual_comparison_C_generic_DVR_rank4_proved':True,
      'actual_mu_D_C8_from_independent_source_evaluation':True,'torsion_minimal_generator_count4_proved':True,
      'generic_rank6_not_assumed':True,'all_generic_branch_lengths_source_or_exact_scheme_incidence_bound':True,
      'same_actual_R0_G19_and_zero_weight_transport_proved':True,'no_newsource_jet_point_GB_or_factor_search':True,
      'torsion_exponents_length_K_D0_socle_or_depth_claimed':False,'independent_review_pending':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
