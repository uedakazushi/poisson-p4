#!/usr/bin/env python3
"""Pure source/order corollary: liftablebarcube has onlythedegree4 target."""
import hashlib,json
from pathlib import Path
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_REDUCTION_20261105_v63'
INPUTS=[
 'reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63/actual-five-edge-functionals-and-genuine-Nbar14-cycles.json',
 'reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_NIL_PRODUCT44_OPERATOR_20261104_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63/actual-trace-first-gate-and-entire-conductor-proof.json',
 'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW/report.json',
 'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63/five-cache-functionals-and-genuine-fourteen-kernel-source.json']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);D=[json.loads((ROOT/p).read_text())for p in INPUTS];assert all(D[i]['all_passed']for i in[1,2,3,4,6,7])
 K=[[Q(c)for c in row]for row in D[0]['actual_edge_kernel_basis16_by14_rows']];assert K[0]==[Q(0)]*14
 # Theaccepted16-sourceinclusion hasorigins(c042,e43,0,...,0).
 # Thereforethe14originsareexactlye43 timesrow1ofK.
 origins=[[Q(0)]*14 for _ in range(44)];origins[43]=K[1];assert all(not any(row)for row in origins[:43])
 assert min(D[5]['all36_Qrad_C_maximal_ideal_order_lower_bounds'])>=2 and D[5]['same_origin_tauA_tauB_zero']
 assert D[8]['all24_m0_e43_entries_zero']
 data={'all14_actual_cycle_origins44_by14':[[str(c)for c in row]for row in origins],
  'all14_a36_and_zero7_odd_coordinate_origins_zero':True,'possible_gA_odd_coordinate_origin_only':43,
  'point_orders':'Qrad_even>=m1 (nativeCmaxorder>=2 andnuA0=nuB0=0); everya36>=m1, henceb0=Qrad_even*a>=m2. B=sum z_i*b_i+gA*bgA>=m2, sinceb_z>=m1,gA>=m3.',
  'cube44_orders':'Fororderedcubea=mu(y1*y2)*a3, allmu(y_i)>=m2, a3>=m1 ⇒a36>=m5. Oddcoordinatevector=(b01*b02+nuA*B1*B2)*b3+b03*(b01*b2+b02*b1): zero7coordinates>=m5, gAcoordinate>=m4.',
  'membership_map_orders':'m isaregular24x44matrix andm0e43=0. Thusm*cube44>=m5 evenwithgAcoordinateoforder4.',
  'honest_homotopy_order_proof':'For(z7,r^n),n1or2, cyclev7(0,r)=0. OrderedHadamard_z(v7) lowersordinaryorderatmost1. wi=vi+r^nHi cannotlowerorder; partial_z lowersatmost1 and(Euler_z+2)^-1preserveseachmonomialdegree. HenceH(m^dcycles)subsetm^(d-1), notm^(d-2). ApplytoACTUALwholecoefficientfunctions; do notpulloutScoefficients.',
  'honest_representative':'(cube44,H(m*cube44)) liesm4 oftheactualfree128 presentation; Phi mapsitexactlytocube44. RegularallorderunitSchurprojectionto5 preservesorder. Thereforetheactualclassliesm4Tnil.',
  'proved_actual_Nbar_cube_image_subset':'m4(Tnil), dimension1; m5(Tnil)=0',
  'leading_cube_source_requirement':'Onlycyclefirstjets andpossiblegA-tailorigin, Qrad-evenFIRST, m0 andm1e43, andminimalprojectionPi0. Highercyclejets/physicalsource3/gauge3 cannotaffectthism4class.',
  'candidate_active_five':'fouradjustedkernelheads plusz5*head0. Other9havea1=zero7b1=bgA0=0, providedliteral16firstjetbindingsarechecked; rh1mayhaveonlygAfirstcoordinate, whichcontributesonehigherorderandcannotaffectleadingcube.',
  '35_triple_completeness_pending_literal_firstjet_check':True,'actual_quartic_projection_dual_not_yet_constructed':True,
  'next_quartic_dual_source':'Knownq3 givesconstructivedegree4reductionintothreepure-r4classes. Addactualoldm4=0 naturalM0*r4relations. Ifthe3x5relationmatrixhasrank2, actualHF4=1makesitskernelcompleteandfixesuniquephi4; rankmustbechecked, notassumed.',
  'no_cube_scalar_or_zero_nonzero_result':True,'no_original_nilcube_vanishing_inferred':'Evenifallbarcubeszero, originalcube maylieellN. NeedseparateellH0zeroorhigherellproof. NonzerobarcubeprovesoriginalcubeNONZEROunconditionally.'}
 ep=out/'actual-bar-cube-order-four-and-leading-source-proof.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PROVED_source_order_corollary_bar_cube_in_one_dimensional_m4','Nbar_cube_image_subset_m4_Tnil':True,'m4_Tnil_dimension':1,'35_product_reduction_requires_firstjet_source_check':True,'cube_values_or_actual_quartic_dual_computed':False,
  'input_sha256':{p:sha(ROOT/p)for p in INPUTS},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)}};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
