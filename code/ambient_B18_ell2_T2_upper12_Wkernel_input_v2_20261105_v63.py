#!/usr/bin/env python3
"""SOURCE-only primitive rows/Hadamard and fixed12x8 kernel job preparation."""
import argparse, hashlib, json, math, signal, time
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_UPPER12_WKERNEL_INPUT_GATE_V2_20261105_v63'
AFF='reports/AMBIENT_B18_ELL2_T2_UPPER12_AFFINE_SOURCE_20261105_v63/remaining12-affine24-second-symbols-factor12-and-cache-ancestor-source.json'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,sort_keys=True)+'\n')
def count(x):
 if isinstance(x,bool) or x is None:return 0
 if isinstance(x,(int,float,F)):return 1
 if isinstance(x,str):
  try:F(x);return 1
  except ValueError:return 0
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return 0
def main(out):
 start=time.monotonic();signal.alarm(55);out.mkdir(parents=True,exist_ok=False);pins={}
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 a=read(AFF);assert a['SOURCE_only_no_rank_adjoint_actual_pairing_or_scalarfunction_values']
 cases=a['cases12'];assert cases==[579,580,593,594,597,598,605,606,611,613,615,616]
 W=[];P=[];Ds=[];gs=[];scales=[];norms=[]
 for j,rec in a['corrected_OLD4_plane_a1_from_full563609_carry']:
  assert j==cases[len(W)];row=[F(0)]*8
  for slot,q in rec:
   p,v=divmod(slot,9);assert p<4 and v in (7,8);row[2*p+v-7]=F(q)
  D=math.lcm(*(q.denominator for q in row));ints=[int(q*D) for q in row];g=math.gcd(*ints) or 1
  prim=[q//g for q in ints];s=F(D,g);assert all(F(x)==s*q for x,q in zip(prim,row))
  W.append(row);P.append(prim);Ds.append(D);gs.append(g);scales.append(s);norms.append(max(1,sum(map(abs,prim))))
 L=math.lcm(*(q.denominator for q in scales));ts=[int(L*q) for q in scales]
 DW=math.lcm(*(q.denominator for row in W for q in row));OW=[[int(DW*q) for q in row] for row in W]
 B=math.prod(sorted(norms,reverse=True)[:8]);T=max(map(abs,ts));O=max(1,max(abs(q) for row in OW for q in row))
 # No pivot, Gaussian elimination, determinant, rank or response calculation occurs here.
 kernel=B*T;factor=8*B*B;originalcheck=12*kernel*O
 bigint=max(B,kernel,abs(L),DW,max(Ds),max(gs),max(norms),T,O)
 temporary=max(2*B*B,factor,originalcheck)
 bigint=max(bigint,temporary) # This worker commits each guarded binary result, including a temporary product.
 assert bigint.bit_length()<=384 and temporary.bit_length()<=1153
 cubicp='reports/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63/literal-existing-cubic-source-input.json'
 cubic=read(cubicp);assert cubic['ell_substitution']=='sA=r+ell/5, ell²=0'
 for p in [
  'reports/AMBIENT_B18_ELL2_RAW45_TRUNCATION_PREFLIGHT_20261105_v63_REVIEW/report.json',
  'reports/AMBIENT_B18_UPPER20_INDEPENDENT_POSTCHECK_V2_20261105_v63/report.json',
  'reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63_REVIEW/report.json',
  'reports/AMBIENT_B18_NATIVE_OPAQUE12_SUPPORT_SOURCE_20261105_v63/report.json']:
  r=read(p);assert r['all_passed']
 data={'SOURCE_only_rank_and_values_UNRUN':True,'V1_prior_scope_clarification':'All guarded binary results are committed by the prepared worker. Thus include positive factor/dot partial-sum bound in committed prior; V1 79 is only input/minor/kernel bound. V2 committed129/temp129, unchanged hard384/1153.','cases12':cases,'W_columns_pivot4_then_r_ell':[[p,v] for p in range(4) for v in [7,8]],'W12_by8_original_QQ':[[str(q) for q in row] for row in W],'P12_by8_primitive_integer':P,'row_denominators_D':Ds,'row_contents_g':gs,'row_scales_s_P_equals_sW':[str(q) for q in scales],'scale_common_denominator_L':L,'scale_integer_t_Ls':ts,'original_W_common_denominator':DW,'original_W_integer12_by8':OW,'row_L1_norms_rounded_up1':norms,'Hadamard_all_minor_abs_upper_B':B,'kernel_cleared_coefficient_abs_upper':kernel,'all_binary_temporary_abs_upper':temporary,'prearith_committed_bits':bigint.bit_length(),'prearith_temporary_bits':temporary.bit_length(),'Cramer_factor_partial_sum_abs_upper':factor,'originalW_kernel_dot_partial_sum_abs_upper':originalcheck,'guaranteed_kernel_dimension_at_least4_NOT_rankvalue':True,'alternative_W9_cases_excluding_trueconductor_third':[j for j in cases if j not in [579,613,616]],'alternative_W9_guaranteed_kernel_at_least1_NOT_rankvalue':True,'full24_opaque_cancellation_rule':'lambda W=0 kills all4pivots*rr/rell*bothchannels exactly by common Gamma; no true SECOND assigned0.','source_combination_recipes':a['source563609_full_correction_for12'],'known_offset_complete_formula':'[Gamma_t*(Xc3-Xp2*a1)]_(r2ell)+D_t*(Xc2-Xp1*a1), fullPsi/nativeL2S0/FIRST/27/99/G3 carry. Uncomputed.','true_conductor_third_required_cases':a['known_TRUE_conductor_THIRD_required_cases'],'true_conductor_cubic_SOURCE_REFERENCE_ONLY':{'path':cubicp,'gA_field':'gA_plane3','gB_field':'gB_own_plane3','linear_CB_chart_field':'YB1_bottom_plane2_by2','rule':'own conductors m3, pointed maps: third uses only ownleadingcubic and knownlinear chart; sA=r+ell/5. No pairing performed.'},'native_actual_kernel_only_source_reference':'reports/AMBIENT_B18_NATIVE_OPAQUE12_SUPPORT_SOURCE_20261105_v63/report.json','rank_or_kernel_execution_requires_ROOT_GO_and_fresh_SOURCE_code_peer':True,'actual_upper_kernel_pairings_require_another_whole_SOURCE_GO':True,'hard_caps':{'committed_bits':384,'temporary_bits':1153,'numeric':12000,'RAW_bytes':33554432,'seconds':60,'memory_bytes':2000000000},'lifetime_ledger_numeric':{'double_parsed_input_max':2400,'rank_matrix96_current_det64_selected64_Cramer64':288,'all_factors96_kernel144_originalproof96':336,'report_metadata_and_pending':1040,'scientific_state_closure_max':2000,'codec_extra_two_copies_of_state':4000,'guard_metadata':1000,'whole_max':9000}}
 assert count(data)<1200
 ep=out/'primitive12x8-Hadamard-originalsource-kernel-input.json';dump(ep,data)
 report={'id':STEM,'all_passed':True,'SOURCE_only_no_rank_determinant_or_values':True,'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'payload_numeric':count(data),'whole_numeric_upper':9000,'private_numeric_cap':12000,'RAW_source_read_upper':2*sum((ROOT/p).stat().st_size for p in pins),'future_RAW_upper':2*sum((ROOT/p).stat().st_size for p in pins)+2000000,'prearith_committed_bits':bigint.bit_length(),'prearith_temporary_bits':temporary.bit_length(),'guaranteed_kernel_dimension_at_least':4,'all_actual_ranks_and_upper_responses_UNRUN':True,'seconds':time.monotonic()-start}
 assert report['future_RAW_upper']<33554432;dump(out/'report.json',report);print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
