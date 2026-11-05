#!/usr/bin/env python3
"""Metadata/source closure for one future shared605 bundle; no responses."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,(F,int)):return 1
 if isinstance(v,str):
  try:F(v);return 1
  except (ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 return 0
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 pins={};held=[];raw=0
 def read(p):
  nonlocal raw
  pp=ROOT/p;pins[p]=sha(pp);raw+=2*pp.stat().st_size;assert raw<32*1024*1024
  v=json.loads(pp.read_text());held.append(v);assert count(held)+15000<150000;return v
 def passed(p):r=read(p);assert r['all_passed'];return r
 jp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_INPUT_GATE_20261105_v63/617-direct-old-known-first-formal-Jacobian-positive-input.json'
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/combined-OLD2-D-and-four-pivot-source-factor-inventory.json'
 fp='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63/617-formal-direct-old-third-and-second-symbol-sensitivities.json'
 op='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json'
 j=read(jp);s=read(sp);jf=read(fp);of=read(op)
 passed('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63_REVIEW/report.json')
 passed('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63_REVIEW/report.json')
 passed('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FRONTIER_SOURCE_20261105_v63_REVIEW_V2/report.json')
 passed('reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63_REVIEW/report.json')
 rawp='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';physical=read(rawp);f164={i for i,v in physical['all_requested164_root_values']};assert len(f164)==164
 typed=dict(s['typed_program']);subtyped=dict(j['typed_source_program']);assert len(typed)==757 and len(subtyped)==629 and all(typed[i]==n for i,n in subtyped.items())
 assert all(n[0]!='product' or n[3] for n in typed.values())
 boundaries={i:k for i,k,v in j['literal_boundary_firsts_point_r_ell_and_membership']};assert len(boundaries)==96
 leaf67={n[1] for n in typed.values() if n[0]=='actual_C0_scalar_function'};assert len(leaf67)==67 and leaf67<=set(boundaries)
 assert len(j['literal_C0_program'])==532
 assert all(n[0] in ['c','n','a','m'] for i,n in j['literal_C0_program'])
 opaque={330164,330289,330524,330594,330184,330309,330544,330614,566400,566401,566402,566403};conductors={299650,505938}
 assert not opaque&f164 and not conductors&f164
 # The absence is from the formal SECOND linear ancestry only; the full
 #757-node program keeps every multiplicative first-cross contribution.
 linear_leaves={}
 for i,n in s['typed_program']:
  if n[0]=='actual_C0_scalar_function':v={n[1]}
  elif n[0]=='sum':v=set().union(*(linear_leaves[k] for k,a,b in n[1] if a)) if n[1] else set()
  else:v=set(linear_leaves[n[3][0]]) if len(n[3])==1 and n[1] else set()
  linear_leaves[i]=v
 Dlinear=set().union(*(linear_leaves[r] for col,roots in s['D_old_root_columns'] for row,r in roots))
 assert not Dlinear&opaque
 outer=s['outer_four_pivot_source'];assert not outer['inverse_nodes']
 assert all(n[0]=='actual_C0_scalar_function' and n[1] in opaque for i,n in outer['typed_program'])
 assert all(n[0] in ['g','a','n'] for i,n in outer['literal_C0_program'])
 assert {i for i,n in outer['literal_C0_program'] if n[0]=='g'}==set(range(566400,566404))
 assert all(n[0]!='product' for i,n in outer['typed_program'])
 recipes=s['kernel617_source_recipes'];assert len(recipes)==617
 fre=[r['free_case619'] for r in recipes];assert set(fre)==set(range(619))-{563,609}
 formal2=dict(of['all617_combined_OLD2_independent_symbol_coefficients']);formal1={r['free_case619']:r['two_formal_J0_Jr_Jell_symbol_rows'] for r in jf['all617_formal_directOLD_J0_Jr_Jell']}
 assert set(formal1)==set(formal2)==set(fre)
 affected=sorted(i for i in fre if any(fn in opaque for row in formal2[i] for pair,fn,q in row))
 assert affected==[579,580,593,594,597,598,605,606,611,613,615,616]
 good=sorted(set(fre)-set(affected));assert len(good)==605
 unknown2={};thirdphysical=set();secondphysical=set();allthird=set();allsecond=set();conductor_third_cases=set()
 for i in fre:
  for ch in range(2):
   for fn,trip in formal1[i][ch]:
    assert len(trip)==3
    for degree_component,q in enumerate(trip):
     if q=='0':continue
     if degree_component==0:allthird.add(fn)
     else:allsecond.add(fn)
     assert fn in f164|conductors,'uncached genuine honest high symbol survived DIRECTOLD J1'
     if degree_component==0 and fn in conductors:conductor_third_cases.add(i)
     if i in good and degree_component==0:assert fn in f164;thirdphysical.add(fn)
   for pair,fn,q in formal2[i][ch]:
    assert fn in f164|conductors|opaque and q!='0'
    if fn in opaque:
     assert pair in [8,9];unknown2.setdefault(fn,set()).add(pair)
    if i in good:
     assert fn not in opaque
     if fn in f164:secondphysical.add(fn)
 assert conductor_third_cases<=set(affected)
 assert set(unknown2)==opaque and all(v=={8,9} for v in unknown2.values())
 # The outer is EXACTLY linear in these12 actual regular functions, not a
 # first-order approximation. No products means no omitted first-cross term.
 # D-side opaque linear support is empty, so zero formal2 opaque weights
 # imply the entire corrected605 OLD4 plane-second outer functional is0.
 nativep='reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/all617-native-constant-sources-factor-through-fixed105-kernel.json';native=read(nativep)
 assert native['all617_literal_full48S0cN_zero_verified'] and native['all617_full128_reconstructions_verified']
 nm={r['free_case619']:r for r in native['all617_source_maps']};assert set(nm)==set(fre)
 native_indices=sorted({k for i in good for k,q in nm[i]['kernel105_coordinates']})
 assert len(native_indices)==72 and sum(len(nm[i]['kernel105_coordinates']) for i in good)==120
 gb='reports/AMBIENT_B18_ELL2_T2_GB_PLANE3_CACHE_SOURCE_20261105_v63_REVIEW/report.json';passed(gb)
 # Global conductor m3 is a separate theorem. It zeros SECOND, never THIRD.
 globalp='reports/AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_SOURCE_20261105_v63/report.json';passed(globalp)
 data={'scope':'SOURCE availability and linearity only; no function jets, scalar response, native basis response, FIRST table or605 upper class evaluated','cases605':good,'opaque_cases12':affected,'all24_opaque_second_symbols_retained':[[fn,p] for fn in sorted(opaque) for p in [8,9]],'one_shared_typed757_contains_Gamma629_exactly':True,'one_shared_literal_C0_program532':True,'all67_typed_leaves_in_shared96_boundary_catalog':True,'shared_typed_and_literal_program_source':{'OLD2_source':sp,'DIRECTOLD_J1_source':jp},'full757_multiproducts_FIRSTcross_kept_not_pruned104':True,'linear_D_opaque_second_support_empty':True,'outerOLD4_only_add_neg_literal12_function_linear_source':True,'corrected605_complete_outerOLD4_SECONDfunctional_ZERO_by_full_formal_certificate_NOT_raw_a1absence':True,'raw_basis609_OLD4plane_carry_explicitly_not_discarded':dict(s['raw619_old4_plane_a1_nonempty28coefficients'])[609],'corrected605_directOLD_third_physical164_ids':sorted(thirdphysical),'corrected605_combinedOLD2_physical164_ids':sorted(secondphysical),'all617_directOLD_third_remaining_functions':sorted(allthird),'all617_directOLD_second_remaining_functions':sorted(allsecond),'conductors299650_505938_thirds_remain_real_cache_inputs_global617':True,'conductor_third_cases_only_in_opaque12':sorted(conductor_third_cases),'gA_third_cache_recipe':'505938=prefix0[7790]; on planezA7=0 its r2ell coefficient is literal7790 raw33 selector30','gB_third_cache_recipe':'accepted negative-boundary tail: fCB28539_r2ell - SUM7(zB_rr*h_ell + zB_rell*h_r); actual zBpoint/plane1=hpoint0; no h2/h3 required','gA_gB_global_m3_SECONDzero_source':globalp,'global617_ONLY_unbound_high_boundary_coefficients_after_complete_formal_certificate':'24 opaque second symbols above; genuine honest thirds/other seconds canceled, physical164 and true conductor thirds remain bound inputs, FIRSTcross remains','future_boundary_extension_rule':'Physical164 use actual B21. Only79 genuine honest higher may receive arbitrary0 extension ONLY inside whole certified functional. This does not assert true function jets0. Conductors SECOND actual0 bym3, thirds real when supported.','B21_raw33_selectors':[0,7,6,5,4,3,2,1,8,9,22,23,20,21,14,15,12,13,24,25,30],'B21_basis':'1; all9firsts; nonzz10 ordinary quadratics; r2ell. Divisorclosed quotient/source maps already independently verified on660NF+1089products.','native105_basis_indices_needed_for605':native_indices,'native_coordinate_nonzeros_for605':120,'native_source_factorization_reference':nativep,'native_raw_auxiliary_responses_not_individually_actual':'only fixed kernel105 restrictions are genuine; do not label extra mask0 auxiliary partial-N columns as actual jets','FIRST_shared_table_recipe':'T[t,p,v]=SUM pairs containingv D[t,pair] dot Xp[first_other]; diagonal rr ONCE. Compute shared27x6x2 contractions once using actual individual nativeL0Se+LeS0. Then apply full corrected a1_j-a*a1_563-b*a1_609, keeping every knownFIRSTcross.','native_outer_table_recipe':'accepted native92 gives23x2x2 rr/rell targets withfullL2S0. Pair with every corrected planea1; OLD4outer handled by complete linear formal theorem above.','corrected_source_basis563609_upper_NOT_set0':True,'future_complete_upper_formula':'shared_DIRECTOLD + shared_native105 - shared_FIRST - accepted_native92outer; all maps on literal corrected617 sources before response','future_positive_height_and_wholeownership_and_immutableworker_and_independentpeer_and_ROOT_GO_REQUIRED':True}
 ep=out/'complete605-shared-source-linearity-and-global24-opaque-frontier.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
 owned=count(held)+count([data,linear_leaves,Dlinear,formal1,formal2,nm])+15000;assert owned<150000
 report={'id':STEM,'all_passed':True,'SOURCE_metadata_ONLY_no_signed_source_response':True,'cases605':605,'opaque_second_symbols24':24,'shared_programs757_and532':True,'native_basis_used72':True,'all_actual_response_values_UNRUN_by_this_job':True,'no_HF2_equality_total_T2_H0_or_all605zero_claim':True,'owned_numeric_upper':owned,'RAW_read_upper':raw,'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:report[k] for k in ['all_passed','cases605','owned_numeric_upper','RAW_read_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
