#!/usr/bin/env python3
"""Pure source corollary: literal gA/gB belong m_R0^3 globally."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_CONDUCTORS_GLOBAL_M3_SOURCE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(exist_ok=False,parents=True);pins={};raw=0
 def read(p):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  return json.loads((ROOT/p).read_text())
 def fields(p,keys):
  nonlocal raw
  pins[p]=H.sha(ROOT/p);raw+=2*(ROOT/p).stat().st_size;assert raw<32*1024*1024
  text=(ROOT/p).read_text();dec=json.JSONDecoder();v={}
  for k in keys:
   pos=text.index(':',text.index('"'+k+'"'))+1
   while text[pos].isspace():pos+=1
   v[k]=dec.raw_decode(text,pos)[0]
  return v
 order='reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63/'
 assert read(order+'report.json')['all_passed']
 o=read(order+'case508-original595-source-order-and-plane-cut-proof.json');assert o['gA_plane_exact_node']==505938 and o['gB_plane_exact_node']==299650 and o['old_f63_distinct_node']==28539 and o['both_zero7_conductor_functions_have_zero_all6_source2_coefficients']
 for p in ['reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json','reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json','reports/AMBIENT_B18_ELL2_V_CONORMAL_HADAMARD_SOURCE_PROFILE_20261105_v63_REVIEW/report.json','reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW/report.json']:
  assert read(p)['all_passed']
 p='reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json';v=fields(p,['source2_monomials','both57_source2'])
 cut=[k for k,m in enumerate(v['source2_monomials']) if set(m)<={7,8}];assert len(cut)==6
 fam=[]
 for name,f in v['both57_source2'].items():
  assert name in ['CA','CB'] and len(f['equations57'])==57 and F(f['determinant'])!=0 and f['all57_equations_all55_coefficients_verified']
  root=f['f63'];assert root==(7790 if name=='CA' else 28539)
  jet=f['pure_function2jets'][str(root)];assert all(F(jet[k])==0 for k in cut)
  fam.append({'family':name,'original_f63':root,'own_F9_nodes':f['F9_nodes'],'unit57_determinant':f['determinant'],'zero7_cut_monomial_indices6':cut,'all6_known_exact_source_coefficients0':['0']*6,'allorder_unit57_implicit_function_defined':True})
 assert len(fam)==2;family={q['family']:q for q in fam};del v
 # The shared physical R0 chart and BOTH own-F9 inverse charts are pointed.
 z=fields('reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',['families'])['families'];center=[]
 for name in ['CA','CB']:
  q=z[name];assert q['own_F9_coordinate_order']==family[name]['own_F9_nodes'] and q['all9_inverse_composition_identities_through2_checked']
  assert F(q['unit9_determinant'])!=0 and all(F(p[0])==0 for p in q['physical_R0_readback_jets_in_own_F9']) and all(F(p[0])==0 for p in q['inverse_F9_in_common_R0_through2'])
  center.append({'family':name,'actual_own_F9_nodes9':q['own_F9_coordinate_order'],'actual_common_physical_R0_rows9':q['physical_R0_rows9'],'all9_readbacks_and_inverse_map_constants0':True,'unit9_determinant':str(q['unit9_determinant']),'local_homomorphism_is_allorder_unique_formal_inverse':True})
 del z
 g=fields('reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json',['CB_coordinate_Y9_nodes','CB_f63_node','CB_g_restriction_node','all16_pointed_unit57_inverse_instances'])
 assert g['CB_coordinate_Y9_nodes']==family['CB']['own_F9_nodes'] and g['CB_f63_node']==28539 and g['CB_g_restriction_node']==299650
 cbprefix=[q for q in g['all16_pointed_unit57_inverse_instances'] if q['retained_Y_indices']==[7,8]];assert len(cbprefix)==1
 cbprefix=cbprefix[0];assert len(cbprefix['all57_pointed_equations'])==57 and len(cbprefix['bindings'])==57 and F(cbprefix['unit_Jacobian'])!=0
 ev={'source_coefficient_ring':'R0=QQ[[physical_x9]]=QQ[[YA0..YA8]] by the accepted pointed unit9 coordinatechange','literal_conductor_ids':[505938,299650],'original_unrestricted_f63_ids_NOT_assigned_order3':[7790,28539],'own_unit57_zero7_function_records':fam,'both_pointed_local_R0_to_own_F9_pullbacks':center,'actual_CB_prefix_retained_indices':[7,8],'actual_CB_prefix57_source_record':cbprefix,'CA_prefix0_allorder_binding_adopted_from_exact_independent273577':'505938=prefix0[7790]; set first7 ownCA coordinates0 before the actual pointed pullback','CB_prefix_allorder_binding_adopted_from_exact_fullrelative_peer':'299650=at(28539,(7,8)); set first7 ownCB coordinates0 BEFORE the actual pointed pullback','formal_proof_steps':['For each family the unique unit57 implicit solution gives a formal function f63(Y0,...,Y8).','The zero7 restriction psi(u,v)=f63(0,...,0,u,v) has all six exact coefficients throughdegree2 equal0, hence liesin(u,v)^3.','The accepted literal prefix constructor identifies gA andgB with these restricted functions pulledback through the corresponding actual own-F9 chart. This is not restrictionofunrestrictedf63 AFTERpullback.','Each actual own-F9 coordinate map is local: allnine constants0 in the same commonphysicalR0 chart, withnonzero unit9Jacobian.','A local homomorphism sends(u,v)^3 into m_R0^3. Therefore gA505938 andgB299650 belongm_R0^3 globally, notonlyontheCA r/ell plane.','Allordinarysecondcoefficients, including zi*r andzi*ell, vanish; reparametrization sA=rA+ell/5 and ell²basechange preserve thiscontainment.'],'global_m3_membership_source_corollary':True,'all_mixed_second_coefficients_vanish_by_order_not_by_zero_fill':True,'no_new_physical_coefficient_quotient_response_gradient_or_upper_value_computed':True,'scope':'Only thetwo literal zero7 conductorfunctions; NOT unrestrictedf63, otherg/h quotients, opaque12 frontier, oranyupper/HF/length claim.'}
 ep=out/'literal-zero7-unit57-and-pointed-pullback-global-m3-proof.json';H.dump(ep,ev)
 owned=count([o,fam,center,g,cbprefix,ev,pins])+15000;assert owned<60000 and raw<32*1024*1024
 report={'id':STEM,'all_passed':True,'pure_source_corollary_from_accepted_exact_unitIFT_and_saved_source2':True,'global_mR0_cubed_for_literal_gA_andgB':True,'allordinarymixed_seconds_zero_by_formal_order':True,'not_unrestricted_f63_or_other_quotients':True,'no_new_function_response_or_sourcejet_values':True,'source_owned_bound':owned,'whole_SOURCE_RAW_upper':raw,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','source_owned_bound','whole_SOURCE_RAW_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
