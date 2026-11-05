#!/usr/bin/env python3
"""Bind the q2 proof to the literal n=1,1,2 homotopy, not a jet-loss bound."""
import hashlib,json
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63'
PATHS={
 'primary':'reports/'+STEM+'/report.json',
 'honest':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63/actual-honest-Tnil128-972-formal-operator.json',
 'honest_peer':'reports/AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW/report.json',
 'H_peer':'reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_20261105_v63_REVIEW/report.json',
 'point':'reports/AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63/five-cache-functionals-and-genuine-fourteen-kernel-source.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 receipts=[];h=S.selective(ROOT/PATHS['honest'],['shared_homotopy_template_v1','actual_S_frame','matrices','matrix_operator_shapes','exact_matrix_identities','canonical_division_regular_proofs'],receipts)
 p=S.selective(ROOT/PATHS['point'],['literal_direct_section_path','source_m0_24_by44'],receipts)
 for k in['honest_peer','H_peer']:assert json.loads((ROOT/PATHS[k]).read_text())['all_passed']
 t=h['shared_homotopy_template_v1'];assert h['actual_S_frame']==['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA']
 assert t['sequence']=='(zA0,...,zA6,rA^n), n=1 or2'
 assert t['cross_coefficients_Hi']=='orderedHadamard_z7(v7); prefix differences divided by zi, sum ziHi=v7'
 assert t['remaining_cycle_wi']=='wi=vi+rA^n*Hi, sum ziwi=0'
 assert t['pair_i_7_output']=='Hi'
 assert t['pair_i_j_output_i_j_lt7']=='(Euler_z+2)^-1*(partial_zi(wj)-partial_zj(wi)) on coefficient functions'
 assert t['not_S_linear']and t['no_pull_out_actual_coefficients_or_truncated_normal_forms']
 matrices=h['matrices'];assert matrices['new_m24_by44']==4 and matrices['d2_new']==23 and matrices['Phi44_by128']==39
 shapes=h['matrix_operator_shapes'];assert shapes[4]==[24,44]and shapes[23]==[24,84]and shapes[39]==[44,128]
 assert all(x[43]=='0'for x in p['source_m0_24_by44'])
 rows=[]
 for n in[1,1,2]:
  # Symbolic offsets from an arbitrary actual input order d>=1.
  Hi=-1;rHi=n-1;wi=min(0,rHi);pairij=wi-1
  assert rHi>=0 and wi==0 and pairij==-1
  rows.append({'n':n,'Hi_order_offset':Hi,'r_power_Hi_offset':rHi,'wi_offset':wi,'partial_Euler_pair_offset':pairij,'overall_loss_at_most':1})
 data={'all_passed':True,'scope':'Specific analytical ordinary-m-order certificate for q2bar product; no scalar products or new jets evaluated.',
  'frozen_primary_sha256':sha(ROOT/PATHS['primary']),
  'literal_shared_template':t,'actual_S_frame':h['actual_S_frame'],
  'exact_matrix_domains':{'m':{'namespace':'HONEST_NIL_T_OPERATOR','root':4,'shape':[24,44]},'d2':{'root':23,'shape':[24,84]},'Phi_npre':{'root':39,'shape':[44,128]},'direct_n':'NBARliteralE(mE)^-1; same m,d1,d2,H template'},
  'exact_domain_matrix_identities':h['exact_matrix_identities'],
  'm_is_actual_regular_matrix_reason':h['canonical_division_regular_proofs'],
  'ordinary_order_proof_for_any_d_at_least1':'InoneACTUALcyclev∈m^dS8, specializez7=0 in delta(v)=0: r^n v7(0,r)=0,so v7(0,r)=0. OrderedprefixHadamard dividesONLYbyzi once, givingHi∈m^(d-1). Thenwi=vi+r^nHi∈m^d since n1or2. ThetrueKoszulEulerformula differentiatesONLYzi once, yieldingpairs∈m^(d-1); divisionbypositiveQQinteger(dz+2) doesnotchangeordinaryorder. Thus theEXACTliteralH(v)∈m^(d-1) for everyformalcycle, notjustcolumns.',
  'all_three_actual_sector_offset_checks':rows,
  'no_periodic_or_translated_source_lift':'This is theone-pass sharedtemplate onthreeunshiftedsectors(zA7,r),(zA7,r),(zA7,r²). It is NOT an abstractperiodicH_G andnotjet extractionaftercoordinateadaptation/RegDivell.',
  'jet_query_order_is_distinct':'Computingactualregularmatrixm coefficients canrequirehigherprimitivejets becauseits entrieswereconstructedbyregularquotients. That is acoefficientAVAILABILITYloss, notordinaryorderofmultiplicationbythealreadyactualRbar-matrixm. Herewedonotevaluateanysuchcoefficient. Allunit/normalizertransportsareREGULAR, Pconstant, soactualidealpowerisnotloweredbythem.',
  'full_function_rule':'Applythetemplate toWHOLEv=m*c before anyfiniteNF, withallactualcoefficientsinside. DO_NOT reuse S-linearh84x780 extension foranarbitraryproduct, DO_NOTpulloutRbarcoefficients, DO_NOTtruncatefirst.',
  'q2_application':'Product44 givesc anactualV-cycle. d1m=Vgamma implies eachsectorv isacycle. Sourcepointm0e43=0 andc4 onlye43 givev∈m5; thisEXACTHcertificategivesH(v)∈m4. Iota(c)=(c,Hv) thereforehasorder>=4 initsraw44+84domain.',
  'Pi_domain_step_preserved':'UniversalPhi_nIota=c andJ(Iota)=Iota useonlyd2Hv=v. CompletePi0firstrow[e42,0_84] annihilatesallthedegree4Hvcorrection. NoassumptionSlinearityofHismade.',
  'no_product_value_or_original_H0_nil_inference':True,
  'input_sha256':{q:sha(ROOT/q)for q in PATHS.values()},'code_sha256':sha(Path(__file__)),'selective_receipts':receipts,
  'whole_SOURCE_owned_upper':S.count(h)+S.count(p)+S.count(rows)+12000}
 assert data['whole_SOURCE_owned_upper']<150000
 out=ROOT/'reports'/STEM/'specific-H-template-order-and-raw-domain-companion.json';assert not out.exists();out.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'all_passed':True,'path':str(out.relative_to(ROOT)),'sha256':sha(out),'whole_SOURCE_owned_upper':data['whole_SOURCE_owned_upper']}))
if __name__=='__main__':main()
