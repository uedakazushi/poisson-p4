#!/usr/bin/env python3
"""Source-bound exact Tnil Loewy theorem and sixteen filtered lifts."""
import hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63'
PATHS={
 'value':'reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63/report.json',
 'values':'reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_VALUES_20261105_v63/actual-head4-z5-map-and-F2-kernel-coefficient.json',
 'peer':'reports/AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_20261105_v63_REVIEW/report.json',
 'filtration':'reports/AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63/report.json',
 'filtrationpeer':'reports/AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63_REVIEW/report.json',
 'oldlength':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/report.json',
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'quad':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_REMAINING_20261105_v63_REVIEW/report.json',
 'rreduction':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_20261105_v63_REVIEW/report.json',
 'cubic':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json',
 'right':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63/report.json',
 'rightpeer':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_20261105_v63_REVIEW/report.json',
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'Q':'reports/AMBIENT_B18_TOLD_TNIL_ZERO3_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def main():
 t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=True);assert not(out/'report.json').exists();D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 for k in ['value','peer','filtration','filtrationpeer','oldlength','quad','rreduction','cubic','right','rightpeer','Q']:assert D[k]['all_passed'],k
 a,b=map(Q,D['values']['actual_w_equals_a_u_plus_b_q2_coefficients']);c=Q(11767,3652);test=1+c*a;assert test==Q(D['values']['actual_qC_coefficient_of_r_kernel_vector']) and test
 assert D['right']['actual_right106_rank']==0
 free=[i for i in range(40)if i not in D['first']['rank_minor_rows']];assert free==[5,7,15,31,39]
 assert D['natural']['old_minimal_generator_free56_rows']==[22,23,29,30,55]
 basis=[]
 for j in range(5):basis.append({'degree':0,'type':'new_minimal_unit','target':j})
 for k in free:basis.append({'degree':1,'type':'new_minimal_monomial','target':k//8,'monomial':[k%8]})
 basis +=[
  {'degree':2,'type':'natural_column_times_monomial','old_minimal_column':3,'old_free56':30,'monomial':[7],'name':'i(v), v=e3*r'},
  {'degree':2,'type':'natural_column_times_monomial','old_minimal_column':4,'old_free56':55,'monomial':[7,7],'name':'i(q2), q2=e4*r^2'},
  {'degree':2,'type':'new_minimal_monomial','target':0,'monomial':[7,7],'name':'newfree42*r^2'},
  {'degree':3,'type':'natural_column_times_monomial','old_minimal_column':4,'old_free56':55,'monomial':[7,7,7],'name':'i(g), g=e4*r^3'},
  {'degree':3,'type':'new_minimal_monomial','target':0,'monomial':[7,7,7],'name':'newfree42*r^3'},
  {'degree':4,'type':'sum_of_natural_column_monomials','terms':[[4,[5,7],Q(1)],[3,[7,7],-a],[4,[7,7,7],-b]],'name':'d=i(qC-a*e3*r^2-b*g), qC=e4*z5*r'}]
 assert len(basis)==16
 hf=[5,5,3,2,1];dims=[sum(hf[j:])for j in range(6)];assert dims==[16,11,6,3,1,0]
 data={'exact_actual_HF':hf,'exact_actual_length':16,'exact_actual_Loewy_length':5,'dimensions_m_power_Tnil':dims,
  'same_S_variables':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA'],'rA_is_NOT_DVR_ell':True,
  'rA_DVR_module_m_squared_Jordan_partition':[3,2,1],'rA_module_statement':'m^2 Tnil is isomorphic to k[[rA]]/(rA^3) plus k[[rA]]/(rA^2) plus k; no R0-module/algebra/product splitting is asserted',
  'complete_actual_cubic917_rank':1,'actual_right106_rank':0,'actual_constant811_rank':1,'no_specific_nonzero_constant811_column_evaluated':True,
  'actual_natural_map_reference':{'path':PATHS['natural'],'matrix':'matrices.actual_natural_map_minimal5_by5','frame':'same source-bound oldfree56[22,23,29,30,55] to actual newfree42,43+three keeper wedges'},
  'explicit_sixteen_filtered_basis_lifts':basis,'actual_a':a,'actual_b':b,'actual_c':c,'nonzero_qC_coefficient':test,
  'proof_exact_HF':'t=w-a*v-b*q2 belongs F3 by all three actual quadratic image identities. r*t=(1+c*a)qC-(a*alpha+b)g is nonzero because qC,g are independent in oldm^2 and1+c*a!=0. Injectivity and r-reduction imply m4Tnil!=0. Accepted m5=0,length16,and first553 leave only5533 or55321, hence55321.',
  'proof_filtered_basis':'degree0/1 are the accepted five minimal units and five fixed first NF monomials. In gr2, i(v),i(q2) are independent in the image F2 plane, while target0*r2 is detected by the literal grQ2 quotient dual. The exact old element v+c*w-alpha*q2 is killed by r; its gr2 image has nonzero i(v) coefficient1+c*a. The one-dimensional r-kernel therefore lies in that plane and avoids i(q2). Thus r*i(q2)=i(g) and target0*r3 are independent and span gr3. The displayed d spans m4 because its qC coefficient is nonzero. No alpha value is needed.',
  'old_alpha_not_evaluated':True,'no_full_actual_action_matrices':True,'not_original_nilideal_Loewy_or_cube_or_socle':True,
  'next_nil_source_gate':'Nil(D0)/ell Nil(D0) is the14-dimensional kernel of the genuine connecting map Tnil16 to Cprime[ell]2; this basis alone does not identify that map or construct actual nil generators.'}
 p=out/'exact-Loewy16-and-explicit-filtered-basis.json';p.write_text(json.dumps(enc(data),separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'status':'PROVED_source_bound_exact_Tnil_HF_and_filtered16_basis','exact_actual_HF':hf,'exact_actual_length':16,'Loewy_length':5,'dimensions_m_power_Tnil':dims,'rA_m_squared_Jordan_partition':[3,2,1],'complete917_rank':1,'constant811_rank':1,'explicit_filtered_basis_count':16,'old_alpha_or_new_source_values_computed':False,'original_nilcube_or_generators_decided':False,'input_sha256':{p:sha(ROOT/p)for p in PATHS.values()},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
