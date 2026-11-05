#!/usr/bin/env python3
"""Pure Loewy-length consequence of the actual injective Told->Tnil map.
No residual scalar coefficient, source jet, or gauge value is evaluated.
"""
import hashlib,itertools,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TNIL_QUADRATIC_DIMENSION_BOUNDS_20261105_v63'
STEMS=['AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63',
 'AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW',
 'AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63',
 'AMBIENT_B18_T_UNISERIAL_DEEP_TAIL_20261104_v63_REVIEW',
 'AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63',
 'AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW',
 'AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63',
 'AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW',
 'AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63',
 'AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW',
 'AMBIENT_B18_HONEST_NIL_T_QUADRATIC_661_OPERATOR_20261104_v63']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 paths=[f'reports/{s}/report.json'for s in STEMS];D=[json.loads((ROOT/p).read_text())for p in paths]
 assert all(d['all_passed']for d in D)
 assert D[0]['actual_cokernel_length']==5 and D[0]['actual_cokernel_mu']==3
 assert D[2]['actual_T_total_length']=='10+s' and D[2]['all_zi_deep_actions_in_rA_squared']
 assert D[4]['actual_HF_through_degree3']==[5,3,2,1]
 assert D[6]['actual_Tnil_length_formula']=='length(T)+5'
 assert D[8]['actual_HF_Tnil_degree0_1']==[5,5]
 assert D[10]['all814_sources_complete'] and D[10]['residual_actual_operator_shape']==[6,661]
 (out/'phase.json').write_text(json.dumps({'phase':'all_actual_module_gates_pinned_before_pure_count','input_sha256':{p:sha(ROOT/p)for p in paths}},indent=2)+'\n')
 # For a graded quotient of k[z0,...,z6,r]^5, exactly one degree2
 # standard monomial permits at most one degree3 standard monomial.
 # All its degree2 divisors must be that SAME monomial, in the SAME target.
 q2=list(itertools.combinations_with_replacement(range(8),2));q3=list(itertools.combinations_with_replacement(range(8),3))
 choices=[]
 for target in range(5):
  for m in q2:
   possible=[]
   for n in q3:
    div={tuple(n[:j]+n[j+1:])for j in range(3)}
    if div=={m}:possible.append(n)
   assert possible==[(m[0],)*3]if m[0]==m[1]else possible==[]
   choices.append([target,list(m),len(possible)])
 assert len(choices)==180
 # Finite range checks only verify the integer consequences. They are not
 # source/module coefficient computations.
 inequalities=[]
 for s in range(1,8):
  old=10+s;new=15+s
  assert old-10==s and new-10-5==s
  # m^(s+6)Tnil=0 allows degrees2,...,s+5: s+4 slots.
  upper_if_h2_one=10+(s+4)
  assert upper_if_h2_one==new-1
  inequalities.append({'s':s,'length_old':old,'length_new':new,'nilpotence_power_bound':s+6,'max_length_if_HF2_is_one':upper_if_h2_one,'strict_gap':1})
 data={'same_actual_coefficient_ring':'S=QQ[[zA0,...,zA6,rA]]=R0/(ell); rA is NOT DVR ell',
  'actual_exact_sequence':'0 -> Told -> Tnil -> Q -> 0, length(Q)=5, mu(Q)=3',
  'actual_Told_HF':'(5,3,2, then s consecutive ones), 1<=s<=7',
  'actual_Tnil_initial_HF':[5,5],
  'proof_upper':'i(m^3 Told) subset m^3 Tnil; lengths s and s+5-h2 imply h2<=5',
  'Q_loewy_power_three_zero':'A finite local module of length5 and mu3 has at most two further nonzero layers; hence m^3 Q=0',
  'Tnil_loewy_power_bound':'m^3 Tnil subset i(Told), m^(s+3)Told=0; hence m^(s+6)Tnil=0',
  'proof_lower':'h2=0 gives length10. If h2=1, a monomial initial module has one degree2 standard monomial; higher standard monomials must be powers of the same variable in the same target, so each higher HF<=1. Total length<=10+(s+4)=s+14 contradicts s+15.',
  'all180_single_degree2_monomial_cases':choices,'all7_integer_length_checks':inequalities,
  'actual_HF2_interval':[2,5],'actual_complete6_by661_rank_interval':[1,4],
  'all_six_residual_quadratic_adjoints_cannot_vanish':True,
  'if_HF2_equals5':'m^3 Tnil=i(m^3 Told) by inclusion and equal finite length; all later graded layers agree with Told',
  'general_truncation_inequality':'length(Tnil/m^k Tnil)<=length(Told/m^k Told)+5 for every k>=0',
  'no_actual_residual_scalar_value_or_HF2_exact_or_total_length_or_nilcube_computed':True}
 p=out/'actual-module-length-and-single-standard-monomial-proof.json';p.write_text(json.dumps(data,separators=(',',':'))+'\n')
 r={'id':STEM,'all_passed':True,'status':'PROVED_actual_HF2_interval_without_new_coefficients','actual_HF2_Tnil_interval':[2,5],'actual_complete6_by661_rank_interval':[1,4],'actual_HF2_exact':None,'native_higher_source_or_residual_scalar_evaluation':False,'input_sha256':{p:sha(ROOT/p)for p in paths},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
