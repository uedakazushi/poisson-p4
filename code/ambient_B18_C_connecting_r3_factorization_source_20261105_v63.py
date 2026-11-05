#!/usr/bin/env python3
"""Exact r^3 map on Told factors through the genuine cyclic top; no values."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_R3_FACTORIZATION_SOURCE_20261105_v63'
PATHS={
 'seven':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63/honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
 'seven_peer':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63_REVIEW/report.json',
 'old40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
 'old40_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
 'exact':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/exact-Loewy11-and-source-basis-corollary.json',
 'exact_peer':'reports/AMBIENT_B18_T_COMBINED619_QUARTIC_WITNESS_20261105_v63_REVIEW/report.json',
 'value':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/phase2-full-coherent-head23-and-low-two-values.json',
 'value_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'nonsplit_peer':'reports/AMBIENT_B18_C_CONNECTING_TOLD_NONSPLIT_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'seven':['actual_r_action7','fixed_I_Told_source_basis','honest_exact_kernel_reduction'],
  'old40':['homogeneous_module_generators40','new_actual_cubic_source'],
  'exact':['exact_Loewy_HF','actual_filtered_monomial_basis11'],
  'value':['actual_cyclic_top_chi']}
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 exp=lambda z,r:[int(j==5)*z+int(j==7)*r for j in range(8)]
 g40=d['old40']['homogeneous_module_generators40'];assert g40[-1]==[[4,exp(1,2),'1']]
 assert g40[-2]==[[3,exp(0,2),'1'],[4,exp(1,1),'11767/3652']]
 assert d['old40']['new_actual_cubic_source']['canonical_initial_relation']==g40[-1]
 assert d['exact']['exact_Loewy_HF']==[5,3,2,1]
 deep=[b for b in d['exact']['actual_filtered_monomial_basis11']if b['degree']==3]
 assert len(deep)==1 and deep[0]['target']==4 and deep[0]['monomial_original8']==exp(0,3)
 I=d['seven']['fixed_I_Told_source_basis'];assert I['first']['monomial_original8']==exp(1,0)
 assert[(b['target'],b['monomial_original8'])for b in I['other_three']]==[(4,exp(1,1)),(4,exp(0,2)),(4,exp(0,3))]
 R=[[Q(q)for q in row]for row in d['seven']['actual_r_action7']]
 assert all(R[i][j]==0 for i in range(7)for j in[0,1])
 assert R[5][2]==Q(1,4)and R[5][3]==1
 chi=list(map(Q,d['value']['actual_cyclic_top_chi']));assert chi[:4]==[0,0,0,0]and chi[4]==Q(-36,1681)
 r3=[Q(0),Q(0),Q(0),Q(0),Q(1)];assert[q/chi[4]for q in chi]==r3
 data={
  'scope':'Exact original bar-comparison source invariant. No higherheadvalue, fullKbar action table or originalH0type is computed.',
  'genuine_source_and_target':'delta:Told11→G=C[ell],sameRbar=R0/(ell),ker=Kbar=K/ellK',
  'actual_ITold_basis':['z5*h4','z5*r*h4','r²*h4','r³*h4'],
  'r2_ITold_exact_zero':'r²(z5h4)=qC=z5r²h4=0ACTUALLY:acceptedcubicinitialsourceqC+O(m4),andoldm4=0. Remainingthreebasisr²imageshaveordinarydegree≥4andvanish. Thusr²ITold=0, notjusta gradedannihilator.',
  'allfive_actual_r3_head_images_in_deep_basis':list(map(str,r3)),
  'allfive_r3_proof':'rh0,rh1∈IToldso theirr³=0. rh2−rh3/4∈IToldso r³h2=r³h3/4. MultiplytheactualG38relationr²h3+c*z5rh4+O(m3)=0byr; O(m4)=0andqC=0give r³h3=0. r³h4=gdeep is the independentactualdegree3basis. Everyheadiscontrolled.',
  'actual_deep_socle_generator':deep[0],
  'deep_socle_m_annihilated':'m*gdeep⊂oldm4=0',
  'actual_cyclic_top_chi':list(map(str,chi)),
  'exact_Rbar_map_factorization':'multiplicationr³:Told→k*gdeep equals(gdeep/d) * chi * delta, withd=−36/1681. EqualitycheckedonALL5minimalheads;bothmapsRbar-linear,mkillsthedeepandcyclictop,soallhighercoefficientspathsagree.',
  'factor_scalar':str(1/chi[4]),'actual_r3_Kbar_zero':True,
  'lifted_same_R0_containment':'SinceKbar=K/ellK HONESTLY,r³Kbar0meansr³K⊂ellK. This is NOT r³K=0 or anellannihilator/torsionexponentbound.',
  'remaining_extension':'Kbar7stillcontainsITold4withquotientJ2+k; itsfullR0moduleextension/actionsandH0DVRexponentsareunresolved.',
  'no_original_H0_length_nilcube_or_nilindex_decision':True,
  'selective_source_receipts':receipts,
 }
 owned=S.count(d)+S.count(data)+10000;assert owned<50000;data['whole_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'exact-r3-factorization-through-cyclic-top-and-Kbar-annihilation.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_EXACT_R3_TOLD_FACTORIZATION_AND_R3_KBAR_ZERO',
  'actual_r3_Kbar_zero':True,'actual_r3_K_zero_claimed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
