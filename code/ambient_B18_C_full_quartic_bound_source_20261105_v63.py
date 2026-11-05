#!/usr/bin/env python3
"""Full-C known HB source gives HF4 in {4,5}; no new coefficient values."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_QUARTIC_BOUND_SOURCE_20261105_v63'
PATHS={
 'reduction':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json',
 'reduction_peer':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/report.json',
 'eta':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63/report.json',
 'eta_peer':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'HB':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json',
 'HB_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json',
 'HB_value':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase2-actual-RIGHT2-quartics-two-factors-and-rank.json',
 'HB_value_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['eta']['actual_full_C_HF_degree3']==6
 f3=[[Q(x)for x in row]for row in d['reduction']['full_no_cut_cubic_NF6_rows12_columns']]
 f4=[[Q(0)]*15 for _ in range(6)]
 for s in range(2):
  a1,a2,a3=[-f3[3*s+k][4*s]for k in range(3)]
  block=[[a1*a1-a2,-a1,1,0,0],[a1*a2-a3,-a2,0,1,0],[a1*a3,-a3,0,0,1]]
  for k,row in enumerate(block):f4[3*s+k][5*s:5*s+5]=list(map(Q,row))
 assert f4[0][:5]==list(map(Q,['44/625','-8/25','1','0','0']))
 assert f4[3][5:10]==list(map(Q,['6252/75625','-96/275','1','0','0']))
 matrix=[[Q(x)for x in row]for row in d['HB_value']['actual_quartic_matrix_rows_A4_B4_columns_HB0_HB1']]
 assert matrix==[[Q(-3,5),Q(6)],[Q(-121,80),Q(121,8)]]
 assert d['HB']['HB_sources']==[{'label':'HB0','formula':'ell*P0-r*P1'},{'label':'HB1','formula':'ell*P1-r*P2'}]
 assert matrix[0][0]!=0
 data={
  'scope':'FULL actual comparison C over SAME R0,D=QQ[[ell]]. This is a known-submodule bound, NOT complete quartic relations.',
  'accepted_full_HF_degrees0_1_2_3':[3,6,6,6],
  'known_full_quartic_NF6_rows15_columns':[[str(x)for x in row]for row in f4],
  'full_NF4_coordinates':['A*r²ell²','A*r ell³','A*ell4','B*r²ell²','B*r ell³','B*ell4'],
  'actual_HB_sources_without_cut':'Both HB0=ellP0−rP1 and HB1=ellP1−rP2 are combinations of actual full-C sources. Their FULL order3 cancellation uses actual Pk secondterms, not ell³H=0. No CUT1 relation is used.',
  'projection_to_old_C3_known_quartic_NF':'Select rows0,3 of fullNF4. They are EXACTLY the old A4/B4 duals including r4/r³ell coefficients. Contributions containing r ell³ or ell4 have zero under this projection; no such full coefficient is assignedzero.',
  'actual_projected_HB_matrix':[[str(x)for x in row]for row in matrix],
  'known_full_quartic_nonzero_witness':{'source':'HB0','full_NF4_row':0,'actual_coefficient':'-3/5'},
  'unconditional_full_HF_degree4_values':[4,5],
  'proof_HF4_upper5':'The full known quartic cover hasdimension6. HB0 is an ACTUAL quartic initial withnonzero fullNFrow0, so at least one dimension isremoved. ThusactualfullHF4<=5.',
  'proof_HF4_lower4':'Reductionm5C=ell*m4C andgenericD-rank4 giveHF4=mu_D(m4C)>=4.',
  'M2_DVR_torsion_shape_proof':'M=m²C is finiteD-rank4 andmu_D M=HF2=6, soTor_D M hastwo cyclicsummands. SinceellM=m³C hasmu6 (actualfullHF3), both summands haveexponent>=2. Sinceell²M=m4C hasmu<=5, atmostone exponentexceeds2. HenceTor_D(m²C)≅D/ell² ⊕D/ell^b forsomeb>=2, noncanonically asD-modules only.',
  'noncanonical_D_type_m2C':{'free_rank':4,'torsion_exponents':[2,'b>=2']},
  'conditional_b_cases':'HF4=4 iff b=2;HF4=5 iff b>=3. This doesnotdetermine the original comparison-C long exponent e.',
  'no_literal_ell3H_cut_in_fullC':True,
  'no_missing_r_ell3_or_ell4_assignedzero':True,
  'no_new_physical_gauge_IFT_jet_or_target_values':True,
  'limits':'HB isnotcompletequarticspace. Original C exponent e>=3, originalH0(D0), nilcube remainopen.',
  'source_profile':{'degree':4,'numeric_envelope':25000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'full-quartic-known-HB-witness-and-m2-torsion-shape.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_FULL_HF4_IN_4_5_AND_M2_TORSION_2_b_KNOWN_SUBMODULE_ONLY',
  'actual_full_C_HF_degree4_in_4_5_proved':True,'complete_full_quartic_space_claimed':False,'longer_C_torsion_exponent_determined':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
