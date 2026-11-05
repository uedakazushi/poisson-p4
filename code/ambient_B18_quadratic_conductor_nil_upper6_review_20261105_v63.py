#!/usr/bin/env python3
"""Independent full-conductor and genuine-barcube composition, no new values."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER6_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER6_20261105_v63/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p in r['input_sha256']:
  if p.endswith('/report.json'):assert read(p)['all_passed']
 source=read('reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER6_20261105_v63/barcube-plus-full-conductor-and-ell-squared-proof.json')
 old=read('reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63/actual-trace-first-gate-and-entire-conductor-proof.json')
 assert old['actual_coefficient_base']=='same physical regular9 R0 throughout'
 assert old['actual_joint_ideal']=='J=IA+IB=(zA7,r²,rs,s²)'
 assert old['m_R0_square_subset_J']and old['K_squared_subset_both_first_normalization_kernels']and old['J_annihilates_K_squared']
 con=read('reports/AMBIENT_B18_NIL_CONNECTING16_SOURCE_OPERATOR_20261105_v63/actual-formal-connecting16-and-ten-source-obstructions.json');assert con['ell_equals']=='5*(sA-rA)'
 finite=read('reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json');assert finite['whole_E4_and_Ered_ell_regularity_all6_Ass_accepted']
 bar=read('reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63_REVIEW/report.json');assert bar['genuine_Nbar_cube_zero']and bar['original_cube_subset_ell_H0']
 # Separate symbolic polynomial multiplication in the two actual plane parameters.
 p={(1,0):-5,(0,1):5};square={}
 for a,c in p.items():
  for b,d in p.items():
   k=tuple(x+y for x,y in zip(a,b));square[k]=square.get(k,0)+c*d
 assert square=={(2,0):25,(1,1):-50,(0,2):25}
 assert source['literal_ell_squared_J_membership']==[['r²',25],['rs',-50],['s²',25]]
 assert r['actual_nil_index_interval']==[3,6]and not r['original_nilcube_vanishing_decided']
 data={'same_full_R0_J_annihilator_checked':True,'explicit_ell_squared_in_J':[['r²',25],['rs',-50],['s²',25]],
 'H0_in_K_proof':'H0=tors_D(D0). The actual target E4 is ell-regular, so every D-torsion element has zero image in E4; hence H0⊂K=ker(D0→E4). This also follows from the accepted D-free quotient N/K.',
 'full_conductor_proof':'Both orders of normalization show K² lies in the corresponding first kernel. A conductor f∈IA gives an AA-linear multiplication map BA→AA; inclusion after this map is multiplication by f. Tensor and completion preserve that composition, so f annihilates the first basechange kernel. Likewise IB. Thus J=IA+IB annihilates K² as a fullR0 scalar ideal, not merely NK².',
 'exact_ideal_power_argument':['N³⊂ellH0⊂ellK','N⁶=(N³)²⊂(ellK)²=ell²K²','ell²∈J, JK²=0 ⇒N⁶=0'],
 'full_germ_transfer':'Inherited accepted finite-origin/mixed-positive nil product theorem and regular G19/projective-scalar transport from UPPER7; no new specialization or reduced-tail assumption.',
 'actual_unconditional_nil_index_interval':[3,6],'known_N_square_nonzero_inherited':True,'original_cube_zero_or_exact_index3_not_proved':True,
 'index7_toy_scope':'Toy model deliberately omitted J K²=0 and has ell²K²≠0, so it does not contradict this stronger actual-source corollary. No old frozen file is changed.',
 'new_source_gauge_jet_or_nil_product_values':False}
 ep=out/'independent-full-conductor-and-actual-nil-sixth-power-proof.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_independent_actual_original_nil_sixth_power_zero','actual_original_and_full_germ_Nil_power6_zero':True,'actual_nil_index_interval':[3,6],'original_cube_vanishing_decided':False,'input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
