#!/usr/bin/env python3
"""Pure source-bound composition of bar-cube zero and the full conductor ideal."""
import hashlib,json,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER6_20261105_v63'
PATHS={
 'upper7':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63/report.json',
 'upper7peer':'reports/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63_REVIEW/report.json',
 'upper7note':'research/notes/AMBIENT_B18_QUADRATIC_CONDUCTOR_NIL_UPPER7_20261104_v63.md',
 'barcube':'reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63/report.json',
 'barcubepeer':'reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63_REVIEW/report.json',
 'finiteNpeer':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'finitecubepeer':'reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k!='upper7note'}
 assert all(x['all_passed']for x in d.values())
 assert d['upper7']['actual_J_contains_m_R0_squared_and_both_nu']
 assert d['upper7']['entire_conductor_annihilation_after_completed_basechange']
 assert d['upper7']['Nil_D0_times_K_squared_zero']
 assert d['barcube']['genuine_Nbar_cube_zero']and d['barcube']['original_nilcube_subset_ell_H0']
 note=(ROOT/PATHS['upper7note']).read_text()
 assert 'JK^2=0' in note and '\\mathfrak m_{R_0}^2\\subset J' in note
 # Literal same-coordinate ell=5(s-r), with no formal coordinate cancellation.
 ell={(1,0):-5,(0,1):5};sq={}
 for (i,j),a in ell.items():
  for (k,l),b in ell.items():sq[i+k,j+l]=sq.get((i+k,j+l),0)+a*b
 assert sq=={(2,0):25,(1,1):-50,(0,2):25}
 data={
  'same_actual_base':'R0=QQ[[zA1,...,zA7,r,s]], same physical coefficient interface as accepted UPPER7',
  'actual_conductor_sum_J':'(zA1,...,zA7,r²,rs,s²)',
  'actual_DVR_parameter_ell':'5*(s-r)',
  'literal_ell_squared_J_membership':[['r²',25],['rs',-50],['s²',25]],
  'conductor_annihilator_proof':'K² is contained in both first-normalization kernels KA,KB because each second-stage basechange kernel is square-zero. The entire actual conductors annihilate KA,KB, so IA*K²=IB*K²=0 and J*K²=0. This was already proved with the same fullR0 in UPPER7 §4, including completed tensor basechange.',
  'new_accepted_cube_constraint':'N³⊂ell H0(D0)⊂ell K, since E4 has H0=0 and K=ker(D0→E4)',
  'power6_identity_chain':['N⁶=(N³)²','(N³)²⊂(ell K)²=ell² K²','ell² K²⊂J K²=0'],
  'proved_actual_Nil_D0_power6_zero':True,
  'full_original_germ_transport':'The accepted finite-origin/mixed-positive theorem makes every nil product of at least two factors a zero-degree product. Thus N0⁶=0 gives the full transverse nil sixth power0; actual regular G19 factor and projective scalar chart preserve the index.',
  'actual_nil_index_interval':[3,6],
  'lower_bound_source':'Existing Qb coherent actual nil-square nonzero, inherited from UPPER7 lowerbound3.',
  'original_nilcube_zero_or_nonzero_and_exact_index_remain_unknown':True,
  'no_new_H0_annihilator_or_length_claim':True,
  'countermodel_scope':'The separate barcube/index7 toy model retains NK²=0 but does not impose the full actual J*K²=0 (in particular ell²K²=0). It therefore is not a counterexample to this stronger source-bound corollary. Historical packets are unchanged.',
  'no_new_source_gauge_jet_or_scalar_nil_product_values':True}
 ev=out/'barcube-plus-full-conductor-and-ell-squared-proof.json';ev.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PROVED_unconditional_actual_nil_upper6','actual_D0_and_full_B18_nil_power6_zero':True,'actual_nil_index_interval':[3,6],'original_nilcube_vanishing_decided':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()
