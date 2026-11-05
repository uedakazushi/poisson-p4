#!/usr/bin/env python3
"""Independent integer-lattice replay of the two frozen raw33 pairings."""
import hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63'
STEM=BASE+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={}
 def read(p):
  pins[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and sha(ROOT/rp)=='849e75ab49a89ef58505bc0c8377e257e6d7ed9b529e79f75a810a8a5086a444'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 ep=f'reports/{BASE}/actual-two29-pairings-and-recovered30-kernel-invariants.json';v=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 g=read('reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_INPUT_GATE_20261105_v63/report.json');d=read('reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_INPUT_GATE_20261105_v63/two-pairing-source-binding-and-positive-height-input.json')
 pr=read('reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_INPUT_GATE_20261105_v63_REVIEW/report.json');assert pr['all_passed'] and pr['two_dot_source_and_worker_condition_closed']
 assert sha(ROOT/'reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_INPUT_GATE_20261105_v63/report.json')==pr['source_gate_sha256']
 # Rebind the primitive input independently to the actual physical63 mask-tauA column.
 p=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')['physical12_coefficients'][7]
 assert p['physical_row']==63 and p['roots']==d['source_raw63_roots']
 vec=[dict(p['delta']),{}, {i:str(-Q(x)) for i,x in p['bB']}]
 assert [dict(x) for x in d['primitive_V29_raw63_mask1']]==vec
 rho=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json')['actual_two_three_component_bar_rho33'];assert rho==d['two_universal_bar_rho33']
 # A single common integer lattice per channel, different from sequential QQ summation.
 receipts=[];answers=[];component=[]
 for ch,rr in enumerate(rho):
  terms=[(target,i,Q(a)*Q(vec[target].get(i,'0'))) for target,col in enumerate(rr) for i,a in col]
  D=math.lcm(*(x.denominator for _,_,x in terms));ints=[x.numerator*(D//x.denominator) for _,_,x in terms];n=sum(ints);ans=Q(n,D)
  comp=[sum((x for target,_,x in terms if target==j),Q(0)) for j in range(3)]
  assert str(ans)==r['actual_alpha29_beta29'][ch]==v['actual_alpha29_beta29'][ch]
  assert list(map(str,comp))==v['two_pairing_component_sums'][ch]
  capD,capU=map(int,d['two_new_dot_positive_bounds'][ch]);assert capD%ans.denominator==0 and abs(ans.numerator)*(capD//ans.denominator)<=capU
  assert max(D.bit_length(),abs(n).bit_length(),*(abs(z).bit_length() for z in ints))<=g['source_committed_bit_bound']==76
  receipts.append({'channel':ch,'integer_common_denominator':str(D),'integer_numerator_sum':str(n),'term_count':len(terms),'all_integer_terms':list(map(str,ints))});answers.append(ans);component.append(list(map(str,comp)))
 recovered=[]
 for ch,z in enumerate(answers):
  # Rearranged exact linear relation: (25/2)u30=uW+(3375/11)u23+50u29.
  z30=(Q(d['known_alpha_beta_on_Wv'][ch])+Q(3375,11)*Q(d['known_alpha23_beta23'][ch])+50*z)*Q(2,25)
  assert str(z30)==r['actual_alpha30_beta30_recovered'][ch]==v['actual_alpha30_beta30_recovered_without_new_pairing'][ch];recovered.append(str(z30))
 assert v['only_TWO_new_pairings_evaluated'] and r['new_pairing_count']==2 and r['actual_upper619_values_or_HF2_not_computed']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 e=out/'independent-integer-lattice-two29-and-linear30-recovery.json';e.write_text(json.dumps({'integer_lattice_receipts':receipts,'actual_alpha29_beta29':list(map(str,answers)),'recovered_alpha30_beta30':recovered,'literal_physical63_mask_tauA_bound':True,'all_component_sums_verified':True,'V30_not_paired':True,'scope':'Two M1 invariant values and algebraic recovery only; no actual upper619 relation, HF2, or total length claim.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'primary_sha256':sha(ROOT/rp),'source_and_value_replay':True,'actual_alpha29_beta29':list(map(str,answers)),'actual_alpha30_beta30_recovered':recovered,'new_response_solve':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
