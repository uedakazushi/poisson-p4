#!/usr/bin/env python3
"""Verify two explicit918-coordinate adjoints; no solve or X2 query."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_T2_FIRST_ADJOINT52_INPUT_GATE_20261105_v63'
P={
 'primary':f'reports/{A}/report.json','cert':f'reports/{A}/actual-two-full918-source-adjoint-certificates-or-obstructions.json',
 'gate_report':f'reports/{G}/report.json','gate':f'reports/{G}/actual52-source-adjoint-input-and-two-full918-functionals.json',
 'first':'reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json',
 'profile':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json',
 'profile_peer':'reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63_REVIEW/report.json',
 'mu_gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
 'mini_peer':'reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'VERIFY_EXPLICIT_FULL918_ADJOINTS_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert all(D[k]['all_passed'] for k in ['primary','gate_report','profile_peer','mini_peer'])
 for name,folder in [('primary',A),('gate_report',G)]:
  for p,h in D[name]['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in D[name]['evidence_sha256'].items():assert sha(ROOT/'reports'/folder/p)==h
 PiY=list(map(sp,D['first']['PiY3_from102_columns']));Pi2=list(map(sp,D['profile']['Pi_remaining_top_linear_columns']));G2=list(map(sp,D['mu_gate']['full_G2_linear234_by102_columns']))
 psi=[{},{}]
 for i in range(102):
  for v in range(9):
   for j,a in PiY[i].items():
    for k,b in Pi2[j*9+v].items():psi[k]=add(psi[k],{i*9+v:a*b})
 assert list(map(enc,psi))==D['gate']['actual_two_PiY_Pi2_composed_functionals918']
 columns=[{} for _ in range(52)];pairs=[(7,7),(7,8)]
 for j,p in enumerate(G2):
  for i,c in p.items():
   row,var=divmod(i,9)
   for v in range(9):
    pair=tuple(sorted((var,v)))
    if pair in pairs:columns[2*row+pairs.index(pair)]=add(columns[2*row+pairs.index(pair)],{j*9+v:c})
 assert list(map(enc,columns))==D['gate']['actual52_adjoint_columns_in918_domain']
 M=list(map(sp,D['mu_gate']['Mnew0_columns']));gammas=[];lambdas=[]
 for k,r in enumerate(D['cert']['results']):
  assert r['channel']==k and r['exact_rr_rell_adjoint_exists']
  lam=sp(r['lambda52']);assert mv(columns,lam)==psi[k]
  assert all(i%2==1 for i in lam) # Verified support: no rr coefficient is used.
  weights={i//2:c for i,c in lam.items()}
  gamma={j:sum((weights.get(i,Q(0))*c for i,c in p.items()),Q(0)) for j,p in enumerate(M)};gamma={j:c for j,c in gamma.items() if c}
  gammas.append(gamma);lambdas.append(lam)
 assert len(lambdas)==2 and sorted(set(gammas[0])|set(gammas[1]))==[1,2,4,17,18,21,22,28,55]
 ep=o/'independent-two-full918-adjoints-and-source-rell-support.json'
 dump(o,ep.name,{'inputs':inputs,'both_true_quotient_functionals918_reconstructed':True,'all52_canonical_ordinary_rr_rell_coefficient_rows_verified':True,'both_explicit_lambda_identities_on_all918_slots_verified':True,'actual_lambda52_certificates':list(map(enc,lambdas)),'actual_Gamma_on_source56_rell':list(map(enc,gammas)),'only_source_pair_required_by_this_certificate':'rA*ell, ordinary coefficient; nofactor2','all_other_source_pairs_zero_in_this_certificate':True,'projected_Y1_from_corrected_source':'psi(Y1)=Gamma*(Xgamma2-Xp1*a1), withGammaonlyrell','linearstrand_condition_closed':'the exact row-span identity replaces any unproved associatedgraded exactness assumption','homotopy_and_unit_choices':'actual lifts differ by completeG3;psi annihilatesitsleadinglinearrelations; correctionusesactualXp, no B-linearity assumption on homotopy','no_actual_X2_M1_N1_physical_or_gauge_value_or_HF1_computed':True})
 report={'id':STEM,'all_passed':True,'both_full918_source_adjoint_certificates_and_rell_only_support_accepted':True,'projected_X2_only_sufficiency_gate_accepted':True,'actual_X2_or_HF1_T2_length_or_higher_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
