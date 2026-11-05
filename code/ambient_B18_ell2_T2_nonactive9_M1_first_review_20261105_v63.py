#!/usr/bin/env python3
"""Independent cached-first postcheck; no M1 or invariant response value."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,add,image
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and r['actual_M1_psi_or_quadratic_values_not_evaluated']
 assert r['code_sha256']==sha(ROOT/'code/ambient_B18_ell2_T2_nonactive9_M1_first_20261105_v63.py')
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 source='reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_M1_FIRST_INPUT_GATE_V2_20261105_v63'
 own=read(source+'_REVIEW/report.json');assert own['all_passed']
 d=read(source+'/pruned-known-firsts-and-positive-M1-four-invariant-test-input.json')
 ep=f'reports/{BASE}/actual-known-first-weighted-M1-operators-and-four-invariant-tests.json';v=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 C={j:sp(x) for j,x in d['selected99_original768']};Y={j:sp(x) for j,x in d['selected99_Y0']};R={j:sp(x) for j,x in d['original_r_first_columns']};E={j:sp(x) for j,x in d['original640_ell_first_columns']}
 L=list(map(sp,d['native_L0_columns']));Se={j:sp(x) for j,x in d['native_Se_columns']};S0=list(map(sp,d['native_S0_columns']));inv=list(map(sp,d['fixed_B30_inverse_columns']));piv=list(map(sp,d['point_Xp27_columns']));free=d['free30_V_indices'];M=list(map(sp,d['M0_columns']));G2=list(map(sp,d['G2linear_columns']));lam=list(map(sp,d['lambda_two']))
 for j,col in Se.items():E[640+j]=image(L,col)
 Z={};checks=0
 for j,c in C.items():
  assert not image(S0,{k-640:x for k,x in c.items() if k>=640});zz=[]
  for var,A in [(7,R),(8,E)]:
   x=image(A,c);a=image(inv,{i:x[row] for i,row in enumerate(free) if row in x});a={i:v for i,v in a.items() if i<27};z=add(x,image(piv,a),Q(-1))
   actual=image(M,z);rhs=image([{i//9:v for i,v in col.items() if i%9==var} for col in G2],Y[j]);assert actual==rhs;checks+=26;zz.append(z)
  Z[j]=zz
 assert len(Z)==34
 weighted={}
 for rec,result in zip(d['formal_eighteen_weights'],v['all18_results'],strict=True):
  assert (rec['original418_source'],rec['contraction'])==(result['original418_source'],result['contraction'])
  w=[[{},{}],[{},{}]]
  for slot,c in rec['weights']:
   j,ch=divmod(slot,2)
   for direction in range(2):w[ch][direction]=add(w[ch][direction],Z[j][direction],Q(c))
  assert w==[[sp(x) for x in row] for row in result['weighted_Xeff_two_channels_r_ell']]
  op={}
  for ch,row in enumerate(lam):
   for target,c in row.items():
    for direction,vec in [(0,w[ch][1]),(1,w[ch][0])]:
     for source,x in vec.items():op=add(op,{direction*26*56+target*56+source:c*x})
  assert not op and all(op.get(i,Q(0))==0 for i in range(2912))
  assert result['represented_by_four_known_invariant_functionals'] and list(map(Q,result['proposed_four_coefficients']))==[Q(0)]*4
  assert result['nonzero_residual_slot_count']==0 and result['first_exact_obstruction'] is None and result['all2912_slots_checked']
  weighted[rec['original418_source'],rec['contraction']]=w
 for source in {j for j,ch in weighted}:
  A,B=weighted[source,'A'],weighted[source,'B'];assert not A[0][0] and not A[1][0] and not A[1][1] and not add(A[0][1],B[0][0]) and B[1]==A[0]
 assert checks==v['all34_full_source_chain_equations_checked']==1768
 assert r['known_M1_operator_representations_succeeded']==18 and r['known_M1_operator_representations_failed']==0
 assert v['no_M1_Phi_native_S2_pivot_G2quadratic_or_actual_invariant_values_queried']
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent34-source-and18-literal-zero-M1-operator-checks.json';e.write_text(json.dumps({'all34_source_r_ell_chain_checks':checks,'all18_times2912_M1_coefficient_slots_literal_zero':True,'all18_proposed_four_coefficients_zero':True,'all9_full56_W2_plane_zero':True,'no_M1_invariant_or_quadratic_response_values':True,'scope':'Exactknownfirstcoefficientidentities.Noactualpsi,old/nativeS2response,HF2orT2lengthcomputed.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all18_M1_operator_rows_literal_zero':True,'no_unknown_response_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'literal_zero_operators':18}))
if __name__=='__main__':main()
