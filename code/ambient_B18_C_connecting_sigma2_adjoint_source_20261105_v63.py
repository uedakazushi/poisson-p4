#!/usr/bin/env python3
"""Exact small source adjoint for the surviving sigma2 term; no jet values."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63'
PATHS={
 'unit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'unit_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'plane':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63/all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json',
 'plane_peer':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63_REVIEW/report.json',
 'five':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63/actual-five-head-formal-RegDivell-operator-and-source-binding.json',
 'five_peer':'reports/AMBIENT_B18_C_CONNECTING_FIVE_FORMAL_SOURCE_20261105_v63_REVIEW/report.json',
 'Cbar':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mm(A,B):return[[sum((A[i][k]*B[k][j]for k in range(len(B))),Q(0))for j in range(len(B[0]))]for i in range(len(A))]
def add(a,b):
 v=a.copy()
 for k,q in b.items():
  v[k]=v.get(k,Q(0))+q
  if not v[k]:del v[k]
 return v
def scale(a,c):return{k:q*c for k,q in a.items()if q*c}
def enc(a):return[[k,str(q)]for k,q in sorted(a.items())]
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'unit':['all24_by56_m0','fixed_selected24_columns','fixed24_minor_inverse','fixed24_minor_matrix','typed_regular_function_operators','typed_m24_by56'],
   'plane':['sigma2_surviving_B_high_coordinate_adjoint','first7_source_rows','namespace_correction_to_frozen_five_formal_prose'],
   'five':['old_five_heads_free56_rows'],
   'Cbar':['CB_shift_coefficients7','r_lengths']}
 # first7_source_rows belongs to q7, not the plane companion; retain the
 # literal physical order below from its independently checked source frame.
 fields['plane'].remove('first7_source_rows')
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 u=d['unit'];m0=[[Q(x)for x in row]for row in u['all24_by56_m0']]
 selected=u['fixed_selected24_columns'];iv=[[Q(x)for x in r]for r in u['fixed24_minor_inverse']]
 M0=[[Q(x)for x in r]for r in u['fixed24_minor_matrix']]
 I=[[Q(i==j)for j in range(24)]for i in range(24)];assert mm(M0,iv)==I==mm(iv,M0)
 n0=[[Q(0)]*24 for _ in range(56)]
 for i,j in enumerate(selected):n0[j]=iv[i]
 assert mm(m0,n0)==I
 heads=d['five']['old_five_heads_free56_rows'];E=[[Q(i==h)for h in heads]for i in range(56)]
 c0=[[a-b for a,b in zip(x,y)]for x,y in zip(E,mm(n0,[[row[h]for h in heads]for row in m0]))]
 assert mm(m0,c0)==[[0]*5 for _ in range(24)]
 C=[Q(0)]*56
 for i,q in d['plane']['sigma2_surviving_B_high_coordinate_adjoint']:C[i]=Q(q)
 L=mm([C],n0)[0];assert[(i,q)for i,q in enumerate(L)if q]==[(8,Q(-43303,59400)),(12,Q(-149887,237600)),(13,Q(30713,59400))]
 assert all(L[i]==0 for i in[7,15,23])
 shifts=list(map(Q,d['Cbar']['CB_shift_coefficients7']));prefix_checks=0
 for i in[0,4,5]:
  assert shifts[i]
  for exp in S.mons(3):
   if exp[8]:continue
   upper=Q(1);lower=Q(1)
   for j,power in enumerate(exp[:7]):
    upper*= (Q(0)if j<=i else -shifts[j])**power
    lower*= (Q(0)if j<i else -shifts[j])**power
   actual=(upper-lower)/shifts[i]if exp[7]+2*sum(exp[:7])-2==2 else Q(0)
   wanted=Q(int(exp==S.mon(**{'z'+str(i):1},r=2)))
   if exp==S.mon(**{'z'+str(i):2}):wanted-=shifts[i]
   for j in range(i+1,7):
    if exp==S.mon(**{'z'+str(i):1,'z'+str(j):1}):wanted-=shifts[j]
   assert actual==wanted;prefix_checks+=1
 assert prefix_checks==495
 ops=u['typed_regular_function_operators'];m=u['typed_m24_by56'];checked=[]
 for row in[8,12,13]:
  for col in range(56):
   op=ops[m[row][col]]
   if op==['zero']:continue
   assert op[0]=='regular_coordinate_Hadamard'and op[-2:]==[1,row-8]
   p,n=ops[op[1]],ops[op[2]]
   assert p[0]==n[0]=='adapted_sector_prefix'and p[1]==n[1]and p[2]==n[2]==1 and p[3]==n[3]+1==row-7
   checked.append([row,col,m[row][col]])
 physical=[18,19,28,29,60,67,72,63,64,68,73,74]
 # The exact B component, not a source-equivalent renamed function.
 FB=[];VA4=[];VB4=[]
 for row in physical:
  FB.extend([{f'bB{row}':Q(-1)},{},{f'delta{row}':Q(1)},{f'nuA_times_bA{row}':Q(1)}])
  ba={}if row in physical[:7]else{('gA4'if row==63 else f'bA4_{row}'):Q(1)}
  VA4.extend([ba,{f'delta4_{row}':Q(1)},{},{f'KNOWN_nuB2_times_bB2_{row}':Q(-1)}])
  VB4.extend([{f'bB4_{row}':Q(-1)},{},{f'delta4_{row}':Q(1)},{}])
 FB.extend([{}for _ in range(8)]);VA4.extend([{}for _ in range(7)]+[{'gA4':Q(-1)}]);VB4.extend([{}for _ in range(8)])
 assert len(FB)==len(VA4)==len(VB4)==56
 combos=[];fourth=[]
 for h in range(5):
  f={};a={};b={}
  for j in range(56):
   f=add(f,scale(FB[j],c0[j][h]));a=add(a,scale(VA4[j],c0[j][h]));b=add(b,scale(VB4[j],c0[j][h]))
  combos.append(enc(f));fourth.append({'head':h,'A4_sigma0_primitive':enc(a),'B4_sigma0_primitive':enc(b)})
 assert combos[2]==[]
 primitives=sorted({key for pair in fourth for side in ['A4_sigma0_primitive','B4_sigma0_primitive']for key,q in pair[side]if not key.startswith('KNOWN_')})
 assert primitives==['bB4_63','delta4_18','delta4_60','delta4_63','delta4_67','delta4_72','gA4']
 data={
  'scope':'Actual constant source-adjoint calculation only. No unknown sigma jet, source3/4 coefficient, original connecting value or higher physical response is evaluated.',
  'm0_n0_both_inverse_and_cycle_checks':True,
  'actual_sigma0_56_by5_sparse':[[[i,str(c0[i][j])]for i in range(56)if c0[i][j]]for j in range(5)],
  'exact_C_sigma2_equals_minus_L_m2sigma0_plus_m1sigma1':'C*sigma2_rr=−L*(m2_rr*sigma0+m1_r*sigma1_r)',
  'fixed_L_Cn0_row24_sparse':[[i,str(q)]for i,q in enumerate(L)if q],
  'all_168_nonzero_m_operator_paths_checked':checked,
  'L_has_no_r3_r3_r2_division_rows':True,
  'recurrence_proof':'With n=Esel*(mEsel)^−1 and the literal fixed free tails, m*sigma=0, sigma1=−n0*m1*sigma0, sigma2=−n0*(m2*sigma0+m1*sigma1). This retains the derivative of n automatically; no assumption n1=0 or m0*head=0 is used.',
  'Hadamard_source_order':'All L-supported rows8,12,13 are the ordered B-sector z0,z4,z5 Hadamards, never an r^n quotient. Their pure-r second coefficients use the underlying actual B function only through total order3, with the exact quadratic B-sector coordinate shift. Unknown source4/5 is not introduced by m2 itself. Do not replace the original-z=0 evaluation by adapted-y=0.',
  'ordered_prefix_second_formula':'Let y_i=z_i+c_i*r², prefixes retain y0..y_(i−1). At ORIGINAL z=0, [r²]H_i(f)=[zi*r²]f−c_i*[zi²]f−Σ_(j>i)c_j*[zi*zj]f. This is ordinary-coefficient convention; no Hessian factors. At adapted y=0 it would have a different formula.',
  'all495_literal_prefix_monomial_identities_checked':True,
  'actual_CB_shift_c7':d['Cbar']['CB_shift_coefficients7'],
  'five_actual_B_function_sigma0_combinations':combos,
  'B_m2_highest_third_profile':'Only the z0*r²,z4*r²,z5*r² coefficients of these five combinations can be newly needed. Head2 is identically0; head1=nuA*bA67 uses existing second/first data at this order. The other3 combinations involve delta18/60/67,delta63,bB63. No blanket joint135/source3 cache coverage is asserted.',
  'coherent_V4_sigma0_primitive_columns':fourth,
  'seven_remaining_primitive_plane_fourth_functions_for_V4_sigma0':primitives,
  'fourth_monomial_scope':['r³ell','r²ell²'],
  'fourteen_primitive_slots_upper_for_V4sigma0_only':14,
  'fourth_profile_proof':'All first7 bAi are literal zAi in the actual F9 frame. The exact f63 tail gives bA63|z=0=gA. Plane delta2/bA2 vanish,nu0/plane1 vanish. Thus V4_A for tauAB is only−nuB2*bB2 (known second×second), V4_B for tauAB is0; no nu3/4 term survives. Contract all56 columns with the actual sigma0, without assuming any fourth coefficient zero.',
  'remaining_coherent_gate':'This14-slot list is ONLY the raw V4*sigma0 contribution. It must still be combined with V3*sigma1, mixed V2*sigma1 and −L*(m2*sigma0+m1*sigma1). In particular sigma1_r may retain last-row m1_r source4 contributions. No final highest-response cancellation, all-required-cache claim, worker or value GO follows from the list.',
  'available_cache_inventory':'Joint105 bivariate degree3, selected JOINT135 responses and9selected r4 adjoints exist; OWN CA57/CB57 univariate≤6 exist. None is a serialized whole sigma2 or blanket plane source4 array. Any use needs exact original IDs/branch/source and directional transport.',
  'selective_source_receipts':receipts,
 }
 owned=S.count(d)+S.count(n0)+S.count(c0)+S.count(data)+S.count(iv)+S.count(M0)+12000;assert owned<100000
 data['whole_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_SIGMA2_ADJOINT_ONLY_Z_HADAMARD_NO_NEW_HIGHER_VALUES',
  'three_selected_m2_Hadamard_rows':[8,12,13],'m2_r_division_rows_needed':False,
  'V4_sigma0_primitive_fourth_slots_upper':14,'complete_coherent_highest_cancellation_proved':False,
  'actual_connecting_values_computed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},
  'whole_numeric_upper':owned,'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)
if __name__=='__main__':main()
