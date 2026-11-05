#!/usr/bin/env python3
"""Positive envelopes only for the finite C72 native70 source program."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_NATIVE70_C72_POSITIVE_SOURCE_20261105_v63'
ZERO=(1,0)
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={};maximum=0;operations=0;cursor={}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 def mark(v):
  nonlocal maximum,operations
  operations+=1;maximum=max(maximum,v[0].bit_length(),v[1].bit_length())
  if maximum>384:
   (out/'positive-guard-frontier.json').write_text(json.dumps({'cursor':cursor,'positive_bound':list(map(str,v)),'bits':maximum,'no_signed_response_or_retry':True})+'\n');raise RuntimeError('positive384guard')
  return v
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else ZERO
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else ZERO
 def vadd(a,b):return [plus(x,y) for x,y in zip(a,b,strict=True)]
 def vscale(a,c):e=ex(c);return [mul(x,e) for x in a]
 base='reports/AMBIENT_B18_ELL2_T2_NATIVE70_C72_FRAME_SOURCE_20261105_v63'
 r=read(base+'/report.json');assert r['all_passed'] and r['regular_G8_program_not_evaluated']
 frame=read(base+'/actual-C72-cache-frame-and-nine-direction-Nfirst-source-bundle.json');assert sha(ROOT/base/'actual-C72-cache-frame-and-nine-direction-Nfirst-source-bundle.json')==r['evidence_sha256']['actual-C72-cache-frame-and-nine-direction-Nfirst-source-bundle.json']
 kb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_A18_SOURCE_20261105_v63';kr=read(kb+'/report.json');assert kr['all_passed'] and kr['all70_fullC4_point_kernel_verified']
 kernel=read(kb+'/all70-C4-kernels-and-full-raw33-A18-ring-map-source.json');recs=kernel['all70_post27_native_point_kernel_records'];assert len(recs)==70 and [a['case619'] for a in recs]==list(range(549,619))
 mons=list(map(tuple,frame['A18_basis_monomials']));mi={m:i for i,m in enumerate(mons)};assert len(mons)==18 and frame['trace_basis_mask_order4']==['1','tauA','tauB','tauA_tauB']
 prod=[[mi.get(tuple(sorted(a+b))) for b in mons] for a in mons]
 nu=[{i:Q(x) for i,x in enumerate(row) if Q(x)} for row in frame['actual_nu_A_B_A18']]
 # Zero by monomial support, without evaluating signed products of cache functions.
 assert all(prod[i][j] is None for a in nu for b in nu for i in a for j in b)
 assert all(not m or len(m)<=2 for m in mons)
 table={};table_records=[]
 for i in range(72):
  for j in range(72):
   u,v=i//4,j//4;p=prod[u][v]
   if p is None:continue
   mask=i%4^j%4;carry=i%4&j%4
   if carry==0:terms=[(4*p+mask,Q(1))]
   elif carry==3:continue
   else:terms=[(4*prod[p][b]+mask,c) for b,c in nu[carry-1].items() if prod[p][b] is not None]
   if terms:table[i,j]=terms;table_records.append([i,j,[[a,str(b)] for a,b in terms]])
 def vmul(a,b):
  v=[ZERO]*72
  for i,x in enumerate(a):
   if not x[1]:continue
   for j,y in enumerate(b):
    if not y[1]:continue
    for k,c in table.get((i,j),[]):v[k]=plus(v[k],mul(mul(x,y),ex(c)))
  return v
 one=[ZERO]*72;one[0]=ex(1)
 kap=[[ex(x) for x in row] for row in frame['actual_kappa9_C72']]
 assert all(Q(row[0])==0 for row in frame['actual_kappa9_C72'])
 prog=frame['regular_G8_source_program'];core=[kap[i] for i in prog['core_kappa_indices4']]
 powers=[]
 for i,a in enumerate(core):
  cursor.update(stage='positive_core_powers',core=i);p=[one,a]
  for e in range(2,7):cursor['power']=e;p.append(vmul(p[-1],a))
  assert not any(x[1] for x in vmul(p[-1],a));powers.append(p)
 def poly(terms):
  z=[ZERO]*72
  for exp,n,d in terms:
   # The maximal ideal J has J^7=0 (base^3=0, two trace generators).
   if sum(exp)>=7:continue
   term=one
   for i,e in enumerate(exp):
    if e:term=vmul(term,powers[i][e])
   z=vadd(z,vscale(term,Q(n,d)))
  return z
 def inv(k,P):
  if k==0:return one
  z=[ZERO]*72
  for e in range(4):
   # Absolute binomial coefficients; the alternating P sign is irrelevant to a positive bound.
   c=Q(math.comb(k+e-1,e),5**(k+e)) if P else Q(math.comb(k+e-1,e))
   z=vadd(z,vscale(powers[1][2*e],c))
  return z
 G=[]
 for k,(name,idx) in enumerate(zip(prog['five_pivot_order'],prog['five_pivot_kappa_indices'],strict=True)):
  cursor.update(stage='positive_five_regular_roots',root=name);q=prog['five_actual_rational_roots'][name]
  root=vmul(vmul(poly(q['numerator']),inv(q['D_power'],False)),inv(q['P_power'],True))
  G.append(vadd(kap[idx],root))
 for k,q in enumerate(prog['regular_CI_generators3']):cursor.update(stage='positive_regular_CI',generator=k);G.append(poly(q))
 assert len(G)==8
 raw_G_bounds=[[[str(D),str(U)] for D,U in row] for row in G]
 # These are the same source functions on the same kappa. The exact old
 # point/first coefficients are already pinned; using their absolute values
 # tightens final bounds, never claims evaluation of an unknown second jet.
 for i in range(8):
  for j,x in enumerate(frame['actual_G8_point_and_first40'][i]):G[i][j]=ex(x)
 sop=read('reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',['native_S12_by32'])['native_S12_by32']
 sdefs={j:[] for j in range(28)}
 for row,col,gi,sign in sop['koszul28_in_first8_rows']:sdefs[col].append((row,gi,sign))
 ga=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json',['Gamma4_on_source56_times8'])
 Gamma=list(map(sp,ga['Gamma4_on_source56_times8']));assert len(Gamma)==4
 needed_g=sorted({(slot//8)//4 for row in Gamma for slot in row});assert needed_g==[0,4,5,7]
 N={}
 for intrinsic in range(8):
  for output in needed_g:
   p=frame['actual_N12_by12_point_and_first40'][intrinsic][output]
   assert max(max(abs(Q(x).numerator).bit_length(),Q(x).denominator.bit_length()) for x in p)<=15
   N[intrinsic,output]=[ex(x) for x in p]+[ZERO]*32
 # At most one case's eight S and four output vectors are live. Keep no
 # all70 dense vector duplicate in the source builder.
 static_owned=count([frame,kernel,ga,sop,table_records])+2*count(table_records)+2*(9*72+28*72+8*72+32*72+12*72)+18000
 assert static_owned<150000
 smax=[];fourbounds=[]
 for rec in recs:
  case=rec['case619'];cursor.update(stage='positive_weighted_S_then_Nfirst',case619=case)
  cn=sp(rec['native128_constant_projection']);assert cn and max(cn)<112
  assert rec['all48_intrinsic_native_S0_components_exact_zero']
  S=[[ZERO]*72 for _ in range(8)]
  for idx,c in cn.items():
   col,mask=divmod(idx,4);tau=[ZERO]*72;tau[mask]=ex(c)
   for row,gi,sign in sdefs[col]:S[row]=vadd(S[row],vmul(tau,G[gi]))
  # Full C4 source kernel, not M0 or old/native cancellation.
  for row in S:row[:4]=[ZERO]*4
  X={output:[ZERO]*72 for output in needed_g}
  for output in needed_g:
   for intrinsic in range(8):X[output]=vadd(X[output],vmul(N[intrinsic,output],S[intrinsic]))
  vals=[]
  for row in Gamma:
   b=ZERO
   for slot,c in row.items():
    source,pair=divmod(slot,8);g,mask=divmod(source,4);b=plus(b,mul(ex(c),X[g][4*(10+pair)+mask]))
   vals.append(list(map(str,b)))
  fourbounds.append({'case619':case,'four_native_partial_positive_bounds':vals})
  smax.append({'case619':case,'maximum_S_coefficient_bound_bits':max(max(D.bit_length(),U.bit_length()) for row in S for D,U in row)})
 owned=max(static_owned,count([frame,kernel,ga,sop,table_records,kap,powers,G,N,fourbounds,smax])+2*count(table_records)+12000)
 assert owned<150000
 future=count([frame['actual_kappa9_C72'],frame['actual_nu_A_B_A18'],prog,frame['actual_G8_point_and_first40'],[[i,j,frame['actual_N12_by12_point_and_first40'][i][j]] for i,j in N],recs,ga,table_records])+2*(8*72+28*72+8*72+4*72)+10000
 ep=out/'positive-C72-program-and-native70-partial-coefficient-envelopes.json'
 ev={'A18_basis_monomials':mons,'C72_nonzero_multiplication_records':table_records,'all_nu_pair_products_zero_by_support':True,'J7_zero_and_binomial_inverse_power3_source_proof':True,'needed_N_intrinsic_rows':list(range(8)),'needed_N_output_columns':needed_g,'all70_free4_source_parts_absent':True,'actual_N_needed_first_entries_maximum_bits':15,'regular_G8_raw_positive_bounds8_by72':raw_G_bounds,'G8_point_first_tightening_uses_same_frozen_functions':True,'all70_native_partial_positive_bounds':fourbounds,'all70_S_positive_height_diagnostics':smax,'source_committed_maximum_positive_bits':maximum,'future_temporary_bound_3B_plus1':3*maximum+1,'positive_operation_count':operations,'source_owned_conservative_bound':owned,'future_native_only_owned_forecast':future,'extension_scope':'Only weighted S(cN) with full C4 point kernel; actual N2 is not set to zero. Any change in base quadratic scalar coefficients annihilates S(cN) because m_base^3=0. Individual pivots retain actual L0Se+LeS0.','remaining_whole_gate':'Old70 source-second and all first-cross responses, actual27 individual first-pivot carries, full combined positive envelopes/owned/raw/height, source peer and ROOT value GO.','no_signed_G8_second_native70_or_pairing_evaluated':True}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'positive_envelopes_only_no_native70_values':True,'source_owned_bound':owned,'future_native_only_owned_forecast':future,'source_committed_bits':maximum,'source_temporary_bits':3*maximum+1,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'positive_bits':maximum,'source_owned':owned,'native_only_forecast':future,'seconds':report['seconds']}))
if __name__=='__main__':main()
