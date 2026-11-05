#!/usr/bin/env python3
"""Independent all224 source reconstruction via transposed native contractions."""
from fractions import Fraction as Q
from pathlib import Path
import json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='5f6924e28aa8f314385ef80eb16d7f5e1adca5ca7a396aa72687b7068f6f0fe1'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 e=rd(b+'actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json')
 pp='reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_INPUT_GATE_20261105_v63_REVIEW/report.json';peer=rd(pp);assert peer['all_passed'] and pins[pp]=='0b45a2f31dd9630d754a88989a84488dd30a505997bdbabf3971b0a4f0606b27'
 cp='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json'
 raw=fld(cp,['Gamma_two','native_L0','native_Le','native_S0','native_Se','Q_dual_two6_rows','pivot27_source_indices'])
 G=list(map(sp,raw['Gamma_two']));L0=list(map(sp,raw['native_L0']));Le=list(map(sp,raw['native_Le']));S0=list(map(sp,raw['native_S0']));Se=list(map(sp,raw['native_Se']));Qleft=list(map(sp,raw['Q_dual_two6_rows']))
 chart='reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json';ch=fld(chart,['free30_V_indices','fixed_bottom30_inverse_columns','original27_X_source_columns'])
 ids=ch['original27_X_source_columns'];assert ids==raw['pivot27_source_indices'] and len(ids)==27
 P=[{} for _ in range(56)];Pi=[{} for _ in range(56)]
 for j,v in zip(ch['free30_V_indices'],map(sp,ch['fixed_bottom30_inverse_columns'])):P[j]={i:x for i,x in v.items() if i<27};Pi[j]={i-27:x for i,x in v.items() if i>=27}
 invp='reports/AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_20261105_v63/actual-two29-pairings-and-recovered30-kernel-invariants.json';inv=fld(invp,['actual_all_three_point_kernel_direction_invariants'])['actual_all_three_point_kernel_direction_invariants'];weights=[[Q(inv[str(j)][k]) for j in [23,29,30]] for k in range(2)]
 old=fld('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',['source_columns640_rows56_sparse_dual9'])['source_columns640_rows56_sparse_dual9']
 native=fld('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',['source_columns128_rows48_sparse_dual9'])['source_columns128_rows48_sparse_dual9']
 oldell=fld('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])['actual_first640_ell_columns']
 Xr=[];oldE={}
 for p,j in enumerate(ids):
  c=old[j] if j<640 else native[j-640];Xr.append({i:sp(v)[8] for i,v in c if sp(v).get(8)})
  if j<640:oldE[p]=sp(oldell[j])
 # Source-only positive envelope for this DIFFERENT, transposed contraction.
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());assert maximum<=147;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def bd(a,b):
  q=(1,0)
  for j,x in a.items():q=add(q,mul(ex(x),ex(b.get(j,Q(0)))))
  return q
 GLb=[[[bd(row,c) for c in L] for L in (L0,Le)] for row in G]
 XPb=[]
 for t0,row in enumerate(G):
  rs=[bd(row,x) for x in Xr];es=[]
  for p,j in enumerate(ids):
   if j<640:es.append(bd(row,oldE[p]))
   else:
    q=(1,0)
    for m,c in Se[j-640].items():q=add(q,mul(GLb[t0][0][m],ex(c)))
    for m,c in S0[j-640].items():q=add(q,mul(GLb[t0][1][m],ex(c)))
    es.append(q)
  XPb.append([rs,es])
 for t0 in range(2):
  for k in range(2):
   for j in range(56):
    q=(1,0)
    for i,x in Pi[j].items():
     q=add(q,mul(ex(Qleft[t0].get(2*i+(1 if k==0 else 0),Q(0))),ex(x)))
     if t0==0:q=add(q,mul(ex(weights[k][i]),ex(x)))
    for p,x in P[j].items():q=add(q,mul(XPb[t0][k][p],ex(x)))
 put('phase.json',{'phase':'independent_transposed_positive_bound_before_signed_values','positive_prior_bits':maximum,'hard_source147':True,'no_K_or_upper_relation_values':True})
 bits=0
 def ck(x):
  nonlocal bits
  bits=max(bits,abs(x.numerator).bit_length(),x.denominator.bit_length());assert bits<=147;return x
 def dot(a,b):
  q=Q(0)
  for j,x in a.items():q=ck(q+ck(x*b.get(j,Q(0))))
  return q
 # Contract Gamma with L first (producer multiplies L with S first).
 GL=[[[dot(row,c) for c in L] for L in (L0,Le)] for row in G];XP=[]
 for t0,row in enumerate(G):
  rs=[dot(row,x) for x in Xr];es=[]
  for p,j in enumerate(ids):
   if j<640:es.append(dot(row,oldE[p]))
   else:
    q=Q(0)
    for m,c in S0[j-640].items():q=ck(q+ck(GL[t0][1][m]*c))
    for m,c in Se[j-640].items():q=ck(q+ck(GL[t0][0][m]*c))
    es.append(q)
  XP.append([rs,es])
 assert [[[str(x) for x in row] for row in cc] for cc in XP]==e['Gamma_Xp1_two_by2_by27']
 C=[]
 for t0 in range(2):
  channel=[]
  for k in range(2):
   row=[]
   for j in range(56):
    # Reversed grouping: all pivot corrections, then invariants and Psi.
    q=Q(0)
    for p,x in reversed(list(P[j].items())):q=ck(q-ck(XP[t0][k][p]*x))
    for i,x in reversed(list(Pi[j].items())):
     if t0==0:q=ck(q+ck(weights[k][i]*x))
     q=ck(q+ck(Qleft[t0].get(2*i+(1 if k==0 else 0),Q(0))*x))
    row.append(str(q))
   channel.append(row)
  C.append(channel)
 assert C==e['Cfull_two_by2_by56']
 proof={'all224_actual_Cfull_coefficients_match':True,'all108_Gamma_Xp_plane_entries_match':True,'original_old640_native128_r_and_full9_oldell_callbacks_reused':True,'native_ell_independent_transposed_contraction':'GammaLe*S0 + GammaL0*Se, both summed; producer used Gamma*(L0Se+LeS0)','point_chart_original27_plus3_and_all_six_invariants':True,'Psi_rell_for_r_Psi_rr_for_ell_reconstructed_from_original_Q_left_inverse':True,'summation_order_independent_reversed_pivots_then_invariants_then_Psi':True,'independent_positive_prior_bits':maximum,'actual_coefficient_max_bits':bits,'source_scope':'coherent existence comparison frame; complementary canonicalM1 not assigned','no_K_Xstar2_Xc3_or_upper_relation_value_computed':True,'Cfull_two_by2_by56':C}
 ep='independent-all224-Cfull-source-and-transposed-native-reconstruction.json';put(ep,proof)
 owned=count([r,e,peer,raw,G,L0,Le,S0,Se,Qleft,ch,P,Pi,inv,weights,old,native,oldell,Xr,oldE,GLb,XPb,GL,XP,C,proof,pins])+8000;assert owned<60000
 report={'id':STEM,'all_passed':True,'all224_Cfull_source_values_independently_checked':True,'no_K_or_upper_relation_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'positive_prior_bits':maximum,'actual_max_bits':bits,'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()
