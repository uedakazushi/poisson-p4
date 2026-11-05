#!/usr/bin/env python3
"""All396 original-cache postcheck, polynomial P5 and direct native vectors."""
from fractions import Fraction as Q
from pathlib import Path
import json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='ac896a778c5e04ee64a1e0501f0bd530712965d1c309c73ca7a18ff971033e7b'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 e=rd(b+'actual396-Kfull-P5-and-native23-existing-cache-values.json')
 pp='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CACHE_INPUT_GATE_V2_20261105_v63_REVIEW/report.json';peer=rd(pp);assert peer['all_passed'] and pins[pp]=='fe07055d348d3ee2e77a9001c2cdcdc33f495ff7224a222030b88f05ed713870'
 gp='reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CACHE_INPUT_GATE_V2_20261105_v63/';g=rd(gp+'report.json');d=rd(gp+'pruned99-Kfull-P5-native23-and-prior396-input.json')
 src=rd('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_CALLBACK_SOURCE_20261105_v63/actual99-selected-source-K-callback-roots-and-literal-ancestry.json')
 raw=fld('reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_CORRECTED_INPUT_GATE_20261105_v63/whole-cache-eight-psi-four-quadratic-classes-source-input.json',['native_L0','native_Lr','native_Le','native_Sr','native_Se'])
 L0,Lr,Le,Sr,Se=[list(map(sp,raw[k])) for k in ['native_L0','native_Lr','native_Le','native_Sr','native_Se']]
 nids=sorted({j-640 for c in src['selected99_original768_constant_combinations'] for j,x in c if j>=640});assert len(nids)==23
 s11=fld('reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',['actual_S_source_columns128_rows48_sparse_B11'])['actual_S_source_columns128_rows48_sparse_B11'];Srr={j:{i:sp(v)[2] for i,v in s11[j] if sp(v).get(2)} for j in nids}
 s16=fld('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',['intrinsic_S2rell_columns'])['intrinsic_S2rell_columns'];Sre={j:sp(s16[j]) for j in nids}
 CF=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_CFULL_20261105_v63/actual224-Cfull-coefficients-with-native-pivot-SUM-and-Psi.json',['Cfull_two_by2_by56'])['Cfull_two_by2_by56'];C=[[{i:Q(x) for i,x in enumerate(row) if Q(x)} for row in ch] for ch in CF];G=list(map(sp,src['Gamma_two']));needed=set().union(*(set(v) for v in G),*(set(row) for ch in C for row in ch))
 original=fld('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',['all_requested164_root_values'])['all_requested164_root_values'];ids={i for i,k,n in src['literal_cache_boundaries'] if k=='source164'};indices={(0,0):0,(1,0):8,(0,1):9,(2,0):24,(1,1):25}
 leaf={i:{m:sp(v).get(k,Q(0)) for m,k in indices.items() if sp(v).get(k)} for i,v in original if i in ids};assert len(leaf)==24
 oldids=[j for j,c in src['needed_original_old_columns']];old=fld('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',['source_columns640_rows56_sparse_dual9'])['source_columns640_rows56_sparse_dual9'];oldr={j:{i:sp(v)[8] for i,v in old[j] if i in needed and sp(v).get(8)} for j in oldids}
 oldeall=fld('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])['actual_first640_ell_columns'];olde={j:{i:x for i,x in sp(oldeall[j]).items() if i in needed} for j in oldids}
 # Own positive direct-vector bound before any signed OLD/native responses.
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());assert maximum<=147;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def ba(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def bm(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 Z=(1,0)
 def wbound(A,c):
  v={}
  for j,x in c.items():
   for i,y in A[j].items():v[i]=ba(v.get(i,Z),bm(ex(x),ex(y)))
  return v
 def mbound(A,v):
  w={}
  for j,x in v.items():
   for i,y in A[j].items():
    if i in needed:w[i]=ba(w.get(i,Z),bm(ex(y),x))
  return w
 def plusbound(a,b):
  v=a.copy()
  for i,x in b.items():v[i]=ba(v.get(i,Z),x)
  return v
 def dbound(a,b):
  q=Z
  for i,x in a.items():q=ba(q,bm(ex(x),b.get(i,Z)))
  return q
 for c0 in src['selected99_original768_constant_combinations']:
  co={j:Q(x) for j,x in c0 if j<640};cn={j-640:Q(x) for j,x in c0 if j>=640};sr=wbound(Sr,cn);se=wbound(Se,cn);srr=wbound(Srr,cn);sre=wbound(Sre,cn);nr=mbound(L0,sr);ne=mbound(L0,se);nrr=plusbound(mbound(L0,srr),mbound(Lr,sr));nre=plusbound(plusbound(mbound(L0,sre),mbound(Lr,se)),mbound(Le,sr));ar=wbound(oldr,co);ae=wbound(olde,co)
  for t0 in range(2):
   ba(dbound(G[t0],nrr),dbound(C[t0][0],plusbound(ar,nr)))
   ba(dbound(G[t0],nre),ba(dbound(C[t0][0],plusbound(ae,ne)),dbound(C[t0][1],plusbound(ar,nr))))
 put('phase.json',{'phase':'independent_direct_native_vector_positive_bound_before_signed_K','prior_bits':maximum,'source147':True,'no_new_source_coefficients':True})
 actual=0
 def ck(x):
  nonlocal actual
  actual=max(actual,abs(x.numerator).bit_length(),x.denominator.bit_length());assert actual<=147;return x
 def add(a,b):
  v=a.copy()
  for i,x in b.items():
   v[i]=ck(v.get(i,Q(0))+x)
   if not v[i]:del v[i]
  return v
 def scale(a,c):return {i:ck(c*x) for i,x in a.items() if c*x}
 def mult(a,b):
  v={}
  for m,x in a.items():
   for n,y in b.items():
    k=(m[0]+n[0],m[1]+n[1])
    if k in indices:v=add(v,{k:ck(x*y)})
  return v
 jets={};skip=set(d['source_order_skip_seven_products'])
 for i,n in src['typed_program']:
  if i==94:jets[i]=None
  elif i==107 or i in skip:jets[i]={}
  elif n[0]=='actual_C0_scalar_function':assert n[1] in leaf;jets[i]=leaf[n[1]]
  elif n[0]=='sum':
   v={}
   for j,a,b0 in reversed(n[1]):assert jets[j] is not None;v=add(v,scale(jets[j],Q(a,b0)))
   jets[i]=v
  else:
   assert n[0]=='product';v={(0,0):Q(n[1],n[2])}
   for j in reversed(n[3]):assert jets[j] is not None;v=mult(v,jets[j])
   jets[i]=v
 old2={}
 for j,col in src['needed_original_old_columns']:
  q=[{},{}]
  for row,(kind,k) in col:
   if kind=='C0':
    if k==0:continue
    assert k in leaf;v=leaf[k]
   else:assert kind=='typed' and jets[k] is not None;v=jets[k]
   for u,m in enumerate([(2,0),(1,1)]):
    if v.get(m):q[u][row]=v[m]
  old2[j]=q
 def weighted(A,c):
  v={}
  for j,x in reversed(list(c.items())):v=add(v,scale(A[j],x))
  return v
 def mv(A,c):
  v={}
  for j,x in reversed(list(c.items())):v=add(v,scale({i:y for i,y in A[j].items() if i in needed},x))
  return v
 def dot(a,b):
  q=Q(0)
  for i,x in reversed(list(a.items())):q=ck(q+ck(x*b.get(i,Q(0))))
  return q
 K=[];pieces=[]
 for n,c0 in enumerate(src['selected99_original768_constant_combinations']):
  co={j:Q(x) for j,x in c0 if j<640};cn={j-640:Q(x) for j,x in c0 if j>=640};sr=weighted(Sr,cn);se=weighted(Se,cn);srr=weighted(Srr,cn);sre=weighted(Sre,cn);ar=weighted(oldr,co);ae=weighted(olde,co)
  nr=mv(L0,sr);ne=mv(L0,se);nr2=mv(L0,srr);ne2=mv(L0,sre);xr=mv(Lr,sr);xe=mv(Lr,se);ex=mv(Le,sr);nrr=add(nr2,xr);nre=add(add(ne2,ex),xe)
  a2=[weighted({j:q[k] for j,q in old2.items()},co) for k in range(2)];ks=[];ps=[]
  for t0 in range(2):
   rr=ck(dot(G[t0],add(a2[0],nrr))+dot(C[t0][0],add(ar,nr)))
   re=ck(ck(dot(G[t0],add(a2[1],nre))+dot(C[t0][1],add(ar,nr)))+dot(C[t0][0],add(ae,ne)))
   vals=[str(rr),str(re)];ks.append(vals)
   pr=[dot(G[t0],a2[0]),dot(C[t0][0],ar),dot(G[t0],nr2),ck(dot(G[t0],xr)+dot(C[t0][0],nr)),Q(0)]
   pe=[dot(G[t0],a2[1]),ck(dot(C[t0][0],ae)+dot(C[t0][1],ar)),dot(G[t0],ne2),ck(dot(G[t0],xe)+dot(C[t0][0],ne)),ck(dot(G[t0],ex)+dot(C[t0][1],nr))]
   ps.append([list(map(str,pr)),list(map(str,pe))])
  assert ks==e['K99_two_channels_rr_rell'][n] and ps==e['all99_two_by2_five_piece_values'][n];K.append(ks);pieces.append(ps)
 assert sum(Q(x)!=0 for c in K for row in c for x in row)==84
 proof={'all396_actual_Kfull_exact':True,'all1980_saved_piece_values_exact':True,'original24physical_P5_polynomial_dictionary_replay':True,'native_whole99_direct_vector_formula':'Xnr=L0Sr,Xnell=L0Sell,Xnrr=L0Srr+LrSr,Xnrell=L0Srell+LrSell+LellSr; thenGammaandCfulldirectcontractions,notproducerR/P/Qgrouping','honest94_opaqueNeverAssigned':True,'only7_order3_product_and_truegBorder3_skips':True,'fullC_Psi_and_individualpivotcarry_history_retained':True,'nonzero_coefficients84':True,'independent_direct_vector_prior_bits':maximum,'actual_max_bits':actual,'no_b1_rho_Z2_upper_or_HF_values':True,'K99_two_channels_rr_rell':K}
 ep='independent396-Kfull-original-P5-and-direct-native-vector-check.json';put(ep,proof)
 owned=count([r,e,peer,g,d,src,raw,L0,Lr,Le,Sr,Se,nids,s11,Srr,s16,Sre,CF,C,G,needed,original,leaf,old,oldr,oldeall,olde,jets,old2,K,pieces,proof,pins])+12000;assert owned<150000
 report={'id':STEM,'all_passed':True,'all396_original_source_values_and1980pieces_checked':True,'no_b1_upper_or_higher_source_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()
