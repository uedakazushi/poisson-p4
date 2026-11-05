#!/usr/bin/env python3
"""Independent shared27 pullback, all1836 alpha slots and836 G3 identities."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63_REVIEW'
def main(expected):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(o/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']==expected
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 e=rd(b+'actual-two-full918-alpha-covectors-and836-G3-identities.json')
 sr='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_CLEARING_INPUT_GATE_20261105_v63_REVIEW/report.json';assert rd(sr)['all_passed'] and pins[sr]=='224a2839ec0840508295ee3fa0e6bdff683041562f38cc443ccbdcceeba4e5d6'
 rows=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63/actual-four-specific-B99-point-dual-rows-STOP-only384.json',['actual_four_rows_channel_pair'])['actual_four_rows_channel_pair'];S={(v['channel'],v['pair_rr_rell']):sp(v['point_row102']) for v in rows}
 pre=[{**{9*i+7:x for i,x in S[t0,1].items()},**{9*i+8:x for i,x in S[t0,0].items()}} for t0 in range(2)]
 h='reports/AMBIENT_B18_ELL2_T2_FIRST_HILBERT_INPUT_GATE_20261105_v63/actual-known-G3linear-and-bottom3-first-Hilbert-inputs.json';z=fld(h,['PiY3_from102_columns','actual_full_G3linear102_by418_columns']);Pi=list(map(sp,z['PiY3_from102_columns']));G=list(map(sp,z['actual_full_G3linear102_by418_columns']))
 sel=fld('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/phase1-known-G3linear-fixed-source-rank.json',['selected_columns'])['selected_columns'];J=[G[i] for i in sel]
 F=list(map(sp,fld('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_FIXED_ROWS_20261105_v63/actual-fixed-specific-dual-rows-and-F25-with-next-clearing-prior.json',['actual_F25_columns27'])['actual_F25_columns27']))
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());assert maximum<=384;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def add(a,b0):
  if not a[1]:return b0
  if not b0[1]:return a
  D=math.lcm(a[0],b0[0]);return mark((D,a[1]*(D//a[0])+b0[1]*(D//b0[0])))
 def mul(a,b0):return mark((a[0]*b0[0],a[1]*b0[1])) if a[1] and b0[1] else (1,0)
 Z=(1,0);Lb=[];Qb=[];Ab=[]
 for t0 in range(2):
  l=[]
  for c in J:
   q=Z
   for i,x in c.items():q=add(q,mul(ex(pre[t0].get(i,Q(0))),ex(x)))
   l.append(q)
  v=[]
  for c in F:
   q=Z
   for i,x in c.items():q=add(q,mul(l[i],ex(x)))
   v.append(q)
  a={}
  for i in range(102):
   for var in range(9):
    q=ex(pre[t0].get(9*i+var,Q(0)))
    for j,x in Pi[i].items():q=add(q,mul(ex(x),v[9*j+var]))
    if q[1]:a[9*i+var]=q
  Lb.append(l);Qb.append(v);Ab.append(a)
 for a in Ab:
  for c in G:
   q=Z
   for i,x in c.items():q=add(q,mul(a.get(i,Z),ex(x)))
 put('phase.json',{'phase':'independent_shared27_positive_before_signed_replay','positive_prior_bits':maximum,'no_rho_b1_Z2_upper_response':True})
 actual=0
 def ck(q):
  nonlocal actual
  actual=max(actual,abs(q.numerator).bit_length(),q.denominator.bit_length());assert actual<=maximum;return q
 def dot(a,b0):
  q=Q(0)
  for i,x in reversed(list(b0.items())):q=ck(q+ck(a.get(i,Q(0))*x))
  return q
 LJ=[];pull27=[];A=[]
 for t0 in range(2):
  lj=[dot(pre[t0],c) for c in J];LJ.append(lj);v=[dot({j:x for j,x in enumerate(lj) if x},c) for c in F];pull27.append(v);a={}
  for i in range(102):
   for var in range(9):
    q=pre[t0].get(9*i+var,Q(0))
    for j,x in reversed(list(Pi[i].items())):q=ck(q-ck(x*v[9*j+var]))
    assert q==sp(e['alpha_two_by918_sparse_rows'][t0]).get(9*i+var,Q(0))
    if q:a[9*i+var]=q
  A.append(a)
 assert [[str(x) for x in v] for v in LJ]==e['actual_preJ25_two_by25']
 gdots=[[dot(a,c) for c in G] for a in A];assert all(not x for row in gdots for x in row) and [[str(x) for x in row] for row in gdots]==e['all836_G3_literal_dots']
 proof={'all1836_alpha_slots_including_zeros':True,'all50_preJ25_coefficients_match':True,'different_grouping':'Computeactualshared27row(preJ25)*F25once,thenpullbackviaPiY3;producerrepeatedF25/PiYproductsper918slot. Reverseexactdotorder.','all836_original_fullG3_annihilation_dots_exact0':True,'ordinaryslot_orientation_Krell_r7_Krr_ell8':True,'full_source_frame_and_fullPsi_K_history':True,'alpha_two_by918_sparse_rows':[[[i,str(x)] for i,x in sorted(a.items())] for a in A],'actual_shared27_correction_rows':[[str(x) for x in v] for v in pull27],'positive_prior_bits':maximum,'own_actual_bits':actual,'scope':'ActualalphasannihilatethefullknownG3imageinthepeeredsourceframe. Noactualrho/rightsolveapplication/RHSZ2/b1/upper/HF/length.'}
 ep='independent1836-alpha-shared27-and836-original-G3-dots.json';put(ep,proof)
 owned=count([r,e,rows,S,pre,z,Pi,G,sel,J,F,Lb,Qb,Ab,LJ,pull27,A,gdots,proof,pins])+12000;assert owned<60000
 report={'id':STEM,'all_passed':True,'all1836_clearing_slots_and836_G3dots_checked':True,'no_rho_RHS_b1_upper_or_higher_response':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'positive_prior_bits':maximum,'actual_bits':actual,'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--expected-primary-sha',required=True);a=p.parse_args();main(a.expected_primary_sha)
