#!/usr/bin/env python3
"""Provided four point rows: original-source memberships and uniqueness only."""
from pathlib import Path
from fractions import Fraction as Q
import json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='c9c4f1863f4dcbf0b941a4cc3cf0049ddbc227bf5e7ae069eb4857c1359fc0e8'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 e=rd(b+'actual-four-specific-B99-point-dual-rows-STOP-only384.json')
 sr='reports/AMBIENT_B18_ELL2_T2_UPPER617_ALPHA_ROWS_GUARDED_INPUT_GATE_20261105_v63_REVIEW/report.json';assert rd(sr)['all_passed'] and pins[sr]=='07b54a35f003463f224f9c48da91cf8ce5895d37ae4d93a41658e8af7fcedd5f'
 mu='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json';z=fld(mu,['normalized_basis_rows','normalized_basis_original_selected_source_representations','actual_selected_Y0_columns102']);N={i:sp(v) for i,v in z['normalized_basis_rows']};R={i:sp(v) for i,v in z['normalized_basis_original_selected_source_representations']};Y=list(map(sp,z['actual_selected_Y0_columns102']))
 K=fld('reports/AMBIENT_B18_ELL2_T2_UPPER617_K99_20261105_v63/actual396-Kfull-P5-and-native23-existing-cache-values.json',['K99_two_channels_rr_rell'])['K99_two_channels_rr_rell'];K=[[[Q(x) for x in row] for row in ch] for ch in K]
 rows=e['actual_four_rows_channel_pair'];assert [(x['channel'],x['pair_rr_rell']) for x in rows]==[(0,0),(0,1),(1,0),(1,1)];S=[sp(v['point_row102']) for v in rows]
 assert [len(v) for v in S]==[3,16,16,18] and all(not(set(v)&{1,36,38}) for v in S)
 assert len(N)==len(R)==len(Y)==len(K)==99 and set(range(102))-set(N)=={1,36,38} and all(v.get(i)==1 and all(j>=i for j in v) for i,v in N.items())
 maximum=0
 def bound(terms):
  nonlocal maximum
  if not terms:return (1,0)
  D=math.lcm(*(x.denominator*y.denominator for x,y in terms));U=sum(abs(x.numerator*y.numerator)*(D//(x.denominator*y.denominator)) for x,y in terms);maximum=max(maximum,D.bit_length(),U.bit_length());assert maximum<=384;return D,U
 for rec,s in zip(rows,S):
  a,b0=rec['channel'],rec['pair_rr_rell']
  for c in Y:bound([(x,c.get(i,Q(0))) for i,x in s.items()])
  for i,c in N.items():
   bound([(x,c.get(j,Q(0))) for j,x in s.items()]);bound([(K[j][a][b0],x) for j,x in R[i].items()])
 put('phase.json',{'phase':'independent_positive_membership_bound_before_exact_dots','positive_prior_bits':maximum,'no_row_solve_or_alpha_query':True})
 actual=0
 def dotterms(terms):
  nonlocal actual
  q=Q(0)
  for x,y in reversed(terms):
   v=x*y;q=q+v;actual=max(actual,abs(v.numerator).bit_length(),v.denominator.bit_length(),abs(q.numerator).bit_length(),q.denominator.bit_length());assert actual<=maximum
  return q
 identities=[]
 for rec,s in zip(rows,S):
  a,b0=rec['channel'],rec['pair_rr_rell']
  for j,c in enumerate(Y):assert dotterms([(x,c.get(i,Q(0))) for i,x in s.items()])==K[j][a][b0]
  for i,c in N.items():assert dotterms([(x,c.get(j,Q(0))) for j,x in s.items()])==dotterms([(K[j][a][b0],x) for j,x in R[i].items()])
  identities.append({'channel':a,'pair_rr_rell':b0,'all99_original_source_memberships':True,'all99_original_normalized_memberships':True,'free3_exact_zero':True})
 proof={'four_actual_point_rows':rows,'all396_original_source_and396_normalized_memberships':identities,'row_uniqueness':'Thenormalizedpivot99submatrixistriangularwithdiagonal1andfreepointaxes1,36,38complementit. ThusrowY=K, rowfree3=0hasatmostonesolution. Providedrowsverifyexistence;noinverse/solvewasrerun.','source_F2_pointframe_and_order102_unchanged':True,'positive_exact_membership_prior_bits':maximum,'own_membership_actual_bits':actual,'no_alpha918_fixed25clearing_G3dots_rho_b1_Z2_upper_or_higher_values':True}
 ep='independent-four-point-rows-original792-memberships-and-uniqueness.json';put(ep,proof)
 owned=count([r,e,z,N,R,Y,K,rows,S,identities,proof,pins])+8000;assert owned<60000
 report={'id':STEM,'all_passed':True,'all_four_provided_rows_and792_memberships_checked':True,'no_new_point_solve_or_alpha_values':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(out/ep)},'positive_prior_bits':maximum,'actual_membership_bits':actual,'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned}))
if __name__=='__main__':main()
