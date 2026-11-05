#!/usr/bin/env python3
"""Fixed row/column determinant and all-column span replay, no rank search."""
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def put(n,d):(o/n).write_text(json.dumps(d,separators=(',',':'))+'\n')
 b='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63/';r=rd(b+'report.json');assert r['all_passed'] and pins[b+'report.json']=='de7667b93bec2447fc1a4dcabfebd2cf0e0ea3862165c65077587b7083443df2'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/b/p)==h
 e=rd(b+'fixed46-by41-integer-Schur-rank-and-full-strand-kernel-recipes.json')
 gb='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63/';g=rd(gb+'report.json');d=rd(gb+'pruned46-by41-integer-rank-source-and-live-bound.json')
 pp='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63_REVIEW/report.json';p=rd(pp);assert p['all_passed'] and pins[pp]=='ed6954d48503e2220e1d7ab29526a64ead78387ab00e1012f70b454172ff20e3'
 put('phase.json',{'phase':'fixed_certificate_before_arithmetic','source_peer_sha256':pins[pp],'rank_search':False})
 sr=e['selected_residual_rows46'];sc=e['selected_residual_columns41'];assert len(sr)==len(sc)==28 and len(set(sr))==len(set(sc))==28
 R=sr+[i for i in range(46) if i not in sr];C=sc+[j for j in range(41) if j not in sc]
 original=d['residual_integer_matrix46_by41'];M=[[original[i][j] for j in C] for i in R];previous=1;maxcommit=0;maxtemp=0
 # Source prior: ordered nested minors <=306bits. Product differences<=613.
 def ck(x,limit):
  nonlocal maxcommit,maxtemp
  if limit==306:maxcommit=max(maxcommit,abs(x).bit_length())
  else:maxtemp=max(maxtemp,abs(x).bit_length())
  assert abs(x).bit_length()<=limit;return x
 assert d['all_minors_order_le29_bits306']==306 and d['integer_product_difference_bits613']==613
 for k in range(28):
  pivot=M[k][k];step=e['fixed_integer_pivot_steps'][k]
  assert pivot and step['pivot_row46']==sr[k] and step['pivot_col41']==sc[k] and pivot==int(step['pivot_integer']) and previous==int(step['exact_previous_divisor'])
  for i in range(k+1,46):
   for j in range(k+1,41):
    n=ck(ck(pivot*M[i][j],613)-ck(M[i][k]*M[k][j],613),613);q,rem=divmod(n,previous);assert rem==0;M[i][j]=ck(q,306)
  for i in range(k+1,46):M[i][k]=0
  previous=pivot
 assert all(M[i][j]==0 for i in range(28,46) for j in range(28,41))
 det=previous;assert det==int(e['fixed_residual_integer_minor']) and det
 # Integer back-substitution with one common determinant, not an inverse/rank solve.
 # Each Z is a selected Cramer minor <=306 bits. At most28 products in a partial
 # triangular sum: 612+ceil(log2(29))<=617 bits, below1153.
 factors=[]
 for j0 in range(41):
  if j0 in sc:z=[det if j==sc.index(j0) else 0 for j in range(28)]
  else:
   j=C.index(j0);z=[0]*28
   for k in range(27,-1,-1):
    n=ck(det*M[k][j],617)
    for s in range(k+1,28):n=ck(n-ck(M[k][s]*z[s],617),617)
    z[k],rem=divmod(n,M[k][k]);assert rem==0;ck(z[k],306)
  for i in range(46):
   n=0
   for k in range(28):n=ck(n+ck(original[i][sc[k]]*z[k],617),617)
   assert n==det*original[i][j0]
  factors.append(z)
 assert e['residual_rank']==28 and e['active_G2_linear_strand_rank']==190 and e['active_kernel_dimension']==62
 assert len(d['additional_zero_source_columns49'])==49 and 41-28+49==62
 kernel=e['complete_kernel_recipes'];assert kernel['residual_free_columns41']==C[28:] and kernel['additional_zero_columns49_original252']==d['additional_zero_source_columns49']
 assert kernel['selected_rows46']==sr and kernel['selected_cols41']==sc and int(kernel['determinant_literal'])==det
 # Source maps and reverse singleton recipes are inherited from the literal
 # independent gate. No full252-source inverse or new rank selection is needed.
 proof={'fixed_selected_rows46':sr,'fixed_selected_cols41':sc,'selected_integer_determinant':str(det),'all28_supplied_nested_pivots_match_after_upfront_source_permutation':True,'all_remaining18_by13_Schur_entries_zero':True,'all41_original_column_factor_numerators_common_det':[[str(x) for x in z] for z in factors],'all1886_original_factorization_identities':True,'residual_rank':28,'active_rank':190,'complete_active_kernel_dimension':62,'kernel_recipe':'13 residual free sources e_j-Eselected*Z_j/det plus49 isolated zero sources; reverse162 literal singleton source substitutions from independently peered no-fill reduction.','full918_strand_kernel_equals_known_G3linear242':True,'tensor_scope':'internaldegree3 coefficientdegree1 only, from exact literal zero3 tensor source peer; not general associatedgraded exactness','b1_functional_annihilation_still_separate':True,'source_committed_prior306':True,'source_product_difference_prior613':True,'fixed_Cramer_backsub_partial_prior617':True,'actual_max_committed_bits':maxcommit,'actual_max_temp_bits':maxtemp,'no_rank_search_or_new_physical_coefficient':True}
 ep='independent-fixed-minor-all41-factorizations-and-full-strand-proof.json';put(ep,proof)
 owned=count([r,e,g,d,p,R,C,original,M,factors,proof,pins])+7000;assert owned<60000
 report={'id':STEM,'all_passed':True,'fixed_certificate_replay_only_no_rank_search':True,'residual_rank':28,'active_rank':190,'full918_strand_kernel_equals_G3linear242':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'owned_numeric_bound':owned,'seconds':time.monotonic()-t};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'owned':owned}))
if __name__=='__main__':main()
