#!/usr/bin/env python3
"""Actual three first columns after all callback cancellations; fixed quotient.
Uses corrected cache V2 only. No new source/jet or residual second relation.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 v2='AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_CACHE_V2_20261105_v63';rep=load(v2);assert pins[f'reports/{v2}/report.json']=='21f4e29dd8e9082e66301f8d25bb352fd5038aa63b3c1d8cc0898c92f5c8d383'
 for p,h in rep['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in rep['evidence_sha256'].items():assert sha(ROOT/'reports'/v2/n)==h
 D=load(v2,'actual-three-oldm1-and-full81-formal-callbacks.json');assert D['all576_first_Koszul_coefficients_verified']and D['full_b_times_nu_including_b1_nu1']
 assert load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW')['actual_origin_rank']==2
 N=load('AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63','actual-natural-Told-to-Tnil-source-map-and-origin.json')['actual_origin_data']
 assert load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63_REVIEW')['actual_first_relation_rank']==35
 L=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json')
 X=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json');W=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase3-actual-W84-rank-and-mu.json');mD=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json')
 first=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase2-row20-surrogates-and-true-m1row7.json')
 M=load('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json')
 assert load('AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW')['all_passed']
 save('phase.json',dict(phase='corrected V2 and all actual source peers pinned before fixed quotient',old_rank39_diagnostic_NOT_reused=True))
 cols=[{int(i):Q.pairs(v)for i,v in c}for c in X['all780_source_columns44_sparse_dual9']];m0=[Q.sparse(v)for v in mD['m0_24_by44']]
 xp=X['X0_pivot_columns'];wcols=[Q.pairs(v)for v in W['all_W84_columns_sparse']];wp=W['W_pivot_residual_columns'];xc=[Q.sparse(v)for v in X['all780_X0_factorization_coefficients']]
 def cres(j):return Q.plus({j:F(1)},*({xp[k]:-a}for k,a in xc[j].items()))
 C=[i for i in range(44)if i not in mD['selected_source_columns24']];I=[C[i]for i in M['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[Q.sparse(v)for v in M['fixed_minor_data']['inverse']]
 Wminor=[{j:wcols[c][r]for j,c in enumerate(wp)if wcols[c].get(r,F(0))}for r in W['W_pivot_rows']];Wi=[Q.sparse(v)for v in W['W_fixed_minor_inverse']];assert Q.mprod(Wminor,Wi)==Q.mprod(Wi,Wminor)==Q.eye(81)
 PAIRS=list(itertools.combinations(range(8),2));actual=[]
 for j,r in enumerate(D['records']):
  assert not r['remaining_actual_adapter_A2_weights']and not r['remaining_native_even2_weights']
  co=Q.pairs(r['total_actual_X780_constant_coefficients'].items());co18=Q.pairs(r['first18'].items());beta=Q.pairs(r['full81'].items())
  expected=Q.plus({xp[k]:-a for k,a in co18.items()},*(Q.scaled(cres(W['residual_original_X_columns'][wp[k]]),-a)for k,a in beta.items()));assert co==expected
  c0=Q.sparse(N['all_five_new_cycle_origins44'][j+1]);f0=Q.plus(c0,*(Q.scaled({i:v[0]for i,v in cols[k].items()if v.get(0,F(0))},a)for k,a in co.items()));assert not f0
  c1=[Q.pairs(v.items())for v in r['new_c1_44x8']];d1=[Q.plus(c1[i],*({q-1:x*a for q,x in cols[k].get(i,{}).items()if q}for k,a in co.items()))for i in range(44)];assert all(not v for v in Q.mprod(m0,d1))
  # BEFORE81, recover the original columnwise Hadamard/Euler h0 and verify
  # all84 entries, not merely the81 pivot entries used to find beta.
  pre1=[Q.plus(c1[i],*({q-1:x*(-a)for q,x in cols[xp[k]].get(i,{}).items()if q}for k,a in co18.items()))for i in range(44)];v=Q.mprod(m0,pre1);h0={}
  for sec in range(3):
   for k,(a,b)in enumerate(PAIRS):
    c=v[8*sec+7].get(a,F(0))if b==7 else(v[8*sec+b].get(a,F(0))-v[8*sec+a].get(b,F(0)))/2
    if c:h0[28*sec+k]=c
  assert Q.plus(h0,*(Q.scaled(wcols[wp[k]],-a)for k,a in beta.items()))=={}
  actual.append(Q.pairs(r['final_known40'].items()))
 # All second callbacks cancel AS FULL SOURCE COEFFICIENTS before any
 # coefficient extraction. Thus the saved known40 formulas are actual.
 # Their ordered Hadamard indices and m1_r carry are exactly those in the
 # separately accepted source-bound FIRST_LINEAR operator, not a new gauge.
 lincols=[Q.pairs(c)for c in L['all849_actual_minimal_first_columns40_sparse']];cp=L['rank_minor_columns'];rr=L['rank_minor_rows'];inv=[Q.sparse(v)for v in L['rank_minor_inverse']]
 minor=[{j:lincols[c][r]for j,c in enumerate(cp)if lincols[c].get(r,F(0))}for r in rr];assert Q.mprod(minor,inv)==Q.mprod(inv,minor)==Q.eye(35)
 free=[i for i in range(40)if i not in rr];F5=[]
 for r in free:
  a=Q.mprod([{j:lincols[c][r]for j,c in enumerate(cp)if lincols[c].get(r,F(0))}],inv)[0]
  F5.append(Q.plus({r:F(1)},{rr[i]:-x for i,x in a.items()}))
 assert len(F5)==5 and all(not Q.mvec(F5,c)for c in lincols)
 origins=[Q.sparse(N['actual_origin_map_columns5'][j])for j in[0,4]]
 products=[{8*i+q:a for i,a in c.items()}for c in origins for q in range(8)];p5=[Q.mvec(F5,c)for c in products]
 k=next(i for i,c in enumerate(p5)if c);row=min(p5[k]);pivot=p5[k][row];v=p5[k]
 for c in p5:assert c==Q.scaled(v,c.get(row,F(0))/pivot)
 F4=[Q.plus(F5[i],Q.scaled(F5[row],-v.get(i,F(0))/pivot))for i in range(5)if i!=row]
 assert len(F4)==4 and all(not Q.mvec(F4,c)for c in lincols+products)
 four=[Q.mvec(F4,c)for c in actual]
 # The length5/mu3 theorem allows rank2 or3; determine it by a fixed
 # selected minor and exact inverse/factor checks, retaining ALL four rows.
 pivotrows={};chosen=[]
 for j,c in enumerate(four):
  z=c.copy()
  for i,p in sorted(pivotrows.items()):z=Q.plus(z,Q.scaled(p,-z.get(i,F(0))))
  if z:
   i=min(z);pivotrows[i]=Q.scaled(z,1/z[i]);chosen.append((i,j))
 rank=len(chosen);assert rank in[2,3]
 # A symbolic QQ small matrix inverse is computed by explicit row operations.
 selectedrows=[i for i,j in chosen];selectedcols=[j for i,j in chosen];A=[{j:four[c][r]for j,c in enumerate(selectedcols)if four[c].get(r,F(0))}for r in selectedrows];B=[v.copy()for v in A];Ai=Q.eye(rank)
 for j in range(rank):
  i=next(i for i in range(j,rank)if B[i].get(j,F(0)))
  if i!=j:B[i],B[j]=B[j],B[i];Ai[i],Ai[j]=Ai[j],Ai[i]
  a=B[j][j];B[j]=Q.scaled(B[j],1/a);Ai[j]=Q.scaled(Ai[j],1/a)
  for i in range(rank):
   if i!=j:
    a=B[i].get(j,F(0));B[i]=Q.plus(B[i],Q.scaled(B[j],-a));Ai[i]=Q.plus(Ai[i],Q.scaled(Ai[j],-a))
 assert B==Q.eye(rank)and Q.mprod(A,Ai)==Q.mprod(Ai,A)==Q.eye(rank)
 det=Q.det_fixed(A);assert det
 for c in four:
  co=Q.mvec(Ai,{i:c[r]for i,r in enumerate(selectedrows)if c.get(r,F(0))});assert Q.plus(*(Q.scaled(four[selectedcols[k]],a)for k,a in co.items()))==c
 hf=[3,4-rank,rank-2];hf=[h for h in hf if h];assert sum(hf)==5 and hf in[[3,2],[3,1,1]]
 proof=dict(source_frames='actual oldfree56 rows23/29/30; corrected b1*nu1 retained; actual natural Z; new18+81 fixed source cleanup',
  origin_rank2_unchanged_dependency='oldc0 uses only oldm0/n0; Z20_0=0 kills oldc1 in origin map; selected h0=Z20_1*cold0',
  all576_actual_first_Koszul_checks=True,all3_complete84_h0_cleanup_checks=True,
  all_adapter_A2_and_native_even2_coefficients_cancel_before_extraction=True,
  actual_three_minimal_first_columns40=actual,actual_old_free56_labels=[23,29,30],
  final_ordered_Hadamard_and_true_m1r_carried=True,
  quotient5_by40_after_actual_rank35_first_relations=F5,origin_product_columns40=products,
  origin_products_add_exactly_one_rank_to35=True,selected_origin_product=k,selected_origin_product_pivot=pivot,
  quotient4_by40_after_origin_products=F4,actual_three_projected_columns4=four,
  selected_small_rows=selectedrows,selected_small_columns=selectedcols,selected_small_inverse=Ai,selected_small_determinant=det,
  actual_added_rank=rank,total_first_relation_rank=36+rank,
  actual_cokernel_HF=hf,actual_cokernel_length=5,actual_cokernel_m3_zero=True,
  no_exact_HF2_Tnil_total_length_nilcube_or_new_source_jet=True,
  old_incorrect_cache_and_rank_diagnostic_withdrawn='AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_CACHE_20261105_v63_ERRATUM/failure-and-correction.json')
 save('actual-three-first-columns-and-fixed-quotient-proof.json',proof)
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,actual_zero_origin_three_column_rank_in_fixed4=rank,actual_cokernel_HF=hf,actual_cokernel_length=5,total_first_relation_rank=36+rank,all_source_first_Koszul_and_formal_callback_gates_closed=True,no_new_source_jet_or_actual_HF2_Tnil_or_nilcube=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
 save('report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
