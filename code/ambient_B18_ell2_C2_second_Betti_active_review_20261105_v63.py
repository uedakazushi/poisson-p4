#!/usr/bin/env python3
"""Independent source split, fixed126 minor and all364 factor postcheck."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_RANK_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63'
P={
 'primary':f'reports/{A}/report.json','rank':f'reports/{A}/active-d3-rank-and-all-source-factorizations.json',
 'gate':f'reports/{G}/report.json','source':f'reports/{G}/active-five-d3-source-split-and-budget.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'first':'reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63/actual-d2-origin-rank-and-source-factorization.json',
 'first_peer':'reports/AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 'resolution_peer':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
ACTIVE=[0,1,4,5,7];ZERO=[2,3,6]
def blocks(n,vs):return [(n-len(w),w) for q in range(min(n,len(vs))+1) for w in itertools.combinations(vs,q)]
def d0(n,vs,actions):
 rows=blocks(n-1,vs);ix={w:i for i,(p,w) in enumerate(rows)};columns=[]
 for p,w in blocks(n,vs):
  for j in range(14):
   col={}
   if p:col=add(col,{14*ix[w]+i:(-1)**p*c for i,c in actions[8][j].items()})
   for pos,k in enumerate(w):
    off=14*ix[w[:pos]+w[pos+1:]];col=add(col,{off+i:(-1)**(p+pos+1)*c for i,c in actions[k][j].items()})
   columns.append(col)
 return columns
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'INDEPENDENT_FIXED126_MINOR_SOURCE_SPLIT','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert all(D[k]['all_passed'] for k in P if k in ['primary','gate'] or k.endswith('_peer'))
 for name in ['primary','gate']:
  for p,h in D[name]['input_sha256'].items():assert sha(ROOT/p)==h
  folder=A if name=='primary' else G
  for p,h in D[name]['evidence_sha256'].items():assert sha(ROOT/'reports'/folder/p)==h
 actions=[[sp(c) for c in M] for M in D['C2']['C2_coordinate_actions9_sparse_columns']]
 assert all(not c for i in ZERO for c in actions[i])
 active={n:d0(n,ACTIVE,actions) for n in [1,2,3]}
 assert [len(active[n]) for n in [1,2,3]]==[84,224,364]
 assert [sum(map(len,active[n])) for n in [1,2,3]]==[37,193,447]
 S=D['source'];assert S['active_variables']==ACTIVE and S['zero_actions']==ZERO
 for n in [1,2,3]:assert list(map(enc,active[n]))==S[f'active_d{n}_columns']
 assert all(not mv(active[1],c) for c in active[2]) and all(not mv(active[2],c) for c in active[3])
 split={}
 for n in [1,2,3]:
  indices={q:{w:i for i,(p,w) in enumerate(blocks(n-q,ACTIVE))} for q in range(min(n,3)+1)};split[n]=[]
  for p,w in blocks(n,range(8)):
   a=tuple(k for k in w if k in ACTIVE);z=tuple(k for k in w if k in ZERO);sgn=(-1)**sum(i>j for i in a for j in z)
   split[n].extend((z,14*indices[len(z)][a]+j,sgn) for j in range(14))
  assert [[list(z),i,s] for z,i,s in split[n]]==S['signed_zero_wedge_decomposition'][str(n)]
 for n in [2,3]:
  full=d0(n,range(8),actions)
  if n==2:assert full==list(map(sp,D['first']['all_d2_constant_columns']))
  for j,col in enumerate(full):
   z,cj,sign=split[n][j];expected={(z,i):sign*c for i,c in active[n-len(z)][cj].items()} if len(z)<n else {}
   got={}
   for i,c in col.items():
    zz,ii,ss=split[n-1][i];got[(zz,ii)]=ss*c
   assert got==expected
 # No rank search: replay supplied exact independent-column operations.
 E=D['rank'];cols=active[3];basis={};comb={};selected=[];pivseq=[]
 for r in E['fixed_pivot_receipts']:
  j=r['source_column'];v=dict(cols[j]);w={j:Q(1)}
  for pivot,c in r['subtract']:
   c=Q(c);assert v.get(pivot)==c;v=add(v,scale(basis[pivot],-c));w=add(w,scale(comb[pivot],-c))
  pivot=r['pivot'];c=Q(r['divide_by']);assert c and min(v)==pivot and v[pivot]==c and pivot not in basis
  basis[pivot]=scale(v,1/c);comb[pivot]=scale(w,1/c);selected.append(j);pivseq.append(pivot)
 assert len(basis)==126 and selected==E['fixed_original_minor_columns']
 rows=sorted(basis);normalized=[basis[p] for p in rows]
 assert list(map(enc,normalized))==E['normalized_independent_columns']
 assert all(mv(normalized,sp(c))==v for c,v in zip(E['all364_column_factorization'],cols))
 assert len(E['all364_column_factorization'])==364
 rowpos={i:j for j,i in enumerate(rows)}
 minor=[{rowpos[i]:c for i,c in cols[j].items() if i in rowpos} for j in selected]
 assert rows==E['fixed_original_minor_rows'] and list(map(enc,minor))==E['fixed_original_minor_sparse_columns']
 # The triangular normalized minor and triangular source change prove the
 # factored determinant, without expanding a huge rational coefficient.
 assert all(basis[p].get(p)==1 and not any(i<p for i in basis[p]) for p in rows)
 sidx={j:k for k,j in enumerate(selected)}
 assert all(comb[p].get(j)==1/Q(E['fixed_pivot_receipts'][k]['divide_by']) and all(sidx[i]<=k for i in comb[p]) for k,(p,j) in enumerate(zip(pivseq,selected)))
 inversions=sum(pivseq[i]>pivseq[j] for i in range(126) for j in range(i+1,126))
 det=E['fixed_minor_exact_nonzero_determinant_product']
 assert det=={'sign':(-1)**inversions,'nonzero_QQ_factors':[r['divide_by'] for r in E['fixed_pivot_receipts']]}
 assert D['first']['rank']==89 and S['active_r1']==11 and S['active_r2']==56
 fullrank=126+3*56+3*11;beta2=518-89-fullrank;assert fullrank==327 and beta2==102
 assert E['full_rank_d3']==fullrank and E['actual_beta2_C2_over_B']==beta2
 ep=o/'independent-active126-minor-and-zero3-source-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all_actual_active_and_full_constant_differentials_reconstructed':True,'all_literal_signed_zero3_wedge_identities_verified':True,'active_d3_shape':[224,364],'active_d3_nonzeros':447,'all126_fixed_receipts_and_original_minor_verified':True,'all364_original_column_factorizations_verified':True,'fixed_nonzero_determinant_kept_factored':det,'full_d3_rank':fullrank,'actual_beta2_C2_over_B':beta2,'no_new_rank_search_or_physical_coefficients':True,'no_actual_T2_or_H0_length_computed':True})
 report={'id':STEM,'all_passed':True,'actual_active_d3_rank126_full327_beta2_102_accepted':True,'source_zero3_wedge_and_all_fixed126_minor_factors_accepted':True,'new_rank_search_or_source_coefficient_or_T2_length_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
