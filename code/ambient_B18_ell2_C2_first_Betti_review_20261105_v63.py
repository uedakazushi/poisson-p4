#!/usr/bin/env python3
"""Independent d2(0) construction and fixed rank89/minor verification."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import add,scale,sp,enc,mv,eye
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_C2_FIRST_BETTI_SOURCE_PROFILE_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','certificate':f'reports/{A}/actual-d2-origin-rank-and-source-factorization.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'resolution_peer':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'ACTUAL_CONSTANT_PERIODIC_KOSZUL_AND_FIXED_MINOR','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert all(D[k]['all_passed'] for k in ['primary','C2_peer','resolution_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 A9=[[sp(p) for p in M] for M in D['C2']['C2_coordinate_actions9_sparse_columns']];E=A9[8];d1=[scale(v,-1) for M in [E]+A9[:8] for v in M]
 blocks1=[(1,())]+[(0,(i,)) for i in range(8)];b1index={p:i for i,p in enumerate(blocks1)}
 blocks2=[(2-q,w) for q in range(3) for w in itertools.combinations(range(8),q)]
 cols=[];labels=[]
 for p,w in blocks2:
  for j in range(14):
   v={}
   if p:
    sign=(-1)**p;off=14*b1index[p-1,w]
    v=add(v,{off+i:sign*c for i,c in E[j].items()})
   for k,var in enumerate(w):
    target=(p,w[:k]+w[k+1:]);off=14*b1index[target];sign=(-1)**(p+k+1)
    v=add(v,{off+i:sign*c for i,c in A9[var][j].items()})
   cols.append(v);labels.append([p,list(w),j])
 C=D['certificate'];assert len(cols)==518 and sum(map(len,cols))==304
 assert list(map(enc,cols))==C['all_d2_constant_columns'] and labels==C['column_labels']
 assert [[p,list(w)] for p,w in blocks1]==C['blocks_F1'] and [[p,list(w)] for p,w in blocks2]==C['blocks_F2']
 assert all(not mv(d1,p) for p in cols)
 piv={};receipts=C['fixed_pivot_receipts'];chosen=[];sequence=[];divides=[]
 for row in receipts:
  v=cols[row['source_column']]
  for k,c in row['subtract']:v=add(v,scale(piv[k],-Q(c)))
  k=row['pivot'];d=Q(row['divide_by']);assert k not in piv and v.get(k)==d and all(i>=k for i in v)
  v=scale(v,1/d);assert enc(v)==row['normalized'];piv[k]=v;chosen.append(row['source_column']);sequence.append(k);divides.append(d)
 assert len(piv)==C['rank']==89
 rows=sorted(piv);basis=[piv[k] for k in rows];factors=list(map(sp,C['all518_column_factorization_in_normalized_basis']))
 assert list(map(enc,basis))==C['normalized_source_basis_columns'] and len(factors)==518
 assert all(mv(basis,p)==v for p,v in zip(factors,cols))
 assert rows==C['fixed_original_minor_rows'] and chosen==C['fixed_original_minor_columns']
 literal=[[str(cols[j].get(i,Q(0))) for j in chosen] for i in rows];assert literal==C['fixed_original_minor']
 # Independent determinant factorization: normalized row-pivot minor is
 # unit lower triangular; selected source factor columns are triangular
 # in receipt order. No determinant or inverse solve is repeated.
 assert all(piv[j].get(i,Q(0))==(Q(1) if i==j else Q(0)) for j in rows for i in rows if i<=j)
 rowposition={p:i for i,p in enumerate(rows)}
 for j,col in enumerate(chosen):
  assert all(factors[col].get(rowposition[sequence[i]],Q(0))==(divides[j] if i==j else 0) for i in range(j,89))
 inversions=sum(sequence[i]>sequence[j] for i in range(89) for j in range(i+1,89));det=Q((-1)**inversions)
 for d in divides:det*=d
 assert det and str(det)==C['fixed_original_minor_determinant']
 assert D['resolution_peer']['actual_C2_Loewy_dimensions'][1]==11
 beta1=126-11-89;assert beta1==C['actual_beta1_C2_over_B']==26
 ep=o/'independent-d2-origin304-fixed89-minor-and-Betti26.json'
 dump(o,ep.name,{'inputs':inputs,'actual_periodic_Koszul_F1_F2_order_and_signs_reconstructed':True,
  'all304_nonzero_entries518columns_verified':True,'all518_literal_d1d2_products_zero':True,
  'all89_fixed_rank_receipts_and518_source_factorizations_verified':True,
  'fixed_original_minor_rows':rows,'fixed_original_minor_columns':chosen,'nonzero_minor_determinant':str(det),
  'determinant_proof':'unit lower triangular normalized pivot minor and receipt-order triangular factor minor; exact pivot permutation and all89 diagonal factors checked',
  'actual_d1_origin_rank':11,'actual_d2_origin_rank':89,'actual_beta1':beta1,
  'Tor_proof':'complete source resolution; Tor1_dimension=ker(d1_tensor_k)-im(d2_tensor_k)=126-11-89',
  'actual_V_Hadamard_source_coefficients_not_computed':True,'no_new_elimination_source_solve_T2_or_H0':True})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'actual_d2_origin_rank':89,'actual_beta1_C2_over_B':26,'all_source_entries_minor_and518_factorizations_accepted':True,
  'actual_V_coefficients_T2_H0_not_computed':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
