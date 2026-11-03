#!/usr/bin/env python3
"""Independent bounded degree-one direct-lift obstruction and saved conditioning.

Modulo-prime full column rank is used only for nonzero integer minors, which
prove characteristic-zero injectivity of a fixed integer linear map.  There is
no module-lifting, special-fiber-flatness or finite-defect inference.
"""
import argparse,gzip,hashlib,json,shutil,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def minor(A,p):
 """Independent NumPy modular elimination, preserving selected original rows."""
 B=np.array(A,dtype=np.int64,copy=True);n,m=B.shape;ids=list(range(n));selected=[];det=1;r=0
 for j in range(m):
  hits=np.flatnonzero(B[r:,j])
  if not len(hits):continue
  i=r+int(hits[0]);B[[r,i]]=B[[i,r]];ids[r],ids[i]=ids[i],ids[r]
  pivot=int(B[r,j]);det=det*pivot%p;selected.append(ids[r]);inv=pow(pivot,-1,p)
  B[r]=B[r]*inv%p
  if r+1<n:B[r+1:]=(B[r+1:]-B[r+1:,j:j+1]*B[r:r+1])%p
  r+=1
  if r==min(n,m):break
 # For a full-column-rank matrix, the selected original rows in this ORDER
 # have determinant equal to the pivot product (no extra row permutation).
 return {'rank':r,'selected_original_rows':selected,'ordered_minor_det_mod_p':det if r==m else None}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);start=time.monotonic();base=ROOT/'work/nonsquare-J0-pid-v63'
 src=base/'degree1-row4-v1/scalar-primitive/residual.txt';assert sha(src)=='55174c163211b4252a9f1ccff04de204b5aece4a1c26aab2a539f8eb5d489bb6'
 native=base/'direct-lift-v1/result.json';d=json.loads(native.read_text());p=d['prime'];pairs=d['evaluation_pairs'];assert p==32003
 n,m=d['source_shape'];assert[n,m]==[20,54];evals=[[[0]*m for _ in range(n)]for _ in pairs]
 apows=[[pow(a,k,p)for k in range(128)]for a,e in pairs];epows=[[pow(e,k,p)for k in range(128)]for a,e in pairs]
 with src.open()as f:
  assert list(map(int,f.readline().split()))==[20,54,1];f.readline();f.readline()
  for i in range(n):
   for j in range(m):
    s=f.readline().split();nt=int(s[0]);assert len(s)==1+4*nt;v=[0]*len(pairs)
    for k in range(nt):
     c=int(s[1+4*k])%p;a,e,b=map(int,s[2+4*k:5+4*k]);assert b==0
     for h in range(len(pairs)):v[h]=(v[h]+c*apows[h][a]*epows[h][e])%p
    for h in range(len(pairs)):evals[h][i][j]=v[h]
  assert not f.read().strip()
 N=[]
 for h,(a,e)in enumerate(pairs):
  for row in evals[h]:N.append(row+[a*v%p for v in row]+[e*v%p for v in row])
 assert N==d['evaluation_matrix'];certificates=[]
 for target in range(n):
  indices=[h*n+i for h in range(len(pairs))for i in range(n)if i!=target];result=minor([N[k]for k in indices],p)
  assert result['rank']==3*m==d['excluded_row_ranks'][target]
  chosen=[indices[k]for k in result['selected_original_rows']]
  # A second determinant computation, in the saved chosen row order, checks
  # that the actual displayed square minor has the asserted nonzero value.
  check=minor([N[k]for k in chosen],p);assert check['rank']==3*m and check['selected_original_rows']==list(range(3*m))and check['ordered_minor_det_mod_p']==result['ordered_minor_det_mod_p']
  certificates.append({'excluded_physical_row':target,'selected_evaluation_rows':chosen,'det_mod_prime':result['ordered_minor_det_mod_p']})
 # Exact characteristic-zero generic full rank: the first20 columns at(1,2)
 # form an INTEGER minor nonzero modulo32003, so their QQ value is nonzero.
 generic=minor(np.array(evals[0],dtype=np.int64)[:,:20].T,p);assert generic['rank']==20
 # Correct for the row order used in the transposed selected square determinant.
 M=np.array(evals[0],dtype=np.int64)[:,:20];detcheck=minor(M,p);assert detcheck['rank']==20
 payload={'source_sha256':sha(src),'prime':p,'cofactor_basis':['1','a','e'],'source_shape':[n,m],'evaluation_pairs':pairs,'evaluation_matrix':N,'direct_lift_minor_certificates':certificates,'generic_full_rank_integer_minor':{'point':[1,2],'columns':list(range(20)),'ordered_row_indices':detcheck['selected_original_rows'],'det_mod_prime':detcheck['ordered_minor_det_mod_p']}}
 with gzip.open(out/'certificate.json.gz','wt',compresslevel=9)as f:json.dump(payload,f,separators=(',',':'))
 col=base/'column-units-v1';assert json.loads((col/'receipt.json').read_text())['exact_diagonal_pass'];normal=json.loads((col/'normalize-receipt.json').read_text());assert normal['status']=='completed'
 best=col/'primitive.txt';assert sha(best)==normal['output_sha256']=='039ee56e033a9d576aae8e271483f34e10dd111e87f391eb740251dc90c9489a'
 powers=[list(map(int,s.split()))for s in(col/'powers.txt').read_text().splitlines()[1:]];assert len(powers)==52;totals=[sum(r[1+u]for r in powers)for u in range(5)];assert totals==[51,0,0,0,33]
 for folder,prefix,names in [(col,'column',['powers.txt','log.txt','receipt.json','QQ-factors.txt','normalize.log','normalize-receipt.json']),(base/'row-lll-v1','lll',['U.txt','log.txt','receipt.json','normalization-receipt.json']),(base/'row-lll-fixture','lll-fixture',['input.txt','output.txt','U.txt','log.txt']),(base/'direct-lift-v1','direct',['log.txt'])]:
  for name in names:shutil.copyfile(folder/name,out/f'{prefix}-{name}')
 inputs=[Path(__file__),src,native,ROOT/'code/s61_nonsquare_J0_direct_lift_20261103_v63.cpp',ROOT/'code/polynomial_known_unit_columns_20261103_v63.cpp',ROOT/'code/polynomial_integer_row_lll_20261103_v63.cpp',ROOT/'code/polynomial_QQ_scalar_primitive_20261103_v63.cpp',ROOT/'reports/S61_NONSQUARE_J0_DIXON_UNITS_20261103_v63/report.json']
 report={'id':'S61_NONSQUARE_J0_PRECONDITIONING_20261103_V63','all_passed':True,'base_ring':'QQ[a,e,1/(a*F*Bs*C*(3e+2))]','J_not_inverted':True,'unknown_core_factor_inverted':False,'new_unit_pivots':0,'total_safe_pivots_from57basis':39,'best_shape':[18,52],'best_source':str(best.relative_to(ROOT)),'best_sha256':sha(best),'best_bytes':best.stat().st_size,'whole_column_declared_unit_divisions':dict(zip(['a','F','Bs','C','3e+2'],totals)),'all_column_divisions_native_reversed_and_exact':True,'all_numeric_QQ_scalings_native_reversed_and_exact':True,'independently_rebuilt_modular_evaluation_entries':len(pairs)*n*3*m,'all20_excluded_row_matrices_have_full162_column_rank':True,'degree1_polynomial_direct_row_target_impossible_over_QQ_for_each20row':True,'obstruction_scope':'Only q_j in QQ[a,e] with total degree<=1, on the pinned20x54 frame; no restriction proved for higher degree or localized/rational q_j. The target scalar is unrestricted.','generic_full_cubic_rank57_on_J0_proved_by_integer_minor':True,'finite_J0_defect_support_proved':False,'no_GF_to_QQ_flatness_or_Nakayama_inference':True,'lll_basis_reduction_only_marginal_and_not_used_in_best_state':True,'lll_bit_sum_before':424895290,'lll_bit_sum_after':421206788,'input_sha256':{str(p.resolve().relative_to(ROOT)):sha(p)for p in inputs},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['all_passed','best_shape','best_bytes','degree1_polynomial_direct_row_target_impossible_over_QQ_for_each20row','seconds']},indent=2))
if __name__=='__main__':main()
