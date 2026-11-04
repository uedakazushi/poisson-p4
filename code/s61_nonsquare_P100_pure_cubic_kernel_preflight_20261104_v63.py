#!/usr/bin/env python3
"""Full-source maps and bounds for10 pure cubic source-kernel vectors. NO values."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
from fractions import Fraction as Q
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
from s61_nonsquare_P100_pure98_source_preflight_20261104_v63 import ZERO,addbound,mb,width,bdict
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE_CUBIC_KERNEL_PREFLIGHT_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'
 bp=base/'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json';matrix=bp.with_name('remaining.txt.gz');sp=base/'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json'
 ap=base/'S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';hp=base/'S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63/certificate.json';rawp=base/'S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63/preflight.json'
 sigp=base/'S61_NONSQUARE_P100_H1_RIGHT_INVERSE_20261104_v63/native.json';etap=base/'S61_NONSQUARE_P100_H1_DUAL_20261104_v63/native.json'
 ci=base/'S61_NONSQUARE_P100_PURE98_SOURCE_PREFLIGHT_20261104_v63/cubic-coefficient-input.json';pivotcache=base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/cubic504.json'
 literalinput=base/'S61_NONSQUARE_P100_PURE21_SOURCE_PREFLIGHT_20261104_v63/source-input.json';g18cache=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63/physical-cubic84.json'
 paths=[Path(__file__),bp,matrix,sp,ap,hp,rawp,sigp,etap,ci,pivotcache,literalinput,g18cache]
 b,s,ad,h,raw,sig,eta,oldci,pvals,lit,g18=map(read,[bp,sp,ap,hp,rawp,sigp,etap,ci,pivotcache,literalinput,g18cache])
 P=[1,6,8,9,11,14];I=[i for i in range(16)if i not in P];physical=[0,1,2,3,4,6,7];W=3
 sigma=sig['linear_lift_coefficients7x16'];assert sig['all7_fixed_source_linear_lifts_verified']and eta['all16_actual_rawH_annihilations_exact']and eta['v3_is_one']
 assert sorted({p for row in sigma for p,v in enumerate(row)if v['numerator']})==P and all(not v['numerator']for v in sigma[W])
 mons=s['normal_monomials'];oldmons=oldci['old_physical_monomials84'];assert len(oldmons)==84
 monindex={tuple(m):j for j,m in enumerate(mons)};purecols=[monindex[tuple(m)]for m in oldmons];hcols=[]
 for v in physical:
  m=[0]*8;m[5]=1;m[v]+=1;hcols.append(monindex[tuple(m)])
 selected=set(range(41,57));rows={};hashraw=hashlib.sha256()
 with gzip.open(matrix,'rb')as f:
  header=f.readline();hashraw.update(header);assert header.strip()==b'57 91 57'
  for row in range(57):
   rr=[]
   for col in range(91):
    ln=f.readline();assert ln;hashraw.update(ln)
    if row in selected:
     t=ln.decode().split();n=int(t[0]);assert len(t)==1+4*n;rr.append([[list(map(int,t[2+4*k:5+4*k])),str(Q(t[1+4*k]))]for k in range(n)])
   if row in selected:rows[row]=rr
  assert not f.read()
 assert hashraw.hexdigest()==b['remaining_matrix_sha256']
 canon=lambda p:sorted((tuple(m),str(Q(c)))for m,c in p if Q(c))
 records=[]
 for i in range(16):
  pos=41+i;eq=h['all112_cubic_source319_definitions'][7*i][2];rr=rows[pos]
  assert ad['basis57_source_definitions'][pos]=={'kind':'genuine_cubic','original_equation':eq}
  assert all(canon(rr[c])==canon(raw['raw16x7'][i][j])for j,c in enumerate(hcols))
  records.append({'genuine_index':i,'actual135_equation':eq,'basis57_position':pos,'original_normalized63_row':b['remaining_original_row_indices'][pos],
                  'oracle_to_saved_row_unit':ad['basis57_encoded_source_units'][pos],'pure84':[rr[c]for c in purecols],'Hlinear7':[rr[c]for c in hcols]})
 literal=[0,3,7];missing=[i for i in I if i not in literal]
 assert all(not p for i in literal for p in records[i]['Hlinear7'])and missing==[2,4,5,10,12,13,15]
 assert [r['genuine_index']for r in oldci['sources6']]==P and len(pvals['values6x84'])==6
 for i in literal:
  lp=next(r for r in lit['three_full_sources']if r['genuine_index']==i)
  for mon,p in zip(lit['physical_pure_cubic_monomials84'],lp['full_pure84_coefficients']):assert canon(p)==canon(records[i]['pure84'][oldmons.index(mon)])
 # No coefficient value is evaluated: propagate only degrees and known denominator powers.
 factors=[(Q(1),Q(0)),(Q(1),Q(20)),(Q(1),Q(60)),(Q(203),Q(2700)),(Q(1),Q(12))];dcache={}
 def denvec(poly):
  key=tuple((m[0],c)for m,c in poly)
  if key in dcache:return dcache[key]
  rem={m[0]:Q(c)for m,c in poly};v=[]
  for a,b0 in factors:
   exponent=0
   while rem and max(rem)>=1:
    r=rem.copy();q={}
    while r and max(r)>=1:
     n=max(r);c=r[n]/a;q[n-1]=q.get(n-1,Q())+c
     for j,t in[(n,c*a),(n-1,c*b0)]:
      r[j]=r.get(j,Q())-t
      if not r[j]:del r[j]
    if r:break
    rem=q;exponent+=1
   v.append(exponent)
  assert rem and max(rem)==0;dcache[key]=tuple(v);return tuple(v)
 def rfbound(r):return ZERO if not r['numerator']else(max(m[0]for m,c in r['numerator']),denvec(r['denominator']))
 def corebound(p):
  if not p:return ZERO
  clear=max(2*m[0]+m[1]for m,c in p);return clear+max(2*m[2]-m[0]for m,c in p),(clear,0,0,0,0)
 sb=[[rfbound(v)for v in row]for row in sigma];pi={i:[rfbound(v)for v in row]for i,row in zip(P,pvals['values6x84'])}
 alpha={};kernelbounds={};recipes=[]
 for i in I:
  hb=[corebound(p)for p in records[i]['Hlinear7']];alpha[i]={p:addbound([mb(hb[j],sb[j][p])for j in range(7)])for p in P}
  kb=[addbound([corebound(records[i]['pure84'][m])]+[mb(alpha[i][p],pi[p][m])for p in P])for m in range(84)];kernelbounds[i]=kb
  recipes.append({'kernel_label':i,'identity_component_at_nonpivot':1,'pivot_coefficient_recipes':[{'pivot_genuine':p,'coefficient':'-sum_j H_i[j]*sigma[j][p]','term_sources':[[i,j,j,p]for j in range(7)if records[i]['Hlinear7'][j]and sigma[j][p]['numerator']]}for p in P],
                  'H_annihilation_proof':'H_i - sum_j H_i[j]*(e_j-eta_j*e_W)=(H_i.eta)*e_W=0; rawH.eta=0 is accepted.','pure_coefficient_values_NOT_evaluated':True})
 perm=P+I;inversions=sum(perm[i]>perm[j]for i in range(16)for j in range(i+1,16));sign=(-1)**inversions
 widths={i:sum(width(b)for b in bs)for i,bs in kernelbounds.items()}
 counts=lambda ii:sum(len(p)for i in ii for p in records[i]['pure84'])
 cache_mons=lit['physical_pure_cubic_monomials84'];g18_order=[cache_mons.index(m)for m in oldmons]
 source={'physical_light_order':physical,'physical_pure84_monomials':oldmons,'full16_normalized_source_polynomials':records,'formal_pivot_genuine_indices6':P,'formal_kernel_nonpivot_indices10':I,
         'ten_kernel_cofactor_recipes':recipes,'six_pivot_curve_value_cache':str(pivotcache.relative_to(ROOT)),'pivot_cache_mon_order_is_this_order':True,
         'G18_curve_value_cache':str(g18cache.relative_to(ROOT)),'G18_cache_to_this_order':g18_order,'seven_missing_genuine_indices':missing,'all112_Hlinear_source_coefficients_match_full57':True}
 (out/'source-input.json').write_text(json.dumps(source,separators=(',',':'))+'\n')
 cert={'source_order':{'pivots6':P,'nonpivots10':I,'permutation16':perm,'permutation_sign':sign},
       'kernel_formula':'K_i=C_i-sum_p alpha_i,p*C_p; alpha_i,p=sum_j H_i[j]*sigma_j,p. Sigma is ORIGINAL16sourceframe.',
       'source_change_determinant':f'{sign}: pivot-first permutation followed by a unit lower block-triangular16x16 matrix; no unknown determinant.',
       'pivot_basis_proof':'For the six physical columns j!=W, accepted sigma_j*H=e_j-eta_j e_W gives a square6x6 leftinverse to the selected H_P submatrix. Hence selected6rows are a basis of ker eta. Kernel vectors K_i with their nonpivot identity block form the WHOLE formal H-source kernel.',
       'pure_image_scope':'Ten formal source-kernel basis vectors produce actual pure cubics; their images need not be independent modulo quadratic products. No full57 inference.',
       'actual_delta7_frame_transport':'For any original cofactor vector k, new6=k6/delta, new1=k1+k6, otherssame, delta=-64000*C^2 is a wholeP100unit. This is separate from the original Sigma frame.',
       'coverage':{'pivot504_curve_values_cached':True,'G18_84_curve_values_cached':True,'literal_G39_G68_source_polynomials_cached_but_not_curve_values':True,'remaining7_pure84_source_polynomials_already_materialized':True,'remaining7_actual135':[records[i]['actual135_equation']for i in missing],'new_oracle_roots_queries_ancestor_nodes':0,'direct_source_missing588_core_terms':counts(missing),'all16_pure1344_core_terms':counts(range(16))},
       'bounds':{'ten_kernel84_output_RF_terms_upper_by_genuine':widths,'all10_retained_output_terms_upper':sum(widths.values()),'max_single_kernel84_output_terms_upper':max(widths.values()),'pivot504_cached_RF_terms':sum(len(v['numerator'])+len(v['denominator'])for row in pvals['values6x84']for v in row),'sigma7x16_cached_RF_terms':sum(len(v['numerator'])+len(v['denominator'])for row in sigma for v in row),'source_parent_array_has_all1344_but_curve_job_can_stream_one_nonpivot':True},
       'potential_product_count70_is_NOT_rank_or_completion':True,'15mixed_Koszul_and4quadH_overlap_vectors_NOT_constructed':True,'no_curve_coefficient_values_or_kernel84_values_evaluated':True,'no_coefficient_job_launched':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 report={'id':ID,'source_only':True,'full112_Hlinear_source_coefficients_verified':True,'formal_kernel_basis_size':10,'source_change_determinant':sign,'seven_remaining_actual135':cert['coverage']['remaining7_actual135'],'new_oracle_queries':0,'coefficient_values_NOT_evaluated':True,'bounds':cert['bounds'],'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'source_input_sha256':sha(out/'source-input.json'),'certificate_sha256':sha(out/'certificate.json'),'seconds':time.monotonic()-start,'frozen':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
