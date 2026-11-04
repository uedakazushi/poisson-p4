#!/usr/bin/env python3
"""Independent integer-convolution check of norm16 source/annihilator evidence.

No author helper/worker import, determinant, solver, point, prime, or GB.
"""
import argparse,gzip,hashlib,json,math,resource,signal,time
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63_REVIEW'
ONE=({0:1},1)
METRICS={'products_checked':0,'integer_coefficient_products':0,'maximum_checked_sum_support':0}
def read(p):
 with (gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def poly(data):
 """Clear only numeric QQ denominators; pack a,e with six e-bits."""
 coefficients=[];den=1;seen=set()
 for m,s in data:
  assert len(m)==3 and m[2]==0 and min(m)>=0 and m[1]<64
  key=(m[0]<<6)+m[1];assert key not in seen;seen.add(key)
  q=Fraction(s);assert q;den=math.lcm(den,q.denominator);coefficients.append((key,q))
 return ({key:q.numerator*(den//q.denominator)for key,q in coefficients},den)
def shift(p,a=0,e=0):
 assert all((m&63)+e<64 for m in p[0]);return ({m+(a<<6)+e:c for m,c in p[0].items()},p[1])
def scale(p,value):
 q=Fraction(value);return ({m:c*q.numerator for m,c in p[0].items()},p[1]*q.denominator)
def check(terms,rhs):
 """Sum exact products after one integer numeric denominator clearing."""
 den=rhs[1]
 for a,b in terms:
  assert max((m&63 for m in a[0]),default=0)+max((m&63 for m in b[0]),default=0)<64
  den=math.lcm(den,a[1]*b[1])
 acc={m:-c*(den//rhs[1])for m,c in rhs[0].items()}
 for a,b in terms:
  factor=den//(a[1]*b[1]);METRICS['products_checked']+=1
  METRICS['integer_coefficient_products']+=len(a[0])*len(b[0])
  for i,c in a[0].items():
   c*=factor
   for j,d in b[0].items():acc[i+j]=acc.get(i+j,0)+c*d
 METRICS['maximum_checked_sum_support']=max(METRICS['maximum_checked_sum_support'],len(acc))
 assert not any(acc.values()), 'Nonzero exact coefficient in asserted identity'
def equal(a,b):check([(a,ONE)],b)
def encode(p):
 return [[[m>>6,m&63,0],str(Fraction(c,p[1]))]for m,c in sorted(p[0].items(),reverse=True)]
def stats(p):
 return {'terms':len(p[0]),'a_degree':max((m>>6 for m in p[0]),default=-1),'e_degree':max((m&63 for m in p[0]),default=-1),'minimum_a_degree':min((m>>6 for m in p[0]),default=-1)}
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 ib=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63';cb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_CANONICAL_20261104_v63';rb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63';eb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_EXCEPTION_SQUAREFREE_20261104_v63';failed=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63'
 ir=read(ib/'report.json');cr=read(cb/'report.json');assert sha(ib/'report.json')=='837fd4a3800fddd470a4fb0c8e7f82100c19613e3355db168ab9fe4211da8c6f';assert sha(cb/'report.json')=='f0920cee69c47a05f5030e6b879256764a7662e37b9568d78629983ab124c2db'
 # Source and all published evidence are bound before arithmetic.
 for b,r in [(ib,ir),(cb,cr),(rb,read(rb/'report.json')),(eb,read(eb/'report.json'))]:
  for name,v in r['evidence_sha256'].items():assert sha(b/name)==v,(b,name)
 c=read(ib/'certificate.json');n=read(cb/'native-certificate.json.gz');original=read(rb/'certificate.json');exc=read(eb/'certificate.json')
 assert c['h1_e16_monic']==exc['old_residual_LC_preflight'][0]['monic_normalized']and c['k1_e25_monic']==exc['old_residual_LC_preflight'][1]['monic_normalized']
 H,K=poly(c['h1_e16_monic']),poly(c['k1_e25_monic']);assert stats(H)['e_degree']==16 and stats(K)['e_degree']==25
 assert {(m>>6):v for m,v in H[0].items()if(m&63)==16}=={0:H[1]}
 assert {(m>>6):v for m,v in K[0].items()if(m&63)==25}=={0:K[1]}
 qs=[poly(p)for p in c['all16_column_quotients_source_order']];M=[[poly(p)for p in row]for row in c['matrix_rows_source_order']];assert len(M)==16 and all(len(row)==16 for row in M)
 assert all((m&63)==0 for row in M for p in row for m in p[0])
 for j in range(16):check([(H,qs[j])]+[(shift(M[i][j],e=i),ONE)for i in range(16)],shift(K,e=j))
 HL,KL=c['source_original_QQleading_units'];Horig=poly(original['oldh_quotient_h1']);Korig=poly(original['newK_quotient_k1']);equal(Horig,scale(H,HL));equal(Korig,scale(K,KL))
 D=int(c['determinant_cost_bounds_only']['matrix_uniform_QQ_denominator']);assert D==int(n['integer_clearing_D'])and D>0
 N=[]
 for row in M:
  z=[]
  for p in row:
   assert D%p[1]==0;z.append(({m:v*(D//p[1])for m,v in p[0].items()},1))
  N.append(z)
 V=[poly(p)for p in n['only16_solution_cofactors']];delta=poly(n['fflu_solver_delta']);assert len(V)==16 and delta[0]and delta[1]==1 and all(p[1]==1 for p in V)
 for i in range(16):check([(N[i][j],V[j])for j in range(16)],delta if i==0 else({},1))
 Ve=poly(n['V_e_polynomial']);qsum=poly(n['Qsum']);check([(shift(V[j],e=j),ONE)for j in range(16)],Ve);check([(V[j],qs[j])for j in range(16)],qsum)
 hc,kc=poly(n['monic_h1_cofactor']),poly(n['monic_k1_cofactor']);equal(hc,scale(qsum,-D));equal(kc,scale(Ve,D));check([(hc,H),(kc,K)],delta)
 hco,kco=poly(n['original_h1_cofactor']),poly(n['original_k1_cofactor']);equal(scale(hco,HL),hc);equal(scale(kco,KL),kc);check([(hco,Horig),(kco,Korig)],delta)
 content=Fraction(n['delta_QQ_content']);primitive=poly(n['delta_primitive']);res=poly(n['remaining_non_a_annihilator']);assert content>0 and content.denominator==1
 assert math.gcd(*(abs(v)for v in delta[0].values()))==content.numerator
 equal(delta,scale(primitive,content));power=n['removed_known_a_power'];assert power==31;equal(primitive,shift(res,a=power))
 assert stats(delta)==dict(terms=118,a_degree=148,e_degree=0,minimum_a_degree=31)
 assert stats(res)==dict(terms=118,a_degree=117,e_degree=0,minimum_a_degree=0)and res[0][0]!=0
 assert res[1]==1 and math.gcd(*(abs(v)for v in res[0].values()))==1
 # The determinant object is only byte/coefficient binding, NOT re-determined.
 equal(poly(n['actual_det16']),delta)
 storage=read(cb/'storage-receipt.json');assert storage['original_aggregate_guard_FAILED']and storage['original_postmath_guard_bytes']==34721772>storage['limit_bytes']==33554432
 restoration=[]
 for e in storage['entries']:
  h=hashlib.sha256();size=0
  with gzip.open(cb/e['compressed_filename'],'rb')as z:
   for data in iter(lambda:z.read(1024*1024),b''):h.update(data);size+=len(data)
  assert h.hexdigest()==e['raw_sha256']and size==e['raw_bytes']and sha(cb/e['compressed_filename'])==e['compressed_sha256']
  p=ROOT/e['raw_workspace_path'];local_matches=sha(p)==e['raw_sha256']if p.exists()else None
  assert local_matches is not False;restoration.append(dict(filename=e['compressed_filename'],raw_bytes=size,raw_SHA=h.hexdigest(),local_raw_SHA_matches=local_matches))
 checkpoint=read(cb/'native-certificate.json.completed-det-solve.json.gz');assert checkpoint['fflu_solver_delta']==n['fflu_solver_delta']and checkpoint['only16_solution_cofactors']==n['only16_solution_cofactors']and checkpoint['actual_det16']==n['actual_det16']
 failure=read(failed/'failure.json');assert failure['status']=='det16_and_NVsolve_passed_but_Bezout_structural_equality_failed'and failure['no_mathematical_norm_conclusion_promoted']
 cert=dict(source_variables=['a','e'],independent_arithmetic='Pythonintegerconvolution, perpolynomialnumericQQclearing, sixe-bits withallproducte-degree<64; noauthorhelper/nativeimport',source_monic_divisions_checked=16,integer_Nv_entries_checked=16,all256_N_entries_bound=True,original_h1_k1_Bezout_fullQQ=True,Ve_Qsum_and_cofactor_normalization_checked=True,only_known_QQ_content_and_a31_removed=True,nonzero_primitive_a117=encode(res),necessary_cover_only='Anyoldh1=k1=0pointwitha!=0hasres(a)=0. ActualBbadnonempty/empty isNOTclaimed.',determinant_not_recomputed_or_needed_for_annihilator=True,canonical_primary_equals_delta_only_checked_as_saved_coefficients=True,storage_restore=restoration,original_guard_failure_preserved=True,first_canonicalization_failure_preserved=True,all_new_artifacts_small=True)
 (out/'independent-certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',ib/'report.json',ib/'certificate.json',cb/'report.json',cb/'native-certificate.json.gz',cb/'storage-receipt.json',cb/'native-certificate.json.completed-det-solve.json.gz',cb/'postcheck-artifact-guard.log',rb/'report.json',rb/'certificate.json',eb/'report.json',eb/'certificate.json',failed/'failure.json',failed/'failed-source.cpp',failed/'canonicalization-fixture.cpp',ROOT/'code/s61_nonsquare_raw_H_B_old_residual_norm16_canonical_20261104_v63.cpp']
 report=dict(id=ID,all_passed=True,scope='Independent fullQQ source divisions/NV/originalBezout andnecessarya117cover ONLY',source_divisions16=True,all16_Nv_entries=True,full_original_h1_k1_Bezout=True,nonzero_a117_necessary_cover=True,no_determinant_recomputation=True,no_author_worker_point_prime_GB_or_solver=True,canonicalization_and_34MB_guard_failures_preserved=True,all_lossless_storage_SHA_verified=True,arithmetic_metrics=METRICS,max_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
