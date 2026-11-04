#!/usr/bin/env python3
"""Independent fixed actual nine-minor and saved QQ gcd certificates."""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, json, resource, subprocess, sys, time

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63_REVIEW'
BASE='reports/'
RSC='S61_NONSQUARE_P100_K7W_FINITE15_RESCUE_20261104_v63'
GCD='S61_NONSQUARE_P100_K7W_FINITE15_QQ_GCD_20261104_v63'
SRC='S61_NONSQUARE_P100_K7W_FINITE15_SOURCE_20261104_v63'
PRE='S61_NONSQUARE_P100_K7W_FINITE15_PREFLIGHT_20261104_v63'
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',default=BASE+STEM)
 ap.add_argument('--resume-final-modular',action='store_true',help='Resume only the saved final reader frontier; do not repeat completed determinant/QQ checks.')
 a=ap.parse_args()
 out=ROOT/a.output;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 hashes={}
 def load(p):hashes[p]=sha(p);return read(p)
 def bind_report(p):
  d=load(p)
  for q,h in d.get('input_sha256',{}).items():assert sha(q)==h,q;hashes[q]=h
  for name,h in d.get('evidence_sha256',{}).items():
   q=str(Path(p).parent/name);assert sha(q)==h,q;hashes[q]=h
  return d
 rr=bind_report(BASE+RSC+'/report.json');gr=bind_report(BASE+GCD+'/report.json')
 assert rr['complete'] and gr['complete'] and gr['exact_QQ_gcd_degree']==2
 assert not rr['whole_P100_original_quartic_pointwise_lower150'] and not gr['whole_P100_original_quartic_pointwise_lower150']
 nr=load(BASE+RSC+'/native.json');ng=load(BASE+GCD+'/native.json')
 source=load(BASE+SRC+'/source.json');pre=load(BASE+PRE+'/certificate.json')
 sp=load(BASE+SRC+'/report.json')
 assert sha(BASE+SRC+'/source-input.txt')==sp['native_input_sha256']
 assert sha(BASE+SRC+'/source.json')==sp['source_sha256']
 for p,h in sp['input_sha256'].items():assert sha(p)==h,p;hashes[p]=h
 actual_source_peer=load('reports/S61_NONSQUARE_P100_K7W_PURE98_20261104_v63_REVIEW/report.json')
 assert actual_source_peer['all_passed'] and actual_source_peer['actual_C68_times_W_row226_source_and_NF98_accepted']
 operator_path='reports/S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63/operator.json'
 operator=load(operator_path)
 matrix=operator['nodes'][operator['exports']['M15x98']]
 typed_peer=load('reports/S61_NONSQUARE_P100_PURE90_TYPED_SPLIT_20261104_v63_REVIEW/report.json')
 assert typed_peer['all_passed'] and typed_peer['input_sha256'][operator_path]==hashes[operator_path]
 assert matrix['shape']==[15,98] and matrix['op']=='external_row_matrix'
 selected_rows=[9,2,3,5,7,6,1,4];targets=[0,1,2,3,4,5,6,9,28]
 assert source['case9_top8_source_ids']==selected_rows==pre['old_source_rows_in_order8']
 assert source['nine_oldlocal98_targets']==targets==pre['nine_oldlocal98_target_columns']
 assert source['original122_targets']==[0,1,2,4,5,6,7,10,33]
 assert pre['old8cover_basis_index']==9
 rows=[];row_receipts=[]
 for i in selected_rows:
  ref=matrix['row_refs'][i];assert ref['json_pointer']=='/values98' and sha(ref['path'])==ref['sha256']
  rows.append(load(ref['path'])['values98']);row_receipts.append({'M15_row':i,**ref})
 bottom=source['NEW_bottom_source'];assert bottom==pre['new_bottom_source']
 assert bottom['original319_source_index']==226 and bottom['literal_normalizedC68_times_physicalv3']
 assert sha(bottom['file'])==bottom['sha256']
 assert actual_source_peer['input_sha256'][bottom['file']]==bottom['sha256']
 rows.append(load(bottom['file'])['values98'])
 old_gcd_path='reports/S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63/native.json'
 old_gcd=load(old_gcd_path)
 old_cover=load('reports/S61_NONSQUARE_P100_NEXT9_COVER_20261104_v63/native.json')
 old_report=load('reports/S61_NONSQUARE_P100_NEXT9_QQ_GCD_20261104_v63/report.json')
 assert old_report['complete'] and old_report['exact_QQ_gcd_degree']==15
 assert old_report['evidence_sha256']['native.json']==hashes[old_gcd_path]
 H=old_gcd['exact_gcd_monic'];assert H==nr['monic_H15']==ng['original_H15']
 assert nr['finite_remainder_R']==ng['original_saved_R14']
 assert source['monicH15_ref']['sha256']==hashes[old_gcd_path]
 # Independently rebuild the exact original flat RF input, not the author
 # builder. This also binds the selected 81 values to their actual row files.
 text=['51034','32003'];punit_count=0;rf_terms=0
 def flatpoly(p):
  nonlocal punit_count
  text.append(str(len(p)))
  for m,c in p:
   assert len(m)==1 and Q(c).denominator%32003
   punit_count+=1;text.append(str(m[0])+' '+c)
 flatpoly(H);text.append('9');text+=[' '.join(map(str,e))for e in pre['nine_LCD_powers5']]
 for row in rows:
  for j in targets:
   for key in ('numerator','denominator'):
    p=row[j][key];rf_terms+=len(p);flatpoly(p)
 assert rf_terms==1057
 assert '\n'.join(text)+'\n'==(ROOT/BASE/SRC/'source-input.txt').read_text()
 del text
 pool=[];cache={}
 def pol(p):
  key=json.dumps(p,separators=(',',':'))
  if key not in cache:cache[key]=len(pool);pool.append(p)
  return cache[key]
 refs=[[(pol(row[j]['numerator']),pol(row[j]['denominator']))for j in targets]for row in rows]
 fields=[nr['exact_determinant_N'],H,nr['division_quotient'],nr['finite_remainder_R'],ng['primitive_ZZ_H'],ng['primitive_ZZ_R'],
  [[[0],ng['H_equals_content_times_primitive_content']]],[[[0],ng['R_equals_content_times_primitive_content']]],ng['exact_gcd_monic'],
  ng['QQ_Bezout_primitiveH_cofactor'],ng['QQ_Bezout_primitiveR_cofactor'],ng['H_div_gcd'],ng['R_div_gcd']]
 ids=list(map(pol,fields));old_ids=[]
 for case,quot,bz in zip(old_cover['determinants5'],old_gcd['all5_quotients'],old_gcd['cofactor_modular_Bezout5']):
  old_ids.append([pol(case['primitive_ZZ_polynomial']),pol(quot),pol([[[e],str(c)]for e,c in bz])])
 work=ROOT/'work/P100-K7W-finite-independent-review';work.mkdir(parents=True,exist_ok=True)
 ip=work/'fixed-Laplace-and-literal-gcd-input.txt'
 with ip.open('w')as f:
  f.write(str(len(pool))+'\n')
  for p in pool:
   f.write(str(len(p))+'\n')
   for m,c in p:f.write(str(m[0])+' '+c+'\n')
  for row in pre['nine_LCD_powers5']:f.write(' '.join(map(str,row))+'\n')
  for row in refs:
   for x in row:f.write(' '.join(map(str,x))+'\n')
  f.write(' '.join(map(str,ids))+'\n')
  for row in old_ids:f.write(' '.join(map(str,row))+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=work/'verify-fixed-K7W'
 hashes[str(Path(__file__).relative_to(ROOT))]=sha(str(Path(__file__).relative_to(ROOT)))
 hashes[str(cpp.relative_to(ROOT))]=sha(str(cpp.relative_to(ROOT)))
 note='research/notes/'+STEM+'.md';hashes[note]=sha(note)
 preflight={'source_matrix':[9,9],'fixed_old_case':9,'top8_rows':selected_rows,'target98_columns':targets,
  'new_bottom_actual319_row':226,'all_source_RF_terms':rf_terms,'original_source_coefficient_punit_checks':punit_count,
  'source_pool_polynomials':len(pool),'source_pool_sparse_terms':sum(map(len,pool)),
  'input_bytes':ip.stat().st_size,'input_sha256':hashes,'ignored_native_input_sha256':sha(str(ip.relative_to(ROOT))),
  'no_new_source_query_point_prime_rank_or_gcd_search':True}
 assert ip.stat().st_size<32<<20
 (out/'preflight-and-source-binding.json').write_text(json.dumps(preflight,indent=2)+'\n')
 if a.resume_final_modular:
  frontier=work/'final-modular-reader-frontier';failure=json.loads((frontier/'failure-and-prefix-proof.json').read_text())
  assert failure['all_fixed_N202_and_QQ_G2_assertions_completed_before_last_modular_assertion']
  assert failure['native_input_sha256']==sha(str(ip.relative_to(ROOT)))
  assert failure['executable_sha256']==sha(str(exe.relative_to(ROOT)))
  failed_source=(frontier/cpp.name).read_text()
  assert 'assert(eq(dp[511],N))' in failed_source
  assert failed_source.index('assert(eq(add(mul(S,PH),mul(T,PR)),G))')<failed_source.index('assert(nmod_poly_is_one(sum.v))')
  assert failed_source.index('assert(eq(mul(G,QH),H)&&eq(mul(G,QR),R))')<failed_source.index('assert(nmod_poly_is_one(sum.v))')
  assert 'nmod_poly_is_one(sum.v)' in (frontier/'native.stderr').read_text()
  for p,h in failure['sha256'].items():
   local=frontier/Path(p).name;assert hashlib.sha256(local.read_bytes()).hexdigest()==h
  (out/'failed-final-modular-reader.cpp').write_text(failed_source)
  (out/'failure-and-prefix-proof.json').write_text(json.dumps(failure,indent=2)+'\n')
  result={'all_passed':True,'all81_actual_source_LCD_checks':81,
   'fixed_subset_Laplace_compared_N202':True,'exact_N202_H15_Q187_R14':True,
   'QQ_gcd_from_literal_Bezout_and_both_divisions':True,
   'prefix_execution_gate':'The saved run reached only its final old-five modular assertion; all preceding unconditional fixed determinant, source and QQ assertions completed.',
   'old_reader_callback_failure_repaired_without_repeating_completed_polynomial_arithmetic':True,
   'peak_live_slots_including_reserve_bound':150000}
 else:
  flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
  with(out/'compile.log').open('w')as f:
   subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=15)
  with(out/'native.log').open('w')as f,(out/'native.stderr').open('w')as e:
   subprocess.run([str(exe),str(ip)],stdout=f,stderr=e,check=True,timeout=max(1,58-(time.monotonic()-started)))
  result=json.loads((out/'native.log').read_text());assert result['all_passed']
  assert result['peak_live_slots_including_reserve']<=150000
 # Independently reduce the original QQ cofactor coefficients using Python
 # integers, compare every saved reduction, and multiply the given Bezout.
 # This avoids the native low-level preinverse callback entirely.
 prime=32003;modsum={};modular_products=0;reduced_coefficients=0
 for quotient,bezout,saved in zip(old_gcd['all5_quotients'],old_gcd['cofactor_modular_Bezout5'],old_gcd['cofactor_reductions5']):
  reduced={}
  for m,c in quotient:
   q=Q(c);assert q.denominator%prime
   v=(q.numerator%prime)*pow(q.denominator%prime,-1,prime)%prime
   if v:reduced[m[0]]=v
   reduced_coefficients+=1
  assert reduced==dict(saved)
  assert max(reduced)==quotient[-1][0][0] # leading degree preserved
  for e,c in reduced.items():
   for f,b in bezout:
    degree=e+f;v=(modsum.get(degree,0)+c*b)%prime
    if v:modsum[degree]=v
    else:modsum.pop(degree,None)
    modular_products+=1
 assert modsum=={0:1}
 result['old_five_QQ_to_modular_coefficients_checked']=reduced_coefficients
 result['old_five_fixed_modular_Bezout_products']=modular_products
 result['old_five_saved_cofactor_Bezout_equals_one']=True
 proof={'all_passed':True,'actual9minor_source_binding':row_receipts+[bottom],
  'independent_determinant_algorithm':'column-subset Laplace expansion; no pivoting or division',
  'all_source81_LCDs_known_five_units':True,'exact_N202_equals_H15_Q187_plus_R14':True,
  'exact_QQ_gcd':'z^2+(200/13)z+22800/13',
  'gcd_proof':'g=S*primitiveH+T*primitiveR, H=g*QH, R=g*QR, nonzero scalar content identities',
  'old_five_H_gate':'saved five source minors accepted by parent; all five primitive/H divisions and old-prime cofactor Bezout independently rechecked',
  'six_fixed_minors_common_frontier_only':True,'whole_P100_pointwise149_preserved':True,
  'pointwise150_on_D_g_only':True,'actual_rank_failure_or_whole150_or_full239_claim':False,
  'failed_author_compile_header_receipt_preserved':True,'no_factorization_or_new_minor':True,
  'old_five_minor_gcd_H15_source_and_literal_certificate_checked':True,
  'native_result':result}
 (out/'independent-fixed-minor-and-QQ-gcd-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 evidence={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()}
 report={'id':STEM,'all_passed':True,'frozen':True,'actual_K7W_fixed9minor_and_exact_QQ_gcd_verified':True,
  'old_five_minors_exact_common_gcd_H15_verified':True,
  'exact_six_chosen_minor_common_frontier_degree':2,'whole150_or_actual_rank_failure_claim':False,
  'input_sha256':hashes,'evidence_sha256':evidence,
  'ignored_native_input_sha256':sha(str(ip.relative_to(ROOT))),'ignored_executable_sha256':sha(str(exe.relative_to(ROOT))),
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32<<20
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
