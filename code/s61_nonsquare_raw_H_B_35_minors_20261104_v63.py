#!/usr/bin/env python3
"""ONE fixed35determinants/content job, no subsequent GB authorization implied."""
import argparse,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 cb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63';cc=read(cb/'certificate.json');cr=read(cb/'report.json');assert sha(cb/'certificate.json')==cr['evidence_sha256']['certificate.json']
 gb=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';g=read(gb/'certificate.json');assert sha(gb/'certificate.json')==read(gb/'report.json')['certificate_sha256']
 # Exact input knownfactors already source-bound; no new factor or sample.
 unitmodel=read(ROOT/'reports/S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63/certificate.json');known=unitmodel['known_factor_polynomials']
 native=out/'native-input.txt';cpp=ROOT/'code/s61_nonsquare_raw_H_B_35_minors_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-35-minors'
 def line(p):return str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in p)+'\n'
 with native.open('w')as f:
  f.write('7 3 6\n')
  for p in known:f.write(line(p))
  for row in g['actual_block7x3']:
   for p in row:f.write(line(p))
  for row in cc['prepared_integer7x3']:
   for p in row:f.write(line(p))
  for u in cc['known_row_unit_divisors']+cc['known_column_unit_divisors']:f.write(line(u['polynomial']))
  for d in cc['integer_column_multipliers_QQ_units']:f.write(line([[[0,0,0],str(d)]]))
  f.write(' '.join(map(str,g['selected_rows7']+g['selected_cols3']))+'\n');f.write(line(g['nonzero_selected_minor_Delta']))
 command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as f:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','B35-native',str(exe),str(native),str(out/'native-certificate.json')],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,56-(time.monotonic()-start)))
 nr=read(out/'native-certificate.json');assert nr['all_passed']and nr['minor_count']==35 and nr['all35_source_scalings_verified']
 bounds={tuple(t['rows7']):t['support_upper']for t in cc['all35_minor_SUPPORT_bounds_only']};primitive=[];seen={};recordmap=[]
 for r in nr['records']:
  assert len(r['prepared_integer_minor'])<=bounds[tuple(r['rows7'])]
  if r['zero']:recordmap.append(None);continue
  p={tuple(m):Q(c)for m,c in r['primitive_integer_polynomial']};lc=p[max(p)];key=tuple((m,str(c/lc))for m,c in sorted(p.items(),reverse=True))
  if key not in seen:seen[key]=len(primitive);primitive.append([[list(m),c]for m,c in key])
  recordmap.append(seen[key])
 assert nr['nonzero_prepared_minor_totalterms']<=7172
 degree=lambda p:max([sum(m)for m,c in p]+[0]);bits=lambda p:max([max(abs(Q(c).numerator).bit_length(),Q(c).denominator.bit_length())for m,c in p]+[0])
 cost={'primitive_distinct_generators':len(primitive),'primitive_distinct_totalterms':sum(map(len,primitive)),'max_primitive_degree':max(map(degree,primitive),default=0),'max_QQ_coefficient_bits':max(map(bits,primitive),default=0),
  'FULLU_product_terms':nr['FULLU_product_terms'],'FULLU_product_degree':nr['FULLU_product_degree'],'one_z_localizer_relation_terms':1+nr['FULLU_product_terms'],'one_z_localizer_relation_degree':1+nr['FULLU_product_degree']}
 cert={'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','all35_minors_in_lexicographic_rowtriple_order':True,'native_certificate_path':'native-certificate.json','native_certificate_sha256':sha(out/'native-certificate.json'),
  'localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling':primitive,'minor35_to_distinct_generator_index':recordmap,'GB_cost':cost,'FULLU_product':nr['FULLU_product'],
  'rankbadlocus_exact_defining_ideal_on_FULLU':'The35 source maximum minors have the same localized ideal as these monic primitive polynomials. No properness/unit assertion until an exact ideal certificate is supplied.',
  'global_B3_proved_by_single_knownunit_minor':nr['single_namedunit_minors']>0,'GB_NOT_executed':True,'no_newpoint_prime_fullH7_or_full239_claim':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');exe.unlink(missing_ok=True)
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),cb/'report.json',cb/'certificate.json',gb/'report.json',gb/'certificate.json',ROOT/'reports/S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63/certificate.json']
 report={'id':ID,'all35_exact_determinants_and_source_unit_scalings_verified':True,'zero_minor_count':nr['zero_minors'],'single_namedunit_minor_count':nr['single_namedunit_minors'],'global_B3_from_single_minor':nr['single_namedunit_minors']>0,
  'native_seconds':nr['native_seconds'],'GB_cost':cost,'GB_jobs':0,'native_executable_sha256_not_committed':exesha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
