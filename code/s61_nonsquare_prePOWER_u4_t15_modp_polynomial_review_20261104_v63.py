#!/usr/bin/env python3
"""Distinct symbolic-modular audit; the earlier QQ timeout stays immutable."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,subprocess,time
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_PREPOWER_U4_T15_MODP_POLYNOMIAL_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 work=ROOT/'work/prePOWER-u4-t15-symbolic-modp-independent-review';work.mkdir(parents=True,exist_ok=True)
 hashes={}
 def load(p):hashes[p]=sha(p);return read(p)
 def get(st,n='certificate.json'):return load('reports/'+st+'/'+n)
 def save(n,d):(out/n).write_text(json.dumps(d,indent=2)+'\n')
 chart=get('S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63_REVIEW','report.json')
 stop=get('S61_NONSQUARE_PREPOWER_DIRECT_SCHUR_ENTRY_20261104_v63_REVIEW','report.json')
 assert chart['all_passed'] and not stop['all_passed'] and stop['status']=='STOPPED_independent_fixed_arithmetic_at48s_guard_no_retry'
 assert hashes['reports/S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63_REVIEW/report.json']=='3a140e24f3cc2d2c97a9a3632ed15e92f87912f9ad5bbc19e7ac32e7f5c0e7dc'
 assert hashes['reports/S61_NONSQUARE_PREPOWER_DIRECT_SCHUR_ENTRY_20261104_v63_REVIEW/report.json']=='3177220753f445d9ac17c31362dcdbb2bc19561892d92e9e77916470b9cacc11'
 for d in (chart,stop):
  for p,h in d['input_sha256'].items():assert sha(p)==h;hashes[p]=h
  for n,h in d['evidence_sha256'].items():p='reports/'+d['id']+'/'+n;assert sha(p)==h;hashes[p]=h
 original=stop['ignored_input']['path'];assert sha(original)==stop['ignored_input']['sha256']=='b8b29b2e9239548811f64e45d828b49d0e3a6e15969d369422e56e1ecd5a3e74'
 hashes[original]=sha(original);body=(ROOT/original).read_text()
 # This is the independently source-bound old input, not a producer stream.
 # Every rational token is checked before the new arithmetic is launched.
 p=32003;fraction_tokens=0
 for word in body.split():
  if '/' in word:assert Q(word).denominator%p;fraction_tokens+=1
 rank=get('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63')
 strip=get('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json')
 entry=get('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63','native.json')
 op=get('S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63','operator-and-cost.json')
 assert len(op['nodes'])==768 and all(n['bound'].get('known_B_denominator_power',0)==0 for n in op['nodes'])
 assert entry['numerator_terms']==7022 and entry['known_a_clearing_power']==0
 appended=[strip['exact_unscaled_pure4_minor'],strip['remaining_polynomial'],rank['p'],rank['nonzero_selected_minor_Delta'],rank['nonzero_minor_of14x5']]
 def reduction(pol):
  a={}
  for m,c in pol:
   q=Q(c);assert q.denominator%p
   k=tuple(m);a[k]=(a.get(k,0)+q.numerator*pow(q.denominator,-1,p))%p
  return {m:c for m,c in a.items() if c}
 assert all(reduction(pol) for pol in appended)
 assert reduction(entry['exact_Saux_Laurent_polynomial'])
 extra=[]
 for pol in appended:
  extra.append(str(len(pol)))
  for m,c in pol:extra.append(' '.join(map(str,m))+' '+str(c))
 ip=work/'input.txt';ip.write_text(body+'\n'.join(extra)+'\n');assert ip.stat().st_size<32<<20
 pre={'method_change':'QQ exact polynomial replay -> symbolic F32003[a,e,d,a^-1] arithmetic; not an unchanged QQ retry.',
  'no_parameter_specialization_or_new_prime':True,'source768_from_independently_bound_frozen_input':True,
  'all768_atom_B_denominator_powers_zero':True,'all_QQ_fraction_tokens_32003_unit_denominators':fraction_tokens,
  'all5_appended_denominator_chart_polynomials_nonzero_mod32003':True,'source_output7022_reduction_nonzero':True,
  'original_QQ_timeout_and_all_frozen_files_unchanged':True,'cap_seconds':60,'cap_memory_bytes':2<<30,'cap_live_coefficients':150000,'input_bytes':ip.stat().st_size,'input_sha256':sha(str(ip.relative_to(ROOT)))}
 save('preflight.json',pre)
 cpp=Path(__file__).with_suffix('.cpp');exe=work/'checker';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-started)))
 with (work/'native.log').open('w') as log:subprocess.run([str(exe),str(ip),str(out/'native.json'),str(out/'stage.log')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,57-(time.monotonic()-started)))
 result=read('reports/'+STEM+'/native.json');assert result['all_passed'] and result['actual_fixed_bordered5_determinant_verified_nonzero_mod_p']
 assert result['peak_live_coefficients_with_reserves']<150000 and not result['QQ7022_exact_independent_verification_claimed']
 for p in [str(Path(__file__).relative_to(ROOT)),str(cpp.relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'actual_one_entry_char0_NONZERO_independently_proved':True,
  'all768_symbolic_mod32003_source_operations_and7022_public_coefficient_reductions_verified':True,
  'independent_fixed_bordered5_determinant_verified_mod32003':True,
  'all_QQ_denominators_punit_and_p_D4_named_chart_denominators_nonzero_polynomials':True,
  'independent_full_QQ7022_exact_identity_claimed':False,
  'scope':'One actual generic residual(0,9) entry only; no residual rank6, full57, other entry, numeric parameter point, or chart-complement removal.',
  'old_QQ_timeout_frontier_unchanged':True,'input_sha256':hashes,
  'ignored_reproducible_input':{'path':str(ip.relative_to(ROOT)),'sha256':sha(str(ip.relative_to(ROOT)))},
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 save('report.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}),flush=True)
if __name__=='__main__':main()
