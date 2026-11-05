#!/usr/bin/env python3
"""Seven fixed old socle contractions of already frozen values only."""
import hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOLD_SOCLE_SEVEN_VALUES_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=True);assert not(out/'report.json').exists();pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(out/n).write_text(json.dumps(v,indent=2)+'\n')
 gate=load('AMBIENT_B18_TOLD_SOCLE_SEVEN_SOURCE_PREFLIGHT_20261105_v63');assert gate['all_passed'] and pins['reports/AMBIENT_B18_TOLD_SOCLE_SEVEN_SOURCE_PREFLIGHT_20261105_v63/report.json']=='9267f4a72212a8d7dd9ccf05061e1db873d69dd525867f3c05476542af9202a3'
 P=load('AMBIENT_B18_TOLD_SOCLE_SEVEN_SOURCE_PREFLIGHT_20261105_v63','seven-old-socle-source-functionals-and-existing-cache-gate.json')
 K=load('AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63','actual-linear-syzygies-and-pure-r3-source-values.json')['actual37_e5_r2_source_coefficients']
 C=load('AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63','phase3-composite-second-source-functional.json');cache={i:({int(j):F(c)for j,c in v},g)for i,v,g in C['all108_composite_values_and_opaque_second_gradients']}
 assert C['B11_monomial_order'][2]==[0,0,0,0,0,0,0,2]
 save('phase.json',{'phase':'only seven fixed contractions; no new source coefficient'})
 roots=P['needed_cached_C0root_ids'];assert all(not cache[i][1]for i in roots);rootvals={i:cache[i][0].get(2,F(0))for i in roots};kv=[F(v['actual_e5_r2_coefficient'])for v in K]
 direct=[];folded=[];certs=[]
 for f,g in zip(P['seven_socle_obstruction_functionals_on_existing37'],P['seven_equivalent_C0root_r_squared_functionals']):
  v=sum((F(c)*kv[int(j)]for j,c in f.items()),F(0));w=sum((F(c)*rootvals[int(j)]for j,c in g.items()),F(0));assert v==w;direct.append(str(v));folded.append(str(w));certs.append({'source37_terms':[[int(j),c,str(kv[int(j)]),str(F(c)*kv[int(j)])]for j,c in f.items()],'folded_root_terms':[[int(j),c,str(rootvals[int(j)]),str(F(c)*rootvals[int(j)])]for j,c in g.items()]})
 positive=any(F(v)for v in direct)
 proof={'seven_actual_QQ_values':direct,'seven_independent_folded_root_values':folded,'literal_product_certificates':certs,'actual_existing_root_values_r2':{str(i):str(c)for i,c in rootvals.items()},'criterion':P['sufficiency'],'any_nonzero':positive,'all_seven_zero':not positive,'no_actual_tail_chosen_if_all_zero':not positive,'new_physical_source_or_jet_values':0,'complete917_rank_or_full_tail_not_inferred_from_zero':not positive}
 save('seven-fixed-old-socle-values-and-source-folds.json',proof)
 for p in[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep={'id':STEM,'all_passed':True,'frozen':True,'all_seven_actual_obstructions_zero':not positive,'FULL_HF55321_forced':positive,'failure_of_sufficient_test_only_if_all_zero':not positive,'no_new_source_jet':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20;save('report.json',rep);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
