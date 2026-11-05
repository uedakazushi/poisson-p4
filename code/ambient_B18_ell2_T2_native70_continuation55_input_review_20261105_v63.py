#!/usr/bin/env python3
"""Independent continuation55 code/custody/positive-integer source gate."""
import argparse,hashlib,itertools,json,math,os,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_INPUT_GATE_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={}
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 base='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_INPUT_GATE_20261105_v63/';r=read(base+'report.json');assert pins[base+'report.json']=='4bcf154e40c19e07b63bb12ee6e9db8694c7a79fd9969dd5fb3f1692b6c81024' and r['all_passed']
 d=read(base+'cached15-plus-new55-and-known-witness-bordered-rank-source-input.json');assert pins[base+'cached15-plus-new55-and-known-witness-bordered-rank-source-input.json']==r['evidence_sha256']['cached15-plus-new55-and-known-witness-bordered-rank-source-input.json']
 for p,z in r['input_sha256'].items():assert sha(ROOT/p)==z['sha256']
 prior='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63_REVIEW/report.json';assert read(prior)['all_passed'] and pins[prior]=='606f7b6e7b30465cae578ed99d5c25ceef74951ce2366823c34cf2f6fcee0671'
 primary='reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63/';prefix=read(primary+'actual-native70-lower-prefix-and-witness-if-any.json');assert d['cached_first15_prefix']==prefix['computed_prefix'] and [v['case619'] for v in d['cached_first15_prefix']]==list(range(549,564))
 assert d['new_case619_indices55']==list(range(564,619)) and d['rank_start_witness_column_index70']==14 and d['rank_start_witness_case619']==563
 whole=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63/complete70-cached-input-source-and-positive-whole-image-bounds.json');DU=[tuple(map(int,v)) for v in whole['all70_four_images_global_D_U']];assert [list(map(str,v)) for v in DU]==d['row_integer_clearing_D_U']
 W=[]
 for (D,U),v in zip(DU,prefix['nonzero_witness']['actual_lower4_values'],strict=True):
  x=Q(v);assert D%x.denominator==0;w=x.numerator*(D//x.denominator);assert 1<=abs(w)<=U and abs(w).bit_length()<=128 and U.bit_length()<=128;W.append(w)
 assert list(map(str,W))==d['known_integer_scaled_witness4']
 bounds=[];maximum=0
 for size in range(1,5):
  for rows in itertools.combinations(range(4),size):
   b=0
   for i in rows:
    p=math.factorial(size-1)*abs(W[i])
    for j in rows:
     if j!=i:p*=DU[j][1]
    b+=p
   maximum=max(maximum,b.bit_length());bounds.append({'rows':list(rows),'permanent_bound_bits':b.bit_length()})
 assert bounds==d['all15_row_subset_permanent_bounds'] and maximum==374<d['rank_temporary_integer_bits']==385
 assert 128+max(D.bit_length() for D,U in DU)<385
 a=(ROOT/'code/ambient_B18_ell2_T2_native70_values_20261105_v63.py').read_text();b=(ROOT/'code/ambient_B18_ell2_T2_native70_continuation55_20261105_v63.py').read_text()
 assert a[a.index(' table={(i,j)'):a.index(' witness=None')]==b[b.index(' table={(i,j)'):b.index(' witness=None')]
 assert a[a.index('  case=rec['):a.index('  computed.append')]==b[b.index('  case=rec['):b.index('  computed.append')]
 assert "zip(d['source70_records'][15:],d['all70_four_component_positive_envelopes'][15:],strict=True)" in b and 'computed=list(cd[\'cached_first15_prefix\'])' in b
 assert "assert cols[0]==known and len(cols)==len(rows)<=4" in b and "selected=[known]" in b and "if rank_temp_bits>385:" in b and 'no_large_minor_or_QQ_factor_committed' in b
 cp='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_PREPARED_CUSTODY_20261105_v63/report.json';c=read(cp);assert c['all_passed'] and c['safe_help_literal_PASS'] and c['worker_sha256']==sha(ROOT/c['worker'])=='501e0a7df82e37ab4d37a6aea1ec65011182c31b423caeecb358f7c8db896b64'
 for p,h in c['worker_and_helper_sha256'].items():assert sha(ROOT/p)==h
 for p,h in c['literal_external_source_sha256'].items():assert sha(ROOT/p)==h
 assert c['continuation_gate_sha256']==pins[base+'report.json']
 h=subprocess.run([sys.executable,str(ROOT/c['worker']),'--help'],cwd=ROOT,capture_output=True,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'},timeout=10);assert h.returncode==0 and not h.stderr and hashlib.sha256(h.stdout).hexdigest()==c['safe_help_stdout_sha256']
 assert 49734+d['cached_prefix_two_owned_copies']+d['extra_rank_working_reserve']+d['extra_profile_and_control_reserve']==d['whole_conservative_owned']==50884<d['unchanged_whole_numeric_cap']==51263
 assert d['cached_prefix_two_owned_copies']==2*count(d['cached_first15_prefix'])==510
 witnessmax=max(17+count(z['actual_a1_27_times9'])+count(z['original768_constant_combination']) for z in whole['source70_records'][15:]);rankoutmax=max(70*(5-k)+7*k+(5+2*k)+8 for k in range(1,5));assert witnessmax==265 and rankoutmax==302 and 2*70*17+witnessmax+rankoutmax==2947<3000
 scratch=4*70+4*4+4*4+4+24*4+4*4+4+4+4+8+4+4;assert scratch<540
 ev={'all15_cached_prefix_read_exact_not_recomputed':True,'new_cases_only_564_through618':True,'finite_cache_and_individual_case_core_byte_identical_to_reviewed_worker':True,'all15_source_permanent_bounds':bounds,'integer_rank_permanent_max374':True,'all_W_and_U_integer_entries_nonzero':True,'bordered_minor_rank_and_final_span_proof':True,'no_evaluated55_images_or_rank':True,'explicit_output_overlap':{'two70_records':2380,'largest_retained_future_witness':265,'rank_certificate_max':302,'total':2947,'reserve':3000},'rank_scratch_bound':scratch,'whole_bound':50884,'safe_help_literal_PASS':True}
 ep=out/'independent55-source-prefix-bordered-rank-height-and-custody.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n');owned=count([d,whole,prefix,c,r,bounds,ev])+10000;assert owned<60000
 report={'id':STEM,'all_passed':True,'source_code_custody_only':True,'no_new55_images_or_rank_queried':True,'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'owned_numeric_bound':owned,'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);x=p.parse_args();main(x.output)
