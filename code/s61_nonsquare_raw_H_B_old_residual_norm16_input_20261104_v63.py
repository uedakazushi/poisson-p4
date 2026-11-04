#!/usr/bin/env python3
"""ONE oldresidual16 normmatrix INPUT only, assignmentdegree/height bounds."""
import hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar);base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_EXCEPTION_SQUAREFREE_20261104_v63';ip=base/'certificate.json';assert sha(ip)==load(base/'report.json')['evidence_sha256'][ip.name];c=load(ip);old=c['old_residual_LC_preflight'];H=ar.poly(old[0]['monic_normalized']);K=ar.poly(old[1]['monic_normalized']);assert {(m[0],m[2]):v for m,v in H.items()if m[1]==16}=={(0,0):Q(1)}
 # Actualmonic h1/k1 domainQQ[a,e]; native reorder makes e leading variable.
 native=lambda p:{(m[1],m[0],m[2]):v for m,v in p.items()};inp=out/'native-input.txt';inp.write_text(ar.line(native(H))+ar.line(native(K)));cpp=ROOT/'code/s61_nonsquare_raw_H_B_old_residual_norm16_input_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-old-residual-norm16-input'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','norm16input-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');un=lambda p:{(m[1],m[0],m[2]):Q(v)for m,v in p if Q(v)};qs=list(map(un,n['column_quotients']));rs=list(map(un,n['column_remainders']));M=[[un(p)for p in row]for row in n['matrix_rows']];one={(0,0,0):Q(1)}
 for j in range(16):
  rebuild={}
  for i in range(16):rebuild=ar.add(rebuild,ar.shift(M[i][j],(0,i,0)))
  assert rebuild==rs[j];ar.check([(qs[j],H),(rs[j],one)],ar.shift(K,(0,j,0)))
 allp=[p for row in M for p in row];terms=sum(map(len,allp));assert terms<150000 and all(m[1]==0 and m[2]==0 for p in allp for m in p)
 den=math.lcm(*(v.denominator for p in allp for v in p.values()));mx=max(abs(int(v*den))for p in allp for v in p.values());support=max(map(len,allp));height=math.factorial(16)*(support*mx)**16;degrees=[[max([m[0]for m in p],default=-100000)for p in row]for row in M]
 # Fixedinteger combinatorialassignment upper, no determinantarithmetic.
 dp=[-10000000]*(1<<16);dp[0]=0
 for mask in range(1<<16):
  i=mask.bit_count()
  if i==16 or dp[mask]<-100000:continue
  for j in range(16):
   if not(mask>>j&1)and degrees[i][j]>=0:
    v=dp[mask]+degrees[i][j]
    if v>dp[mask|1<<j]:dp[mask|1<<j]=v
 dega=dp[-1];assert dega>=0;bounds={'matrix_uniform_QQ_denominator':str(den),'matrix_clear_integer_max_coeff_bits':mx.bit_length(),'matrix_max_entry_a_degree':max(max(r)for r in degrees),'norm16_a_degree_upper_assignment':dega,'norm16_support_upper':dega+1,'cleared_integer_norm16_coefficient_bits_upper':height.bit_length(),'conservative_norm16_serialized_bytes_upper':(dega+1)*(100+(height.bit_length()*30103+99999)//100000),'norm16_not_computed':True};assert bounds['conservative_norm16_serialized_bytes_upper']<32*1024*1024
 cert={'source_variable_order':['a','e','d_unused'],'h1_e16_monic':old[0]['monic_normalized'],'k1_e25_monic':old[1]['monic_normalized'],'source_original_QQleading_units':[old[0]['LC_is_QQconstant_global_unit'],old[1]['LC_is_QQconstant_global_unit']],'all16_column_quotients_source_order':[[[list(m),str(v)]for m,v in sorted(p.items(),reverse=True)]for p in qs],'matrix_rows_source_order':[[[[list(m),str(v)]for m,v in sorted(p.items(),reverse=True)]for p in row]for row in M],'all16_source_identities_independent_integer_verified':True,'matrix_a_degrees':degrees,'determinant_cost_bounds_only':bounds,'source_algebra':'QQ[a][e]/(h1monic) finitefree16; multiplyingk1. H4 wouldgive64base directions beforethek1constraint.','no_unknown_LC_inverted':True,'norm_nonzero_by_saved_h1_k1_coprimality_and_Gauss':True,'no_determinant_resultant_GB_point_prime_job':True};(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip];report={'id':ID,'matrix_shape':[16,16],'matrix_total_terms':terms,'max_entry_terms':support,'bounds':bounds,'native_seconds':n['native_seconds'],'peak_native_support':n['peak_single_polynomial_support'],'all16_independent_integer_identities':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
