#!/usr/bin/env python3
"""ONE norm16 det+oneRHSfflu+originalh1/k1 Bezout; onlyQQ/a removal."""
import hashlib,importlib.util,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_CANONICAL_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(54);out.mkdir(exist_ok=False,parents=True)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 def pin(stem,file='certificate.json'):
  b=ROOT/'reports'/stem;p=b/file;assert sha(p)==load(b/'report.json')['evidence_sha256'][file];return b,p,load(p)
 ib,ip,c=pin('S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63');M=[[ar.poly(p)for p in row]for row in c['matrix_rows_source_order']];H=ar.poly(c['h1_e16_monic']);K=ar.poly(c['k1_e25_monic']);qs=list(map(ar.poly,c['all16_column_quotients_source_order']));D=c['determinant_cost_bounds_only']['matrix_uniform_QQ_denominator'];HL,KL=c['source_original_QQleading_units'];rb,rp,rc=pin('S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63');Horig=ar.poly(rc['oldh_quotient_h1']);Korig=ar.poly(rc['newK_quotient_k1']);assert Horig==ar.scalar(H,Q(HL))and Korig==ar.scalar(K,Q(KL))
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for row in M:
   for p in row:f.write(ar.line(p))
  for p in[H,K]+qs:f.write(ar.line(p))
  f.write(D+' '+HL+' '+KL+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_old_residual_norm16_canonical_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-old-residual-norm16'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','norm16-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');det=ar.poly(n['actual_det16']);ann=ar.poly(n['fflu_solver_delta']);assert det and ann;v=list(map(ar.poly,n['only16_solution_cofactors']));assert len(v)==16
 for i in range(16):ar.check([(ar.scalar(M[i][j],Q(D)),v[j])for j in range(16)],ann if i==0 else{})
 hc=ar.poly(n['monic_h1_cofactor']);kc=ar.poly(n['monic_k1_cofactor']);hco=ar.poly(n['original_h1_cofactor']);kco=ar.poly(n['original_k1_cofactor']);assert ar.scalar(hco,Q(HL))==hc and ar.scalar(kco,Q(KL))==kc;ar.check([(hco,Horig),(kco,Korig)],ann)
 prim=ar.poly(n['delta_primitive']);res=ar.poly(n['remaining_non_a_annihilator']);ar.check([(prim,{(0,0,0):Q(n['delta_QQ_content'])})],ann);ar.check([(res,{(n['removed_known_a_power'],0,0):Q(1)})],prim);assert all(m[1]==m[2]==0 for m in res)and any(m[0]==0 for m in res)
 classification='old_finite_projection_residual_EMPTY_on_aunit'if len(res)==1 and(0,0,0)in res else'nonzero_univariate_a_NECESSARY_cover_for_old_finite_residual'
 proof={'classification':classification,'actual_det16_and_solver_delta_distinguished':True,'solver_delta_equals_det':ann==det,'solver_delta_equals_negative_det':ann==ar.scalar(det,-1),'all16_exact_integer_Nv_delta_e0_checks':True,'full_original_h1_k1_Bezout_independent_integer_checked':True,'delta_in_original_residual_ideal':'delta=original_h1_cofactor*h1+original_k1_cofactor*k1','only16_sourcecofs_not256adjugate':True,'known_unit_removal_only_QQcontent_and_a_power':True,'non_a_annihilator_untouched_factorization_NOT_run':True,'no_unknown_a_factor_or_sourceLC_inverted':True,'no_point_prime_GB_41Sylvester_or_largefinitefiber_job':True,'meaning':'Anyoldfinitebadprojectionwitha!=0mustannihilateres(a). Nonconstantresdoesnotassertactualbadexists. IfresisQQunit, oldprojectionresidualisemptyonFULLU.'}
 (out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63/failure.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63/canonicalization-fixture.cpp',ROOT/'research/notes/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63.md',ib/'report.json',ip,rb/'report.json',rp]
 report={'id':ID,'classification':classification,'actual_det_stats':ar.stats(det),'solver_delta_stats':ar.stats(ann),'solver_delta_equals_actual_det':ann==det,'removed_known_a_power':n['removed_known_a_power'],'remaining_non_a_annihilator_stats':ar.stats(res),'source_Bezout_cofactor_terms':[len(hco),len(kco)],'only16_cofactor_terms':list(map(len,v)),'native_seconds':n['native_seconds'],'peak_native_support':n['peak_single_polynomial_support'],'total_solution_terms':n['total_solution_polynomial_terms'],'all_independent_integer_Nv_and_original_Bezout_checks':True,'GB_point_prime_41Sylvester_factorization_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
