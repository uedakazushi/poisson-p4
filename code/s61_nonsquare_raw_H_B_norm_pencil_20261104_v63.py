#!/usr/bin/env python3
"""ONE normpencil mixed3 determinant+exact source GCD; no GB or parameterpoints."""
import hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def pin(base,filename):
 p=base/filename;r=load(base/'report.json');assert sha(p)==r['evidence_sha256'][filename];return p,load(p)
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 gb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_INPUT_20261104_v63';gp,gc=pin(gb,'certificate.json');kb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_INPUT_20261104_v63';kp,kc=pin(kb,'certificate.json');M0=[[ar.poly(p)for p in row]for row in gc['matrix_rows']];M1=[[ar.poly(p)for p in row]for row in kc['matrix_rows']]
 D0=Q(gc['determinant_cost_bounds_NOT_value']['matrix_uniform_QQ_denominator']);D1=Q(kc['cost_bounds_only']['matrix_uniform_QQ_denominator']);a0=gc['determinant_cost_bounds_NOT_value']['matrix_uniform_a_clear_power'];a1=kc['cost_bounds_only']['matrix_uniform_a_clear_power'];assert a0==10 and a1==11 and D1/D0==5000
 A=[[ar.scalar(ar.shift(p,(a1,0,0)),D1)for p in row]for row in M0];B=[[ar.scalar(ar.shift(p,(a1,0,0)),D1)for p in row]for row in M1];allp=[p for row in A+B for p in row];assert all(v.denominator==1 and m[0]>=0 for p in allp for m,v in p.items());maxa=max(m[0]for p in allp for m in p);maxe=max(m[1]for p in allp for m in p);mx=max(abs(int(c))for p in allp for c in p.values());sp=max(map(len,allp));height=24*6*(sp*mx)**4;box=(4*maxa+1)*(4*maxe+1);bytesbound=5*box*(100+(height.bit_length()*30103+99999)//100000);assert 5*box<150000 and bytesbound<32*1024*1024
 nb0=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_20261104_v63';p0,n0=pin(nb0,'native-certificate.json');nb1=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_20261104_v63';p1,n1=pin(nb1,'native-certificate.json');sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';sup,sc=pin(sb,'certificate.json');U=ar.poly(sc['FULLU_product']);G=ar.poly(n1['gcd_newK_oldh']);R0=ar.poly(n0['det_cleared_matrix_R']);R4=ar.poly(n1['det_cleared_matrix_R'])
 preflight={'source_matrix_terms':[sum(len(p)for row in M0 for p in row),sum(len(p)for row in M1 for p in row)],'common_scalar_s':{'QQ':str(D1),'a_power':a1},'oldNormg_to_common_ratio':'5000*a','lambda0_det_expected':'(5000*a)^4*oldR_g','lambda4_det_expected':'oldR_k','max_matrix_degree_a_e':[maxa,maxe],'all5_support_bound':5*box,'each_coefficient_support_bound':box,'integer_coefficient_bits_upper':height.bit_length(),'all5_serialization_bytes_upper':bytesbound,'variable_order':['a','e','d'],'lambda_formal_not_corepoint':True}
 (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n');inp=out/'native-input.txt'
 with inp.open('w')as f:
  for mat in[A,B]:
   for row in mat:
    for p in row:f.write(ar.line(p))
  for p in[U,G,R0,R4]:f.write(ar.line(p))
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_norm_pencil_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-norm-pencil'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','normpencil-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');C=list(map(ar.poly,n['pencil_coefficients_lambda0_to4']));assert len(C)==5;ar.check([(R0,{(4,0,0):Q(5000**4)})],C[0]);assert C[4]==R4
 RR,a,e,d=ring('a,e,d',s.QQ);rp=lambda p:RR.from_dict({m:s.QQ(str(v))for m,v in p.items()});prev=G;mixed=[]
 for i,rec in enumerate(n['mixed3_records']):
  primitive=ar.poly(rec['primitive']);ar.check([(primitive,{(0,0,0):Q(rec['QQ_content'])})],C[i+1]);current=primitive
  for step in rec['localizer_removal_steps']:
   t=ar.poly(step['removed_gcd_with_U']);after=ar.poly(step['remaining_after']);ar.check([(t,after)],current);ar.check([(t,ar.poly(step['U_quotient']))],U);current=after
  h=ar.poly(rec['remaining_h_coprime_U']);assert current==h;mixed.append(h)
  if h:cg=rp(U).gcd(rp(h));assert cg/cg.LC==RR.one
  common=ar.poly(rec['curve_gcd_after_this']);ar.check([(common,ar.poly(rec['previous_curve_quotient_by_gcd']))],prev);ar.check([(common,ar.poly(rec['mixed_quotient_by_gcd']))],h);other=rp(h).gcd(rp(prev));assert other/other.LC==rp(common);prev=common
 final=ar.poly(n['final_common_curve_gcd']);assert prev==final;ar.check([(final,ar.poly(n['original_G_quotient_by_final']))],G)
 for h,q in zip(mixed,n['mixed3_quotients_by_final']):ar.check([(final,ar.poly(q))],h)
 finite=final=={(0,0,0):Q(1)};classification='entire_B_bad_finite_or_empty'if finite else'only_smaller_necessary_common_projection_curve_plus_finite_residual'
 proof={'classification':classification,'norm_pencil_root_sharing_statement':'For each algebraicallyclosed fieldpoint with a nonzero, H is monic4. det(multiplicationby(g+lambda*k))=product_r(g(r)+lambda*k(r)), countingrootswithmultiplicity. All5coefficients vanish iff at leastoneroot hasg(r)=k(r)=0. AllactualBbad satisfies this, but converse needs allotherminors.','full_cofactors_not_dropped':'Bothsource monicdivision4column identities are retained; commonclear scalar uses onlyaandQQunits.','independent_integer_endpoint_and_ALL_unit_strip_gcd_quotient_reconstructions':True,'independent_reverse_SymPy_ALL3_curve_gcd_and_U_coprimality_matches':True,'all5_coefficients_source_bound_exact_native':True,'if_gcd1_finiteness':'Oldnecessaryprojection V(h,K)=V(G) union finiteZ. OnV(G), allthree mixedcoefficients vanish atBbad. gcd(G,mixed1,mixed2,mixed3)=1 inQQ[a,e] excludesallheight1primes, so thissupportis0D/empty. MonicH makes everyBbad finite/empty. No existence/emptinessinference.','if_gcd_nonconstant':'Retainnewcommoncurve AND the finite residual ofdividedG/mixedcoeffs; also retain oldfinite residual. A survivingnecessarycurve doesnotassertactualbadexists.','no_GB_nested_resultant_corepoint_prime_or_globalrank3_inference':True}
 (out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),gb/'report.json',gp,kb/'report.json',kp,nb0/'report.json',p0,nb1/'report.json',p1,sb/'report.json',sup]
 report={'id':ID,'classification':classification,'raw5_coefficient_stats':list(map(ar.stats,C)),'postUstrip_mixed3_stats':list(map(ar.stats,mixed)),'curve_gcd_chain_stats':[ar.stats(ar.poly(x['curve_gcd_after_this']))for x in n['mixed3_records']],'final_common_curve_stats':ar.stats(final),'all3_unit_removal_counts':[len(x['localizer_removal_steps'])for x in n['mixed3_records']],'native_seconds':n['native_seconds'],'peak_single_native_support':n['peak_single_polynomial_support'],'all_native_and_independent_integer_gcd_checks':True,'GB_nested_resultant_newpoint_prime_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
