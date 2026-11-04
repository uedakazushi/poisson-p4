#!/usr/bin/env python3
"""ONE additional fixed4 norm+adj+H/k membership and exact oldcurve GCD."""
import hashlib,importlib.util,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integernormhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 ib=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_INPUT_20261104_v63';ip=ib/'certificate.json';assert sha(ip)==load(ib/'report.json')['evidence_sha256'][ip.name];c=load(ip);assert c['selected_generator_index']==20;H=ar.poly(c['source_H']);k=ar.poly(c['source_k']);M=[[ar.poly(p)for p in row]for row in c['matrix_rows']];Qs=list(map(ar.poly,c['column_quotients_by_H']));b=c['cost_bounds_only'];aclear=b['matrix_uniform_a_clear_power'];D=Q(b['matrix_uniform_QQ_denominator']);assert aclear==11;N=[[ar.scalar(ar.shift(p,(aclear,0,0)),D)for p in row]for row in M]
 sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';sp=sb/'certificate.json';assert sha(sp)==load(sb/'report.json')['evidence_sha256'][sp.name];U=ar.poly(load(sp)['FULLU_product']);ob=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_20261104_v63';op=ob/'native-certificate.json';assert sha(op)==load(ob/'report.json')['evidence_sha256'][op.name];oldh=ar.poly(load(op)['remaining_h_coprime_U'])
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for row in N:
   for p in row:f.write(ar.line(p))
  for p in [H,k]+Qs:f.write(ar.lline(p))
  f.write(ar.line(U));f.write(ar.line(oldh));f.write(str(aclear)+' '+str(D)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_additional_minor_norm_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-additional-minor-norm'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','additional-norm-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');R=ar.poly(n['det_cleared_matrix_R']);adj=list(map(ar.poly,n['adjugate_first_column']));assert R and max(m[2]for m in R)==0
 for i in range(4):ar.check([(N[i][j],adj[j])for j in range(4)],R if i==0 else{})
 vh=ar.poly(n['literal_membership_H_cofactor']);vk=ar.poly(n['literal_membership_k_cofactor']);ar.check([(vh,H),(vk,k)],R)
 primitive=ar.poly(n['R_primitive']);one={(0,0,0):Q(1)};ar.check([(primitive,{(0,0,0):Q(n['R_QQ_content'])})],R);cur=primitive
 for step in n['localizer_removal_steps']:
  f=ar.poly(step['removed_gcd_with_U']);after=ar.poly(step['remaining_after']);ar.check([(f,after)],cur);ar.check([(f,ar.poly(step['U_quotient']))],U);cur=after
 K=ar.poly(n['remaining_h_coprime_U']);G=ar.poly(n['gcd_newK_oldh']);assert cur==K and ar.poly(n['final_gcd_h_U'])==one;ar.check([(G,ar.poly(n['newK_quotient_by_gcd']))],K);ar.check([(G,ar.poly(n['oldh_quotient_by_gcd']))],oldh)
 RR,a,e,d=ring('a,e,d',s.QQ);rp=lambda p:RR.from_dict({m:s.QQ(str(v))for m,v in p.items()});cg=rp(U).gcd(rp(K));assert cg/cg.LC==RR.one;other=rp(oldh).gcd(rp(K));assert other/other.LC==rp(G)
 common_constant=G==one;classification='entire_B_bad_finite_or_empty'if common_constant else'common_projection_curve_plus_residual_intersections_only_necessary'
 proof={'classification':classification,'source_k_generator_index':20,'independent_integer_all4_matrix_adj_and_H_k_memberships':True,'independent_integer_ALL_removed_factor_Uquotients_and_reconstructions':True,'independent_reverse_SymPy_gcd_newK_U_is1_and_gcd_newK_oldh_matches':True,'K_is_nonzero':True,'gcd_oldh_K_is1':common_constant,'all_bad_projection_in_V_oldh_K':True,'if_gcd1':'Two coprime nonzero QQ[a,e]polynomials define a0D/empty locus; H monicd4 makes allFULLU B-bad points finite/empty. No point existence or emptiness inference.','if_gcd_nonconstant':'As sets V(oldh,K)=V(G) union V(oldh/G,K/G); retain both. The latter finite/empty; V(G) onlynecessary curve, no actualbadexistence.','no_nested_resultant_GB_newpoint_prime_primarydec_or_globalB3_claim':True,'K_raw_stats':ar.stats(R),'K_after_U_strip_stats':ar.stats(K),'gcd_oldh_K_stats':ar.stats(G)}
 (out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),ib/'report.json',ip,ob/'report.json',op,sb/'report.json',sp]
 report={'id':ID,'classification':classification,'K_raw_stats':ar.stats(R),'K_after_U_strip_stats':ar.stats(K),'gcd_oldh_K_stats':ar.stats(G),'adj4_terms':list(map(len,adj)),'membership_H_k_terms':[len(vh),len(vk)],'unit_factor_removal_steps':len(n['localizer_removal_steps']),'native_seconds':n['native_seconds'],'peak_single_native_support':n['peak_single_polynomial_support'],'all_native_and_independent_integer_identities_checked':True,'nested_resultant_GB_newpoint_prime_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
