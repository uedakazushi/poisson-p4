#!/usr/bin/env python3
"""ONE exact squarefree exception job; oldh1 LC-only preflight, no fiber solve."""
import hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_EXCEPTION_SQUAREFREE_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 def pin(stem,file='certificate.json'):
  b=ROOT/'reports'/stem;p=b/file;assert sha(p)==load(b/'report.json')['evidence_sha256'][file];return b,p,load(p)
 cb,cp,c=pin('S61_NONSQUARE_RAW_H_B_FINITE_BAD_LOCUS_20261104_v63');groups=[r['e_polynomial_factors']for r in c['all_B_bad_contained_in_three_covers'][1:]];pre=[]
 for factors in groups:
  bit=0;de=0
  for p in factors:
   co=[Q(v)for m,v in p];den=math.lcm(*(v.denominator for v in co));z=[int(v*den)for v in co];con=math.gcd(*z);z=[x//con for x in z];bit+=sum(abs(x)for x in z).bit_length();de+=max(m for m,v in p)
  pre.append({'factor_count':len(factors),'product_degree':de,'primitive_integer_product_coeff_bits_upper':bit,'single_product_serialization_bound':(de+1)*(150+(bit*30103+99999)//100000)})
 assert [x['product_degree']for x in pre]==[65,361] and sum(x['single_product_serialization_bound']for x in pre)<32*1024*1024;(out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  f.write('2\n')
  for factors in groups:
   f.write(str(len(factors))+'\n')
   for p in factors:f.write(str(len(p))+''.join(' '+str(m)+' '+str(v)for m,v in p)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_exception_squarefree_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-exception-squarefree'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','squarefree-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');up=lambda p:{(0,int(m),0):Q(v)for m,v in p if Q(v)};one={(0,0,0):Q(1)}
 def multiply(p,q):
  a,da=ar.intp(p);b,db=ar.intp(q);r={}
  for m,c in a.items():
   for nn,d in b.items():
    z=tuple(x+y for x,y in zip(m,nn));r[z]=r.get(z,0)+c*d
  return{m:Q(v,da*db)for m,v in r.items()if v}
 squarefree=[]
 for source,rec in zip(groups,n['groups']):
  P=one
  for src,mon,lc in zip(source,rec['normalized_input_factors'],rec['original_input_leading_constants']):
   ar.check([(up(mon),{(0,0,0):Q(lc)})],up(src));P=multiply(P,up(mon))
  assert P==up(rec['monic_product_P']);der={(0,m[1]-1,0):m[1]*v for m,v in P.items()if m[1]};assert der==up(rec['derivative_Pprime']);G=up(rec['derivative_gcd_G']);S=up(rec['squarefree_S_P_over_G']);ar.check([(up(rec['P_Pprime_Bezout'][0]),P),(up(rec['P_Pprime_Bezout'][1]),der)],G);ar.check([(G,S)],P);ar.check([(G,up(rec['Pprime_quotient_by_G']))],der);sd={(0,m[1]-1,0):m[1]*v for m,v in S.items()if m[1]};assert sd==up(rec['S_derivative']);ar.check([(up(rec['S_Sprime_Bezout1'][0]),S),(up(rec['S_Sprime_Bezout1'][1]),sd)],one)
  remaining=P
  for step in rec['root_support_removal_steps']:
   factor=up(step['removed_factor']);after=up(step['remaining_after']);ar.check([(factor,after)],remaining);ar.check([(factor,up(step['S_quotient']))],S);remaining=after
  assert remaining==one;squarefree.append(S)
 rb,rp,rc=pin('S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63');h1=ar.poly(rc['oldh_quotient_h1']);k1=ar.poly(rc['newK_quotient_k1'])
 old=[]
 for label,p in [('h1',h1),('k1',k1)]:
  ed=max(m[1]for m in p);lc={m:v for m,v in p.items()if m[1]==ed};assert len(lc)==1 and (0,ed,0)in lc;v=lc[(0,ed,0)];mon=ar.scalar(p,1/v);ar.check([(mon,{(0,0,0):v})],p);old.append({'label':label,'e_degree':ed,'literal_e_leading_piece':[[list(m),str(c)]for m,c in lc.items()],'LC_is_QQconstant_global_unit':str(v),'monic_normalized':[[list(m),str(c)]for m,c in sorted(mon.items(),reverse=True)]})
 fb,fp,fc=pin('S61_NONSQUARE_RAW_H_B_COMMON_CURVE_FACTOR_20261104_v63');degs=[max(m[1]for m in p)for p in squarefree];bounds=[4*degs[0],8*degs[1]];pointbound=1600+sum(bounds)
 proof={'squarefree_L':n['groups'][0]['squarefree_S_P_over_G'],'squarefree_Q':n['groups'][1]['squarefree_S_P_over_G'],'all_source_products_and_both_xgcd_identities_independent_integer_verified':True,'all_iterated_support_factors_divide_S_and_reconstruct_P':True,'zero_sets_unchanged_exact_support_certificate':True,'old_residual_LC_preflight':old,'old_h1_e16_monic_and_H4_give_free64_over_QQ_aunit':True,'source_finite_cover_inputs':{'L':fc['factors'][1]['factor'],'Q':fc['factors'][0]['factor'],'H':c['H_monic_d4_on_aunit'],'squarefree_e_polynomials':[n['groups'][0]['squarefree_S_P_over_G'],n['groups'][1]['squarefree_S_P_over_G']]},'finite_QQ_cover_dimension_upper_after_a_localization_and_H4':bounds,'total_geometric_B_bad_bound':pointbound,'geometric_point_bound_NOT_scheme_length':True,'no_fiber_unit_test_resultant_GB_newpoint_or_prime_job':True,'localization_scope':'OverQQ finitealgebralocalization is a directfactor and staysfinite; thisdoesNOT supplyanintegralpropermodel/noescape forfuturemodpNakayama.'}
 (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),cb/'report.json',cp,rb/'report.json',rp,fb/'report.json',fp]
 report={'id':ID,'original_exception_degrees':[65,361],'squarefree_exception_degrees':degs,'root_support_removal_step_counts':[len(x['root_support_removal_steps'])for x in n['groups']],'L_Q_finite_QQ_algebra_dimension_upper':bounds,'old_residual_point_bound':1600,'updated_geometric_B_bad_point_bound':pointbound,'old_h1_LC_global_QQconstant_confirmed':True,'native_seconds':n['native_seconds'],'all_independent_integer_reconstructions_checked':True,'fiber_tests_GB_resultant_newpoint_prime_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
