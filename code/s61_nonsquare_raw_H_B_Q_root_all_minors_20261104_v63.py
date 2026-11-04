#!/usr/bin/env python3
"""ONE exact ALL22 source-minor evaluation on Q commonroot, no fieldcancel/GB/points."""
import hashlib,importlib.util,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_Q_ROOT_ALL_MINORS_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 def pin(base,file):
  p=base/file;assert sha(p)==load(base/'report.json')['evidence_sha256'][file];return p,load(p)
 sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';sp,c=pin(sb,'certificate.json');A=list(map(ar.poly,c['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']));U=ar.poly(c['FULLU_product'])
 fb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_COMMON_CURVE_FACTOR_20261104_v63';fp,fc=pin(fb,'certificate.json');F=ar.poly(fc['factors'][0]['factor']);assert len(F)==9 and max(m[0]for m in F)==2
 hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';hp,hc=pin(hb,'certificate.json');Hclear=ar.shift(ar.poly(hc['monic_H_laurent']),(5,0,0))
 qb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_COMPONENT_EXACT_QUOTIENT_Q_20261104_v63';qp,qc=pin(qb,'certificate.json');root=qc['final']['common_gcd'];assert len(root)==2 and root[0]['d_power']==1 and root[1]['d_power']==0
 N={(1,0,0):Q(18000),(0,2,0):Q(6480),(0,1,0):Q(3420),(0,0,0):Q(-980)};D={(0,3,0):Q(486),(0,2,0):Q(243),(0,1,0):Q(-297),(0,0,0):Q(-54)}
 for rec in root[1]['coefficient']:
  p=rec['fraction'];nu={(0,m,0):Q(v)for m,v in p['numerator_e']};den={(0,m,0):Q(v)for m,v in p['denominator_e']};Nc={(0,m[1],0):v for m,v in N.items()if m[0]==rec['a_power']};ar.check([(D,nu)],{})if not Nc else ar.check([(D,nu),(ar.scalar(Nc,-1),den)],{})
 maxAi=max(max(m[0]for m in p)for p in A);maxEi=max(max(m[1]for m in p)for p in A);maxdi=max(max(m[2]for m in p)for p in A);pre={'all22_source_total_terms':sum(map(len,A)),'max_source_degrees_a_e_d':[maxAi,maxEi,maxdi],'root_minus_N_over_D':True,'cleared_eval_term_box_upper':(maxAi+maxdi+1)*(maxEi+4*maxdi+2*maxAi+1),'reduced_term_box_upper':2*(maxEi+4*maxdi+2*maxAi+1),'all22_reduced_support_upper':22*2*(maxEi+4*maxdi+2*maxAi+1)};assert pre['cleared_eval_term_box_upper']<150000;(out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  for p in[F,N,D,U,Hclear]:f.write(ar.line(p))
  f.write('22\n')
  for p in A:f.write(ar.line(p))
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_Q_root_all_minors_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-Q-root-all-minors'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','Q-root-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json')
 def multiply(p,q):
  pi,dp=ar.intp(p);qi,dq=ar.intp(q);r={}
  for m,c in pi.items():
   for nn,d in qi.items():
    z=tuple(x+y for x,y in zip(m,nn));r[z]=r.get(z,0)+c*d
  return{m:Q(v,dp*dq)for m,v in r.items()if v}
 one={(0,0,0):Q(1)};NP=[one];DP=[one]
 for j in range(7):NP.append(multiply(NP[-1],ar.scalar(N,-1)));DP.append(multiply(DP[-1],D))
 def evaluate(p):
  m=max(x[2]for x in p);v={}
  for j in range(m+1):pj={(x[0],x[1],0):c for x,c in p.items()if x[2]==j};v=ar.add(v,multiply(pj,multiply(NP[j],DP[m-j])))
  return v
 def verify(p,rec):
  value=ar.poly(rec['cleared_evaluation']);assert value==evaluate(p);q=ar.poly(rec['Q_quotient']);r=ar.poly(rec['remainder']);ar.check([(q,F),(r,one)],value);norm=ar.poly(rec['remainder_quadratic_norm_e']);conj=ar.poly(rec['remainder_conjugate']);normq=ar.poly(rec['norm_membership_Q_quotient']);ar.check([(r,conj),(ar.scalar(normq,-1),F)],norm);assert all(m[0]==0 and m[2]==0 for m in norm);return r,norm
 rows=[verify(p,rec)for p,rec in zip(A,n['all22_evaluations'])];assert len(rows)==22;ur,unorm=verify(U,n['localizer_U_evaluation']);hr,hnorm=verify(Hclear,n['Hclear_evaluation']);assert not hr and ur and unorm
 nonzero=[i for i,(p,z)in enumerate(rows)if p];assert nonzero==n['nonzero_remainder_indices'];allzero=not nonzero;classification='actual_B_rankbad_curve_on_FULLU_proved'if allzero else'generic_Q_candidate_excluded_B_bad_entirely_finite_or_empty'
 cert={'classification':classification,'Q':fc['factors'][0]['factor'],'root_numerator_N':[[list(m),str(v)]for m,v in N.items()],'root_denominator_D':[[list(m),str(v)]for m,v in D.items()],'source_Q_field_commonroot_matches_minus_N_D':True,'all22_actual_source_minors_indexed_by_B35_certificate':True,'independent_integer_all22_cleared_evaluations_and_Q_memberships':True,'independent_integer_ALL_quadratic_norm_conjugate_identities':True,'Hclear_remainder_zero':True,'U_remainder_and_quadratic_norm_nonzero':True,'nonzero_minor_remainder_indices':nonzero,'finite_e_exceptions':['D(e)=0','a=0 excludedbyFULLU','zerosofNorm(Uclear)(e)','previousQfielddenominatorpolynomials'],'if_allzero':'The irreduciblecurveQ=0 has a denseopen where a,D,NormU nonzero. Map(a,e)->(a,e,d=-N/D) hasH=0 andALL22minors0, U nonzero. B35sourceknownunitrelations then giveALL35actualBmaximalminors0, so a genuineFULLUB-rankbadcurve. No fullcubicrankdrop follows.','if_nonzero':'OnQoffthe savedfielddenominators the commonroot isunique. A nonzeroremainder hasnonzeroquadraticnorm ine, so anyactualbadpoint there lies overitsfinitelymanyzeros. TogetherwithLboundaryandoldfiniteprojectionresidual,this makesALLBbadfinite/empty. No empty/globalrank3claim.','no_new_GB_point_prime_nestedresultant_or_fullcubicquartic_rank_inference':True}
 (out/'independent-verification.json').write_text(json.dumps(cert,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),sb/'report.json',sp,fb/'report.json',fp,hb/'report.json',hp,qb/'report.json',qp];report={'id':ID,'classification':classification,'all22_remainders_zero':allzero,'nonzero_indices':nonzero,'all22_remainder_term_counts':[len(p)for p,z in rows],'U_remainder_stats':ar.stats(ur),'U_norm_stats':ar.stats(unorm),'native_seconds':n['native_seconds'],'peak_native_support':n['peak_single_polynomial_support'],'all_independent_integer_evaluations_memberships_norms_checked':True,'GB_point_prime_nestedresultant_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
