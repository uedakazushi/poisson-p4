#!/usr/bin/env python3
"""ONE small QQ4 ideal GB/defaultM2 with cached lift; no larger retry or new point."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_SMALL_IDEAL_GB_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def matrixfile(p):
 ls=p.read_text().splitlines();nr,nc=map(int,ls[0].split());out=[]
 for line in ls[1:]:
  t=line.split();n=int(t[0]);assert len(t)==1+5*n;out.append({tuple(map(int,t[k+1:k+5])):Q(t[k])for k in range(1,len(t),5)})
 assert len(out)==nr*nc;return nr,nc,out
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';ip=base/'certificate.json';pr=read(base/'report.json');assert sha(ip)==pr['evidence_sha256'][ip.name];c=read(ip)
 gens=[[(tuple(m)+(0,),Q(v))for m,v in p]for p in c['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']]
 loc=[(tuple(m)+(1,),Q(v))for m,v in c['FULLU_product']]+[((0,0,0,0),Q(-1))];gens.append(loc);assert len(gens)==23 and sum(len(p)for p in gens)==3013
 names=['a','e','d','z']
 def expr(p):
  return '+'.join('('+str(v)+')'+''.join('*'+names[i]+('^'+str(e)if e!=1 else'')for i,e in enumerate(m)if e)for m,v in p)or'0'
 expected=[dict(p)for p in gens];source_text='matrix{{'+','.join(expr(p)for p in gens)+'}}'
 version=subprocess.check_output(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec M2 --version'],text=True,timeout=5);(out/'M2-version.txt').write_text(version)
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_small_ideal_identity_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-small-ideal-identity';command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=10)
 exesha=sha(exe)
 # Proven listForm native export, with actual-string byte counting and only default GB.
 helper='''outputBytes=0;
emitPoly=(stream,pol)->(
 ts:=listForm pol; piece:=toString(#ts);outputBytes=outputBytes+#piece+1;
 assert(outputBytes<30000000);stream<<piece;
 scan(ts,term->(piece=concatenate(join({" ",toString(last term)},flatten apply(first term,ee->{" ",toString ee})));outputBytes=outputBytes+#piece;assert(outputBytes<30000000);stream<<piece;));stream<<endl;);
emitMatrix=(mat,filename)->(stream:=openOut filename;stream<<numRows mat<<" "<<numColumns mat<<endl;scan(numRows mat,i->scan(numColumns mat,j->emitPoly(stream,mat_(i,j))));close stream;);
phase=(text)->(stream:=openOut "'''+str(out/'phase.txt')+'''";stream<<text<<endl;close stream;print text;<<flush;);
'''
 script='R=QQ[a,e,d,z,MonomialOrder=>GRevLex];\n'+helper+'M='+source_text+';\nassert(numRows M==1 and numColumns M==23);\nemitMatrix(M,"'+str(out/'source23.native')+'");phase("GB_START");\nG=gb(M,ChangeMatrix=>true);B=gens G;phase("GB_COMPLETE");\nemitMatrix(B,"'+str(out/'basis.native')+'");phase("BASIS_SAVED");\nassert(M%G==0*M);\nidealIsUnit=(1_R%G==0_R);\nif idealIsUnit then (\nL=matrix{{1_R}}//G;assert(numRows L==23 and numColumns L==1);assert(M*L==matrix{{1_R}});emitMatrix(L,"'+str(out/'unit-lift.native')+'");phase("UNIT_LIFT_SAVED_VERIFIED");\n) else (\nT=B//G;assert(numRows T==23);assert(M*T==B);emitMatrix(T,"'+str(out/'basis-change.native')+'");phase("PROPER_BASIS_AND_CHANGE_SAVED_VERIFIED");\n);\n'
 m2path=out/'input.m2';m2path.write_text(script);timedout=False
 with(out/'M2.log').open('w')as f:
  try:run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec M2 --script "$1"','B-small-GB',str(m2path)],stdout=f,stderr=subprocess.STDOUT,timeout=max(1,52-(time.monotonic()-start)));exitcode=run.returncode
  except subprocess.TimeoutExpired:timedout=True;exitcode=124
 phase=(out/'phase.txt').read_text().strip()if(out/'phase.txt').exists()else'not_started';basis=None;unit_candidate=False;verified=False;native_audit=None;dimension=None
 if(out/'source23.native').exists():
  nr,nc,got=matrixfile(out/'source23.native');assert(nr,nc)==(1,23)and got==expected
 if(out/'basis.native').exists()and phase in['BASIS_SAVED','UNIT_LIFT_SAVED_VERIFIED','PROPER_BASIS_AND_CHANGE_SAVED_VERIFIED']:
  nr,nc,basis=matrixfile(out/'basis.native');assert nr==1;unit_candidate=any(set(p)=={(0,0,0,0)}and p[(0,0,0,0)]for p in basis)
  if not unit_candidate:
   leads=[max(p,key=lambda m:(sum(m),tuple(-m[i]for i in reversed(range(4)))))for p in basis if p]
   height=min(mask.bit_count()for mask in range(16)if all(any(mask>>i&1 and m[i]for i in range(4))for m in leads));dimension=4-height
 if phase in['UNIT_LIFT_SAVED_VERIFIED','PROPER_BASIS_AND_CHANGE_SAVED_VERIFIED']and time.monotonic()-start<55:
  lift=out/('unit-lift.native'if unit_candidate else'basis-change.native');target=out/'basis.native'
  if unit_candidate:target=out/'target-one.native';target.write_text('1 1\n1 1 0 0 0 0\n')
  with(out/'identity-check.log').open('w')as f:verify=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','B-small-identity',str(exe),str(out/'source23.native'),str(lift),str(target),str(out/'identity-check.json')],stdout=f,stderr=subprocess.STDOUT,timeout=max(1,57-(time.monotonic()-start)))
  assert verify.returncode==0;native_audit=read(out/'identity-check.json');assert native_audit['all_passed'];verified=True
 if timedout:classification='blocked_timecap'
 elif unit_candidate and verified:classification='unit_ideal_B_rank3_on_ALL_FULLU'
 elif not unit_candidate and basis is not None and verified:classification='proper_exact_QQ_rankbadlocus_GB_primary'
 else:classification='incomplete_or_unverified_solver_result'
 recipe={'ring':'QQ[a,e,d,z]','order':'GRevLex, variables in displayed order, degrees all1','algorithm':'Macaulay2 DEFAULTgb(M,ChangeMatrix=>true), no AlgorithmLinearAlgebra','one_GB_object_reused':'G; unit1//G or basisB//G uses the cached change matrix, not a fresh GB call',
  'independent_GB_verification':'Parse basis.native; recompute leading monomials by total degree then reverselex. Reduce every S-pair and each of the23 original generators to zero using the fixed basis (no GB regeneration). Check source23.native*saved_change=basis.native, or source23*unitlift=1, with exactQQ4 arithmetic; the latter is independently checked here by FLINT.',
  'localization_interpretation':'z*U-1 identifies the quotient with the actual22-minor ideal localized at U=adFBCJ. Unitlift gives1 in that ideal, hence B rank3 everywhereFULLU. A proper completebasis describes the exact rankbad scheme overQQ; no primary decomposition or point selected.',
  'no_fullH7_full239_or_kernel_count_change_from_Balone':True}
 (out/'verification-recipe.json').write_text(json.dumps(recipe,indent=2)+'\n');exe.unlink(missing_ok=True)
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip]
 report={'id':ID,'classification':classification,'phase':phase,'M2_exit_code':exitcode,'timed_out':timedout,'M2_version':version.strip(),'precise_ring_order':'QQ[a,e,d,z], GRevLex, standarddegree1','source_generator_count':23,'source_terms':3013,
  'basis_polynomial_count':None if basis is None else len(basis),'basis_total_terms':None if basis is None else sum(map(len,basis)),'proper_quotient_dimension_from_leading_monomial_ideal':dimension,'independent_FLINT_membership_verified':verified,'identity_audit':native_audit,
  'GB_calls_in_script':1,'no_newpoint_prime_or_automatic_retry':True,'native_executable_sha256_not_committed':exesha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
