#!/usr/bin/env python3
"""ONE finiteL QQ ALL22ideal membership U^N, independent algebra identitycheck."""
import hashlib,importlib.util,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_L_FINITE_UNIT_FIXED_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(58);out.mkdir(exist_ok=False,parents=True)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 def pin(stem,file='certificate.json'):
  b=ROOT/'reports'/stem;p=b/file;assert sha(p)==load(b/'report.json')['evidence_sha256'][file];return b,p,load(p)
 sf,sp,sc=pin('S61_NONSQUARE_RAW_H_B_EXCEPTION_SQUAREFREE_20261104_v63');S=sc['squarefree_L'];assert max(m for m,v in S)==9
 hb,hp,hc=pin('S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63');Hclear=ar.shift(ar.poly(hc['monic_H_laurent']),(5,0,0));bb,bp,bc=pin('S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63');Pi=list(map(ar.poly,bc['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']));U=ar.poly(bc['FULLU_product']);aL=[[2,'-9/50'],[1,'-3/50'],[0,'3/25']]
 pre={'source_e_degree':9,'finite_algebra_dimension_before_a_clopen_upper':36,'ideal_column_count_upper':22*36,'matrix_entry_count_upper':22*36*36,'total_original_minor_terms':2963,'only_clopen_removed_a_zero_known_FULLU_excluded':True,'native_coefficient_height_and_actual_serialization_guards_before_solve':True};assert pre['matrix_entry_count_upper']<150000;(out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 def lineE(p):return str(len(p))+''.join(' '+str(m)+' '+str(v)for m,v in p)+'\n'
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  f.write(lineE(S));f.write(lineE(aL));f.write(ar.line(Hclear));f.write(ar.line(U));f.write('22\n');f.writelines(ar.line(p)for p in Pi)
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_L_finite_unit_fixed_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-L-finite-unit'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','Lfinite-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');en=n['base_e_degree'];N=n['dimension_N'];assert N==4*en;up=lambda p:{(0,m,0):Q(v)for m,v in p if Q(v)};one3={(0,0,0):Q(1)};current=up(S)
 for step in n['a_zero_clopen_removal_steps']:
  factor=up(step['removed_factor']);after=up(step['remaining_after']);ar.check([(factor,after)],current);ar.check([(factor,up(step['aL_quotient']))],up(aL));current=after
 assert current==up(n['remaining_e_modulus']);ar.check([(up(n['modulus_aL_Bezout1'][0]),current),(up(n['modulus_aL_Bezout1'][1]),up(aL))],one3)
 # Separate QQ finitealgebra implementation; no FLINT arithmetic or sourcebuilderimport.
 ep=lambda p:{int(m):Q(v)for m,v in p if Q(v)};mod=ep(n['remaining_e_modulus']);assert mod[en]==1
 def ea(p,q):
  r=dict(p)
  for m,v in q.items():
   z=r.get(m,Q())+v
   if z:r[m]=z
   else:r.pop(m,None)
  return r
 def es(p,c):return{m:v*c for m,v in p.items()if v*c}
 def er(p):
  r=dict(p)
  while r and max(r)>=en:
   m=max(r);v=r[m];r=ea(r,{m-en+j:-v*c for j,c in mod.items()})
  return r
 def em(p,q):
  if not p or not q:return{}
  denp=__import__('math').lcm(*(c.denominator for c in p.values()));denq=__import__('math').lcm(*(c.denominator for c in q.values()));r={}
  for m,c in p.items():
   for j,v in q.items():r[m+j]=r.get(m+j,0)+int(c*denp)*int(v*denq)
  return er({m:Q(v,denp*denq)for m,v in r.items()if v})
 o={0:Q(1)};al=er(ep(aL));ainv=ep(n['aL_inverse']);assert em(ainv,al)==o
 ap=[o];ee=[o]
 for j in range(16):ap.append(em(ap[-1],al));ee.append(em(ee[-1],{1:Q(1)}))
 def hco(i):
  r={}
  for m,c in Hclear.items():
   if m[2]==i:r=ea(r,es(em(ap[m[0]],ee[m[1]]),c))
  return r
 hc0=[hco(i)for i in range(5)];assert hc0==list(map(ep,n['Hclear_coefficients']));inv5=o
 for j in range(5):inv5=em(inv5,ainv)
 h=[em(p,inv5)for p in hc0];assert h==list(map(ep,n['H_normalized_monic_coefficients']))and h[4]==o
 az=lambda:[{}for j in range(4)]
 def ao():r=az();r[0]=o;return r
 def aa(p,q):return[ea(p[j],q[j])for j in range(4)]
 def am(p,q):
  r=[{}for j in range(7)]
  for i in range(4):
   for j in range(4):r[i+j]=ea(r[i+j],em(p[i],q[j]))
  for j in range(6,3,-1):
   v=r[j]
   for i in range(4):r[j-4+i]=ea(r[j-4+i],es(em(v,h[i]),-1))
   r[j]={}
  return r[:4]
 dp=[ao()];d=az();d[1]=o
 for j in range(16):dp.append(am(dp[-1],d))
 def evaluate(p):
  r=az()
  for m,c in p.items():
   ev=es(em(ap[m[0]],ee[m[1]]),c);r=aa(r,[em(v,ev)for v in dp[m[2]]])
  return r
 images=[evaluate(p)for p in Pi];assert images==[[ep(p)for p in F]for F in n['all22_original_source_images']];ui=evaluate(U);assert ui==list(map(ep,n['U_image']));target=ao();pw=ui;power=N
 while power:
  if power&1:target=am(target,pw)
  power>>=1
  if power:pw=am(pw,pw)
 assert target==list(map(ep,n['target_U_power_N']));solved=n['membership_solved']
 if solved:
  cof=[[ep(p)for p in F]for F in n['source22_cofactor_witness']];assert len(cof)==22;dot=az()
  for i in range(22):dot=aa(dot,am(cof[i],images[i]))
  assert dot==target
 proof={'membership_solved':solved,'independent_clopen_a_knownunit_reconstruction_and_global_Bezout1':True,'independent_finite_algebra_sourceH_ALL22_sourceimages_and_U_power_checked':True,'independent_original_source22_cofactor_dot_U_power_checked':solved,'if_solved':'U^N in ALL22sourceideal onentirefiniteLcover; afterUinversion unitideal, so NOactualFULLUB-badpoint onLcover, includingALLsavedLdenominatorexceptions.','if_not_solved':'NativeexactfiniteQQmembership returnedfalse; retainfullmatrix/target forindependentnonmembershipreview. No globalB3/fullcubicnegative inferred.','scope_Lcover_only':'Qcover andoldh1/k1finite residual untouched. RawH7p0/full57/full239 unchanged.','basis':'e^i*d^j, indexj*en+i, i<en,j<4','N':N,'no_fiber_point_prime_GB_nestedresultant':True}
 (out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True);paths=[Path(__file__),cpp,helper,ROOT/'reports/S61_NONSQUARE_RAW_H_B_L_FINITE_UNIT_20261104_v63/failure.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_L_FINITE_UNIT_20261104_v63/failed-source.cpp',ROOT/'research/notes/S61_NONSQUARE_RAW_H_B_L_FINITE_UNIT_20261104_v63.md',sf/'report.json',sp,hb/'report.json',hp,bb/'report.json',bp]
 report={'id':ID,'L_cover_FULLU_empty_proved':solved,'remaining_e_degree_after_a_zero_peel':en,'finite_QQ_algebra_dimension_N':N,'matrix_shape':[N,22*N],'matrix_nonzero':n['matrix_nonzero_entries'],'matrix_max_numerator_plus_denominator_bits':n['matrix_max_numerator_plus_denominator_bits'],'original_source_cofactor_count':len(n['source22_cofactor_witness']),'native_seconds':n['native_seconds'],'all_independent_source_and_witness_checks':True,'Q_old_residual_full57_239_NOT_solved':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
