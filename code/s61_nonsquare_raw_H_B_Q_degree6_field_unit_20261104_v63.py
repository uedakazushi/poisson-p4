#!/usr/bin/env python3
"""ONEdegree6component unit-certified K-algebra8x176, independent fullcofactor proof."""
import hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_Q_DEGREE6_FIELD_UNIT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True)
 helper=ROOT/'code/s61_nonsquare_raw_H_B_pair_norm_20261104_v63.py';spec=importlib.util.spec_from_file_location('integerhelpers',helper);ar=importlib.util.module_from_spec(spec);spec.loader.exec_module(ar)
 def pin(stem,file='certificate.json'):
  b=ROOT/'reports'/stem;p=b/file;assert sha(p)==load(b/'report.json')['evidence_sha256'][file];return b,p,load(p)
 fb,fp,fc=pin('S61_NONSQUARE_RAW_H_B_Q_EXCEPTION_FACTOR_20261104_v63');factor=fc['factor_inventory'][5]['factor'];field=[[m[1],v]for m,v in factor];assert max(m for m,v in field)==6
 qb,qp,qc=pin('S61_NONSQUARE_RAW_H_B_COMMON_CURVE_FACTOR_20261104_v63');Qpoly=ar.poly(qc['factors'][0]['factor']);hb,hp,hc=pin('S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63');Hclear=ar.shift(ar.poly(hc['monic_H_laurent']),(5,0,0));sb,sp,sc=pin('S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63');Pi=list(map(ar.poly,sc['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling']));U=ar.poly(sc['FULLU_product'])
 pre={'selected_factor_index':5,'field_degree':6,'K_algebra_dimension':8,'full_ideal_matrix_shape':[8,176],'rational_coefficient_count_upper':8*176*6,'q0_degree_e':4,'q0_unit_check':'QQ[e]xgcdagainstf6withliteralBezout1','all_pivot_inverses':'sameQQ[e]xgcdunitcertificates; nonunitSTOP','actual_model_serialized_height_guard_before_solve':8*1024*1024,'no_QQ48x1056_matrix':True};(out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 inp=out/'native-input.txt'
 with inp.open('w')as f:
  f.write(str(len(field))+''.join(' '+str(m)+' '+str(v)for m,v in field)+'\n')
  for p in[Qpoly,Hclear,U]:f.write(ar.line(p))
  f.write('22\n');f.writelines(ar.line(p)for p in Pi)
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_Q_degree6_field_unit_20261104_v63.cpp';exe=ROOT/'work/raw-H-B-Q-degree6-field-unit'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','Q6field-native',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 n=load(out/'native-certificate.json');ep=lambda p:{int(m):Q(v)for m,v in p if Q(v)};up=lambda p:{(0,int(m),0):Q(v)for m,v in p if Q(v)};one3={(0,0,0):Q(1)};mod=ep(field);en=6
 # Independent polynomial quotient arithmetic, no field inverse/GCD calls.
 def ea(p,q):
  r=dict(p)
  for m,c in q.items():
   v=r.get(m,Q())+c
   if v:r[m]=v
   else:r.pop(m,None)
  return r
 def es(p,c):return{m:v*c for m,v in p.items()if v*c}
 def er(p):
  r=dict(p)
  while r and max(r)>=6:
   j=max(r);c=r[j];r=ea(r,{j-6+k:-c*v for k,v in mod.items()})
  return r
 def em(p,q):
  if not p or not q:return{}
  da=math.lcm(*(x.denominator for x in p.values()));db=math.lcm(*(x.denominator for x in q.values()));r={}
  for m,c in p.items():
   for j,v in q.items():r[m+j]=r.get(m+j,0)+int(c*da)*int(v*db)
  return er({m:Q(v,da*db)for m,v in r.items()if v})
 o={0:Q(1)}
 for cert in n['all_field_inverse_certificates']:
  element=ep(cert['element']);iv=ep(cert['inverse']);ar.check([(up(cert['Bezout_modulus']),up(field)),(up(cert['Bezout_element']),up(cert['element']))],one3);assert em(element,iv)==o
 ee=[o]
 for j in range(20):ee.append(em(ee[-1],{1:Q(1)}))
 q0={};q1={}
 for m,c in Qpoly.items():
  if m[0]==0:q0=ea(q0,es(ee[m[1]],c))
  if m[0]==1:q1=ea(q1,es(ee[m[1]],c))
 assert q0==ep(n['q0_reduced'])and q1==ep(n['q1_reduced'])
 bz=lambda:[{},{}]
 def bo():return[o,{}]
 def ba(p,q):return[ea(p[0],q[0]),ea(p[1],q[1])]
 def bn(p):return[es(p[0],-1),es(p[1],-1)]
 def bs(p,c):return[em(p[0],c),em(p[1],c)]
 def bm(p,q):
  vv=em(p[1],q[1]);return[ea(em(p[0],q[0]),es(em(q0,vv),-1)),ea(ea(em(p[0],q[1]),em(p[1],q[0])),es(em(q1,vv),-1))]
 ag=[{},o];ainv=list(map(ep,n['a_inverse']));assert bm(ag,ainv)==bo();ap=[bo()]
 for j in range(20):ap.append(bm(ap[-1],ag))
 def hcoef(j):
  r=bz()
  for m,c in Hclear.items():
   if m[2]==j:r=ba(r,bs(ap[m[0]],es(ee[m[1]],c)))
  return r
 hclear=[hcoef(j)for j in range(5)];assert hclear==[[ep(p)for p in b]for b in n['Hclear_coefficients']];iv5=bo()
 for j in range(5):iv5=bm(iv5,ainv)
 h=[bm(b,iv5)for b in hclear];assert h==[[ep(p)for p in b]for b in n['H_monic_coefficients']]and h[4]==bo()
 az=lambda:[bz()for j in range(4)]
 def ao():r=az();r[0]=bo();return r
 def aa(p,q):return[ba(p[j],q[j])for j in range(4)]
 def am(p,q):
  r=[bz()for j in range(7)]
  for i in range(4):
   for j in range(4):r[i+j]=ba(r[i+j],bm(p[i],q[j]))
  for j in range(6,3,-1):
   v=r[j]
   for i in range(4):r[j-4+i]=ba(r[j-4+i],bn(bm(v,h[i])))
   r[j]=bz()
  return r[:4]
 def ab(p,b):return[bm(v,b)for v in p]
 dp=[ao()];dg=az();dg[1]=bo()
 for j in range(20):dp.append(am(dp[-1],dg))
 def evaluate(p):
  r=az()
  for m,c in p.items():r=aa(r,ab(dp[m[2]],bs(ap[m[0]],es(ee[m[1]],c))))
  return r
 images=[evaluate(p)for p in Pi];nativeF=lambda f:[[ep(p)for p in b]for b in f];assert images==list(map(nativeF,n['all22_original_source_images']));ui=evaluate(U);assert ui==nativeF(n['U_image']);target=ao();pw=ui;power=8
 while power:
  if power&1:target=am(target,pw)
  power>>=1
  if power:pw=am(pw,pw)
 assert target==nativeF(n['target_U_power8']);solved=n['membership_solved']
 if solved:
  cofs=list(map(nativeF,n['source22_cofactor_witness']));assert len(cofs)==22;dot=az()
  for i in range(22):dot=aa(dot,am(cofs[i],images[i]))
  assert dot==target
 proof={'factor_index':5,'field_factor':factor,'unit_inverse_all_Bezout1_independent_integer_checked':True,'independent_sourceQ_H_and_ALL22_images_checked':True,'independent_U8_and_original22cofactor_dot_checked':solved,'all_used_inverses_are_actual_units_even_without_field_assumption':True,'membership_solved':solved,'basis':'a^i*d^j overQQ[e]/f6, index2*j+i, i<2,j<4','closed_component_scope':'Entiref6=Q=H=0 withaunit; q0inversemakesaunitautomatically. Sixconjugateeembeddingsallcovered, notsampling.','no_other_Q_component_or_old_residual_solved':True,'no_point_prime_GB_resultant':True};(out/'independent-verification.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True)
 paths=[Path(__file__),cpp,helper,ROOT/'research/notes'/(ID+'.md'),fb/'report.json',fp,qb/'report.json',qp,hb/'report.json',hp,sb/'report.json',sp]
 report={'id':ID,'selected_factor_index':5,'field_degree':6,'K_algebra_dimension':8,'QQ_dimension48_matrix_NOT_built':True,'field_ideal_matrix_shape':[8,176],'matrix_rational_coefficient_count':n['matrix_rational_coefficient_count'],'matrix_max_numerator_plus_denominator_bits':n['matrix_max_numerator_plus_denominator_bits'],'rank_over_K':n['rank'],'whole_selected_component_FULLU_empty_proved':solved,'original_source_cofactor_count':len(n['source22_cofactor_witness']),'all_field_inverse_certificate_count':len(n['all_field_inverse_certificates']),'all_independent_source_and_inverse_and_cofactor_checks':True,'native_seconds':n['native_seconds'],'other_component_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
