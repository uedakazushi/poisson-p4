#!/usr/bin/env python3
"""Independent fixednative det + owninteger source/division/gcd/Bezout review; no authorimport."""
import hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def pd(p):return{tuple(m):Q(v)for m,v in p if Q(v)}
def up(p):return{(int(m),0,0):int(v)for m,v in p if int(v)}
def integerize(p):
 den=math.lcm(*(Q(v).denominator for v in p.values()));return {m:int(Q(v)*den)for m,v in p.items()},den
def check(products,expected):
 terms=[]
 for p,q in products:
  a,da=integerize(p);b,db=integerize(q);terms.append((a,b,da*db))
 ep,de=integerize(expected);den=math.lcm(de,*(x[2]for x in terms));out={m:-v*(den//de)for m,v in ep.items()}
 for p,q,d in terms:
  scale=den//d
  for m,c in p.items():
   for n,v in q.items():
    z=tuple(x+y for x,y in zip(m,n));r=out.get(z,0)+c*v*scale
    if r:out[z]=r
    else:out.pop(z,None)
 assert not out
 return True
def shift(p,v):return{tuple(x+y for x,y in zip(m,v)):c for m,c in p.items()}
def add(p,q):
 r=dict(p)
 for m,c in q.items():
  v=r.get(m,Q())+c
  if v:r[m]=v
  else:r.pop(m,None)
 return r
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True)
 b=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_20261104_v63';br=load(b/'report.json');assert br['all_passed']
 for section in['input_sha256','evidence_sha256']:
  for rel,h in br[section].items():assert sha((ROOT/rel)if section=='input_sha256'else(b/rel))==h,rel
 n=load(b/'native-certificate.json');assert not n['determinant_zero'];ib=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MIXED_NORM16_INPUT_20261104_v63';ir=load(ib/'report.json');c=load(ib/'certificate.json');cost=load(ib/'cost-preflight.json');assert sha(ib/'certificate.json')==ir['evidence_sha256']['certificate.json'];assert sha(ib/'cost-preflight.json')==ir['evidence_sha256']['cost-preflight.json']
 pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63';pp=pb/'native-certificate.json';assert sha(pp)==load(pb/'report.json')['evidence_sha256'][pp.name];pencil=load(pp);G=pd(c['G']);assert G==pd(pencil['final_common_curve_gcd'])
 for i in range(3):
  source=pd(pencil['mixed3_records'][i]['remaining_h_coprime_U']);mi=pd(c['all3_quotients_by_G'][i]);assert source==pd(c['mixed3_stripped'][i])and mi==pd(pencil['mixed3_quotients_by_final'][i]);check([(G,mi)],source)
 selected=c['chosen_mixed_lambda_coefficient_index'];assert selected==1;LC=Q(c['chosen_mi_QQ_leading_unit']);K=pd(c['chosen_mi_monic']);check([(K,{(0,0,0):LC})],pd(c['all3_quotients_by_G'][0]));assert {(m[0],m[2]):v for m,v in K.items()if m[1]==19}=={(0,0):Q(1)}
 H=pd(c['H_h1_monic16']);oldib=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63';oldinput=load(oldib/'certificate.json');assert H==pd(oldinput['h1_e16_monic']);assert {(m[0],m[2]):v for m,v in H.items()if m[1]==16}=={(0,0):Q(1)}
 M=[[pd(p)for p in row]for row in c['matrix_rows_source_order']];q=list(map(pd,c['all16_column_quotients_source_order']));r=list(map(pd,c['all16_column_remainders_source_order']));D=int(cost['matrix_uniform_QQ_denominator']);allN=[]
 for j in range(16):
  rebuilt={}
  for i in range(16):
   assert all(m[1]==m[2]==0 for m in M[i][j]);rebuilt=add(rebuilt,shift(M[i][j],(0,i,0)))
  assert rebuilt==r[j];check([(H,q[j]),(r[j],{(0,0,0):Q(1)})],shift(K,(0,j,0)))
 for row in M:
  for p in row:
   v={}
   for m,c0 in p.items():c1=c0*D;assert c1.denominator==1;v[m[0]]=c1.numerator
   allN.append(v)
 det=up(n['determinant']);nativeinput=out/'independent-native-input.txt'
 with nativeinput.open('w')as f:
  for p in allN+[{m[0]:v for m,v in det.items()}]:
   f.write(str(len(p))+'\n')
   for m,v in sorted(p.items(),reverse=True):f.write(str(m)+' '+str(v)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_mixed_norm16_review_20261104_v63.cpp';exe=ROOT/'work/mixed-norm16-independent-review'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 exesha=sha(exe)
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','indepfixed-det',str(exe),str(nativeinput),str(out/'independent-native-det.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-start)))
 own=load(out/'independent-native-det.json');assert own['independent_exact_det_equals_author']and own['degree']==98 and own['a_valuation']==35 and up(own['determinant'])==det
 prim=up(n['primitive']);new=up(n['new_primitive_without_a']);old=up(n['old_primitive_a117']);g=up(n['gcd']);qn=up(n['quotient_new']);qo=up(n['quotient_old']);one={(0,0,0):1};check([(prim,{(0,0,0):int(n['determinant_QQ_content'])})],det);check([(new,{(35,0,0):1})],prim);assert min(m[0]for m in new)==0
 prior=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63_REVIEW';assert load(prior/'report.json')['all_passed'];assert old==pd(load(prior/'independent-certificate.json')['nonzero_primitive_a117']);check([(g,qn)],new);check([(g,qo)],old);assert n['Bezout_saved']and int(n['Bezout_constant'])!=0;check([(up(n['Bezout_S']),new),(up(n['Bezout_T']),old)],{m:int(n['Bezout_constant'])*v for m,v in g.items()});assert max(m[0]for m in g)==13 and min(m[0]for m in g)==0
 proof={'all_passed':True,'native_algorithm':'independentFLINTfmpz_poly_mat_det_fflu; authorownscalarBareissnotimported','own_all3_Mi_G_source_products':True,'own_all16_monic_divisions_and_sameH1_basis':True,'own_all256_integer_matrix_clearing':True,'independent_exact_fixed_det_equals_author':True,'det_degree':98,'det_a_order':35,'primitive_non_a_degree':63,'old_non_a_degree':117,'own_integer_allcontent_a_removal_gcdquotient_Bezout_identities':True,'actual_gcd_degree':13,'triple_quotient':'T=QQ[a,a^-1,e]/(h1,k1,m1) = coker_R([M_k1 | M_m1]), R=QQ[a,a^-1]','Fitting_bound_proof':'Fitt0(T) contains both16x16determinants, hence theirgcdg13 byexactQQBezout afterknowna/QQunitstrips. OverPIDR, Fitt0(T)=(productSmithinvariants), whichthereforedividesg13; T finiteQQandlength<=13.','triple_quotient_QQ_length_upper':13,'H4_containing_cover_QQ_length_upper':52,'H4_bound_requires_actual_monicH_source':True,'actual_Bbad_transfer_scope':'ForD(G), allbadneedsm1 afterdivisionbyG; G0Q/Lcasesmustremainunlessseparateallpointunitgatesaccepted. ThisreviewaloneNOTglobalBempty.','nonconstant_g13_NOT_badexistence':True,'no_solve_newnormsearch_GB_point_prime':True};(out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n');exe.unlink(missing_ok=True)
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),b/'report.json',b/'native-certificate.json',ib/'report.json',ib/'certificate.json',ib/'cost-preflight.json',pb/'report.json',pp,oldib/'report.json',oldib/'certificate.json',prior/'report.json',prior/'independent-certificate.json'];report={'id':ID,'all_passed':True,'independent_fixed_det_degree98_aorder35':True,'independent_gcd13_certificate_verified':True,'triple_coker_length_upper':13,'H4_cover_length_upper':52,'G0_scope_kept_not_global_empty':True,'native_seconds':own['seconds'],'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'native_executable_SHA_not_committed':exesha,'seconds':time.monotonic()-start,'frozen':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
