#!/usr/bin/env python3
"""One fixed raw12 cleared Schur defect; known-unit stripping only if nonzero."""
import argparse,hashlib,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_RAW12_DEFECT_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def dec(p):return{tuple(m):Q(v)for m,v in p if Q(v)}
def enc(p):return[[list(m),str(v)]for m,v in sorted(p.items())]
def add(p,q,c=Q(1)):
 r=p.copy()
 for m,v in q.items():
  r[m]=r.get(m,Q())+c*v
  if not r[m]:del r[m]
 assert len(r)<=150000
 return r
def mul(p,q):
 r={}
 for m,v in p.items():
  for n,w in q.items():
   s=tuple(x+y for x,y in zip(m,n));r[s]=r.get(s,Q())+v*w
 assert len(r)<=150000
 return{m:v for m,v in r.items()if v}
def shift(p,s,c=Q(1)):return{tuple(x+y for x,y in zip(m,s)):c*v for m,v in p.items()if c*v}
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';cp=ROOT/'reports/S61_NONSQUARE_THETA_RATIONAL_COORDINATES_20261104_v63';tp=ROOT/'reports/S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63';dp=ROOT/'reports/S61_NONSQUARE_THETA_DELTA_LAMBDA_PROFILES_20261104_v63'
 cpp=Path(__file__).with_suffix('.cpp');paths=[Path(__file__),cpp,gp/'report.json',gp/'certificate.json',cp/'report.json',cp/'certificate.json',tp/'report.json',tp/'certificate.json',dp/'report.json',dp/'certificate.json']
 rd=lambda p:json.loads(p.read_text());gr,g,cr,c,tr,t,dr,d=[rd(p)for p in paths[2:]]
 for r,f in[(gr,paths[3]),(cr,paths[5]),(tr,paths[7]),(dr,paths[9])]:assert sha(f)==r['certificate_sha256']
 B=[[dec(p)for p in row]for row in g['actual_block7x3']];Delta=dec(g['nonzero_selected_minor_Delta']);lam=[dec(p)for p in g['cleared_upper_row_coordinates7xr'][4]];piv=g['selected_rows7'];assert piv==[0,5,2]
 assert not lam[1]and lam[2]==Delta
 # Literal full original three-column coordinate identity, before any chart.
 for j in range(3):
  rhs={}
  for k in range(3):rhs=add(rhs,mul(lam[k],B[piv[k]][j]))
  assert mul(Delta,B[4][j])==rhs
 A=dec(c['target_coordinates']['a']);E=dec(c['target_coordinates']['e']);An=shift(A,(0,2),Q(500));En=shift(E,(0,1),Q(60))
 source_indices=[4,0,2];names=['S_raw12','S_raw7','S_raw10'];text=[]
 for q in[An,En]:text.append(str(len(q)));text.extend(f'{v} {m[0]} {m[1]}'for m,v in sorted(q.items()))
 text.append('3')
 for name,i in zip(names,source_indices):
  p=t['fixed_7x4_matrix'][i][3];text.append(f'{name} {len(p)}');text.extend(f'{v} {m[0]} {m[1]} {m[2]}'for m,v in p)
 nativein=out/'native-input.txt';nativein.write_text('\n'.join(text)+'\n');flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 exe=out/'substitute-native';subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20);subprocess.run([str(exe),str(nativein),str(out/'native.json')],check=True,timeout=40)
 nv=rd(out/'native.json');S={v['name']:dec(v['polynomial'])for v in nv['images']};assert len(S)==3
 images={v['name']:dec(v['polynomial'])for v in d['actual_images']};delta=images['Delta'];l0=images['lambda_4_0'];l2=images['lambda_4_2'];assert delta and l2==delta
 lhs=mul(delta,S['S_raw12']);r0=mul(l0,S['S_raw7']);r2=mul(l2,S['S_raw10']);defect=add(add(lhs,r0,Q(-1)),r2,Q(-1))
 stripping=[];remaining=defect.copy();dunit=0;content=Q()
 if remaining:
  dunit=min(m[1]for m in remaining);remaining=shift(remaining,(0,-dunit));num=0;den=1
  for v in remaining.values():num=math.gcd(num,abs(v.numerator));den=math.lcm(den,v.denominator)
  content=Q(num,den);remaining={m:v/content for m,v in remaining.items()}
  factors=[('zeta',{(1,0):Q(1)}),('zeta+20',{(1,0):Q(1),(0,0):Q(20)})]+[(name,dec(c['target_coordinates'][name]))for name in['L','K','Fnum','Jnum']]
  def divide(f,h):
   degree=max(m[0]for m in h);assert {m:v for m,v in h.items()if m[0]==degree}=={(degree,0):Q(1)}
   r=f.copy();q={};steps=0
   while r and max(m[0]for m in r)>=degree:
    high=max(m[0]for m in r);b={(m[0]-degree,m[1]):v for m,v in r.items()if m[0]==high};q=add(q,b);r=add(r,mul(b,h),Q(-1));steps+=1;assert steps<100
   assert add(mul(q,h),r)==f
   return q,r
  for name,factor in factors:
   exp=0
   while remaining:
    q,r=divide(remaining,factor)
    if r:break
    remaining=q;exp+=1;assert exp<100
   stripping.append({'known_unit':name,'exponent':exp})
  rebuilt=remaining.copy()
  for item,(_,factor)in zip(stripping,factors):
   for _ in range(item['exponent']):rebuilt=mul(rebuilt,factor)
  assert shift(rebuilt,(0,dunit),content)==defect
 stored=sum(len(q)for q in S.values())+sum(map(len,[delta,l0,l2,lhs,r0,r2,defect,remaining]));assert stored<=150000
 cert={'original_B_coordinate_group_raw12':{'Delta':enc(Delta),'lambda':[enc(q)for q in lam],'pivot_raw_rows':[7,13,10],'all3_full_literal_B_identities_verified':True},
       'actual_transformed_S':nv['images'],'transformed_Delta':enc(delta),'transformed_lambda0':enc(l0),'transformed_lambda2':enc(l2),
       'fixed_defect_formula':'Delta*S_raw12-lambda0*S_raw7-Delta*S_raw10',
       'three_product_polynomials':[enc(lhs),enc(r0),enc(r2)],'defect':enc(defect),'defect_zero':not defect,
       'known_unit_strip':{'d_exponent':dunit,'rational_content':str(content),'named_factor_exponents':stripping,'remaining_unknown_factor':enc(remaining)},
       'Delta_remains_nonunit_no_inverse':True,
       'scope':'Zero closes only this fixed raw12 source relation. Nonzero is an exact extra-column obstruction; any unknown residual factor and Delta=0 must retain their actual scope. No raw9/raw14 test or solver.',
       'no_rank_minor_GB_newpoint_prime_or_automatic_other_defect':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'all3_original_cleared_B_coordinate_identities_verified':True,'source_S_rows':[12,7,10],
         'S_image_termcounts':{name:len(p)for name,p in S.items()},'defect_terms':len(defect),'defect_identically_zero':not defect,
         'known_unit_factors':stripping,'unknown_residual_terms':len(remaining),'stored_polynomial_terms':stored,
         'Delta_not_inverted':True,'no_other_defect_evaluated':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'native_input_sha256':sha(nativein),'native_output_sha256':sha(out/'native.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
