#!/usr/bin/env python3
"""Only Delta and three lambda groups in the rational Theta chart, plus supports."""
import argparse,hashlib,json,resource,shlex,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_DELTA_LAMBDA_PROFILES_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def dec(p):return{tuple(m):Q(v)for m,v in p if Q(v)}
def profile(p):return{'terms':len(p),'zeta_maxdegree':max([m[0]for m in p]+[0]),'d_maxdegree':max([m[1]for m in p]+[0]),'d_min_exponent':min([m[1]for m in p]+[0])}
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';cp=ROOT/'reports/S61_NONSQUARE_THETA_RATIONAL_COORDINATES_20261104_v63';tp=ROOT/'reports/S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63'
 cpp=Path(__file__).with_suffix('.cpp');paths=[Path(__file__),cpp,gp/'report.json',gp/'certificate.json',cp/'report.json',cp/'certificate.json',tp/'report.json',tp/'certificate.json']
 rd=lambda p:json.loads(p.read_text());gr,g,cr,c,tr,t=[rd(p)for p in paths[2:]]
 for r,source in[(gr,paths[3]),(cr,paths[5]),(tr,paths[7])]:assert sha(source)==r['certificate_sha256']
 targets=c['target_coordinates'];A=dec(targets['a']);E=dec(targets['e']);An={(m[0],m[1]+2):v*500 for m,v in A.items()};En={(m[0],m[1]+1):v*60 for m,v in E.items()}
 group_rows=[1,4,6];items=[('Delta',g['nonzero_selected_minor_Delta'])]+[(f'lambda_{i}_{j}',g['cleared_upper_row_coordinates7xr'][i][j])for i in group_rows for j in range(3)]
 assert len(items)==10 and sum(len(p)for _,p in items)<10000
 text=[]
 for p in[An,En]:
  text.append(str(len(p)));text.extend(f'{v} {m[0]} {m[1]}'for m,v in sorted(p.items()))
 text.append(str(len(items)))
 for name,p in items:
  text.append(f'{name} {len(p)}');text.extend(f'{v} {m[0]} {m[1]} {m[2]}'for m,v in p)
 nativein=out/'native-input.txt';nativein.write_text('\n'.join(text)+'\n')
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 exe=out/'substitute-native';subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20)
 subprocess.run([str(exe),str(nativein),str(out/'native.json')],check=True,timeout=45)
 native=rd(out/'native.json');images={x['name']:dec(x['polynomial'])for x in native['images']};assert sum(len(p)for p in images.values())<=150000
 # Source input support estimates only. No coefficient substitution of S and
 # no Delta*S-lambda*S defect is evaluated in this packet.
 sourceS=[dec(row[3])for row in t['fixed_7x4_matrix']];sp=[set(A),set(E)];cache={}
 def mink(p,q):return{(x[0]+y[0],x[1]+y[1])for x in p for y in q}
 def suppower(j,k):
  if(j,k)not in cache:
   r={(0,0)}
   for _ in range(k):r=mink(r,sp[j])
   cache[j,k]=r
  return cache[j,k]
 def bound(p):
  ans=set()
  for (i,j,k)in p:
   r=mink(suppower(0,i),suppower(1,j));ans.update((x,y+k)for x,y in r)
  assert len(ans)<=150000
  return ans
 supports=[];scores=[]
 if not native['stopped_Delta_zero']:
  Sbounds=[bound(p)for p in sourceS];supports=[{'block_row':i,'raw_row':t['full_source_rows'][i],'source_terms':len(sourceS[i]),'substituted_support_upper':len(s),'zeta_degree_upper':max([x[0]for x in s]+[0]),'d_range_upper':[min([x[1]for x in s]+[0]),max([x[1]for x in s]+[0])]}for i,s in enumerate(Sbounds)]
  piv=g['selected_rows7'];assert piv==[0,5,2]
  for i in group_rows:
   pairupper=len(images['Delta'])*len(Sbounds[i])+sum(len(images[f'lambda_{i}_{j}'])*len(Sbounds[piv[j]])for j in range(3))
   sourcesets=[mink(set(images['Delta']),Sbounds[i])]+[mink(set(images[f'lambda_{i}_{j}']),Sbounds[piv[j]])for j in range(3)]
   outputbound=len(set().union(*sourcesets))
   scores.append({'block_row':i,'raw_row':t['full_source_rows'][i],'uncollected_product_pair_upper':pairupper,'stored_defect_support_upper':outputbound})
  scores.sort(key=lambda x:(x['uncollected_product_pair_upper'],x['raw_row']))
 cert={'source_Delta_and_lambda_names':[name for name,_ in items],'actual_images':native['images'],
       'Delta_nonzero_in_rational_Theta_chart':bool(images.get('Delta')),
       'source_input_binding':'Exact frozen Delta407 and raw9/raw12/raw14 row-coordinate groups, selected B rows raw7/raw13/raw10. No new minor computed.',
       'S_substitution_support_only_no_coefficients':supports,'candidate_costs_only':scores,
       'no_defect_rank_inverse_GCD_GB_point_prime_or_other_source_substitution':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'Delta_nonzero':bool(images.get('Delta')),'stopped_Delta_zero':native['stopped_Delta_zero'],
         'image_profiles':{name:profile(p)for name,p in images.items()},'candidate_support_only':scores,
         'stored_image_terms':native['stored_image_terms'],'peak_native_terms':native['peak_single_polynomial_terms'],
         'no_defect_polynomial_computed':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},
         'native_input_sha256':sha(nativein),'native_output_sha256':sha(out/'native.json'),'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
