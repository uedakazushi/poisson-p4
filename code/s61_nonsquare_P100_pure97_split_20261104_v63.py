#!/usr/bin/env python3
"""Explicit saved q5 GL2 split and seven saved G18 rows; no Gaussian."""
import argparse,gzip,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_P100_H1_right_inverse_20261104_v63 import read,sha
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_P100_PURE97_SPLIT_20261104_v63'
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports';rp=base/'S61_NONSQUARE_P100_Q5_PURE98_ROW_20261104_v63/pure98.json';gp=base/'S61_NONSQUARE_P100_PURE98_ROW_GCD_20261104_v63/native.json';gr=gp.with_name('report.json')
 b=base/'S61_NONSQUARE_P100_G18_PURE7_20261104_v63';br=read(b/'report.json');cpp=Path(__file__).with_suffix('.cpp');utility=ROOT/'code/s61_nonsquare_P100_q5_pure98_row_20261104_v63.cpp'
 paths=[Path(__file__),cpp,utility,rp,gp,gr,b/'report.json']+[b/f'pure98-row{i}.json'for i in range(7)]
 assert sha(gp)==read(gr)['native_output_sha256']and read(gr)['row_unimodular_on_whole_P100']and br['complete']
 r=read(rp)['values98'];dual=read(gp)['global_regular_row_dual'];assert [d['position']for d in dual]==[0,1]and all(d['denominator']==[[[0],'1']]for d in dual)
 rows=[]
 for i in range(7):
  p=b/f'pure98-row{i}.json';assert sha(p)==br['seven_row_sha256'][p.name];d=read(p);assert d.get('source319_index',d.get('actual_source319_index'))==174+i
  assert d['target122_columns']==read(rp)['target122_columns'];rows.append(d['values98'])
 # Source-only denominator/degree bounds, using only the five certified linears.
 factors=[(Q(1),Q(0)),(Q(1),Q(20)),(Q(1),Q(60)),(Q(203),Q(2700)),(Q(1),Q(12))]
 cache={}
 def deg(p):return max((m[0]for m,c in p),default=-1)
 def denpowers(p):
  key=tuple((m[0],c)for m,c in p)
  if key in cache:return cache[key]
  rem={m[0]:Q(c)for m,c in p};powers=[]
  for a,b0 in factors:
   exponent=0
   while rem and max(rem)>=1:
    q={};v=rem.copy()
    while v and max(v)>=1:
     n=max(v);c=v[n]/a;q[n-1]=q.get(n-1,Q())+c
     for j,k in[(n,c*a),(n-1,c*b0)]:
      v[j]=v.get(j,Q())-k
      if not v[j]:del v[j]
    if v:break
    rem=q;exponent+=1
   powers.append(exponent)
  assert rem and max(rem)==0;cache[key]=powers;return powers
 def common(row):
  powers=[denpowers(v['denominator'])for v in row];mx=[max(v[j]for v in powers)for j in range(5)];dd=sum(mx)
  numer=[-1 if not v['numerator']else deg(v['numerator'])+dd-sum(p)for v,p in zip(row,powers)]
  return dd,numer,mx
 dr,nr,er=common(r);dc=[deg(v['numerator'])for v in dual];profiles=[];outupper=0
 for row in rows:
  dg,ng,eg=common(row);an=max(ng[0]+dc[0]if ng[0]>=0 else-1,ng[1]+dc[1]if ng[1]>=0 else-1)
  ns=[max(ng[0]+nr[1],ng[1]+nr[0])]+[max(ng[j]+dr,an+nr[j])for j in range(2,98)];ds=dg+dr
  w=sum(max(0,n+1)+ds+1 for n in ns);outupper+=w;profiles.append({'common_denominator_powers5':eg,'source_common_denominator_degree':dg,'first_source_cofactor_numerator_degree_upper':an,'residual97_common_denominator_degree':ds,'residual97_numerator_degrees_upper':ns,'output_RF_terms_upper':w})
 allrf=r+dual+[v for row in rows for v in row];reserve=sum(len(v['numerator'])+len(v['denominator'])for v in allrf)
 inputdense=sum(max(0,deg(v['numerator'])+1)+deg(v['denominator'])+1 for v in allrf)
 livebound=reserve+inputdense+max(p['output_RF_terms_upper']for p in profiles)+15000
 pre={'q5_common_denominator_degree':dr,'q5_common_denominator_powers5':er,'q5_cleared_numerator_degrees98':nr,'dual_numerator_degrees':dc,'seven_profiles':profiles,'literal_input_terms':reserve,'input_RF_dense_terms_upper':inputdense,'peak_live_RF_terms_upper':livebound,'serialized_output_terms_upper':outupper+7*200,'separate_rows_streamed_gzip':True,'no_rank_Gaussian_or_unknown_inverses':True}
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n');print(json.dumps({'preflight_live_bound':livebound,'serialized_bound':outupper+1400,'max_one_output_row_terms_upper':max(p['output_RF_terms_upper']for p in profiles)},indent=2),flush=True);assert livebound<150000 and outupper+1400<150000
 text=[str(reserve)]
 def rf(v):
  for k in['numerator','denominator']:
   text.append(str(len(v[k])));text.extend(f'{m[0]} {c}'for m,c in v[k])
 for v in r+dual+[v for row in rows for v in row]:rf(v)
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');assert ip.stat().st_size<3*1024*1024
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();exe=out/'pure97-native'
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags,'-lz'],check=True,timeout=max(1,52-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:run=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-start)))
 a=read(out/'audit.json')if(out/'audit.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else{'all_passed':False,'returncode':run.returncode}
 complete=run.returncode==0 and a['all_passed']
 if complete:
  assert a['all686_G18_source_reconstruction_identities_verified']and a['peak_live_RF_terms_with_inputreserve']<=150000
  for i in range(7):
   with gzip.open(out/f'pure97-row{i}.json.gz','rt')as f:d=json.load(f)
   assert len(d['values97'])==97 and d['source319_index']==174+i
 cert={'target_GL2':[['c0','-r1'],['c1','r0']],'inverse_GL2':[['r0','r1'],['-c1','c0']],'column_j_clear':'j>=2: col_j -> col_j-r_j*col0','source_cofactor':'a=v0*c0+v1*c1','residual97':['-v0*r1+v1*r0','v_j-a*r_j for j>=2'],'full_source_reconstruction':['v0=a*r0-w1*c1','v1=a*r1+w1*c0','v_j=w_j+a*r_j for j>=2'],'target_order':'first residual column is original col1 afterGL2, then originalcols2..97','actual_pure_quotient_source':'G18*vj minus a_j times actual accepted q5H-cleared pure source; no pure tails discarded','same_coker_with_q5_unit_direction_removed':complete,'no_generic_or_uniform_rank_increment_claimed':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 report={'id':ID,'complete':complete,'preflight_live_bound':livebound,'native_audit':a,'samecoker_pure98_with_one_unit_split_and_seven_actual_pure97_rows':complete,'no_rank_Gaussian_or_automatic_nextjob':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
