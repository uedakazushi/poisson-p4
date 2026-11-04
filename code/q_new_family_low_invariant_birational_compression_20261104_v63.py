#!/usr/bin/env python3
"""Exact sparse rational compression of seven accepted low invariants."""
import argparse, gzip, hashlib, json, math, resource, signal, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def clean(p):return {k:v for k,v in p.items()if v}
def mul(p,q):
 d={}
 for e,c in p.items():
  for f,b in q.items():
   k=tuple(x+y for x,y in zip(e,f));d[k]=d.get(k,0)+c*b
 return clean(d)
def pw(p,n):
 q={(0,0,0,0):Q(1)}
 for _ in range(n):q=mul(q,p)
 return q
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s cap')));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63'
 rev=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63_REVIEW/report.json'
 assert json.loads(rev.read_text())['all_passed']
 inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),rev];maps={};results=[]
 for name in ['F4','p0','p1','p2','p3','p4','N12']:
  path=base/(name+'.json.gz');inputs.append(path);d=json.loads(gzip.decompress(path.read_bytes()))
  cp=d['source_c_power'];hp=d['denominator_H_power'];bp=d['denominator_B_power'];assert cp%2==0
  poly={}
  for row in d['terms']:
   num,den,i,j,k,l=row;i,j,k,l=map(int,(i,j,k,l))
   exr=i-j+cp//2+hp;exs=2*j-i;exb=l+cp-bp
   assert min(exr,exs,exb)>=0,(name,row)
   c=Q(int(num),int(den))*Q(76800**i,3686400**j*3072**k)*Q(25**hp,int(d['denominator_integer']))
   for t in range(k+1):
    ex=(exr,exs,t,exb);poly[ex]=poly.get(ex,0)+c*math.comb(k,t)
  poly=clean(poly);maps[name]=poly
  lcm=math.lcm(*(c.denominator for c in poly.values()));ints={e:int(c*lcm)for e,c in poly.items()}
  g=math.gcd(lcm,math.gcd(*(abs(c)for c in ints.values())));lcm//=g;ints={e:c//g for e,c in ints.items()}
  data=dict(name=name,variables=['r=u^2/1600','s=48*u*X','q=3072*W-1','B'],
   xi_power=cp,denominator_integer=str(lcm),denominator_h_power=hp,h='s^2-25*r*q^2',
   terms=[[str(c),*e]for e,c in sorted(ints.items(),reverse=True)])
  dest=out/(name+'.json.gz');dest.write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0))
  results.append(dict(name=name,terms=len(ints),xi_power=cp,denominator_h_power=hp,
   denominator_integer=str(lcm),multidegree=[max(e[i]for e in ints)for i in range(4)]))
 F={(2,0,2,0):Q(100),(1,0,2,0):Q(-100),(1,1,0,0):Q(-4),(0,2,0,0):Q(1),(1,0,0,0):Q(4)}
 assert maps['F4']=={e:c*Q(2,405)for e,c in F.items()}
 # Q0=p0/F4^2. Its exact B-linear recovery requires the dense D(p0) chart.
 l={(0,1,0,0):Q(-1),(0,0,1,0):Q(5),(0,0,0,0):Q(1)}
 rminus={(1,0,0,0):Q(1),(0,0,0,0):Q(-1)}
 expected=mul({(4,0,2,1):Q(131072000000,2278125)},mul(rminus,l))
 assert maps['p0']==expected
 recipe=dict(
  old_coordinates=['a=u/X','Y=X^2','W','B','c/X'],
  new_coordinates=['xi=40*c/(B*u)','r=u^2/1600','s=48*u*X','q=3072*W-1','B'],
  inverse=dict(a='76800*r/s',Y='s^2/(3686400*r)',W='(q+1)/3072',c_over_X='1920*xi*B*r/s'),
  forward_in_old=dict(xi='40*(c/X)/(B*a)',r='Y*a^2/1600',s='48*Y*a',q='3072*W-1'),
  h='s^2-25*r*q^2',H='h/(25*r)',F='100*r*(r-1)*q^2-4*r*s+s^2+4*r',
  F4='(2/405)*xi^2*F',Q0='p0/F4^2=1/400-4*p4/F4^2',
  D='2359296000*r^4*(r-1)*q^2*(-s+5*q+1)',
  B_recovery='B=Q0*h*F^2/D',
  recovery_scope='Dense rational source chart r*s*B*h*F*D!=0; no global localization assertion',
  known_tau_fixed_source_field='QQ(xi,r,s,q,B)',
  actual_NEW_target_field_equality_claimed=False)
 (out/'coordinate-recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
 report=dict(id=ID,all_passed=True,results=results,exact_rational_coordinate_substitution=True,
  F4_five_term_numerator_and_p0_factorization_verified=True,B_rationally_recovered_from_normalized_p0=True,
  new_beta_point_relation_or_higher_degree_job=False,NEW_equation_constructed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);a=p.parse_args();main(a.output.resolve())
