#!/usr/bin/env python3
"""Known-source full covectors, with a prior guard before signed construction."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,add,image
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_COVECTORS_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 kp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_KNOWN_INPUT_GATE_20261105_v63/pruned-known619-lower4-input-and-prior-resource-gate.json';d=read(kp)
 prefix='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_FIXED_ROWS_20261105_v63';r=read(prefix+'/report.json');assert r['all_passed'] and r['next_clearing_ready']
 p=prefix+'/actual-fixed-specific-dual-rows-and-F25-with-next-clearing-prior.json';f=read(p);assert sha(ROOT/p)==r['evidence_sha256'][Path(p).name]
 FF=list(map(sp,f['actual_F25_columns27']));S=[list(map(sp,M)) for M in f['actual_specific_dual_rows2_by2']]
 J=list(map(sp,d['J25_full918_columns']));Pi=list(map(sp,d['PiY0_source102_columns']));raw=[]
 for sr,se in S:
  a={9*i+7:x for i,x in sr.items()};b={9*i+7:x for i,x in se.items()}
  for i,x in sr.items():b[9*i+8]=x
  raw.extend([a,b])
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def plus(a,b):
  if not a[1]:return mark(b)
  if not b[1]:return mark(a)
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def accum(v,k,b):
  if b[1]:v[k]=plus(v.get(k,(1,0)),b)
 bounds=[]
 for L in raw:
  w=[]
  for col in J:
   b=(1,0)
   for i,c in col.items():b=plus(b,mul(ex(c),ex(L.get(i,Q(0)))))
   w.append(b)
  v={}
  for i,col in enumerate(FF):
   for j,c in col.items():accum(v,i,mul(w[j],ex(c)))
  a={i:ex(c) for i,c in L.items()}
  for row,col in enumerate(Pi):
   for target,c in col.items():
    for var in range(9):accum(a,9*row+var,mul(ex(c),v.get(9*target+var,(1,0))))
  bounds.append(a)
 assert maximum==f['next_clearing_prior_bits']==75
 # Include the fresh direct annihilation checks before signed construction.
 for a in bounds:
  for col in J:
   b=(1,0)
   for i,c in col.items():b=plus(b,mul(a.get(i,(1,0)),ex(c)))
 assert maximum<=384
 prior=maximum;owned=count(d)+count([f,raw,J,Pi,FF,bounds])+2*4*918+12000;assert owned<150000
 # No Gaussian elimination or source evaluation: these are known operator rows.
 alpha=[];actualbits=0
 def check(x):
  nonlocal actualbits
  actualbits=max(actualbits,abs(x.numerator).bit_length(),x.denominator.bit_length());assert actualbits<=prior;return x
 for L,bound in zip(raw,bounds,strict=True):
  w=[check(sum((c*L.get(i,Q(0)) for i,c in col.items()),Q(0))) for col in J]
  projected={i:check(sum((w[j]*c for j,c in col.items()),Q(0))) for i,col in enumerate(FF)}
  a=dict(L)
  for row,col in enumerate(Pi):
   for target,c in col.items():
    for var in range(9):
     k=9*row+var;x=check(a.get(k,Q(0))-check(c*projected.get(9*target+var,Q(0))))
     if x:a[k]=x
     else:a.pop(k,None)
  for i,x in a.items():
   D,U=bound[i];assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
  for col in J:assert check(sum((c*a.get(i,Q(0)) for i,c in col.items()),Q(0)))==0
  alpha.append(a)
 assert not alpha[0] and not alpha[1]
 # The52 chosen degree2 output coordinate covectors of G2_1*Y1.
 G2=list(map(sp,d['G2linear_source102_columns']));A=[{} for _ in range(52)]
 for j,col in enumerate(G2):
  for k,c in col.items():
   target,var=divmod(k,9)
   if var==7:
    A[2*target]=add(A[2*target],{9*j+7:c})
    A[2*target+1]=add(A[2*target+1],{9*j+8:c})
   elif var==8:A[2*target+1]=add(A[2*target+1],{9*j+7:c})
 assert all(all(k%9 in (7,8) for k in a) for a in A)
 zsupport=[[i for i in a if i%9<7] for a in alpha]
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-four-full918-source-covectors-and52-plane-adjoint-input.json'
 ev={'full918_covectors4':[[[i,str(x)] for i,x in sorted(a.items())] for a in alpha],
 'plane52_output_covectors_on918_slots':[[[i,str(x)] for i,x in sorted(a.items())] for a in A],
 'all25_G3linear_annihilations_checked':True,'all714_z_slots_zero':not any(zsupport),'z_support4':zsupport,
 'positive_build_and_annihilation_prior_bits':prior,'actual_known_covector_bits':actualbits,'owned_numeric_bound':owned,
 'no_adjoint_solve_or_X2_619_values':True,'next_scope':'Only if all z slots vanish may a52-output adjoint be sought. Every coefficient is known-source algebra; no actual619 relation image has been evaluated.',
 'source_formula':'alpha=L-(LJ25)F25PiY0; L=chiHstar1B99 onplaneY1, withall25actualG3firstclearingsretained.'}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'z_slots_zero':not any(zsupport),'known_covector_prior_bits':prior,'actual_bits':actualbits,'no_adjoint_solve_or_619_values':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'zslots_zero':not any(zsupport),'support':[len(a) for a in alpha],'zsupport':[len(a) for a in zsupport],'prior':prior,'actual':actualbits,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
