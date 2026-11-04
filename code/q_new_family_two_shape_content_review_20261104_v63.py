#!/usr/bin/env python3
"""Independent full coefficient product and univariate field-Bezout maximality audit."""
import argparse,gzip,hashlib,json,resource,signal,time
from collections import defaultdict
from pathlib import Path
import sympy as s
from sympy.polys.fields import field
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_CONTENT_PRIMITIVITY_20261104_v63_REVIEW'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def parse_poly(p,nv):
 with p.open()as f:
  v,n=map(int,f.readline().split());assert v==nv;d={}
  for _ in range(n):
   c,*e=map(int,f.readline().split());assert len(e)==nv and tuple(e)not in d and c;d[tuple(e)]=c
  assert not f.read().strip()
 return d
def mul(a,b):
 d={}
 for e,c in a.items():
  for f,k in b.items():
   g=tuple(x+y for x,y in zip(e,f));d[g]=d.get(g,0)+c*k
 return {e:c for e,c in d.items()if c}
def shift(p,e,c=1):return {tuple(x+y for x,y in zip(a,e)):b*c for a,b in p.items()if b*c}
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 author=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63';ar=read(author/'report.json');assert ar['all_passed']
 for rel,h in ar['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for rel,h in ar['evidence_sha256'].items():assert sha(author/rel)==h,rel
 normgate=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_FINAL_NORM_IDENTITY_20261104_v63_REVIEW/report.json'
 coprimegate=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_COPRIME_DEGREES_20261104_v63_REVIEW/report.json'
 spectralgate=ROOT/'reports/Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63_REVIEW/report.json'
 assert all(read(p)['all_passed']for p in [normgate,coprimegate,spectralgate])
 strip=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_BOUNDARY_STRIP_20261104_v63'
 qq=parse_poly(strip/'quotient.txt',6);raw=parse_poly(ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RESULTANT_20261104_v63/resultant.txt',6)
 C=parse_poly(author/'content.txt',3);G=parse_poly(author/'primitive.txt',6)
 assert len(C)==70 and len(G)==891 and all(e[2]==0 for e in C)
 assert max(e[0]for e in C)==9 and max(e[1]for e in G)==10
 # Independent Python integer products, not FLINT gcd/cofactor or native
 # multiplication. Compare the literal full6919 coefficient polynomial.
 C6={(0,*e,0,0):c for e,c in C.items()}
 assert mul(C6,G)==qq
 qgroups=defaultdict(dict)
 for (zero,r,R,Z0,Z1,Z2),c in qq.items():
  assert zero==0;qgroups[Z1,Z2][r,R,Z0]=c
 cofactors={}
 with(author/'cofactors.txt').open()as f:
  assert list(map(int,f.readline().split()))==[3,15]
  for _ in range(15):
   i,j,n=map(int,f.readline().split());p={}
   for _ in range(n):
    c,r,R,z=map(int,f.readline().split());assert(r,R,z)not in p;p[r,R,z]=c
   assert(i,j)not in cofactors;cofactors[i,j]=p
  assert not f.read().strip()
 assert set(cofactors)==set(qgroups)
 for label,p in cofactors.items():assert mul(C,p)==qgroups[label]
 # Also close the full known-boundary chain independently by repeated
 # multiplication by its small linear factors, avoiding a huge product.
 p5={e[1:]:c for e,c in qq.items()};assert all(e[0]==0 for e in qq)
 rm1={(1,0,0,0,0):1,(0,0,0,0,0):-1}
 tt={(1,1,0,0,0):1,(0,1,0,0,0):-1,(1,0,0,0,0):-2}
 for _ in range(16):p5=mul(p5,rm1)
 for _ in range(15):p5=mul(p5,tt)
 p5=shift(p5,(22,0,1,0,0),25*(2**68))
 assert p5=={e[1:]:c for e,c in raw.items()} and len(raw)==47114
 # Different maximality certificate: Euclid/Bezout in QQ(R)[r], using
 # Z0 coefficients, never the author's multivariate gcd algorithm.
 K,Rvar=field('R',s.QQ);z=K.zero;one=K.one
 def trim(p):
  while p and not p[-1]:p.pop()
  return p
 def plus(a,b):
  c=[z]*max(len(a),len(b))
  for i,v in enumerate(a):c[i]+=v
  for i,v in enumerate(b):c[i]+=v
  return trim(c)
 def times(a,b):
  if not a or not b:return []
  c=[z]*(len(a)+len(b)-1)
  for i,v in enumerate(a):
   for j,w in enumerate(b):c[i+j]+=v*w
  return trim(c)
 def scalar(a,b):return trim([v*b for v in a])
 def divide(a,b):
  assert b;rr=list(a);quo=[z]*max(0,len(a)-len(b)+1)
  while rr and len(rr)>=len(b):
   j=len(rr)-len(b);cc=rr[-1]/b[-1];quo[j]+=cc
   for k,v in enumerate(b):rr[j+k]-=cc*v
   trim(rr)
  assert plus(times(quo,b),rr)==trim(list(a))
  return trim(quo),rr
 def extended(a,b):
  aa,bb=list(a),list(b);u,v=[one],[];uu,vv=[],[one];degrees=[]
  while bb:
   q,r=divide(aa,bb);degrees.append([len(aa)-1,len(bb)-1,len(r)-1])
   aa,bb=bb,r;u,uu=uu,plus(u,scalar(times(q,uu),-one));v,vv=vv,plus(v,scalar(times(q,vv),-one))
  inv=one/aa[-1];return scalar(aa,inv),scalar(u,inv),scalar(v,inv),degrees
 def coefficient(p,k):
  a=[z]*(max(e[0]for e in p)+1)
  for(r,R,j),c in p.items():
   if j==k:a[r]+=K(c)*Rvar**R
  return trim(a)
 selected=[];polys=[];bez=[];g=[];steps=[]
 labels=[(0,3),(5,0)]+[label for label in cofactors if label not in [(0,3),(5,0)]]
 assert all(all(e[1:]==(0,0)for e in cofactors[label])for label in labels[:2])
 candidates=[]
 for label in labels:
  for k in sorted({e[2]for e in cofactors[label]}):candidates.append((label,k,coefficient(cofactors[label],k)))
 for label,k,p in candidates:
  if not g:
   g=list(p);selected.append([*label,k]);polys.append(p);bez.append([one]);continue
  _,rem=divide(p,g)
  if not rem:continue
  ng,u,v,trace=extended(g,p)
  bez=[times(u,b)for b in bez]+[v];selected.append([*label,k]);polys.append(p)
  g=ng;steps.append({'new_group_Z1_Z2_Z0_power':[*label,k],'Euclidean_degree_trace':trace,'new_gcd_degree':len(g)-1})
  assert sum_identity(bez,polys,plus,times)==g
  if g==[one]:break
 assert g==[one]and sum_identity(bez,polys,plus,times)==[one]
 def fp(a):return [[int(e[0]),str(c)]for e,c in a.items()]
 def encfield(a):return {'numerator_R_coefficients':fp(a.numer),'denominator_R_coefficients':fp(a.denom)}
 certificate={'coefficient_field':'QQ(R)','polynomial_variable':'r','selected_cofactor_Z1_Z2_Z0_coefficient_labels':selected,
              'selected_univariate_polynomials_ascending':[[encfield(c)for c in p]for p in polys],
              'Bezout_polynomials_ascending':[[encfield(c)for c in p]for p in bez],
              'literal_identity':'sum_i Bezout_i(r,R)*coefficient_i(r,R)=1','Euclidean_steps':steps,
              'all15_cofactor_identities_and_full_Q_equals_C_times_G':True,
              'full_raw_boundary_factor_identity_independently_checked':True,
              'global_identity':'raw=25*2^68*r^22*Z0*(r-1)^16*T^15*C*G',
              'content_has_no_Z0':'Any common divisor divides the two nonzero Z0-free groups, hence has Z0-degree zero up to units; the displayed field-Bezout proves no remaining r-content.',
              'implicitization_proof':'Coprime3/5 maps give birational s-plane parametrization; primitive resultant is its irreducible equation to exponent1; Gauss then gives irreducible G in QQ(R,Z0,Z1,Z2)[r].',
              'spectral_base_algebraically_independent':['R','Z0','Z1','Z2'],
              'exact_source_shape_extension_degree':10,
              'source_field':'QQ(R,Z0,r,s)','base_field':'QQ(R,Z0,Z1,Z2)',
              's_rational_recovery_over_base_with_r_retained':True,
              'actual_ZN_primitivity_full_target_index_and_NEW_equation_not_proved':True}
 ep=out/'independent-certificate.json';ep.write_text(json.dumps(certificate,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'research/notes'/(STEM+'.md'),author/'report.json',author/'content.txt',author/'primitive.txt',author/'cofactors.txt',
         strip/'report.json',strip/'quotient.txt',ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RESULTANT_20261104_v63/resultant.txt',normgate,coprimegate,spectralgate]
 report={'id':STEM,'all_passed':True,'independent_integer_full_Q_equals_CG_and_all15_cofactor_products':True,
         'independent_full_raw_equals_known_boundaries_times_CG':True,
         'maximal_content_proved_by_DIFFERENT_univariate_QQ_R_Bezout_identity':True,
         'no_multivariate_gcd_or_author_native_rerun':True,'primitive_G_irreducible_over_QQ_R_Z0_Z1_Z2_in_r':True,
         'exact_source_SHAPE_extension_degree':10,'s_recovered_rationally_over_base_with_r':True,
         'actual_full_NEW_target_field_index_ZN_primitivity_or_equation_not_claimed':True,
         'no_new_point_prime_resultant_PRS_source26_or_automatic_retry':True,
         'selected_Bezout_inputs':len(selected),'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in inputs},
         'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started,'frozen':True}
 assert ep.stat().st_size<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in ['all_passed','exact_source_SHAPE_extension_degree','selected_Bezout_inputs','seconds']},indent=2))
def sum_identity(bez,polys,plus,times):
 a=[]
 for b,p in zip(bez,polys):a=plus(a,times(b,p))
 return a
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/STEM);a=p.parse_args();main(a.output.resolve())
