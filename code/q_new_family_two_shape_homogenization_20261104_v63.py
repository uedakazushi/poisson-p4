#!/usr/bin/env python3
"""ONE bounded polynomial homogenization/tracked-division job; NO generic gcd/GB."""
import argparse,gzip,hashlib,json,resource,signal,time
from fractions import Fraction
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60second wall cap')));signal.alarm(60)
 r,s,q,B,R=S.symbols('r s q B R');variables=(r,s,R);P=lambda x:S.Poly(x,*variables,domain=S.QQ)
 base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63';review=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_BIRATIONAL_COMPRESSION_20261104_v63_REVIEW/report.json'
 recovery=ROOT/'reports/Q_NEW_FAMILY_MOBIUS_SHAPE_RECOVERY_20261104_v63/report.json';assert read(review)['all_passed']and read(recovery)['all_passed']
 paths=[];poly={};front=[];homogeneous={};pairs={};constraints={};error=None
 def save(stage):
  front.append(dict(stage=stage,seconds=time.monotonic()-started));(out/'frontier.json').write_text(json.dumps(front,indent=2)+'\n')
 def dump_poly(p):return [[str(c.p),str(c.q),*m]for m,c in p.terms()]
 def record(p):return dict(terms=len(p.terms()),multidegree=list(p.degree_list()),coefficients=dump_poly(p))
 try:
  for name in['F4','p0','p1','p2','N12']:
   path=base/(name+'.json.gz');paths.append(path);d=read(path)
   poly[name]=S.Poly.from_dict({tuple(e):S.Rational(c,int(d['denominator_integer']))for c,*e in d['terms']},(r,s,q,B),domain=S.QQ)
  T=(r-1)*R-2*r;L=8*R*r*(r-1)-6*r*r-11*r+1;Nq=(s-1)*L-(r-1);J=(10*r-1)*s-10*r;alpha=S.Rational(2,405)
  Hn=64*r*T*T*s*s-Nq*Nq;Fn=(r-1)*Nq*Nq+16*r*T*T*(s*s-4*r*s+4*r)
  factors={k:P(v)for k,v in dict(r=r,rm1=r-1,T=T,Nq=Nq,J=J,Hn=Hn,Fn=Fn).items()}
  numerator=P(Nq);denominator=P(40*r*T);powersN=[P(1)];powersD=[P(1)]
  for i in range(1,9):powersN.append(powersN[-1]*numerator);powersD.append(powersD[-1]*denominator)
  def split(name,power):
   return S.Poly.from_dict({e[:3]:c for e,c in poly[name].terms()if e[3]==power},(r,s,q),domain=S.QQ)
  def hom(name,pp):
   degree=pp.degree(q);ans=P(0)
   for (i,j,k),c in pp.terms():ans+=c*P(r**i*s**j)*powersN[k]*powersD[degree-k]
   assert len(ans.terms())<=20000
   homogeneous[name]=dict(q_degree=int(degree),**record(ans));save('homogenized_'+name);return ans,int(degree)
  def normalize(name,num,constant,den):
   den={k:v for k,v in den.items()if v};num=P(num);constant=S.Rational(constant);canceled={}
   for k in list(den):
    if den[k]<0:num*=factors[k]**(-den[k]);del den[k]
   for k in list(den):
    while den[k]>0:
     quotient,remainder=num.div(factors[k])
     if not remainder.is_zero:break
     assert quotient*factors[k]==num;num=quotient;den[k]-=1;canceled[k]=canceled.get(k,0)+1
   den={k:v for k,v in den.items()if v};content,num=num.primitive();constant*=content
   ans=dict(numerator=num,constant=constant,den=den,canceled=canceled)
   pairs[name]=dict(numerator=record(num),constant=str(constant),denominator_factors=den,verified_exact_factor_divisions=canceled)
   save('tracked_reduction_'+name);return ans
  # h(qback)=Hn/(64rT²), f(qback)=alpha Fn/(16rT²),
  # L0(qback)=alpha²184320r(r−1)²Nq²J/T³.
  for i in[1,2]:
   A,m=hom('A'+str(i),split('p'+str(i),0));normalize('A'+str(i),A,S.Rational(16384,1)/(alpha**2*40**m),{'r':m-3,'T':m-6,'Hn':1,'Fn':2})
   C,m=hom('L'+str(i),split('p'+str(i),1));normalize('C'+str(i),C,1/(alpha**2*184320*40**m),{'r':m+1,'T':m-3,'rm1':2,'Nq':2,'J':1})
  M,m=hom('M0',split('N12',1));normalize('E1',M,S.Rational(1024,1)/(alpha**3*184320*40**m),{'r':m-1,'T':m-7,'Hn':1,'Fn':1,'rm1':2,'Nq':2,'J':1})
  M,m=hom('M1',split('N12',2));normalize('E2',P(Fn)*M,1/(16*alpha**3*184320**2*40**m),{'r':m+3,'T':m-4,'rm1':4,'Nq':4,'J':2})
  def pair_obj(name):
   z=pairs[name];return dict(num=S.Poly.from_dict({tuple(e):S.Rational(a,b)for a,b,*e in z['numerator']['coefficients']},variables,domain=S.QQ),constant=S.Rational(z['constant']),den=z['denominator_factors'])
  def cleared(name,first,second,powers):
   a,b=pair_obj(first),pair_obj(second);lcm={k:max(a['den'].get(k,0),b['den'].get(k,0))for k in factors};lcm={k:v for k,v in lcm.items()if v};den=P(1)
   for k,e in lcm.items():den*=factors[k]**e
   nums=[]
   for x in[a,b]:
    num=x['constant']*x['num']
    for k,e in lcm.items():num*=factors[k]**(e-x['den'].get(k,0))
    nums.append(num)
   assert len(den.terms())<=20000 and all(len(z.terms())<=20000 for z in nums)
   degrees=[z.degree_list()for z in nums+[den]];bideg=[max(d[j]for d in degrees)for j in[0,1]]
   constraints[name]=dict(equation=f'coefficient1*Z0^{powers[0]}+coefficient2*Z0^{powers[1]}-target*denominator=0',
    coefficient1=record(nums[0]),coefficient2=record(nums[1]),denominator=record(den),
    denominator_factors=lcm,shape_bidegree=bideg,R_degree=max(d[2]for d in degrees),
    terms_with_independent_target_parameter=sum(len(z.terms())for z in nums+[den]),
    Z0_powers=powers,global_gcd_not_computed=True)
   save('localized_constraint_'+name+'_saved')
  cleared('Z1','A1','C1',[0,1]);cleared('Z2','A2','C2',[0,1]);cleared('ZN','E1','E2',[1,2])
  cert=dict(variables=list(map(str,variables)),T=str(T),Nq=str(Nq),J=str(J),Hn=str(Hn),Fn=str(Fn),
   factors={k:record(v)for k,v in factors.items()},homogenized_input_polynomials=homogeneous,rational_coefficient_pairs=pairs,localized_constraints=constraints,
   scope='xi*r*s*Z0*(r-1)*T*Nq*J*Hn*Fn!=0. Only exact localized equivalence; no extension on zeros, properness, targetindex1, or actual degree assertion',
   no_multivariate_gcd_cancel_factorization_or_GB_used=True)
  (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
  assert (out/'certificate.json').stat().st_size<=32*1024*1024
 except Exception as exc:error=repr(exc)
 signal.alarm(0)
 # Even if the cap fires, preserve every completed polynomial/pair/constraint.
 if error:(out/'partial-certificate.json').write_text(json.dumps(dict(homogenized=homogeneous,pairs=pairs,constraints=constraints),indent=2)+'\n')
 report=dict(id=STEM,all_passed=error is None,error=error,one_polynomial_homogenization_job=True,
  exact_division_only_by_named_nonzero_chart_factors=True,no_generic_gcd_cancel_factorization_GB_elimination_source26_newpoints=True,
  completed_homogenized_inputs=list(homogeneous),completed_rational_pairs=list(pairs),completed_constraints=list(constraints),
  constraint_summary={k:{z:v[z]for z in['shape_bidegree','R_degree','terms_with_independent_target_parameter','denominator_factors']}for k,v in constraints.items()},
  limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,certificate_bytes=32*1024*1024),
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),review,recovery,*paths]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','error','constraint_summary','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
