#!/usr/bin/env python3
"""ONE component exact field Euclid with whole identities, denominator exceptions and U-root removal."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as FQ
from pathlib import Path
import sympy as s
from sympy.polys.agca.extensions import MonogenicFiniteExtension
ROOT=Path(__file__).resolve().parents[1];PREFIX='S61_NONSQUARE_RAW_H_B_COMPONENT_EXACT_QUOTIENT_';TAIL='_20261104_v63';CAP=150000

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main(component,out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(exist_ok=False,parents=True);ID=PREFIX+component+TAIL
 assert component=='Q'
 base=ROOT/'reports/S61_NONSQUARE_RAW_H_B_CURVE_COMPONENTS_INPUT_20261104_v63';ip=base/'certificate.json';assert sha(ip)==load(base/'report.json')['evidence_sha256'][ip.name];c=load(ip);rec=next(x for x in c['component_records']if x['factor_name']==component)
 a,e,d=s.symbols('a e d');K0=s.QQ.frac_field(e);expr=lambda p:sum(s.Rational(v)*a**m[0]*e**m[1]*d**m[2]for m,v in p);factor_expr=expr(rec['factor'])
 if component=='Q':
  K=MonogenicFiniteExtension(s.Poly(factor_expr,a,domain=K0));assert K.is_Field and K.rank==2
  def exactfieldquot(f,g):
   if not g:raise ZeroDivisionError('zeroextensionelement')
   h=f*g.inverse();assert h*g==f
   return h
  # Instance-onlyadapter keepsSymPy'snativeDomainconversiontypename unchanged.
  K.exquo=exactfieldquot;K.quo=exactfieldquot
 else:K=K0
 def toP(p):
  by={}
  for m,v in p:by[m[2]]=by.get(m[2],s.S.Zero)+s.Rational(v)*a**m[0]*e**m[1]
  return s.Poly.from_dict({(j,):K.from_sympy(v)for j,v in by.items()},(d,),domain=K)
 P=[toP(x['remainder'])for x in rec['polynomials']];assert all(p.degree()==4 for p in P)
 denoms={};peak=0;encoded_terms=0
 def rat(c):
  nonlocal encoded_terms
  nu,de=c.numer,c.denom
  nd=[[int(m[0]),str(v)]for m,v in sorted(nu.items(),reverse=True)if v];dd=[[int(m[0]),str(v)]for m,v in sorted(de.items(),reverse=True)if v];assert dd
  if de.degree()>0:
   lc=de.LC;key=tuple((int(m[0]),str(v/lc))for m,v in sorted(de.items(),reverse=True));denoms[key]=list(map(list,key))
  encoded_terms+=len(nd)+len(dd);assert encoded_terms<CAP
  return {'numerator_e':nd,'denominator_e':dd}
 def coeff(c):
  if component=='L':return [{'a_power':0,'fraction':rat(c)}]
  return [{'a_power':int(m[0]),'fraction':rat(v)}for m,v in sorted(c.rep.to_dict().items(),reverse=True)if v]
 def enc(p):
  nonlocal peak,encoded_terms
  before=encoded_terms;z=[{'d_power':int(m[0]),'coefficient':coeff(v)}for m,v in sorted(p.rep.to_dict().items(),reverse=True)if v];peak=max(peak,encoded_terms-before);return z
 def scalar_poly(c):return s.Poly.from_dict({(0,):c},(d,),domain=K)
 one=s.Poly(1,d,domain=K);zero=s.Poly(0,d,domain=K);trace=[]
 def gcdex(p,q,label):
  r0,r1=p,q;s0,s1=one,zero;t0,t1=zero,one;stage=[]
  while not r1.is_zero:
   lc=r1.rep.LC();inv=K.exquo(K.one,lc);assert K.mul(lc,inv)==K.one
   quotient,remainder=r0.div(r1);assert r0==quotient*r1+remainder
   stage.append({'divisor_degree':int(r1.degree()),'divisor_LC':coeff(lc),'divisor_LC_inverse':coeff(inv),'quotient':enc(quotient),'remainder':enc(remainder)})
   r0,r1=r1,remainder;s0,s1=s1,s0-quotient*s1;t0,t1=t1,t0-quotient*t1
   assert r0==s0*p+t0*q
   assert len(stage)<=6
  lc=r0.rep.LC();inv=K.exquo(K.one,lc);gg=r0.mul_ground(inv);ss=s0.mul_ground(inv);tt=t0.mul_ground(inv);assert ss*p+tt*q==gg
  trace.append({'label':label,'steps':stage,'gcd_monic':enc(gg),'source2_Bezout':[enc(ss),enc(tt)]});return ss,tt,gg
 s0,t0,g01=gcdex(P[0],P[1],'a5H_vs_g');u,v,G=gcdex(g01,P[2],'gcd01_vs_k');C=[u*s0,u*t0,v];assert sum((C[i]*P[i]for i in range(3)),zero)==G
 quot=[]
 for p in P:
  q,r=p.div(G);assert r.is_zero and q*G==p;quot.append(q)
 # Original localizer U reduced by the same monic component factor in a.
 sb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';sup=sb/'certificate.json';assert sha(sup)==load(sb/'report.json')['evidence_sha256'][sup.name];Up=s.Poly(expr(load(sup)['FULLU_product']),a,domain=s.QQ.poly_ring(e,d));fp=s.Poly(factor_expr,a,domain=s.QQ.poly_ring(e,d));uq,ur=Up.div(fp);assert uq*fp+ur==Up
 U=toP([[[int(m[0]),int(m[1]),int(m[2])],str(v)]for m,v in s.Poly(ur.as_expr(),a,e,d,domain=s.QQ).terms()]);assert not U.is_zero
 remaining=G;removed=[]
 while remaining.degree()>0:
  t=remaining.gcd(U).monic()
  if t.degree()==0:break
  q,rr=remaining.div(t);assert rr.is_zero and q*t==remaining;uquot,rr=U.div(t);assert rr.is_zero and uquot*t==U
  removed.append({'removed_generic_gcd_with_U':enc(t),'remaining_after':enc(q),'U_quotient':enc(uquot)});remaining=q
  assert len(removed)<=4
 boundary_only=remaining.degree()==0
 # Serialize all final identities before collecting the complete denominator list.
 final={'common_gcd':enc(G),'Bezout3':list(map(enc,C)),'all3_source_quotients_by_gcd':list(map(enc,quot)),'U_reduced_on_component':enc(U),'remaining_common_factor_after_U_support':enc(remaining)}
 proof={'component':component,'coefficient_domain':'QQ(e)'if component=='L'else'QQ(e)[a]/Q, irreduciblequadratic','factor_source':rec['factor'],'source_reduced_polynomials':[x['remainder']for x in rec['polynomials']],'exact_Euclid_traces':trace,'final':final,'U_support_removal_steps':removed,'generic_common_roots_all_U_boundary':boundary_only,'finite_e_exception_polynomials':list(denoms.values()),'all3_Bezout_and_divisibility_and_U_removal_identities_checked':True,'full_U_generic_nonzero':True,'LC_nonunit_exceptions_not_discarded':True,'generic_vs_global':'Fieldinverses are valid onlyoffthe displayedfinitee exceptions; no globalemptyclaim. OnL a ispolynomialine; onQ a isdegree2 overe, and actuala!=0 keepsHmonicd4, so displayedexceptions supportonlyfinitelymany actualBbad points onthiscomponent.','no_other_minor_or_globalB3_inference':True,'exact_Q_field_division_adapter':'instanceK.exquo/quo=f*g.inverse(); everyproductchecked; botholdAPIpreflightfailurespreserved','no_GB_point_prime_or_nested_resultant_jobs':True}
 (out/'certificate.json').write_text(json.dumps(proof,indent=2)+'\n');paths=[Path(__file__),ROOT/'reports/S61_NONSQUARE_RAW_H_B_COMPONENT_COMMON_ROOT_Q_20261104_v63/failure.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_COMPONENT_FIELD_REPAIR_Q_20261104_v63/failure.json',ROOT/'research/notes'/'S61_NONSQUARE_RAW_H_B_COMPONENT_COMMON_ROOT_20261104_v63.md',base/'report.json',ip,sb/'report.json',sup];report={'id':ID,'component':component,'source_terms':[len(x['remainder'])for x in rec['polynomials']],'common_gcd_d_degree':int(G.degree()),'remaining_after_U_support_d_degree':int(remaining.degree()),'generic_roots_all_on_U_boundary':boundary_only,'number_exception_polynomials':len(denoms),'total_exception_degree_upper':sum(max(x[0]for x in p)for p in denoms.values()),'serialized_field_polynomial_terms':encoded_terms,'peak_field_polynomial_terms':peak,'all_whole_field_identities_checked':True,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--component',choices=['L','Q'],required=True);p.add_argument('--output',type=Path);args=p.parse_args();main(args.component,(args.output or(ROOT/'reports'/(PREFIX+args.component+TAIL))).resolve())
