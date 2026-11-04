#!/usr/bin/env python3
"""H4/g/k source preparation only, no coefficient inverse or d-Euclid."""
import argparse,hashlib,json,math,resource,signal,sys,time
from fractions import Fraction as Q
from pathlib import Path
from s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63 import Ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_GK_INPUT_20261104_v63';sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(a):int(v)for a,v in p}
def pa(a,b,c=1):
 z=dict(a)
 for i,v in b.items():
  z[i]=z.get(i,0)+c*v
  if not z[i]:del z[i]
 return z
def pm(a,b):
 z={}
 for i,c in a.items():
  for j,d in b.items():z[i+j]=z.get(i+j,0)+c*d
 return{k:v for k,v in z.items()if v}
def powers(p,n):
 z=[{0:1}]
 for i in range(n):z.append(pm(z[-1],p))
 return z
def exactdiv(p,g):
 z=dict(p);q={};lc=g[13]
 while z and max(z)>=13:
  k=max(z)-13;v=z[max(z)];assert v%lc==0;c=v//lc;q[k]=c;z=pa(z,{i+k:c*v for i,v in g.items()},-1)
 assert not z;return q
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63';eb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_E_RECONSTRUCTION_20261104_v63';bb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';pg=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_INPUT_20261104_v63';pk=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_MINOR_NORM_INPUT_20261104_v63'
 h=read(hb/'H4-monic-field-coefficients.json');g13={i:int(v)for i,v in h['g13']};ring=Ring(g13);tn=ip(h['e_numerator']);td=int(h['e_denominator']);t=ring.make(tn,td);a=ring.make({1:1});one=ring.c(1)
 def load(p):
  z={i:Q(v)for i,v in p};D=math.lcm(*(v.denominator for v in z.values()));return ring.make({i:int(v*D)for i,v in z.items()},D)
 H=list(map(load,h['H4_coefficients_d0_to4']));assert H[4]==one
 polys=read(bb/'certificate.json')['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'];sources=[polys[15],polys[20]]
 def dp(p):return {tuple(m):Q(v)for m,v in p}
 assert dp(sources[0])==dp(read(pg/'certificate.json')['g'])and dp(sources[1])==dp(read(pk/'certificate.json')['source_k']);assert list(map(len,sources))==[59,88]
 maxa=max(m[0]for p in sources for m,v in p);maxe=max(m[1]for p in sources for m,v in p);A=[one];T=[one]
 for i in range(maxa):A.append(ring.mul(A[-1],a))
 for i in range(maxe):T.append(ring.mul(T[-1],t))
 mon={(i,j):ring.mul(A[i],T[j])for i in range(maxa+1)for j in range(maxe+1)};tp=powers(tn,maxe);records=[];profiles=[]
 for name,p in zip(['g_Pi15','k_Pi20'],sources):
  adeg=max(m[0]for m,v in p);edeg=max(m[1]for m,v in p);ddeg=max(m[2]for m,v in p);D=math.lcm(*(Q(v).denominator for m,v in p));raw=[{}for _ in range(ddeg+1)];field=[ring.c(0)for _ in range(ddeg+1)]
  for m,v in p:
   q=Q(v);field[m[2]]=ring.add(field[m[2]],ring.scale(mon[m[0],m[1]],q.numerator,q.denominator));raw[m[2]]=pa(raw[m[2]],{i+m[0]:c for i,c in tp[m[1]].items()},int(q*D)*td**(edeg-m[1]))
  common=D*td**edeg;proof=[]
  for rr,x in zip(raw,field):
   diff=pa({i:c*x[1]for i,c in rr.items()},{i:c*common for i,c in x[0].items()},-1);quot=exactdiv(diff,g13);proof.append(dict(raw_numerator=[[i,str(v)]for i,v in sorted(rr.items(),reverse=True)],field=ring.enc(x),source_divisibility_quotient=[[i,str(v)]for i,v in sorted(quot.items(),reverse=True)]))
  before=list(field);quotient=[]
  while len(field)>4:
   c=field[-1];power=len(field)-5
   while len(quotient)<=power:quotient.append(ring.c(0))
   quotient[power]=c
   for j in range(5):field[power+j]=ring.sub(field[power+j],ring.mul(c,H[j]))
   assert field[-1]==ring.c(0);field.pop()
  while len(field)>1 and not field[-1][0]:field.pop()
  rebuilt=[ring.c(0)for _ in before]
  for i,c in enumerate(quotient):
   for j,v in enumerate(H):rebuilt[i+j]=ring.add(rebuilt[i+j],ring.mul(c,v))
  for i,c in enumerate(field):rebuilt[i]=ring.add(rebuilt[i],c)
  assert rebuilt==before
  terms=sum(len(x[0])for x in field);bitsn=max((abs(v).bit_length()for x in field for v in x[0].values()),default=0);bitsd=max(x[1].bit_length()for x in field)
  profiles.append(dict(name=name,source_terms=len(p),source_degrees=[adeg,edeg,ddeg],reduced_d_degree=len(field)-1,QQ_coefficient_terms=terms,numerator_bits=bitsn,denominator_bits=bitsd))
  records.append(dict(name=name,common_source_denominator=str(common),source_polynomial=p,coefficient_source_proofs=proof,raw_field_coefficients=[ring.enc(x)for x in before],H4_quotient=[ring.enc(x)for x in quotient],H4_remainder=[ring.enc(x)for x in field],whole_monic_H4_reduction_identity_checked=True))
 cert=dict(coefficient_algebra='QQ[a]/g13',g13=h['g13'],H4=[ring.enc(x)for x in H],gk=records,profiles=profiles,no_d_Euclid_or_coefficient_inverse_computed=True,source_projection='Existingallnormpencilcoeffs0 implyatleastonecommonH4/g/krootovereachgeometricbasefield; notallBconditions.',no_actual_bad_existence=True)
 payload=json.dumps(cert,indent=2)+'\n';assert len(payload.encode())<32*1024*1024;(out/'source-input.json').write_text(payload)
 plan=dict(input_d_degrees=[4]+[p['reduced_d_degree']for p in profiles],input_QQ_terms=53+sum(p['QQ_coefficient_terms']for p in profiles),initial_max_coefficient_bits=max([1986,1922]+[p[k]for p in profiles for k in['numerator_bits','denominator_bits']]),max_live_field_coefficients_under300=True,unit_inverse_query_bound9=9,coefficient_inverse_rule='AnyneededLCp(a): sourceprimitiveQQcontentonly thenxgcd(p,g13)=1 andliteralBezoutverified. Propergcdretain/splitbasecomponents; zeroLCskipped.',future_completion='Ifmonicgcd_d(H4,g,k)=d-s, wholeidealidentity+Ka13unitpivots makescommonrootcoverKa13 via d=s. ThenallB21ats mayrank3unitoremptiness, orallmaxminors0 givesactualcandidateonlyafterFULLU/sourcegate. No suchcomputationhere.',height_growth_not_assumed_small=True,proposed_budget='ONE60s/2GB/32MiB checkpointedddegree4Euclidwithunitcertificate-sizeguard, onlyafterparentallocation.',no_extra_minor_or_Norm4_planned=True)
 (out/'next-unit-Euclid-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
 paths=[Path(__file__),ROOT/'code/s61_nonsquare_raw_H_B_K13_named_saturation_20261104_v63.py',ROOT/f'research/notes/{ID}.md',hb/'report.json',hb/'H4-monic-field-coefficients.json',eb/'report.json',bb/'report.json',bb/'certificate.json',pg/'report.json',pg/'certificate.json',pk/'report.json',pk/'certificate.json']
 report=dict(id=ID,all_passed=True,exact_source_g15_k20_binding=True,all147_source_terms_and_H4_reductions_checked=True,profiles=profiles,no_d_Euclid_unitinverse_Gaussian_Norm4_newpoint_prime=True,field_integer_operations=ring.operations,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(json.dumps(plan,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
