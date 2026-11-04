#!/usr/bin/env python3
"""Exact monic4 norm INPUT only; no determinant/GB/factorization/evaluation."""
import hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_PAIR_NORM_INPUT_20261104_v63';CAP=150000
peak=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def poly(p):return{tuple(m):Q(c)for m,c in p if Q(c)}
def ser(p):return[[list(m),str(c)]for m,c in sorted(p.items(),reverse=True)]
def guard(p):
 global peak
 peak=max(peak,len(p));assert len(p)<=CAP;return p
def add(*ps):
 r={}
 for p in ps:
  for m,c in p.items():
   v=r.get(m,Q())+c
   if v:r[m]=v
   else:r.pop(m,None)
  guard(r)
 return r
def mul(p,q):
 r={}
 for m,c in p.items():
  for n,d in q.items():
   k=tuple(x+y for x,y in zip(m,n));v=r.get(k,Q())+c*d
   if v:r[k]=v
   else:r.pop(k,None)
  guard(r)
 return r
def scale(p,c):return{m:x*c for m,x in p.items()if x*c}
def shift(p,v):return{tuple(m[i]+v[i]for i in range(3)):c for m,c in p.items()}
def divide_monic(p,H):
 r=dict(p);q={};steps=[]
 while r and max(m[2]for m in r)>=4:
  n=max(m[2]for m in r);t={(m[0],m[1],m[2]-4):c for m,c in r.items()if m[2]==n};q=add(q,t);r=add(r,scale(mul(t,H),-1));steps.append({'old_d_degree':n,'pivot_terms':len(t),'remaining_terms':len(r)})
 assert add(mul(q,H),r)==p
 return q,r,steps
def stats(p):
 return {'terms':len(p),'a_exponent_min':min([m[0]for m in p],default=0),'a_exponent_max':max([m[0]for m in p],default=0),'e_degree':max([m[1]for m in p],default=0),'d_degree':max([m[2]for m in p],default=-1),'coefficient_numerator_bits':max([abs(c.numerator).bit_length()for c in p.values()],default=0),'coefficient_denominator_bits':max([c.denominator.bit_length()for c in p.values()],default=0)}
def main():
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 mb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_MONIC_COMPLETION_20261104_v63';mp=mb/'certificate.json';mr=load(mb/'report.json');assert sha(mp)==mr['evidence_sha256'][mp.name];mc=load(mp);H=poly(mc['monic_H_laurent']);assert {(m[0],m[1],0):c for m,c in H.items()if m[2]==4}=={(0,0,0):Q(1)}
 pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_COPRIME_PAIR_PREFLIGHT_20261104_v63';pp=pb/'certificate.json';pr=load(pb/'report.json');assert sha(pp)==pr['evidence_sha256'][pp.name];pc=load(pp);f,g=map(poly,pc['selected_pair_source']['polynomials']);assert pc['selected_pair']==[10,15]
 assert add(scale(shift(H,(5,0,0)),Q(-11)),scale(add(shift(f,(1,0,0)),scale(mul(poly(mc['q']),g),-1)),Q(-9)))=={}
 cols=[];quot=[];trace=[];M=[[{}for j in range(4)]for i in range(4)]
 for j in range(4):
  q,r,t=divide_monic(shift(g,(0,0,j)),H);cols.append(r);quot.append(q);trace.append(t)
  for i in range(4):M[i][j]={(m[0],m[1],0):c for m,c in r.items()if m[2]==i}
 assert sum(len(x)for row in M for x in row)<=CAP
 for j in range(4):assert add(*(shift(M[i][j],(0,0,i))for i in range(4)),mul(quot[j],H))==shift(g,(0,0,j))
 allp=[p for row in M for p in row];nonzero=[p for p in allp if p];mina=min(m[0]for p in nonzero for m in p);maxa=max(m[0]for p in nonzero for m in p);maxe=max(m[1]for p in nonzero for m in p);aclear=max(0,-mina);den=math.lcm(*(c.denominator for p in nonzero for c in p.values()));maxint=max(abs(int(c*den))for p in nonzero for c in p.values());support=max(map(len,nonzero));height=math.factorial(4)*(support*maxint)**4
 # Matrix uniformly cleared by D*a^aclear; determinant has (D*a^aclear)^4 scaling.
 bounds={'matrix_uniform_a_clear_power':aclear,'matrix_uniform_QQ_denominator':str(den),'cleared_matrix_max_a_degree':maxa+aclear,'cleared_matrix_max_e_degree':maxe,'cleared_matrix_max_integer_coefficient_bits':maxint.bit_length(),'det_cleared_a_degree_upper':4*(maxa+aclear),'det_e_degree_upper':4*maxe,'det_cleared_support_box_upper':(4*(maxa+aclear)+1)*(4*maxe+1),'det_cleared_integer_abs_coefficient_bits_upper':height.bit_length(),'det_laurent_a_exponent_lower':4*mina,'det_laurent_a_exponent_upper':4*maxa,'det_permutation_count':24}
 bounds['det_serialization_conservative_bytes_upper']=bounds['det_cleared_support_box_upper']*(100+(height.bit_length()*30103+99999)//100000+len(str(den**4)))
 assert bounds['det_serialization_conservative_bytes_upper']<32*1024*1024
 cert={'variable_order':['a','e','d'],'quotient_coefficient_ring':'QQ[a,a^-1,e]','basis':['1','d','d^2','d^3'],'H':ser(H),'g':ser(g),'matrix_rows':[[ser(p)for p in row]for row in M],'column_quotients_by_H':[ser(p)for p in quot],'monic_division_steps_by_column':trace,'all4_full_polynomial_identities':'g*d^j=H*column_quotient[j]+sum_i M[i,j]*d^i, j=0..3','all16_entry_stats':[[stats(p)for p in row]for row in M],'matrix_exact_total_terms':sum(map(len,allp)),'matrix_exact_nonzero_entries':len(nonzero),'stage_peak_single_polynomial_support':peak,'determinant_cost_bounds_NOT_value':bounds,'determinant_resultant_GB_point_factorization_jobs':0,'pair_ideal_equality':'(f,g)=(H,g) in QQ[a,a^-1,e,d], because H=(-9/(11*a^4))*f+(9*q/(11*a^5))*g and f=q*g/a-(11*a^4/9)*H.','norm_nonzero_justification':'QQ[a,e,d]gcd(f,g)=1 implies QQ(a,e)[d]gcd(H,g)=1, so multiplication by g is invertible overQQ(a,e)[d]/H. Its4norm is nonzero, but no determinant computed.','next_exact_certificate_recipe':'If R=detM and v=adj(M)*e0 then g*sum(v_j*d^j)-R is divisible by the monicH. Substitute the savedH combination to get a literalR membership in(f,g); all division is monic or knowna-unit.','no_unknown_localizer_or_leading_coefficient_inverted':True}
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),mb/'report.json',mp,pb/'report.json',pp]
 report={'id':ID,'input_H_g_terms':[len(H),len(g)],'exact_matrix_size':[4,4],'exact_matrix_terms':cert['matrix_exact_total_terms'],'entry_term_counts':[[len(p)for p in row]for row in M],'cost_bounds':bounds,'all4_full_polynomial_division_identities_checked':True,'resultant_or_determinant_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
