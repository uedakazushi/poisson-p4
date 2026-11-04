#!/usr/bin/env python3
"""Fixed p factor/comaximal identity and28+28 monic remainder profiles only."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def add(p,q,c=Q(1)):
 z=p.copy()
 for m,v in q.items():
  z[m]=z.get(m,Q())+c*v
  if not z[m]:del z[m]
 return z
def mul(p,q):
 z={}
 for m,c in p.items():
  for n,d in q.items():
   k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,Q())+c*d
 return{k:v for k,v in z.items()if v}
def monshift(p,m,c=Q(1)):return{tuple(x+y for x,y in zip(k,m)):c*v for k,v in p.items()if c*v}
def scale(p,c):return{k:c*v for k,v in p.items()if c*v}
def encode(p):return[[list(m),str(c)]for m,c in sorted(p.items())]
def profile(p):return{'terms':len(p),'a_e_d_maxdegrees':[max([m[j]for m in p]+[0])for j in range(3)],'a_e_d_minexponents':[min([m[j]for m in p]+[0])for j in range(3)],'max_coefficient_bits':max([max(abs(c.numerator).bit_length(),c.denominator.bit_length())for c in p.values()]+[0])}
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False);read=lambda p:json.loads(p.read_text());dec=lambda p:{tuple(m):Q(c)for m,c in p if Q(c)}
 raw=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63';gp=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';bind=ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63'
 paths=[Path(__file__),raw/'report.json',raw/'checkpoint.json',gp/'report.json',gp/'certificate.json',bind/'report.json',bind/'certificate.json',ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63_REVIEW/report.json']
 rr,rc,gr,gc,br,bc=[read(p)for p in paths[1:7]]
 assert sha(paths[2])==rr['evidence_sha256']['checkpoint.json']and sha(paths[4])==gr['certificate_sha256']and sha(paths[6])==br['certificate_sha256'];assert read(paths[7])['all_passed']
 matrix=[[dec(p)for p in row]for row in rc['matrix15x6']];rows=rc['remaining_rows'];cols=[0,5]
 a={(1,0,0):Q(1)};e={(0,1,0):Q(1)};d={(0,0,1):Q(1)};one={(0,0,0):Q(1)}
 g1=add(add(scale(mul(e,d),3),scale(d,-2)),scale(one,20));g2=add(add(scale(mul(e,d),6),d),scale(one,40))
 C=add(add(add(scale(a,100),scale(mul(e,e),18)),scale(e,-9)),scale(one,-2))
 B=add(add(mul(add(add(add(scale(a,50),scale(mul(e,e),18)),scale(e,-9)),scale(one,-2)),d),scale(e,60)),scale(one,-40))
 Theta=add(mul(B,B),mul(C,mul(g1,g2)),Q(-1));p=matrix[2][0];assert p==dec(gc['p'])==matrix[4][0]
 assert scale(mul(g1,Theta),Q(1,750000))==p
 # Cleared exact comaximal certificate B^2=Theta+C*g1*g2.
 assert add(Theta,mul(C,mul(g1,g2)))==mul(B,B)
 assert{m:v for m,v in Theta.items()if m[0]==2}=={(2,0,2):Q(2500)}
 # Exact monic branch models invert only d and nonzero rational constants.
 G=monshift(g1,(0,0,-1),Q(1,3));T=monshift(Theta,(0,0,-2),Q(1,2500))
 assert{m:v for m,v in G.items()if m[1]==1}=={(0,1,0):Q(1)}
 assert{m:v for m,v in T.items()if m[0]==2}=={(2,0,0):Q(1)}
 def divide(f,H,var,degree):
  rem=f.copy();quo={};steps=0
  while rem and max(m[var]for m in rem)>=degree:
   highest=max(m[var]for m in rem)
   block={tuple(v-(degree if j==var else 0)for j,v in enumerate(m)):c for m,c in rem.items()if m[var]==highest}
   quo=add(quo,block);rem=add(rem,mul(block,H),Q(-1));steps+=1
   assert steps<20 and len(quo)<20000 and len(rem)<20000
  assert add(mul(quo,H),rem)==f
  return quo,rem,steps
 branches=[]
 for name,H,var,degree,scaleMon,scaleQ in[('g1',G,1,1,(0,0,-1),Q(1,3)),('Theta',T,0,2,(0,0,-2),Q(1,2500))]:
  records=[];remainders={}
  for row in rows:
   remainders[row]={}
   for col in cols:
    original=matrix[row][col];quot,rem,steps=divide(original,H,var,degree)
    actualquot=monshift(quot,scaleMon,scaleQ)
    assert add(mul(actualquot,g1 if name=='g1'else Theta),rem)==original
    assert all(m[0]>=0 and m[1]>=0 for m in rem)
    remainders[row][col]=rem
    records.append({'raw_checkpoint_row':row,'column':col,'quotient_by_original_branch_polynomial':encode(actualquot),'remainder':encode(rem),'remainder_profile':profile(rem),'monic_reduction_steps':steps})
  assert not remainders[2][0]and not remainders[4][0]
  assert remainders[2][5]==remainders[4][5]
  sparse=[row for row in rows if all(not matrix[row][j]for j in[2,3,4])]
  branches.append({'branch':name,'monic_polynomial':encode(H),'polynomial_leading_variable':['a','e','d'][var],
   'monic_degree':degree,'coefficient_ring':'QQ[a,d,d^-1]'if name=='g1'else'QQ[e,d,d^-1]',
   'full14x2_records':records,'nonzero_first_column_rows':[row for row in rows if remainders[row][0]],
   'nonzero_last_column_rows':[row for row in rows if remainders[row][5]],
   'total_remainder_terms':sum(len(remainders[row][col])for row in rows for col in cols),
   'max_single_remainder_terms':max(len(remainders[row][col])for row in rows for col in cols),
   'sparse_other3_zero_rows':sparse,'known_pair_first_column_zero_and_last_equal':True,
   'no_rank_or_alternate_pair_minor_computed':True})
 cert={'selected_p_source':'reviewed rawcompactcheckpoint row2,column0',
  'source_p':encode(p),'g1':encode(g1),'g2':encode(g2),'Theta':encode(Theta),
  'literal_factor_identity':'750000*p=g1*Theta','full_coefficient_factor_check':True,
  'comaximal_identity_cleared':'B^2=Theta+C*g1*g2',
  'comaximal_identity_on_FULLU':'1=Theta/B^2+(C*g2/B^2)*g1',
  'p0_exact_scheme_decomposition':'OnFULLU theideals(g1,Theta)arecomaximal; R/(p)=R/(g1)xR/(Theta). No intersectionpiece discarded without this certificate.',
  'known_unit_factors_inverted':['nonzeroQQ','d','B (comaximalidentityonly)'],
  'branches':branches,
  'all56_quotient_remainder_identities_verified':True,
  'all_five_heavy_columns_scope':'Only first/last are profiled. Otherthree retain accepted7x3B dependence; alternatepair/globalheavy rank remains open.',
  'no_14x5_allminors_GB_prime_point_or_pureprojection_claim':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'actual_p_factor_and_FULLU_comaximal_identity_verified':True,
  'p0_branches':['g1','Theta'],'all56_full_quotient_remainder_identities_passed':True,
  'branch_profiles':[{k:b[k]for k in['branch','nonzero_first_column_rows','nonzero_last_column_rows','total_remainder_terms','max_single_remainder_terms','sparse_other3_zero_rows']}for b in branches],
  'no_rank_minor_solver_newpoint_prime_GB_or_automatic_extra_operation':True,
  'no_Bbad_cubicbad_or_globalheavy_rank_claim':True,
  'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),
  'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
