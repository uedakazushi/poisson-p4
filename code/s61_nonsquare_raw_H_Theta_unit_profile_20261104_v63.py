#!/usr/bin/env python3
"""One fixed Theta row8 unit and seven-by-four Schur profiles; no solver."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
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
def scale(p,c):return{m:c*v for m,v in p.items()if c*v}
def shift(p,m,c=Q(1)):return{tuple(x+y for x,y in zip(k,m)):c*v for k,v in p.items()if c*v}
def power(p,k):
 z={(0,0,0):Q(1)}
 for _ in range(k):z=mul(z,p)
 return z
def enc(p):return[[list(m),str(v)]for m,v in sorted(p.items())]
def dec(p):return{tuple(m):Q(v)for m,v in p if Q(v)}
def stat(p):return{'terms':len(p),'a_e_d_maxdegrees':[max([m[j]for m in p]+[0])for j in range(3)],'min_d_exponent':min([m[2]for m in p]+[0]),'max_coefficient_bits':max([max(abs(v.numerator).bit_length(),v.denominator.bit_length())for v in p.values()]+[0])}
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 prior=ROOT/'reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63';raw=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63'
 paths=[Path(__file__),prior/'report.json',prior/'certificate.json',raw/'report.json',raw/'checkpoint.json',ROOT/'reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63_REVIEW/report.json']
 rd=lambda p:json.loads(p.read_text());pr,pc,rr,rc=[rd(p)for p in paths[1:5]]
 assert sha(paths[2])==pr['certificate_sha256']and sha(paths[4])==rr['evidence_sha256']['checkpoint.json'];assert rd(paths[5])['all_passed']
 M=[[dec(p)for p in row]for row in rc['matrix15x6']]
 branch=next(b for b in pc['branches']if b['branch']=='Theta');records={(r['raw_checkpoint_row'],r['column']):r for r in branch['full14x2_records']}
 Theta=dec(pc['Theta']);g1=dec(pc['g1']);g2=dec(pc['g2']);monic=dec(branch['monic_polynomial'])
 a={(1,0,0):Q(1)};e={(0,1,0):Q(1)};d={(0,0,1):Q(1)};one={(0,0,0):Q(1)}
 C=add(add(add(scale(a,100),scale(power(e,2),18)),scale(e,-9)),scale(one,-2))
 B=add(mul(add(add(add(scale(a,50),scale(power(e,2),18)),scale(e,-9)),scale(one,-2)),d),add(scale(e,60),scale(one,-40)))
 assert add(Theta,mul(C,mul(g1,g2)))==power(B,2)
 assert all(not M[8][j]for j in[2,3,4])
 f=scale(mul(g1,power(g2,2)),Q(1,37500));rec=records[8,0]
 qf=dec(rec['quotient_by_original_branch_polynomial']);assert f==dec(rec['remainder'])and add(mul(qf,Theta),f)==M[8][0]
 invnum=scale(mul(power(C,2),g1),Q(37500));B4=power(B,4)
 # f*invnum=(B^2-Theta)^2; the full original row entry has a saved quotient.
 unitquot=add(Theta,scale(power(B,2),Q(-2)))
 assert add(mul(f,invnum),B4,Q(-1))==mul(Theta,unitquot)
 fullunitquot=add(unitquot,mul(qf,invnum));assert add(mul(M[8][0],invnum),B4,Q(-1))==mul(Theta,fullunitquot)
 def divide(f):
  rem=f.copy();quo={};steps=0
  while rem and max(m[0]for m in rem)>=2:
   highest=max(m[0]for m in rem);block={(m[0]-2,m[1],m[2]):v for m,v in rem.items()if m[0]==highest}
   quo=add(quo,block);rem=add(rem,mul(block,monic),Q(-1));steps+=1
   assert steps<25 and len(rem)<50000 and len(quo)<50000
  actualq=shift(quo,(0,0,-2),Q(1,2500));assert add(mul(actualq,Theta),rem)==f
  return actualq,rem,steps
 # Only the seven B-support rows survive after deleting the row8 unit.
 support=[7,9,10,11,12,13,14];deleted=[r for r in rc['remaining_rows']if r not in support+[8]]
 for r in deleted:
  assert all(not M[r][j]for j in[2,3,4])
  assert all(not dec(records[r,j]['remainder'])for j in[0,5])
 lastfactor=mul(invnum,M[8][5]);profiles=[];newmatrix=[]
 for row in support:
  vals=[]
  for col in[2,3,4,5]:
   source=M[row][col]if col!=5 else add(mul(B4,M[row][5]),mul(M[row][0],lastfactor),Q(-1))
   q,r,steps=divide(source);vals.append(r)
   profiles.append({'raw_checkpoint_row':row,'new_column':col,'source_kind':'unchanged_middle'if col!=5 else'B4*original_last-original_first*inverse_numerator*original_pivot_last',
                    'quotient_by_Theta':enc(q),'remainder':enc(r),'remainder_profile':stat(r),'division_steps':steps})
  newmatrix.append(vals)
 cert={'branch':'FULLU/(Theta)','source_pivot_row':8,'source_pivot_column':0,'all_middle_three_pivot_entries_literally_zero':True,
       'pivot_remainder':enc(f),'pivot_remainder_formula':'g1*g2^2/37500','pivot_original_quotient_by_Theta':enc(qf),
       'inverse_numerator':enc(invnum),'inverse_denominator':enc(B4),'inverse_formula':'37500*C^2*g1/B^4',
       'full_pivot_times_inverse_minus_one_cleared_quotient_by_Theta':enc(fullunitquot),
       'full_source_rows':support,'new_target_columns':[2,3,4,'B4*last-Schur'],
       'zero_rows_after_specialization':deleted,'fixed_7x4_matrix':[[enc(p)for p in row]for row in newmatrix],
       'all28_full_original_quotient_identities':profiles,
       'block_rank_identity':'At every residue field on FULLU/(Theta), rank(heavy14x5)=1+rank(this7x4). Its first3 columns are unchanged B7x3 moduloTheta.',
       'known_unit_denominators':['nonzeroQQ','d','B'],'scheme_validity':'All identities hold moduloTheta without irreducibility/reducedness assumptions.',
       'no_7x4_rank_minor_kernel_solver_or_full_cubic_claim':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'row8_factor_and_full_source_quotient_verified':True,'known_unit_inverse_and_source_quotient_verified':True,
         'all28_full_source_Theta_remainder_identities_verified':True,'new_shape':[7,4],
         'column_nonzero_counts':[sum(bool(row[j])for row in newmatrix)for j in range(4)],
         'column_remainder_terms':[sum(len(row[j])for row in newmatrix)for j in range(4)],
         'max_single_remainder_terms':max(len(p)for row in newmatrix for p in row),
         'all_last_Schur_entries_zero':all(not row[3]for row in newmatrix),
         'no_rank_minors_GB_new_point_prime_or_automatic_operation':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
