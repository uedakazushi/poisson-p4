#!/usr/bin/env python3
"""One fixed raw11/raw7 B-row relation and Theta Schur-column defect."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_THETA_ONE_DEFECT_20261104_v63'
PAIRS=0;CAP=150000
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
 global PAIRS
 PAIRS+=len(p)*len(q);assert PAIRS<=CAP,('fixed_total_product_support_budget',PAIRS,CAP)
 z={}
 for m,v in p.items():
  for n,w in q.items():
   k=tuple(x+y for x,y in zip(m,n));z[k]=z.get(k,Q())+v*w
 return{k:v for k,v in z.items()if v}
def scale(p,c):return{m:c*v for m,v in p.items()if c*v}
def shift(p,s,c=Q(1)):return{tuple(x+y for x,y in zip(m,s)):c*v for m,v in p.items()if c*v}
def dec(p):return{tuple(m):Q(v)for m,v in p if Q(v)}
def enc(p):return[[list(m),str(v)]for m,v in sorted(p.items())]
def stats(p):return{'terms':len(p),'a_e_d_maxdegree':[max([m[j]for m in p]+[0])for j in range(3)],'d_min_exponent':min([m[2]for m in p]+[0])}
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(57);out.mkdir(parents=True,exist_ok=False)
 gp=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';tp=ROOT/'reports/S61_NONSQUARE_RAW_H_THETA_UNIT_PROFILE_20261104_v63';pp=ROOT/'reports/S61_NONSQUARE_RAW_H_P_BRANCH_PROFILES_20261104_v63';raw=ROOT/'reports/S61_NONSQUARE_RAW_H_LINEAR_UNITS_20261104_v63'
 paths=[Path(__file__),gp/'report.json',gp/'certificate.json',tp/'report.json',tp/'certificate.json',pp/'report.json',pp/'certificate.json',raw/'report.json',raw/'checkpoint.json']
 rd=lambda p:json.loads(p.read_text());gr,g,tr,t,pr,p,rr,rc=[rd(x)for x in paths[1:]]
 for r,cpath in[(gr,paths[2]),(tr,paths[4]),(pr,paths[6])]:assert sha(cpath)==r['certificate_sha256']
 assert sha(paths[8])==rr['evidence_sha256']['checkpoint.json']
 B=[[dec(x)for x in row]for row in g['actual_block7x3']];M=[[dec(x)for x in row]for row in t['fixed_7x4_matrix']]
 R=[[dec(x)for x in row]for row in rc['matrix15x6']]
 assert t['full_source_rows']==[7,9,10,11,12,13,14]
 assert all(B[i][j]==R[t['full_source_rows'][i]][[2,3,4][j]]for i in range(7)for j in range(3))
 # Fixed small-coordinate clearings. No ratio is assumed to be a known unit.
 U=B[3][0];V=B[0][0];u=M[3][0];v=M[0][0]
 estimate=3*(len(U)*len(V)+len(V)*len(U))+len(U)*len(dec(g['nonzero_selected_minor_Delta']))+len(V)*len(dec(g['cleared_upper_row_coordinates7xr'][3][0]))+len(v)*len(M[3][3])+len(u)*len(M[0][3])
 assert estimate<100000
 (out/'preflight.json').write_text(json.dumps({'fixed_rows':[11,7],'clearing_U_V_terms':[len(U),len(V)],'reduced_u_v_terms':[len(u),len(v)],'initial_product_pairs_bound':estimate,'total_product_pairs_guard':CAP},indent=2)+'\n')
 assert all(add(mul(V,B[3][j]),mul(U,B[0][j]),Q(-1))=={}for j in range(3))
 Delta=dec(g['nonzero_selected_minor_Delta']);lam=dec(g['cleared_upper_row_coordinates7xr'][3][0])
 assert all(not g['cleared_upper_row_coordinates7xr'][3][j]for j in[1,2])
 assert add(mul(Delta,U),mul(lam,V),Q(-1))=={}
 Theta=dec(p['Theta']);monic=shift(Theta,(0,0,-2),Q(1,2500))
 def divtheta(f):
  r=f.copy();q={};steps=0
  while r and max(m[0]for m in r)>=2:
   h=max(m[0]for m in r);b={(m[0]-2,m[1],m[2]):c for m,c in r.items()if m[0]==h};q=add(q,b);r=add(r,mul(b,monic),Q(-1));steps+=1;assert steps<20
  q=shift(q,(0,0,-2),Q(1,2500));assert add(mul(q,Theta),r)==f
  return q,r
 middleproof=[]
 for j in range(3):
  source=add(mul(v,M[3][j]),mul(u,M[0][j]),Q(-1));q,r=divtheta(source);assert not r;middleproof.append(enc(q))
 defect_before=add(mul(v,M[3][3]),mul(u,M[0][3]),Q(-1));q,defect=divtheta(defect_before)
 # The nonunit Delta is only the old B-rank-three chart. No inverse is added.
 cert={'fixed_raw_source_rows':[11,7],'block_rows':[3,0],'U':enc(U),'V':enc(V),
       'literal_all3_source_relation':'V*B[3,:]-U*B[0,:]=0','literal_relation_verified':True,
       'Delta':enc(Delta),'old_lambda_first':enc(lam),'Delta_U_minus_lambda_V_zero':True,
       'Theta_reduced_U':enc(u),'Theta_reduced_V':enc(v),'all3_middle_cleared_relations_Theta_quotients':middleproof,
       'cleared_extra_column_formula':'v*S[3]-u*S[0] moduloTheta (S is the accepted B^4-cleared last Schur column)',
       'defect_before_final_Theta_reduction':enc(defect_before),'quotient_by_Theta':enc(q),'defect_remainder':enc(defect),
       'defect_profile':stats(defect),'defect_identically_zero':not defect,
       'scope':'U,V,Delta remain nonunits. A nonzero defect gives a generic additional direction on the recorded B-rank-three chart after also checking nonvanishing there; a zero defect closes only this fixed pair route. No global rank/minor/inverse or other row test.',
       'total_actual_sparse_product_pairs':PAIRS,'no_other_defect_solver_minors_GB_points_or_primes':True}
 (out/'certificate.json').write_text(json.dumps(cert,separators=(',',':'))+'\n')
 result={'id':ID,'all3_literal_B_row_relation_verified':True,'old_Delta_lambda_compatibility_verified':True,
         'Theta_middle_relations_verified':True,'single_defect_profile':stats(defect),'single_defect_identically_zero':not defect,
         'total_actual_sparse_product_pairs':PAIRS,'support_guard':CAP,
         'input_sha256':{str(x.relative_to(ROOT)):sha(x)for x in paths},'certificate_sha256':sha(out/'certificate.json'),
         'seconds':time.monotonic()-started,'frozen':True,'independent_review_pending':True}
 assert sum(x.stat().st_size for x in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
