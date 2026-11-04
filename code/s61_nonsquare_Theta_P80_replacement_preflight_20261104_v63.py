#!/usr/bin/env python3
"""Knownchartunit row/column preconditioning ofP80matrix; cofactorcostONLY, no solve/rank."""
import hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_THETA_P80_REPLACEMENT_PREFLIGHT_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def poly(p):return {int(m):Q(v)for m,v in p if Q(v)}
def ser(p):return [[m,str(v)]for m,v in sorted(p.items())if v]
def multiply(p,q):
 if not p or not q:return{}
 dp=math.lcm(*(c.denominator for c in p.values()));dq=math.lcm(*(c.denominator for c in q.values()));r={}
 for m,c in p.items():
  for n,v in q.items():r[m+n]=r.get(m+n,0)+int(c*dp)*int(v*dq)
 return {m:Q(v,dp*dq)for m,v in r.items()if v}
def add(p,q):
 r=dict(p)
 for m,c in q.items():
  v=r.get(m,Q())+c
  if v:r[m]=v
  else:r.pop(m,None)
 return r
def divide(p,q):
 r=dict(p);s={};n=max(q);lc=q[n]
 while r and max(r)>=n:
  m=max(r);c=r[m]/lc;s[m-n]=s.get(m-n,Q())+c;r=add(r,{m-n+j:-c*v for j,v in q.items()})
 return s,r
def powpoly(p,n):
 r={0:Q(1)}
 for i in range(n):r=multiply(r,p)
 return r
def linval(p,L):
 if not p:return 10**9
 r=p;n=0
 while r and max(r)>=1:
  q,rem=divide(r,L)
  if rem:break
  r=q;n+=1
 return n
def costs(M,rows,cols,extras):
 integer=[[{m:int(v)for m,v in p.items()}for p in row]for row in M];assert all(v.denominator==1 for row in M for p in row for v in p.values());rowdeg=[max([max(p)for p in integer[i]if p],default=0)for i in rows];bounds=[]
 for c in range(17):
  for extra in extras:
   used=cols[:c]+[extra]+cols[c+1:];bounds.append(sum(max([max(integer[i][j])for j in used if integer[i][j]],default=0)for i in rows))
 rownorm=[max(sum(abs(v)for v in p.values())for p in integer[i])for i in rows];height=math.factorial(17)*math.prod(rownorm);degree=max(bounds);bits=height.bit_length();return {'row_max_a_or_zeta_degree':rowdeg,'max_all24_row_degree_sum':sum(rowdeg),'max119_replacement_minor_degree_upper':degree,'all119_replacement_coeff_slot_upper':119*(degree+1),'single17vector_coeff_slot_upper':17*(degree+1),'replacement_integer_coefficient_bits_upper':bits,'single17vector_raw_byte_upper_without_additional_Bezout':17*(degree+1)*(100+(bits*30103+99999)//100000),'ALL119rawcofactor_dump_NOT_authorized':True,'future_gcd_Bezout_total_size_NOT_inferred_from_cofactor_bound':True,'matrix_terms':sum(len(p)for row in M for p in row),'matrix_max_integer_coefficient_bits':max(abs(v.numerator).bit_length()for row in M for p in row for v in p.values())}
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(exist_ok=False,parents=True)
 base=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_GENERIC_20261104_v63';rep=load(base/'report.json');ip=base/'substituted17x24-before-rank.json';np=base/'native-rank.json';assert sha(ip)==rep['evidence_sha256'][ip.name]and sha(np)==rep['evidence_sha256'][np.name];matrix=load(ip)['integer_matrix17x24'];M=[[poly(p)for p in row]for row in matrix];raw=[[dict(p)for p in row]for row in M];nr=load(np);rows=nr['selected_local_rows'];cols=nr['selected_local_columns'];extras=[j for j in range(24)if j not in cols];assert cols==list(range(17))and extras==list(range(17,24))and set(rows)==set(range(17));detsource=poly(nr['nonzero_selected_minor']);assert max(detsource)==447
 ub=ROOT/'reports/S61_NONSQUARE_THETA_P80_CUBIC_MINOR_UNITS_20261104_v63';up=ub/'exact-known-unit-strip.json';assert sha(up)==load(ub/'report.json')['evidence_sha256'][up.name];u=load(up);linear=[{1:Q(x['linear_a_b'][0]),0:Q(x['linear_a_b'][1])}for x in u['all6_known_unit_power_divisions']];linear=[{m:v for m,v in p.items()if v}for p in linear];rf=[{0:Q(1)}for i in range(17)];cf=[{0:Q(1)}for j in range(24)];rexp=[[0]*6 for i in range(17)];cexp=[[0]*6 for j in range(24)]
 for i in range(17):
  for k,L in enumerate(linear):
   v=min(linval(p,L)for p in M[i]);assert v<10**9
   if v:
    factor=powpoly(L,v);rf[i]=multiply(rf[i],factor);rexp[i][k]+=v
    for j in range(24):M[i][j],rem=divide(M[i][j],factor);assert not rem
 for j in range(24):
  for k,L in enumerate(linear):
   v=min(linval(M[i][j],L)for i in range(17));assert v<10**9
   if v:
    factor=powpoly(L,v);cf[j]=multiply(cf[j],factor);cexp[j][k]+=v
    for i in range(17):M[i][j],rem=divide(M[i][j],factor);assert not rem
 rowqq=[]
 for i in range(17):
  den=math.lcm(*(v.denominator for p in M[i]for v in p.values()));content=math.gcd(*(int(v*den)for p in M[i]for v in p.values()));assert content>0;scale=Q(content,den);rowqq.append(scale);rf[i]={m:v*scale for m,v in rf[i].items()};M[i]=[{m:v/scale for m,v in p.items()}for p in M[i]]
 for i in range(17):
  for j in range(24):assert multiply(multiply(rf[i],cf[j]),M[i][j])==raw[i][j]
 unitdet={0:Q(1)}
 for i in rows:unitdet=multiply(unitdet,rf[i])
 for j in cols:unitdet=multiply(unitdet,cf[j])
 expected,rem=divide(detsource,unitdet);assert not rem and all(v.denominator==1 for v in expected.values());assert multiply(unitdet,expected)==detsource
 removed=[sum(rexp[i][k]for i in rows)+sum(cexp[j][k]for j in cols)for k in range(6)];oldpowers=[x['power']for x in u['all6_known_unit_power_divisions']];remaining=[a-b for a,b in zip(oldpowers,removed)];assert all(x>=0 for x in remaining)
 rawcost=costs(raw,rows,cols,extras);preparedcost=costs(M,rows,cols,extras);assert preparedcost['all119_replacement_coeff_slot_upper']<150000 and preparedcost['single17vector_raw_byte_upper_without_additional_Bezout']<32*1024*1024
 cert={'variable':'zeta','source_row_order_for_A_B':rows,'A_source_columns':cols,'B_source_columns':extras,'source_matrix17x24':str(ip.relative_to(ROOT)),'six_known_chart_linear_units':[ser(p)for p in linear],'row_knownunit_exponents':rexp,'column_knownunit_exponents':cexp,'row_QQ_units':list(map(str,rowqq)),'full_expanded_row_units':[ser(p)for p in rf],'full_expanded_column_units':[ser(p)for p in cf],'all408_source_entry_reconstruction_checked':True,'prepared_integer_matrix17x24':[[ser(p)for p in row]for row in M],'prepared_A17x17':[[ser(M[i][j])for j in cols]for i in rows],'prepared_B17x7':[[ser(M[i][j])for j in extras]for i in rows],'expected_prepared_det_from_saved_source_NO_newdet':ser(expected),'whole_det_unit_multiplier':ser(unitdet),'source_det_equals_unitdet_times_expected_checked':True,'old_known_unit_det_powers':oldpowers,'removed_row_column_det_unit_powers':removed,'remaining_known_unit_det_powers':remaining,'raw_cost':rawcost,'prepared_cost':preparedcost,'R46_irreducible_exponent1_preserved_by_knownunit_scalings':True,'future_solver_gate':'RequireA*X=delta*B,deltaequals(plus/minus)expected; FIRSTXijcoprimeR46plus17vectorsourceproof andQQ[z]BezoutwouldcoverR46; doNOTdumpall119redundantly.','all_X_zero_mod_R46_only_sufficient_projection_criterion_fails_not_full57_negative':True,'no_solve_rank_rightinverse_point_prime_GB':True}
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),base/'report.json',ip,np,ub/'report.json',up,ROOT/'reports/S61_NONSQUARE_THETA_P80_FACTOR_INPUT_20261104_v63/report.json'];r={'id':ID,'source_all408_entries_and_det_unit_identity_checked':True,'raw_cost':rawcost,'prepared_cost':preparedcost,'prepared_expected_det_degree':max(expected),'known_det_powers_removed':removed,'known_det_powers_remaining':remaining,'solve_rank_point_prime_GB_jobs':0,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-start,'frozen':True,'independent_review_pending':True};assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
