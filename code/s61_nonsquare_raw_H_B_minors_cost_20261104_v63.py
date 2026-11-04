#!/usr/bin/env python3
"""Read-only coefficient preparation and35 max-minor SUPPORT bounds; no determinants/GB."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 gb=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63';gp=gb/'certificate.json';assert sha(gp)==read(gb/'report.json')['certificate_sha256'];g=read(gp)
 review=ROOT/'reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63_REVIEW/report.json';peer=read(review);assert peer['all_passed']and peer['input_sha256'][str(gp.relative_to(ROOT))]==sha(gp)
 R,a,e,d=ring('a,e,d',s.QQ);dec=lambda p:R.from_dict({tuple(m):s.QQ(c)for m,c in p});enc=lambda p:[[list(m),str(c)]for m,c in p.items()]
 known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
 original=[[dec(q)for q in row]for row in g['actual_block7x3']];A=[row[:]for row in original];assert sum(len(q)for row in A for q in row)==908
 def clean(values):
  first=next(q for q in values if q);c=first.LC;v=[q/c for q in values];powers=[]
  for f in known:
   n=0
   while all(not q.rem(f)for q in v):v=[q.exquo(f)for q in v];n+=1
   powers.append(n)
  unit=R.ground_new(c)
  for f,n in zip(known,powers):unit*=f**n
  assert all(unit*x==y for x,y in zip(v,values));return v,{'constant':str(c),'a_d_F_B_C_J_powers':powers,'polynomial':enc(unit)}
 rowunits=[]
 for i in range(7):A[i],u=clean(A[i]);rowunits.append(u)
 colunits=[]
 for j in range(3):
  vals,u=clean([A[i][j]for i in range(7)]);colunits.append(u)
  for i in range(7):A[i][j]=vals[i]
 # Only QQ column scaling clears rational coefficient denominators; no unknown inverse.
 D=[]
 for j in range(3):
  q=1
  for i in range(7):
   for c in A[i][j].values():q=math.lcm(q,Q(str(c)).denominator)
  D.append(q)
 integer=[[A[i][j]*D[j]for j in range(3)]for i in range(7)];assert all(Q(str(c)).denominator==1 for row in integer for p in row for c in p.values())
 for i in range(7):
  for j in range(3):assert dec(rowunits[i]['polynomial'])*integer[i][j]*dec(colunits[j]['polynomial'])==original[i][j]*D[j]
 supports=[[set(p)for p in row]for row in integer];paircache={}
 def pairset(i,j,k,l):
  key=(i,j,k,l)
  if key not in paircache:paircache[key]={tuple(a+b for a,b in zip(m,n))for m in supports[i][j]for n in supports[k][l]}
  return paircache[key]
 triples=[]
 for I in itertools.combinations(range(7),3):
  union=set();rawproductbound=0
  for perm in itertools.permutations(range(3)):
   p=pairset(I[0],perm[0],I[1],perm[1]);ss={tuple(a+b for a,b in zip(m,n))for m in p for n in supports[I[2]][perm[2]]};union.update(ss)
   rawproductbound+=len(supports[I[0]][perm[0]])*len(supports[I[1]][perm[1]])*len(supports[I[2]][perm[2]])
  triples.append({'rows7':list(I),'support_upper':len(union),'max_totaldegree_upper':max(map(sum,union),default=0),'raw_multiply_contributions_bound':rawproductbound,
   'monomial_support_upper':sorted(map(list,union),reverse=True)})
 coefficientbits=max(abs(Q(str(c)).numerator).bit_length()for row in integer for p in row for c in p.values());contributionmax=max(t['raw_multiply_contributions_bound']for t in triples)
 heightbound=3*coefficientbits+(contributionmax-1).bit_length();ground=[[i,j,str(A[i][j])]for i in range(7)for j in range(3)if A[i][j]and A[i][j].is_ground]
 native=out/'prepared-integer7x3.txt'
 with native.open('w')as f:
  f.write('7 3\n')
  for row in integer:
   for q in row:f.write(str(len(q))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in q.items())+'\n')
 certificate={'coefficient_domain':'FULLU=QQ[a,e,d,1/(adFBCJ)]','variable_order':['a','e','d'],'actual_original_block_rows14':g['all_other3_column_support_exact_rows'],'actual_original_block_columns14':g['other3_original_columns'],
  'known_row_unit_divisors':rowunits,'known_column_unit_divisors':colunits,'integer_column_multipliers_QQ_units':D,'prepared_integer7x3':[[enc(q)for q in row]for row in integer],
  'literal_all21_unit_preconditioning_reconstructions_verified':True,'all35_minor_SUPPORT_bounds_only':triples,'max_prepared_coefficient_integer_bits':coefficientbits,'conservative_minor_integer_coefficient_bit_upper':heightbound,
  'literal_ground_entries_after_knownunit_preconditioning':ground,
  'rank_and_localized_maxminor_ideal_unchanged_proof':'Row and column factors are only QQ constants and a,d,F,B,C,J, allunits onFULLU. Each maxminor changes by its3 rowunit product and common columnunit product; the localized ideal is unchanged.',
  'no_new_minor_coefficients_determinants_GB_fieldpoint_or_factorization':True,
  'candidate_next_method_NOT_executed':'Evaluate fixed35 determinants of this prepared integer7x3 using6 triple products and exactFLINT polynomial arithmetic, cap150k total support. Strip onlynamedunits; inspect actualsupport/gcd. If a small35-generator ideal remains, ONE separately authorized exactQQ sequential-knownunit saturation/Bezout may decide allFULLU rank3. Nonunit outcome needs explicit locus, not a conclusion about fullquartic239.'}
 (out/'certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n');paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),gb/'report.json',gp,review]
 report={'id':ID,'original_block_terms':908,'prepared_block_terms':sum(len(q)for row in integer for q in row),'prepared_max_entry_terms':max(len(q)for row in integer for q in row),'literal_ground_entry_count':len(ground),
  'all35_minor_support_upper_sum':sum(t['support_upper']for t in triples),'max_minor_support_upper':max(t['support_upper']for t in triples),'max_minor_degree_upper':max(t['max_totaldegree_upper']for t in triples),'max_integer_input_bits':coefficientbits,'max_minor_integer_coefficient_bits_upper':heightbound,
  'row_named_factor_powers':[u['a_d_F_B_C_J_powers']for u in rowunits],'column_named_factor_powers':[u['a_d_F_B_C_J_powers']for u in colunits],
  'minor_coefficients_computed':0,'GB_jobs':0,'newpoint_prime_or_rank_claim':False,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-start,'frozen':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
