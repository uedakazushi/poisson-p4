#!/usr/bin/env python3
"""Bounds only: no F4, pencil, N12 coefficient expansion and no point evaluation."""
import argparse,gzip,hashlib,itertools,json,math,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_LOW_INVARIANT_PREFLIGHT_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 sourcepath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
 sectionpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json'
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 src,section,spin=map(read,[sourcepath,sectionpath,spinpath]);u,c,X,W,B=s.symbols('u c X W B');v=s.Symbol('v')
 D=147456*X**2-(3072*W-1)**2;vs=s.sympify(section['section_v']);assert s.cancel(vs+98304*X/D)==0
 polys=[s.Poly(s.cancel(B*D*s.sympify(z).subs({c:1,v:vs},simultaneous=True)),u,X,W,B)for z in src['primitive_beta16']]
 denominator=s.ilcm(1,1,*(s.denom(z)for p in polys for z in p.coeffs()))
 integer=[s.Poly(p.as_expr()*denominator,u,X,W,B)for p in polys]
 norms=[sum(abs(int(z))for z in p.coeffs())for p in integer]
 maxes=[max((m[k]for p in integer for m,z in p.terms()if z),default=0)for k in range(4)]
 total=max(sum(m)for p in integer for m,z in p.terms()if z)
 ux=max(m[0]+m[1]for p in integer for m,z in p.terms()if z)
 target2=max(3*m[0]+m[1]+2*m[2]+2*m[3]for p in integer for m,z in p.terms()if z)
 assert maxes==[1,3,3,1]and total==4 and ux==4 and target2==9
 mon=spin['quadratic_monomials']
 def adj_norm(key):
  values=[]
  for row in spin[key]:
   z=sum(abs(4*s.Rational(a))*norms[mon[k][0]]*norms[mon[k][1]]for k,a in enumerate(row));assert z.is_Integer;values.append(int(z))
  return [values[4*i:4*i+4]for i in range(4)]
 Cn,Dn=map(adj_norm,['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136'])
 def mm(A,B):return [[sum(A[i][k]*B[k][j]for k in range(4))for j in range(4)]for i in range(4)]
 def tr(A):return sum(A[i][i]for i in range(4))
 CC,DD=mm(Cn,Cn),mm(Dn,Dn);M=[[Cn[i][j]+2*Dn[i][j]for j in range(4)]for i in range(4)]
 Fbound=tr(mm(M,M));pbounds=[0]*5
 for perm in itertools.permutations(range(4)):
  coeff=[1]
  for i in range(4):
   old=coeff;coeff=[0]*(len(old)+1)
   for j,z in enumerate(old):coeff[j]+=z*Cn[i][perm[i]];coeff[j+1]+=z*Dn[i][perm[i]]
  pbounds=[a+b for a,b in zip(pbounds,coeff)]
 CD2=mm(Cn,DD);qC,qD,D4,T4=tr(CC),tr(DD),tr(mm(DD,DD)),tr(mm(CC,DD))
 Nbound=8*tr(mm(CD2,CD2))+4*qD*T4+qD*qD*qC+2*D4*qC
 def monomial_bound(degree):
  count=0
  for aa in range(degree+1):
   for YY in range(2*degree+1):
    if not 0<=2*YY-aa<=3*degree:continue
    for ww in range(3*degree+1):
     for bb in range(degree+1):
      if 2*YY+ww+bb<=4*degree and 2*(aa+YY+ww+bb)<=9*degree:count+=1
  return count
 counts={d:monomial_bound(d)for d in[4,8,12]}
 bounds=[Fbound,*pbounds,Nbound];degrees=[4,*([8]*5),12]
 records=[]
 for name,d,b in zip(['F4',*[f'p{i}'for i in range(5)],'N12'],degrees,bounds):
  # 40 extra bytes/term covers quoted integer, commas, four exponents, whitespace.
  digits=len(str(b));size=counts[d]*(digits+40)
  records.append(dict(output=name,beta_degree=d,known_denominator_power_B_and_D=d,
   numerator_max_total_degree_in_u_X_W_B=4*d,numerator_max_total_degree_after_parity=9*d//2,
   monomial_count_rigorous_polytope_upper_bound=counts[d],max_integer_coefficient_bits=b.bit_length(),
   max_integer_coefficient_decimal_digits=digits,compact_sparse_JSON_byte_upper_bound=size))
 size=sum(z['compact_sparse_JSON_byte_upper_bound']for z in records)
 report=dict(id=STEM,all_passed=True,only_sixteen_small_beta_numerators_expanded=True,
  no_invariant_output_or_point_expanded_or_evaluated=True,clearing_beta_factor=f'{denominator}*B*D',
  D=str(D),shape_variables=['a=u/X','Y=X^2','W','B'],
  input_beta_numerator_term_counts=[len(p.terms())if not p.is_zero else 0 for p in integer],
  global_beta_numerator_exponent_box=maxes,polytope_count_by_beta_degree=counts,output_preflight=records,
  total_compact_sparse_output_byte_upper_bound=size,requested_output_hard_cap_bytes=32*1024*1024,
  predicted_bound_within32MiB=size<=32*1024*1024,
  proposed_hard_runtime_seconds=60,proposed_hard_AS_bytes=2_000_000_000,
  proposed_term_cap_each_output=100000,proposed_serialize_stop_bytes=32*1024*1024,
  no_actual_output_expansion_launched=True,no_automatic_budget_enlargement=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),sourcepath,sectionpath,spinpath]},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','polytope_count_by_beta_degree','total_compact_sparse_output_byte_upper_bound','predicted_bound_within32MiB','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
