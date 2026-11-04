#!/usr/bin/env python3
"""Saved-plane trace identities and integer-scale/source-liveness preflight only."""
import argparse,gzip,hashlib,json,math,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
from weight1_structure import J,sp_basis
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_HESSIAN_SOURCE_VECTOR_PREFLIGHT_20261104_v63'
P=32003
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dictionary(v,d):return {(i,j):v[k]for k,(i,j)in enumerate((i,j)for i in range(d+1)for j in range(d-i+1))if v[k]}
def add(*terms):
 ans={}
 for scalar,a in terms:
  for m,z in a.items():ans[m]=(ans.get(m,0)+scalar*z)%P
 return {m:z for m,z in ans.items()if z}
def mul(a,b):
 ans={}
 for(i,j),x in a.items():
  for(k,l),y in b.items():ans[i+k,j+l]=(ans.get((i+k,j+l),0)+x*y)%P
 return {m:z for m,z in ans.items()if z}
def power(a,n):
 b={(0,0):1}
 for _ in range(n):b=mul(b,a)
 return b
def mm(a,b):return [[sum(a[i][k]*b[k][j]for k in range(len(b)))for j in range(len(b[0]))]for i in range(len(a))]
def trace_word_norm(C,D,w):
 a=[[int(i==j)for j in range(4)]for i in range(4)]
 for x in w:a=mm(a,[C,D][x])
 return sum(a[i][i]for i in range(4))
def main(out):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 basispath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/native-certificate.json.gz'
 cdpath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_UNIVERSAL_BASIS_20261104_v63/native-certificate.json.gz'
 prepath=ROOT/'reports/Q_NEW_FAMILY_DEGREE24_THETA_H_REPAIR_PREFLIGHT_20261104_v63/report.json'
 betapath=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63/native-input.txt.gz'
 spin,inv,basis,cd,pre=map(read,[spinpath,invpath,basispath,cdpath,prepath]);assert basis['status']=='complete_Inv24_basis_certified'
 # Universal weight identities, compared over QQ in TWO CARTAN variables only.
 a,b=s.symbols('a b');s1=a*a+b*b;s2=a*a*b*b;weights=inv['Cartan_weights']
 assert s.Poly(sum((u*a+v*b)**6 for u,v in weights)-(264*s1**3+408*s1*s2),a,b).is_zero
 assert s.Poly(sum((u*a+v*b)**8 for u,v in weights)-(1032*s1**4+3488*s1**2*s2-4208*s2**2),a,b).is_zero
 # SAME saved formal plane; no fresh point or invariant matrix evaluation.
 F=dictionary(cd['F4'],4);E=[dictionary(c,8)for c in cd['E5']];p0=add((pow(400,-1,P),E[0]),(-4,E[4]));S1=add((pow(2,-1,P),F))
 S2=add((1,p0),(-2,E[1]),(4,E[2]),(-8,E[3]),(16,E[4]))
 tr6=add((264,power(S1,3)),(408,mul(S1,S2)))
 tr8=add((1032,power(S1,4)),(3488,mul(power(S1,2),S2)),(-4208,power(S2,2)))
 assert tr6==dictionary(basis['old_trace6_columns'][0],12)
 assert tr8==dictionary(basis['old_trace8_columns'][0],16)
 # Moment tensor denominators and exact rho map are bound to accepted quadratic coefficients.
 O=s.Matrix(inv['symplectic_form16']);rho=[s.Matrix(r)for r in inv['all_ten_representation_matrices']];gs=sp_basis();Gi=s.Matrix([[s.trace(x*y)for y in gs]for x in gs]).inv()
 assert all(s.denom(x)==1 for x in O)and all(s.denom(x)==1 for r in rho for x in r)and all(s.denom(2*x)==1 for x in Gi)
 moment=spin['moment_transpose_coefficients16x136'];assert all((2*Q(str(x))).denominator==1 for r in moment for x in r)
 # For mu=J*S, coefficients in sp_basis are exactly S_ij; no extra division.
 for k,A in enumerate(gs):
  S=J.inv()*A;assert S==S.T
  ij=list((i,j)for i in range(4)for j in range(i,4))[k]
  assert S[ij[0],ij[1]]==1
 tok=iter(gzip.decompress(betapath.read_bytes()).decode().split());assert[int(next(tok))for _ in range(4)]==[4,16,136,1080];norm=[]
 for _ in range(16):
  n=int(next(tok));v=0
  for _ in range(n):v+=abs(int(next(tok)));[next(tok)for _ in range(4)]
  norm.append(v)
 mon=spin['quadratic_monomials'];mat=[]
 for key in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  vals=[sum(abs(4*Q(str(z)))*norm[mon[k][0]]*norm[mon[k][1]]for k,z in enumerate(row))for row in spin[key]]
  assert all(v.denominator==1 for v in vals);mat.append([list(map(int,vals[i:i+4]))for i in range(0,16,4)])
 bounds=[]
 for w in[(0,0,0,1,0,1),(0,0,0,1,1,1)]:
  z=trace_word_norm(*mat,w);support=27418;bounds.append(dict(kind='CD_trace6',word=w,coefficient_l1_upper_bits=z.bit_length(),coefficient_decimal_digits=len(str(z)),support_upper=support,raw_JSON_bytes_upper=support*(len(str(z))+40)))
 for n in[6,8]:
  old=pre['future_source_old_word_bounds'][0 if n==6 else 1];bits=old['coefficient_l1_upper_bit_bounds'][1]
  # The old bound used scale20. The ACTUAL new integer scale2 divides all coefficients by10^n.
  z=((1<<bits)-1+10**n-1)//10**n;support=old['output_support_bound']
  bounds.append(dict(kind='theta_H_trace',length=n,word=[0]*(n-1)+[1],integer_scale=2,coefficient_l1_upper_bits=z.bit_length(),coefficient_decimal_digits=len(str(z)),support_upper=support,raw_JSON_bytes_upper=support*(len(str(z))+40)))
 total=sum(v['raw_JSON_bytes_upper']for v in bounds)
 report=dict(id=STEM,all_passed=True,metadata_and_saved_coefficient_preflight_only=True,no_new_invariant_source_expansion=True,no_NEW25_points=True,
  universal_character_trace6_and8_identities_verified=True,same_saved_formal_plane_all_coefficients_trace6_and8_verified=True,
  integral_Omega_and_rho=True,twice_Gram_inverse_integral=True,twice_moment_quadratic_coefficients_integral=True,
  exact_integer_recipe='Theta=2theta(b), Z=rho(Aj)b,V=Omega Z,G=2GramInv; numerator(Tr theta^k H)=Tr Theta^(k+1)+Tr(G V^T Theta^k Z)',
  exact_word_denominators={'theta5H':'2^6*(1080BH)^12','theta7H':'2^8*(1080BH)^16','CD_trace6':'4^6*(1080BH)^12'},
  source_vector_recurrence_polynomial_products_upper=7*16*16*10,vector_live_levels=2,
  selected_output_bounds=bounds,total_compact_JSON_bytes_upper=total,rigorous_raw32MiB_completion_bound=total<=32*1024*1024,
  proposed_source_job_limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,output_terms=100000,intermediate_terms=150000,total_raw_coefficient_bytes=32*1024*1024),
  not_a_runtime_or_memory_completion_guarantee=True,no67_relation_test_authorized_by_this_preflight=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),spinpath,invpath,basispath,cdpath,prepath,betapath]},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in['all_passed','selected_output_bounds','total_compact_JSON_bytes_upper','rigorous_raw32MiB_completion_bound','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
