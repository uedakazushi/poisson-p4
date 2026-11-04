#!/usr/bin/env python3
"""ONE capped exact four-output expansion with T^2=B torus preprocessing."""
import argparse,gzip,hashlib,json,resource,shlex,subprocess,time
from fractions import Fraction
from pathlib import Path
import sympy as s
from weight1_structure import J,sp_basis
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63'
NAMES=['CD000101','CD000111','theta5H','theta7H']
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def main(out,work):
 started=time.monotonic();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);raw=work/'raw';raw.mkdir()
 prepath=ROOT/'reports/Q_NEW_FAMILY_HESSIAN_SOURCE_VECTOR_PREFLIGHT_20261104_v63/report.json'
 basispath=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/report.json'
 oldpath=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63/native-input.txt.gz'
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 pre,basis,inv,spin=map(read,[prepath,basispath,invpath,spinpath]);assert pre['all_passed']and pre['same_saved_formal_plane_all_coefficients_trace6_and8_verified']and basis['universal_complete_QQ_Inv24_basis']
 O=s.Matrix(inv['symplectic_form16']);rho=[s.Matrix(r)for r in inv['all_ten_representation_matrices']];gs=sp_basis();Gi=s.Matrix([[s.trace(a*b)for b in gs]for a in gs]).inv()
 assert all(z.is_Integer for z in O)and all(z.is_Integer for r in rho for z in r)and all((2*z).is_Integer for z in Gi)
 # Verify exact torus cocharacter from the original physical representation.
 h2=s.diag(0,0,1,-1);coef=J.inv()*h2;Hrep=s.zeros(16);k=0
 for i in range(4):
  for j in range(i,4):Hrep+=coef[i,j]*rho[k];k+=1
 weights=[int(w[1])for w in inv['Cartan_weights']];assert Hrep==s.diag(*weights)
 # The accepted mu^T=-C+2D tensor is half-integral, hence2mu and2theta integral.
 assert all((2*Fraction(str(z))).denominator==1 for row in spin['moment_transpose_coefficients16x136']for z in row)
 tokens=iter(gzip.decompress(oldpath.read_bytes()).decode().split());assert [int(next(tokens))for _ in range(4)]==[4,16,136,1080]
 polys=[];transforms=[]
 for i in range(16):
  n=int(next(tokens));terms=[]
  for _ in range(n):
   z=int(next(tokens));exponents=[int(next(tokens))for _ in range(4)];oldB=exponents[3];exponents[3]=2*oldB+weights[i]-1
   assert exponents[3]>=0;terms.append((z,exponents))
  if terms:assert max(m[3]for z,m in terms)<=1
  polys.append(terms);transforms.append(dict(index=i,second_Cartan_weight=weights[i],input_terms=n,T_exponents=sorted({m[3]for z,m in terms})))
 mon=[[int(next(tokens)),int(next(tokens))]for _ in range(136)];assert mon==spin['quadratic_monomials']
 coeff=[[int(next(tokens))for _ in range(136)]for _ in range(32)];assert next(tokens,None)is None
 expected=[]
 for key in['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
  for row in spin[key]:
   q=[4*Fraction(str(z))for z in row];assert all(z.denominator==1 for z in q);expected.append(list(map(int,q)))
 assert coeff==expected
 inp=work/'native-input.txt'
 with inp.open('w')as f:
  f.write('4 16 136 1080\n')
  for terms in polys:
   f.write(str(len(terms))+'\n')
   for z,m in terms:f.write(str(z)+' '+' '.join(map(str,m))+'\n')
  for row in mon+coeff:f.write(' '.join(map(str,row))+'\n')
  for mat in[O,2*Gi,*rho]:f.write(' '.join(str(int(z))for z in mat)+'\n')
 with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
 native=ROOT/'code/q_new_family_hessian_low_invariant_expand_20261104_v63.cpp';helper=ROOT/'code/q_new_family_low_invariant_expand_20261104_v63.cpp';binary=work/'source-native'
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 error=None
 try:
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(raw)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=60,preexec_fn=cap)
 except Exception as e:error=repr(e)
 completed=[]
 for name in NAMES:
  path=raw/(name+'.json')
  if not path.is_file():continue
  try:d=read(path)
  except Exception:continue
  degree=16 if name=='theta7H'else 12;seen=set();terms=d['terms']
  assert d['name']==name and d['T_squared_equals_B_preconditioned']and d['denominator_D_power']==degree
  assert len(terms)<=100000 and d['source_c_power']==degree//2 and 0<=d['denominator_B_power']<=degree//2
  assert int(d['denominator_integer'])>0 and int(d['denominator_integer'])%32003
  bits=0
  for z,*m in terms:
   assert len(m)==4 and tuple(m)not in seen;seen.add(tuple(m));assert int(z);bits=max(bits,abs(int(z)).bit_length())
   assert all(isinstance(k,int)and k>=0 for k in m)and m[3]<=degree//2
  with gzip.open(out/(name+'.json.gz'),'wb')as f:f.write(path.read_bytes())
  completed.append(dict(name=name,terms=len(terms),max_coefficient_bits=bits,raw_bytes=path.stat().st_size,compressed_bytes=(out/(name+'.json.gz')).stat().st_size,
   denominator_integer=d['denominator_integer'],denominator_B_power=d['denominator_B_power'],denominator_H_power=degree))
 nr=read(raw/'native-report.json')if error is None else None
 if nr:
  assert nr['all_passed']and len(completed)==4 and nr['max_intermediate_terms']<=150000
  assert sum(v['raw_bytes']for v in completed)==nr['coefficient_output_bytes']<=32*1024*1024
  assert nr['all_T_parity_and_u_X_parity_verified']and nr['integer_scale_theta']==2
  (out/'native-report.json').write_text(json.dumps(nr,indent=2)+'\n')
 recipe=dict(old_integer_beta='b=1080*B*H*beta at c=1',source_H='147456*X^2-(3072W-1)^2',preprocessed_beta='bbar_i=T^(w_i-1)*b_i(B=T^2)',
  invariant_transformation='I(bbar)=B^(-d/2)*I(b); I(beta)=I(bbar)/((1080H)^d*B^(d/2))',
  exact_torus_cocharacter_verified=True,per_coordinate_transform=transforms,
  integer_vector_formula='Theta=2theta(bbar);Z=rho(Aj)bbar;V=Omega Z;G=2GramInv; numerator=TrTheta^(k+1)+Tr(G V^T Theta^k Z)',
  pure_trace6='264s1^3+408s1s2',pure_trace8='1032s1^4+3488s1^2s2-4208s2^2',s1='tr((2mu)^2)/2',s2='det(2mu)',
  output_conversion='u^a X^x T^t -> (u/X)^a*(X^2)^((a+x)/2)*B^(t/2); all exponents even where required')
 (out/'recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
 inputs=[Path(__file__),native,helper,prepath,basispath,oldpath,invpath,spinpath,ROOT/'code/weight1_structure.py']
 report=dict(id=STEM,all_passed=error is None,one_native_expansion_run=True,native_run_count=1,completed_outputs=completed,error=error,
  source_c_algebraically_set1_and_restored_by_homogeneity=True,T_squared_equals_B_preconditioning_exact=True,
  same_saved_formal_plane_word_normalizations_preverified=True,source_vector_recursion_two_live_levels=True,
  native_hard_limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,coefficient_bytes=32*1024*1024,output_terms=100000,intermediate_terms=150000),
  coefficient_output_bytes=sum(v['raw_bytes']for v in completed),compressed_coefficient_bytes=sum(v['compressed_bytes']for v in completed),native_report=nr,
  no_new_actual_point_or_degree24_coefficient_test=True,no_automatic_retry_or_budget_growth=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','completed_outputs','coefficient_output_bytes','seconds','error']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();main(args.output.resolve(),args.work.resolve())
