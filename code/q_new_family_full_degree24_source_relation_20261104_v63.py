#!/usr/bin/env python3
"""ONE full67 coefficient test after BOTH independent basis/source gates PASS."""
import argparse,gzip,hashlib,itertools,json,math,resource,shlex,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_FULL_DEGREE24_SOURCE_RELATION_20261104_v63'
P=32003
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def detmod(A):
 A=[r[:]for r in A];value=1
 for k in range(len(A)):
  j=next(i for i in range(k,len(A))if A[i][k]%P)
  if j!=k:A[j],A[k]=A[k],A[j];value=-value
  pivot=A[k][k]%P;value=value*pivot%P;inverse=pow(pivot,-1,P)
  for i in range(k+1,len(A)):
   q=A[i][k]*inverse%P
   if q:A[i]=[(a-q*b)%P for a,b in zip(A[i],A[k])]
 return value%P
def main(out,work):
 sys.set_int_max_str_digits(0);started=time.monotonic();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 gatepaths=[ROOT/'reports'/stem/'report.json'for stem in[
  'Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW',
  'Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW',
  'Q_NEW_FAMILY_FULL_DEGREE24_LOGIC_20261104_v63_REVIEW']]
 for p in gatepaths:
  d=read(p);assert d['all_passed']
  for n,h in d['input_sha256'].items():assert sha(ROOT/n)==h,n
  for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n
 basisdir=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63'
 basis=read(basisdir/'report.json');cert=read(basisdir/'native-certificate.json.gz');assert basis['universal_complete_QQ_Inv24_basis']
 olddir=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63';newdir=ROOT/'reports/Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63'
 names=['F4','p1','p2','p3','p4','N12','CD000101','CD000111','theta5H','theta7H'];degrees=[4,8,8,8,8,12,12,12,12,16]
 paths=[olddir/(n+'.json.gz')for n in names[:6]]+[newdir/(n+'.json.gz')for n in names[6:]];polys=list(map(read,paths))
 all70=[]
 for r in range(4):
  for ij in itertools.combinations_with_replacement(range(4),r):
   e=[6-2*r]+[0]*9
   for j in ij:e[1+j]+=1
   all70.append(e)
 for i in range(4):
  e=[3]+[0]*9;e[5+i]=1;all70.append(e)
  for j in range(4):
   e=[1]+[0]*9;e[1+j]=1;e[5+i]=1;all70.append(e)
 for i,j in itertools.combinations_with_replacement(range(4),2):
  e=[0]*10;e[5+i]+=1;e[5+j]+=1;all70.append(e)
 for j in range(5):
  e=[0]*10;e[9]=1;e[0 if j==0 else j]=2 if j==0 else 1;all70.append(e)
 assert len(all70)==70
 selected=cert['degree24_selected_columns'];assert selected==[j for j in range(70)if j not in[62,68,69]]
 exponents=[all70[j]for j in selected];constden=[]
 for e in exponents:
  assert sum(a*d for a,d in zip(e,degrees))==24
  assert sum(a*p['denominator_D_power']for a,p in zip(e,polys))==24
  assert sum(a*p['denominator_B_power']for a,p in zip(e,polys))<=12
  constden.append(math.prod(int(p['denominator_integer'])**a for a,p in zip(e,polys)))
 inp=work/'native-input.txt'
 with inp.open('w')as f:
  f.write('4 10 67 32003\n')
  for p in polys:
   assert int(p['denominator_integer'])%P
   f.write(f"{p['denominator_integer']} {p['denominator_B_power']} {p['denominator_D_power']}\n{len(p['terms'])}\n")
   for z,*m in p['terms']:f.write(z+' '+' '.join(map(str,m))+'\n')
  for e in exponents:f.write(' '.join(map(str,e))+'\n')
 with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
 generators=dict(generator_order=names,weights_beta=degrees,selected_columns_from70=selected,all70_generator_exponents=all70,
  selected67_generator_exponents=exponents,column_constant_denominators=list(map(str,constden)),common_denominator='B^12*H^24',
  individual_generator_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
 (out/'generators.json').write_text(json.dumps(generators,indent=2)+'\n')
 native=ROOT/'code/q_new_family_full_degree24_source_relation_20261104_v63.cpp';binary=work/'full24-source-native'
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 raw=work/'native-certificate.json';frontier=out/'frontier.json';error=None;data=None
 try:
  with(out/'native.log').open('w')as log:subprocess.run([str(binary),str(inp),str(raw),str(frontier)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60,preexec_fn=cap)
  assert raw.stat().st_size<=32*1024*1024;data=read(raw);assert data['all_passed']
 except Exception as e:error=repr(e)
 verified=False;injective=False;relation=[];rank=None
 if data:
  rank=data['coefficient_rank_modp'];assert 0<=rank<=67 and data['max_intermediate_terms']<=150000
  assert len(data['column_term_counts'])==67 and max(data['column_term_counts'])<=100000
  A=[[int(z)for z in row]for row in data['selected_integer_matrix']];pc=data['selected_columns'];assert len(A)==len(pc)==rank and all(len(r)==67 for r in A)
  minor=detmod([[row[j]%P for j in pc]for row in A]);assert minor==data['selected_minor_modp']and minor
  injective=rank==67
  if rank==66:
   kernel=list(map(int,data['integer_column_kernel']));G=list(map(int,data['universal_primitive_integer_relation']));assert len(kernel)==len(G)==67 and any(kernel)and any(G)
   assert all(sum(a*b for a,b in zip(row,kernel))==0 for row in A)
   common=next(k*d//g for k,d,g in zip(kernel,constden,G)if g);assert common and all(k*d==common*g for k,d,g in zip(kernel,constden,G))and math.gcd(*G)==1
   verified=data['exact_integer_coefficient_identity_verified']
   if not verified:
    counter=list(map(int,data['counterrow_integer_coefficients']));assert len(counter)==67
    dot=sum(k*z for k,z in zip(kernel,counter));assert dot==int(data['countercoefficient'])and dot and data['augmented_QQ_rank67_by_kernel_pairing']
    injective=True
   relation=[dict(coefficient=str(g),generator_exponents=e,universal_basis_column=j)for j,(g,e)in enumerate(zip(G,exponents))if g]
  with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
 lowerpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63/report.json'
 imagepath=ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63/report.json';sourcepath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/report.json'
 assert all(read(p)['all_passed']for p in[lowerpath,imagepath,sourcepath])
 (out/'relation.json').write_text(json.dumps(dict(exact_source_identity_verified=verified,terms=relation,
  universal_nonzero_reason='fixed independent67 basis and nonzero relation coefficient vector'if verified else None),indent=2)+'\n')
 inputs=[Path(__file__),native,*gatepaths,basisdir/'report.json',basisdir/'native-certificate.json.gz',*paths,lowerpath,imagepath,sourcepath]
 report=dict(id=STEM,all_passed=error is None,status=data['status']if data else 'bounded_failure',prime=P,native_run_count=1,
  both_independent_basis_and_source_gates_passed_before_run=True,one_GF_polynomial_coefficient_rank_no_source_points=True,
  coefficient_rank_modp=rank,all_Inv24_QQ_restrictions_injective=injective,NEW25_defining_degree_lower_bound=26 if injective else 24,
  exact_all_integer_polynomial_coefficients_identity_verified=verified,relation_support=len(relation),NEW25_degree24_defining_equation_proved=verified,
  odd_degree25_absent_by_central_minus_identity=injective,no_higher_degree_search_or_other_prime_or_automatic_retry=True,
  native_limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,coefficient_certificate_bytes=32*1024*1024,column_terms=100000,intermediate_terms=150000),
  error=error,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','status','coefficient_rank_modp','all_Inv24_QQ_restrictions_injective','NEW25_defining_degree_lower_bound','exact_all_integer_polynomial_coefficients_identity_verified','relation_support','NEW25_degree24_defining_equation_proved','seconds','error']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();main(args.output.resolve(),args.work.resolve())
