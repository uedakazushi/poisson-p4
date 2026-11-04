#!/usr/bin/env python3
"""One point-free 41-column weighted-degree24 coefficient test and exact lift."""
import argparse,gzip,hashlib,itertools,json,math,resource,shlex,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_DEGREE24_COEFFICIENT_RELATION_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def main(out,work):
 started=time.monotonic();sys.set_int_max_str_digits(0);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 expdir=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63';expreport=read(expdir/'report.json');assert expreport['all_passed']
 for name,digest in expreport['evidence_sha256'].items():assert sha(expdir/name)==digest,name
 names=['F4','p1','p2','p3','p4','N12'];polys=[read(expdir/(name+'.json.gz'))for name in names]
 exponents=[]
 for r in range(4):
  for choices in itertools.combinations_with_replacement(range(4),r):
   row=[6-2*r,0,0,0,0,0]
   for j in choices:row[1+j]+=1
   exponents.append(row)
 exponents.append([3,0,0,0,0,1])
 for j in range(4):
  row=[1,0,0,0,0,1];row[j+1]=1;exponents.append(row)
 exponents.append([0,0,0,0,0,2]);assert len(exponents)==41
 constden=[]
 for e in exponents:
  assert e[0]+2*sum(e[1:5])+3*e[5]==6
  assert sum(x*p['denominator_D_power']for x,p in zip(e,polys))==24
  assert sum(x*p['denominator_B_power']for x,p in zip(e,polys))<=12
  constden.append(math.prod(int(p['denominator_integer'])**x for x,p in zip(e,polys)))
 assert max(t[4]for t in polys[-1]['terms'])==1
 bounds=[sum(e[i]*max(t[4]for t in p['terms'])for i,p in enumerate(polys))+12-sum(e[i]*p['denominator_B_power']for i,p in enumerate(polys))for e in exponents]
 assert bounds[-1]==4 and max(bounds[:-1])<=3
 inp=work/'input.txt'
 with inp.open('w')as f:
  f.write('4 6 41 32003\n')
  for p in polys:
   assert int(p['denominator_integer'])%32003
   f.write(f"{p['denominator_integer']} {p['denominator_B_power']} {p['denominator_D_power']}\n{len(p['terms'])}\n")
   for coefficient,*mon in p['terms']:f.write(coefficient+' '+' '.join(map(str,mon))+'\n')
  for e in exponents:f.write(' '.join(map(str,e))+'\n')
 with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
 native=ROOT/'code/q_new_family_degree24_coefficient_relation_20261104_v63.cpp';binary=work/'relation-native'
 flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True));subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 raw=work/'native-certificate.json';error=None
 try:
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(raw)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=60,preexec_fn=cap)
 except Exception as exc:error=repr(exc)
 data=read(raw)if error is None else None;verified=False;relation=[]
 if data is not None:
  assert data['all_passed']and 0<=data['coefficient_rank_modp']<=41 and data['max_polynomial_terms']<=100000
  assert len(data['selected_rows'])==data['coefficient_rank_modp']
  if data['coefficient_rank_modp']==40:
   kernel=list(map(int,data['integer_column_kernel']));G=list(map(int,data['universal_primitive_integer_relation']))
   assert len(kernel)==len(G)==41 and any(kernel)and any(G)
   assert all(sum(int(z)*k for z,k in zip(row,kernel))==0 for row in data['selected_integer_matrix'])
   common=next((k*d//g for k,d,g in zip(kernel,constden,G)if g),None)
   assert common and all(k*d==common*g for k,d,g in zip(kernel,constden,G))
   assert math.gcd(*G)==1
   verified=data['exact_integer_coefficient_identity_verified']
   if verified:assert G[-1]==0
   relation=[dict(coefficient=str(g),exponents=e)for g,e in zip(G,exponents)if g]
  with gzip.open(out/'native-certificate.json.gz','wb')as f:f.write(raw.read_bytes())
  assert raw.stat().st_size<=32*1024*1024
 sample=dict(generator_order=names,weights_beta=[4,8,8,8,8,12],weighted_monomial_exponents=exponents,
  column_constant_denominators=list(map(str,constden)),common_denominator='B^12*D^24',
  polynomial_column_B_degree_upper_bounds=bounds,N12_squared_coefficient_forced_zero_by_unique_B_degree4=True,
  universal_relation_terms=relation,coefficient_identity_verified=verified)
 (out/'relation.json').write_text(json.dumps(sample,indent=2)+'\n')
 gatepaths=[ROOT/'reports/Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63/report.json',
  ROOT/'reports/Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE_20261103_v63/report.json',
  ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/report.json',
  ROOT/'reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63/report.json']
 assert all(read(p)['all_passed']for p in gatepaths)
 inputs=[Path(__file__),native,expdir/'report.json',*[expdir/(name+'.json.gz')for name in names],*gatepaths]
 report=dict(id=STEM,all_passed=error is None,status=data['status']if data else 'bounded_failure',
  prime=32003,one_GF_coefficient_rank_only=True,native_run_count=1,no_source_point_evaluations=True,
  coefficient_rank_modp=data['coefficient_rank_modp']if data else None,
  exact_all_polynomial_coefficients_identity_verified=verified,relation_support=len(relation),
  generator_six_algebraically_independent_universally=True,
  verified_relation_is_nonzero_universal_weight24=verified,
  exact_NEW25_firstjet_defining_equation_degree24_proved=verified,
  N12_squared_forced_zero_by_B_degrees=True,no_higher_degree_search_or_automatic_retry=True,
  limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,terms_per_poly=100000,certificate_bytes=32*1024*1024),
  error=error,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','status','coefficient_rank_modp','relation_support','exact_NEW25_firstjet_defining_equation_degree24_proved','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
 a=p.parse_args();main(a.output.resolve(),a.work.resolve())
