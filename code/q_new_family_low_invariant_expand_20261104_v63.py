#!/usr/bin/env python3
"""One capped exact low-invariant expansion after the accepted preflight."""
import argparse,gzip,hashlib,json,resource,shlex,subprocess,sys,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cap():
 resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def main(out,work):
 started=time.monotonic();sys.set_int_max_str_digits(0);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False);raw=work/'raw';raw.mkdir()
 prepath=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_PREFLIGHT_20261104_v63/report.json'
 sourcepath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
 sectionpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json'
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 pre,src,section,spin=map(read,[prepath,sourcepath,sectionpath,spinpath]);assert pre['all_passed']and pre['predicted_bound_within32MiB']
 for name,digest in pre['input_sha256'].items():assert sha(ROOT/name)==digest,name
 u,c,X,W,B,v=s.symbols('u c X W B v');D=147456*X**2-(3072*W-1)**2;vs=s.sympify(section['section_v'])
 polys=[s.Poly(s.cancel(B*D*s.sympify(z).subs({c:1,v:vs},simultaneous=True)),u,X,W,B)for z in src['primitive_beta16']]
 L=s.ilcm(1,1,*(s.denom(z)for p in polys for z in p.coeffs()));assert L==1080
 polys=[s.Poly(p.as_expr()*L,u,X,W,B)for p in polys];native=ROOT/'code/q_new_family_low_invariant_expand_20261104_v63.cpp'
 inp=work/'input.txt'
 with inp.open('w')as f:
  f.write('4 16 136 1080\n')
  for p in polys:
   terms=[(m,int(z))for m,z in p.terms()if z];f.write(str(len(terms))+'\n')
   for m,z in terms:f.write(str(z)+' '+' '.join(map(str,m))+'\n')
  for i,j in spin['quadratic_monomials']:f.write(f'{i} {j}\n')
  for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
   for row in spin[key]:
    co=[4*s.Rational(x)for x in row];assert all(x.is_Integer for x in co);f.write(' '.join(map(str,co))+'\n')
 with gzip.open(out/'native-input.txt.gz','wb')as f:f.write(inp.read_bytes())
 binary=work/'expand-native';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
 subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
 error=None
 try:
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(raw)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=60,preexec_fn=cap)
 except Exception as exc:error=repr(exc)
 completed=[]
 for name in ['F4','p0','p1','p2','p3','p4','N12']:
  path=raw/(name+'.json')
  if not path.is_file():continue
  try:data=read(path)
  except Exception:continue
  terms=data['terms'];d={'F4':4,'N12':12}.get(name,8)
  assert data['name']==name and len(terms)<=100000 and int(data['denominator_integer'])>0
  assert data['denominator_D_power']==d and 0<=data['denominator_B_power']<=d
  seen=set();bits=0
  for coefficient,*mon in terms:
   assert len(mon)==4 and tuple(mon)not in seen;seen.add(tuple(mon));z=int(coefficient);assert z;bits=max(bits,abs(z).bit_length())
   aa,YY,ww,bb=mon;assert 0<=aa<=d and 0<=YY<=2*d and 0<=ww<=3*d and 0<=bb<=d
  with gzip.open(out/(name+'.json.gz'),'wb')as f:f.write(path.read_bytes())
  completed.append(dict(name=name,terms=len(terms),max_coefficient_bits=bits,raw_bytes=path.stat().st_size,
   compressed_bytes=(out/(name+'.json.gz')).stat().st_size,denominator_integer=data['denominator_integer'],
   denominator_B_power=data['denominator_B_power'],denominator_D_power=data['denominator_D_power']))
 success=error is None
 nr=read(raw/'native-report.json')if success else None
 if success:
  assert nr['all_passed']and len(completed)==7 and nr['max_intermediate_terms']<=150000
  assert sum(z['raw_bytes']for z in completed)==nr['coefficient_output_bytes']<=32*1024*1024
  assert nr['quartic_relations_and_F4_square_pencil_identity_verified']
  (out/'native-report.json').write_text(json.dumps(nr,indent=2)+'\n')
 inputs=[Path(__file__),native,prepath,sourcepath,sectionpath,spinpath]
 report=dict(id=STEM,all_passed=success,one_native_expansion_run=True,native_run_count=1,
  native_hard_limits=dict(cpu_and_wall_seconds=60,address_space_bytes=2_000_000_000,coefficient_bytes=32*1024*1024,output_terms=100000,intermediate_terms=150000),
  shape_variables=['a=u/X','Y=X^2','W','B'],D='147456*Y-(3072*W-1)^2',source_c_algebraically_set1_and_restored_by_homogeneity=True,
  exact_native_word_scales=dict(F4='16*(1080BD)^4',p_i='256*(1080BD)^8',N12='32768*(1080BD)^12'),
  B_monomial_and_integer_content_cancelled=True,D_factors_not_cancelled=True,completed_outputs=completed,
  coefficient_output_bytes=sum(z['raw_bytes']for z in completed),compressed_coefficient_bytes=sum(z['compressed_bytes']for z in completed),
  native_report=nr,error=error,no_new_beta_point_no_Inv24_evaluation_matrix=True,
  no_41_relation_search_launched=True,no_automatic_retry_or_budget_enlargement=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:report[k]for k in['all_passed','completed_outputs','coefficient_output_bytes','seconds']},indent=2))
 if not success:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
 a=p.parse_args();main(a.output.resolve(),a.work.resolve())
