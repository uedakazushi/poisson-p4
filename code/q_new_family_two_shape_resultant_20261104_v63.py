#!/usr/bin/env python3
"""ONE shared120s exact resultant job, after a LIGHT named-content preflight."""
import argparse,gzip,hashlib,json,math,resource,shlex,signal,subprocess,time
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_RESULTANT_20261104_v63';BYTECAP=32*1024*1024
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for z in iter(lambda:f.read(1<<20),b''):h.update(z)
 return h.hexdigest()
def cap(seconds):
 resource.setrlimit(resource.RLIMIT_CPU,(seconds,seconds));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def main(out,work):
 started=time.monotonic();deadline=started+120;resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 path=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_HOMOGENIZATION_20261104_v63/certificate.json';review=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_RECOVERY_20261104_v63_REVIEW/report.json'
 data=read(path);gate=read(review);assert gate['all_passed']
 for n,h in gate['input_sha256'].items():assert sha(ROOT/n)==h,n
 r,s,R,Z0,Z1,Z2=S.symbols('r s R Z0 Z1 Z2');v3=(r,s,R);v6=(s,r,R,Z0,Z1,Z2);P=lambda z:S.Poly(z,*v3,domain=S.QQ)
 def load(rec):return S.Poly.from_dict({tuple(e):S.Rational(a,b)for a,b,*e in rec['coefficients']},v3,domain=S.QQ)
 original=[]
 for name in['Z1','Z2']:
  c=data['localized_constraints'][name];original.append([load(c[k])for k in['coefficient1','coefficient2','denominator']])
 Fn=P(S.sympify(data['Fn']));J=P(S.sympify(data['J']));preflight={};native_data=None;error=None;calls=0;front=[]
 def checkpoint(stage):
  front.append(dict(stage=stage,seconds=time.monotonic()-started));(out/'frontier.json').write_text(json.dumps(front,indent=2)+'\n')
 try:
  # Factor ONLY the small known denominator Fn; no input-equation factorization.
  signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('light Fn factorization5s cap')));signal.alarm(5)
  constant,ff=S.factor_list(Fn);signal.alarm(0);product=P(constant)
  for f,e in ff:product*=f**e
  assert product==Fn
  factors=[('r',P(r)),('r-1',P(r-1)),('J',J)]+[('Fn_factor'+str(i),f)for i,(f,e)in enumerate(ff)]
  # Every Fn factor is a unit on the PREVIOUS Fn!=0 chart, including repeated factors.
  cleared=[];metadata=[]
  for i,group in enumerate(original):
   group=list(group);removed=[];deleted=P(1)
   for label,f in factors:
    count=0
    while True:
     divisions=[p.div(f)for p in group]
     if not all(rem.is_zero for quo,rem in divisions):break
     old=group;group=[quo for quo,rem in divisions]
     assert all(a*f==b for a,b in zip(group,old));deleted*=f;count+=1
    if count:removed.append(dict(factor=label,power=count,expression=str(f.as_expr())))
   # Independent target variables mean ALL3 coefficient polynomials must share a deletedfactor.
   for label,f in factors:assert not all(p.rem(f).is_zero for p in group)
   expr=group[0].as_expr()+Z0*group[1].as_expr()-[Z1,Z2][i]*group[2].as_expr();pp=S.Poly(expr,*v6,domain=S.QQ)
   lcm=math.lcm(*(int(S.denom(c))for c in pp.coeffs()));integer={m:int(c*lcm)for m,c in pp.terms()};g=math.gcd(*(abs(c)for c in integer.values()));integer={m:c//g for m,c in integer.items()};scale=S.Rational(lcm,g)
   integerpoly=S.Poly.from_dict(integer,v6,domain=S.ZZ);before=original[i][0].as_expr()+Z0*original[i][1].as_expr()-[Z1,Z2][i]*original[i][2].as_expr()
   assert S.Poly(integerpoly.as_expr()*deleted.as_expr()-scale*before,*v6,domain=S.QQ).is_zero
   assert integerpoly.degree(s)>0
   cleared.append(integerpoly);metadata.append(dict(input_index=i,removed_parameter_content=removed,deleted_product=str(deleted.as_expr()),
    clearing_rational=str(scale),terms=len(integer),multidegree=list(integerpoly.degree_list()),
    exact_input_identity='integer_input*deleted_product=clearing_rational*original_constraint',
    no_known_denominator_factor_common_to_all3_remaining_coefficients=True))
  preflight=dict(Fn_factorization_constant=str(constant),Fn_factors=[dict(expression=str(f.as_expr()),multiplicity=int(e))for f,e in ff],
   Fn_squarefree=all(e==1 for f,e in ff),Fn_factorization_exact_remultiplication=True,inputs=metadata,
   retained_localizers=data['scope'],generic_finiteness_dependency='accepted octic dominance; independent Z0,R,Z1,Z2 gives finite rational source-field extension in(r,s)',
   no_generic_chart_component_left_reason='every irreducible factor of knownDenr(r−1)JFn tested in ALL independenttarget coefficients; nongeneric denominator loci remain explicitly excluded')
  (out/'content-preflight.json').write_text(json.dumps(preflight,indent=2)+'\n');checkpoint('light_content_preflight_passed')
  inp=work/'integer-input.txt'
  with inp.open('w')as f:
   f.write('6 0\n')
   for p in cleared:
    f.write(str(len(p.terms()))+'\n')
    for e,c in p.terms():f.write(str(c)+' '+' '.join(map(str,e))+'\n')
  with gzip.open(out/'integer-input.txt.gz','wb')as f:f.write(inp.read_bytes())
  m,n=[p.degree(s)for p in cleared];coeffs=[]
  for p in cleared:
   table=[]
   for k in range(p.degree(s)+1):
    terms=[[str(c),*e[1:]]for e,c in p.terms()if e[0]==k];table.append(terms)
   coeffs.append(table)
  matrix=[[None]*(m+n)for _ in range(m+n)]
  for row in range(n):
   for k in range(m+1):matrix[row][row+k]=dict(input_index=0,s_power=m-k)
  for row in range(m):
   for k in range(n+1):matrix[n+row][row+k]=dict(input_index=1,s_power=n-k)
  semantic=dict(variable_order=list(map(str,v6)),eliminate_variable='s',integer_clearing=metadata,s_degrees=[m,n],
   coefficient_variable_order=['r','R','Z0','Z1','Z2'],input_s_coefficients_ascending=coeffs,Sylvester_matrix_entry_references=matrix,
   determinant_definition='n shifted rows of descendingE1 coefficients, thenm shifted rows of descendingE2; det=Res_s(E1,E2)',
   source_constraint_sha256=sha(path),raw_resultant_is_NOT_NEW25_equation=True)
  (out/'determinant-recipe.json').write_text(json.dumps(semantic,indent=2)+'\n')
  native=ROOT/'code/q_new_family_two_shape_resultant_20261104_v63.cpp';binary=work/'resultant-native'
  flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint','openssl'],text=True));subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(.1,deadline-time.monotonic()-3))
  remaining=deadline-time.monotonic()-3;assert remaining>0
  existing=sum(p.stat().st_size for p in out.iterdir()if p.is_file());budget=BYTECAP-existing-1_000_000;assert budget>0
  calls=1;checkpoint('starting_one_exact_FLINT_resultant_shared_deadline')
  with(out/'native.log').open('w')as log:
   subprocess.run([str(binary),str(inp),str(out),str(budget)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=remaining,preexec_fn=lambda:cap(max(1,math.ceil(remaining))))
  native_data=read(out/'native-report.json');assert native_data['resultant_call_count']==1
  if native_data['full_coefficients_saved']:
   result=out/'resultant.txt';hsh=hashlib.sha256();count=0
   with result.open('rb')as f:
    header=f.readline().split();assert header==[b'6',str(native_data['terms']).encode()]
    for line in f:hsh.update(line);count+=1
   assert hsh.hexdigest()==native_data['canonical_term_SHA256']and count==native_data['terms']
   assert native_data['multidegree'][0]==0
   # Degree bounds of the defining Sylvester determinant, no identity rerun.
   for j in range(1,6):assert native_data['multidegree'][j]<=n*cleared[0].degree_list()[j]+m*cleared[1].degree_list()[j]
  checkpoint(native_data['status'])
 except Exception as exc:error=repr(exc);signal.alarm(0);checkpoint('bounded_failure_or_preflight_stop')
 inputs=[Path(__file__),ROOT/'code/q_new_family_two_shape_resultant_20261104_v63.cpp',path,review]
 final=dict(id=STEM,all_passed=error is None and bool(native_data and native_data['all_passed']),error=error,
  status=native_data['status']if native_data else 'bounded_failure_or_preflight_stop',resultant_call_count=calls,
  preflight_completed=bool(preflight),native_result=native_data,one_shared_deadline_seconds=120,address_space_bytes=2_000_000_000,
  no_automatic_retry_alternate_algorithm_newpoint_prime_source26_or_r_elimination=True,
  FLINT_library_return_status_distinguished_from_zero_polynomial=True,library_output_requires_later_independent_identity_audit=True,
  raw_resultant_NOT_NEW25_defining_equation=True,actual_target_index_and_NEW_equation_unresolved=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<BYTECAP-10000
 (out/'report.json').write_text(json.dumps(final,indent=2)+'\n')
 print(json.dumps({k:final[k]for k in['all_passed','status','resultant_call_count','error','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();main(args.output.resolve(),args.work.resolve())
