#!/usr/bin/env python3
"""ONE60s point-free content gcd of15 coefficient groups; no degree claim before peer."""
import argparse,gzip,hashlib,json,math,resource,shlex,subprocess,time
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_TWO_SHAPE_CONTENT_GCD_20261104_v63';CAP=32*1024*1024
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def cap(seconds):
 resource.setrlimit(resource.RLIMIT_CPU,(seconds,seconds));resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
def main(out,work):
 start=time.monotonic();deadline=start+60;resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 qdir=ROOT/'reports/Q_NEW_FAMILY_TWO_SHAPE_BOUNDARY_STRIP_20261104_v63';qpath=qdir/'quotient.txt';qreport=qdir/'report.json';q=read(qreport);assert q['all_passed']and sha(qpath)==q['evidence_sha256']['quotient.txt']
 gates=[ROOT/'reports'/stem/'report.json'for stem in['Q_NEW_FAMILY_TWO_SHAPE_FINAL_NORM_IDENTITY_20261104_v63_REVIEW','Q_NEW_FAMILY_TWO_SHAPE_COPRIME_DEGREES_20261104_v63_REVIEW']]
 assert all(read(p)['all_passed']for p in gates)
 groups=defaultdict(dict)
 with qpath.open()as f:
  assert f.readline().split()==['6','6919']
  for line in f:
   c,*e=map(int,line.split());assert e[0]==0;label=(e[4],e[5]);mon=tuple(e[1:4]);assert mon not in groups[label];groups[label][mon]=c
 assert len(groups)==15 and sum(map(len,groups.values()))==6919
 order=[(0,3),(5,0)]+sorted((k for k in groups if k not in[(0,3),(5,0)]),key=lambda k:(len(groups[k]),k))
 assert all(all(m[2]==0 for m in groups[k])for k in order[:2])and len(groups[order[0]])!=len(groups[order[1]])
 metadata=[]
 for k in order:
  g=groups[k];metadata.append(dict(Z1_Z2_powers=list(k),terms=len(g),multidegree_r_R_Z0=[max(e[i]for e in g)for i in range(3)],min_exponents=[min(e[i]for e in g)for i in range(3)]))
 prep=dict(variable_order=['r','R','Z0'],coefficient_group_order=metadata,first_two_intrinsically_Z0_free_and_nonproportional=True,
  field='QQ(R,Z0)[r]',Gauss_lemma='a primitive3variablegcd equals functionfield r-content up to units; first twoZ0-free groups force content independentZ0 up to fieldunits',
  global_Q_identity_requested='Q=content*G using all15 exact coefficient cofactors and literal6variable reconstruction',
  maximality_gate='native quotientcoefficientgcd1 plus later independent fieldEuclid/Bezout before exactrdegree statement',
  field_irreducibility_gate='accepted3/5 coprime s-map degrees give birational planeimage; primitiveRes irreducible overQQ(r,R,Z0)[Z1,Z2], then Gauss',
  generic_source_field='QQ(R,Z0,r,s)',spectral_base_field='QQ(R,Z0,Z1,Z2)',
  no_global_boundary_removal_properness_actual_NEW_target_index_or_equation_claim=True)
 (out/'method-preflight.json').write_text(json.dumps(prep,indent=2)+'\n')
 inp=work/'groups.txt'
 with inp.open('w')as f:
  f.write('3 15\n')
  for k in order:
   f.write(f'{k[0]} {k[1]} {len(groups[k])}\n')
   for m,c in sorted(groups[k].items(),reverse=True):f.write(str(c)+' '+' '.join(map(str,m))+'\n')
 with gzip.open(out/'coefficient-groups.txt.gz','wb')as f:f.write(inp.read_bytes())
 native=ROOT/'code/q_new_family_two_shape_content_gcd_20261104_v63.cpp';binary=work/'content-native';flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint','openssl'],text=True));error=None;nr=None
 try:
  subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True,timeout=max(.1,deadline-time.monotonic()-2))
  remaining=deadline-time.monotonic()-2;assert remaining>0
  with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(qpath),str(out)],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=remaining,preexec_fn=lambda:cap(max(1,math.ceil(remaining))))
  nr=read(out/'native-report.json');assert nr['all_passed']and nr['quotient_coefficient_gcd_is_one']
  assert all(s['gcd_API_return_status']==1 for s in nr['API_status_log'])
  assert nr['all15_exact_cofactor_identities_verified']and nr['source_grouping_all6919_terms_verified']and nr['full_Q_equals_content_times_primitive_verified']and nr['content_independent_of_Z0']
  for name,stat in[('content',nr['content']),('primitive',nr['primitive']),('quotient-coefficient-gcd',nr['quotient_coefficient_gcd'])]:
   digest=hashlib.sha256();count=0
   with(out/(name+'.txt')).open('rb')as f:
    header=f.readline().split();assert header[1]==str(stat['terms']).encode()
    for line in f:digest.update(line);count+=1
   assert count==stat['terms']and digest.hexdigest()==stat['canonical_SHA256']
 except Exception as exc:error=repr(exc)
 final=dict(id=STEM,all_passed=error is None,error=error,status='full_content_and_all_cofactors_native_verified_pending_independent_maximality'if error is None else 'bounded_failure',
  one_shared_deadline_seconds=60,address_space_bytes=2_000_000_000,coefficient_groups=15,
  content=nr['content']if nr else None,primitive=nr['primitive']if nr else None,
  all15_cofactors_and_exact_FULL_Q_identity_verified=bool(nr and nr['full_Q_equals_content_times_primitive_verified']),
  native_quotient_coefficient_gcd_is_one=bool(nr and nr['quotient_coefficient_gcd_is_one']),
  independent_maximality_gate_completed=False,exact_generic_source_r_field_degree_claimed=False,
  raw_determinant_identification_independently_accepted=True,actual_NEW_target_field_index_and_equation_still_open=True,
  no_general_factorization_resultant_PRS_point_prime_source26_or_automatic_retry=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),native,qpath,qreport,*gates]},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<CAP-10000
 (out/'report.json').write_text(json.dumps(final,indent=2)+'\n');print(json.dumps({k:final[k]for k in['all_passed','error','content','primitive','native_quotient_coefficient_gcd_is_one','seconds']},indent=2))
 if error:raise SystemExit(2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);args=p.parse_args();main(args.output.resolve(),args.work.resolve())
