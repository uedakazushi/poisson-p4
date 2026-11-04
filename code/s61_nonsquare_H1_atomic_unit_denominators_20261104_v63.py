#!/usr/bin/env python3
"""ONE exact228-denominator unit model plus graph-only native-cost estimate; no coefficient queries."""
import argparse,collections,hashlib,json,os,resource,signal,subprocess,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 base=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';ip=base/'certificate.json';assert sha(ip)==read(base/'report.json')['certificate_sha256'];x=read(ip);atoms=x['scalar_circuit']['atoms_num_den'];nodes=x['scalar_circuit']['nodes'];block=x['complete_H1_image_block_scalar_roots']
 unique=[];ids={};atom_to_den=[]
 for n,d in atoms:
  key=json.dumps(d,separators=(',',':'))
  if key not in ids:ids[key]=len(unique);unique.append(d)
  atom_to_den.append(ids[key])
 assert len(unique)==228 and len(atoms)==2378
 R,a,e,d=ring('a,e,d',s.QQ);known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
 enc=lambda p:[[list(m),str(c)]for m,c in p.items()]
 def line(p):return str(len(p))+''.join(' '+str(c)+' '+' '.join(map(str,m))for m,c in p)+'\n'
 native=out/'native-input.txt';cpp=ROOT/'code/s61_nonsquare_H1_atomic_unit_denominators_20261104_v63.cpp';exe=ROOT/'work/H1-atomic-unit-denominators-native'
 with native.open('w')as f:
  f.write('228 2378 6\n')
  for p in known:f.write(line(enc(p)))
  for p in unique:f.write(line(p))
  for i,(n,d)in enumerate(atoms):f.write(str(atom_to_den[i])+'\n');f.write(line(n))
 # Existing activated toolchain only; no installation or coefficient DAG execution.
 command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-began)))
 executable_sha=sha(exe)
 with(out/'native.log').open('w')as f:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','H1-native',str(exe),str(native),str(out/'native-certificate.json')],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-began)))
 nr=read(out/'native-certificate.json');assert nr['all_passed']and nr['denominator_reconstructions']==228 and nr['atomic_fraction_crossproducts']==2378
 # Graph-only cost estimates: exact reachability/counts, no scalar value evaluation.
 def profile(roots):
  need=set();todo=list(roots)
  while todo:
   i=todo.pop()
   if i in need:continue
   need.add(i);nd=nodes[i]
   if nd[0]in['a','s','m']:todo.extend(nd[1:])
  refs=collections.Counter(roots)
  for i in need:
   if nodes[i][0]in['a','s','m']:refs.update(nodes[i][1:])
  live=set();maxlive=0;depth={0:0}
  for i in sorted(need):
   live.add(i);nd=nodes[i];depth[i]=0 if nd[0]in['z','c']else 1+max(depth[j]for j in nd[1:]);maxlive=max(maxlive,len(live))
   if nd[0]in['a','s','m']:
    for dep in nd[1:]:
     refs[dep]-=1
     if not refs[dep]:live.discard(dep)
  ais={nodes[i][1]for i in need if nodes[i][0]=='c'}
  return{'reachable_nodes':len(need),'node_operations':dict(collections.Counter(nodes[i][0]for i in need)),'atomic_inputs':len(ais),'atomic_numerator_terms':sum(len(atoms[i][0])for i in ais),'atomic_denominator_terms':sum(len(atoms[i][1])for i in ais),'max_simultaneously_live_nodes':maxlive,'max_DAG_depth':max(depth.values())}
 rowprofiles=[[i,profile(block[i])]for i in range(113,134)];whole=profile([r for row in block[113:]for r in row])
 unit_model={'source_input_sha256':sha(ip),'known_factor_order':['a','d','F','B','C','J'],'known_factor_polynomials':[enc(p)for p in known],
  'atomic_denominator_id_by_atom':atom_to_den,'unique_atomic_denominators':unique,'literal_constant_and_powers_by_denominator':nr['factors'],
  'finite_unit_model_all2378_source_atoms':True,'native_reconstructs228_and_checks2378_literal_crossproducts':True,
  'quartic_rows113_to133_graph_profiles':rowprofiles,'first_psi113_graph_profile':rowprofiles[0][1],'all21_quartics_shared_graph_profile':whole,
  'graph_only_no_quartic_coefficient_value_computed':True,
  'exact_native_arithmetic_recipe':'Represent scalar as QQ polynomial numerator / product known_factors^denvector. Add by componentwise maxden(lcm) alignment; multiply numerator polynomials and adddenvectors; strip ONLY named knownfactors by exact divides to control growth. Generic polynomialGCD is unnecessary for correctness.',
  'cost_limit':'Reachability/liveness counts are exact; runtime and term growth cannot be inferred without a separately budgeted native coefficient pilot. Guards remain150k intermediate/100k output terms perentry/32MiB artifacts.',
  'old113_rows_and7_pivots_unchanged':True,'no_Hquartic_retry_Gaussian_DAGvalues_point_prime_GB_or_new_kernel_claim':True}
 (out/'certificate.json').write_text(json.dumps(unit_model,indent=2)+'\n');exe.unlink()
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),ip,base/'report.json',ROOT/'reports/S61_NONSQUARE_H1_NATIVE_ARITHMETIC_PREFLIGHT_20261104_v63/report.json',ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz']
 result={'id':ID,'all228_denominators_exactly_known_FULLU_units':True,'all2378_source_atom_fraction_crossproducts_checked':True,'exact_named_divisions':nr['exact_named_divisions'],'first_psi113_graph_profile':rowprofiles[0][1],'all21_shared_graph_profile':whole,
  'quartic_DAG_value_queries':0,'Gaussian_jobs':0,'native_executable_sha256_not_committed':executable_sha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
