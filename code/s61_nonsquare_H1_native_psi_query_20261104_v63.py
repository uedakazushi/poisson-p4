#!/usr/bin/env python3
"""ONE native factored-unit psi113 query; no other quartic, old-row reevaluation or Gaussian."""
import argparse,gzip,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_H1_NATIVE_PSI_QUERY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 ub=ROOT/'reports/S61_NONSQUARE_H1_ATOMIC_UNIT_DENOMINATORS_20261104_v63';ur=read(ub/'report.json');up=ub/'certificate.json';assert sha(up)==ur['evidence_sha256'][up.name];units=read(up);assert units['finite_unit_model_all2378_source_atoms']
 base=ROOT/'reports/S61_NONSQUARE_H_SQUARE_H1_PREFLIGHT_20261104_v63';ip=base/'certificate.json';assert sha(ip)==units['source_input_sha256'];x=read(ip);nodes=x['scalar_circuit']['nodes'];atoms=x['scalar_circuit']['atoms_num_den'];roots=x['complete_H1_image_block_scalar_roots'][113];assert len(roots)==23
 cp=ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_CONTINUATION_20261104_v63/checkpoint.json.gz';cr=read(cp.parent/'report.json');assert sha(cp)==cr['evidence_sha256'][cp.name];cache=json.loads(gzip.decompress(cp.read_bytes()));assert len(cache['materialized_raw_rows'])==113 and cache['completed_pivots']==7
 rv=ROOT/'reports/S61_NONSQUARE_H1_KNOWN_UNIT_ELIMINATION_20261104_v63_REVIEW/report.json';assert read(rv)['all_passed']
 needed=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if nodes[i][0]in['a','s','m']:todo.extend(nodes[i][1:])
 aids=sorted({nodes[i][1]for i in needed if nodes[i][0]=='c'});assert len(needed)==31860 and len(aids)==631
 native=out/'native-input.txt';cpp=ROOT/'code/s61_nonsquare_H1_native_psi_query_20261104_v63.cpp';exe=ROOT/'work/H1-native-psi-query'
 def polyline(p):return str(len(p))+''.join(' '+c+' '+' '.join(map(str,m))for m,c in p)+'\n'
 with native.open('w')as f:
  f.write(f'6 {len(aids)} {len(needed)} 23 7 {len(nodes)}\n')
  for p in units['known_factor_polynomials']:f.write(polyline(p))
  for ai in aids:
   proof=units['literal_constant_and_powers_by_denominator'][units['atomic_denominator_id_by_atom'][ai]];f.write(str(ai)+' '+proof['constant']+' '+' '.join(map(str,proof['a_d_F_B_C_J_powers']))+'\n');f.write(polyline(atoms[ai][0]))
  for i in sorted(needed):f.write(str(i)+' '+' '.join(map(str,nodes[i]))+'\n')
  f.write(' '.join(map(str,roots))+'\n')
  for st in cache['unit_pivots']:
   f.write(str(st['pivot_column'])+'\n')
   for n,d in st['normalized23']:f.write(polyline(n));f.write(polyline(d))
 command='source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 with(out/'compile.log').open('w')as f:subprocess.run(['bash','-lc',command],stdout=f,stderr=subprocess.STDOUT,check=True,timeout=max(1,55-(time.monotonic()-began)))
 executable_sha=sha(exe)
 with(out/'native.log').open('w')as f:
  run=subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','H1-psi-native',str(exe),str(native),str(out)],stdout=f,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-began)))
 if(out/'native-audit.json').exists():audit=read(out/'native-audit.json')
 else:audit={'all_passed':False,'stop_reason':'native_abnormal_exit_without_complete_audit','native_returncode':run.returncode}
 complete=run.returncode==0 and audit['all_passed'];outputs={}
 def parse_factored(p):
  ls=p.read_text().splitlines();count=int(ls[0]);records=[]
  for line in ls[1:]:
   t=line.split();den=list(map(int,t[:6]));n=int(t[6]);assert len(t)==7+4*n;records.append({'signed_a_d_F_B_C_J_denominator_powers':den,'numerator':[[list(map(int,t[k+1:k+4])),t[k]]for k in range(7,len(t),4)]})
  assert len(records)==count;return records
 if complete:
  assert audit['evaluated_DAG_nodes']==31860 and audit['old7_input_fraction_crossproducts']==161 and audit['all23_final_schur_identities_checked']
  for name,count in[('raw23-factored.txt',23),('schur23-factored.txt',23),('schur7-coefficients-factored.txt',7)]:
   outputs[name]=parse_factored(out/name);assert len(outputs[name])==count
  for st in cache['unit_pivots']:assert not outputs['schur23-factored.txt'][st['pivot_column']]['numerator']
 certificate={'coefficient_domain':units['known_factor_order'],'source_H1_raw134_row':113,'source_row_label':x['all_remaining21_quartic_H2free_labels'][0],
  'literal_actual_full_source':'psi=-b131*Qhat77+b77*Qhat131, in the reviewed joint frame and NF117; after H2 deletion its H1 projection is scalar row113.',
  'completed_all23_raw_and_old7_schur_coefficients':complete,'factored_scalar_representation':'n / (a^v0*d^v1*F^v2*B^v3*C^v4*J^v5), signed integers permitted',
  'arithmetic_input_scalar_roots23':roots,'needed_node_original_ids':sorted(needed),'needed_atom_original_ids':aids,'native_unit_arithmetic_audit':audit,
  'all23_H1_coefficients_and_schur7_lift_if_complete':outputs,'old7_pivot_source_rows':[st['raw134_source_row']for st in cache['unit_pivots']],
  'full_source_schur_recipe':'psi_reduced = psi - sum_i c_i * previous_normalized_full_pivot_i, c_i from saved schur7 file; the old normalized full source pivots have their exact saved recursive lifts. This retains the actual pure98 coefficients lazily and applies the same target change to all134 sources.',
  'old113_source_DAG_reevaluations':0,'other20_quartic_queries':0,'Gaussian_jobs':0,'pure98_queries':0,'known38_missing42_and7_H1_units_unchanged':True,
  'no_full121_full239_or_new_unit_pivot_claim':True}
 (out/'certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n');exe.unlink(missing_ok=True)
 paths=[Path(__file__),cpp,ROOT/'research/notes'/(ID+'.md'),ub/'report.json',up,ub/'native-certificate.json',base/'report.json',ip,cp,rv]
 result={'id':ID,'completed_psi113_query_only':complete,'native_audit':audit,'old113_source_DAG_reevaluations':0,'other20_quartic_queries':0,'Gaussian_jobs':0,'pure98_queries':0,
  'existing7_units_and_global125_lower_unchanged':True,'native_executable_sha256_not_committed':executable_sha,'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in paths},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-began,'frozen':True,'independent_review_pending':True}
 assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())<32*1024*1024;signal.alarm(0)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
