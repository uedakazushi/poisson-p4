#!/usr/bin/env python3
"""ONE cleared Schur numerator u4/target15, under the refined frozen resource gate."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_ENTRY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);assert not(out/'report.json').exists();phase=read(out/'phase.json');paths=[];cache={}
 def source(rel):
  p=ROOT/rel
  if p not in cache:cache[p]=read(p);paths.append(p)
  return cache[p]
 oprel='reports/S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63/operator-and-cost.json';op=source(oprel);N=op['nodes'];assert len(N)==768
 gatep=ROOT/'reports/S61_NONSQUARE_PREPOWER_U4_T15_DEGREE_SUPPORT_PREFLIGHT_20261104_v63/report.json';gate=read(gatep);paths.append(gatep);assert sha(gatep)=='e1aa471a137491230806be414f11a6cb89f0c313472d6563a6598669c2d20f07'and gate['resource_gate_passes']
 dc=source('reports/S61_NONSQUARE_PREPOWER_U4_T15_DEGREE_SUPPORT_PREFLIGHT_20261104_v63/certificate.json');direct=source('reports/S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63/operator.json');assert sha(ROOT/'reports/S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63/operator.json')=='4bbcd0b3aa6dcb70c991ed75ce0f213cd8253a3e060736ce63b363d1a3123397'
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json')
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if r['path']=='accepted_Q2L_source_atom':
   x=r['definition'];nd=cp['matrix25x32'][x[1]][x[2]]if x[0]=='checkpoint'else pd['Lframe_dual6x2'][x[1]][1]if x[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][x[1]]
   assert len(nd[1])==1 and nd[1][0][0][1:]==[0,0]
   return nd[0]
  x=source(r['path']);k=r['key']
  if k=='node_values':return next(y['Laurent_polynomial']for y in x[k]if y['old_DAG_root']==r['old_DAG_root'])
  if k=='source_h9':return x[k][x['original_target_labels'].index(r['target'])]
  x=x[k]
  for i in r.get('index',[]):x=x[i]
  return x
 # This parser holds all literal atom polynomials until their use. Charge
 # one EXTRA literal-size reserve compared with the abstract last-use proof.
 actual_bound=gate['new_total_live_upper']+op['literal_atom_terms'];assert actual_bound==149178 and actual_bound<=150000
 preflight=dict(refined_DAG_live_upper=gate['new_total_live_upper'],additional_preloaded_atom_upper=op['literal_atom_terms'],implementation_live_upper=actual_bound,
  literal_reserve=op['literal_atom_terms'],generic_reserve=15000,output_scratch_included=True,
  no_shifted_input_polynomial_copies=True,output_support_upper=gate['output_support_upper'],coefficient_bits=346,QQ_denominator_bits=141,raw_artifact_bytes_upper=gate['raw_artifact_bytes_upper']+1048576)
 (out/'preflight.json').write_text(json.dumps(preflight,indent=2)+'\n');lines=[f"{op['literal_atom_terms']} 15000 768"];atoms_manifest=[]
 for i,n in enumerate(N):
  assert n['bound'].get('known_B_denominator_power',0)==0
  if n['op']=='atom':
   p=atom(n);m=min((ex[0]for ex,c in p),default=0);shift=n['a_shift']+m;assert n['B_denominator_power']==0
   lines.append(f"p {shift} {n['QQ_denominator_unit']}");lines.append(str(len(p)))
   for ex,c in p:lines.append(f'{ex[0]-m} {ex[1]} {ex[2]} {c}')
   atoms_manifest.append(dict(node=i,reference=n['reference'],normalization=dict(a_shift=shift,QQ_unit=n['QQ_denominator_unit']),terms=len(p)))
  else:assert n['op']in['add','sub','mul'];lines.append(f"{dict(add='a',sub='s',mul='m')[n['op']]} "+' '.join(map(str,n['args'])))
 specials=op['source_coefficients3'];lines.append(f"{op['output']} {len(specials)} "+' '.join(map(str,specials)))
 ip=out/'native-input.txt';ip.write_text('\n'.join(lines)+'\n')
 identity=dict(Mark_hybrid_export=direct['exports']['residual6x13'],Mark_entry=[0,9],current15_source5_and_original32target15=True,
  scalar_relation=op['source_to_Mark_hybrid_relation'],rho_LCM_powers=op['known_rho_LCM_powers'],
  all768_operation_identities_bound_to_source=True,literal_atoms=atoms_manifest,
  source_Hpair_identity='N[u,remaining3]=0; (Nupair*Adj(A2))*A2=delta2*Nupair. Thus lower3 cofactorneverenters thisrow.',
  source_P4_identity='Adj(P4)*P4=P4*Adj(P4)=D4 I4; multiplyingout gives the direct4pivotSchur and the stated HYBRIDscalar relation.',
  all13_Q9_remaining_tails_kept_by_parent_operator=True,
  no_other_entry_or_minor=True)
 (out/'source-binding.json').write_text(json.dumps(identity,indent=2)+'\n');phase.update(phase='all_source_atom_and_149178_resource_gates_passed_before_worker');(out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
 # The parent process releases parsed source coefficient caches before native.
 cache.clear();del q,cp,pd,N,lines
 cpp=Path(__file__).with_suffix('.cpp');exe=out/'single-entry-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out)],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-start)))
 n=read(out/'native.json')if(out/'native.json').exists()else read(out/'frontier.json')if(out/'frontier.json').exists()else dict(all_passed=False,returncode=r.returncode);ok=r.returncode==0 and n['all_passed']
 if ok:assert n['peak_live_terms_with_input_work_and_output_scratch']<=150000 and n['all_inplace_exponent_denominator_identities_checked']
 paths += [Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md'];result=dict(id=ID,all_passed=ok,frozen=True,independent_review_pending=True,
  numerator_is_zero=n.get('numerator_is_zero'),numerator_terms=n.get('numerator_terms'),known_a_clearing_power=n.get('known_a_clearing_power'),
  peak_live_terms_with_all_reserves=n.get('peak_live_terms_with_input_work_and_output_scratch'),implementation_forecast149178=True,
  source_only_HYBRID_scalar_comparison_saved=True,no_full6_rank_or_global57_or_chart_complement_removal=True,
  no_other_entry_minor_factorsearch_point_prime=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
