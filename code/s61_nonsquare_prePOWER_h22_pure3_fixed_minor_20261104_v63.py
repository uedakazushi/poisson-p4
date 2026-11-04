#!/usr/bin/env python3
"""ONE fixed actualpure3 minor, no alternative columns or factor search."""
import argparse,hashlib,importlib.util,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);assert not(out/'report.json').exists();phase=read(out/'phase.json');paths=[]
 def load(stem,file='certificate.json'):
  p=ROOT/'reports'/stem/file;paths.append(p);return read(p)
 height=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FRAGMENT_PREFLIGHT_20261104_v63');hr=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FRAGMENT_PREFLIGHT_20261104_v63','report.json');assert sha(paths[-1])=='a9459e904628708b2b1f01094c304792f3bcd665531e58f2cbc5d94a3830dd06'
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json')
 nodes=q['scalar_DAG']['nodes'];needed=height['needed_scalar_DAG_roots'];assert len(needed)==108;index={old:i for i,old in enumerate(needed)}
 helper=ROOT/'code/s61_nonsquare_prePOWER_h22_kernel9_support_preflight_20261104_v63.py';paths.append(helper);sp=importlib.util.spec_from_file_location('frozen_h22_bound',helper);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
 bound={};text=['1120 108'];literal=0;binding=[]
 def poly(p):
  nonlocal literal
  literal+=len(p);text.append(str(len(p)));text.extend(' '.join(map(str,ex))+' '+c for ex,c in p)
 for old in needed:
  n=nodes[old];op=n[0]
  if op=='z':b=m.Z;text.append(f'{old} z')
  elif op=='one':b=m.ONE;text.append(f'{old} 1')
  elif op in ['checkpoint','accepted_pair_second','new_unit_atom']:
   nd=cp['matrix25x32'][n[1]][n[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
   num,den=nd;assert len(den)==1 and den[0][0][1:]==[0,0];text.append(f'{old} p {den[0][0][0]} {den[0][1]}');poly(num);b=m.shift(m.atom(num),-den[0][0][0]);binding.append(dict(old_root=old,source_node=n,denominator=den,numerator_terms=len(num)))
  else:
   assert op in ['a','s','m'];text.append(f'{old} {op} {index[n[1]]} {index[n[2]]}');b=m.mm(bound[n[1]],bound[n[2]])if op=='m'else m.Z if op=='s'and n[1]==n[2]else m.aa(bound[n[1]],bound[n[2]])
  bound[old]=b
 cols=q['new22_target_original32_columns'][:17];rows=[1,3,6];chosen=[11,13,25]
 for i in rows:text.append(' '.join(str(index[q['remaining15x22_scalar_roots'][i][cols.index(j)]])for j in chosen+[22]))
 for label in chosen+[22]:poly(h['source_h9'][h['original_target_labels'].index(label)])
 inputcache=1120;cache=sum(m.size(b)for b in bound.values());forecast=dict(input_literal_reserve=inputcache,actual_literal_atom_polynomial_terms=literal,all108_scalar_node_value_support_upper=cache,
  nine_output_support_upper=4630,determinant_support_upper=10080,
  two_determinant_copies_plus_leibniz_workspace_upper=4*10080,cofactor_and_extra_multiply_workspace_upper=6*10080,
  total_live_upper=inputcache+cache+4630+10*10080,strict150k_guard=True)
 assert forecast['total_live_upper']<=150000
 (out/'preflight.json').write_text(json.dumps(forecast,indent=2)+'\n');(out/'source-binding.json').write_text(json.dumps(dict(atom_bindings=binding,selected_source15=rows,source_actual135=[18,39,68],fixed_columns=chosen,h22_alias='accepted actual14termQ9 coefficient',all8_Q9_tails_retained_in_parent_operator=True,source_global_chart_reference='S61_NONSQUARE_PREPOWER_HEAVY_Q_CHART_20261104_v63',no_Hadj_Bpole_or_other_minor=True),indent=2)+'\n')
 ip=out/'native-input.txt';ip.write_text('\n'.join(text)+'\n');del q,cp,pd,h,bound,text
 phase.update(phase='input_and_resource_gate_passed_before_single_native_job',preflight=forecast,height_report_SHA256=sha(paths[1]));(out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=out/'fixed3-native';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=max(1,53-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:r=subprocess.run([str(exe),str(ip),str(out/'native.json')],stdout=log,stderr=subprocess.STDOUT,timeout=max(1,56-(time.monotonic()-start)))
 n=read(out/'native.json')if(out/'native.json').exists()else dict(all_passed=False,returncode=r.returncode);good=r.returncode==0 and n['all_passed']
 if good:assert n['peak_live_terms_with_input_reserve']<=150000 and n['cofactor_and_six_permutation_identity_verified']
 paths += [Path(__file__),cpp,ROOT/f'research/notes/{ID}.md'];report=dict(id=ID,all_passed=good,frozen=True,independent_review_pending=True,
  fixed_actual3_minor_zero=n.get('determinant_is_zero'),determinant_terms=n.get('determinant_terms'),determinant_common_a_shift=n.get('determinant_common_a_shift'),
  peak_live_terms_with_input_reserve=n.get('peak_live_terms_with_input_reserve'),forecast_total_live_upper=forecast['total_live_upper'],
  no_alternative_minor_or_factor_search_point_prime_GB=True,no_full_kernel9_or_globalFULLU57_claim=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 try:main(a.output.resolve())
 except Exception as e:
  if a.output.exists():(a.output/'driver-failure.json').write_text(json.dumps(dict(all_passed=False,error=repr(e),no_auto_retry=True),indent=2)+'\n')
  raise
