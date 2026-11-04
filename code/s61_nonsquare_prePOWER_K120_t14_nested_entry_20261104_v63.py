#!/usr/bin/env python3
"""ONE nestedK120/t14 coefficient worker after the independent source gate."""
import argparse,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K120_T14_NESTED_ENTRY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={};cache={};dump(out/'phase.json',dict(phase='source/peer/partial-serialization gates beforeONEworker'))
 def source(rel):
  if rel not in cache:cache[rel]=read(ROOT/rel);pins[rel]=sha(ROOT/rel)
  return cache[rel]
 base='reports/S61_NONSQUARE_PREPOWER_K120_T14_NESTED_SOURCE_PREFLIGHT_20261104_v63';gate=source(base+'/report.json');plan=source(base+'/nested-plan.json');cost=source(base+'/cost.json');assert gate['resource_gate_passes']and gate['total_live_upper']==144877 and pins[base+'/report.json']=='10ff652bd664c2227399e7fb809917fb6e8b868d05294e8f5f29b712bcaff42c'
 peer=source(base+'_REVIEW/report.json');assert peer['all_passed']and pins[base+'_REVIEW/report.json']=='2c48439b98402610c26b5d9c22c3e34b5b06b36bfd27bacca901841c70eda2bd'
 assert source('reports/S61_NONSQUARE_PREPOWER_K120_WEDGE_SIGNS_20261104_v63/report.json')['all_passed']
 old=source(plan['external_scalar_nodes']['path']);assert pins[plan['external_scalar_nodes']['path']]==plan['external_scalar_nodes']['sha256'];N=old['nodes']+plan['appended_scalar_nodes'];assert len(N)==11523 and len(plan['outer21'])==21
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json')
 need=set();atoms=set(plan['all_atom_nodes'])
 for e in plan['outer21']:
  need.update(e['weighted2_source_nodes'])
  for t in e['inner10']:need.update(t['source_nodes']);assert t['left2minor']!=t['right3minor']
 assert atoms=={i for i in need if N[i]['op']=='atom'}
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if 'named_FULLU_unit'in r:return [[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
  if r.get('derived_existing_coefficient_difference'):
   hc=source('reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json');p={tuple(m):Q(c)for m,c in hc['raw_compact14x5'][11][1]}
   for m,c in hc['raw_compact14x5'][9][1]:p[tuple(m)]=p.get(tuple(m),Q(0))-Q(c)
   answer=[[list(m),str(c)]for m,c in sorted(p.items())if c];assert answer==r['polynomial'];return answer
  if r['path']=='accepted_Q2L_source_atom':
   a=r['definition'];return cp['matrix25x32'][a[1]][a[2]][0]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1][0]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]][0]
  x=source(r['path'])
  if 'old_DAG_root'in r:return next(v['Laurent_polynomial']for v in x['node_values']if v['old_DAG_root']==r['old_DAG_root'])
  if 'h9_target'in r:return x['source_h9'][x['original_target_labels'].index(r['h9_target'])]
  x=x[r['key']]
  for i in r.get('index',[]):x=x[i]
  return x
 ip=out/'native-input.txt';records=[];atomterms=0
 def spec(b,d):return ' '.join(map(str,[v for p in b['box']for v in p]+b['degree']+[d]))
 with ip.open('w')as f:
  f.write('11990 11523\n')
  for i,n in enumerate(N):
   if i not in need:f.write('x\n');continue
   pref=f"{n['bound']['support']if n['bound']else 0} {n['QQ_denominator']}"
   if n['op']=='atom':
    p=atom(n);atomterms+=len(p);m=min((e[0]for e,c in p),default=0);s=n['a_shift']+m;f.write(f"p {pref} {s} {n['QQ_unit']}\n{len(p)}\n")
    for e,c in p:f.write(f'{e[0]-m} {e[1]} {e[2]} {c}\n')
    records.append(dict(node=i,reference=n['reference'],QQ_unit=n['QQ_unit'],a_shift=s,terms=len(p)))
   elif n['op']=='known_a_shift':f.write(f"k {pref} {n['args'][0]} {n['a_exponent_shift']}\n")
   else:f.write(dict(add='a',sub='s',mul='m')[n['op']]+' '+pref+' '+' '.join(map(str,n['args']))+'\n')
  assert atomterms==11990;f.write(spec(plan['output_bound'],plan['global_output_integer_denominator'])+'\n21\n')
  for e in plan['outer21']:
   f.write(str(e['weighted2_root'])+'\n');f.write(str(len(e['weighted2_source_nodes']))+' '+' '.join(map(str,e['weighted2_source_nodes']))+'\n');f.write(spec(e['top5_bound'],e['top5_integer_denominator'])+'\n10\n')
   for t in e['inner10']:
    f.write(' '.join(map(str,[t['left2minor'],t['right3minor'],t['sign'],t['integer_QQ_scale']]))+'\n');f.write(str(len(t['source_nodes']))+' '+' '.join(map(str,t['source_nodes']))+'\n')
   f.write(str(e['outer_sign'])+' '+e['outer_integer_QQ_scale']+'\n')
 # Partial-inner serialization must bound all prefixes, not just its final sum.
 innerbytes=0;innerheights=[]
 for e in plan['outer21']:
  TD=int(e['top5_integer_denominator']);bound=0
  for t in e['inner10']:
   a=N[t['left2minor']];b=N[t['right3minor']];bound+=int(a['integer_l1_upper'])*int(b['integer_l1_upper'])*int(t['integer_QQ_scale'])
  ib=e['top5_dense_slots']*(len(str(bound))+45)+len(str(TD))+1000;innerbytes=max(innerbytes,ib);innerheights.append(dict(outer=e['index'],prefix_integer_l1_bits=bound.bit_length(),byte_bound=ib))
 prefix_bound=cost['whole_raw_artifact_byte_bound']+innerbytes+1000000;assert prefix_bound<32*1024*1024
 dump(out/'preflight.json',dict(total_owned_coefficient_slots=144877,complete_artifact_raw_bound=cost['whole_raw_artifact_byte_bound'],max_extra_partial_inner_bytes=innerbytes,whole_prefix_stop_artifact_bound=prefix_bound,inner_prefix_bounds=innerheights,no_hidden_global_value_DAG=True,clock_guard_seconds=50))
 dump(out/'source-binding.json',dict(source_plan_path=base+'/nested-plan.json',source_plan_sha256=pins[base+'/nested-plan.json'],independent_gate_report=base+'_REVIEW/report.json',independent_gate_sha256=pins[base+'_REVIEW/report.json'],literal_source_atoms=records,
  all_source_operations_pinned_by_phase_lists=True,source_current15='P4[0,1,3,6],G20current2,A=current11/c10,B=(current12-current10)/c11,C=current7/c6',original32_columns=[22,11,13,25,15,19,14],
  represented_integer_output='normalized(D7/g):raw(D7/g)=a²*output; actualnextentry=c11*output/N6',all3_multilinear_terms_retained=True,all_tail2_effects_retained=True,no_numeric_g_inverse_or_other_entry=True))
 dump(out/'phase.json',dict(phase='allpeer/source/prefixbytegatespassed beforeONE nativeworker'));cache.clear();del q,cp,pd,N,old,records,plan
 cpp=Path(__file__).with_suffix('.cpp');hpp=ROOT/'code/s61_nonsquare_prePOWER_K120_t14_nested_entry_main_20261104_v63.hpp';exe=out/'native-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)';run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,56-(time.monotonic()-began)));(out/'compile.log').write_text(run.stdout);assert run.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(ip)+' '+str(out);run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(run.stdout)
 np=out/'native.json'if(out/'native.json').exists()else out/'frontier.json';n=read(np)if np.exists()else dict(all_passed=False,returncode=run.returncode);ok=run.returncode==0 and n['all_passed'];dump(out/'phase.json',dict(phase='FINAL completeN7'if ok else'FINAL exactprefixonly, noautoretry'))
 for p in[Path(__file__),cpp,hpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 size=sum(p.stat().st_size for p in out.iterdir()if p.is_file()and p!=exe);assert size<32*1024*1024
 result=dict(id=ID,all_passed=ok,frozen=True,independent_peer_pending=True,normalized_N7_is_zero=n.get('normalized_N7_is_zero'),normalized_N7_terms=n.get('terms'),completed_inner=n.get('completed_inner210',n.get('completed_inner')),completed_outer=n.get('completed_outer21',n.get('completed_outer')),source_evaluations=n.get('source_operation_evaluations',n.get('source_evaluations')),actual_scalar_updates=n.get('actual_scalar_updates',n.get('scalar_updates')),max_owned_coefficient_slots=n.get('max_owned_live_coefficient_slots',n.get('max_slots')),whole_scalar_claim_ONLY_if_complete=ok,no_other_entry_or_rank_or_global57=True,artifact_bytes_excluding_binary=size,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-began)
 dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
