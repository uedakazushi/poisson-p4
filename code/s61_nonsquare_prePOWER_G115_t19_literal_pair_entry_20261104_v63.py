#!/usr/bin/env python3
"""ONE actual G115/target19 normalized bordered6 under its frozen source gate."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_ENTRY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={};cache={};dump(out/'phase.json',dict(phase='source chart checks before ONE scalar worker'))
 def source(rel):
  if rel not in cache:cache[rel]=read(ROOT/rel);pins[rel]=sha(ROOT/rel)
  return cache[rel]
 base='reports/S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_PREFLIGHT_20261104_v63';gate=source(base+'/report.json');plan=source(base+'/operator-and-cost.json');assert pins[base+'/report.json']=='eb7495b6e5692bb84c9ca2b9df40922fc0a3132146a065e2651e0118776f98ca'and gate['resource_gate_passes']and gate['total_live_upper']==108513
 bot=source('reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json');hc=source('reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json');anchor=bot['bottom4x6'][0][0];g=hc['raw_compact14x5'][6][1];assert len(anchor)==156 and len(g)==26 and all(c!='0'for m,c in anchor+g)
 assert source('reports/S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW/report.json')['all_passed'];assert source('reports/S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63_REVIEW/report.json')['all_passed']
 dump(out/'nonempty-chart-proof.json',dict(literal_G20_anchor22=anchor,nonzero_anchor_terms=156,literal_g_N6_1=g,nonzero_g_terms=26,D5_G20_nonzero='Pair pure6 rank<=1 with nonzero G20 anchor; existing nonzero Saux is a linear combination of pair5minors, so all nonzero pair rows includingG20 have nonzero5minor in fractionfield.',no_D5_expansion=True,additional_chart='D(g*D5_G20)',complement_retained=True))
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json');N=plan['nodes'];need=set(plan['needed_fused_nodes']);assert len(N)==1357 and len(need)==1325 and plan['row_a_factors']==[0,1,1,0,0,0]and plan['row_B_clearings']==[0]*6
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if 'named_FULLU_unit'in r:return [[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
  if r['path']=='accepted_Q2L_source_atom':
   a=r['definition'];return cp['matrix25x32'][a[1]][a[2]][0]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1][0]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]][0]
  x=source(r['path'])
  if 'old_DAG_root'in r:return next(v['Laurent_polynomial']for v in x['node_values']if v['old_DAG_root']==r['old_DAG_root'])
  if 'h9_target'in r:return x['source_h9'][x['original_target_labels'].index(r['h9_target'])]
  x=x[r['key']]
  for i in r.get('index',[]):x=x[i]
  return x
 ip=out/'native-input.txt';atoms=[]
 with ip.open('w')as f:
  f.write('7174 1357\n')
  for i,n in enumerate(N):
   if i not in need:f.write('x\n');continue
   pref=f"{n['bound']['support']if n['bound']else 0} {n['QQ_denominator']}"
   if n['op']=='atom':
    p=atom(n);m=min((e[0]for e,c in p),default=0);shift=n['a_shift']+m;f.write(f"p {pref} {shift} {n['QQ_unit']}\n{len(p)}\n")
    for e,c in p:f.write(f'{e[0]-m} {e[1]} {e[2]} {c}\n')
    atoms.append(dict(node=i,reference=n['reference'],QQ_unit=n['QQ_unit'],a_shift=shift,terms=len(p)))
   elif n['op']=='known_a_shift':f.write(f"k {pref} {n['args'][0]} {n['a_exponent_shift']}\n")
   else:f.write(dict(add='a',sub='s',mul='m')[n['op']]+' '+pref+' '+' '.join(map(str,n['args']))+'\n')
  ob=N[plan['output']]['bound'];D=int(N[plan['output']]['QQ_denominator']);f.write(' '.join(map(str,[v for pp in ob['box']for v in pp]+ob['degree']+[D]))+'\n15\n');events=[]
  for k,e in enumerate(plan['laplace15']):
   a,b=e['top4_root'],e['bottom2_root'];den=int(N[a]['QQ_denominator'])*int(N[b]['QQ_denominator']);assert D%den==0;events.append([max(a,b),k,a,b,e['sign'],D//den])
  for e in sorted(events):f.write(' '.join(map(str,e))+'\n')
 dump(out/'source-binding.json',dict(current15_source_rows=[0,1,3,6,2,11,7],source_combinations=plan['source_rows'],original32_columns=plan['original32_columns'],literal_atoms=atoms,all1325_operations=sorted(need),all15_laplace_events=events,row_a_factors=plan['row_a_factors'],raw_D6_equals='a^2*normalized_D6',row_B_clearings=plan['row_B_clearings'],residual='c10*D6/(g*D5_G20)',no_D5_expansion=True,extra_chart_complement_retained=True,full57_source_lift_peer_bound=True,all7_remaining_tail_effects_retained=True))
 dump(out/'preflight.json',dict(total_whole_worker_live_upper=108513,output_lattice_slots=38931,raw_artifact_bytes_upper=gate['raw_artifact_bytes_upper'],hard_native_seconds=50,no_other_candidate=True))
 dump(out/'phase.json',dict(phase='all source/chart/denominator checks passed before ONE native worker'));cache.clear();del q,cp,pd,N,plan,atoms,bot,hc,anchor,g
 cpp=Path(__file__).with_suffix('.cpp');exe=out/'native-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)';run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,56-(time.monotonic()-began)));(out/'compile.log').write_text(run.stdout);assert run.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(ip)+' '+str(out);run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(run.stdout)
 nf=out/'native.json'if(out/'native.json').exists()else out/'frontier.json';n=read(nf)if nf.exists()else dict(all_passed=False,returncode=run.returncode);ok=run.returncode==0 and n['all_passed'];dump(out/'phase.json',dict(phase='FINAL complete normalized D6'if ok else'FINAL stopped exact prefix only; no retry'))
 for p in[Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 size=sum(p.stat().st_size for p in out.iterdir()if p.is_file()and p!=exe);assert size<32*1024*1024
 report=dict(id=ID,all_passed=ok,frozen=True,independent_peer_pending=True,normalized_D6_is_zero=n.get('normalized_D6_is_zero'),normalized_D6_terms=n.get('terms'),max_total_live_coefficient_slots=n.get('max_total_live_coefficient_slots',n.get('max_slots')),actual_scalar_updates=n.get('actual_scalar_updates',n.get('scalar_updates')),all15_completed=n.get('completed_laplace15',n.get('completed_laplace')),raw_D6_equals_a_squared_times_normalized=True,whole_D6_claim_ONLY_if_complete=ok,additional_chart_g_D5_retained=True,no_D5_expansion_or_other_entry_or_rank_or_global57=True,artifact_bytes_excluding_binary=size,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-began)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
