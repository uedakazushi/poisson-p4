#!/usr/bin/env python3
"""ONE stronger rawcross18-only sourcecheck; no kappa18 evaluation."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={};cache={};dump(out/'phase.json',dict(phase='sourcegates beforeONE3rawcross18onlyworker'))
 def source(rel):
  if rel not in cache:cache[rel]=read(ROOT/rel);pins[rel]=sha(ROOT/rel)
  return cache[rel]
 base='reports/S61_NONSQUARE_PREPOWER_K92_T18_SOURCE_PREFLIGHT_20261104_v63';gate=source(base+'/report.json');op=source(base+'/operator-and-cost.json');profiles=source(base+'/profiles.json');p=profiles[1];assert pins[base+'/report.json']=='875b0fc35f2ab3cd5da262d3b689295e1e850cd1d5923ff9be15706f280ca659'and p['resource_gate_passes']and p['total_live_upper']==37236
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json');N=op['nodes'];need=set(p['needed_nodes']);assert len(N)==76 and len(need)==62 and op['common_B_clearing_power']==0
 assert not(set(op['kappa4_roots'])&need);assert op['kappa18_output']not in need
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if 'named_FULLU_unit'in r:return [[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
  if r['path']=='accepted_Q2L_source_atom':
   a=r['definition'];return cp['matrix25x32'][a[1]][a[2]][0]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1][0]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]][0]
  x=source(r['path'])[r['key']]
  for i in r.get('index',[]):x=x[i]
  return x
 cross=op['raw_cross3_outputs'];anchors=[N[N[cross[0]]['args'][0]]['args'][0]]+[N[N[i]['args'][1]]['args'][0]for i in cross];columns=op['column18_roots'];ip=out/'native-input.txt';atoms=[]
 with ip.open('w')as f:
  f.write('3564 76\n')
  for i,n in enumerate(N):
   if i not in need:f.write('x\n');continue
   pref=f"{n['bound']['support']if n['bound']else 0} {n['QQ_denominator']}"
   if n['op']=='atom':
    x=atom(n);m=min((e[0]for e,c in x),default=0);s=n['a_shift']+m;f.write(f"p {pref} {s} {n['QQ_unit']}\n{len(x)}\n")
    for e,c in x:f.write(f'{e[0]-m} {e[1]} {e[2]} {c}\n')
    atoms.append(dict(node=i,reference=n['reference'],QQ_unit=n['QQ_unit'],a_shift=s,terms=len(x)))
   else:f.write(dict(add='a',sub='s',mul='m')[n['op']]+' '+pref+' '+' '.join(map(str,n['args']))+'\n')
  f.write(' '.join(map(str,cross+columns+anchors))+'\n')
 dump(out/'source-binding.json',dict(current15_rows=[2,4,5,8],original32_target=18,anchor_original32=22,column18_roots=columns,anchor_roots=anchors,rawcross3_outputs=cross,source_ops62=sorted(need),literal_atoms=atoms,B_clear_power=0,NO_kappa18_nodes_evaluated=True,previous_raw7_scope=[22,11,13,25,15,19,14],all8tail_effects_retained=True))
 dump(out/'preflight.json',dict(whole_job_live_bound=37236,post_readback_bound=31012,raw_artifact_bound_with_JSON_ASCII=2*p['raw_artifact_bytes_upper']+1000000,source_gate_sha256=pins[base+'/report.json'],no_det11_or_othercolumn=True));dump(out/'phase.json',dict(phase='source/readback gatespassed beforeONE worker'));cache.clear();del q,cp,pd,N,atoms,op
 cpp=Path(__file__).with_suffix('.cpp');hpp=ROOT/'code/s61_nonsquare_prePOWER_rawcross18_entry_main_20261104_v63.hpp';exe=out/'native-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)';run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,56-(time.monotonic()-began)));(out/'compile.log').write_text(run.stdout);assert run.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(ip)+' '+str(out);run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(run.stdout)
 nf=out/'native.json'if(out/'native.json').exists()else out/'frontier.json';n=read(nf)if nf.exists()else dict(all_passed=False,returncode=run.returncode);ok=run.returncode==0 and n['all_passed'];rz=ok and all(n['raw_cross3_zero']);dump(out/'phase.json',dict(phase='FINAL sourcecross18complete'if ok else'FINAL stopped; noautoretry'))
 for x in[Path(__file__),cpp,hpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']:pins[str(x.relative_to(ROOT))]=sha(x)
 size=sum(x.stat().st_size for x in out.iterdir()if x.is_file()and x!=exe);assert size<32*1024*1024
 report=dict(id=ID,all_passed=ok,frozen=True,independent_peer_pending=True,raw_cross3_zero=n.get('raw_cross3_zero'),output_terms=n.get('output_terms'),rawpair_fixed8_rank_le1=rz,rawpair_fixed8_original32=[22,11,13,25,15,19,14,18]if rz else None,NO_kappa18_evaluation=True,independent_readback_reproducts=n.get('independent3_rawcross_readback_reproducts_zero'),max_total_live_slots=n.get('max_total_live_coefficient_slots',n.get('max_slots')),no_det11_or_othercolumn_or_rankclaim=True,artifact_bytes_excluding_binary=size,input_sha256=pins,evidence_sha256={x.name:sha(x)for x in out.iterdir()if x.is_file()and x!=exe},seconds=time.monotonic()-began)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
