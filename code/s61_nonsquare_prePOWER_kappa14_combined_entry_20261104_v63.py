#!/usr/bin/env python3
"""ONE combined kappa14/rawcross3 job under the frozen75k source gate."""
import argparse,hashlib,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={};cache={};dump(out/'phase.json',dict(phase='beforeONE combinedkappa14/rawcross3 source worker'))
 def source(rel):
  if rel not in cache:cache[rel]=read(ROOT/rel);pins[rel]=sha(ROOT/rel)
  return cache[rel]
 base='reports/S61_NONSQUARE_PREPOWER_KAPPA14_PREFLIGHT_20261104_v63';gate=source(base+'/report.json');op=source(base+'/operator-and-cost.json');profiles=source(base+'/profiles.json');p=profiles[2];assert pins[base+'/report.json']=='b9de7b3ac05ae9ecb7585ce3630576e8d00f6ed7134b047106433aff63f65075'and p['resource_gate_passes']and p['total_live_upper']==75284
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json');N=op['nodes'];need=set(p['needed_nodes']);assert len(N)==76 and len(need)==73 and op['common_B_clearing_power']==0
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if 'named_FULLU_unit'in r:return [[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
  if r['path']=='accepted_Q2L_source_atom':
   a=r['definition'];return cp['matrix25x32'][a[1]][a[2]][0]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1][0]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]][0]
  x=source(r['path'])[r['key']]
  for i in r.get('index',[]):x=x[i]
  return x
 cross=op['raw_cross3_outputs'];anchors=[N[N[cross[0]]['args'][0]]['args'][0]]+[N[N[i]['args'][1]]['args'][0]for i in cross];assert len(set(anchors))==4;outputs=[op['kappa14_output']]+cross;cofs=op['kappa4_roots'];columns=op['column14_roots'];ip=out/'native-input.txt';atoms=[]
 with ip.open('w')as f:
  f.write('8484 76\n')
  for i,n in enumerate(N):
   if i not in need:f.write('x\n');continue
   pref=f"{n['bound']['support']if n['bound']else 0} {n['QQ_denominator']}"
   if n['op']=='atom':
    x=atom(n);m=min((e[0]for e,c in x),default=0);shift=n['a_shift']+m;f.write(f"p {pref} {shift} {n['QQ_unit']}\n{len(x)}\n")
    for e,c in x:f.write(f'{e[0]-m} {e[1]} {e[2]} {c}\n')
    atoms.append(dict(node=i,reference=n['reference'],QQ_unit=n['QQ_unit'],a_shift=shift,terms=len(x)))
   else:f.write(dict(add='a',sub='s',mul='m')[n['op']]+' '+pref+' '+' '.join(map(str,n['args']))+'\n')
  f.write(' '.join(map(str,outputs+cofs+columns+anchors))+'\n')
 dump(out/'source-binding.json',dict(current15_rows=op['current15_source_rows'],target_original32=14,anchor_original32=22,B_clear_power=0,kappa_coefficients_saved_source_peer='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63_REVIEW',all73_source_ops=sorted(need),literal_atoms=atoms,kappa14_output=outputs[0],rawcross3_outputs=cross,cofactor_roots=cofs,column14_roots=columns,anchor_roots=anchors,
  exact_rawcross_formulas='M[2,22]*M[i,14]-M[i,22]*M[2,14],i=4,5,8',all8tail_effects_retained=True,full57_lift_accepted=True,zero_kappa_NOT_equal_raw_proportionality=True))
 dump(out/'preflight.json',dict(core_total_live_bound=75284,post_reverse_reproduct_bound=69504,full_job_live_bound=75284,raw_artifact_bound_including_double_JSON_ASCII=2*p['raw_artifact_bytes_upper']+1000000,source_gate_report_sha256=pins[base+'/report.json'],no_other_candidate=True))
 dump(out/'phase.json',dict(phase='allsource/scale/provenance gatespassed beforeONE nativeworker'));cache.clear();del q,cp,pd,N,atoms,op
 cpp=Path(__file__).with_suffix('.cpp');exe=out/'native-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)';run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,56-(time.monotonic()-began)));(out/'compile.log').write_text(run.stdout);assert run.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(ip)+' '+str(out);run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-began)));(out/'native.log').write_text(run.stdout)
 np=out/'native.json'if(out/'native.json').exists()else out/'frontier.json';n=read(np)if np.exists()else dict(all_passed=False,returncode=run.returncode);ok=run.returncode==0 and n['all_passed'];dump(out/'phase.json',dict(phase='FINAL nativecompleted'if ok else'FINAL guardedstop; noautoretry'))
 rawzero=ok and all(n['raw_cross3_zero']);kzero=ok and n['kappa14_zero'];scope=dict(no_raw_proportionality_from_kappa_alone=True,raw_proportionality_on_seven_columns_proved=rawzero,
  kappa_quotient_span_proved=kzero,raw_proportionality_obstructed=ok and not rawzero,all_other_columns_retained=True,
  consequence='Threeexactrawcrosszeros plus nonzero savedG20anchor provefullpair4x7rank<=1 oncolumns22,11,13,25,15,19,14; uncorrectedliteralpair representationvalidwithchart.'if rawzero else'Onlykappazero: requireexplicitthree-columnspan corrections onALL15rows beforebare7 inference.'if kzero else'Exactobstructiontonextcolumn recorded; donotextrapolateanyothercolumn.')
 dump(out/'scope.json',scope)
 for f in[Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']:pins[str(f.relative_to(ROOT))]=sha(f)
 size=sum(f.stat().st_size for f in out.iterdir()if f.is_file()and f!=exe);assert size<32*1024*1024
 report=dict(id=ID,all_passed=ok,frozen=True,independent_peer_pending=True,kappa14_zero=n.get('kappa14_zero'),raw_cross3_zero=n.get('raw_cross3_zero'),output_terms=n.get('output_term_counts'),raw_proportionality_on_seven_columns_proved=rawzero,logical_kappa_vs_rawdistinction_preserved=True,all_source_columns_cofactors_and_anchor_values_saved=n.get('saved_columns_cofactors_anchors_and_all4_outputs'),reverse_reproducts_passed=ok,max_total_live_slots=n.get('max_total_live_coefficient_slots',n.get('max_slots')),no_further_determinant_column_family_or_rankclaim=True,artifact_bytes_excluding_binary=size,input_sha256=pins,evidence_sha256={f.name:sha(f)for f in out.iterdir()if f.is_file()and f!=exe},seconds=time.monotonic()-began)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
