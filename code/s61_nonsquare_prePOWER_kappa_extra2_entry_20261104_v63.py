#!/usr/bin/env python3
"""ONE B*kappa*M(extra-tail2) exact coefficient worker; projected-source only."""
import argparse,hashlib,itertools,json,resource,signal,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def universal_cofactor_identity():
 def mon(v):return tuple(sorted(v))
 def add(a,b,s=1):
  r=a.copy()
  for m,c in b.items():
   r[m]=r.get(m,0)+s*c
   if not r[m]:del r[m]
  return r
 def mul(a,b):
  r={}
  for m,c in a.items():
   for n,d in b.items():
    k=mon(m+n);r[k]=r.get(k,0)+c*d
    if not r[k]:del r[k]
  return r
 A=[[{(3*i+j,):1}for j in range(3)]for i in range(4)];co=[]
 for i in range(4):
  rows=[r for r in range(4)if r!=i];q={}
  for p in itertools.permutations(range(3)):
   t={():1}
   for k in range(3):t=mul(t,A[rows[k]][p[k]])
   parity=i+sum(p[k]>p[l]for k in range(3)for l in range(k+1,3));q=add(q,t,-1 if parity%2 else 1)
  co.append(q)
 for j in range(3):
  r={}
  for i in range(4):r=add(r,mul(co[i],A[i][j]))
  assert not r
 return dict(all_three_universal_kappa_A_columns_zero=True,formal_independent_entries=12,signed_cofactor_terms=[len(c)for c in co],known6_proof='Accepted90 two-minors +nonzero anchor polynomial imply all6 columns proportional in fractionfield; cofactor kills anchor; domain descent gives all6 identities, no pointwise inverse.')
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);assert out.is_dir()and not(out/'report.json').exists();phase=read(out/'phase.json');pins={};cache={}
 def source(rel):
  p=ROOT/rel
  if rel not in cache:cache[rel]=read(p);pins[rel]=sha(p)
  return cache[rel]
 base='reports/S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_PREFLIGHT_20261104_v63';gate=source(base+'/report.json');cert=source(base+'/certificate.json');op=source(base+'/operator-and-cost.json');assert gate['resource_gate_passes']and gate['total_live_upper']==119119 and sha(ROOT/(base+'/report.json'))=='73d433857d810e49ce32b74756c14928b8ed7033d7b94a5eed83d0fcfc9d64fc'
 peer=source('reports/S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW/report.json');assert peer['all_passed']
 N=op['nodes'];needed=set(cert['needed_nodes']);assert len(N)==240 and len(needed)==238
 q=source('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json');cp=source('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');pd=source('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json')
 def atom(n):
  r=n['reference']
  if 'constant'in r:return [[[0,0,0],r['constant']]]
  if 'named_FULLU_unit'in r:
   assert r['named_FULLU_unit']=='B';return [[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
  if r['path']=='accepted_Q2L_source_atom':
   d=r['definition'];return cp['matrix25x32'][d[1]][d[2]][0]if d[0]=='checkpoint'else pd['Lframe_dual6x2'][d[1]][1][0]if d[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][d[1]][0]
  x=source(r['path'])[r['key']]
  for i in r.get('index',[]):x=x[i]
  return x
 ip=out/'native-input.txt';atoms=[]
 with ip.open('w')as f:
  f.write('10456 240\n')
  for i,n in enumerate(N):
   if i not in needed:f.write('x\n');continue
   prefix=f"{n['bound']['support']if n['bound']else 0} {n['QQ_denominator']}"
   if n['op']=='atom':
    p=atom(n);m=min((e[0]for e,c in p),default=0);s=n['a_shift']+m;f.write(f"p {prefix} {s} {n['QQ_unit']}\n{len(p)}\n")
    for e,c in p:f.write(f'{e[0]-m} {e[1]} {e[2]} {c}\n')
    atoms.append(dict(node=i,reference=n['reference'],QQ_unit=n['QQ_unit'],a_shift=s,terms=len(p)))
   else:f.write(dict(add='a',sub='s',mul='m')[n['op']]+' '+prefix+' '+' '.join(map(str,n['args']))+'\n')
  f.write(' '.join(map(str,[op['output']]+op['cofactor_roots']+op['cleared_column_roots']))+'\n')
 formal=universal_cofactor_identity();dump(out/'universal-cofactor-proof.json',formal)
 dump(out/'source-binding.json',dict(source_current15_rows=op['source_current15_rows'],original32_target=2,common_B_clearing_power=op['common_B_clearing_power'],cofactor_roots=op['cofactor_roots'],cleared_column_roots=op['cleared_column_roots'],output=op['output'],literal_atoms=atoms,all_4_row_composition_retained=True,all8_extra_tails_retained_in_parent=True,projected_source_ONLY=True,fullrawH7_eta_Cdelta_HZ_HW_lift_SEPARATE=True,known6_zero_peer_path='reports/S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW/report.json'))
 dump(out/'preflight.json',dict(core_and_parser_bound=119119,post_reverse_subtraction_bound=112416,full_job_bound=119119,all_output_and_input_atom_copy_bound=2*gate['raw_artifact_bytes_upper']+1000000,post_reproduct='Save/release core DAG, reread answer, reverse-subtract four products, release inputs before subtraction; every partial difference is in same21740 lattice.'))
 phase.update(phase='source/universal identities/preflight pinned; before one native column2');dump(out/'phase.json',phase);cache.clear();del q,cp,pd,N,atoms
 cpp=Path(__file__).with_suffix('.cpp');exe=out/'native-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)'
 run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,56-(time.monotonic()-started)));(out/'compile.log').write_text(run.stdout);assert run.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(ip)+' '+str(out);run=subprocess.run(['bash','-lc',cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=max(1,58-(time.monotonic()-started)));(out/'native.log').write_text(run.stdout)
 p=out/'native.json'if(out/'native.json').exists()else out/'frontier.json';n=read(p)if p.exists()else dict(all_passed=False,returncode=run.returncode);ok=run.returncode==0 and n['all_passed'];phase.update(phase='native completed FINAL'if ok else'one native stopped; no automaticretry');dump(out/'phase.json',phase)
 for p in [Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_h22_pure3_fixed_minor_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 size=sum(p.stat().st_size for p in out.iterdir()if p.is_file()and p!=exe);assert size<32*1024*1024
 result=dict(id=ID,all_passed=ok,frozen=True,independent_review_pending=True,output_is_zero=n.get('output_is_zero'),output_terms=n.get('output_terms'),max_live_slots=n.get('max_total_live_coefficient_slots'),all4_signed_cofactors_and_source_columns_saved=n.get('all4_signed_cofactors_and_all4_cleared_columns_saved'),reverse_product_verification=n.get('independent_reverse_four_product_reproduct_zero'),universal_kappa_A_and_known6_zero_proof_pinned=True,projected_source_ONLY_fullrawH7_gap_explicit=True,no_fullrank_or_global57_or_othercolumn_claim=True,raw_artifact_bytes_before_report_excluding_binary=size,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p!=exe},seconds=time.monotonic()-started)
 dump(out/'report.json',result);signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
