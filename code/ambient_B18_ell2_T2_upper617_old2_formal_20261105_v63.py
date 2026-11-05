#!/usr/bin/env python3
"""FORMAL combined OLD2 symbol coefficients only; no source-response values."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63'
GATE='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_INPUT_GATE_20261105_v63/'
GATE_SHA='e193894d22d0da56a2420c3441c7aecdfd523405c84149316104713cde48c0ce'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,F)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out,peer,custody):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);pins={};cursor={};bits=0;records=[];live_state={}
 def read(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def dump(name,x):
  (out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
  assert sum(p.stat().st_size for p in out.iterdir())<d['future_formal_RAW_cap']
 r=read(GATE+'report.json');assert r['all_passed'] and pins[GATE+'report.json']==GATE_SHA
 ip=GATE+'617-combined-OLD2-independent-symbol-positive-input.json';d=read(ip)
 assert pins[ip]==r['evidence_sha256'][Path(ip).name]
 p=read(peer);c=read(custody)
 assert p['all_passed'] and GATE_SHA in json.dumps(p) and c['all_passed'] and c['safe_help_literal_PASS']
 assert c['worker_sha256']==sha(Path(__file__)) and c['gate_sha256']==GATE_SHA
 for path,h in c['literal_input_sha256'].items():assert sha(ROOT/path)==h
 assert count([r,p,c,pins])<1000 and count(d)<=d['prior_input_count']+1000
 assert d['source_committed_bits']==120 and d['source_temporary_bits']==361 and d['future_formal_numeric_cap']==49698
 assert d['source_proofs']['no_leaf_second_zero_filled']
 dump('phase.json',{'phase':'authorized_formal_combined_OLD2_before_signed_symbol_fold','no_actual_second_firstcross_native_or_upper_response':True,'gate_sha256':GATE_SHA})
 def ck(q):
  nonlocal bits
  b=max(abs(q.numerator).bit_length(),q.denominator.bit_length());bits=max(bits,b)
  if b>120:
   dump('guard-frontier.json',{'cursor':cursor,'uncommitted_formal_QQ':str(q),'bits':b,'live_state':live_state,'completed_formal_records':records,'source_replay_pins':pins,'no_cap_raise_retry_or_response_values':True})
   raise RuntimeError('prior120 STOP before use')
  return q
 def plus(a,b):return ck(a+b)
 def times(a,b):return ck(a*b)
 def madd(z,k,q):
  q=plus(z.get(k,F(0)),q)
  if q:z[k]=q
  else:z.pop(k,None)
 def axpy(z,a,k):
  for j,q in a.items():madd(z,j,times(k,q))
 def ser(z):return [[list(k) if isinstance(k,tuple) else k,str(q)] for k,q in sorted(z.items())]
 grad={}
 for i,n in d['linear_typed_program104']:
  cursor={'stage':'signed_linear_projection_only','node':i};live_state={'gradient_nodes_done':len(grad)}
  if n[0]=='actual_C0_scalar_function':z={n[1]:F(1)}
  elif n[0]=='sum':
   z={}
   for j,a,b in n[1]:axpy(z,grad[j],F(a,b))
  else:assert n[0]=='product' and len(n[3])==1;z={};axpy(z,grad[n[3][0]],F(n[1],n[2]))
  grad[i]=z
 for i,v in d['positive_linear_node_gradients']:
  env={f:(int(den),int(u)) for f,den,u in v}
  assert set(grad[i])<=set(env)
  for f,q in grad[i].items():den,u=env[f];assert den%q.denominator==0 and abs(q.numerator)*(den//q.denominator)<=u
 outer=d['literal_outer_four_pivot_source'];G=[{i:F(q) for i,q in v} for v in outer['Gamma_two']]
 cg={i:{i:F(1)} for i,kind,point in d['literal_function_origins71'] if F(point)==0}
 assert len(cg)==71
 for i,n in outer['literal_C0_program']:
  cursor={'stage':'literal_outer_old_minus_new_symbol_route','node':i}
  if n[0]=='g':assert i in cg
  elif n[0]=='n':cg[i]={f:times(F(-1),q) for f,q in cg[n[1]].items()}
  else:
   assert n[0]=='a';cg[i]=dict(cg[n[1]]);axpy(cg[i],cg[n[2]],F(1))
 D=[[[F(q) for q in v] for v in ch] for ch in d['D_two_by10_by56']]
 cols=dict(d['linear_D_old_root_columns45']);c619=[{j:F(q) for j,q in v if j<640} for v in d['original619_c0_source_combinations']]
 a619={j:{i:F(q) for i,q in v} for j,v in d['raw619_old4_plane_a1_nonempty28coefficients']}
 pcols=dict(outer['four_source_root_columns']);old4=outer['original4_old_pivots']
 bounds={(case,ch):{(p,f):(int(den),int(u)) for p,f,den,u in v} for case,ch,v in d['positive617_combined_symbol_coefficients']}
 nonzero=0;symbolset=set();actualslots=0
 for rec in d['kernel617_source_recipes']:
  j=rec['free_case619'];cursor={'stage':'617_source_subtraction_before_rootfold','case':j}
  cc=dict(c619[j]);aa=dict(a619.get(j,{}))
  for basis,coef in zip([563,609],rec['subtract_actual_basis563609'],strict=True):
   q=F(coef)
   if q:axpy(cc,c619[basis],-q);axpy(aa,a619.get(basis,{}),-q)
  channels=[]
  for ch in range(2):
   roots={};live_state={'case':j,'channel':ch,'corrected_old_source':ser(cc),'corrected_old4_a1':ser(aa)}
   cursor={'stage':'identical_root_fold_before_atom_gradient','case':j,'channel':ch}
   for old,x in sorted(cc.items()):
    for row,node in cols.get(old,[]):
     for pair in range(10):
      if D[ch][pair][row]:madd(roots,('typed',node,pair),times(x,D[ch][pair][row]))
   for slot,x in sorted(aa.items()):
    pi,var=divmod(slot,9);assert pi<4 and var in [7,8];pair=8 if var==8 else 9
    for row,(kind,node) in pcols[old4[pi]]:
     if row in G[ch]:madd(roots,(kind,node,pair),times(-x,G[ch][row]))
   assert len(roots)<=d['largest_per_case_rootweight_support']
   live_state['combined_identical_roots']=ser(roots)
   z={};cursor={'stage':'independent_function_second_symbol_fold','case':j,'channel':ch}
   for (kind,node,pair),q in sorted(roots.items()):
    for f,a in (grad[node] if kind=='typed' else cg[node]).items():madd(z,(pair,f),times(q,a))
   env=bounds.get((j,ch),{});assert set(z)<=set(env)
   for k,q in z.items():
    den,u=env[k];assert den%q.denominator==0 and abs(q.numerator)*(den//q.denominator)<=u
   v=[[pair,f,str(q)] for (pair,f),q in sorted(z.items())];channels.append(v)
   actualslots+=len(v);nonzero+=bool(v);symbolset.update(f for pair,f,q in v)
  records.append([j,channels])
  # The record reserve also covers the JSON codec copy; the 2*input reserve
  # covers all independently parsed maps, including the bound table.
  actual_owned=2*count(d)+2*count(records)+d['future_scratch_upper']+2000
  assert actual_owned<d['future_formal_numeric_cap']
 ev={'all617_combined_OLD2_independent_symbol_coefficients':records,'ordinary_pair_order10':d['pair_order10'],'all71_literal_function_classifications_and_zero_origins':d['literal_function_origins71'],'actual_symbol_function_ids':sorted(symbolset),'possible1926_slots_all_checked':True,'identical_root_fold_after617_source_subtraction_before_atom_gradient':True,'scope':'FORMAL D Xc2 minus OLD4 GammaXp2a1 only. No actual second, firstcross, native, b1 or upper response.','all_physical_and_only79_newg_symbols_independent':True,'native23_L2S0_and_remaining_source_terms_not_zero_filled':True}
 dump('actual617-formal-combined-OLD2-second-symbol-coefficients.json',ev)
 owned=2*count(d)+2*count(ev)+d['future_scratch_upper']+2000;assert owned<d['future_formal_numeric_cap']
 report={'id':STEM,'all_passed':True,'formal_only_no_actual_second_or_upper_response':True,'cases':len(records),'nonzero_case_channels':nonzero,'actual_formal_slots':actualslots,'actual_symbol_function_ids':sorted(symbolset),'actual_max_bits':bits,'owned_upper':owned,'source_positive_bits120':True,'ordinary_no_factor2':True,'input_sha256':pins,'evidence_sha256':{'actual617-formal-combined-OLD2-second-symbol-coefficients.json':sha(out/'actual617-formal-combined-OLD2-second-symbol-coefficients.json')},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','cases','nonzero_case_channels','actual_formal_slots','actual_symbol_function_ids','actual_max_bits','owned_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--source-peer',required=True);p.add_argument('--custody',required=True);a=p.parse_args();main(a.output,a.source_peer,a.custody)
