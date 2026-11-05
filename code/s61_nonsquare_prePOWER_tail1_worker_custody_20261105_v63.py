#!/usr/bin/env python3
"""Prepare ONLY copied literal atoms/native instructions; never run the tail1 value."""
import argparse, hashlib, json, math, mmap, resource, signal, subprocess, time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_TAIL1_WORKER_CUSTODY_20261105_v63'
CAP=9517864
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def numeric(x):
 if type(x)is int:return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except ValueError:return 0
 if isinstance(x,list):return sum(map(numeric,x))
 if isinstance(x,dict):return sum(map(numeric,x.values()))
 return 0
class Spans:
 """Lexical byte spans only: skip old checkpoint values without decoding them."""
 def __init__(self,p):self.f=p.open('rb');self.b=mmap.mmap(self.f.fileno(),0,access=mmap.ACCESS_READ);self.cache={}
 def ws(self,i):
  while self.b[i]in b' \r\n\t':i+=1
  return i
 def end(self,i):
  i=self.ws(i);c=self.b[i]
  if c==34:
   j=i+1
   while True:
    if self.b[j]==92:j+=2
    elif self.b[j]==34:return j+1
    else:j+=1
  if c in [91,123]:
   depth=1;j=i+1
   while depth:
    c=self.b[j]
    if c==34:j=self.end(j);continue
    if c in [91,123]:depth+=1
    elif c in [93,125]:depth-=1
    j+=1
   return j
  j=i
  while j<len(self.b)and self.b[j]not in b',]} \r\n\t':j+=1
  return j
 def obj(self,key,span=None):
  i=self.ws(0 if span is None else span[0]);assert self.b[i]==123;i+=1
  while True:
   i=self.ws(i)
   if self.b[i]==125:raise KeyError(key)
   e=self.end(i);k=json.loads(self.b[i:e]);i=self.ws(e);assert self.b[i]==58;i=self.ws(i+1);e=self.end(i)
   if k==key:return i,e
   i=self.ws(e)
   if self.b[i]==44:i+=1
 def arr(self,span):
  if span in self.cache:return self.cache[span]
  i=self.ws(span[0]);assert self.b[i]==91;i+=1;r=[]
  while True:
   i=self.ws(i)
   if self.b[i]==93:break
   e=self.end(i);r.append((i,e));i=self.ws(e)
   if self.b[i]==44:i+=1
  self.cache[span]=r;return r
 def value(self,span):return json.loads(self.b[span[0]:span[1]])
 def index(self,span,*ij):
  for i in ij:span=self.arr(span)[i]
  return span
 def close(self):self.b.close();self.f.close()
def main(out):
 began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False);pins={}
 def pin(rel):p=ROOT/rel;pins[rel]=sha(p);return p
 def dump(name,x):(out/name).write_text(json.dumps(x,indent=2)+'\n')
 base='reports/S61_NONSQUARE_PREPOWER_TAIL8_COMPLEMENT_SOURCE_20261105_v63/'
 r=json.loads(pin(base+'report.json').read_text());assert r['all_passed']and r['selected_source_candidate']==1
 family=Spans(pin(base+'typed-family-and-source-plans.json'))
 plan=family.value(family.index(family.obj('source_plans'),0));family.close();assert plan['target']==1
 NN=plan['nodes'];cols=plan['column_roots'];ks=plan['cached_kappa_nodes']
 assert cols==[66,94,122,150]and ks==[151,152,153,154]
 needed=set();todo=cols[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);todo+=NN[i]['args']
 order=sorted(needed);assert len(order)==149
 atoms=sorted([i for i in needed if NN[i]['op']=='external_atom']+ks)
 assert len(set(atoms))==len(atoms)
 cp=Spans(pin('reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json'))
 pd=Spans(pin('reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json'))
 q=Spans(pin('reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63/certificate.json'))
 cps=cp.obj('matrix25x32');pds=pd.obj('Lframe_dual6x2');qds=q.obj('scalar_DAG');qns=q.obj('nodes',qds);qnew=q.obj('new_unit_atoms_num_den',qds)
 B=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
 finalD=int(NN[plan['output_root']]['QQ_D']);assert finalD.bit_length()==197
 allcoef=all_den=colcoef=kcoef=0;maxatom=0;prepbits=0;atom_witness=[]
 native=out/'native-input.txt';literal=out/'literal-atoms58.jsonl'
 with native.open('w')as f,literal.open('w')as g:
  f.write('J1K001 155 149\n'+' '.join(map(str,cols))+'\n'+' '.join(map(str,ks))+'\n'+str(finalD)+'\n'+' '.join(map(str,order))+'\n')
  for i,n in enumerate(NN[:155]):
   op={'external_atom':0,'add_bound_only':1,'sub_bound_only':2,'mul_bound_only':3}[n['op']]
   args=n['args']if op else [-1,-1];sup=n['bound']['support']if n['bound']else 0
   bits=int(n['integer_l1_upper']).bit_length()
   f.write(f'{i} {op} {args[0]} {args[1]} {sup} {bits} {n["QQ_D"]}\n')
  f.write(str(len(atoms))+'\n')
  for i in atoms:
   n=NN[i];ref=n['reference'];den=[[[0,0,0],'1']]
   if 'kind'in ref:
    definition=ref['definition'];assert q.value(q.index(qns,ref['old_node']))==definition
    if ref['kind']=='checkpoint':p=cp.value(cp.index(cps,*definition[1:]))
    elif ref['kind']=='accepted_pair_second':p=pd.value(pd.index(pds,definition[1],1))
    else:assert ref['kind']=='new_unit_atom';p=q.value(q.index(qnew,definition[1]))
    raw,den=p
   elif 'cached_kappa_root'in ref:
    kpath=pin(ref['path']);raw=json.loads(kpath.read_text())['polynomial']
   elif 'named_FULLU_unit'in ref:assert ref['named_FULLU_unit']=='B';raw=B
   else:raw=[]if ref['QQ']=='0'else [[[0,0,0],ref['QQ']]]
   allcoef+=len(raw);all_den+=len(den)
   if i in ks:kcoef+=len(raw)
   else:colcoef+=len(raw)
   maxatom=max(maxatom,numeric(raw)+numeric(den));unit=Q(n['QQ_divisor']);D=int(n['QQ_D']);L=0
   record=dict(node=i,reference=ref,original_num_den=[raw,den],QQ_divisor=n['QQ_divisor'],a_shift=n['a_shift'])
   g.write(json.dumps(record,separators=(',',':'))+'\n')
   f.write(f'{i} {len(raw)}\n')
   for m,c in raw:
    cc=Q(c);numb=cc.numerator*unit.denominator*D;denb=cc.denominator*unit.numerator
    assert denb and numb%denb==0;v=numb//denb
    # Both the exact division and its reverse use literal integer bounds.
    h=max(abs(cc.numerator).bit_length()+unit.denominator.bit_length()+D.bit_length(),
          cc.denominator.bit_length()+abs(unit.numerator).bit_length(),
          abs(v).bit_length()+cc.denominator.bit_length()+abs(unit.numerator).bit_length())+2
    prepbits=max(prepbits,h);assert h<2048
    assert v*denb==numb
    a,e,d=m[0]+n['a_shift'],m[1],m[2];assert e>=0 and d>=0
    assert -128<=a<=128 and e<=128 and d<=128
    f.write(f'{a*1048576+e*1024+d} {v}\n');L+=abs(v)
   assert L==int(n['integer_l1_upper'])
   atom_witness.append(dict(node=i,terms=len(raw),integer_l1=L,denominator=D,normalization_reverse_exact=True))
 assert colcoef==9268 and kcoef==4920 and allcoef==14188
 spans=2*sum(len(v)for ss in [cp,pd,q]for v in ss.cache.values())
 for ss in [cp,pd,q]:ss.close()
 maxNN=max(int(n['integer_l1_upper']).bit_length()for n in NN);assert maxNN==359
 maxDD=max([finalD.bit_length()]+[int(n['QQ_D']).bit_length()for n in NN[:155]])
 for i in order:
  n=NN[i]
  if not n['args']:continue
  x,y=(NN[j]for j in n['args']);di,dj=int(x['QQ_D']),int(y['QQ_D']);li,lj=int(x['integer_l1_upper']),int(y['integer_l1_upper'])
  D=int(n['QQ_D']);L=int(n['integer_l1_upper'])
  if n['op']=='mul_bound_only':assert D==di*dj and L==li*lj
  else:assert D%di==D%dj==0 and L==li*(D//di)+lj*(D//dj)
 srcpairs=sum(NN[n['args'][0]]['bound']['support']*NN[n['args'][1]]['bound']['support']for n in NN[:151]if n['op']=='mul_bound_only')
 dotpairs=r['chosen_candidate_profile']['direct_accumulator_pair_upper'];assert dotpairs==14283624
 dump('selected-plan.json',plan)
 dump('literal-normalization-witness.json',dict(atom_witnesses=atom_witness,all_copied_numerator_coefficients=14188,
  column_atoms=9268,cached_kappa_atoms=4920,denominator_record_terms=all_den,
  all_original_records_preserved=True,no_new_column_or_contraction_values=True,
  source_scalar_order_and_B_clearing_unchanged=True))
 cpp=pin('code/s61_nonsquare_tail1_cached_kappa_20261105_v63.cpp');driver=pin('code/s61_nonsquare_prePOWER_tail1_cached_kappa_entry_20261105_v63.py')
 me=pin(str(Path(__file__).relative_to(ROOT)))
 cmd=['bash','-c','source /workspace/cloud-setup/activate.sh && c++ -std=c++17 -O2 "$1" -o "$2" -lflint','tail1-build',str(cpp),str(out/'native-worker')]
 with (out/'compile.log').open('w')as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=30)
 helpcmd=[str(out/'native-worker'),'--help']
 h=subprocess.run(helpcmd,text=True,capture_output=True,check=True,timeout=5);(out/'worker-help.log').write_text(h.stdout+h.stderr)
 h=subprocess.run([str(ROOT/'.venv/bin/python'),str(driver),'--help'],text=True,capture_output=True,check=True,timeout=5);(out/'driver-help.log').write_text(h.stdout+h.stderr)
 prepowned=numeric(plan)+spans+maxatom+15000+2000
 assert prepowned<=91676
 own=dict(limit91676=91676,coefficient_and_packed_key_cells_per_term=2,native_ALL_metadata_and_scalar_scratch_reserve=15100,
  native_live_terms_limit=(91676-15100)//2,column_source_profile_bound=r['chosen_candidate_profile']['phase_owned']['column'],
  direct_dot_exact_profile_bound=91676,preparation_owned_upper=prepowned,
  prepared_numeric_plan_cells=numeric(plan),byte_span_index_cells=spans,max_one_literal_atom_numeric_cells=maxatom,
  preparation_integer_normalization_temp_bits_upper=prepbits,native_coefficient_bits_upper=359,
  native_denominator_temp_bits_upper=2*maxDD+1,hard_integer_temp_ceiling_for_preparation=2048,
  one_source_column_closure_nodes=149,source_integer_mul_pairs_upper=srcpairs,
  one_dot_pairs_upper=dotpairs,three_dot_passes_pairs_upper=3*dotpairs,
  three_passes='forward; reverse four-term subtraction to literal zero; forward reproduction. One accumulator, no product/second result map.',
  soft_seconds45=45,hard_seconds60=60,hard_RAM_bytes=2000000000,
  source_order='all149 column nodes; four cached kappa leaves; three fixed dot passes',
  no_value_execution=True)
 dump('ownership-and-worker-proof.json',own)
 raw=sum(p.stat().st_size for p in out.iterdir())+cpp.stat().st_size+driver.stat().st_size+me.stat().st_size
 future=raw+4_700_000+100_000
 assert future<CAP,(raw,future)
 dump('RAW-closure.json',dict(prepared_actual_bytes_plus_external3_codes=raw,prefix_OR_success_file_upper=4700000,
  all_future_receipts_report_logs_reserve=100000,total_upper=future,unchanged_cap=CAP,
  one_state_file_no_atomic_doublebuffer=True,no_success_output_before_all3_passes=True,
  no_Q2_checkpoint_or_old572_read_by_scalar_worker_or_driver=True))
 report=dict(id=ID,all_passed=True,frozen=True,prepared_only=True,source_gate_sha256=sha(ROOT/(base+'report.json')),
  literal_column_terms=9268,cached_kappa_terms=4920,column_DAG_nodes=149,source_only_no_values=True,
  worker_and_driver_help_passed=True,worker_compiled=True,value_worker_NOT_run=True,
  strict_owned_limit=91676,native_coefficient_bits=359,preparation_temp_bits=prepbits,
  all_phase_RAW_upper=future,strict_RAW_limit=CAP,independent_source_code_review_pending=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-began)
 dump('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
