#!/usr/bin/env python3
"""ONE literal-added-source M12 axis+fixed determinant; frozen source/guard scope."""
import argparse,gc,hashlib,importlib.util,json,math,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_M12_AXIS_ENTRY_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
class SourceStop(Exception):pass
def main(gatepath,out):
 started=time.monotonic();signal.alarm(59);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));g=read(gatepath);assert g['all_passed']and g['mathematical_axis_job_NOT_run'];pins={str(gatepath.relative_to(ROOT)):sha(gatepath)}
 for rel,h in g['input_sha256'].items():assert sha(ROOT/rel)==h;pins[rel]=h
 for fn,h in g['source_input_sha256'].items():assert sha(gatepath.parent/fn)==h
 exe=ROOT/'work/d-zero-v63/M12axis/worker';assert exe.exists()and sha(exe)==g['executable_sha256']
 recipe=read(gatepath.parent/'recipe.json');old=read(gatepath.parent/'old-axis579.json');out.mkdir(parents=True,exist_ok=False)
 atoms={};values={};extras=[];C=[[{}for _ in range(12)]for _ in range(12)];cur={};phase='literal_axis_atom';cursor={};literalcells=0;sourcebits=0;sourceproducts=0;peak=0
 def guard(scratch=128):
  nonlocal peak
  n=11000+literalcells+3*(sum(map(len,atoms.values()))+sum(map(len,values.values()))+sum(map(len,extras)))+2*sum(len(p)for row in C for p in row)+scratch
  peak=max(peak,n)
  if n>33275:raise SourceStop('source_numeric_metadata_guard')
  if time.monotonic()-started>20:raise SourceStop('source_time_guard')
 def bit(v):
  nonlocal sourcebits
  sourcebits=max(sourcebits,abs(v).bit_length())
  if abs(v).bit_length()>2048:raise SourceStop('source_integer_bit_guard')
  return v
 def qm(a,b):
  nonlocal sourceproducts
  if abs(a.numerator).bit_length()+abs(b.numerator).bit_length()>2048 or a.denominator.bit_length()+b.denominator.bit_length()>2048:raise SourceStop('source_before_fraction_product_guard')
  sourceproducts+=1;return Q(bit(a.numerator*b.numerator),bit(a.denominator*b.denominator))
 def qa(a,b):
  if max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1>2048:raise SourceStop('source_before_fraction_sum_guard')
  n=bit(bit(a.numerator*b.denominator)+bit(b.numerator*a.denominator));d=bit(a.denominator*b.denominator);return Q(n,d)
 def scale(p,q=Q(1),shift=0):
  nonlocal cur
  cur={}
  for a,v in p.items():
   guard()
   z=qm(v,q)
   if z:cur[a+shift]=z
  return cur
 def add(a,b,sign=1):
  nonlocal cur
  cur=a.copy()
  for d,v in b.items():
   guard();z=qa(cur.get(d,Q(0)),v if sign>0 else -v)
   if z:cur[d]=z
   else:cur.pop(d,None)
  return cur
 def mul(a,b):
  nonlocal cur
  cur={}
  for i,v in a.items():
   for j,w in b.items():
    cursor['degrees']=[i,j];guard();z=qa(cur.get(i+j,Q(0)),qm(v,w))
    if z:cur[i+j]=z
    else:cur.pop(i+j,None)
  return cur
 def check(p,b):
  if b is None:assert not p;return
  assert len(p)<=b['terms']
  for a,v in p.items():assert b['lo']<=a<=b['hi']and abs(v.numerator)<=b['L1']and b['D']%v.denominator==0
 def witness(p):return [[a,str(v)]for a,v in sorted(p.items())]
 def prefix(reason):
  data=dict(version='M12SP001',reason=reason,phase=phase,cursor=cursor,recipe_sha256=sha(gatepath.parent/'recipe.json'),literal_atoms_sha256=sha(gatepath.parent/'literal-atoms67.jsonl'),
   atoms={str(i):witness(p)for i,p in atoms.items()},completed_source_values={str(i):witness(p)for i,p in values.items()},extras=[witness(p)for p in extras],current_uncommitted_polynomial=witness(cur),
   integer_matrix_prefix=[[[a,str(v)]for a,v in sorted(p.items())]for row in C for p in row],no_axis_or_determinant_result_claim=True)
  dump(out/'source-prefix-M12SP001.json',data)
 try:
  BP={(1,0,1):Q(50),(0,2,1):Q(18),(0,1,1):Q(-9),(0,0,1):Q(-2),(0,1,0):Q(60),(0,0,0):Q(-40)}
  with(gatepath.parent/'literal-atoms67.jsonl').open()as f:
   for line in f:
    z=json.loads(line);idx=z['atom'];cursor={'atom':idx};num,den=z['literal_num_den'];literalcells=4*(len(num)+len(den));guard();meta=z['expected'];aa=meta['denominator_a_power'];bb=meta['denominator_B_power'];unit=Q(meta['denominator_QQ_unit']);dd={tuple(m):Q(c)for m,c in den}
    if bb:assert aa==0 and set(dd)==set(BP)and all(dd[m]==unit*BP[m]for m in BP)
    else:assert dd=={(aa,0,0):unit}
    cur={};uq=qm(unit,Q((-40)**bb));inv=Q(uq.denominator,uq.numerator)
    for k,(m,c)in enumerate(num):
     cursor['literal_term']=k;assert m[1]>=0 and m[2]>=0
     if m[1]or m[2]:continue
     guard();v=qm(Q(c),inv);a=m[0]-aa;assert a not in cur
     if v:cur[a]=v
    check(cur,meta['bound']['axis']);atoms[idx]=cur;literalcells=0;del z,num,den;guard()
  assert len(atoms)==67
  phase='source214_axis_DAG'
  with(out/'source214-axis-witness.jsonl').open('w')as w:
   for n in recipe['nodes']:
    i=n['node'];d=n['definition'];cursor={'source_node':i,'opcode':d[0]};guard()
    if n['atom']is not None:p=atoms[n['atom']]
    elif d[0]=='z':p={}
    elif d[0]=='one':p={0:Q(1)}
    elif d[0]=='m':
     p=mul(values[d[1]],values[d[2]]);saved=p.copy();reverse=mul(values[d[2]],values[d[1]]);assert saved==reverse;p=saved
    else:
     sign=-1 if d[0]=='s'else 1;p=add(values[d[1]],values[d[2]],sign);saved=p.copy();assert add(add(saved,values[d[1]],-1),values[d[2]],-sign)=={};p=saved
    check(p,n['axis_bound']);values[i]=p;w.write(json.dumps(dict(node=i,definition=d,axis_polynomial=witness(p),bounds_verified=True,source_reverse_products_backadds=True))+'\n')
  phase='literal_M12_added_adapter';cols=recipe['current22_original_labels'];roots=recipe['current15_roots'];const={int(i):Q(c)for i,c in recipe['constants'].items()}
  def s(r,t):return values[roots[r][cols.index(t)]]
  def keep(p):extras.append(p);return p
  col=[keep(atoms[recipe['h9_atom']]),keep(scale(s(1,17),shift=-1)),keep(scale(s(3,17),shift=-1)),keep(s(6,17)),keep(s(2,17)),keep(scale(s(11,17),1/const[10])),keep(scale(add(s(12,17),s(10,17),-1),1/const[11])),keep(scale(s(9,17),1/const[8])),keep(scale(s(7,17),1/const[6])),keep(scale(s(13,17),1/const[12],shift=1)),keep(scale(s(10,17),1/const[9]))]
  col[1]=add(col[1],scale(col[2],Q(2)));extras[1]=col[1]
  bottom=[keep(s(14,t))for t in[22,11,13,25,15,19,14,18]]+[keep(atoms[i])for i in recipe['H3_atoms_in_order']];corner=keep(s(14,17))
  for p,b in zip(col,recipe['source_and_extra_axis_bounds']['addedcol']):check(p,b)
  for p,b in zip(bottom,recipe['source_and_extra_axis_bounds']['bottom']):check(p,b)
  check(corner,recipe['source_and_extra_axis_bounds']['corner']);guard()
  dump(out/'added23-axis-values.json',dict(column17=[witness(p)for p in col],G134row=[witness(p)for p in bottom],corner=witness(corner),old_precondition_on_newcolumn=True,H3carry='rho13/L=8000000/16119'))
  rowD=list(map(int,old['old_row_LCD']));E=recipe['extra_row_QQ_integer_scales'];D12=int(recipe['new_row_LCD']);rs=recipe['rowdual'];cs=recipe['coldual']
  for a in old['axis_coefficients']:
   i,j=a['row'],a['column'];C[i][j]={int(k):bit(int(v)*E[i])for k,v in a['polynomial']}
  assert sum(len(C[i][j])for i in range(11)for j in range(11))==579
  def integerize(p,factor,shift):
   ans={}
   for a,v in p.items():
    z=qm(v,Q(factor));assert z.denominator==1 and a+shift>=0 and a+shift<=9 and abs(z.numerator).bit_length()<=118
    if z:ans[a+shift]=z.numerator
   return ans
  for i,p in enumerate(col):C[i][11]=integerize(p,rowD[i]*E[i],rs[i])
  for j,p in enumerate(bottom):C[11][j]=integerize(p,D12,cs[j])
  C[11][11]=integerize(corner,D12,0);guard();assert sum(len(p)for row in C for p in row)<=747 and sourcebits<=238
  dump(out/'source-normalization-and-guards.json',dict(all67literal_e_d_nonnegative_and_aBunit_denominators=True,B_axis_image=-40,all214_source_axis_bounds_reverse_products_backadds=True,
   old579_literal_scaled_only=True,all23_rowLCD_and_a_shifts_integer_exact=True,H3carry='rho13/L=8000000/16119, compact13[1,3,2]',no_kappa3_denominator=True,axis_is_coefficient_restriction_NOT_FULLU_point=True,
   source_Fraction_observed_bits=sourcebits,source_numeric_peak=peak,source_products=sourceproducts))
  sec=min(45,54-(time.monotonic()-started));assert sec>1;nativein=out/'native-input.txt'
  with nativein.open('w')as w:
   w.write(f'M12AXIS001 12 108 1096 2048 33275 {sec}\n'+recipe['all_minor_H']+'\n')
   for i,row in enumerate(C):
    for j,p in enumerate(row):w.write(f'{i} {j} {len(p)} '+' '.join(f'{a} {v}'for a,v in sorted(p.items()))+'\n')
  Dtotal=D12*math.prod(d*e for d,e in zip(rowD,E));del atoms,values,extras,C,col,bottom,corner,cur;gc.collect()
 except SourceStop as e:
  prefix(str(e));dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,source_guard_prefix_only=True,no_axis_determinant_result=True,input_sha256=pins));return
 except Exception as e:
  prefix(type(e).__name__+': '+str(e));dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,source_arithmetic_started=True,source_failure_prefix_saved=True,no_axis_determinant_result=True,error=type(e).__name__+': '+str(e),input_sha256=pins));return
 dump(out/'phase.json',dict(phase='source214 and literal579 integer adapter COMPLETE; ONE fixed12 determinant starts',native_input_sha256=sha(nativein),exe_sha256=sha(exe)))
 def limits():resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));resource.setrlimit(resource.RLIMIT_CPU,(52,52))
 rr=subprocess.run([str(exe),str(nativein),str(out),sha(nativein)],capture_output=True,text=True,timeout=min(sec+3,58-(time.monotonic()-started)),preexec_fn=limits);(out/'native.stdout').write_text(rr.stdout);(out/'native.stderr').write_text(rr.stderr)
 parser=ROOT/'code/s61_nonsquare_M12_binary_inspect_20261105_v63.py';ps=importlib.util.spec_from_file_location('m12_binary_reader',parser);pm=importlib.util.module_from_spec(ps);ps.loader.exec_module(pm)
 if rr.returncode:
  state=out/'guard-M12AS001.bin'
  if state.exists():dump(out/'guard-parse.json',pm.parse(state))
  dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,all_binary_patches_and_guard_preserved=True,no_generic_or_axis_result_claim=True,input_sha256=pins));return
 n=read(out/'native.json');assert n['all_passed']and n['owned_with_ALL_metadata_and_input']<=33275 and n['observed_integer_bits']<=2048
 parsed=[pm.parse(out/f'stage-{i}-M12BP001.bin')for i in range(12)];assert all(p['input_sha256']==sha(nativein)for p in parsed);dump(out/'all12-binary-witness-parse.json',parsed)
 qaxis=[[k+5,str(Q(int(v),Dtotal))]for k,v in n['determinant_C12']];assert bool(qaxis)==n['axis_determinant_NONZERO']
 dump(out/'axis-determinant-and-scope.json',dict(detM12_axis_QQ_coefficients=qaxis,exact_normalization='detM12_axis=a^5*detC12/prod(oldrowLCD*extraD)/newrowLCD',axis_rank=n['axis_rank'],
  generic_D12_NONZERO_proved_PRIMARY=bool(qaxis),axis_zero_IS_NOT_generic_zero=True,independent_peer_required_if_positive=True,
  sufficient_chart_only='oldOmega intersect D(D12), via A*D15_norm=-T*D12, no extraPchart; all complements/fullU57/global239/exhaustion OPEN'))
 new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in gatepath.parent.iterdir());assert new+external<=22239243
 rep=dict(id=ID,all_passed=True,frozen=True,ONE_added_axis_source_and_fixed12_test=True,old579_reused=True,added214_source_verified=True,axis_determinant_NONZERO=n['axis_determinant_NONZERO'],axis_rank=n['axis_rank'],
  generic_D12_nonzero_PRIMARY=bool(qaxis),axis_zero_inconclusive_for_generic=True,source_fraction_bits=sourcebits,native_integer_bits=n['observed_integer_bits'],whole_numeric_peak=max(peak,n['owned_with_ALL_metadata_and_input']),
  all12_binary_witnesses_complete_and_parsed=True,no_new_GCD_fieldrank_point_prime_or_source572=True,no_global57_global239_or_exhaustion_claim=True,
  integer_products=n['integer_products'],new_RAW=new,external_input_RAW=external,TOTAL_RAW=new+external,positive_independent_peer_pending=bool(qaxis),input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--gate',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.gate.resolve(),a.output.resolve())
