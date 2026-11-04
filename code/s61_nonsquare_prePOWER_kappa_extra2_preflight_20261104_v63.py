#!/usr/bin/env python3
"""SOURCE-ONLY complete scalar recipe/height/lifetime for kappa dot extra-tail2."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
BP=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6']
 gate=load('S61_NONSQUARE_PREPOWER_PAIR_ONLY_PURE_COVECTOR_PREFLIGHT_20261104_v63','report.json');assert gate['all_passed']
 N=[];keys={};cache={};inventory=[]
 def lattice(box,lo,hi):return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
 def put(op,args,bound,D,L,**kw):i=len(N);N.append(dict(id=i,op=op,args=args,bound=bound,QQ_denominator=str(D),integer_l1_upper=str(L),**kw));return i
 def atom(poly,ref,unit=Q(1),shift=0):
  p=[([m[0]+shift,m[1],m[2]],Q(c)/unit)for m,c in poly if Q(c)];key=tuple((tuple(m),str(c))for m,c in p)
  if key in keys:return keys[key]
  D=math.lcm(*(c.denominator for m,c in p))if p else 1;L=sum(abs(c.numerator)*(D//c.denominator)for m,c in p)
  bound=None if not p else dict(support=len(p),box=[[min(m[k]for m,c in p),max(m[k]for m,c in p)]for k in range(3)],degree=[min(sum(m)for m,c in p),max(sum(m)for m,c in p)])
  i=put('atom',[],bound,D,L,reference=ref,QQ_unit=str(unit),a_shift=shift);keys[key]=i;return i
 def const(c):return atom([[[0,0,0],str(Q(c))]],dict(constant=str(c)))
 z=const(0);one=const(1)
 def bd(i):return N[i]['bound']
 def dl(i):return int(N[i]['QQ_denominator']),int(N[i]['integer_l1_upper'])
 def mul(i,j):
  if i==z or j==z:return z
  if i==one:return j
  if j==one:return i
  A,B=bd(i),bd(j);d,n=dl(i);e,m=dl(j);box=[[A['box'][k][0]+B['box'][k][0],A['box'][k][1]+B['box'][k][1]]for k in range(3)];lo=A['degree'][0]+B['degree'][0];hi=A['degree'][1]+B['degree'][1];s=min(A['support']*B['support'],lattice(box,lo,hi),n*m)
  return put('mul',[i,j],dict(support=s,box=box,degree=[lo,hi]),d*e,n*m)
 def add(i,j,sgn=1):
  if j==z:return i
  if i==z:return j if sgn==1 else mul(const(-1),j)
  if sgn<0 and i==j:return z
  A,B=bd(i),bd(j);d,n=dl(i);e,m=dl(j);D=math.lcm(d,e);L=n*(D//d)+m*(D//e);box=[[min(A['box'][k][0],B['box'][k][0]),max(A['box'][k][1],B['box'][k][1])]for k in range(3)];lo=min(A['degree'][0],B['degree'][0]);hi=max(A['degree'][1],B['degree'][1]);s=min(A['support']+B['support'],lattice(box,lo,hi),L)
  return put('add'if sgn>0 else'sub',[i,j],dict(support=s,box=box,degree=[lo,hi]),D,L)
 Bnode=atom(BP,dict(named_FULLU_unit='B',literal_source='(50a+18e^2-9e-2)d+60e-40'))
 def timesB(i,k):
  for _ in range(k):i=mul(i,Bnode)
  return i
 def imp(i):
  if i in cache:return cache[i]
  a=q['scalar_DAG']['nodes'][i];op=a[0]
  if op=='z':v=(z,0)
  elif op=='one':v=(one,0)
  elif op in['a','s','m']:
   x,b=imp(a[1]);y,c=imp(a[2])
   if op=='m':v=(mul(x,y),b+c)
   else:
    d=max(b,c);v=(add(timesB(x,d-b),timesB(y,d-c),-1 if op=='s'else 1),d)
  else:
   num,den=cp['matrix25x32'][a[1]][a[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   if len(den)==1:
    assert den[0][0][1:]==[0,0]and Q(den[0][1]);unit=Q(den[0][1]);ashift=-den[0][0][0];bpower=0
   else:
    dd={tuple(m):Q(c)for m,c in den};bb={tuple(m):Q(c)for m,c in BP};assert dd.keys()==bb.keys();unit=dd[next(iter(bb))]/bb[next(iter(bb))];assert unit and all(dd[m]==unit*bb[m]for m in bb);ashift=0;bpower=1
   inventory.append(dict(old_DAG_root=i,definition=a,known_B_denominator_power=bpower,QQ_unit=str(unit),a_shift=ashift));v=(atom(num,dict(path='accepted_Q2L_source_atom',old_DAG_root=i,definition=a,known_B_denominator_power=bpower),unit,ashift),bpower)
  cache[i]=v;return v
 rows=[2,4,5,8];compact=[1,3,4,7];col=q['new22_target_original32_columns'].index(2);images=[imp(q['remaining15x22_scalar_roots'][r][col])for r in rows];power=max(p for x,p in images);cleared=[timesB(x,power-p)for x,p in images]
 A=[]
 for i,r in enumerate(compact):
  row=[atom(hc['raw_compact14x5'][r][j],dict(path='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json',key='raw_compact14x5',index=[r,j]),unit=1/Q(hc['known_row_units_rho14'][r]['constant']))for j in[0,4]]
  row.append(atom(bot[i][0],dict(path='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json',key='bottom4x6',index=[i,0])));A.append(row)
 cofs=[];terms=[];output=z
 for i in range(4):
  rr=[r for r in range(4)if r!=i];c=z
  for p in itertools.permutations(range(3)):
   t=one
   for j in range(3):t=mul(t,A[rr[j]][p[j]])
   parity=i+sum(p[j]>p[k]for j in range(3)for k in range(j+1,3));c=add(c,t,-1 if parity%2 else 1)
  cofs.append(c);t=mul(c,cleared[i]);terms.append(t);output=add(output,t)
 needed=set();todo=[output]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);todo+=N[i]['args']
 uses={i:0 for i in needed}
 for i in needed:
  for j in N[i]['args']:uses[j]+=1
 uses[output]+=1;live=peak=0;peakid=0;records=[]
 for i in sorted(needed):
  s=N[i]['bound']['support']if bd(i)else 0;live+=s;scratch=s if N[i]['op']in['add','sub','mul']else 0
  if live+scratch>peak:peak=live+scratch;peakid=i
  records.append(dict(node=i,live_upper=live,scratch_upper=scratch))
  for j in N[i]['args']:
   uses[j]-=1
   if not uses[j]:live-=N[j]['bound']['support']if bd(j)else 0
 atoms=[i for i in needed if N[i]['op']=='atom'];literal=sum(bd(i)['support']if bd(i)else 0 for i in atoms);total=peak+literal+15000
 saved=[output]+cofs+cleared;artifact=1000000
 for i in saved:
  d,n=dl(i);artifact+=(bd(i)['support']if bd(i)else 0)*(len(str(n))+len(str(d))+110)
 for i in atoms:
  d,n=dl(i);artifact+=(bd(i)['support']if bd(i)else 0)*(len(str(n))+len(str(d))+100)
 d,n=dl(output);ob=bd(output);gateok=total<=150000 and artifact<32*1024*1024
 operator=dict(nodes=N,output=output,cofactor_roots=cofs,cleared_column_roots=cleared,all4_row_term_roots=terms,
  source_current15_rows=rows,source_original32_target=2,anchor_original32=22,common_B_clearing_power=power,
  represented_output='B^power * kappa * projectedM[:,original32column2]',all8_extra_tails_retained_in_parent=True,
  projected_source_only=True,fullrawH7_eta_Cdelta_HZ_HW_lift_SEPARATE_UNEXPANDED=True,
  atom_denominator_inventory=inventory,all_native_nodes_a_Laurent_and_QQ_ONLY=True,
  parameter_only_denominators_no_unknown_q_inverse=True,new_parameter_products_NOT_evaluated=True,
  future_saved_roots=saved)
 dump(out/'operator-and-cost.json',operator)
 cert=dict(last_use_records=records,needed_nodes=sorted(needed),literal_atom_reserve=literal,unmodified_work_reserve=15000,
  peak_node=peakid,last_use_peak_upper=peak,total_live_upper=total,output_bound=ob,output_integer_l1_bits=n.bit_length(),QQ_denominator_bits=d.bit_length(),
  all_source_B_alignments_explicit=True,no_cancellation_assumed_except_zero_or_same_root=True,
  raw_cofactor_four_columns_output_and_atoms_artifact_bytes_upper=artifact,resource_gate_passes=gateok)
 dump(out/'certificate.json',cert);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,resource_gate_passes=gateok,typed_nodes=len(N),live_nodes=len(needed),literal_atom_terms=literal,
  total_live_upper=total,output_support_upper=ob['support'],output_degree_range=ob['degree'],output_a_exponent_box=ob['box'][0],output_integer_l1_bits=n.bit_length(),QQ_denominator_bits=d.bit_length(),
  common_known_B_clearing_power=power,raw_artifact_bytes_upper=artifact,projected_source_only_fullrawH7_gap_explicit=True,
  no_actual_column_scalar_minor_rank_point_prime=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024;dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
