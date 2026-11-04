#!/usr/bin/env python3
"""SOURCE ONLY one kappa14 gate, separately profiled raw-proportionality tests."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_KAPPA14_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
BP=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
@lru_cache(None)
def lattice(box,lo,hi):return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];ke='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63';kb=load(ke,'source-binding.json');assert load(ke+'_REVIEW','report.json')['all_passed'];assert load('S61_NONSQUARE_PREPOWER_KAPPA11_SCHUR_20261104_v63_REVIEW','report.json')['all_passed']
 N=[];keys={};cache={};inventory=[]
 def put(op,args,b,D,L,**kw):i=len(N);N.append(dict(id=i,op=op,args=args,bound=b,QQ_denominator=str(D),integer_l1_upper=str(L),**kw));return i
 def atom(poly,ref,unit=Q(1),shift=0):
  p=[([m[0]+shift,m[1],m[2]],Q(c)/unit)for m,c in poly if Q(c)];key=tuple((tuple(m),str(c))for m,c in p)
  if key in keys:return keys[key]
  D=math.lcm(*(c.denominator for m,c in p))if p else 1;L=sum(abs(c.numerator)*(D//c.denominator)for m,c in p);b=None if not p else dict(support=len(p),box=[[min(m[k]for m,c in p),max(m[k]for m,c in p)]for k in range(3)],degree=[min(sum(m)for m,c in p),max(sum(m)for m,c in p)])
  i=put('atom',[],b,D,L,reference=ref,QQ_unit=str(unit),a_shift=shift);keys[key]=i;return i
 def const(c):return atom([[[0,0,0],str(Q(c))]],dict(constant=str(c)))
 z=const(0);one=const(1)
 def bd(i):return N[i]['bound']
 def size(i):return bd(i)['support']if bd(i)else 0
 def dl(i):return int(N[i]['QQ_denominator']),int(N[i]['integer_l1_upper'])
 def mul(i,j):
  if i==z or j==z:return z
  if i==one:return j
  if j==one:return i
  A,B=bd(i),bd(j);d,n=dl(i);e,m=dl(j);box=[[A['box'][k][0]+B['box'][k][0],A['box'][k][1]+B['box'][k][1]]for k in range(3)];lo=A['degree'][0]+B['degree'][0];hi=A['degree'][1]+B['degree'][1];b=dict(support=min(A['support']*B['support'],lattice(tuple(map(tuple,box)),lo,hi),n*m),box=box,degree=[lo,hi]);return put('mul',[i,j],b,d*e,n*m)
 def add(i,j,sgn=1):
  if j==z:return i
  if i==z:return j if sgn==1 else mul(const(-1),j)
  if sgn<0 and i==j:return z
  A,B=bd(i),bd(j);d,n=dl(i);e,m=dl(j);D=math.lcm(d,e);L=n*(D//d)+m*(D//e);box=[[min(A['box'][k][0],B['box'][k][0]),max(A['box'][k][1],B['box'][k][1])]for k in range(3)];lo=min(A['degree'][0],B['degree'][0]);hi=max(A['degree'][1],B['degree'][1]);b=dict(support=min(A['support']+B['support'],lattice(tuple(map(tuple,box)),lo,hi),L),box=box,degree=[lo,hi]);return put('add'if sgn>0 else'sub',[i,j],b,D,L)
 Bnode=atom(BP,dict(named_FULLU_unit='B'))
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
   else:d=max(b,c);v=(add(timesB(x,d-b),timesB(y,d-c),-1 if op=='s'else 1),d)
  else:
   num,den=cp['matrix25x32'][a[1]][a[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   if len(den)==1:assert den[0][0][1:]==[0,0]and Q(den[0][1]);unit=Q(den[0][1]);ashift=-den[0][0][0];bp=0
   else:
    dd={tuple(m):Q(c)for m,c in den};bb={tuple(m):Q(c)for m,c in BP};assert dd.keys()==bb.keys();unit=dd[next(iter(bb))]/bb[next(iter(bb))];assert unit and all(dd[m]==unit*bb[m]for m in bb);ashift=0;bp=1
   inventory.append(dict(old_DAG_root=i,definition=a,known_B_denominator_power=bp,QQ_unit=str(unit),a_shift=ashift));v=(atom(num,dict(path='accepted_Q2L_source_atom',old_DAG_root=i,definition=a),unit,ashift),bp)
  cache[i]=v;return v
 rows=[2,4,5,8];col=q['new22_target_original32_columns'].index(14);images=[imp(q['remaining15x22_scalar_roots'][r][col])for r in rows];bp=max(p for i,p in images);vec=[timesB(i,bp-p)for i,p in images]
 cofs=[atom(load(ke,f'root{i}.json')['polynomial'],dict(path=f'reports/{ke}/root{i}.json',key='polynomial'))for i in kb['cofactor_roots']];anchor=[atom(bot[k][0],dict(path='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json',key='bottom4x6',index=[k,0]))for k in range(4)]
 k14=z;terms=[]
 for k in range(4):t=mul(cofs[k],vec[k]);terms.append(t);k14=add(k14,t)
 cross=[add(mul(anchor[0],vec[k]),mul(anchor[k],vec[0]),-1)for k in range(1,4)]
 def profile(name,outputs,saved):
  need=set();todo=outputs[:]
  while todo:
   i=todo.pop()
   if i in need:continue
   need.add(i);todo+=N[i]['args']
  use={i:0 for i in need}
  for i in need:
   for j in N[i]['args']:use[j]+=1
  for i in outputs:use[i]+=1
  live=peak=0;records=[]
  for i in sorted(need):
   live+=size(i);scratch=size(i)if N[i]['op']in['mul','add','sub']else 0;peak=max(peak,live+scratch);records.append(dict(node=i,live_upper=live,scratch=scratch))
   for j in N[i]['args']:
    use[j]-=1
    if not use[j]:live-=size(j)
  atoms=[i for i in need if N[i]['op']=='atom'];literal=sum(size(i)for i in atoms);total=peak+literal+15000;artifact=1000000
  for i in set(saved)|set(atoms):
   d,n=dl(i);artifact+=size(i)*(len(str(n))+len(str(d))+115)
  return dict(name=name,output_roots=outputs,future_saved_roots=saved,needed_nodes=sorted(need),needed_count=len(need),output_bounds=[bd(i)for i in outputs],output_numerator_l1_bits=[dl(i)[1].bit_length()for i in outputs],output_QQ_denominator_bits=[dl(i)[0].bit_length()for i in outputs],literal_atom_reserve=literal,lastuse_peak_with_scratch=peak,total_live_upper=total,raw_artifact_bytes_upper=artifact,resource_gate_passes=total<=150000 and artifact<32*1024*1024,lastuse_records=records)
 profiles=[profile('kappa14_only',[k14],[k14]+vec),profile('three_raw_cross_minors_only',cross,cross+vec),profile('kappa14_and_rawcross3',[k14]+cross,[k14]+cross+vec)]
 op=dict(nodes=N,current15_source_rows=rows,original32_target=14,anchor_original32=22,common_B_clearing_power=bp,column14_roots=vec,kappa4_roots=cofs,kappa14_output=k14,raw_cross3_outputs=cross,
  raw_cross_formulas=['M[2,22]*M[i,14]-M[i,22]*M[2,14] for i=4,5,8'],atom_denominator_inventory=inventory,all14_source_values_unevaluated=True,all8_Q9_tail_effects_retained=True)
 dump(out/'operator-and-cost.json',op);dump(out/'profiles.json',profiles)
 proof=dict(logical_distinction_preserved=True,zero_kappa14_means='Onrank3A=(Hbar27,Hbar31,M22), vectorM14 liesinspanA. Equivalentlydet[A,M14]=0; notrawproportionality.',
  normalized_chart='OnexistingD(Saux*T), kappa3isunit so A first3rows invertible; coefficientsalpha,beta,gamma have finite3adjugate recipe, no newchart.',
  zero_kappa14_permitted_column_change='Replace column14 bycolumn14-alpha*Hbar27-beta*Hbar31, retaining these heavycorrections onALLcurrent15rows. Pairfour-vector then gamma*M22. DoNOT useuncorrectedbare7border.',
  nonzero_kappa14_means='M14isindependentofA generically; no newrankbeyondalreadyknownA4sinceextra2wasused. Givesanalternatepairbasis onlyonitsownchart.',
  raw_cross3_zero_means='SinceM[2,22]isasavednonzeropolynomial, the threeidentities proveM14 israwproportional toM22 inthefractionfield, andallcross-product identitiesdescendtothedomain. Thisstrongergatepermitsrawpairreplacement ontheextendedsevenpurecolumns.',
  raw_cross3_nonzero_means='Rawproportionalityisfalse; doesnotdecidekappa14zero/nonzero. Heavycorrectionsmustremain.',
  no_new_scalar_computation=True,previous_raw10_support_failure_preserved=True,no_global_rank_or_complement_claim=True)
 dump(out/'scope-proof.json',proof);dump(out/'phase.json',dict(phase='FINAL source-only; no kappa14 or rawcross scalar job launched'));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,logical_kappa_vs_raw_proportionality_gate_explicit=True,typed_nodes=len(N),common_B_clear_power=bp,profiles=[{k:v for k,v in p.items()if k not in['lastuse_records','needed_nodes']}for p in profiles],no_scalar_value_rank_point_prime=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
