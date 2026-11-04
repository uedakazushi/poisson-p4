#!/usr/bin/env python3
"""Source-only full-worker bounds for ONE G115/target19 bordered6 recipe."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
BP=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
@lru_cache(None)
def lattice(box,lo,hi):
 return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 dump(out/'phase.json',dict(phase='source-only same G115 target19; no scalar coefficient work'))
 def load(st,name='certificate.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return read(p)
 st443=load('S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63','report.json');assert st443['all_passed']
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');h9=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');fx=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');fixed={r['old_DAG_root']:r['Laurent_polynomial']for r in fx['node_values']};strip=load('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json');bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];assert load('S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW','report.json')['all_passed']
 N=[];keys={};cache={};inventory=[]
 def put(op,args,bound,D,L,**kw):i=len(N);N.append(dict(id=i,op=op,args=args,bound=bound,QQ_denominator=str(D),integer_l1_upper=str(L),**kw));return i
 def atom(poly,ref,unit=Q(1),shift=0):
  p=[([m[0]+shift,m[1],m[2]],Q(c)/unit)for m,c in poly if Q(c)];key=tuple((tuple(m),str(c))for m,c in p)
  if key in keys:return keys[key]
  D=math.lcm(*(c.denominator for m,c in p))if p else 1;L=sum(abs(c.numerator)*(D//c.denominator)for m,c in p)
  b=None if not p else dict(support=len(p),box=[[min(m[k]for m,c in p),max(m[k]for m,c in p)]for k in range(3)],degree=[min(sum(m)for m,c in p),max(sum(m)for m,c in p)])
  i=put('atom',[],b,D,L,reference=ref,QQ_unit=str(unit),a_shift=shift);keys[key]=i;return i
 def const(c):return atom([[[0,0,0],str(Q(c))]],dict(constant=str(c)))
 z=const(0);one=const(1)
 def bd(i):return N[i]['bound']
 def dl(i):return int(N[i]['QQ_denominator']),int(N[i]['integer_l1_upper'])
 def support(i):return bd(i)['support']if bd(i)else 0
 def makebd(A,B,additive):
  if additive:box=[[min(A['box'][k][0],B['box'][k][0]),max(A['box'][k][1],B['box'][k][1])]for k in range(3)];lo=min(A['degree'][0],B['degree'][0]);hi=max(A['degree'][1],B['degree'][1]);s=A['support']+B['support']
  else:box=[[A['box'][k][0]+B['box'][k][0],A['box'][k][1]+B['box'][k][1]]for k in range(3)];lo=A['degree'][0]+B['degree'][0];hi=A['degree'][1]+B['degree'][1];s=A['support']*B['support']
  return dict(support=min(s,lattice(tuple(map(tuple,box)),lo,hi)),box=box,degree=[lo,hi])
 def mul(i,j):
  if i==z or j==z:return z
  if i==one:return j
  if j==one:return i
  d,n=dl(i);e,m=dl(j);b=makebd(bd(i),bd(j),False);b['support']=min(b['support'],n*m);return put('mul',[i,j],b,d*e,n*m)
 def add(i,j,sgn=1):
  if j==z:return i
  if i==z:return j if sgn==1 else mul(const(-1),j)
  if sgn<0 and i==j:return z
  d,n=dl(i);e,m=dl(j);D=math.lcm(d,e);L=n*(D//d)+m*(D//e);b=makebd(bd(i),bd(j),True);b['support']=min(b['support'],L);return put('add'if sgn>0 else'sub',[i,j],b,D,L)
 def shift(i,k):
  if i==z or not k:return i
  A=bd(i);b=dict(support=A['support'],box=[[x+k for x in A['box'][0]]]+A['box'][1:],degree=[x+k for x in A['degree']]);d,n=dl(i);return put('known_a_shift',[i],b,d,n,a_exponent_shift=k)
 Bnode=atom(BP,dict(named_FULLU_unit='B'))
 def timesB(i,k):
  for _ in range(k):i=mul(i,Bnode)
  return i
 def imp(i):
  if i in cache:return cache[i]
  a=q['scalar_DAG']['nodes'][i];op=a[0]
  if i in fixed:v=(atom(fixed[i],dict(path='reports/S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63/native.json',old_DAG_root=i)),0)
  elif op=='z':v=(z,0)
  elif op=='one':v=(one,0)
  elif op in['a','s','m']:
   x,b=imp(a[1]);y,c=imp(a[2])
   if op=='m':v=(mul(x,y),b+c)
   else:d=max(b,c);v=(add(timesB(x,d-b),timesB(y,d-c),-1 if op=='s'else 1),d)
  else:
   num,den=cp['matrix25x32'][a[1]][a[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   if len(den)==1:assert den[0][0][1:]==[0,0]and Q(den[0][1]);unit=Q(den[0][1]);ashift=-den[0][0][0];bpower=0
   else:
    dd={tuple(m):Q(c)for m,c in den};bb={tuple(m):Q(c)for m,c in BP};assert dd.keys()==bb.keys();unit=dd[next(iter(bb))]/bb[next(iter(bb))];assert unit and all(dd[m]==unit*bb[m]for m in bb);ashift=0;bpower=1
   inventory.append(dict(old_DAG_root=i,definition=a,known_B_denominator_power=bpower,QQ_unit=str(unit),a_shift=ashift));v=(atom(num,dict(path='accepted_Q2L_source_atom',old_DAG_root=i,definition=a,known_B_denominator_power=bpower),unit,ashift),bpower)
  cache[i]=v;return v
 cols=q['new22_target_original32_columns'];labels=[22,11,13,25,15,19];matrix=[];rowBp=[]
 for r in[0,1,3,6]:
  rr=[]
  for j in labels:
   if r==0:x=(atom(h9['source_h9'][h9['original_target_labels'].index(j)],dict(path='reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json',h9_target=j)),0)
   else:x=imp(q['remaining15x22_scalar_roots'][r][cols.index(j)])
   rr.append(x)
  bp=max(p for i,p in rr);matrix.append([timesB(i,bp-p)for i,p in rr]);rowBp.append(bp)
 matrix.append([atom(bot[0][j],dict(path='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json',key='bottom4x6',index=[0,j]))for j in range(6)]);rowBp.append(0)
 g=atom(hc['raw_compact14x5'][6][1],dict(path='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json',key='raw_compact14x5',index=[6,1]));hh=atom(hc['raw_compact14x5'][10][1],dict(path='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json',key='raw_compact14x5',index=[10,1]));c6=Q(hc['known_row_units_rho14'][6]['constant']);c10=Q(hc['known_row_units_rho14'][10]['constant']);Kr=[]
 for j in labels:
  x,b=imp(q['remaining15x22_scalar_roots'][11][cols.index(j)]);y,c=imp(q['remaining15x22_scalar_roots'][7][cols.index(j)]);bp=max(b,c);Kr.append((add(mul(mul(g,const(1/c10)),timesB(x,bp-b)),mul(mul(hh,const(1/c6)),timesB(y,bp-c)),-1),bp))
 bp=max(p for i,p in Kr);matrix.append([timesB(i,bp-p)for i,p in Kr]);rowBp.append(bp)
 rowa=[min(bd(i)['box'][0][0]for i in rr if i!=z)for rr in matrix];matrix=[[shift(i,-a)for i in rr]for rr,a in zip(matrix,rowa)]
 def det3(m):
  ans=z
  for p in itertools.permutations(range(3)):
   t=one
   for i in range(3):t=mul(t,m[i][p[i]])
   ans=add(ans,t,-1 if sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))%2 else 1)
  return ans
 def det4(chosen):
  if tuple(chosen)==(0,1,2,3):return atom(strip['exact_unscaled_pure4_minor'],dict(path='reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63/native.json',key='exact_unscaled_pure4_minor'),shift=-sum(rowa[:4]))
  # Expand along sparse Q9 row, rather than recomputing all24 products.
  ans=z
  for c in range(4):
   mm=[[matrix[i][chosen[j]]for j in range(4)if j!=c]for i in[1,2,3]];t=mul(matrix[0][chosen[c]],det3(mm));ans=add(ans,t,-1 if c%2 else 1)
  return ans
 events=[];output=z
 for chosen in itertools.combinations(range(6),4):
  rem=[i for i in range(6)if i not in chosen];a=det4(chosen);b=add(mul(matrix[4][rem[0]],matrix[5][rem[1]]),mul(matrix[4][rem[1]],matrix[5][rem[0]]),-1);events.append(dict(top4_columns=list(chosen),top4_root=a,bottom2_root=b,sign=-1 if sum(chosen)%2 else 1));output=add(output,mul(a,b),events[-1]['sign'])
 roots=[x for e in events for x in[e['top4_root'],e['bottom2_root']]];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);todo+=N[i]['args']
 use={i:0 for i in need}
 for i in need:
  for j in N[i]['args']:use[j]+=1
 schedule={}
 for e in events:
  a,b=e['top4_root'],e['bottom2_root'];use[a]+=1;use[b]+=1;schedule.setdefault(max(a,b),[]).append(e)
 live=peak=0;records=[]
 for i in sorted(need):
  live+=support(i);scratch=support(i)if N[i]['op']in['mul','add','sub']else 0;peak=max(peak,live+scratch);records.append(dict(node=i,live=live,scratch=scratch,events_after=[e['top4_columns']for e in schedule.get(i,[])]))
  for j in N[i]['args']:
   use[j]-=1
   if not use[j]:live-=support(j)
  for e in schedule.get(i,[]):
   for j in[e['top4_root'],e['bottom2_root']]:
    use[j]-=1
    if not use[j]:live-=support(j)
 atoms=[i for i in need if N[i]['op']=='atom'];literal=sum(support(i)for i in atoms);ob=bd(output);D,L=dl(output);dense=lattice(tuple(map(tuple,ob['box'])),*ob['degree']);total=peak+literal+dense+15000
 scalarpairs=sum(support(e['top4_root'])*support(e['bottom2_root'])for e in events);maxL=max(L,*[dl(i)[1]for i in atoms]);maxD=max(D,*[dl(i)[0]for i in atoms]);artifact=(dense+literal)*(len(str(maxL))+len(str(maxD))+115)+2000000
 assert all(D%(dl(e['top4_root'])[0]*dl(e['bottom2_root'])[0])==0 for e in events)
 op=dict(nodes=N,matrix6x6=matrix,output=output,laplace15=events,needed_fused_nodes=sorted(need),row_a_factors=rowa,row_B_clearings=rowBp,original32_columns=labels,source_rows='P4(Q9,G18,G39,G68),literalG20,Ksmall(g*G115/c10-h*G85/c6)',actual_input_roots_current15=[0,1,3,6,2,11,7],literal_pair_current15=2,atom_B_denominator_inventory=inventory,new_parameter_products_NOT_evaluated=True)
 dump(out/'operator-and-cost.json',op)
 proof=dict(same_one_source_row_and_target=True,chosen_remaining_current15=11,actual135=115,target_original32=19,Ksmall_middleheavy_zero_from_published_Cramer=True,known_pair_pure6_rank1=True,
  D5_literal='det5(P4,literalG20;columns[22,11,13,25,15])',D6_literal='det6(P4,literalG20,Ksmall;columns[22,11,13,25,15,19])',residual_formula='c10*D6_literal/(g*D5_literal)',
  field_identity='det6(P4,Bhat4,Ksmall)*D5_literal=Saux*D6_literal by pair-pure6 rank1 and determinant multilinearity',
  explicit_extra_chart='D(g*D5_literal)',extra_chart_complement_retained=True,D5_literal_nonzero_generic='G20 anchor22 is a saved nonzero polynomial; all pair pure6 proportional; Saux nonzero forces its5minor nonzero in fractionfield',
  no_delta2_p_R_h22_or_kappa_extra_inverse=True,full57_lift='S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63',all_remaining7_tail_effects_retained=True,no_scalar_D5_D6_or_rank_computation=True)
 dump(out/'proof.json',proof)
 cost=dict(output_bound=ob,output_integer_l1_bits=L.bit_length(),output_QQ_denominator_bits=D.bit_length(),matrix_input_support_upper=sum(support(i)for rr in matrix for i in rr),literal_atom_reserve=literal,lastuse_source_peak_with_scratch=peak,dense_output_lattice_slots=dense,work_reserve=15000,total_whole_worker_coefficient_slots_upper=total,raw_artifact_bytes_upper=artifact,
  scalar_convolution_pairs_upper=scalarpairs,runtime60s_not_guaranteed=True,all15_common_QQ_scales_integer=True,resource_gate_passes=total<=150000 and artifact<32*1024*1024,lastuse_records=records,all_B_and_a_normalization_accounted=True)
 dump(out/'cost.json',cost);note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,typed_scalar_nodes=len(N),needed_fused_nodes=len(need),literal_atom_terms=literal,output_support_upper=ob['support'],output_degree_range=ob['degree'],output_box=ob['box'],total_live_upper=total,raw_artifact_bytes_upper=artifact,output_integer_l1_bits=L.bit_length(),output_QQ_denominator_bits=D.bit_length(),resource_gate_passes=cost['resource_gate_passes'],scalar_convolution_pairs_upper=scalarpairs,no_actual_value_minor_rank_point_prime=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);dump(out/'phase.json',dict(phase='FINAL source-only; scalar job NOT launched'));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
