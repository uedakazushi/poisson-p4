#!/usr/bin/env python3
"""Source-only integer recipe for ONE K120/target14 numerator after symbolic g cancellation."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K120_T14_CANCELLED_SOURCE_PREFLIGHT_20261104_v63'
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
 dump(out/'phase.json',dict(phase='source-only K120 target14 with three retained multilinear terms; no scalar work'))
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
 # Fixed source7 proportionality gate; no substitution of unknown columns.
 assert load('S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63','report.json')['raw_proportionality_on_seven_columns_proved']
 cols=q['new22_target_original32_columns'];labels=[22,11,13,25,15,19,14];top=[];bottom=[];allrows=[]
 for r in[0,1,3,6]:
  rr=[]
  for j in labels:
   if r==0:x=(atom(h9['source_h9'][h9['original_target_labels'].index(j)],dict(path='reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json',h9_target=j)),0)
   else:x=imp(q['remaining15x22_scalar_roots'][r][cols.index(j)])
   rr.append(x)
  top.append(rr)
 pair=[]
 for j in labels:
  if j==14:x=atom(load('S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63','root21.json')['polynomial'],dict(path='reports/S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63/root21.json',key='polynomial'))
  else:k=[22,11,13,25,15,19].index(j);x=atom(bot[0][k],dict(path='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json',key='bottom4x6',index=[0,k]))
  pair.append((x,0))
 top.append(pair)
 c6=Q(hc['known_row_units_rho14'][6]['constant']);c10=Q(hc['known_row_units_rho14'][10]['constant']);c11=Q(hc['known_row_units_rho14'][11]['constant']);assert hc['known_row_units_rho14'][11]==hc['known_row_units_rho14'][9]
 for kind in['A','B','C']:
  rr=[]
  for j in labels:
   r=11 if kind=='A'else 12 if kind=='B'else 7;x,b=imp(q['remaining15x22_scalar_roots'][r][cols.index(j)])
   if kind=='B':
    y,c=imp(q['remaining15x22_scalar_roots'][10][cols.index(j)]);d=max(b,c);x=add(timesB(x,d-b),timesB(y,d-c),-1);b=d
   x=mul(x,const(1/(c10 if kind=='A'else c11 if kind=='B'else c6)));rr.append((x,b))
  bottom.append(rr)
 allrows=top+bottom;columnBp=[max(bp for rr in allrows for i,bp in[rr[j]])for j in range(7)];allrows=[[timesB(i,columnBp[j]-bp)for j,(i,bp)in enumerate(rr)]for rr in allrows];top=allrows[:5];A,Bb,C=allrows[5:]
 g=atom(hc['raw_compact14x5'][6][1],dict(path='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json',key='raw_compact14x5',index=[6,1]));h=atom(hc['raw_compact14x5'][10][1],dict(path='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json',key='raw_compact14x5',index=[10,1]));dif={tuple(m):Q(c)for m,c in hc['raw_compact14x5'][11][1]}
 for m,c in hc['raw_compact14x5'][9][1]:dif[tuple(m)]=dif.get(tuple(m),Q(0))-Q(c)
 hpPoly=[[list(m),str(c)]for m,c in sorted(dif.items())if c];assert len(hpPoly)==35;hp=atom(hpPoly,dict(derived_existing_coefficient_difference=True,source='N11,1-N9,1',polynomial=hpPoly))
 # Remove only the five definite known a-row factors; bottom weights unchanged.
 rowa=[min(bd(i)['box'][0][0]for i in rr if i!=z)for rr in top];top=[[shift(i,-a)for i in rr]for rr,a in zip(top,rowa)]
 def minor(m):
  k=len(m)
  if k==1:return m[0][0]
  ans=z
  for j in range(k):ans=add(ans,mul(m[0][j],minor([[m[i][c]for c in range(k)if c!=j]for i in range(1,k)])),-1 if j%2 else 1)
  return ans
 def top5(chosen):
  ans=z
  for j in range(5):
   cc=[chosen[k]for k in range(5)if k!=j]
   if cc==[0,1,2,3]:d4=atom(strip['exact_unscaled_pure4_minor'],dict(path='reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63/native.json',key='exact_unscaled_pure4_minor'),shift=-sum(rowa[:4]))
   else:d4=minor([[top[i][c]for c in cc]for i in range(4)])
   ans=add(ans,mul(top[4][chosen[j]],d4),-1 if j%2 else 1)
  return ans
 def det2(U,V,j,k):return add(mul(U[j],V[k]),mul(U[k],V[j]),-1)
 events=[];output=z
 for rem in itertools.combinations(range(7),2):
  chosen=[j for j in range(7)if j not in rem];tp=top5(chosen);j,k=rem
  w=add(add(mul(g,det2(A,Bb,j,k)),mul(hp,det2(A,C,j,k)),-1),mul(h,det2(C,Bb,j,k)),-1)
  # Expand along finaltwo rows: sum(rowindices5,6)=11.
  sgn=-1 if(11+j+k)%2 else 1;events.append(dict(top5_columns=chosen,bottom2_columns=list(rem),top5_root=tp,weighted_bottom2_root=w,sign=sgn));output=add(output,mul(tp,w),sgn)
 roots=[r for e in events for r in[e['top5_root'],e['weighted_bottom2_root']]];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);todo+=N[i]['args']
 use={i:0 for i in need}
 for i in need:
  for j in N[i]['args']:use[j]+=1
 schedule={}
 for e in events:
  a,b=e['top5_root'],e['weighted_bottom2_root'];use[a]+=1;use[b]+=1;schedule.setdefault(max(a,b),[]).append(e)
 live=peak=0;records=[]
 for i in sorted(need):
  live+=support(i);scratch=support(i)if N[i]['op']in['mul','add','sub']else 0;peak=max(peak,live+scratch);records.append(dict(node=i,live=live,scratch=scratch,events_after=[e['bottom2_columns']for e in schedule.get(i,[])]))
  for j in N[i]['args']:
   use[j]-=1
   if not use[j]:live-=support(j)
  for e in schedule.get(i,[]):
   for j in[e['top5_root'],e['weighted_bottom2_root']]:
    use[j]-=1
    if not use[j]:live-=support(j)
 atoms=[i for i in need if N[i]['op']=='atom'];literal=sum(support(i)for i in atoms);ob=bd(output);D,L=dl(output);dense=lattice(tuple(map(tuple,ob['box'])),*ob['degree']);total=peak+literal+dense+15000;scalarpairs=sum(support(e['top5_root'])*support(e['weighted_bottom2_root'])for e in events);maxL=max(L,*[dl(i)[1]for i in atoms]);maxD=max(D,*[dl(i)[0]for i in atoms]);artifact=(dense+literal)*(len(str(maxL))+len(str(maxD))+115)+3000000
 assert all(D%(dl(e['top5_root'])[0]*dl(e['weighted_bottom2_root'])[0])==0 for e in events)
 op=dict(nodes=N,top5_matrix=top,bottomABC_matrix=[A,Bb,C],weights_g_h120_h=[g,hp,h],output=output,laplace21=events,needed_fused_nodes=sorted(need),known_a_top5_row_factors=rowa,known_B_column_clearings=columnBp,source_pure_columns_original32=labels,
  exact_universal_cancellation='D7/g = g det(P4,G20,A,B) - h120 det(P4,G20,A,C) - h det(P4,G20,C,B)',
  weighted_bottom2='g*(A_j B_k-A_k B_j)-h120*(A_j C_k-A_k C_j)-h*(C_j B_k-C_k B_j)',all_three_terms_retained=True,no_numeric_g_division=True,new_parameter_products_NOT_evaluated=True)
 dump(out/'operator-and-cost.json',op)
 # Tiny independent exterior-algebra check of the cancellation and all signs.
 wedge={('A','B'):1,('A','C'):1,('C','B'):1}
 cancellation=dict(formal_coefficient_of_AB='g^2',formal_coefficient_of_AC='-g*h120',formal_coefficient_of_CB='-g*h',duplicate_CC_coefficient_cancels=True,divide_by_g_is_symbolic_multilinearity_NOT_parameter_division=True,
  normalized_output_relation=f"raw(D7/g)=a^{sum(rowa)}*B^(-{sum(columnBp)})*output",next_residual='c11*(raw D7/g)/(a^2*N6)',all7_raw_pair_columns_rank1_gate_bound=True,tail2_effects_retained_in_prior_operator=True,
  original_current15_mnemonics_not_bare135=True,no_next_scalar_job=True)
 dump(out/'source-proof.json',cancellation)
 cost=dict(output_bound=ob,output_integer_l1_bits=L.bit_length(),output_QQ_denominator_bits=D.bit_length(),matrix_input_support_upper=sum(support(i)for rr in allrows for i in rr),literal_atom_reserve=literal,lastuse_source_peak_with_scratch=peak,dense_output_lattice_slots=dense,work_reserve=15000,total_whole_worker_slots_upper=total,
  scalar_convolution_pairs_upper=scalarpairs,raw_artifact_bytes_upper=artifact,all21_integer_denominator_scales_verified=True,lastuse_records=records,
  resource_gate_passes=total<=150000 and artifact<32*1024*1024,runtime60s_not_guaranteed=True)
 dump(out/'cost.json',cost);dump(out/'phase.json',dict(phase='FINAL source-only; no N7 coefficient calculation'));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,one_fixed_K120_target14=True,three_multilinear_terms_retained=True,no_numeric_g_division=True,typed_nodes=len(N),needed_nodes=len(need),source_atom_terms=literal,output_support_upper=ob['support'],output_degree_range=ob['degree'],total_live_upper=total,raw_artifact_bytes_upper=artifact,scalar_convolution_pairs_upper=scalarpairs,output_integer_l1_bits=L.bit_length(),output_QQ_denominator_bits=D.bit_length(),resource_gate_passes=cost['resource_gate_passes'],no_scalar_value_minor_rank_point_prime=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
