#!/usr/bin/env python3
"""SOURCE ONLY fixedK92/18 det11: rankoneH3 Laplace weights and denominator/support bounds."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_T18_DET11_SOURCE_PREFLIGHT_20261104_v63'
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
 dump(out/'phase.json',dict(phase='source-only fixedK92target18 det11 profile; no scalar work'))
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
 assert load('S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63','report.json')['rawpair_fixed8_rank_le1']
 rk=load('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63');cols=q['new22_target_original32_columns'];pure=[22,11,13,25,15,19,14,18];top=[]
 for r in[0,1,3,6]:
  rr=[]
  for j in pure:
   if r==0:x=(atom(h9['source_h9'][h9['original_target_labels'].index(j)],dict(path='reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json',h9_target=j)),0)
   else:x=imp(q['remaining15x22_scalar_roots'][r][cols.index(j)])
   rr.append(x)
  top.append(rr)
 pair=[]
 for j in pure:
  if j in[14,18]:st='S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63'if j==14 else'S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63';x=atom(load(st,'root21.json')['polynomial'],dict(path=f'reports/{st}/root21.json',key='polynomial'))
  else:k=[22,11,13,25,15,19].index(j);x=atom(bot[0][k],dict(path='reports/S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63/native.json',key='bottom4x6',index=[0,k]))
  pair.append((x,0))
 top.append(pair);rho=hc['known_row_units_rho14'];c={i:Q(rho[i]['constant'])for i in[6,8,9,10,11,12]};assert rho[11]==rho[9];lower=[]
 for kind,r,ci in[('A',11,10),('B',12,11),('D',9,8),('C',7,6),('E',13,12),('F',10,9)]:
  rr=[]
  for j in pure:
   x,bp=imp(q['remaining15x22_scalar_roots'][r][cols.index(j)])
   if kind=='B':y,cpower=imp(q['remaining15x22_scalar_roots'][10][cols.index(j)]);d=max(bp,cpower);x=add(timesB(x,d-bp),timesB(y,d-cpower),-1);bp=d
   x=mul(x,const(1/c[ci]))
   if kind=='E':x=shift(x,1)
   rr.append((x,bp))
  lower.append(rr)
 rows=top+lower;Bclear=[max(b for r in rows for x,b in[r[j]])for j in range(8)];rows=[[timesB(x,Bclear[j]-b)for j,(x,b)in enumerate(r)]for r in rows];rowa=[min(bd(i)['box'][0][0]for i in r if i!=z)for r in rows];rows=[[shift(i,-a)for i in r]for r,a in zip(rows,rowa)]
 assert Bclear==[0]*8
 # Detached majorants; do not store/evaluate determinant coefficient nodes.
 def value(i):return None if bd(i)is None else dict(bd(i),D=dl(i)[0],L=dl(i)[1])
 def bm(a,b):
  if a is None or b is None:return None
  box=[[a['box'][j][0]+b['box'][j][0],a['box'][j][1]+b['box'][j][1]]for j in range(3)];lo=a['degree'][0]+b['degree'][0];hi=a['degree'][1]+b['degree'][1];return dict(support=min(a['support']*b['support'],lattice(tuple(map(tuple,box)),lo,hi),a['L']*b['L']),box=box,degree=[lo,hi],D=a['D']*b['D'],L=a['L']*b['L'])
 def ba(a,b):
  if a is None:return b
  if b is None:return a
  box=[[min(a['box'][j][0],b['box'][j][0]),max(a['box'][j][1],b['box'][j][1])]for j in range(3)];lo=min(a['degree'][0],b['degree'][0]);hi=max(a['degree'][1],b['degree'][1]);D=math.lcm(a['D'],b['D']);L=a['L']*(D//a['D'])+b['L']*(D//b['D']);return dict(support=min(a['support']+b['support'],lattice(tuple(map(tuple,box)),lo,hi),L),box=box,degree=[lo,hi],D=D,L=L)
 def bdet(m):
  state={0:value(one)}
  for r in m:
   nxt={}
   for mask,a in state.items():
    for j,b in enumerate(r):
     if not(mask>>j&1)and b:nxt[mask|1<<j]=ba(nxt.get(mask|1<<j),bm(a,b))
   state=nxt
  return state[(1<<len(m))-1]
 def va(poly):
  if not poly:return None
  D=math.lcm(*(Q(c).denominator for m,c in poly));L=sum(abs(Q(c).numerator)*(D//Q(c).denominator)for m,c in poly);return dict(support=len(poly),box=[[min(m[j]for m,c in poly),max(m[j]for m,c in poly)]for j in range(3)],degree=[min(sum(m)for m,c in poly),max(sum(m)for m,c in poly)],D=D,L=L)
 g=va(hc['raw_compact14x5'][6][1]);h=va(hc['raw_compact14x5'][10][1]);dif={tuple(m):Q(c)for m,c in hc['raw_compact14x5'][11][1]}
 for m,v in hc['raw_compact14x5'][9][1]:dif[tuple(m)]=dif.get(tuple(m),Q(0))-Q(v)
 hp=va([[list(m),str(v)]for m,v in dif.items()if v]);Delta=va(rk['nonzero_selected_minor_Delta']);l1,l2,l3=[va(p)for p in rk['cleared_upper_row_coordinates7xr'][1]]
 # g timesH3 minors: ten exact surviving coefficient profiles.
 data=[((0,2,4),-1,h,l3,'-h*lambda428'),((0,2,5),1,h,l2,'h*lambda419'),((0,4,5),1,h,Delta,'h*Delta'),((1,2,4),-1,hp,l3,'-h120*lambda428'),((1,2,5),1,hp,l2,'h120*lambda419'),((1,4,5),1,hp,Delta,'h120*Delta'),((2,3,4),1,g,l3,'g*lambda428'),((2,3,5),-1,g,l2,'-g*lambda419'),((2,4,5),1,g,l1,'g*lambda505'),((3,4,5),1,g,Delta,'g*Delta')]
 terms=[];answer=None
 for selected,hsign,a,b,label in data:
  remaining=[i for i in range(6)if i not in selected];m=[[value(j)for j in r]for r in rows[:5]+[rows[5+i]for i in remaining]];p8=bdet(m);w=bm(a,b);term=bm(w,p8);answer=ba(answer,term);terms.append(dict(H3_selected_lower6=list(selected),pure_remaining_lower6=remaining,H3_weight=label,H3_weight_sign=hsign,outer_sign=-1 if sum(selected)%2 else 1,pure8_bound={k:v for k,v in p8.items()if k not in['D','L']},weight_bound={k:v for k,v in w.items()if k not in['D','L']}))
 latticeout=lattice(tuple(map(tuple,answer['box'])),*answer['degree']);matrixsupport=sum(support(i)for r in rows for i in r);integerbytes=latticeout*(len(str(answer['L']))+45)+1000
 cost=dict(proposed_numerator='g*det11',pure8_column_known_B_clears=Bclear,row_known_a_factors=rowa,source_matrix11x8_pure_support_bound=matrixsupport,
  ten_nonzero_H3_weighted_terms=terms,ten_other_H3_minors_zero_from_three_proportional_rows=True,output_bound={k:v for k,v in answer.items()if k not in['D','L']},dense_output_slots=latticeout,global_integer_l1_bits=answer['L'].bit_length(),global_QQ_denominator_bits=answer['D'].bit_length(),integer_output_byte_bound=integerbytes,
  dense_output_alone_exceeds150k=latticeout>150000,NO_actual_det11_or_weight_or_polynomial_products=True,whole_worker_gate_NOT_passed=True,inner_lifetime_NOT_evaluated_if_output_alone_fails=True)
 dump(out/'cost.json',cost)
 proof=dict(source_rows='P4[0,1,3,6],G20current2,A=current11/c10,B=(current12-current10)/c11,D=current9/c8,C=current7/c6,E=a*current13/c12,F=current10/c9',
  source_pure8_columns=pure,H3_columns_normalized_by_known_L=[28,30,29],all_three_B3_cofactors_retained=True,factor_a_lambda419_retained_inside_E=True,
  source_H3_row_order=['N10','N11-N9','N8','N6','N12','N9'],only_A_B_C_proportional=True,no_other_H3_vanishing_assumed=True,
  raw8_proportionality_primary_bound=True,conditional_formula_now_primary='actualnextentry=c8*(g*det11)/(Delta407*a²*N7)',
  no_g_division_in_polynomial_recipe=True,all_existing_unknown_chart_complements_and_tail2_retained=True,
  next_method_required='Outputlatticefailuremustbeaddressedbyfactor/sourcecompressionorblockedcoefficientstreamingwithprovenRAWbound,notcapraise.'if latticeout>150000 else'Needcompletephase/height/RAWgatebeforeanyscalar.')
 dump(out/'source-proof.json',proof);dump(out/'phase.json',dict(phase='FINAL source-only det11 bound; no scalar worker'));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,fixed_K92_target18=True,source_metadata_nodes=len(N),pure11x8_input_support_upper=matrixsupport,weighted_H3_terms10=True,dense_output_slots=latticeout,output_degree=answer['degree'],integer_l1_bits=answer['L'].bit_length(),QQ_denominator_bits=answer['D'].bit_length(),integer_output_byte_bound=integerbytes,resource_gate_passes=False,actual_coefficient_job_NOT_launched=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
