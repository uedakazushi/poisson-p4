#!/usr/bin/env python3
"""SOURCE ONLY weightedNewton and GaussQQcontent for SAME K92/18 ten-term numerator."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_T18_NEWTON_CONTENT_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
BP=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
WEIGHTS=[(1,1,1),(2,1,0),(2,1,1),(1,1,2),(2,1,2),(1,2,1),(1,0,1),(0,1,1),(1,1,0)]
def termcount(x):
 if isinstance(x,list):
  if len(x)==2 and isinstance(x[0],list)and len(x[0])==3 and all(isinstance(i,int)for i in x[0])and isinstance(x[1],str):return 1
  return sum(termcount(y)for y in x)
 if isinstance(x,dict):return sum(termcount(y)for y in x.values())
 return 0
def qgcd(a,b):return Q(math.gcd(a.numerator,b.numerator),math.lcm(a.denominator,b.denominator))
@lru_cache(None)
def weighted_lattice(box,linear):
 count=0
 for a in range(box[0][0],box[0][1]+1):
  for e in range(box[1][0],box[1][1]+1):
   lo,hi=box[2]
   for w,(lower,upper)in zip(WEIGHTS,linear):
    base=w[0]*a+w[1]*e
    if w[2]:lo=max(lo,(lower-base+w[2]-1)//w[2]);hi=min(hi,(upper-base)//w[2])
    elif not lower<=base<=upper:hi=lo-1;break
    if hi<lo:break
   count+=max(0,hi-lo+1)
 return count
@lru_cache(None)
def lattice(box,lo,hi):
 return sum(max(0,min(box[2][1],hi-a-e)-max(box[2][0],lo-a-e)+1)for a in range(box[0][0],box[0][1]+1)for e in range(box[1][0],box[1][1]+1))
def main(out):
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False);pins={}
 dump(out/'phase.json',dict(phase='source-only fixedK92target18 det11 profile; no scalar work'))
 retained_terms=0;load_peak=0;load_events=[];atom_max=0;atom_coefficients_checked=0
 def load(st,name='certificate.json'):
  nonlocal retained_terms,load_peak
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);x=read(p);full=termcount(x);load_peak=max(load_peak,retained_terms+full);assert load_peak+15000<=150000
  if st=='S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63':x={k:x[k]for k in['raw_compact14x5','known_row_units_rho14']}
  elif st=='S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63':x={'node_values':x['node_values']}
  elif st=='S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63':x={'Lframe_dual6x2':x['Lframe_dual6x2']}
  elif st=='S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63':x={'nonzero_selected_minor_Delta':x['nonzero_selected_minor_Delta'],'cleared_upper_row_coordinates7xr':[[],x['cleared_upper_row_coordinates7xr'][1]]}
  elif st=='S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63':x={'bottom4x6':x['bottom4x6']}
  kept=termcount(x);load_events.append(dict(path=str(p.relative_to(ROOT)),full_file_terms=full,retained_before=retained_terms,kept_terms=kept,parse_peak=retained_terms+full));retained_terms+=kept;return x
 st443=load('S61_NONSQUARE_PREPOWER_STAGED443_PREFLIGHT_20261104_v63','report.json');assert st443['all_passed']
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');h9=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');fx=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');fixed={r['old_DAG_root']:r['Laurent_polynomial']for r in fx['node_values']};bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];assert load('S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW','report.json')['all_passed']
 N=[];keys={};cache={};inventory=[]
 def put(op,args,bound,D,L,**kw):
  if op=='atom':c=kw.pop('content_lower');lp=kw.pop('content_l1');linear=kw.pop('linear_bounds')
  elif op=='mul':
   a,b=[N[j]for j in args];c=a['content_lower']*b['content_lower'];lp=a['content_l1']*b['content_l1'];linear=[[x[0]+y[0],x[1]+y[1]]for x,y in zip(a['linear_bounds'],b['linear_bounds'])]
  elif op in['add','sub']:
   a,b=[N[j]for j in args];c=qgcd(a['content_lower'],b['content_lower']);r=a['content_lower']/c;s=b['content_lower']/c;assert r.denominator==s.denominator==1;lp=a['content_l1']*r.numerator+b['content_l1']*s.numerator;linear=[[min(x[0],y[0]),max(x[1],y[1])]for x,y in zip(a['linear_bounds'],b['linear_bounds'])]
  elif op=='known_a_shift':
   a=N[args[0]];c=a['content_lower'];lp=a['content_l1'];linear=[[x+k*w[0]for x in bounds]for bounds,w in zip(a['linear_bounds'],WEIGHTS)for k in[kw['a_exponent_shift']]]
  else:raise AssertionError(op)
  if bound:bound=dict(bound,support=min(bound['support'],lp,weighted_lattice(tuple(map(tuple,bound['box'])),tuple(map(tuple,linear)))))
  i=len(N);N.append(dict(id=i,op=op,args=args,bound=bound,QQ_denominator=str(D),integer_l1_upper=str(L),content_lower=c,content_l1=lp,linear_bounds=linear,**kw));return i
 def atom(poly,ref,unit=Q(1),shift=0):
  nonlocal atom_max,atom_coefficients_checked
  p=[([m[0]+shift,m[1],m[2]],Q(c)/unit)for m,c in poly if Q(c)];key=('zero',)if not p else('literalconstant',str(p[0][1]))if len(p)==1 and p[0][0]==[0,0,0]else(json.dumps(ref,sort_keys=True),str(unit),shift);atom_max=max(atom_max,len(p));assert 3*len(p)<=15000
  if key in keys:return keys[key]
  D=math.lcm(*(c.denominator for m,c in p))if p else 1;L=sum(abs(c.numerator)*(D//c.denominator)for m,c in p)
  b=None if not p else dict(support=len(p),box=[[min(m[k]for m,c in p),max(m[k]for m,c in p)]for k in range(3)],degree=[min(sum(m)for m,c in p),max(sum(m)for m,c in p)])
  G=__import__('functools').reduce(math.gcd,(abs(c.numerator*(D//c.denominator))for m,c in p),0);content=Q(G,D)if G else Q(1);lp=L//G if G else 0
  for m,c0 in p:assert(c0/content).denominator==1
  atom_coefficients_checked+=len(p);linear=[[min(sum(w[k]*m[k]for k in range(3))for m,c in p),max(sum(w[k]*m[k]for k in range(3))for m,c in p)]for w in WEIGHTS]if p else[[0,0]for w in WEIGHTS]
  i=put('atom',[],b,D,L,reference=ref,QQ_unit=str(unit),a_shift=shift,content_lower=content,content_l1=lp,linear_bounds=linear);keys[key]=i;return i
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
 # Actualsourcefiles cannowbereleased; allfollowingobjects are support/content metadata.
 del cp,pd,fx,fixed,bot,h9
 def value(i):
  x=N[i]
  return None if x['bound']is None else dict(x['bound'],linear=x['linear_bounds'],content=x['content_lower'],Lprim=x['content_l1'])
 def bm(a,b):
  if a is None or b is None:return None
  box=[[a['box'][j][0]+b['box'][j][0],a['box'][j][1]+b['box'][j][1]]for j in range(3)];linear=[[x[0]+y[0],x[1]+y[1]]for x,y in zip(a['linear'],b['linear'])];lo=a['degree'][0]+b['degree'][0];hi=a['degree'][1]+b['degree'][1];lp=a['Lprim']*b['Lprim']
  return dict(support=min(a['support']*b['support'],lp,weighted_lattice(tuple(map(tuple,box)),tuple(map(tuple,linear)))),box=box,degree=[lo,hi],linear=linear,content=a['content']*b['content'],Lprim=lp)
 def ba(a,b):
  if a is None:return b
  if b is None:return a
  box=[[min(a['box'][j][0],b['box'][j][0]),max(a['box'][j][1],b['box'][j][1])]for j in range(3)];linear=[[min(x[0],y[0]),max(x[1],y[1])]for x,y in zip(a['linear'],b['linear'])];lo=min(a['degree'][0],b['degree'][0]);hi=max(a['degree'][1],b['degree'][1]);c=qgcd(a['content'],b['content']);r=a['content']/c;s=b['content']/c;assert r.denominator==s.denominator==1;lp=a['Lprim']*r.numerator+b['Lprim']*s.numerator
  return dict(support=min(a['support']+b['support'],lp,weighted_lattice(tuple(map(tuple,box)),tuple(map(tuple,linear)))),box=box,degree=[lo,hi],linear=linear,content=c,Lprim=lp)
 def bdet(m):
  state={0:value(one)}
  for r in m:
   nex={}
   for mask,a in state.items():
    for j,b in enumerate(r):
     if not(mask>>j&1)and b:nex[mask|1<<j]=ba(nex.get(mask|1<<j),bm(a,b))
   state=nex
  return state[(1<<len(m))-1]
 def va(poly):
  if not poly:return None
  D=math.lcm(*(Q(c).denominator for m,c in poly));ints=[Q(c).numerator*(D//Q(c).denominator)for m,c in poly];G=__import__('functools').reduce(math.gcd,map(abs,ints));assert G>0
  box=[[min(m[j]for m,c in poly),max(m[j]for m,c in poly)]for j in range(3)];linear=[[min(sum(w[j]*m[j]for j in range(3))for m,c in poly),max(sum(w[j]*m[j]for j in range(3))for m,c in poly)]for w in WEIGHTS];lp=sum(map(abs,ints))//G
  return dict(support=min(len(poly),lp,weighted_lattice(tuple(map(tuple,box)),tuple(map(tuple,linear)))),box=box,degree=[min(sum(m)for m,c in poly),max(sum(m)for m,c in poly)],linear=linear,content=Q(G,D),Lprim=lp)
 g=va(hc['raw_compact14x5'][6][1]);h=va(hc['raw_compact14x5'][10][1]);dif={tuple(m):Q(c)for m,c in hc['raw_compact14x5'][11][1]}
 for m,v in hc['raw_compact14x5'][9][1]:dif[tuple(m)]=dif.get(tuple(m),Q(0))-Q(v)
 hp=va([[list(m),str(v)]for m,v in dif.items()if v]);Delta=va(rk['nonzero_selected_minor_Delta']);l1,l2,l3=[va(p)for p in rk['cleared_upper_row_coordinates7xr'][1]]
 data=[((0,2,4),-1,h,l3,'-h*lambda428'),((0,2,5),1,h,l2,'h*lambda419'),((0,4,5),1,h,Delta,'h*Delta'),((1,2,4),-1,hp,l3,'-h120*lambda428'),((1,2,5),1,hp,l2,'h120*lambda419'),((1,4,5),1,hp,Delta,'h120*Delta'),((2,3,4),1,g,l3,'g*lambda428'),((2,3,5),-1,g,l2,'-g*lambda419'),((2,4,5),1,g,l1,'g*lambda505'),((3,4,5),1,g,Delta,'g*Delta')]
 terms=[];answer=None
 for k,(selected,hsign,a,b,label)in enumerate(data):
  remaining=[i for i in range(6)if i not in selected];m=[[value(j)for j in r]for r in rows[:5]+[rows[5+i]for i in remaining]];p8=bdet(m);w=bm(a,b);answer=ba(answer,bm(w,p8));terms.append(dict(H3_selected_lower6=list(selected),pure_remaining_lower6=remaining,H3_weight=label,H3_weight_sign=hsign,outer_sign=-1 if sum(selected)%2 else 1,pure8_bound=p8,weight_bound=w));dump(out/'phase.json',dict(phase='source-only weightedNewton/content metadata',completed_H3_terms=k+1))
 weighted=weighted_lattice(tuple(map(tuple,answer['box'])),tuple(map(tuple,answer['linear'])));loose=lattice(tuple(map(tuple,answer['box'])),*answer['degree']);assert weighted<=loose
 cuts=[]
 for i,w in enumerate(WEIGHTS):
  broad=[[sum(v[j]*answer['box'][j][0]for j in range(3)),sum(v[j]*answer['box'][j][1]for j in range(3))]for v in WEIGHTS];broad[0]=answer['linear'][0];broad[i]=answer['linear'][i];cuts.append(dict(weight=w,interval=answer['linear'][i],count_with_totaldegree=weighted_lattice(tuple(map(tuple,answer['box'])),tuple(map(tuple,broad)))))
 amin=answer['box'][0][0];L=answer['Lprim'];c0=answer['content'];matrixsupport=sum(N[i]['bound']['support']if N[i]['bound']else 0 for r in rows for i in r);uniqueatoms={i for i,x in enumerate(N)if x['op']=='atom'};atomterms=sum(N[i]['bound']['support']if N[i]['bound']else 0 for i in uniqueatoms)
 outbytes=weighted*(len(str(L))+45)+len(str(c0))+1000;naive_resident=weighted+matrixsupport+atomterms+15000
 # Check intervalcounter on artificialfiniteboxes, notsource/parameterpoints.
 for limit in[1,2,3]:
  box=((0,limit),)*3;linear=tuple((0,2*limit+(i%2))for i in range(len(WEIGHTS)));expected=sum(all(lo<=sum(w[j]*p[j]for j in range(3))<=hi for w,(lo,hi)in zip(WEIGHTS,linear))for p in itertools.product(range(limit+1),repeat=3));assert weighted_lattice(box,linear)==expected
 norm=dict(content_lower_divisor=str(c0),known_output_a_power=amin,normalization='actualnormalizedgdet11 = content_lower*a^amin*I, I integral; content_lowerneednotbemaximal',integer_L1_upper=str(L),integer_coefficient_bits_upper=L.bit_length(),primitive_coefficient_digits_upper=len(str(L)),no_exact_primitive_content_claim=True)
 dump(out/'normalization.json',norm);dump(out/'weighted-output.json',dict(weights=WEIGHTS,output_metadata=answer,all_weight_count=weighted,old_style_degree_box_count=loose,single_weight_comparisons=cuts,matrix_input_support_upper=matrixsupport,source_atom_terms=atomterms,ten_H3_surviving_terms=terms,row_known_a_factors=rowa,B_column_clearings=Bclear))
 sourcebudget=dict(file_load_events=load_events,max_input_parameter_term_peak=load_peak,transient_and_metadata_work_reserve=15000,full_source_only_owned_term_upper=load_peak+15000,largest_atom_terms=atom_max,all_input_atom_coefficients_content_checked=atom_coefficients_checked,never_materialized_output_lattice_points=True,weighted_counter='for(a,e),intersectalldintervals; d0weights arebooleanconstraints',all_integer_content_ratios_verified=True)
 dump(out/'source-resource.json',sourcebudget);dump(out/'source-circuit-metadata.json',dict(nodes=N,weights=WEIGHTS,matrix11x8_input_nodes=rows,no_actual_polynomial_products=True))
 proof=dict(QQcontent_atom='c=gcd(integercoefficients)/LCD; storeintegerL1afterdividingknownpositivec',QQcontent_product='knownfactorsmultiply(Gauss), integerL1boundsmultiply',QQcontent_sum='c=gcd_Q(c1,c2); both c1/c and c2/c integral, sointegerL1 <=L1*(c1/c)+L2*(c2/c)',known_content_only_not_maximal=True,
  all10_H3_terms_and_a_lambda419_retained=True,weighted_support_functions_propagated_through_every_source_operation=True,all_existing_source_unit_and_prefix_conventions_unchanged=True,no_parameter_samples_or_polynomial_gcd_factorization=True,
  scope='ThisisSOURCE/METADATAgate. Scalarwholeworkerrequiresseparateactualphaselifetime/input/outputprotocol; no determinantcoefficientscomputed.')
 dump(out/'proof.json',proof);dump(out/'phase.json',dict(phase='FINAL source-only Newton/content bound; no scalar'));note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<32*1024*1024 and load_peak+15000<=150000
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,source_metadata_gate_passes=True,source_nodes=len(N),weight_count=len(WEIGHTS),weighted_output_lattice=weighted,totaldegree_box_lattice=loose,known_content_lower=str(c0),integer_primitive_coefficient_bits_bound=L.bit_length(),primitive_integer_output_byte_bound=outbytes,source_input_peak_plus_reserve=load_peak+15000,matrix_input_support_upper=matrixsupport,scalar_all_matrix_input_output_reserve=naive_resident,
  scalar_whole_worker_gate_NOT_passed=True,scalar_dense_output_alone_exceeds150k=weighted>150000,scalar_integer_output_alone_exceeds32MiB=outbytes>=32*1024*1024,old238909_64506430_failure_preserved=True,no_scalar_or_new_coefficient_values=True,artifact_bytes_before_report=size,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2,default=str))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
