#!/usr/bin/env python3
"""SOURCE ONLY det11withgfactored; boundedadaptiveintegerarithmetic/prefixserialization gate."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_PREPOWER_K92_T18_G_FACTORED_GUARDED_PREFLIGHT_20261104_v63'
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
  elif st=='S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63':x={'nonzero_selected_minor_Delta':x['nonzero_selected_minor_Delta'],'cleared_upper_row_coordinates7xr':[[],x['cleared_upper_row_coordinates7xr'][1],[],x['cleared_upper_row_coordinates7xr'][3],x['cleared_upper_row_coordinates7xr'][4]]}
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
 # Compute the four unmaterialized H3 minors by support/content only.
 H=hc['raw_compact14x5'];hcrows={i:[va(H[i][j])for j in[1,3,2]]for i in[6,8,9,10,11,12]}
 HA=hcrows[10];HB=[ba(x,y)for x,y in zip(hcrows[11],hcrows[9])];HD=hcrows[8];HC=hcrows[6];HE=hcrows[12];HF=hcrows[9]
 cofactor_proofs=[]
 def quotient_bound(mat,p,q,label):
  direct=bdet(mat);exactc=p['content']*q['content']/g['content'];ratio=exactc/direct['content'];assert ratio.denominator==1 and ratio>0
  box=[[max(direct['box'][j][0],p['box'][j][0]+q['box'][j][0]-g['box'][j][0]),min(direct['box'][j][1],p['box'][j][1]+q['box'][j][1]-g['box'][j][1])]for j in range(3)]
  linear=[[max(a[0],b[0]+c[0]-d[0]),min(a[1],b[1]+c[1]-d[1])]for a,b,c,d in zip(direct['linear'],p['linear'],q['linear'],g['linear'])];assert all(lo<=hi for lo,hi in box+linear)
  lp=direct['Lprim']//ratio.numerator;assert lp>0;result=dict(support=min(direct['support'],lp,weighted_lattice(tuple(map(tuple,box)),tuple(map(tuple,linear)))),box=box,degree=linear[0],linear=linear,content=exactc,Lprim=lp)
  cofactor_proofs.append(dict(label=label,known_identity='g*H3minor=sign*atomproduct',exact_content_from_Gauss=str(exactc),direct_lower_content=str(direct['content']),integer_ratio=str(ratio),direct_integerL1_bits=direct['Lprim'].bit_length(),new_integerL1_bits=lp.bit_length(),quotient_support_only_no_coefficientdivision=True))
  return result
 ade=quotient_bound([HA,HD,HE],h,l3,'ADE');adf=quotient_bound([HA,HD,HF],h,l2,'ADF');bde=quotient_bound([HB,HD,HE],hp,l3,'BDE');bdf=quotient_bound([HB,HD,HF],hp,l2,'BDF')
 la=va(rk['cleared_upper_row_coordinates7xr'][3][0]);lb=va(rk['cleared_upper_row_coordinates7xr'][4][0])
 data=[((0,2,4),-1,ade,'minor(A,D,E)'),((0,2,5),1,adf,'minor(A,D,F)'),((0,4,5),1,la,'lambda415'),((1,2,4),-1,bde,'minor(B,D,E)'),((1,2,5),1,bdf,'minor(B,D,F)'),((1,4,5),1,lb,'lambda442'),((2,3,4),1,l3,'lambda428'),((2,3,5),-1,l2,'-lambda419'),((2,4,5),1,l1,'lambda505'),((3,4,5),1,Delta,'Delta')]
 terms=[];answer=None
 for k,(selected,hsign,w,label)in enumerate(data):
  remaining=[i for i in range(6)if i not in selected];m=[[value(j)for j in r]for r in rows[:5]+[rows[5+i]for i in remaining]];p8=bdet(m);answer=ba(answer,bm(w,p8));terms.append(dict(H3_selected_lower6=list(selected),pure_remaining_lower6=remaining,H3_weight=label,H3_weight_sign=hsign,outer_sign=-1 if sum(selected)%2 else 1,pure8_bound=p8,weight_bound=w));dump(out/'phase.json',dict(phase='source-only weightedNewton/content metadata',completed_H3_terms=k+1))
 weighted=weighted_lattice(tuple(map(tuple,answer['box'])),tuple(map(tuple,answer['linear'])));loose=lattice(tuple(map(tuple,answer['box'])),*answer['degree']);assert weighted<=loose
 cuts=[]
 for i,w in enumerate(WEIGHTS):
  broad=[[sum(v[j]*answer['box'][j][0]for j in range(3)),sum(v[j]*answer['box'][j][1]for j in range(3))]for v in WEIGHTS];broad[0]=answer['linear'][0];broad[i]=answer['linear'][i];cuts.append(dict(weight=w,interval=answer['linear'][i],count_with_totaldegree=weighted_lattice(tuple(map(tuple,answer['box'])),tuple(map(tuple,broad)))))
 amin=answer['box'][0][0];L=answer['Lprim'];c0=answer['content'];matrixsupport=sum(N[i]['bound']['support']if N[i]['bound']else 0 for r in rows for i in r);uniqueatoms={i for i,x in enumerate(N)if x['op']=='atom'};atomterms=sum(N[i]['bound']['support']if N[i]['bound']else 0 for i in uniqueatoms)
 outbytes=weighted*(len(str(L))+45)+len(str(c0))+1000;naive_resident=weighted+matrixsupport+atomterms+15000
 # Check intervalcounter on artificialfiniteboxes, notsource/parameterpoints.
 for limit in[1,2,3]:
  box=((0,limit),)*3;linear=tuple((0,2*limit+(i%2))for i in range(len(WEIGHTS)));expected=sum(all(lo<=sum(w[j]*p[j]for j in range(3))<=hi for w,(lo,hi)in zip(WEIGHTS,linear))for p in itertools.product(range(limit+1),repeat=3));assert weighted_lattice(box,linear)==expected
 # The scalarattempt below is only a finite guarded operatorplan, not a run.
 # These caps are BEFORE-COMMIT limits and never assert deterministiccompletion.
 CELL_CAP=135000;WORK=15000;POLY_CAP=1024;COEFF_BITS=1024;CONTENT_BITS=4096;PENDING_BITS=16384
 assert CELL_CAP+WORK==150000
 atomInput=0;atomMaxBits=0
 for x in N:
  if x['op']=='atom'and x['bound']:
   atomInput+=x['bound']['support']*(len(str(x['content_l1']))+45)+len(str(x['content_lower']))+128;atomMaxBits=max(atomMaxBits,x['content_l1'].bit_length())
 extraPolys=[p for row in hc['raw_compact14x5']for p in row]+[rk['nonzero_selected_minor_Delta']]+rk['cleared_upper_row_coordinates7xr'][1]+rk['cleared_upper_row_coordinates7xr'][3]+rk['cleared_upper_row_coordinates7xr'][4]
 for poly in extraPolys:
  v=va(poly)
  if v:atomInput+=v['support']*(len(str(v['Lprim']))+45)+len(str(v['content']))+128;atomMaxBits=max(atomMaxBits,v['Lprim'].bit_length())
 # Operator/table text allowance; actual preparation must stop before this cap.
 nativeInputUpper=atomInput+500000;INPUT_CAP=3*1024*1024;METADATA_CAP=2*1024*1024;LOG_CAP=512*1024;MANIFEST_CAP=512*1024
 assert nativeInputUpper<INPUT_CAP and atomMaxBits<=COEFF_BITS
 # 3 signed32 exponents +uint16 bytecount+sign+128 magnitudebytes <144.
 termBytes=144;polyBytes=2*(CONTENT_BITS//8)+192
 prefixUpper=CELL_CAP*termBytes+POLY_CAP*polyBytes+4096+INPUT_CAP+METADATA_CAP+LOG_CAP+MANIFEST_CAP
 assert prefixUpper<32*1024*1024
 norm=dict(output_is_det11_with_g_kept_factored=True,content_lower_divisor=str(c0),known_output_a_power=amin,normalization='normalizeddet11=positiveQQcontent_lower*a^amin*IwithIintegral;contentneednotbemaximal',integer_L1_upper=str(L),integer_coefficient_bits_upper=L.bit_length(),g_only_factored_in_actualentry='c8*g*det11/(Delta*a²*N7)',NO_polynomialdivisionby_g=True,cofactor_content_Newton_proofs=cofactor_proofs)
 dump(out/'normalization.json',norm);dump(out/'weighted-output.json',dict(weights=WEIGHTS,output_metadata=answer,weighted_count=weighted,old_style_box_totaldegree_count=loose,single_weight_counts=cuts,all10_H3_terms=terms,row_known_a_factors=rowa,B_column_clearings=Bclear,matrix_input_support_upper=matrixsupport))
 procedure=dict(actual_source='current15prefix,notbare135;all7remainingtailsandtail2effectsretained',source_matrix11x8_nodes=rows,source_nodes_external_or_included='source-circuit-metadata.json',
  source_H3_rows=['N10','N11-N9','N8','N6','N12','N9'],ten_H3_nonzero_coefficients='sixfrozencoefs+fourliteral3x3Hminors;all10signsfromweighted-output.json',
  H3_cofactor_exactcontent='FournewcofscomputedonlyfromsmallH3polynomials,thenrescaledbyprovenexactQQcontentusingwholecoefficientintegerdivisionchecks.No polynomialgdivision.',
  determinant_expression='Sumall10 H3coef * signedpure8det. Eachpure8detusesrowstreamedsubsetDPwithTWOadjacentlayers; NOunapprovedpolynomialgcd/factorization.',
  source_evaluation='Sourceatomstreamreadintoownedcells;sourceDAGlast-use. KnowncontentfactorsperoperationassignafixedQQscalar;integercontentratiosarecheckedbeforeanyaddmul.',
  sparse_addmul='Eachproductcoefficientisgeneratedinlexicographicinputkeyorder. Computeonependingintegerupdate,checkBEFOREcommitnewcell/bitgrowth. Zeroresultsreleasecells. Nofullproductbuffer.',
  source_and_DP_content='Allsourceadditions/productsandDPcandidatecontentrulesareprecomputablegcdQ/multiplication.Metadatafactorensuresintegerquotientscalesforeverypartialsum.',
  release='AfteraDPtransitionhasconsumedapreviousstateitscellsarereleased. Aftereachpure8detisaccumulatedreleaseallitsDPstates. StreamH3cofactors ifneeded;neverdropamathterm.',
  stop='Beforecommitstopifownedcellpoolwouldexceed135000,activepolys1024,coefficient1024bits,QQcontent4096bits,metadata/input/logcaps,or45sclock. Saveexactcurrentstate/pendingoperationcursor. Noautomaticretry.',
  completion='OnlycompletionofALL10terms/ALLsourceandDPloopspermitsfullscalar/nonzeroclaim. A nonzeroprefixisNOTnonzerodet11.',
  no_scalar_job_executed=True)
 dump(out/'guarded-operator-recipe.json',procedure)
 binary=dict(version=1,ring_order=['a','e','d'],term_record='3*signedint32exponent,uint16magnitudebytecount,uint8sign,magnitudebytes<=128;pad144byteupper',poly_record='id/commonQQcontentnumden<=4096bits/rowshift/sourceopandDPcursors<=192headerbytes',
  term_cap=CELL_CAP,term_record_upper=termBytes,active_poly_cap=POLY_CAP,poly_record_upper=polyBytes,pending_update_cap_bits=PENDING_BITS,pending_header_and_integer_bytes_upper=4096,
  input_cap_bytes=INPUT_CAP,operator_metadata_cap_bytes=METADATA_CAP,log_cap_bytes=LOG_CAP,report_manifest_cap_bytes=MANIFEST_CAP,total_safe_prefix_raw_bytes=prefixUpper,
  source_input_preparation_byte_upper=nativeInputUpper,actual_source_atom_primitive_bit_upper=atomMaxBits,
  checkpoint_once='WriteONCEatguard/completion,neverretainoldlargecheckpointandnewtemporarycopy. Streambinaryrecords;noJSON/base64duplicationofcoefficients.',
  checkpoint_contents='Allactivepolys/commonQQscalars,sparsekeys/coeffs,layerstates/globalaccumulator,operationnumber,row/mask/left/righttermindices,allinput/sourcehashes,pendinguncommittedupdate.',
  pending_scale_bound='Atmostthree4096bitcontentfactorsplus2*1024bitcoefficients<16384pendingbits;checkbeforemultiplyingiflarger.',
  actual_RAW_guard='Beforeeachfilewritecheckcurrenttotal+nextrecord<=32MiB;savecomplete-or-partial-file statushonestly. Noallocation/outputcapincrease.')
 dump(out/'safe-prefix-format.json',binary)
 sourcebudget=dict(file_load_events=load_events,input_peak=load_peak,source_input_plus_reserve=load_peak+15000,max_atom_terms=atom_max,all_input_atom_coefficients_content_checked=atom_coefficients_checked,output_lattice_points_NEVER_materialized=True,source_only_metadata_gate_passes=load_peak+15000<=150000)
 dump(out/'source-resource.json',sourcebudget);dump(out/'source-circuit-metadata.json',dict(nodes=N,weights=WEIGHTS,matrix11x8_input_nodes=rows,no_new_polynomial_coefficients=True))
 dump(out/'phase.json',dict(phase='FINAL source-only g-factored/guarded-prefix gate; no scalar'))
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));size=sum(p.stat().st_size for p in out.iterdir());assert size<METADATA_CAP and load_peak+15000<=150000
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,g_kept_factored_and_knownunit_on_N6_chart=True,source_nodes=len(N),weighted_det11_output_lattice=weighted,known_content_lower=str(c0),primitive_integer_coefficient_bits=L.bit_length(),primitive_integer_full_output_byte_bound=outbytes,ordinary_full_output_gate_passes=weighted<=150000 and outbytes<32*1024*1024,
  guarded_attempt_gate_passes=True,completion_NOT_guaranteed=True,allocated_parameter_coefficient_cell_cap=CELL_CAP,work_reserve=WORK,total_owned_hardcap=CELL_CAP+WORK,coefficient_bit_stop_cap=COEFF_BITS,QQcontent_bit_stop_cap=CONTENT_BITS,safe_prefix_total_raw_bytes=prefixUpper,native_input_byte_upper=nativeInputUpper,source_input_peak_plus_reserve=load_peak+15000,source_artifacts_before_report=size,
  all10_H3_terms_and_factor_a_lambda419_retained=True,NO_actual_scalar_or_new_coefficient_values=True,older_FAIL_reports_unchanged=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-started)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2,default=str))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
