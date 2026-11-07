#!/usr/bin/env python3
"""Prepared fixed20 postcheck: covector-L90/denseOLD; no new case or rank search."""
import argparse,gc,hashlib,json,resource,signal,time,traceback
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_UPPER20_INDEPENDENT_POSTCHECK_20261105_v63'
DATA='reports/AMBIENT_B18_ELL2_T2_UPPER20_WHOLE_INPUT_GATE_V4_20261105_v63/whole20-coupled-pruned-cache-input-and-literal-owner-height.json'
PRIMARY='reports/AMBIENT_B18_ELL2_T2_UPPER20_VALUES_V4_20261105_v63/'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple,set)):return sum(map(count,v))
 if isinstance(v,F):return 1
 if isinstance(v,str):return int(all(x.lstrip('-').isdigit() for x in v.split('/')) and len(v.split('/'))<=2)
 return int(isinstance(v,(int,float)))
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return list(map(enc,v))
 return v
def main(out,gate,custody):
 def timeout(signum,frame):raise TimeoutError('POSTCHECK55s')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False)
 start=time.monotonic();pins={};raw=0;cursor={};pending={};work={};maximum=0;R={};Sc=[];aux={};points={};kernels={};C={};T={};first={};old3={};old2={};rawparts={};FIRSTvals={};results=[]
 def read(p):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size;assert raw<32*1024*1024;return json.loads(q.read_text())
 def dump(name,v):(out/name).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 def ck(q):
  nonlocal maximum
  maximum=max(maximum,abs(q.numerator).bit_length(),q.denominator.bit_length())
  if maximum>354:cursor['offending']=str(q);raise RuntimeError('POST354')
  return q
 def op(a,b,kind):
  nonlocal pending
  pending={'kind':kind,'a':str(a),'b':str(b),'cursor':dict(cursor)}
  assert max(abs(a.numerator).bit_length(),a.denominator.bit_length())+max(abs(b.numerator).bit_length(),b.denominator.bit_length())+1<=1063
  z=ck(a+b if kind=='add' else a*b);work['last_binary']=z;pending={};return z
 def add(a,b):return op(a,b,'add')
 def mul(a,b):return op(a,b,'mul')
 def parse(q):return ck(F(q))
 def env(q,v):
  work['env']=[q,v];den,u=map(int,v);assert den%q.denominator==0 and abs(q.numerator)*(den//q.denominator)<=u
 def dot(a,b):
  z=F(0)
  for i in reversed(range(len(a))):
   if a[i] and b[i]:z=add(z,mul(a[i],b[i]));work['dot']=z
  return z
 def put(v,k,q):
  z=add(v.get(k,F(0)),q)
  if z:v[k]=z
  else:v.pop(k,None)
 try:
  sg=read(gate);cu=read(custody);assert sg['all_passed'] and cu['all_passed'] and cu['worker_sha256']==sha(Path(__file__))
  for p,h in cu['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
  d=read(DATA);pr=read(PRIMARY+'report.json');sv=read(PRIMARY+'actual20-full-components-upper2-six-and-shared-tables.json');rank=read(PRIMARY+'actual20-fixed-integer-rank-all-factors-and-kernel.json')
  assert pins[DATA]=='dd8c9d761d6954cee49150166e11fa76c50071702f6450857c1523decf7e9342' and pins[PRIMARY+'report.json']=='253aff768ad00dfe9bdbfb5771bc2151fa2006ea472567697dbbdf8912a943a0'
  assert pr['all_passed'] and rank['rank']==0 and len(sv['all20'])==20 and sg['private_numeric_cap']==60000
  o=d['OLD'];n=d['NATIVE'];p=n['primitive'];Gamma=[{i:parse(q) for i,q in row} for row in o['Gamma_two']];D=[[list(map(parse,row)) for row in ch] for ch in o['D_two_by10_by56']]
  needed=sorted({i//4 for row in Gamma for i in row if i<48}|{i//4 for ch in D for row in ch for i,q in enumerate(row) if q and i<48});assert len(needed)==7
  N=[{i:[[parse(q) for q in cell] for cell in src[i]] for i in needed} for src in p['N_known_first8_by12_by9_by4']];G=[[parse(q) for q in row] for row in p['G8_C72']];nu=[[parse(q) for q in row] for row in p['nu_two_J18']]
  basis=[tuple(v) for v in p['basis18']];bi={v:i for i,v in enumerate(basis)};table=[[bi.get(tuple(a+b for a,b in zip(x,y))) for y in basis] for x in basis];known=p['N_known_positions'];assert all(table[i][j] is None for a in nu for b in nu for i,x in enumerate(a) if x for j,y in enumerate(b) if y);pairs=o['pair_order10'];vf={0:1,1:2,4:3,5:4,7:5,8:6}
  keys=set()
  for t in range(2):
   keys.update((0,t,-1,b) for b in [0,5,6,15,16,17])
   for k,(a,b) in enumerate(pairs):keys.update([(1,t,k,0),(1,t,k,7+k),(1,t,k,vf[a]),(1,t,k,vf[b])])
  assert len(keys)==90
  for key in sorted(keys,reverse=True):
   typ,t,k,target=key;cov=Gamma[t] if typ==0 else {i:q for i,q in enumerate(D[t][k]) if q};v=[F(0)]*32;R[key]=v
   for src in range(7,-1,-1):
    for im in range(3,-1,-1):
     z=F(0);cursor={'stage':'covector_L90','key':key,'source':src,'mask':im}
     for phys in reversed(needed):
      for kp,nb in reversed(list(enumerate(known))):
       for mask,x in reversed(list(enumerate(N[src][phys][kp]))):
        a=cov.get(4*phys+(mask^im),F(0));co=mask&im
        if not x or not a:continue
        if not co and nb==target:z=add(z,mul(a,x))
        elif co in [1,2]:
         for vb,y in enumerate(nu[co-1]):
          if y and table[nb][vb]==target:z=add(z,mul(mul(a,x),y))
     v[4*src+im]=z
  for rec in n['raw_program']:
   src=rec['native128'];im=rec['inputmask'];Sc=[[F(0)]*32 for _ in basis];cursor={'stage':'intrinsic_raw_S','raw':src}
   for row,g,sign in reversed(rec['Koszul_entries']):
    for flat in range(71,-1,-1):
     x=G[g][flat]
     if not x:continue
     b,mask=divmod(flat,4);co=mask&im;col=4*row+(mask^im);q=mul(F(sign),x)
     if not co:Sc[b][col]=add(Sc[b][col],q)
     elif co in [1,2]:
      for vb,y in enumerate(nu[co-1]):
       dst=table[b][vb]
       if y and dst is not None:Sc[dst][col]=add(Sc[dst][col],mul(q,y))
   points[src]={i:q for i,q in enumerate(Sc[0]) if q};assert points[src]=={i:F(q) for i,q in dict(n['raw_point_S0_columns48'])[src]};aux[src]=[]
   for t in range(2):
    a=dot(R[0,t,-1,17],Sc[0]);b=F(0)
    for rb,sb in [(0,17),(5,16),(6,15),(15,6),(16,5)]:a=add(a,dot(R[0,t,-1,rb],Sc[sb]))
    for k,(x,y) in enumerate(pairs):
     b=add(b,dot(R[1,t,k,7+k],Sc[0]));b=add(b,dot(R[1,t,k,0],Sc[7+k]))
     for xx,yy in ([(x,y)] if x==y else[(x,y),(y,x)]):b=add(b,dot(R[1,t,k,vf[xx]],Sc[vf[yy]]))
    aux[src].append([a,b]);assert [a,b]==list(map(F,sv['native_raw_auxiliaries_FORMAL_partialN_extensions'][str(src)][t]))
  fines={r['kernel105']:r['two_channels_Gamma_D_sum_DU'] for r in n['kernel_fine']}
  for rec in n['restrictions']:
   k=rec['kernel105'];pt={};kernels[k]=[]
   for j,q in rec['original_native_source_combination']:
    for row,x in points[j].items():put(pt,row,mul(F(q),x))
   assert not pt
   for t in range(2):
    v=[F(0),F(0)]
    for j,q in reversed(rec['original_native_source_combination']):
     for c in range(2):v[c]=add(v[c],mul(F(q),aux[j][t][c]))
    z=add(*v)
    for x,eb in zip(v+[z],fines[k][t]):env(x,eb)
    kernels[k].append(z)
  for k in n['source_zero_kernel_ids']:assert all(int(v[1])==0 for row in fines[k] for v in row);kernels[k]=[F(0),F(0)]
  assert {str(k):list(map(str,v)) for k,v in kernels.items()}==sv['native_actual_kernel_restrictions']
  dump('phase1-native-L90-raw16-and-restricted25.json',{'FORMAL_raw16':aux,'actual_kernel25_including10_SOURCEzero':kernels,'points':points,'samechosenpartialN_extension_verified_notactualrawjets':True})
  N=None;G=None;nu=None;R.clear();Sc=[];work.clear()
  # Dense ordinary B21, different multiplication/summation order.
  b=[tuple(v) for v in o['ordinary_B21_basis']];idx={v:i for i,v in enumerate(b)};tab=[[idx.get(tuple(x+y for x,y in zip(a,c))) for c in b] for a in b];leaves={i:list(map(parse,v)) for i,mode,v in o['boundary_cache']};C=dict(leaves)
  def pm(a,b):
   z=[F(0)]*21;work['poly']=z
   for i in range(20,-1,-1):
    if a[i]:
     for j in range(20,-1,-1):
      k=tab[i][j]
      if b[j] and k is not None:z[k]=add(z[k],mul(a[i],b[j]))
   return z
  for i,r in o['literal_C0_program']:
   cursor={'stage':'C0_FORMAL','node':i}
   if r[0]=='c':v=[F(0)]*21;v[0]=parse(F(r[1],r[2]))
   elif r[0]=='n':v=[-q for q in C[r[1]]]
   elif r[0]=='a':v=[add(x,y) for x,y in zip(C[r[1]],C[r[2]])]
   else:assert r[0]=='m';v=pm(C[r[1]],C[r[2]])
   C[i]=v
  gc=dict(o['OLD_Gamma_columns']);dc=dict(o['OLD_D_columns']);old3={j:[F(0)]*2 for j in gc};old2={j:[F(0)]*2 for j in dc}
  for j,roots in gc.items():
   for row,(kind,i) in roots:
    if kind=='C0':
     for t in range(2):old3[j][t]=add(old3[j][t],mul(Gamma[t].get(row,F(0)),C[i][20]))
  C.clear()
  for i,r in o['typed_program']:
   cursor={'stage':'typed_FORMAL','node':i}
   if r[0]=='actual_C0_scalar_function':v=leaves[r[1]][:]
   elif r[0]=='sum':
    v=[F(0)]*21
    for j,a,c in reversed(r[1]):v=[add(x,mul(F(a,c),y)) for x,y in zip(v,T[j])]
   else:
    assert r[0]=='product';v=[F(0)]*21;v[0]=F(r[1],r[2])
    for j in reversed(r[3]):v=pm(v,T[j])
   T[i]=v
  for j,roots in gc.items():
   for row,(kind,i) in roots:
    if kind=='typed':
     for t in range(2):old3[j][t]=add(old3[j][t],mul(Gamma[t].get(row,F(0)),T[i][20]))
  for j,roots in dc.items():
   for row,i in roots:
    for t in range(2):
     for k in range(10):old2[j][t]=add(old2[j][t],mul(D[t][k][row],T[i][10+k]))
  dump('phase2-shared-OLD21-column-coefficients-FORMAL.json',{'OLD_Gamma':old3,'OLD_D':old2,'all_raw_helpers_FORMAL_untilcorrected':True});T.clear();C.clear();leaves.clear();work.clear()
  L0={j:dict((i,F(q)) for i,q in v) for j,v in o['native_L0_selected_columns']};Le={j:dict((i,F(q)) for i,q in v) for j,v in o['native_Le_selected_columns']}
  for r in o['FIRST_column_recipes']:
   v=[F(0)]*56;rr=r['first_column_recipe']
   if rr['kind']=='literal56':
    for i,q in rr['column']:v[i]=parse(q)
   else:
    for label,mat in [('S0',Le),('Se',L0)]:
     for j,q in reversed(rr[label]):
      for i,x in mat[j].items():v[i]=add(v[i],mul(F(q),x))
   first[r['pivot27'],r['variable9']]=v
  c0=dict(o['original_source_allc0_raw21']);a1=dict(o['original_full_a1_raw21']);piv=o['point27_source_indices'];v92=dict(o['accepted_native_Gamma_plane_targets'])
  for case in o['required_raw_cases21']:
   rawparts[case]=[];FIRSTvals[case]=[]
   for t in range(2):
    cursor={'stage':'raw21_FORMAL','case':case,'channel':t};a=F(0);bb=F(0);carry=F(0);outer=F(0);fv=[]
    for j,q in reversed(c0[case]):
     if j<640:a=add(a,mul(F(q),old3[j][t]));bb=add(bb,mul(F(q),old2[j][t]))
    for slot,q in a1[case]:
     pv,var=divmod(slot,9)
     for k,(x,y) in enumerate(pairs):
      for oth in ([y] if var==x else[])+([x] if var==y and x!=y else[]):z=mul(F(q),dot(D[t][k],first[pv,oth]));fv.append(z);carry=add(carry,z)
     if pv>=4 and var in [7,8]:outer=add(outer,mul(F(q),F(v92[piv[pv]][t][0 if var==8 else 1])))
    rawparts[case].append([a,bb,carry,outer]);FIRSTvals[case].append(fv)
    assert [a,bb,carry,outer]==list(map(F,sv['raw21_partial_components_FORMAL_until605correction'][str(case)][t])) and fv==list(map(F,sv['all_raw21_FIRST_factor_values'][str(case)][t]))
  first.clear();work.clear();coords=dict(n['case_coordinates']);savedcases={r['case619']:r for r in sv['all20']};Q=[dict(row) for row in o['Q_two6']]
  for case,weights in o['corrected20_source_recipes']:
   up=[];pieces=[]
   for t in range(2):
    v=list(rawparts[case][t])
    for base,q in zip([563,609],weights):
     if F(q):v=[add(a,-mul(F(q),bb)) for a,bb in zip(v,rawparts[base][t])]
    native=F(0)
    for k,q in reversed(coords[case]):native=add(native,mul(F(q),kernels[k][t]))
    z=add(add(v[0],v[1]),native);z=add(z,-add(v[2],v[3]));up.append(z);piece=[v[0],v[1],native,v[2],v[3]];pieces.append(piece)
    sk=savedcases[case]['two_source_component_pieces'][t];assert piece==[F(sk[k]) for k in ['OLD_Gamma','OLD_D','native_Gamma_plus_D','FIRST_to_subtract','native92_outer_to_subtract']]
   full=[add(mul(F(Q[1].get(i,'0')),up[1]),mul(F(Q[0].get(i,'0')),up[0])) for i in range(6)];assert up==[F(0)]*2==list(map(F,savedcases[case]['actual_upper_two'])) and full==[F(0)]*6==list(map(F,savedcases[case]['actual_full_six']));results.append({'case619':case,'actual_upper':up,'actual_six':full,'pieces':pieces})
  assert rank['rank']==0 and rank['selected_original20_columns']==[] and rank['fixed_minor']==0 and rank['actual_integer_matrix2_by20']==[[0,0]]*20 and rank['all20_rational_factors']==[[]]*20
  assert rank['cleared_kernel_vectors']==[[[j,1]] for j in range(20)]
  live=2*count(d)+2*count(sv)+count([Gamma,D,aux,points,kernels,old3,old2,rawparts,FIRSTvals,results,L0,Le,rank,pr,pins])+6000;assert live<60000
  dump('fixed20-independent-nativeOLD-FIRST-fullcomponents-and-zero-rank.json',{'all20':results,'all1742_FIRST_matches':True,'all_raw21_FORMAL_partial_matches':True,'all16_FORMALaux_and25_actual_orSOURCEzero_kernel_matches':True,'complete605_upperzero_after584SOURCE_and554prior':True,'remaining12opaque24seconds_OPEN':True,'no_HF2_length_H0_claim':True,'fixedrank0_factor_kernel_certificate_verified':True})
  assert raw+sum(x.stat().st_size for x in out.iterdir())<sg['RAW_cap']
  ep='fixed20-independent-nativeOLD-FIRST-fullcomponents-and-zero-rank.json';report={'id':STEM,'all_passed':True,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'actual_QQ_bits':maximum,'live_upper':live,'private_cap60000':True,'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start};dump('report.json',report);print(json.dumps(report));signal.alarm(0)
 except BaseException as error:
  signal.alarm(0);message=repr(error);tb=error.__traceback__;traceback.clear_frames(tb);error.__traceback__=None;error.__context__=None;error.__cause__=None;del tb;gc.collect();signal.signal(signal.SIGALRM,signal.SIG_DFL);signal.alarm(5)
  dump('FAILURE_PREFIX.json',{'error':message,'cursor':cursor,'pending':pending,'work':work,'R':R,'Sc':Sc,'aux':aux,'points':points,'kernels':kernels,'C':C,'T':T,'first':first,'old3':old3,'old2':old2,'rawparts':rawparts,'FIRSTvals':FIRSTvals,'results':results,'pins':pins});raise
if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--output',required=True,type=Path);pa.add_argument('--gate',required=True);pa.add_argument('--custody',required=True);a=pa.parse_args();main(a.output,a.gate,a.custody)
