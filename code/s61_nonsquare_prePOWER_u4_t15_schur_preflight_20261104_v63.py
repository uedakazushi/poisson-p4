#!/usr/bin/env python3
"""ONE selected Schur entry: typed source recipe and rigorous support/l1/lifetime bounds."""
import argparse,hashlib,importlib.util,itertools,json,math,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_U4_T15_SCHUR_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();assert not(out/'report.json').exists();paths=[]
 def load(st,file='certificate.json'):
  p=ROOT/'reports'/st/file;paths.append(p);return read(p)
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63');rank=load('S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63');h=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json');fixed=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');strip=load('S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63','native.json')
 direct=load('S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63','operator.json');assert sha(paths[-1])=='4bbcd0b3aa6dcb70c991ed75ce0f213cd8253a3e060736ce63b363d1a3123397'
 dr=load('S61_NONSQUARE_PREPOWER_DIRECT_P4_H5_SCHUR_20261104_v63','report.json');assert sha(paths[-1])=='264049acfbbc7abeae10dfd5994b03967b599d5a0ea18f0a10c07efb4eb57b11'
 assert direct['context']['P4_current15_rows']==[0,1,3,6]and direct['context']['remaining6_current15'][0]==5 and direct['context']['remaining13_original32_seed_labels'][9]==15
 hp=ROOT/'code/s61_nonsquare_prePOWER_h22_kernel9_support_preflight_20261104_v63.py';paths.append(hp);sp=importlib.util.spec_from_file_location('bounds',hp);b=importlib.util.module_from_spec(sp);sp.loader.exec_module(b)
 DAG=[];cache={};literal={};cached={r['old_DAG_root']:r['Laurent_polynomial']for r in fixed['node_values']}
 def put(op,bb,dd,nn,args=[],**meta):
  n=len(DAG);DAG.append(dict(id=n,op=op,args=args,bound=b.enc(bb),common_QQ_denominator=str(dd),integer_l1_upper=str(nn),**meta));return n
 def B(i):
  x=DAG[i]['bound'];return None if x['zero']else(x['support_upper'],tuple(tuple(t)for t in x['Laurent_exponent_box']),x['known_B_denominator_power'])
 def DN(i):return int(DAG[i]['common_QQ_denominator']),int(DAG[i]['integer_l1_upper'])
 def atom(p,reference,den_unit=Q(1),a_shift=0,Bpower=0):
  cc=[Q(c)/den_unit for m,c in p];D=math.lcm(*(c.denominator for c in cc))if cc else 1;N=sum(abs(c.numerator)*(D//c.denominator)for c in cc);bb=b.shift(b.atom(p),a_shift)
  if bb is not None:bb=(bb[0],bb[1],Bpower)
  key=json.dumps([p,str(den_unit),a_shift,Bpower],sort_keys=True,separators=(',',':'));literal.setdefault(key,len(p));return put('atom',bb,D,N,reference=reference,QQ_denominator_unit=str(den_unit),a_shift=a_shift,B_denominator_power=Bpower)
 def const(c):return atom([[[0,0,0],str(Q(c))]],dict(constant=str(c)))
 def add(x,y,sign=1):
  bx,by=B(x),B(y);bp=max(0 if bx is None else bx[2],0 if by is None else by[2]);d,n=DN(x);e,m=DN(y);D=math.lcm(d,e);N=n*(D//d)*179**(bp-(0 if bx is None else bx[2]))+m*(D//e)*179**(bp-(0 if by is None else by[2]));return put('add'if sign==1 else'sub',b.aa(bx,by),D,N,[x,y])
 def mul(x,y):
  d,n=DN(x);e,m=DN(y);return put('mul',b.mm(B(x),B(y)),d*e,n*m,[x,y])
 def summ(xs):
  z=const(0)
  for x in xs:z=add(z,x)
  return z
 def imp(old):
  if old in cache:return cache[old]
  n=q['scalar_DAG']['nodes'][old];op=n[0]
  if old in cached:x=atom(cached[old],dict(path='reports/S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63/native.json',key='node_values',old_DAG_root=old))
  elif op=='z':x=const(0)
  elif op=='one':x=const(1)
  elif op in['a','s','m']:
   a,c=imp(n[1]),imp(n[2]);x=mul(a,c)if op=='m'else add(a,c,-1 if op=='s'else 1)
  else:
   num,den=cp['matrix25x32'][n[1]][n[2]]if op=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if op=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
   if len(den)==1:assert den[0][0][1:]==[0,0];unit=Q(den[0][1]);ashift=-den[0][0][0];bp=0
   else:
    dd={tuple(m):Q(c)for m,c in den};bb={tuple(m):Q(c)for m,c in b.BP};assert dd.keys()==bb.keys();unit=dd[next(iter(bb))]/bb[next(iter(bb))];assert all(dd[m]==unit*bb[m]for m in bb);ashift=0;bp=1
   x=atom(num,dict(path='accepted_Q2L_source_atom',old_DAG_root=old,definition=n),unit,ashift,bp)
  cache[old]=x;return x
 cols=q['new22_target_original32_columns'][:17];piv=[22,11,13,25];pure_rows=[0,1,3,6]
 def P(i,j):
  if i==0 and j in h['original_target_labels']:return atom(h['source_h9'][h['original_target_labels'].index(j)],dict(path='reports/S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63/prePOWER-h9-source.json',key='source_h9',target=j))
  return imp(q['remaining15x22_scalar_roots'][i][cols.index(j)])
 A=[[P(i,j)for j in piv]for i in pure_rows];D4=atom(strip['exact_unscaled_pure4_minor'],dict(path='reports/S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63/native.json',key='exact_unscaled_pure4_minor'))
 # Pair-only u=4: no lower3-H cofactor appears in this selected row.
 N=hc['raw_compact14x5'];assert all(not N[4][j]for j in[1,2,3]);na=lambda i,j:atom(N[i][j],dict(path='reports/S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63/certificate.json',key='raw_compact14x5',index=[i,j]))
 det2=atom(rank['literal_pair_determinant'],dict(path='reports/S61_NONSQUARE_RAW_H_GENERIC_RANK_20261104_v63/certificate.json',key='literal_pair_determinant'))
 mu0=add(mul(na(4,0),na(3,4)),mul(na(4,4),na(3,0)),-1);mu1=add(mul(na(4,4),na(1,0)),mul(na(4,0),na(1,4)),-1)
 indices=[4,1,3];pw=[hc['known_row_units_rho14'][i]['a_d_F_B_C_J_powers']for i in indices];lcm=[max(r[j]for r in pw)for j in range(6)]
 names=[[[[1,0,0],'1']],[[[0,0,1],'1']],[[[1,0,0],'500'],[[0,2,0],'90'],[[0,1,0],'3'],[[0,0,0],'-24']],b.BP,[[[1,0,0],'100'],[[0,2,0],'18'],[[0,1,0],'-9'],[[0,0,0],'-2']],[[[1,0,1],'5'],[[0,1,0],'-6'],[[0,0,0],'-4']]]
 src=[]
 for i,co,sgn in zip(indices,[det2,mu0,mu1],[1,-1,-1]):
  row=hc['known_row_units_rho14'][i];x=const(Q(sgn)/Q(row['constant']))
  for j,p in enumerate(lcm):
   for _ in range(p-row['a_d_F_B_C_J_powers'][j]):x=mul(x,atom(names[j],dict(namedFULLUunit=['a','d','F','B','C','J'][j])))
  src.append(mul(co,x))
 # Form only5 needed Bhat components, not a full17-column row.
 Bu=[summ([mul(c,P(i+1,j))for i,c in zip(indices,src)])for j in piv+[15]]
 aremain=[P(i,15)for i in pure_rows]
 # Adj4 cofactor computed-and-consumed per position, never all16 retained.
 Z=[]
 for k in range(4):
  terms=[]
  for l in range(4):
   mat=[r[:k]+r[k+1:]for r in A[:l]+A[l+1:]];minor=const(0)
   for perm in itertools.permutations(range(3)):
    z=const(1)
    for i,j in enumerate(perm):z=mul(z,mat[i][j])
    parity=sum(perm[i]>perm[j]for i in range(3)for j in range(i+1,3))+(k+l);minor=add(minor,z,-1 if parity%2 else 1)
   terms.append(mul(minor,aremain[l]))
  Z.append(summ(terms))
 S=mul(D4,Bu[4])
 for k in range(4):S=add(S,mul(Bu[k],Z[k]),-1)
 # Exact last-use upper bound over typed scalar-operation DAG.
 need=set();todo=[S]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i);todo+=DAG[i]['args']
 use={i:0 for i in need}
 for i in need:
  for j in DAG[i]['args']:use[j]+=1
 use[S]+=1;live=peak=0;peakid=0
 for i in sorted(need):
  n=DAG[i]['bound']['support_upper'];live+=n
  # A multiply/add may require an additional output-size scratch polynomial.
  temp=n if DAG[i]['op']in['add','sub','mul']else 0
  if live+temp>peak:peak=live+temp;peakid=i
  for j in DAG[i]['args']:
   use[j]-=1
   if not use[j]:live-=DAG[j]['bound']['support_upper']
 reserve=sum(literal.values());DD,NN=DN(S);bound=B(S);raw=(bound[0]+sum(DAG[i]['bound']['support_upper']for i in need if DAG[i]['op']=='atom'))*(len(str(NN))+len(str(DD))+100)
 relation=dict(expression='Saux=D4*Bhat[15]-Bhat[piv4]*Adj(P4)*P4[15]',
  Bhat='rhoLCM*(det2*P_u/rho_u-(N_u,pair*Adj(A2))*(P_pair/rho_pair))',
  hybrid_Schur_entry='rho_u/(det2*rhoLCM*D4)*Saux',
  polynomial_clearing=dict(a_power=max(0,-bound[1][0][0]),B_power=bound[2]),
  full6_chart='FULLU intersectD(Delta5*R434); chosenoneentryneedsno h22 orDelta3inverse. Allcomplementsretained.',
  Mark_final_DIRECT_P4_H5_operator_binding=True,Mark_operator_export=direct['exports']['residual6x13'],Mark_residual_entry=[0,9])
 data=dict(nodes=DAG,output=S,source_to_Mark_hybrid_relation=relation,selected_compact_row=4,current15_row=5,target_original32=15,
  known_rho_LCM_powers=lcm,source_coefficients3=src,source_current15_rows=[5,2,4],Bhat_pivot_and_target_roots=Bu,
  output_bound=b.enc(bound),common_QQ_denominator=str(DD),integer_l1_upper=str(NN),
  last_use_live_numerator_terms_upper=peak,last_use_peak_node=peakid,literal_atom_terms=reserve,
  total_live_with_literal_and_15000reserve=peak+reserve+15000,
  output_integer_coefficient_bits=NN.bit_length(),QQ_denominator_bits=DD.bit_length(),raw_artifact_bytes_upper=raw,
  no_new_coefficient_value_or_parameter_product_evaluated=True)
 (out/'operator-and-cost.json').write_text(json.dumps(data,indent=2)+'\n')
 paths += [Path(__file__),ROOT/f'research/notes/{ID}.md'];report=dict(id=ID,all_passed=True,source_only=True,frozen=True,typed_nodes=len(DAG),live_nodes=len(need),
  output_support_upper=bound[0],output_integer_coefficient_bits=NN.bit_length(),QQ_denominator_bits=DD.bit_length(),
  total_live_with_literal_and_15000reserve=peak+reserve+15000,raw_artifact_bytes_upper=raw,
  resource_gate_passes=(peak+reserve+15000<=150000 and raw<32*1024*1024),
  no_actual_Schur_coefficient_minor_rank_point_or_prime=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},operator_sha256=sha(out/'operator-and-cost.json'),seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
