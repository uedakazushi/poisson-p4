#!/usr/bin/env python3
"""Independent ONE source563 image replay; sparse algebra, no source solve."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_WITNESS563_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(a):return {int(i):Q(x) for i,x in a if Q(x)}
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};bits=0
 def read(p,keys=None):pins[p]=sha(ROOT/p);return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 primary='reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63/';r=read(primary+'report.json');assert pins[primary+'report.json']=='def393dc9f744be2835f744ff68c3f0f509b5fd012c3540e552335a32b6a46f1' and r['all_passed']
 result=read(primary+'actual-native70-lower-prefix-and-witness-if-any.json');w=result['nonzero_witness'];assert w['case619']==563 and r['computed_prefix_count']==15
 peer='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63_REVIEW/report.json';assert read(peer)['all_passed'] and pins[peer]=='606f7b6e7b30465cae578ed99d5c25ceef74951ce2366823c34cf2f6fcee0671'
 gp='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63/';d=read(gp+'complete70-cached-input-source-and-positive-whole-image-bounds.json');assert d['whole_source_committed_bound']==128
 rec=next(x for x in d['source70_records'] if x['case619']==563);assert w['original768_constant_combination']==rec['original768_constant_combination'] and w['actual_a1_27_times9']==rec['actual_a1_27_times9']
 def ck(x):
  nonlocal bits
  bits=max(bits,abs(x.numerator).bit_length(),x.denominator.bit_length());assert bits<=128;return x
 def add(a,b,c=Q(1)):
  z=dict(a)
  for i,x in b.items():
   q=ck(z.get(i,Q(0))+ck(c*x))
   if q:z[i]=q
   else:z.pop(i,None)
  return z
 mons=list(map(tuple,d['A18_basis_monomials']));mi={v:i for i,v in enumerate(mons)};nu=[{i:Q(x) for i,x in enumerate(v) if Q(x)} for v in d['actual_nu_A_B_A18']]
 # Compile monic trace products directly, without the producer's table.
 table={}
 for i in range(72):
  for j in range(72):
   p=mi.get(tuple(sorted(mons[i//4]+mons[j//4])))
   if p is None:continue
   common=(i%4)&(j%4);mask=(i%4)^(j%4);v={4*p+mask:Q(1)}
   for fam in range(2):
    if common&(1<<fam):
     z={}
     for k,c in v.items():
      for b,x in nu[fam].items():
       q=mi.get(tuple(sorted(mons[k//4]+mons[b])))
       if q is not None:z=add(z,{4*q+k%4:ck(c*x)})
     v=z
   if v:table[i,j]=v
 def mult(a,b):
  z={}
  for i,x in a.items():
   for j,y in b.items():
    for k,c in table.get((i,j),{}).items():z=add(z,{k:ck(ck(x*y)*c)})
  return z
 kap=[{i:Q(x) for i,x in enumerate(v) if Q(x)} for v in d['actual_kappa9_C72']];prog=d['regular_G8_source_program'];one={0:Q(1)};powers=[]
 for idx in prog['core_kappa_indices4']:
  p=[one,kap[idx]]
  for _ in range(2,7):p.append(mult(p[-1],kap[idx]))
  assert not mult(p[-1],kap[idx]);powers.append(p)
 def poly(terms):
  z={}
  for exp,n,den in terms:
   if sum(exp)>=7:continue
   v=one
   for i,e in enumerate(exp):
    if e:v=mult(v,powers[i][e])
   z=add(z,v,Q(n,den))
  return z
 def inv(k,p):
  if not k:return one
  z={}
  for e in range(4):z=add(z,powers[1][2*e],Q((-1)**e*math.comb(k+e-1,e),5**(k+e)) if p else Q(math.comb(k+e-1,e)))
  return z
 G=[]
 for name,idx in zip(prog['five_pivot_order'],prog['five_pivot_kappa_indices'],strict=True):
  q=prog['five_actual_rational_roots'][name];G.append(add(kap[idx],mult(mult(poly(q['numerator']),inv(q['D_power'],False)),inv(q['P_power'],True)),Q(-1)))
 G.extend(poly(p) for p in prog['regular_CI_generators3'])
 for g,v in zip(G,d['actual_G8_point_and_first40'],strict=True):assert [g.get(i,Q(0)) for i in range(40)]==list(map(Q,v))
 ell=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json',['actualG8'])['actualG8']
 for g,v in zip(G,ell,strict=True):assert [g.get(4*i+k,Q(0)) for i in [0,8,9,11] for k in range(4)]==list(map(Q,v))
 op=read('reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',['native_S12_by32'])['native_S12_by32'];defs={j:[] for j in range(28)}
 for row,col,gi,sign in op['koszul28_in_first8_rows']:defs[col].append((row,gi,Q(sign)))
 S=[{} for _ in range(8)]
 for j,c in sp(rec['native128_constant_projection']).items():
  for row,g,sign in defs[j//4]:S[row]=add(S[row],mult({j%4:c},G[g]),sign)
 assert all(not any(k<4 for k in v) for v in S)
 N={(i,j):{k:Q(x) for k,x in enumerate(v) if Q(x)} for i,j,v in d['actual_needed_N_point_and_first40']};X={j:{} for j in [0,4,5,7]}
 for j in X:
  for i in range(8):X[j]=add(X[j],mult(N[i,j],S[i]))
 # Ordinary sparse monomial convolution instead of the producer's dense18 jet.
 first={i:sp(v) for i,v in d['actual_C0_leaf_first10']};sec={i:list(map(Q,v)) for i,v in d['actual_13_honest_second8']};ops=dict(d['typed_closure_nodes']);jets={}
 def pmul(a,b):
  z={}
  for i,x in a.items():
   for j,y in b.items():
    k=mi.get(tuple(sorted(mons[i]+mons[j])))
    if k is not None:z=add(z,{k:ck(x*y)})
  return z
 for k,n in sorted(ops.items()):
  if n[0]=='actual_C0_scalar_function':
   v=dict(first[n[1]])
   for j,x in enumerate(sec.get(n[1],[Q(0)]*8)):
    if x:v[10+j]=x
  elif n[0]=='c':v={0:Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='sum':
   v={}
   for child,a,b in n[1]:v=add(v,jets[child],Q(a,b))
  else:
   assert n[0]=='product';v={0:Q(n[1],n[2])}
   for child in n[3]:v=pmul(v,jets[child])
  jets[k]=v
 gamma=list(map(sp,d['Gamma4_source56_times8']));B={i:dict(v) for i,v in d['oldB_source_rows_pruned']};old={i:x for i,x in sp(rec['original768_constant_combination']).items() if i<640};a1=sp(rec['actual_a1_27_times9']);Xp={i:[sp(v) for v in a] for i,a in d['actual_pivot_old8_first_pruned']};E={i:sp(v) for i,v in d['actual_old_pivot_ell_first_pruned']};L0=list(map(sp,d['native_L0_rows_pruned']));Le=list(map(sp,d['native_L1ell_rows_pruned']));S0={i:sp(v) for i,v in d['native_pivot_S0_columns']};Se={i:sp(v) for i,v in d['native_pivot_S1ell_columns']}
 for p in S0:
  v={}
  for mat,s in [(L0,Se[p]),(Le,S0[p])]:
   for j,c in s.items():v=add(v,mat[j],c)
  E[p]=v
 vals=[];pieces=[];pairs=list(map(tuple,d['ordinary_pair_order8']))
 for row in gamma:
  ov=nv=cv=Q(0)
  for slot,c in row.items():
   source,pair=divmod(slot,8);u,v=pairs[pair];g,mask=divmod(source,4);nv=ck(nv+ck(c*X[g].get(4*(10+pair)+mask,Q(0))))
   for j,c0 in old.items():
    k=B[j].get(source)
    if k is not None:ov=ck(ov+ck(c*ck(c0*jets[k].get(10+pair,Q(0)))))
   for p in d['needed_pivot_indices']:
    au,av=a1.get(9*p+u,Q(0)),a1.get(9*p+v,Q(0))
    if au:cv=ck(cv+ck(c*ck(au*(Xp[p][v].get(source,Q(0)) if v<8 else E[p].get(source,Q(0))))))
    if u!=v and av:cv=ck(cv+ck(c*ck(av*(Xp[p][u].get(source,Q(0)) if u<8 else E[p].get(source,Q(0))))))
  vals.append(ck(ck(ov+nv)-cv));pieces.append(list(map(str,[ov,nv,cv])))
 assert list(map(str,vals))==w['actual_lower4_values']==['-1000/3','-212/3','-19000','-3968'] and pieces==w['old_native_and_subtracted_pivot_pieces'] and any(vals)
 owned=count([d,table,kap,powers,G,N,S,X,first,sec,jets,gamma,B,a1,Xp,E,L0,Le,op,defs,w,result,r])+10000;assert owned<60000
 ep=out/'independent-one563-full-source-and-four-piece-witness.json';e={'case619':563,'original768_constant_combination':rec['original768_constant_combination'],'actual_a1_27_times9':rec['actual_a1_27_times9'],'independent_old_native_subtracted_pivot_pieces':pieces,'actual_lower4_values':list(map(str,vals)),'nonzero_actual_quadratic_initial_relation':True,'all70_values_or_HF2_or_total_length_claimed':False,'no_new_source_gauge_implicit_solve':True,'actual_bits':bits,'source_bound':128,'owned_numeric_bound':owned};ep.write_text(json.dumps(e,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'one563_lower_witness_independently_verified':True,'four_lower_values':list(map(str,vals)),'no_other70_images_queried':True,'actual_HF2_or_total_length_computed':False,'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'actual_bits':bits,'owned_numeric_bound':owned,'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'values':list(map(str,vals)),'bits':bits,'seconds':report['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
