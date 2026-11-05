#!/usr/bin/env python3
"""Independent remaining55 replay and fixed70 lower certificate; no source solve."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_V2_20261105_v63_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(a):return {int(i):Q(x) for i,x in a if Q(x)}
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 return int(isinstance(x,(int,float,Q)) or isinstance(x,str) and x.lstrip('-').replace('/','').isdigit())
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(exist_ok=False,parents=True);pins={};bits=0
 def read(p,keys=None):pins[p]=sha(ROOT/p);return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 primary='reports/AMBIENT_B18_ELL2_T2_NATIVE70_CONTINUATION55_V2_20261105_v63/';r=read(primary+'report.json');assert pins[primary+'report.json']=='9de0afb961f09ca111d8cccdf2b882299ae0b6e9bd0b75860c63872b554e2f91' and r['all_passed']
 result=read(primary+'actual-native70-complete-lower-images-from-cached15-and-new55.json');combined=result['computed_prefix'];assert len(combined)==70
 oldprefix=read('reports/AMBIENT_B18_ELL2_T2_NATIVE70_VALUES_20261105_v63/actual-native70-lower-prefix-and-witness-if-any.json')['computed_prefix'];assert combined[:15]==oldprefix
 fixed=read(primary+'actual-complete70-lower-fixed-rank-and-span-certificate.json');assert fixed['selected_case619_indices']==[563,609] and fixed['selected_original_column70_indices']==[14,60] and fixed['fixed_minor_row_indices']==[0,1]
 peer='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63_REVIEW/report.json';assert read(peer)['all_passed'] and pins[peer]=='606f7b6e7b30465cae578ed99d5c25ceef74951ce2366823c34cf2f6fcee0671'
 gp='reports/AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63/';d=read(gp+'complete70-cached-input-source-and-positive-whole-image-bounds.json');assert d['whole_source_committed_bound']==128
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
 N={(i,j):{k:Q(x) for k,x in enumerate(v) if Q(x)} for i,j,v in d['actual_needed_N_point_and_first40']}
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
 gamma=list(map(sp,d['Gamma4_source56_times8']));B={i:dict(v) for i,v in d['oldB_source_rows_pruned']};Xp={i:[sp(v) for v in a] for i,a in d['actual_pivot_old8_first_pruned']};E={i:sp(v) for i,v in d['actual_old_pivot_ell_first_pruned']};L0=list(map(sp,d['native_L0_rows_pruned']));Le=list(map(sp,d['native_L1ell_rows_pruned']));S0={i:sp(v) for i,v in d['native_pivot_S0_columns']};Se={i:sp(v) for i,v in d['native_pivot_S1ell_columns']}
 for p in S0:
  v={}
  for mat,s in [(L0,Se[p]),(Le,S0[p])]:
   for j,c in s.items():v=add(v,mat[j],c)
  E[p]=v
 checked=[]
 for rec,provided in zip(d['source70_records'][15:],combined[15:],strict=True):
  assert rec['case619']==provided['case619']
  S=[{} for _ in range(8)]
  for j,c in sp(rec['native128_constant_projection']).items():
   for row,g,sign in defs[j//4]:S[row]=add(S[row],mult({j%4:c},G[g]),sign)
  assert all(not any(k<4 for k in v) for v in S)
  X={j:{} for j in [0,4,5,7]}
  for j in X:
   for i in range(8):X[j]=add(X[j],mult(N[i,j],S[i]))
  old={i:x for i,x in sp(rec['original768_constant_combination']).items() if i<640};a1=sp(rec['actual_a1_27_times9'])
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
  assert list(map(str,vals))==provided['actual_lower4_values'] and pieces==provided['old_native_and_subtracted_pivot_pieces']
  checked.append({'case619':rec['case619'],'actual_lower4_values':list(map(str,vals)),'old_native_and_subtracted_pivot_pieces':pieces})
 # Read the original15; no duplicate image arithmetic is performed.
 allrecords=oldprefix+checked;assert allrecords==combined and len(checked)==55
 D=list(map(int,fixed['fixed_minor_literal_recipe']['integer_scaling_D_byrow']));ints=[]
 for v in allrecords:
  col=[]
  for den,x in zip(D,map(Q,v['actual_lower4_values']),strict=True):
   assert den%x.denominator==0;n=x.numerator*(den//x.denominator);assert abs(n).bit_length()<=128;col.append(n)
  ints.append(col)
 ibits=0
 def ig(x):
  nonlocal ibits
  ibits=max(ibits,abs(x).bit_length());assert ibits<=385;return x
 def det(cols,rows):
  assert cols[0]==14
  ans=0
  for perm in itertools.permutations(range(len(rows))):
   t=(-1)**sum(perm[i]>perm[j] for i in range(len(rows)) for j in range(i+1,len(rows)))
   for i,j in enumerate(perm):t=ig(t*ints[cols[j]][rows[i]])
   ans=ig(ans+t)
  return ans
 v=det([14,60],[0,1]);assert v and (1 if v>0 else -1)==fixed['fixed_scaled_minor_sign'] and abs(v).bit_length()==fixed['fixed_scaled_minor_bits'] and hashlib.sha256(str(v).encode()).hexdigest()==fixed['fixed_scaled_minor_integer_sha256']
 for j,receipt in enumerate(fixed['all70_fixed_span_checks']):
  assert receipt['case619']==549+j and receipt['fixed_bordered_minor_zero_checks']==[True,True]
  assert det([14,60,j],[0,1,2])==0 and det([14,60,j],[0,1,3])==0
 assert fixed['actual_lower4_by70_rank']==2 and fixed['selected_original_QQ_columns']==[allrecords[j]['actual_lower4_values'] for j in [14,60]]
 owned=count([d,table,kap,powers,G,N,S,X,first,sec,jets,gamma,B,a1,Xp,E,L0,Le,op,defs,result,r,combined,checked,ints,fixed])+10000;assert owned<60000
 ep=out/'independent55-source-images-and-all70-fixed-lower-rank.json';e={'all55_new_source_images_and_all12_piece_values_verified':True,'actual55_images':checked,'cached15_data_read_not_recomputed':True,'provided_original_selected_cases':[563,609],'provided_minor_rows':[0,1],'all140_fixed_border_minors_zero':True,'nonzero_selected_minor_sign_bits_and_hash_verified':True,'actual_lower4_by70_rank':2,'actual_quadratic_initial_relation_rank_atleast':2,'actual_HF2_upper':4,'actual_HF2_equality_or_total_length_claimed':False,'no_rank_search_or_source_solve':True,'actual_coefficient_bits':bits,'actual_temporary_integer_bits':ibits};ep.write_text(json.dumps(e,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all55_actual_images_and_fixed70_rank_independently_verified':True,'cached15_not_recomputed':True,'actual_lower_rank':2,'actual_HF2_upper_only':4,'no_higher_source_or_gauge_coefficients':True,'input_sha256':pins,'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'actual_bits':bits,'integer_temporary_bits':ibits,'owned_numeric_bound':owned,'seconds':time.monotonic()-start};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'rank':2,'bits':bits,'integer_bits':ibits,'seconds':report['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
