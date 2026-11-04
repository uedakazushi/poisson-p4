#!/usr/bin/env python3
"""Independent sparse C36 source/finite-gauge certificate verification."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import product
import hashlib,json,resource,time
import ambient_B18_comparison_X_origin_review_20261104_v63 as V
import ambient_B18_comparison_X_first640_linear_review_20261104_v63 as D9

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_NATIVE128_FIRST_20261104_v63';STEM=AUTHOR+'_REVIEW'

class Ring:
 def __init__(self,nus):self.nus=nus;self.Z={};self.ONE={(0,0):F(1)};self.inverse_count=0
 def c(self,q):return{(0,0):F(q)}if q else{}
 def add(self,*xs):
  out={}
  for x in xs:
   for k,c in x.items():
    q=out.get(k,F(0))+c
    if q:out[k]=q
    else:out.pop(k,None)
  return out
 def neg(self,x):return{k:-c for k,c in x.items()}
 def scale(self,x,q):return{k:c*q for k,c in x.items()if c*q}
 def mul(self,a,b):
  out={}
  for(e,m),x in a.items():
   for(f,n),y in b.items():
    if e and f:continue
    carry=m&n;mask=m^n
    if not carry:contrib=[((e+f,mask),x*y)]
    elif e or f or carry==3:continue
    else:contrib=[((i,mask),x*y*c)for i,c in enumerate(self.nus[0 if carry==1 else 1])if c]
    for k,c in contrib:
     q=out.get(k,F(0))+c
     if q:out[k]=q
     else:out.pop(k,None)
  return out
 def tm(self,*xs):
  out=self.ONE
  for x in xs:out=self.mul(out,x)
  return out
 def power(self,x,n):
  if n<0:return self.power(self.inv(x),-n)
  out=self.ONE
  while n:
   if n&1:out=self.mul(out,x)
   x=self.mul(x,x);n//=2
  return out
 def inv(self,x):
  c=x.get((0,0),F(0));assert c,'nonunit inversion prohibited'
  nil=self.scale(self.add(x,self.c(-c)),-1/c);term=self.ONE;ans=self.ONE
  for _ in range(4):term=self.mul(term,nil);ans=self.add(ans,term)
  assert not self.mul(term,nil)
  ans=self.scale(ans,1/c);assert self.mul(x,ans)==self.ONE;self.inverse_count+=1;return ans
 def dual(self,raw):return{(i,0):F(c)for i,c in enumerate(raw)if F(c)}
 def tau(self,i):return{(0,1<<i):F(1)}
 def flat(self,raw):return{divmod(i,4):F(c)for i,c in enumerate(raw)if F(c)}
 def install(self):
  V.Z={};V.ONE=self.ONE;V.TA=self.tau(0);V.TB=self.tau(1)
  for name,f in [('scalar',self.c),('plus',self.add),('neg',self.neg),('sc',self.scale),('times',self.mul),('tm',self.tm),('inverse',self.inv),('power',self.power)]:setattr(V,name,f)
  def mi(A):
   n=len(A);base=[[x.get((0,0),F(0))for x in row]for row in A]
   C=V.constmat(V.rational_inverse(base));E=V.matmul(C,[[self.add(x,self.c(-x.get((0,0),F(0))))for x in row]for row in A])
   term=V.eye(n);S=V.eye(n)
   for k in range(1,5):
    term=V.matmul(term,E)
    S=[[self.add(S[i][j],self.scale(term[i][j],(-1)**k))for j in range(n)]for i in range(n)]
   assert all(not x for row in V.matmul(term,E)for x in row)
   answer=V.matmul(S,C);assert V.matmul(A,answer)==V.matmul(answer,A)==V.eye(n);return answer
  V.matrix_inverse=mi

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def unpack(A,ring):return[[ring.flat(x)for x in row]for row in A]
def mmQ(A,B):return[[sum((x*y for x,y in zip(row,col)),F(0))for col in zip(*B)]for row in A]

def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 rp='reports/'+AUTHOR+'/report.json';R=read(rp);assert R['all_passed'];hashes={rp:sha(rp)}
 for p,h in R['input_sha256'].items():assert sha(p)==h,p;hashes[p]=h
 for e in R['phase_evidence']+[R['data_evidence']]:assert sha(e['path'])==e['sha256'];hashes[e['path']]=e['sha256']
 code='code/ambient_B18_native128_first_20261104_v63.py';assert sha(code)==R['code_sha256'];hashes[code]=sha(code)
 for st in ['AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW','AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW','AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW']:
  p='reports/'+st+'/report.json';a=read(p);assert a['all_passed'];hashes[p]=sha(p)
  for kind in ['input_sha256','checked_dependency_sha256']:
   for p,h in a.get(kind,{}).items():
    if p in R['input_sha256']:assert R['input_sha256'][p]==h
 def src(st,name):return read('reports/'+st+'_20261104_v63/'+name)
 P=src('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT','literal-proper-source14-preflight.json')
 I=src('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET','actual-source-physical-firstjet-interface.json')
 S=src('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT','actual-source-coefficients.json')
 A20=src('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION','actual-source105-Artin20-values.json')
 PF=src('AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT','actual-native128-first-input-and-mixed-index-profile.json')
 Q=src('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL','actual-proper-source14-equations17.json')
 H=src('AMBIENT_B18_PROPER_RATIONAL_ELIMINATION','actual12-flag-divided-polynomials.json')
 K=src('AMBIENT_B18_PROPER_CORE_ELIMINATION','actual4variable7core-polynomials.json')
 M=src('AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP','membership-and-source-cofactors.json')
 prep=src('AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP','preparation.json')
 Op=src('AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR','actual-complete-Res4-syzygy-operator.json')
 E=src('AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL','actual-normalizer-tensor-free4-DAG.json')
 origin=src('AMBIENT_B18_COMPARISON_X_ORIGIN','phase2-native-adapter.json')
 X0=src('AMBIENT_B18_COMPARISON_X_ORIGIN','actual-X0-and-unit-cancellation.json')
 phases=[read(e['path'])for e in R['phase_evidence']];p1,p2,p3,p4,p5=phases;final=read(R['data_evidence']['path'])
 for name in ['code/ambient_B18_comparison_X_origin_review_20261104_v63.py','code/ambient_B18_comparison_X_first640_linear_review_20261104_v63.py']:
  hashes[name]=sha(name)
 first={i:tuple(map(F,p))for i,p in PF['all105_existing_unknown_first9']}
 basis=A20['A20_basis'];positions=[0]+[basis.index([int(i==j)for i in range(8)])for j in range(8)]
 rawvalues=dict(zip(A20['joint105_unknown_order'],A20['all105_unknown_final_values']))
 assert len(rawvalues)==105 and set(first)==set(rawvalues)
 for node,jet in first.items():
  raw={i:F(c)for i,c in rawvalues[node]};assert jet==tuple(raw.get(i,F(0))for i in positions)
 nus=[list(map(F,p))for p in PF['nu_A_B_first9']];assert all(not p[0]for p in nus)
 for nu,raw in zip(nus,A20['nu_CA_CB_values']):
  values={i:F(c)for i,c in raw};assert nu==[values.get(i,F(0))for i in positions]
 roots=[node for family in S['families'].values()for pair in family['complete_normalized_physical_w_coefficients'][0]for node in pair]
 dag=S['DAG_scalar_nodes'];needed=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);n=dag[i]
  if n[0]in['a','m','n','i']:todo.extend(n[1:])
 order=sorted(needed);assert order==PF['dual9_source_evaluation']['closure_order']and len(order)==17926
 (out/'preflight.json').write_text(json.dumps({'source_DAG_closure':len(order),'fixed_unknown105':len(first),'native_dimension':36,
  'no_new_source_or_implicit_solve':True,'input_sha256':hashes},indent=2)+'\n')
 uses=Counter()
 for i in order:
  if dag[i][0]in['a','m','n','i']:uses.update(dag[i][1:])
 uses.update(roots);vals={};peak=0
 def dual(raw):return{i:F(c)for i,c in enumerate(raw)if F(c)}
 for i in order:
  n=dag[i];tag=n[0]
  if tag=='v':value=dual(first[i])
  elif tag=='c':value=D9.const(F(n[1],n[2]))
  elif tag=='a':value=D9.plus(vals[n[1]],vals[n[2]])
  elif tag=='m':value=D9.times(vals[n[1]],vals[n[2]])
  elif tag=='n':value={k:-q for k,q in vals[n[1]].items()}
  elif tag=='i':
   p=vals[n[1]];c=p.get(0,F(0));assert c;value={0:1/c};value.update({k:-q/c**2 for k,q in p.items()if k});assert D9.times(p,value)=={0:F(1)}
  else:raise AssertionError(n)
  vals[i]=value;peak=max(peak,9*len(vals));assert peak<150000
  if tag in['a','m','n','i']:
   for k in n[1:]:
    uses[k]-=1
    if not uses[k]:del vals[k]
 pairs={f:[[tuple(vals[a].get(i,F(0))for i in range(9)),tuple(vals[b].get(i,F(0))for i in range(9))]for a,b in spec['complete_normalized_physical_w_coefficients'][0]]for f,spec in S['families'].items()}
 assert pairs=={f:[[tuple(map(F,a)),tuple(map(F,b))]for a,b in spec]for f,spec in p1['physical150_tracepairs_first9'].items()}
 assert nus==[[F(x)for x in row]for row in p1['nu_first9']]
 positionsG=I['actual_group15_positions'];GL=PF['GL18_source_parameter_positions'];odd=list(map(F,Op['group30_to120']['group30_odd_origin']))
 rows=P['zero_coordinate_rows36'];selected=P['CA_constant_left_inverse_rows19'];fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices'];center=list(map(F,P['actual_CA_source19_center']))
 def frame_family(ring,family):
  spec=S['families'][family];tau=ring.tau(0 if family=='CA'else 1)
  p=ring.add(tau,ring.scale(ring.dual(first[spec['p_quadratic_relation']['s']]),F(1,2)))
  G=V.eye(5)
  for name,(aa,bb)in spec['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);G[i][j]=ring.add(G[i][j],ring.dual(first[aa]),ring.mul(ring.dual(first[bb]),p))
  C=V.eye(5);C[0][3]=ring.scale(p,-2 if family=='CA'else F(1,2));C[1][0]=p;C[1][3]=ring.scale(ring.mul(p,p),-1 if family=='CA'else F(1,4))
  return V.matmul(G,C)
 def physical(ring,family):return[ring.add(ring.dual(a),ring.mul(ring.dual(b),ring.tau(0 if family=='CA'else 1)))for a,b in pairs[family]]
 def gauge(ring,ev,od):
  G=V.eye(5)
  for(i,j),a,b in zip(positionsG,ev,od):G[i][j]=ring.add(G[i][j],ring.dual(a),ring.mul(ring.dual(b),ring.tau(0)))
  return G
 def constraints(ring,G,HA,PA):
  Fm=V.normalized_flag_frame(V.matmul(G,HA));moved=V.tensor_push(PA,V.matmul(Fm,G))
  own=[row[0]for row in V.matmul(V.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in selected])]
  return[Fm[1][0],Fm[1][3],Fm[1][4],Fm[0][3],Fm[0][4]]+[ring.add(own[j],ring.c(-center[j]))for j in fixed]
 # Literal mixed225 replay in length4. This verifies a fixed finite source
 # coefficient certificate, not a new physical jet or a source solve.
 r4=Ring([[F(0)],[F(0)]]);r4.install();HA=V.eye(5)
 sourceMA=next(f for f in I['families']if f['family']=='CA')['source_frame5_by5_by9_first'][7]
 for i,j in product(range(5),repeat=2):HA[i][j]=r4.add(HA[i][j],r4.scale(r4.tau(0),F(sourceMA[i][j])))
 base=V.eye(5)
 for(i,j),c in zip(positionsG,odd):base[i][j]=r4.add(base[i][j],r4.scale(r4.tau(0),c))
 A0=[[F(0)]*15 for _ in range(15)];A1=[[F(0)]*15 for _ in range(15)]
 for k,(i,j)in enumerate(positionsG):
  G=[row[:]for row in base];G[i][j]=r4.add(G[i][j],r4.tau(1))
  actual=constraints(r4,G,HA,[r4.c(x)for x in P['actual_P0_physical150']])
  for row,p in enumerate(actual):
   assert not p.get((0,0),F(0))and not p.get((0,1),F(0))
   A0[row][k]=p.get((0,2),F(0));A1[row][k]=p.get((0,3),F(0))
 assert A0==[[F(x)for x in row]for row in p2['A0']]==[[F(x)for x in row]for row in I['group15_unit_Jacobian']]
 assert A1==[[F(x)for x in row]for row in p2['A1_eta_delta']]
 JJ=[[F(0)]*30 for _ in range(30)]
 for i,j in product(range(15),repeat=2):JJ[i][j]=A0[i][j];JJ[15+i][j]=A1[i][j];JJ[15+i][15+j]=A0[i][j]
 Ji=[[F(x)for x in row]for row in p2['J30_inverse']];identity=[[F(i==j)for j in range(30)]for i in range(30)]
 assert mmQ(JJ,Ji)==mmQ(Ji,JJ)==identity and F(p2['det_J30'])==F(4,9)
 ring=Ring(nus);ring.install();HA=frame_family(ring,'CA');PA=physical(ring,'CA')
 ev0=[(F(0),)*9 for _ in range(15)];od0=[(c,)+(F(0),)*8 for c in odd]
 residual=constraints(ring,gauge(ring,ev0,od0),HA,PA)
 rhs=[[p.get((k+1,b),F(0))for k in range(8)]for b in range(2)for p in residual]
 assert rhs==[[F(x)for x in row]for row in p3['RHS30_by8']]
 correction=[[F(x)for x in row]for row in p3['correction30_by8']]
 assert mmQ(JJ,correction)==[[-x for x in row]for row in rhs]
 ev=[tuple(map(F,row))for row in p3['even15_first9']];od=[tuple(map(F,row))for row in p3['odd15_first9']]
 assert ev==[(F(0),)+tuple(row)for row in correction[:15]]
 assert od==[(odd[i],)+tuple(row)for i,row in enumerate(correction[15:])]
 G=gauge(ring,ev,od);assert all(not p for p in constraints(ring,G,HA,PA))
 assert V.normalized_flag_frame(V.matmul(G,HA))==V.eye(5)
 FB=V.normalized_flag_frame(V.matmul(G,frame_family(ring,'CB')))
 moved=V.tensor_push(PA,G);own=[row[0]for row in V.matmul(V.constmat(P['CA_left_inverse19']),[[moved[rows[j]]]for j in selected])]
 assert all(own[j]==ring.c(center[j])for j in fixed)
 actual=[ring.add(own[j],ring.c(-center[j]))for j in free]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]]
 kappa=[actual[j]for j in I['retained_source9_indices']]
 assert actual==[ring.flat(p)for p in p4['actual14']]and kappa==[ring.flat(p)for p in p4['kappa9']]
 targets=[18,19,28,29,60,67,72,63,64,68,73,74];PB=physical(ring,'CB')
 F12=[ring.add(PA[j],ring.neg(PB[j]))for j in targets];assert F12==[ring.flat(p)for p in p4['physicalF12']]and any(F12)
 R17=[[{}]*12 for _ in range(17)];transform=V.matmul(FB,G);ann=V.constmat(P['CB_annihilator17_by36'])
 for j,target in enumerate(targets):
  inp=[{}]*150;inp[target]=ring.ONE
  for k,row in enumerate(E['eliminated_zero_rows4']):inp[row]=ring.c(E['eliminated_four_coordinate_reconstruction'][k][rows.index(target)])
  pushed=V.tensor_push(inp,transform);column=V.matmul(ann,[[pushed[i]]for i in rows])
  for i in range(17):R17[i][j]=column[i][0]
 assert R17==unpack(p4['R17'],ring)
 F17=[[V.polynomial(p,actual)]for p in Q['equations14variable_sparse']];assert V.matmul(R17,[[p]for p in F12])==F17
 D=ring.add(ring.ONE,ring.neg(ring.mul(kappa[1],kappa[1])));PP=ring.add(ring.c(5),ring.mul(kappa[1],kappa[1]))
 axes=[0,1,10,11,13];firstroots=[]
 for name in['v0','v1','b','c','e']:
  p=H['five_explicit_rational_roots'][name];firstroots.append(ring.tm(V.polynomial(p['numerator'],kappa),ring.power(D,-p['D_denominator_power'])))
 co=[V.synthetic_prefix_coefficients(p,actual,axes,firstroots)for p in Q['equations14variable_sparse']]
 BA=V.matmul(co[5:],V.matrix_inverse(co[:5]));R12=[[ring.add(R17[i+5][j],ring.neg(ring.add(*(ring.mul(BA[i][k],R17[k][j])for k in range(5)))))for j in range(12)]for i in range(12)]
 for i,p in enumerate(H['all12_actual_remaining_and_divided_records']):R12[i]=[ring.mul(ring.scale(ring.power(D,p['D_denominator_power']),F(p['constant_denominator'],p['primitive_content_unit'])),x)for x in R12[i]]
 assert R12==unpack(p4['R12'],ring)
 hpolys=[p['divided_primitive_polynomial']for p in H['all12_actual_remaining_and_divided_records']];hv=[[V.polynomial(p,kappa)]for p in hpolys]
 assert V.matmul(R12,[[p]for p in F12])==[[ring.mul(kappa[1],p[0])]for p in hv]
 core=[kappa[j]for j in[0,1,2,6]];axes=[K['source9_variable_order'].index(n)for n in K['selected5_pivot_order']];roots2=[]
 for name in K['selected5_pivot_order']:
  p=K['five_actual_rational_roots'][name];roots2.append(ring.tm(V.polynomial(p['numerator'],core),ring.power(D,-p['D_power']),ring.power(PP,-p['P_power'])))
 co=[V.synthetic_prefix_coefficients(p,kappa,axes,roots2)for p in hpolys];perm=K['selected5_source_H_rows']+[i for i in range(12)if i not in K['selected5_source_H_rows']]
 Ai=V.matrix_inverse([co[i]for i in perm[:5]]);BA=V.matmul([co[i]for i in perm[5:]],Ai);T=[[{}]*12 for _ in range(12)]
 for i,j in product(range(5),repeat=2):T[i][perm[j]]=Ai[i][j]
 for i in range(7):
  T[i+5][perm[i+5]]=ring.ONE
  for j in range(5):T[i+5][perm[j]]=ring.neg(BA[i][j])
 for i,p in enumerate(K['remaining7_actual_core_polynomials']):
  lam=ring.scale(ring.tm(ring.power(D,-p['D_denominator_power']),ring.power(PP,-p['P_denominator_power'])),F(p['QQ_unit_content_numerator'],p['QQ_unit_content_denominator']))
  T[i+5]=[ring.mul(ring.inv(lam),x)for x in T[i+5]]
 r=ring.add(core[2],ring.neg(core[0]),ring.neg(core[1]));E3=V.eye(12);E3[5][5]=ring.inv(D);E3[6][6]=ring.inv(D)
 E3[7][5]=ring.neg(ring.tm(r,ring.inv(ring.mul(PP,D))));E3[7][7]=ring.scale(ring.inv(ring.mul(PP,D)),-2);T=V.matmul(E3,T);E4=V.eye(12)
 for j in range(4):
  for i in range(3):E4[8+j][5+i]=ring.neg(ring.tm(ring.power(D,[2,2,2,3][j]),V.polynomial(M['cofactor3_by4'][i][j],core)))
 T=V.matmul(E4,T);N=V.matmul(T,R12);G8=[ring.add(kappa[j],ring.neg(v))for j,v in zip(axes,roots2)]+[V.polynomial(p,core)for p in prep['generators3']]
 assert T==unpack(p5['T'],ring)and N==unpack(p5['N'],ring)and G8==[ring.flat(p)for p in p5['G8']]
 assert V.matmul(T,hv)==[[p]for p in G8]+[[{}]for _ in range(4)]
 assert V.matmul(N,[[p]for p in F12])==[[ring.mul(kappa[1],p)]for p in G8]+[[{}]for _ in range(4)]
 for name,A in [('R12',R12),('T',T),('N',N)]:assert[[[str(x.get((0,k),F(0)))for k in range(4)]for x in row]for row in A]==origin[name]
 assert[[str(x.get((0,k),F(0)))for k in range(4)]for x in G8]==origin['G8_C0']
 Sn=[[{}]*32 for _ in range(12)]
 for i,j,g,sign in Op['native_S12_by32']['koszul28_in_first8_rows']:Sn[i][j]=ring.scale(G8[g],F(sign))
 for i,j,c in Op['native_S12_by32']['free4_zero_source_columns']:Sn[i][j]=ring.c(c)
 native=V.matmul(list(map(list,zip(*N))),Sn);output=[];nnz=0
 for j in range(128):
  col,b=divmod(j,4);entries=[]
  for i in range(12):
   value=ring.mul(native[i][col],{(0,b):F(1)})
   for mask in range(4):
    row=4*i+mask;terms=[[e,str(value.get((e,mask),F(0)))]for e in range(9)if value.get((e,mask),F(0))]
    assert value.get((0,mask),F(0))==F(X0['all48_by128_SyzQ_origins'][row][j])
    if terms:entries.append([row,terms]);nnz+=len(terms)
  output.append(entries)
 assert output==final['source_columns128_rows48_sparse_dual9']
 proof={'all_passed':True,'all105_firsts_original_Artin_values_bound':True,'all17926_source_DAG_nodes_dual9_replayed':True,
  'peak_source9_slots':peak,'all225_mixed_gauge_coefficients_and240_RHS_checked':True,
  'saved30_fixed_inverse_both_products_and240_correction_equations':True,
  'saved_gauge_checked_in_actual_C36_all30_constraints_zero':True,
  'full_F12_nonzero_and_all17_offintersection_identities_checked':True,
  'full_first5_ordered_telescopes_not_origin_derivatives':True,
  'all_R12_T_N_G8_and_dH_TH_NF_identities_checked':True,
  'all_constant_N_G8_S_entries_match_accepted_C4_peer':True,
  'all128_Res4_NtransposeS_columns_and6144_origins_checked':True,'output_nonzero_coefficients':nnz,
  'C36_square_reduction_retains_tau3_contributions':True,
  'no_new_source_implicit_solve_or_physical_jet_GB_ranksearch_muT_depth':True}
 (out/'independent-source-gauge-offshell-and-native128-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,'actual_native128_X1_independently_source_verified':True,
  'all_C36_source_gauge_R12_T_N_and_origin_interfaces_checked':True,'all128_X1_columns_checked':True,
  'X1_rank_muT_length_K_torsion_or_depth_claim':False,'input_sha256':hashes,
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
