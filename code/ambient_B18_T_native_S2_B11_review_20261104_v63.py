#!/usr/bin/env python3
"""Independent B11/C44 source gate and restricted effective-S2 audit."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
import hashlib,json,resource,time
import ambient_B18_comparison_X_origin_review_20261104_v63 as V
from ambient_B18_native128_first_review_20261104_v63 import Ring as DualRing

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63';STEM=AUTHOR+'_REVIEW'
SURVIVE=[0,1,4,5]
BASIS=[(0,)*8,(0,)*7+(1,),(0,)*7+(2,)]+[tuple(int(j==i)for j in range(8))for i in SURVIVE]+[tuple(int(j==i or j==7)for j in range(8))for i in SURVIVE]
INDEX={m:i for i,m in enumerate(BASIS)};LINEAR=[(0,0),(1,8),(3,1),(4,2),(5,5),(6,6)];QUAD=[2,7,8,9,10]
def sp(raw):return{int(i):F(q)for i,q in raw if F(q)}
def add(a,b,c=1):
 a=dict(a)
 for i,q in b.items():
  x=a.get(i,F(0))+c*q
  if x:a[i]=x
  else:a.pop(i,None)
 return a
def bm(a,b):
 out={}
 for i,c in a.items():
  for j,d in b.items():
   k=INDEX.get(tuple(x+y for x,y in zip(BASIS[i],BASIS[j])))
   if k is not None:out=add(out,{k:c*d})
 return out
def embed_first(raw):return{bi:F(raw[di])for bi,di in LINEAR if F(raw[di])}

class B44:
 def __init__(self,nus):self.nus=nus;self.ONE={(0,0):F(1)};self.inverse_count=0
 def c(self,x):return{(0,0):F(x)}if x else{}
 def add(self,*xs):
  out={}
  for x in xs:
   for k,q in x.items():
    v=out.get(k,F(0))+q
    if v:out[k]=v
    else:out.pop(k,None)
  return out
 def neg(self,p):return{k:-q for k,q in p.items()}
 def scale(self,p,q):return{k:v*q for k,v in p.items()if v*q}
 def mul(self,a,b):
  out={}
  for(i,u),c in a.items():
   for(j,v),d in b.items():
    part=bm({i:c},{j:d});carry=u&v;mask=u^v
    if carry&1:part=bm(part,self.nus[0])
    if carry&2:part=bm(part,self.nus[1])
    out=self.add(out,{(k,mask):q for k,q in part.items()})
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
  c=x.get((0,0),F(0));assert c
  nil=self.scale(self.add(x,self.c(-c)),-1/c);term=self.ONE;value=self.ONE
  for _ in range(8):term=self.mul(term,nil);value=self.add(value,term)
  assert not self.mul(term,nil);value=self.scale(value,1/c);assert self.mul(x,value)==self.ONE;self.inverse_count+=1;return value
 def base(self,p):return{(i,0):q for i,q in p.items()}
 def tau(self,k):return{(0,1<<k):F(1)}
 def flat(self,p):return{divmod(i,4):F(c)for i,c in enumerate(p)if F(c)}
 def install(self):
  V.Z={};V.ONE=self.ONE;V.TA=self.tau(0);V.TB=self.tau(1)
  for n,f in[('scalar',self.c),('plus',self.add),('neg',self.neg),('sc',self.scale),('times',self.mul),('tm',self.tm),('inverse',self.inv),('power',self.power)]:setattr(V,n,f)
  def mi(A):
   n=len(A);C=V.constmat(V.rational_inverse([[x.get((0,0),F(0))for x in row]for row in A]))
   E=V.matmul(C,[[self.add(x,self.c(-x.get((0,0),F(0))))for x in row]for row in A]);term=V.eye(n);value=V.eye(n)
   for k in range(1,9):
    term=V.matmul(term,E);value=[[self.add(value[i][j],self.scale(term[i][j],(-1)**k))for j in range(n)]for i in range(n)]
   assert all(not x for row in V.matmul(term,E)for x in row)
   answer=V.matmul(value,C);assert V.matmul(A,answer)==V.matmul(answer,A)==V.eye(n);return answer
  V.matrix_inverse=mi

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def mmQ(A,B):return[[sum((x*y for x,y in zip(r,c)),F(0))for c in zip(*B)]for r in A]

def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 rp='reports/'+AUTHOR+'/report.json';R=read(rp);assert R['all_passed'];hashes={rp:sha(rp)}
 for p,h in R['input_hashes'].items():assert sha(p)==h,p;hashes[p]=h
 for e in R['phase_evidence']+[R['data_evidence']]:assert sha(e['path'])==e['sha256'];hashes[e['path']]=e['sha256']
 code='code/ambient_B18_T_native_S2_B11_20261104_v63.py';assert sha(code)==R['code_sha256'];hashes[code]=sha(code)
 for st in['AMBIENT_B18_NATIVE128_FIRST_20261104_v63_REVIEW','AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW','AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW']:
  p='reports/'+st+'/report.json';d=read(p);assert d['all_passed'];hashes[p]=sha(p)
 def get(st,f):
  p='reports/'+st+'_20261104_v63/'+f;hashes[p]=sha(p);return read(p)
 H=get('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20','actual-special-fiber20-coefficients.json')
 Art=get('AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION','actual-source105-Artin20-values.json')
 B=get('AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT','actual-B11-source-and-gauge-preflight.json')
 S=get('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT','actual-source-coefficients.json')
 I=get('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET','actual-source-physical-firstjet-interface.json')
 P=get('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT','literal-proper-source14-preflight.json')
 K=get('AMBIENT_B18_PROPER_CORE_ELIMINATION','actual4variable7core-polynomials.json')
 Prep=get('AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP','preparation.json')
 Op=get('AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR','actual-complete-Res4-syzygy-operator.json')
 PF=get('AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT','actual-native128-first-input-and-mixed-index-profile.json')
 Old3=get('AMBIENT_B18_NATIVE128_FIRST','phase3-gauge30-firsts.json');Old2=get('AMBIENT_B18_NATIVE128_FIRST','phase2-fixed30-Jacobian.json')
 Old4=get('AMBIENT_B18_NATIVE128_FIRST','phase4-native-source-and-R12.json');Old5=get('AMBIENT_B18_NATIVE128_FIRST','phase5-native-adapter-firsts.json')
 old_array=get('AMBIENT_B18_NATIVE128_FIRST','actual-native128-X-mod-m2.json')
 M1=get('AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT','actual-projected-three-class-X2-operator-preflight.json')
 L2=get('AMBIENT_B18_T_NATIVE_ADAPTER2_CANCELLATION_PREFLIGHT','actual-native-L2-cancellation-source-preflight.json')
 for st in['AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT','AMBIENT_B18_T_NATIVE_ADAPTER2_CANCELLATION_PREFLIGHT','AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT']:
  report=get(st,'report.json');assert report['all_passed']
  for val in report['input_hashes'].values():assert sha(val['path'])==val['sha256'];hashes[val['path']]=val['sha256']
  for name,h in report['evidence_sha256'].items():
   p='reports/'+st+'_20261104_v63/'+name;assert sha(p)==h;hashes[p]=h
 X0=get('AMBIENT_B18_COMPARISON_X_ORIGIN','actual-X0-and-unit-cancellation.json')
 X1=get('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR','actual-first640-X-mod-m2.json')['source_columns640_rows56_sparse_dual9']+old_array['source_columns128_rows48_sparse_dual9']
 low=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC','phase2-certified-quadratic-input.json')
 lowfull=get('AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC','actual-linear-quadratic-T-certificate.json')
 schur=get('AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR','actual-minimal5-by825-Schur-operator.json')
 phases=[read(e['path'])for e in R['phase_evidence']];p1,p2,p3=phases;final=read(R['data_evidence']['path'])
 assert[list(m)for m in BASIS]==B['B11_basis_monomials8']
 hb=H['quotient_basis_monomials'];maps=[{INDEX[tuple(m)]:F(1)}if tuple(m)in INDEX else{}for m in hb]
 assert maps==[sp(p)for p in B['H20_basis_to_B11_columns']]
 def project(v):
  ans={}
  for i,c in v.items():ans=add(ans,maps[i],c)
  return ans
 actions=[[sp(p)for p in rows]for rows in H['multiplication8_column_sparse']]
 def hp(i,j):
  val={i:F(1)}
  for k,n in enumerate(hb[j]):
   for _ in range(n):
    outval={}
    for a,c in val.items():outval=add(outval,actions[k][a],c)
    val=outval
  return val
 for i,j in product(range(20),repeat=2):assert project(hp(i,j))==bm(maps[i],maps[j])
 assert B['B11_basis_lifts_in_H20']==[hb.index(list(m))for m in BASIS]
 first={node:project(sp(v))for node,v in zip(Art['joint105_unknown_order'],Art['all105_unknown_final_values'])}
 assert first=={i:sp(v)for i,v in B['all105_actual_B11_root_values']}
 nus=[project(sp(v))for v in Art['nu_CA_CB_values']];assert nus==[sp(v)for v in B['nu_A_B_actual_B11']]
 assert all(not bm(a,b)for a,b in product(nus,repeat=2))
 ring=B44(nus)
 # Direct source quotient construction checks the full44x44 table; its
 # associativity follows from the standard monomial quotient law.
 for i,j in product(range(44),repeat=2):
  assert ring.mul({divmod(i,4):F(1)},{divmod(j,4):F(1)})=={divmod(k,4):q for k,q in sp(B['C_B44_multiplication_basis_products'][i][j]).items()}
 (out/'preflight.json').write_text(json.dumps({'H20_products':400,'C44_basis_products':1936,'source105_map':105,'effective_only_no_individual_X2':True,'input_sha256':hashes},indent=2)+'\n')
 # Source-specific L2 cancellation: all actual zero-origin combinations
 # vanish SEPARATELY in the old640 and native128 origin blocks.
 rest=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations']if r.get('zero')]
 combos=[sp(r['constant_original_source_combination'])for r in rest]+[sp(c)for c in M1['620_actual_source_combinations_in_original768']]
 assert len(combos)==1361
 for c in combos:
  for row in X0['full_X0_56_by768']:
   assert sum((F(row[j])*q for j,q in c.items()if j<640),F(0))==0
   assert sum((F(row[j])*q for j,q in c.items()if j>=640),F(0))==0
 assert L2['actual_old_B_origin_rank']==4 and L2['actual_joint_origin_rank']==27 and L2['rank_native_origin_proof']['sum']==23
 oldp=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations']if'pivot'in r and r['source_column']<640]
 assert[r['source_column']for r in oldp]==[604,613,622,631]
 assert all(all(j<640 for j,q in sp(r['constant_original_source_combination']).items())for r in X0['all768_fixed_elimination_receipts_and_constant_combinations']if r['source_column']<640)
 # Recheck the M1 full56 first span witnesses and the actual Pi3 annihilator.
 xp0=[sp(p)for p in M1['actual_pivot_X0_columns27_sparse56']]
 xlin=[{8*r+k-1:F(q)for r,p in col for k,q in p if k and F(q)}for col in X1]
 pivots=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations']if'pivot'in r]
 for i,r in enumerate(pivots):
  a0={};a1={}
  for j,c in sp(r['constant_original_source_combination']).items():
   a0=add(a0,{k:F(row[j])for k,row in enumerate(X0['full_X0_56_by768'])if F(row[j])},c);a1=add(a1,xlin[j],c)
  assert a0==xp0[i]and a1==sp(M1['actual_pivot_X1_by27_sparse56_times8'][i])
 # Reconstruct the620 null combinations from the frozen rank37 full-column
 # factorization, then transport both earlier constant source changes.
 cert=lowfull['complete_linear_rank_certificate'];fac=[sp(p)for p in cert['all_column_factorization_in_normalized_basis']];selectedL=cert['selected_original_columns'];nb=[]
 for k,j in enumerate(selectedL):
  p={j:F(1)}
  for i,q in fac[j].items():
   if i<k:p=add(p,nb[i],-q)
   else:assert i==k
  nb.append({i:q/fac[j][k]for i,q in p.items()})
 raw741=[sp(p)for p in schur['source741_zero_constant_combinations']];raw657=[]
 for c in schur['remaining657_source741_combinations']:
  p={}
  for j,q in sp(c).items():p=add(p,raw741[j],q)
  raw657.append(p)
 built620=[]
 for j in range(657):
  if j in selectedL:continue
  c={j:F(1)}
  for i,q in fac[j].items():c=add(c,nb[i],-q)
  p={}
  for k,q in c.items():p=add(p,raw657[k],q)
  built620.append(p)
 assert built620==[sp(c)for c in M1['620_actual_source_combinations_in_original768']]
 for combo,coeff in zip(M1['620_actual_source_combinations_in_original768'],M1['620_actual_linear_pivot_coefficients_a1_by27_times8']):
  actual={}
  for j,c in sp(combo).items():actual=add(actual,xlin[j],c)
  expected={}
  for ik,c in sp(coeff).items():
   j,k=divmod(ik,8);expected=add(expected,{8*i+k:q for i,q in xp0[j].items()},c)
  assert expected==actual
 Pi=[sp(p)for p in M1['Pi3_180_dual_rows']]
 assert len(Pi)==3 and all(not sum((p.get(i,F(0))*q for i,q in sp(col).items()),F(0))for p in Pi for col in low['all464_certified_quadratic_columns_sparse'])
 assert[[p.get(k,F(0))for k in[143,176,179]]for p in Pi]==[[F(i==j)for j in range(3)]for i in range(3)]
 # Evaluate only the original source DAG, using already solved Artin roots.
 dag=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values()for pair in f['complete_normalized_physical_w_coefficients'][0]for i in pair]
 order=PF['dual9_source_evaluation']['closure_order'];uses=Counter()
 for i in order:
  if dag[i][0]in['a','m','n','i']:uses.update(dag[i][1:])
 uses.update(roots);vals={}
 def binv(p):
  c=p.get(0,F(0));assert c;q=add(p,{0:-c});value={0:1/c};value=add(value,q,-1/c**2);value=add(value,bm(q,q),1/c**3);assert bm(p,value)=={0:F(1)};return value
 for i in order:
  n=dag[i];tag=n[0]
  if tag=='v':v=first[i]
  elif tag=='c':v={0:F(n[1],n[2])}if n[1]else{}
  elif tag=='a':v=add(vals[n[1]],vals[n[2]])
  elif tag=='m':v=bm(vals[n[1]],vals[n[2]])
  elif tag=='n':v={k:-q for k,q in vals[n[1]].items()}
  elif tag=='i':v=binv(vals[n[1]])
  else:raise AssertionError(n)
  vals[i]=v
  if tag in['a','m','n','i']:
   for k in n[1:]:
    uses[k]-=1
    if not uses[k]:del vals[k]
 pairs={f:[[vals[a],vals[b]]for a,b in spec['complete_normalized_physical_w_coefficients'][0]]for f,spec in S['families'].items()}
 assert pairs=={f:[[sp([[i,q]for i,q in enumerate(a)if F(q)]),sp([[i,q]for i,q in enumerate(b)if F(q)])]for a,b in ps]for f,ps in p1['physical150_tracepairs_B11'].items()}
 ring.install();positions=I['actual_group15_positions'];GL=PF['GL18_source_parameter_positions'];rows=P['zero_coordinate_rows36'];pick=P['CA_constant_left_inverse_rows19'];fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices'];center=list(map(F,P['actual_CA_source19_center']))
 def frame(family):
  spec=S['families'][family];p=ring.add(ring.tau(0 if family=='CA'else 1),ring.scale(ring.base(first[spec['p_quadratic_relation']['s']]),F(1,2)));G=V.eye(5)
  for name,(aa,bb)in spec['coefficient_unknown_pairs'].items():
   if name.startswith('gamma')and name!='gamma0':
    i,j=divmod(GL[int(name[5:])-1],5);G[i][j]=ring.add(G[i][j],ring.base(first[aa]),ring.mul(ring.base(first[bb]),p))
  C=V.eye(5);C[0][3]=ring.scale(p,-2 if family=='CA'else F(1,2));C[1][0]=p;C[1][3]=ring.scale(ring.mul(p,p),-1 if family=='CA'else F(1,4));return V.matmul(G,C)
 def physical(family):return[ring.add(ring.base(a),ring.mul(ring.base(b),ring.tau(0 if family=='CA'else 1)))for a,b in pairs[family]]
 def gauge(ev,od):
  G=V.eye(5)
  for(i,j),a,b in zip(positions,ev,od):G[i][j]=ring.add(G[i][j],ring.base(a),ring.mul(ring.base(b),ring.tau(0)))
  return G
 HA=frame('CA');PA=physical('CA')
 def constraints(G):
  Fm=V.normalized_flag_frame(V.matmul(G,HA));pushed=V.tensor_push(PA,V.matmul(Fm,G));own=[r[0]for r in V.matmul(V.constmat(P['CA_left_inverse19']),[[pushed[rows[j]]]for j in pick])]
  return[Fm[1][0],Fm[1][3],Fm[1][4],Fm[0][3],Fm[0][4]]+[ring.add(own[j],ring.c(-center[j]))for j in fixed]
 ev0=[embed_first(p)for p in Old3['even15_first9']];od0=[embed_first(p)for p in Old3['odd15_first9']];res=constraints(gauge(ev0,od0))
 assert all(not c for p in res for(b,mask),c in p.items()if b not in QUAD)
 rhs=[[p.get((b,mask),F(0))for b in QUAD]for mask in range(2)for p in res]
 assert rhs==[[F(c)for c in row]for row in p2['RHS30_by5']]
 corr=[[F(c)for c in row]for row in p2['correction30_by5']]
 A0=[[F(c)for c in row]for row in Old2['A0']];A1=[[F(c)for c in row]for row in Old2['A1_eta_delta']];J=[[F(0)]*30 for _ in range(30)]
 for i,j in product(range(15),repeat=2):J[i][j]=A0[i][j];J[15+i][j]=A1[i][j];J[15+i][15+j]=A0[i][j]
 assert mmQ(J,corr)==[[-c for c in row]for row in rhs]
 ev=[sp([[j,c]for j,c in enumerate(p)if F(c)])for p in p2['even15_B11']];od=[sp([[j,c]for j,c in enumerate(p)if F(c)])for p in p2['odd15_B11']]
 for i in range(15):assert ev[i]==add(ev0[i],dict(zip(QUAD,corr[i])))and od[i]==add(od0[i],dict(zip(QUAD,corr[15+i])))
 G=gauge(ev,od);assert all(not p for p in constraints(G));assert V.normalized_flag_frame(V.matmul(G,HA))==V.eye(5)
 FB=V.normalized_flag_frame(V.matmul(G,frame('CB')));pushed=V.tensor_push(PA,G);own=[p[0]for p in V.matmul(V.constmat(P['CA_left_inverse19']),[[pushed[rows[j]]]for j in pick])]
 actual=[ring.add(own[j],ring.c(-center[j]))for j in free]+[FB[1][0],FB[1][3],FB[1][4],FB[0][3],FB[0][4]];kappa=[actual[j]for j in I['retained_source9_indices']]
 PB=physical('CB');F12=[ring.add(PA[j],ring.neg(PB[j]))for j in[18,19,28,29,60,67,72,63,64,68,73,74]]
 assert actual==[ring.flat(p)for p in p3['actual14_C_B44']]and kappa==[ring.flat(p)for p in p3['kappa9_C_B44']]and F12==[ring.flat(p)for p in p3['physicalF12_C_B44']]
 def match_first(new,old):return all(new.get((bi,mask),F(0))==F(old[4*di+mask])for bi,di in LINEAR for mask in range(4))
 assert all(match_first(p,q)for p,q in zip(actual,Old4['actual14']))and all(match_first(p,q)for p,q in zip(kappa,Old4['kappa9']))
 assert all(match_first(p,q)for p,q in zip(F12,Old4['physicalF12']))
 core=[kappa[j]for j in[0,1,2,6]];D=ring.add(ring.ONE,ring.neg(ring.mul(kappa[1],kappa[1])));PP=ring.add(ring.c(5),ring.mul(kappa[1],kappa[1]));axes=[K['source9_variable_order'].index(n)for n in K['selected5_pivot_order']];rr=[]
 for name in K['selected5_pivot_order']:
  q=K['five_actual_rational_roots'][name];rr.append(ring.tm(V.polynomial(q['numerator'],core),ring.power(D,-q['D_power']),ring.power(PP,-q['P_power'])))
 gs=[ring.add(kappa[j],ring.neg(v))for j,v in zip(axes,rr)]+[V.polynomial(p,core)for p in Prep['generators3']]
 assert gs==[ring.flat(p)for p in p3['actualG8_C_B44']]
 assert all(match_first(p,q)for p,q in zip(gs,Old5['G8']))
 Sn=[[{}]*32 for _ in range(12)]
 for i,j,g,c in Op['native_S12_by32']['koszul28_in_first8_rows']:Sn[i][j]=ring.scale(gs[g],F(c))
 for i,j,c in Op['native_S12_by32']['free4_zero_source_columns']:Sn[i][j]=ring.c(c)
 # Build scalar L0+L1 in C36 FIRST, then embed its scalar base coefficients.
 # No native-N truncation with second monics is used.
 oldnu=[[F(c)for c in row]for row in PF['nu_A_B_first9']];oldring=DualRing(oldnu);L={}
 for j in range(48):
  jn,b=divmod(j,4)
  for i in range(12):
   value=oldring.mul(oldring.flat(Old5['N'][jn][i]),{(0,b):F(1)})
   for mask in range(4):
    p={bi:value.get((di,mask),F(0))for bi,di in LINEAR if value.get((di,mask),F(0))}
    if p:L[4*i+mask,j]=p
 effective=[];sourceS=[];nzE=nzS=0
 for j in range(128):
  col,b=divmod(j,4);scol={}
  for i in range(12):
   value=ring.mul(Sn[i][col],{(0,b):F(1)})
   for mask in range(4):
    p={k:q for(k,m),q in value.items()if m==mask}
    if p:scol[4*i+mask]=p
  ecol={}
  for(i,k),p in L.items():
   if k in scol:ecol[i]=add(ecol.get(i,{}),bm(p,scol[k]))
  ecol={i:p for i,p in ecol.items()if p}
  old={i:sp(p)for i,p in old_array['source_columns128_rows48_sparse_dual9'][j]}
  assert all(ecol.get(i,{}).get(bi,F(0))==old.get(i,{}).get(di,F(0))for i in range(48)for bi,di in LINEAR)
  encode=lambda c:[[i,[[k,str(q)]for k,q in sorted(p.items())]]for i,p in sorted(c.items())]
  sourceS.append(encode(scol));effective.append(encode(ecol));nzS+=sum(map(len,scol.values()));nzE+=sum(map(len,ecol.values()))
 assert sourceS==final['actual_S_source_columns128_rows48_sparse_B11']and effective==final['source_columns128_rows48_sparse_B11_effective_second']
 assert(nzS,nzE)==(2109,6876)
 proof={'all_passed':True,'all400_H20_to_B11_products_and105_root_maps_checked':True,'all1936_C44_monomial_basis_products_checked':True,
  'all17926_source_nodes_replayed_in_B11':True,'fixed150_gauge_coefficients_and_all30_B_constraints_checked':True,
  'actual14_kappa9_G8_C44_source_values_and_actual_monic_carries_checked':True,
  'actual_S128_columns_all2109coefficients_checked':True,'scalar_L0_L1_built_in_C36_before_B11_embedding':True,
  'effective128_columns_all6876coefficients_checked':True,'all741_and620_origin_combinations_separately_old_native_zero':True,
  'all620_full56_firstspan_and_Pi3_known_quadratic_annihilation_checked':True,
  'M1_cancellation_proof':'For x1=Xp0*a1, z=x-Xp*a1 starts degree2; m1 and unknown a2 vanish under m0/pi. Euler commutator has zero boundary; its first part is d3_linear because H2 K(z7,0)=0. Pi3 kills its transported known d3 quadratics.',
  'L2_cancellation_proof':'rank old640=4 plus native23 equals joint27, so origin images are direct. Thus S0*c_native=0 separately and L2*S0*c_native=0 before projection.',
  'effective_not_individual_actual_X2':True,'valid_only_for_proved_zero_origin_source_combinations':True,
  'no_new_source105_or_physical_jet_N2_R12_T2_or_T_length_depth':True}
 (out/'independent-B11-C44-gauge-and-effective-S2-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[str(Path(__file__).relative_to(ROOT)),'code/ambient_B18_native128_first_review_20261104_v63.py','research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,'actual_B11_S2_source_and_fixed_gauge_verified':True,
  'M1_and_source_specific_L2_cancellation_accepted':True,'effective_native128_B11_coefficients_verified':6876,
  'individual_actual_X2_or_T_length_depth_claim':False,'input_sha256':hashes,
  'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
