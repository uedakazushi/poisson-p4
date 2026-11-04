#!/usr/bin/env python3
"""Independent exact source/formal-algebra review; no author imports or huge-core expansion."""
import argparse,hashlib,itertools,json,math,resource,signal,time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63_REVIEW'
BASE='S61_NONSQUARE_CUBIC_LAST_Q9_POWER_BEZOUT_20261104_v63';DEADLINE=0;MAX_TERMS=0
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class P:
 __slots__=('n','d')
 def __init__(self,n,d=None):self.n=n;self.d={tuple(m):F(c)for m,c in(d or{}).items()if c}
 def coerce(self,p):return p if isinstance(p,P)else const(self.n,p)
 def __bool__(self):return bool(self.d)
 def __eq__(self,p):p=self.coerce(p);return self.n==p.n and self.d==p.d
 def __neg__(self):return P(self.n,{m:-c for m,c in self.d.items()})
 def __add__(self,p):
  p=self.coerce(p);assert p.n==self.n;z=self.d.copy()
  for m,c in p.d.items():
   q=z.get(m,F())+c
   if q:z[m]=q
   else:z.pop(m,None)
  return P(self.n,z)
 __radd__=__add__
 def __sub__(self,p):return self+-self.coerce(p)
 def __rsub__(self,p):return self.coerce(p)+-self
 def __mul__(self,p):
  global MAX_TERMS
  p=self.coerce(p);assert self.n==p.n;z={}
  for m,c in self.d.items():
   for k,t in p.d.items():
    e=tuple(a+b for a,b in zip(m,k));q=z.get(e,F())+c*t
    if q:z[e]=q
    else:z.pop(e,None)
  MAX_TERMS=max(MAX_TERMS,len(z));assert len(z)<20000,'unexpected polynomial expansion';assert time.monotonic()<DEADLINE,'50s independent arithmetic boundary';return P(self.n,z)
 __rmul__=__mul__
 def __pow__(self,n):
  assert n>=0;z=const(self.n,1);b=self
  while n:
   if n&1:z=z*b
   n//=2
   if n:b=b*b
  return z
def const(n,c):return P(n,{(0,)*n:F(c)})
def variable(n,i):return P(n,{tuple(int(k==i)for k in range(n)):F(1)})
def dec(ts,n):
 assert all(len(m)==n for m,c in ts);return P(n,{tuple(m):F(c)for m,c in ts})
def frac3(q):
 n,d=dec(q[0],3),dec(q[1],3);assert len(d.d)==1;exp,c=next(iter(d.d.items()));assert exp[1:]==(0,0)and c
 return P(3,{tuple(a-b for a,b in zip(m,exp)):v/c for m,v in n.d.items()})
def inv_monomial(p):
 assert len(p.d)==1;m,c=next(iter(p.d.items()));assert m[1:]==(0,0)and c;return P(3,{tuple(-a for a in m):1/c})
def lift(p):return P(6,{m+(0,0,0):c for m,c in p.d.items()})
def ym(m):return P(6,{(0,0,0)+tuple(m):F(1)})
def coeff(p,m):return P(3,{k[:3]:c for k,c in p.d.items()if k[3:]==tuple(m)})
def encode(p):return[[list(m),str(c)]for m,c in sorted(p.d.items(),reverse=True)]
def core_polynomial(p):assert all(min(m)>=0 for m in p.d);return p
def main(out):
 global DEADLINE
 began=time.monotonic();DEADLINE=began+50;signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 primary=ROOT/'reports'/BASE;pr=read(primary/'report.json');pc=read(primary/'certificate.json');assert sha(primary/'certificate.json')==pr['certificate_sha256'];hashes={}
 for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;hashes[p]=h
 cp=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63/checkpoint.json');data=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63/input.json')
 pair=read(ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63/certificate.json');bc=read(ROOT/'reports/S61_NONSQUARE_CUBIC_LAST_Q9_PREFLIGHT_20261104_v63/certificate.json')
 qb=ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63';qr=read(qb/'report.json');qc=read(qb/'certificate.json');assert sha(qb/'certificate.json')==qr['certificate_sha256'];gate=ROOT/'reports/S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63_REVIEW/report.json';qgate=read(gate);assert qgate['all_passed']and qgate['input_sha256'][str((qb/'certificate.json').relative_to(ROOT))]==sha(qb/'certificate.json')
 assert qc['new15_source_original25_labels'][0]==9 and qr['actual_new_shape_same_cokernel']==[15,22]
 for p in[primary/'report.json',primary/'certificate.json',qb/'report.json',qb/'certificate.json',gate]:hashes[str(p.relative_to(ROOT))]=sha(p)
 a,e,d=[variable(3,i)for i in range(3)];one3=const(3,1);zero3=const(3,0);ys=[variable(6,i)for i in range(3,6)];pairmons=[(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
 datum=lambda i,j:frac3(cp['matrix25x32'][i][j]);Lcols=[21,22,25];b=[datum(5,j)for j in Lcols];sv=[datum(9,j)for j in Lcols];v=[(108*e+42)*inv_monomial(5*a),const(3,-18),const(3,-9)];cc=sum((x*y for x,y in zip(sv,v)),zero3)
 assert sum((x*y for x,y in zip(b,v)),zero3)==one3 and cc==frac3(bc['c_s_dot_v'])
 selected=[11,13,14,15,17,18,19,22,25];assert selected==pc['actual_h9_original_target_labels']==bc['query_only_original32_target_columns']
 roots=[qc['remaining15x22_scalar_roots'][0][qc['new22_target_original32_columns'].index(j)]for j in selected];assert roots==bc['queried_roots'];nodes=qc['scalar_DAG']['nodes'];needed=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i)
  if nodes[i][0]in['a','s','m']:todo.extend(nodes[i][1:])
 actual={};refs=set();allowed={(i,j)for i in[5,9]for j in selected[:7]+Lcols}
 for i in sorted(needed):
  n=nodes[i];tag=n[0]
  if tag=='z':q=zero3
  elif tag=='one':q=one3
  elif tag=='checkpoint':ref=tuple(n[1:]);assert ref in allowed;refs.add(ref);q=datum(*ref)
  elif tag=='new_unit_atom':q=frac3(qc['scalar_DAG']['new_unit_atoms_num_den'][n[1]])
  else:assert tag in['a','s','m'];q=actual[n[1]]+actual[n[2]]if tag=='a'else actual[n[1]]-actual[n[2]]if tag=='s'else actual[n[1]]*actual[n[2]]
  actual[i]=q
 hvals=[core_polynomial(actual[i])for i in roots];direct=[datum(9,j)-cc*datum(5,j)for j in selected[:7]]+[sv[i]-cc*b[i]for i in[0,1]]
 assert hvals==direct==[dec(p,3)for p in pc['actual_h9']]==[frac3(z['actual_coefficient'])for z in bc['actual_nine_coefficients']]
 assert len(needed)==114 and len(refs)==19
 # Original quadratic sources and all four pair-dual identities.
 initial=[[dec(q,3)for q in row]for row in data['matrix25x32'][:10]];pairs=[[initial[i][j]for j in range(20,26)]for i in[8,9]];assert pairs==[[frac3(q)for q in row]for row in pair['Lframe_pair2x6']]
 dual=[[frac3(q)for q in row]for row in pair['Lframe_dual6x2']]
 for i in range(2):
  for j in range(2):assert sum((pairs[i][k]*dual[k][j]for k in range(6)),zero3)==const(3,int(i==j))
 phi=[row[1]for row in dual]
 P0=sum((lift(c)*ym(m)for c,m in zip(pairs[0],pairmons)),const(6,0));RR=sum((lift(c)*ym(m)for c,m in zip(pairs[1],pairmons)),const(6,0));ell=sum((lift(initial[5][20+i])*ys[i]for i in range(3)),const(6,0))
 assert coeff(P0,(2,0,0))==const(3,F(-5,4))and coeff(P0,(1,0,1))==const(3,F(-5,2))and not coeff(P0,(0,0,2));assert -F(5,4)**2==F(-25,16)
 for row,var in zip([5,6,7],ys):assert[initial[row][j]for j in range(20,26)]==[coeff(ell*var,m)for m in pairmons]
 encoded_y=[sum((lift(frac3(c))*ym(m)for m,c in rec),const(6,0))for rec in pc['P_and_R_and_ell_y']];assert encoded_y==[P0,RR,ell]
 # Only the fixed quadratic-source transcript; skip the nonquadratic eighth row after checking no feedback.
 V=[[const(3,int(i==j))for j in range(10)]for i in range(10)]
 for st in cp['all_completed_unit_steps']:
  ri=st['pivot_row']
  if ri>=10:assert not any(i<10 for i,c in st['eliminations']);continue
  inverse=inv_monomial(frac3(st['pivot']));V[ri]=[q*inverse for q in V[ri]]
  for i,c in st['eliminations']:
   if i<10:co=frac3(c);V[i]=[x-co*y for x,y in zip(V[i],V[ri])]
 assert all(not V[4][i]for i in range(10)if i not in[0,4]);assert[data['source25_labels'][i][2]for i in[0,4]]==['q0*x0','scaled_q5*x0']
 vf=[x-cc*y for x,y in zip(V[9],V[5])];assert vf[9]==one3
 m=sum((lift(vf[i])*ys[i-5]for i in[5,6,7]),const(6,0));n=sum((lift(vf[i])*ys[i-1]for i in[1,2,3]),const(6,0));lam=lift(vf[8]);puremons=[(2,1,0),(1,2,0),(1,1,1),(1,0,2),(0,2,1),(0,1,2),(0,0,3)]
 AA=sum((lift(hvals[i])*ym(mon)for i,mon in enumerate(puremons)),const(6,0));BB=lift(hvals[7])*ys[0]*ys[1]+lift(hvals[8])*ys[0]*ys[2]+lift((v[0]*hvals[7])*F(1,9)-2*hvals[8])*ys[2]**2
 assert AA==RR*m+P0*n and BB==RR+ell*m+lam*P0
 # Both polynomial coordinate maps and monic reduction; q is never inverted.
 S=-F(5,2)*ys[0]+lift((9*e+4)*F(1,5))*ys[1];T=ys[1];W=ys[2]-F(1,5)*S-lift((3*e+8)*F(1,25))*T;q=(625*a+108*e**2+156*e+48)*F(1,125);assert q==dec(pc['q'],3);assert P0==S*W+lift(q)*T*T
 s0,t0,w0=ys;inverse=[-F(2,5)*(s0-lift((9*e+4)*F(1,5))*t0),t0,w0+F(1,5)*s0+lift((3*e+8)*F(1,25))*t0]
 def sub(f):
  z=const(6,0)
  for ex,c in f.d.items():
   term=P(6,{ex[:3]+(0,0,0):c})
   for im,k in zip(inverse,ex[3:]):term=term*im**k
   z=z+term
  return z
 assert sub(S)==s0 and sub(T)==t0 and sub(W)==w0
 qpow={0:const(6,1)}
 def nf(f):
  z=const(6,0)
  for ex,c in f.d.items():
   k=min(ex[3],ex[5])
   if k not in qpow:qpow[k]=lift(-q)**k
   mon=ex[:3]+(ex[3]-k,ex[4]+2*k,ex[5]-k);z=z+P(6,{mon:c})*qpow[k]
  return z
 assert sub(P0)==s0*w0+lift(q)*t0*t0;Rbar=nf(sub(RR));mons2=[(2,0,0),(1,1,0),(0,2,0),(0,1,1),(0,0,2)];rv=[core_polynomial(coeff(Rbar,m))for m in mons2];assert rv==[dec(x,3)for x in pc['r5']]
 Tbasis=[S*S,S*T,T*T,T*W,W*W];tv=[sum((coeff(f,m)*phi[i]for i,m in enumerate(pairmons)),zero3)for f in Tbasis];Tpoly=[core_polynomial(a*a*x)for x in tv];assert tv==[frac3(x)for x in pc['transported_pair_functional_t5']]and Tpoly==[dec(x,3)for x in pc['polynomial_T5_a2_times_t']];assert sum((x*y for x,y in zip(Tpoly,rv)),zero3)==a*a
 # Universal coefficient and power identities by independent integer/Fraction arithmetic.
 uv=[variable(6,i)for i in range(6)];u,r0,r1,r2,r3,r4=uv;rs=uv[1:];J=[r0*r0,2*r0*r1,r1*r1+2*r0*r2,2*r1*r2-2*u*r0*r3,r2*r2-2*u*r1*r3+2*u*u*r0*r4,2*r2*r3-2*u*r1*r4,r3*r3+2*r2*r4,2*r3*r4,r4*r4]
 assert J==[dec(x,6)for x in pc['universal_q_r_generator9']];bounds=[2,4,6,4,2];assert pc['universal_power_bounds']==bounds;witnesses=[[dec(x,6)for x in row]for row in pc['universal_power_witness5x9']]
 for i in range(5):assert sum((x*y for x,y in zip(witnesses[i],J)),const(6,0))==rs[i]**bounds[i]
 def evaluni(f):
  z=zero3
  for ex,c in f.d.items():
   term=const(3,c)*q**ex[0]
   for v,k in zip(rv,ex[1:]):term=term*v**k
   z=z+term
  return z
 Jactual=[evaluni(x)for x in J];basis4=[(4,0,0),(3,1,0),(2,2,0),(1,3,0),(0,4,0),(0,3,1),(0,2,2),(0,1,3),(0,0,4)];assert pc['degree4_normal_basis']==[list(x)for x in basis4]
 square=nf(Rbar*Rbar);assert[coeff(square,x)for x in basis4]==Jactual;assert nf(sub(RR*BB-ell*AA))==square
 A9=[[zero3 for j in range(9)]for i in range(9)]
 for j in range(9):
  if j<7:expr=-lift(a)*ell*ym(puremons[j])
  elif j==7:expr=RR*(lift(a)*ys[0]*ys[1]+lift((36*e+14)*F(1,15))*ys[2]**2)
  else:expr=RR*(lift(a)*ys[0]*ys[2]-lift(2*a)*ys[2]**2)
  f=nf(sub(expr))
  for i,mon in enumerate(basis4):A9[i][j]=core_polynomial(coeff(f,mon))
 assert A9==[[dec(x,3)for x in row]for row in pc['a_J_equals_matrix_times_actual_h9']]
 for i in range(9):assert a*Jactual[i]==sum((A9[i][j]*hvals[j]for j in range(9)),zero3)
 # Formal commutative signatures prove the finite coefficient circuit without expanding corecofactors.
 dag=pc['scalar_DAG'];dn=dag['nodes'];ap=[core_polynomial(dec(x,3))for x in dag['atoms_polynomial_QQ_a_e_d']];assert len(dn)==pr['scalar_nodes']==25572 and len(ap)==pr['core_polynomial_atoms']==145;canonical=lambda p:tuple(sorted(p.d.items()));aliases={};atomkey=[]
 for p in ap:
  c=canonical(p)
  if c not in aliases:aliases[c]=len(aliases)
  atomkey.append((0,aliases[c]))
 def psig(p):
  if not p:return(F(0),())
  if set(p.d)=={(0,0,0)}:return(next(iter(p.d.values())),())
  return(F(1),(((0,aliases[canonical(p)]),1),))
 def sigmul(x,y):
  c=x[0]*y[0]
  if not c:return(F(0),())
  ex=dict(x[1])
  for k,v in y[1]:ex[k]=ex.get(k,0)+v
  return c,tuple(sorted((k,v)for k,v in ex.items()if v))
 def monsig(i,cut):
  stack=[i];v=(F(1),())
  while stack:
   j=stack.pop()
   if j in cut:v=sigmul(v,(F(1),(((1,j),1),)));continue
   n=dn[j]
   if n[0]=='z':return(F(0),())
   if n[0]=='one':continue
   if n[0]=='c':v=sigmul(v,psig(ap[n[1]]))
   elif n[0]=='m':stack.extend(n[1:])
   else:raise AssertionError('unexpected uncut additive factor')
  return v
 def signature(root,cut=frozenset()):
  todo=[root];ans={}
  while todo:
   i=todo.pop()
   if i not in cut and dn[i][0]=='a':todo.extend(dn[i][1:]);continue
   c,m=monsig(i,cut)
   if c:
    z=ans.get(m,F())+c
    if z:ans[m]=z
    else:ans.pop(m,None)
  return ans
 def expected_uni(f):
  ans={}
  for ex,c in f.d.items():
   v=(c,())
   for p,k in zip([q]+rv,ex):
    for _ in range(k):v=sigmul(v,psig(p))
   ans[v[1]]=ans.get(v[1],F())+v[0]
  return{k:v for k,v in ans.items()if v}
 degrees=[]
 for i,n in enumerate(dn):
  tag=n[0];assert tag in['z','one','c','a','m']
  if tag in['z','one']:degree=0
  elif tag=='c':assert 0<=n[1]<len(ap);degree=max([sum(m)for m in ap[n[1]].d]+[0])
  else:assert all(0<=j<i for j in n[1:]);degree=max(degrees[j]for j in n[1:])if tag=='a'else sum(degrees[j]for j in n[1:])
  degrees.append(degree)
 wr=pc['universal_power_witness_scalar_roots']
 for i in range(5):
  for j in range(9):assert signature(wr[i][j])==expected_uni(witnesses[i][j])
 terms=pc['multinomial3060_terms_exp_choice_integer_multiplier_root'];bars=set()
 for b0,b1,b2,b3 in itertools.combinations(range(18),4):bars.add((b0,b1-b0-1,b2-b1-1,b3-b2-1,17-b3))
 assert len(bars)==3060 and len(terms)==3060;seen=set();coefsum=0;termroots=[]
 for ex,i,c,root in terms:
  ex=tuple(ex);assert ex in bars and ex not in seen;seen.add(ex);assert i==next(k for k in range(5)if ex[k]>=bounds[k])
  count=1;left=14
  for k in ex:count*=math.comb(left,k);left-=k
  assert c==count;coefsum+=c;v=(F(c),())
  for p,k in zip(Tpoly,ex):
   for _ in range(k):v=sigmul(v,psig(p))
  for j,(p,k)in enumerate(zip(rv,ex)):
   for _ in range(k-(bounds[i]if j==i else 0)):v=sigmul(v,psig(p))
  assert monsig(root,frozenset())==v;termroots.append(root)
 assert seen==bars and coefsum==5**14 and len(set(termroots))==3060 and sum(k-1 for k in bounds)==13
 cut=frozenset(termroots+[r for row in wr for r in row if r>1]);acc=pc['universal_J_cofactor_roots']
 for j in range(9):
  expected={}
  for (_,i,c,root)in terms:
   v=sigmul(monsig(root,cut),monsig(wr[i][j],cut))
   if v[0]:expected[v[1]]=expected.get(v[1],F())+v[0]
  assert signature(acc[j],cut)=={k:v for k,v in expected.items()if v}
 Ccut=frozenset(i for i in acc if i>1);Hroots=pc['actual_h9_cofactor_roots']
 for j in range(9):
  expected={}
  for k in range(9):
   v=sigmul(monsig(acc[k],Ccut),psig(A9[k][j]))
   if v[0]:expected[v[1]]=expected.get(v[1],F())+v[0]
  assert signature(Hroots[j],Ccut)=={k:v for k,v in expected.items()if v}
 Hcut=frozenset(i for i in Hroots if i>1);expected={}
 for root,h in zip(Hroots,hvals):
  v=sigmul(monsig(root,Hcut),psig(h))
  if v[0]:expected[v[1]]=expected.get(v[1],F())+v[0]
 assert signature(pc['identity_right_root_sum_Hi_hi'],Hcut)==expected;lhs=dn[pc['identity_left_root_a29']];assert lhs[0]=='c'and ap[lhs[1]]==a**29
 assert[degrees[r]for r in Hroots]==pc['polynomial_cofactor_degree_upper_bounds']and max(degrees[r]for r in Hroots)==207
 proof={'own_Fraction_Laurent_engine':True,'actual9_coefficients_from_reviewed_GL3_roots_and114_ancestors':True,'old_quadratic_source_transcript_and_x0_vanishing_row4_checked':True,'all4_pairdual_identities_checked':True,
  'all_source_P_R_ell_A_B_identities_and_inverse_coordinate_map_checked':True,'actual_r5_T5_dot_equals_a2_checked':True,'all81_source_A9_coefficients_and9_membership_identities_checked':True,
  'all5_universal_power_witnesses_literal_checked':True,'all3060_terms_starsbars_and_binomial_multinomials_and_factor_roots_checked':True,'all9_J_cofactor_accumulations_and9_actual_H_cofactors_checked_formally':True,
  'all25572_DAG_nodes_typed_acyclic_polynomial_and_degree_bounds_checked':len(dn)==25572,'actual_polynomial_identity_sum_Hi_hi_equals_a29_accepted':True,'actual_dual_denominator_only_a29':True,
  'no_degree207_corecofactor_expansion_authorbuilder_GB_point_or_degree_search':True,'no_GL10_target_completion_or_full57_full239_physical_C0_application_claim':True,'max_expanded_small_polynomial_terms':MAX_TERMS}
 (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[Path(__file__),ROOT/'research/notes'/(ID+'.md')]:hashes[str(p.relative_to(ROOT))]=sha(p)
 result={'id':ID,'all_passed':True,'actual9row_polynomial_Bezout_a29_independently_verified':True,'all5_universal_and3060_structural_witnesses_verified':True,'input_sha256':hashes,'evidence_sha256':{'independent-proof.json':sha(out/'independent-proof.json')},'seconds':time.monotonic()-began,'max_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
 signal.alarm(0);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_passed':True,'seconds':result['seconds'],'max_RSS_KiB':result['max_RSS_KiB'],'report_sha256':sha(out/'report.json')},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'reports'/ID);main(p.parse_args().output.resolve())
