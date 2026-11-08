#!/usr/bin/env python3
"""Independent entire constant-root nonreduced-m3 review over Fraction(tau).
All producer files are data. Only previous independent rational/sparse helpers.
"""
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import ast,json,hashlib,time,resource,signal,argparse
import classification_plane_m4_special_review as H
from classification_unimodular_logarithmic_review import RF,det,rank,roots_supported
ROOT=Path(__file__).resolve().parents[1]
BASE='reports/CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_'
Y=(0,1,3,4);T=RF((F(0),F(1)))

def load(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def matrix(M,R):return[[R.scalar(c)for c in r]for r in M]
def mm(A,B):return[[sum((a*b for a,b in zip(r,c)if a and b),RF())for c in zip(*B)]for r in A]

class Ring:
 def __init__(self,params=(),tau=T):self.names=('a','b','c','z','w')+tuple(params);self.zero=(0,)*len(self.names);self.tau=RF(tau)
 def atom(self,c=1,i=None):
  e=list(self.zero)
  if i is not None:e[i]=1
  return{tuple(e):RF(c)}if c else{}
 def parse(self,s):
  def walk(v):
   if isinstance(v,ast.Constant)and type(v.value)is int:return self.atom(v.value)
   if isinstance(v,ast.Name):return self.atom(self.tau)if v.id=='tau'else self.atom(1,self.names.index(v.id))
   if isinstance(v,ast.UnaryOp):return H.scale(walk(v.operand),-1 if isinstance(v.op,ast.USub) else 1)
   if isinstance(v,ast.BinOp):
    a=walk(v.left)
    if isinstance(v.op,ast.Pow):
     k=int(v.right.value)if isinstance(v.right,ast.Constant)else-int(v.right.operand.value);o=self.atom()
     if k<0:
      assert set(a)<=set([self.zero]);v=RF(1)
      for _ in range(-k):v=v/a[self.zero]
      return self.atom(v)
     for _ in range(k):o=H.mul(o,a)
     return o
    b=walk(v.right)
    if isinstance(v.op,ast.Add):return H.add(a,b)
    if isinstance(v.op,ast.Sub):return H.add(a,H.scale(b,-1))
    if isinstance(v.op,ast.Mult):return H.mul(a,b)
    if isinstance(v.op,ast.Div):assert set(b)<=set([self.zero]);return H.scale(a,RF(1)/b[self.zero])
   raise ValueError(ast.dump(v))
  return walk(ast.parse(str(s),mode='eval').body)
 def scalar(self,s):
  p=self.parse(s);assert set(p)<=set([self.zero]);return p.get(self.zero,RF())
 def skew(self,data):
  P={(i,j):{}for i in range(5)for j in range(5)}
  for(i,j),p in data.items():P[i,j]=p;P[j,i]=H.scale(p,-1)
  return P
 def plus(self,*Ps):return self.skew({ij:H.add(*(P.get(ij,{})for P in Ps))for ij in H.PAIRS})
 def tensor(self,v):
  out={}
  for(ij,e),c in zip(H.COLS,v):
   if c:out.setdefault(ij,{})[e+(0,)*(len(self.names)-5)]=RF(c)
  return self.skew(out)
 def vector(self,P):return[P.get(ij,{}).get(e+(0,)*(len(self.names)-5),RF())for ij,e in H.COLS]
 def strict(self,P):return{tri:H.add(*(H.mul(P[a,l],H.der(P[b,c],l))for a,b,c in(tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1]))for l in range(5)))for tri in combinations(range(5),3)}
 def jac(self,P):
  J=self.strict(P);return{quad:H.add(*(H.scale(H.mul(self.atom(1,quad[k]),J[quad[:k]+quad[k+1:]]),(-1)**k)for k in range(4)))for quad in combinations(range(5),4)}
 def jvec(self,P):
  J=self.jac(P);return[J[q].get(e+(0,)*(len(self.names)-5),RF())for q,e in H.ROWS]
 def actual(self,P):
  D={}
  for at in range(5):
   i,j,k,l=[v for v in range(5)if v!=at];pf=H.add(H.mul(P[i,j],P[k,l]),H.scale(H.mul(P[i,k],P[j,l]),-1),H.mul(P[i,l],P[j,k]));D=H.add(D,H.scale(H.mul(self.atom(1,at),pf),(-1)**at))
  return D
 def basis(self):
  b1,b2=H.sources();convert=lambda P:self.skew({ij:{e+(0,)*(len(self.names)-5):RF(v)for e,v in P.get(ij,{}).items()}for ij in H.PAIRS});b1=b1[:18]+b1[24:27]+b1[18:21]+b1[27:30]+b1[21:24];return list(map(convert,b1)),list(map(convert,b2))
 def decode(self,d,local=False):return self.skew({tuple(Y[i]for i in ast.literal_eval(k))if local else tuple(ast.literal_eval(k)):self.parse(v)for k,v in d.items()})


def cert(M,C,R):
 rr,cc=C['minor_rows'],C['minor_columns'];n=len(M[0]);rk=C['rank'];assert len(rr)==len(cc)==rk;d=det([[M[i][j]for j in cc]for i in rr]);assert d and d==R.scalar(C['minor_determinant']);K=matrix(C['complete_kernel'],R);null=C['nullity'];assert len(K)==n and all(len(row)==null for row in K)and null==n-rk
 if null:assert not any(any(r)for r in mm(M,K))and rank(H.transpose(K))==null
 return K,d

def polynomial_divides(R,f,g):
 assert f and any(any(e[:5])for e in f)
 r=g.copy();lf=max(f);lc=f[lf]
 while r:
  m=max(r)
  if not all(a>=b for a,b in zip(m,lf)):return False
  c=r[m]/lc;delta=tuple(a-b for a,b in zip(m,lf));multiple={tuple(a+b for a,b in zip(e,delta)):v*c for e,v in f.items()};r=H.add(r,H.scale(multiple,-1))
 return True

def point_poly(P,t):return{e:evalrf(c,t)for e,c in P.items()if evalrf(c,t)}
def evalrf(x,t):
 x=RF(x);p=sum((c*t**i for i,c in enumerate(x.p)),F(0));q=sum((c*t**i for i,c in enumerate(x.q)),F(0));assert q;return RF(p/q)

def Wbasis(R):
 out=[]
 for i in(0,1):
  for s in('a*z','b*z','a*w','b*w','z**2','z*w','w**2'):
   v=[{}for _ in range(4)];v[i]=R.parse(s);out.append(v)
 for i,s in[(2,'z*w'),(2,'w**2'),(3,'z**2'),(3,'z*w')]:
  v=[{}for _ in range(4)];v[i]=R.parse(s);out.append(v)
 return out

def pad(R,Q,W):return R.skew({ij:Q.get(ij,{})for ij in H.PAIRS if 2 not in ij}|{(min(2,j),max(2,j)):H.scale(p,1 if 2<j else-1)for j,p in zip(Y,W)})

def module(row,Qspec=None,Wspec=None,tau=T):
 R=Ring(tau=tau);Q=R.decode(Qspec or row['full_Q'],True);W=list(map(R.parse,Wspec or row['W0']));assert not any(R.strict(Q).values());P=pad(R,Q,W);C=[2*v for v in R.jvec(P)];wb=Wbasis(R);cols=[]
 for V in wb:
  P1=pad(R,Q,[H.add(a,b)for a,b in zip(W,V)]);cols.append([2*a-b for a,b in zip(R.jvec(P1),C)])
 M=H.transpose(cols);assert M==matrix(row['full350_affine_matrix'],R)and C==[R.scalar(s[0]if isinstance(s,list)else s)for s in row['full350_affine_constant']];K,dm=cert(M,row['linear_certificate'],R);KA,da=cert([r+[-c]for r,c in zip(M,C)],row['augmented_certificate'],R)
 if row['status']=='inconsistent':
  assert row['augmented_certificate']['rank']>row['linear_certificate']['rank'];return dict(type=row.get('type','wall'),inconsistent=True,G_empty=True),[dm,da],K+KA
 assert row['augmented_certificate']['rank']==row['linear_certificate']['rank'];part=[R.scalar(s)for s in row['particular_solution']];assert all(not sum((a*b for a,b in zip(r,part)),c)for r,c in zip(M,C))
 D=R.actual(P);DL=[H.add(R.actual(pad(R,Q,[H.add(a,b)for a,b in zip(W,V)])),H.scale(D,-1))for V in wb];DP=H.add(D,*(H.scale(d,p)for d,p in zip(DL,part)));DK=[H.add(*(H.scale(d,K[i][j])for i,d in enumerate(DL)))for j in range(row['linear_certificate']['nullity'])];assert DP==R.parse(row['full_actual_D_particular'])and DK==list(map(R.parse,row['full_actual_D_kernel_coefficients']))
 if row['whole_actual_D_zero']:assert not DP and not any(DK)
 else:
  fac=R.parse(row['proper_D_factor']);assert all(polynomial_divides(R,fac,p)for p in[DP]+DK)
 return dict(type=row.get('type','wall'),inconsistent=False,kernel_dimension=len(K[0]),all350_affine_and_full_kernel_actualD=True,proper_factor=row['proper_D_factor'],G_empty=True),[dm,da],K+KA+[part]

def cover_roots(row,dets,K,R,key):
 walls=row[key];roots={F(v['alpha'])for v in walls};polys=[]
 for x in dets:polys.extend([RF(x).p,RF(x).q])
 for r in K:
  for v in r:polys.append(RF(v).q)
 if row.get('all_exception_polynomial'):
  ex=R.scalar(row['all_exception_polynomial']);polys.extend([ex.p,ex.q])
 seen=set()
 for p in polys:seen|=roots_supported(p,roots)
 assert seen==roots
 return walls

def modules_main(raw):
 records=[]
 for row in raw['cases']:
  res,dets,K=module(row);records.append(res)
  if row['whole_alpha_line']:
   for wall in cover_roots(row,dets,K,Ring(),'all_exceptional_fields'):
    if wall.get('source_excluded_here'):
     assert row['type'].startswith('LOG_B_')and wall['alpha']in['-2/3','1/3']
     assert wall['reason']==('coveredbyallfourSCALAR_Bfixedmodules'if wall['alpha']=='-2/3'else'leadingzero')
     assert sum(r['type'].startswith('SCALAR_B_')for r in raw['cases'])==4
     records.append(dict(type=row['type'],alpha=wall['alpha'],source_excluded_to_named_full_complement=wall['reason']));continue
    rr,_,_=module(wall,row['full_Q'],row['W0'],F(wall['alpha']));rr['alpha']=wall['alpha'];records.append(rr)
 assert len(raw['cases'])==22
 return dict(all22_and_all_complex_walls_passed=True,cases=records)

# ROOTS and full Ore-source stages are added before the finite review is run.

def monbideg(b,n):return[e[:3]+(k,n-k)for e in H.monomials(b)if e[3:]==(0,0)for k in range(n+1)]

def initial_model(R,name):
 p=R.parse
 if name.startswith('DOUBLE_LINE'):
  e=int(name.split('_')[2][1]);delta=int(name.endswith('EXACT_C'));gamma=int(name.endswith('EXACT_B'))
  return R.skew({(0,2):p('-a**2'),(1,2):p('2*a*b'),(2,3):p(str(4-3*e)+'*a*z'),(0,4):p(str(-delta)+'*a*w'),(1,4):p(str(2*delta)+'*b*w'),(2,4):p(str(-delta)+'*c*w+3*tau*a*w-'+str(3*gamma)+'*b*w'),(3,4):p(str(delta*(3*e-4))+'*z*w')})
 V={'TRIPLE_LINE_ZERO':['0','0','0'],'TRIPLE_LINE_LINEAR':['0','0','3*a'],'TRIPLE_LINE_SQUARE':['0','0','6*b'],'TRIPLE_LINE_SQUARE_PLUS_LINEAR':['0','-3*a','6*b'],'TRIPLE_LINE_RANK_TWO':['0','-3*b','3*c']}[name]
 return R.skew({(1,2):p('3*a**2')}|{(i,4):H.scale(H.mul(p('w'),p(v)),-1)for i,v in enumerate(V)})

def operator(R,P,basis):
 J=R.jvec(P);return H.transpose([[a-b-c for a,b,c in zip(R.jvec(R.plus(P,Q)),J,R.jvec(Q))]for Q in basis])

def shears(R,P,B):
 labels=[next((ij,e,c)for ij in H.PAIRS for e,c in Q[ij].items()if c)for Q in B];out=[];E=[R.atom(1,i)for i in range(5)]
 for normal in(3,4):
  for idx in range(3):
   O=R.skew({(i,j):H.add(H.scale(H.mul(E[normal],H.der(P[i,j],idx)),-1),P[normal,j]if i==idx else{},P[i,normal]if j==idx else{})for i,j in H.PAIRS});Z=[{}for _ in range(5)]
   for i in range(3):
    ez=[0]*5;ez[3]=2;ew=[0]*5;ew[4]=2;Z[i]=H.add(H.scale(E[3],O[i,3].get(tuple(ez),RF())),H.scale(E[4],O[i,4].get(tuple(ew),RF())))
   O=R.plus(O,R.skew({(i,j):H.add(H.mul(E[i],Z[j]),H.scale(H.mul(E[j],Z[i]),-1))for i,j in H.PAIRS}));v=[O[ij].get(e,RF())/c for ij,e,c in labels];assert R.vector(O)==[sum((a*b for a,b in zip(r,v)),RF())for r in H.transpose([R.vector(Q)for Q in B])];out.append(v)
 return H.transpose(out)

def roots_case(row,P,R,base=None):
 B1,B2=R.basis();base=row if base is None else base
 if 'full150_P0'in row:assert R.vector(P)==R.vector(R.decode(row['full150_P0']))
 assert not R.actual(P)and not any(R.jac(P).values())
 if 'full30_K_consistency_matrix'in row:
  raw=R.skew({ij:P[ij]for ij in H.PAIRS if ij!=(3,4)});D0=R.actual(raw);km=monbideg(3,2);nm=[(0,0,0,2,0),(0,0,0,1,1),(0,0,0,0,2)];C=H.transpose([[H.add(R.actual(R.plus(raw,R.skew({(3,4):{ex:RF(1)}}))),H.scale(D0,-1)).get(e,RF())for e in km]for ex in nm]);constant=[D0.get(e,RF())for e in km]
  assert C==matrix(row['full30_K_consistency_matrix'],R);assert[[v]for v in constant]==matrix(row['full30_K_consistency_constant'],R);cert(C,row['K_certificate'],R);cert([r+[-c]for r,c in zip(C,constant)],row['K_augmented_certificate'],R);kk=[P[3,4].get(e,RF())for e in nm];assert all(not sum((a*b for a,b in zip(r,kk)),c)for r,c in zip(C,constant))and P[3,4]==R.parse(row['unique_global_K'])
 J1,J2=operator(R,P,B1),operator(R,P,B2);mons=monbideg(2,3);D3=H.transpose([[H.add(R.actual(R.plus(P,Q)),H.scale(R.actual(Q),-1)).get(e,RF())for e in mons]for Q in B1]);assert J2==matrix(row['full350_weight2_matrix'],R)
 if 'full350_weight1_operator'in row:assert J1==matrix(row['full350_weight1_operator'],R)and D3==matrix(row['full24by30_actualD3_map'],R)
 target=R.parse(base['exact_nonzero_leading_target']);tv=[target.get(e,RF())for e in mons];Ann=matrix(base['full_line_annihilator'],R);assert len(Ann)==23 and rank(Ann)==23 and not any(any(v)for v in mm(Ann,[[v]for v in tv]));M=J1+mm(Ann,D3);assert M==matrix(row['full_weight1_plus_line_matrix'],R)
 K1,d1=cert(M,row['weight1_certificate'],R);K2,d2=cert(J2,row['weight2_certificate'],R);S=shears(R,P,B1);assert S==matrix(row['full_shear_matrix'],R);KS,ds=cert(S,row['actual_shear_certificate'],R);assert not any(any(r)for r in mm(M,S));image=mm(D3,K1);assert image==matrix(row['actual_D3_complete_kernel_image'],R)
 empty=not any(any(r)for r in image);assert empty==row['whole_m3_source_empty'];assert int(not empty)==row['actual_D3_line_image_rank']
 if not empty:
  at=next(i for i,v in enumerate(tv)if v);j=row['actual_D3_nonzero_certificate']['kernel_column'];assert image[at][j]/tv[at]==R.scalar(row['actual_D3_nonzero_certificate']['coefficient'])
 name=base['type'];B=(10 if 'EXACT_C'in name else 9 if 'EXACT_B'in name else 8)if name.startswith('DOUBLE')else{'TRIPLE_LINE_ZERO':5,'TRIPLE_LINE_LINEAR':7,'TRIPLE_LINE_SQUARE':8,'TRIPLE_LINE_SQUARE_PLUS_LINEAR':9,'TRIPLE_LINE_RANK_TWO':10}[name]
 if base.get('whole_alpha_line')and not row.get('whole_alpha_line',False):B-=1
 k1=row['weight1_certificate']['nullity'];k2=row['weight2_certificate']['nullity'];assert B==row['whole_fixed_plane_symbol_dimension_upper'];upper=B+k1+k2+5;assert upper==row['whole_projective_image_upper']and row['G_dimension_criterion']==(empty or upper<=24)
 return dict(type=name,alpha=row.get('alpha','whole'),globalK_and_P0_exact=True,whole350_30_9_and_D3line=True,K1=k1,K2=k2,whole_first_symbol_upper=B,whole_image_upper=upper,entire_m3_empty=empty,G_closed=empty or upper<=24),[d1,d2,ds],K1+K2

def roots_main(raw):
 R=Ring();mons=monbideg(2,3);assert raw['all24_leading_monomials']==[list(e)for e in mons]and raw['all30_K_monomials']==[list(e)for e in monbideg(3,2)];B1,B2=R.basis()
 RP=Ring(tuple('l'+str(i)for i in range(9))+tuple('n'+str(i)for i in range(9))+tuple(s+str(i)for s in('d','b','e','c')for i in range(3)))
 for key,expr in raw['full_source1'].items():
  pair=ast.literal_eval(key);expect={}
  for j,Q in enumerate(B1):
   for ex,c in Q[pair].items():ee=list(ex)+( [0]*30);ee[5+j]=1;expect=H.add(expect,{tuple(ee):c})
  assert RP.parse(expr)==expect
 records=[]
 for row in raw['cases']:
  P=initial_model(R,row['type']);res,dets,K=roots_case(row,P,R);records.append(res)
  if row['whole_alpha_line']:
   for wall in cover_roots(row,dets,K,R,'all_exceptional_roots'):
    RW=Ring(tau=F(wall['alpha']));p=initial_model(RW,row['type']);rr,_,_=roots_case(wall,p,RW,row);records.append(rr)
 assert len(raw['cases'])==14
 return dict(all14_and_every_alpha_wall_passed=True,cases=records,old_six_coarse_false_not_rewritten=True,Cexact_only_G_not_wholeE=True)

def keep_spatial(P):return{e[:5]:c for e,c in P.items()if not any(e[5:])}
def coefficient_source(P,spatial,n):return[P.get(spatial+tuple(int(i==j)for i in range(n)),RF())for j in range(n)]

def ore_case(row,oldroot,tau=T):
 nr=len(row['original_full_M1_certificate']['complete_kernel'][0]);names=tuple('p'+str(i)for i in range(nr))+tuple('r'+str(i)for i in range(9));R=Ring(names,tau);K=matrix(row['original_full_M1_certificate']['complete_kernel'],R);oldK=matrix(oldroot['weight1_certificate']['complete_kernel'],R);assert K==oldK
 B1,B2=R.basis();P0=initial_model(R,row.get('type',oldroot.get('type')));data={ij:P0[ij].copy()for ij in H.PAIRS}
 for ij in H.PAIRS:
  for j in range(nr):
   V=H.add(*(H.scale(Q[ij],K[i][j])for i,Q in enumerate(B1)));data[ij]=H.add(data[ij],H.mul(R.atom(1,5+j),V))
  for j,Q in enumerate(B2):data[ij]=H.add(data[ij],H.mul(R.atom(1,5+nr+j),Q[ij]))
 P=R.skew(data);E=[R.atom(1,i)for i in range(5)];L=[H.scale(H.der(P[0,2],2),-1),H.scale(H.der(P[1,2],2),-1),{},{},{}];assert L==list(map(R.parse,row['literal_first_Euler_gauge']));G=R.skew({(i,j):H.add(P[i,j],H.scale(H.add(H.mul(E[i],L[j]),H.scale(H.mul(E[j],L[i]),-1)),-1))for i,j in H.PAIRS});Lc=R.parse(row['second_Euler_c_normal_component']);L=[{},{},Lc,{},{}];G=R.skew({(i,j):H.add(G[i,j],H.scale(H.add(H.mul(E[i],L[j]),H.scale(H.mul(E[j],L[i]),-1)),-1))for i,j in H.PAIRS});assert all(not H.der(G[ij],2)for ij in H.PAIRS)
 mons=[e for e in H.monomials(2)if e[2]==0];labels=[(i,j,e)for i,j in combinations(Y,2)for e in mons];QM=[coefficient_source(G[i,j],e,len(names))for i,j,e in labels];assert QM==matrix(row['full60by_source_Q_coefficient_matrix'],R);QK,dq=cert(QM,row['Q_space_certificate'],R);qnum=row['Q_space_certificate']['rank'];assert qnum==row['Q_source_dimension']
 QR=Ring(tuple('Q'+str(i)for i in range(qnum)),tau);Qfam=QR.decode(row['full_affine_Q'],True);qorigin=QR.decode(row['Q0'],True);cols=row['Q_space_certificate']['minor_columns']
 for i,j,e in labels:
  assert Qfam[i,j].get(e+(0,)*qnum,RF())==keep_spatial(G[i,j]).get(e,RF())==qorigin[i,j].get(e+(0,)*qnum,RF())
  for k,col in enumerate(cols):assert coefficient_source(Qfam[i,j],e,qnum)[k]==QM[labels.index((i,j,e))][col]
 W=[G[2,i]for i in Y];W0=list(map(R.parse,row['full_W0']));assert all(keep_spatial(v)==keep_spatial(w)for v,w in zip(W,W0));allowed=[set()for _ in Y]
 for V in Wbasis(R):
  for i,p in enumerate(V):allowed[i]|={e[:5]for e in p}
 for i in range(4):assert all(e[:5]in allowed[i]for e in H.add(W[i],H.scale(W0[i],-1)))
 J=QR.strict(Qfam)
 for key,v in row['full_strict_Q_Jacobi'].items():assert J[tuple(Y[i]for i in ast.literal_eval(key))]==QR.parse(v)
 name=row.get('type',oldroot.get('type'));constantfactor=None
 if name in('DOUBLE_LINE_E0_EXACT_B','DOUBLE_LINE_E2_EXACT_B'):
  vv=0 if 'E0'in name else 1;assert all(not p for ij,p in Qfam.items()if ij[0]<ij[1]and vv not in ij);assert all(not H.der(p,vv)for p in Qfam.values());constantfactor=vv
 return dict(type=name,alpha=row.get('alpha','whole'),whole_complete_K1_plus_full9P2_Euler_source=True,all_Q_coefficients_and_complete_image=True,Q_dimension=qnum,full_W18_is_covering_supersource=True,all_full_strict_Q_Jacobi_reconstructed=True,constant_factor_direction=constantfactor,negative_rowspan_diagnostic_not_used_in_proof=True),[dq],QK

def ore_main(raw,roots):
 rootmap={r['type']:r for r in roots['cases']};records=[]
 for row in raw['cases']:
  rr=ore_case(row,rootmap[row['type']]);records.append(rr[0])
  if row['whole_alpha_line']:
   oldwalls={F(w['alpha']):w for w in rootmap[row['type']]['all_exceptional_roots']}
   for wall in row['all_exceptional_sources']:
    tt=F(wall['alpha']);assert tt in oldwalls;wall=dict(wall,type=row['type']);result=ore_case(wall,oldwalls[tt],tt);records.append(result[0])
 assert len(raw['cases'])==6
 return dict(all6_and_all_source_walls_Euler_Q_W_covered=True,cases=records,negative_factor_rowspan_false_preserved_not_reproved_or_used=True)


def main():
 parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['all','sources','roots','ore','modules'],default='all');parser.add_argument('--job',required=True);parser.add_argument('--output',required=True);args=parser.parse_args();signal.alarm(420);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);t=time.monotonic();job=load(args.job);assert all(sha(b['path'])==b['sha256']for b in job['exact_inputs_and_hashes']);stages={}
 roots=load(BASE+'ROOTS_20261107.json');ore=load(BASE+'ORE_20261107.json');mods=load(BASE+'MODULES_20261107.json')
 for r in(roots,ore,mods):assert hashlib.sha256(r['producer_bytes'].encode()).hexdigest()==r['producer_sha256']
 if args.mode in('all','sources','roots'):stages['roots']=roots_main(roots);print('PASS independent14 source+allwalls',flush=True)
 if args.mode in('all','sources','ore'):stages['ore']=ore_main(ore,roots);print('PASS independent6 Euler Q/W+allwalls',flush=True)
 if args.mode in('all','modules'):stages['modules']=modules_main(mods)
 out=dict(task_id='G',all_requested_independent_arithmetic_passed=True,mode=args.mode,stages=stages,no_producer_or_CAS_imports=True,all_whole_constant_root_G_closure_requires_all_stages_and_geometry=args.mode=='all',Cexact_only_weak_G_bound_not_E=True,other_root_degrees_and_all_G_not_claimed=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=time.monotonic()-t,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss);p=ROOT/args.output;assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in('all_requested_independent_arithmetic_passed','mode','seconds','maximum_RSS_KiB')}))
if __name__=='__main__':main()
