#!/usr/bin/env python3
"""Exact two-embedding main-M1 exception at 7*c^2=1."""
import time,itertools
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from s61_model import S61
from s61_minimal_baseline_thickenings_20261101 import tensor,coords
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_main_jordan_removal_20261101 import orbit

def number_field_case(m,ranks,main):
 st=time.time();a=s.Symbol('a');root=s.sqrt(7)/7;K=s.QQ.algebraic_field(root);alpha=K.from_sympy(root);zero=K.zero;one=K.one

 def cv(e):
  n,d=s.fraction(s.cancel(e));pp=lambda q:sum((K.from_sympy(c)*alpha**mon[0]for mon,c in s.Poly(q,a).terms()),zero);return pp(n)/pp(d)
 def dm(M):return DM([[cv(v)for v in row]for row in M.tolist()],M.shape,K)
 def null(D):
  R,p=D.rref();R=R.to_list();f=[i for i in range(D.shape[1])if i not in p];N=[[zero]*len(f)for i in range(D.shape[1])]
  for j,i in enumerate(f):
   N[i][j]=one
   for r,c in enumerate(p):N[c][j]=-R[r][i]
  return DM(N,(D.shape[1],len(f)),K)
 def sex(v):
  cc=v.to_list();return sum(s.QQ.to_sympy(x)*a**(len(cc)-i-1)for i,x in enumerate(cc))
 def en(D):return [[str(sex(v))for v in row]for row in D.to_list()]
 source=ranks['cases']['M1'];c=s.Symbol('c');b=s.symbols('b0:9');BK=s.Matrix(source['beta_kernel']).applyfunc(s.sympify);bv=(BK*s.Matrix(b)).subs({b[7]:(5*c+3)*b[7],b[8]:(5*c-3)*b[8]},simultaneous=True).applyfunc(s.cancel).subs(c,a);beta=bv.subs(dict(zip(b,range(1,10))));oi=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];two,F=specialize(m,dict(zip(oi,beta)),True);Q=tensor('M1',0,0,(a-1)/2);q0=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q0)if v);co=[i for i in range(45)if i!=fi];Jc=s.zeros(135,45)
 for r,f in enumerate(F):
  val=0
  for mon,z in f.items():
   if len(mon)==2:
    i,j=mon;val+=z*q0[i]*q0[j]
    if i!=fi:Jc[r,co.index(i)]+=z*q0[j]
    if j!=fi:Jc[r,co.index(j)]+=z*q0[i]
   elif len(mon)==1:Jc[r,44]+=z*q0[mon[0]]
  assert cv(val)==zero
 D=dm(Jc);T=null(D);print('rank',D.rank(),'tan',T.shape[1],time.time()-st,flush=True)
 pars=list(map(s.sympify,main['source_variables']));sb=s.Matrix(main['source_beta']).applyfunc(s.sympify);mq=s.Matrix(main['cases']['M1']['q_coordinates']).applyfunc(s.sympify);sp=s.solve(list(sb-beta),pars,dict=True)[0];assert all(cv(x)==zero for x in mq.subs(sp)-q0);u,act=beta_actions(m);beta_coordinates=u;BD=dm(sb.jacobian(pars).subs(sp).row_join(act.subs(dict(zip(u,beta)))));dq=mq.jacobian(pars).subs(sp);dq-=q0*dq[fi,:]/q0[fi];O=orbit(m,Q,two,fi);QD=dm(dq[co,:].row_join(O[co,:]));C=QD*null(BD);assert C.rank()==1;j=next(j for j in range(C.shape[1])if any(row[j]for row in C.to_list()));tv=C.extract(range(44),[j]).vstack(DM.zeros((1,1),K));rr=list(T.transpose().rref()[1]);vt=T.extract(rr,range(T.shape[1])).inv()*tv.extract(rr,[0]);assert T*vt==tv;first=next(i for i,row in enumerate(vt.to_list())if row[0]);perm=[first]+[i for i in range(T.shape[1])if i!=first];T=T.extract(range(45),perm);vt=vt.extract(perm,[0]);print('main tangent',en(vt),'beta',BD.rank(),time.time()-st,flush=True)
 n=T.shape[1];rows=list(D.transpose().rref()[1]);piv=list(D.rref()[1]);UU=[[zero]*len(piv)for i in range(45)];iv=D.extract(rows,piv).inv().to_list()
 for j,i in enumerate(piv):UU[i]=iv[j]
 UD=DM(UU,(45,len(piv)),K);JU=D*UD;Tlist=T.to_list();mon2=list(itertools.combinations_with_replacement(range(n),2));mon3=list(itertools.combinations_with_replacement(range(n),3));ix={v:i for i,v in enumerate(mon3)};raw=[[zero]*len(mon2)for _ in range(135)];terms=[]
 for row,f in enumerate(F):
  ts=[]
  for mon,z in f.items():
   if fi in mon:continue
   i,j=[co.index(k)for k in mon]+[44]*(2-len(mon));z=cv(z);ts.append((i,j,z))
   for h,(u,v)in enumerate(mon2):raw[row][h]+=z*(Tlist[i][u]*Tlist[j][v]+(Tlist[i][v]*Tlist[j][u]if u!=v else zero))
  terms.append(ts)
 Raw=DM(raw,(135,len(mon2)),K);w2=-(UD*Raw.extract(rows,range(len(mon2))));A2=Raw-JU*Raw.extract(rows,range(len(mon2)));print('Qrank',A2.rank(),time.time()-st,flush=True);wl=w2.to_list();raw3=[[zero]*len(mon3)for _ in range(135)]
 for row,ts in enumerate(terms):
  for i,j,z in ts:
   for k in range(n):
    for h,mon in enumerate(mon2):raw3[row][ix[tuple(sorted((k,)+mon))]]+=z*(Tlist[i][k]*wl[j][h]+Tlist[j][k]*wl[i][h])
 R3=DM(raw3,(135,len(mon3)),K);A3=R3-JU*R3.extract(rows,range(len(mon3)));ar=list(A2.transpose().rref()[1]);ac=list(A2.rref()[1]);C=A2.extract(range(135),ac)*A2.extract(ar,ac).inv();B3=A3-C*A3.extract(ar,range(len(mon3)));br=list(B3.transpose().rref()[1]);v=s.symbols('v0:'+str(n));forms=[]
 for MM,mons,rs in [(A2,mon2,ar),(B3,mon3,br)]:
  for row in MM.extract(rs,range(len(mons))).to_list():forms.append(sum(K.to_sympy(z)*s.prod(v[i]for i in mon)for z,mon in zip(row,mons)))
 def length(g):
  leads=[p.LM(order=g.order).exponents for p in g.polys];counts=[]
  for degree in range(15):
   ct=0
   for mon in itertools.combinations_with_replacement(range(n-1),degree):
    e=tuple(mon.count(i)for i in range(n-1));ct+=not any(all(x<=y for x,y in zip(l,e))for l in leads)
   counts.append(ct)
   if not ct:break
  return sum(counts),counts
 sec=s.groebner([f.subs(v[0],0)for f in forms],*v[1:],order='grevlex',domain=K);gb=s.groebner([f.subs(v[0],1)for f in forms],*v[1:],order='grevlex',domain=K);assert sec.is_zero_dimensional and gb.is_zero_dimensional;print('lengths',length(sec),length(gb),'elapsed',time.time()-st,flush=True);out=dict(field='QQ[a]/(7*a^2-1)',rank=D.rank(),tangent=en(T),main_tangent=en(vt),quadratic_rank=A2.rank(),pure_cubic_rank=len(br),quadratic=en(A2.extract(ar,range(len(mon2)))),pure_cubic=en(B3.extract(br,range(len(mon3)))),section_length=length(sec),localized_length=length(gb),source_beta_minor_rank=BD.rank(),covers_both_embeddings=True);fixedBD=dm(bv.jacobian(b).row_join(act.subs(dict(zip(beta_coordinates,beta)))));assert fixedBD.rank()==16;assert length(gb)[0]==2
 out.update(fixed_c_beta_rank=16,fixed_c_beta_differential=en(fixedBD),main_beta_differential=en(BD),main_q_differential=en(QD),jacobian=en(D),source_point={str(k):str(v)for k,v in sp.items()},section_groebner=[str(g).replace('sqrt(7)','7*a')for g in sec],localized_groebner=[str(g).replace('sqrt(7)','7*a')for g in gb],seconds=time.time()-st)
 return out
