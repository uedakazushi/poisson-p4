#!/usr/bin/env python3
"""Exact universal main-curve full-Jacobian and normal-quadratic certificates."""
import itertools,time
import sympy as s
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_baseline_thickenings_20261101 import enc
def coeffs(F,T,cols):
 mons=list(itertools.combinations_with_replacement(range(T.cols),2));H=s.zeros(135,len(mons));fixed=next(i for i in range(45)if i not in cols);tf={i:T[cols.index(i),:]for i in cols};tf[fixed]=s.zeros(1,T.cols)
 for row,f in enumerate(F):
  for mon,v in f.items():
   aa,bb=[tf[i]for i in mon]+[T[-1,:]]*(2-len(mon))
   for j,(u,w) in enumerate(mons):H[row,j]+=v*(aa[u]*bb[w]+(aa[w]*bb[u]if u!=w else 0))
 return H,mons
def uniform(m,source,kind):
 c=s.Symbol('c');b=s.symbols('b0:9');one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];st=time.time()
 r=source['cases'][kind];sub={b[7]:(5*c-6)*b[7],b[8]:(5*c+6)*b[8]}if kind=='M0' else{b[7]:(5*c+3)*b[7],b[8]:(5*c-3)*b[8]};K=s.Matrix(r['beta_kernel']).applyfunc(s.sympify);bv=(K*s.Matrix(b)).subs(sub,simultaneous=True).applyfunc(s.cancel);two,F=specialize(m,dict(zip(one,bv)),True);q=s.Matrix(r['canonical_q_coordinates']).applyfunc(s.sympify);fi=next(i for i,v in enumerate(q)if v);cols=[i for i in range(45)if i!=fi];A=s.zeros(135,45);tc=s.zeros(135,1)
 for row,f in enumerate(F):
  for mon,v in f.items():
   if len(mon)==2:
    i,j=mon;A[row,i]+=v*q[j];A[row,j]+=v*q[i]
   elif len(mon)==1:tc[row]+=v*q[mon[0]]
 assert all(s.expand(sum(v*s.prod(q[i]for i in mon)for mon,v in f.items()if len(mon)==2))==0 for f in F);A=A[:,cols];pr=r['first_pivot_rows'];pc=r['first_pivot_columns'];U=s.zeros(44,len(pc));iv=A[pr,pc].inv(method='DM')
 for j,i in enumerate(pc):U[i,:]=iv[j,:]
 C=s.Matrix(r['first_kernel']).applyfunc(s.sympify);N=s.Matrix(r['rank_certificates']['135']['residual_kernel']).applyfunc(s.sympify).subs(sub,simultaneous=True).applyfunc(s.cancel);T=(C*N[:-1,:]-U*tc[pr,:]*N[-1,:]).col_join(N[-1,:]).applyfunc(s.cancel);assert(A.row_join(tc).to_DM()*T.to_DM()).is_zero_matrix;E=s.Matrix(r['residual_matrix']).applyfunc(s.sympify).row_join(s.Matrix(r['residual_t_column']).applyfunc(s.sympify)).subs(sub,simultaneous=True).applyfunc(s.cancel);# E has rows70..134 only
 sample={c:5,**dict(zip(b,range(1,10)))};En=E.subs(sample);ep=list(En.to_DM().rref()[1]);er=list(En[:,ep].T.to_DM().rref()[1]);normdim=2 if kind=='M0'else 1;chosen=None
 for inds in itertools.combinations(range(T.cols),normdim):
  H,mons=coeffs([{mon:v.subs(sample)for mon,v in f.items()}for f in F],T[:,list(inds)].subs(sample),cols);Yn=H-A.subs(sample)*U.subs(sample)*H[pr,:];Yn=Yn[70:,:];res=Yn-En[:,ep]*En[er,ep].inv(method='DM')*Yn[er,:]
  rows=list(res.T.to_DM().rref()[1])
  if len(rows)==len(mons):chosen=list(inds);break
 assert chosen is not None;fullframe=T.copy();ff=r['rank_certificates']['135']['free_columns'];fc=[i for i in range(44)if i not in pc];frame_rows=[fc[j]if j<len(fc)else 44 for j in ff];assert fullframe[frame_rows,:]==s.eye(fullframe.cols);T=T[:,chosen];den=[]
 for j in range(T.cols):
  dd=s.factor(s.lcm([s.denom(v)for v in T[:,j]]));den.append(dd)
  for i in range(T.rows):T[i,j]=s.cancel(T[i,j]*dd)
 from sympy.polys.rings import ring
 RR,*gens=ring(','.join(map(str,(c,*b))),s.QQ);TP=[[RR.from_expr(T[i,j])for j in range(T.cols)]for i in range(T.rows)];HP=[[RR.zero for mon in mons]for row in range(135)];forms={i:TP[cols.index(i)]for i in cols};forms[fi]=[RR.zero]*T.cols
 for row,f in enumerate(F):
  for mon,vv in f.items():
   aa,bb=[forms[i]for i in mon]+[TP[-1]]*(2-len(mon));vp=RR.from_expr(vv)
   for j,(uu,ww)in enumerate(mons):HP[row][j]+=vp*(aa[uu]*bb[ww]+(aa[ww]*bb[uu]if uu!=ww else RR.zero))
 AU=(A*U).applyfunc(s.cancel);YP=[];yden=[]
 for i in range(70,135):
  yd=s.factor(s.lcm([s.denom(v)for v in AU[i,:]]));yden.append(yd);yp=RR.from_expr(yd);ap=[RR.from_expr(s.cancel(v*yd))for v in AU[i,:]];YP.append([yp*HP[i][j]-sum((a*HP[r][j]for a,r in zip(ap,pr)),RR.zero)for j in range(len(mons))])
 ad,bd=E[er,ep].T.to_DM().unify(E[rows,ep].T.to_DM());ad=ad.to_field();bd=bd.convert_to(ad.domain);Cr=ad.lu_solve(bd).to_Matrix().T.applyfunc(s.cancel);result=s.zeros(len(rows),len(mons))
 for i,row in enumerate(rows):
  dd=s.factor(s.lcm([yden[row]]+[s.denom(s.cancel(v/yden[r]))for v,r in zip(Cr[i,:],er)]));lead=RR.from_expr(s.cancel(dd/yden[row]));cc=[RR.from_expr(s.cancel(v*dd/yden[r]))for v,r in zip(Cr[i,:],er)]
  for j,(aa,bb)in enumerate(mons):
   pn=lead*YP[row][j]-sum((a*YP[r][j]for a,r in zip(cc,er)),RR.zero);denom=s.factor(dd*den[aa]*den[bb]);co,fac=s.factor_list(denom);left=s.sympify(co)
   for ff,ee in fac:
    fp=RR.from_expr(ff);unused=0
    for h in range(ee):
     qq,rem=pn.div(fp)
     if not rem:pn=qq
     else:unused=ee-h;break
    left*=ff**unused
   result[i,j]=s.factor(pn.as_expr()/left)
 dm=result.to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));jm=E[er,ep].to_DM();jdet=s.factor(jm.domain.to_sympy(jm.det()));A0=A[:70,:];dm=A0[pr,pc].to_DM();firstdet=s.factor(dm.domain.to_sympy(dm.det()));assert (A0.to_DM()*C.to_DM()).is_zero_matrix;ct=(tc-A*U*tc[pr,:]).applyfunc(s.cancel);assert ct[:70,:]==s.zeros(70,1);directE=(A*C).row_join(ct)[70:,:].applyfunc(s.cancel);assert (directE-E).applyfunc(s.cancel)==s.zeros(*E.shape)
 L=s.zeros(len(rows),135)
 for i,row in enumerate(rows):
  L[i,70+row]=1
  for j,rr in enumerate(er):L[i,70+rr]-=Cr[i,j]
 correction=(L*A*U).applyfunc(s.cancel)
 for j,col in enumerate(pr):L[:,col]-=correction[:,j]
 assert (L.to_DM()*A.row_join(tc).to_DM()).is_zero_matrix
 all_den=s.factor(s.lcm([s.denom(v)for v in fullframe]+[s.denom(v)for v in C]+[s.denom(v)for v in Cr]+[s.denom(v)for v in L]));print(kind,'UNIFORM',time.time()-st,flush=True)
 return dict(kind=kind,beta=list(map(str,bv)),frame_substitution={str(k):str(v)for k,v in sub.items()},q_point=list(map(str,q)),chart_columns=cols,full_kernel_frame=enc(fullframe),frame_identity_rows=frame_rows,normal_columns=chosen,normal_denominators=list(map(str,den)),normal_cokernel=enc(L),normal_minor=enc(result),normal_determinant=str(det),jacobian_residual_rows=er,jacobian_residual_columns=ep,jacobian_residual_determinant=str(jdet),first_minor=str(firstdet),all_frame_denominator_factors=str(all_den),first_pivot_rows=pr,first_pivot_columns=pc,normal_rows=rows,all135_highest_verified=True)
