"""Exact quartic initial-ideal bound in a specified formal tangent chart."""
import time
from itertools import combinations_with_replacement
import sympy as s
from s61_flag_transverse_obstruction_20261101 import specialize
def quartic_bound(m,local,branches,expected):
    r=local;order_bound=4
    start=time.time();p=list(map(s.Rational,r['point61']));beta={i:p[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};_,F=specialize(m,beta,True);cols=r['chart_columns'];fixed=next(i for i in range(45) if i not in cols);J=s.Matrix(r['jacobian']).applyfunc(s.Rational);T=s.Matrix(r['tangent_basis']).applyfunc(s.Rational);n=T.cols;rows=list(J.T.to_DM().rref()[1]);piv=list(J.to_DM().rref()[1]);U=s.zeros(45,len(piv));iv=J.extract(rows,piv).inv(method='DM')
    for i,c in enumerate(piv):U[c,:]=iv[i,:]
    jac_piv=list(piv);JU=J*U;mon2=list(combinations_with_replacement(range(n),2));mon3=list(combinations_with_replacement(range(n),3));i3={x:i for i,x in enumerate(mon3)};Q=s.zeros(135,len(mon2));terms=[]
    for row,f in enumerate(F):
     ts=[]
     for mon,c in f.items():
      if fixed in mon:continue
      i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
      for h,(a,b) in enumerate(mon2):Q[row,h]+=c*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
     terms.append(ts)
    w2=-U*Q[rows,:];A2=Q-JU*Q[rows,:];Q3=s.zeros(135,len(mon3))
    for row,ts in enumerate(terms):
     for i,j,c in ts:
      for a in range(n):
       for h,mon in enumerate(mon2):Q3[row,i3[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w2[j,h]+T[j,a]*w2[i,h])
    if order_bound==4:
        A3=Q3-JU*Q3[rows,:];w3=-U*Q3[rows,:];mon4=list(combinations_with_replacement(range(n),4));i4={x:i for i,x in enumerate(mon4)};Q4=s.zeros(135,len(mon4))
        for row,ts in enumerate(terms):
         for i,j,c in ts:
          for h,ma in enumerate(mon2):
           for l,mb in enumerate(mon2):Q4[row,i4[tuple(sorted(ma+mb))]]+=c*w2[i,h]*w2[j,l]
          for a in range(n):
           for h,mon in enumerate(mon3):Q4[row,i4[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w3[j,h]+T[j,a]*w3[i,h])
        A4=Q4-JU*Q4[rows,:];mons=mon2+mon3+mon4;ix={mon:i for i,mon in enumerate(mons)};rr,rp=s.Matrix.hstack(A2,A3,A4).to_DM().rref();R=rr.to_Matrix()[:len(rp),:];mult=[]
        for row in R.tolist():
         for fac in [(i,) for i in range(n)]+mon2:
          z=[s.S.Zero]*len(mons)
          for j,mon in enumerate(mons):
           if len(mon)+len(fac)<=4 and row[j]:z[ix[tuple(sorted(mon+fac))]]+=row[j]
          if any(z):mult.append(z)
        full=R.col_join(s.Matrix(mult));print('macaulay',full.shape,flush=True);rr,piv=full.to_DM().rref();rr=rr.to_Matrix();v=s.symbols(f'v0:{n}');forms=[];ranks={2:0,3:0,4:0}
        for row,pc in enumerate(piv):
         degree=len(mons[pc]);ranks[degree]+=1;forms.append(s.expand(sum(rr[row,j]*s.prod(v[i] for i in mon) for j,mon in enumerate(mons) if len(mon)==degree)))
        ar=list(range(ranks[2]));br=list(range(ranks[3]));print('initial ranks',ranks,flush=True)
    else:
        A3=Q3-JU*Q3[rows,:];ar=list(A2.T.to_DM().rref()[1]);ac=list(A2.to_DM().rref()[1]);C=A2[:,ac]*A2.extract(ar,ac).inv(method='DM') if ar else s.zeros(135,0);B=A3-C*A3[ar,:];br=list(B.T.to_DM().rref()[1]);v=s.symbols(f'v0:{n}');forms=[s.expand(sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,mons))) for M,mons in [(A2[ar,:],mon2),(B[br,:],mon3)] for row in M.tolist()]
    choices=[[1]+[0]*(n-1),list(range(1,n+1)),[1]*(n)]
    result=None
    for coeff in choices:
     substitute=-sum(coeff[i]*v[i] for i in range(1,n));section=s.groebner([f.subs(v[0],substitute) for f in forms],*v[1:],order='grevlex')
     if not section.is_zero_dimensional:continue
     gb=s.groebner([f.subs(v[0],1+substitute) for f in forms],*v[1:],order='grevlex');assert gb.is_zero_dimensional;leads=[p.LM(order=gb.order).exponents for p in gb.polys];counts=[]
     for degree in range(20):
      count=0
      for inds in combinations_with_replacement(range(n-1),degree):
       e=tuple(inds.count(i) for i in range(n-1))
       if not any(all(a<=b for a,b in zip(l,e)) for l in leads):count+=1
      counts.append(count)
      if count==0:break
     result={'parameter_coefficients':coeff,'multiplicity_upper_bound':sum(counts),'localized_basis':list(map(str,gb)),'standard_degrees':counts,'initial_forms':list(map(str,forms)),'section_groebner':list(map(str,section))};break
    assert result is not None and result['multiplicity_upper_bound']==expected
    coeff=result['parameter_coefficients'];values=[]
    for branch in branches:
        value=sum(s.Rational(a)*b for a,b in zip(branch['tangent'],coeff));assert value!=0;values.append(str(value))
    result.update(initial_degree_bound=4,parameter_on_branch_tangents=values,jacobian_pivot_rows=rows,jacobian_pivot_columns=jac_piv,seconds=time.time()-start)
    return result
