#!/usr/bin/env python3
"""A binary/mixed-tail junction: actual infinity line and transverse length upper bound67."""
import sys,json,time,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import sympy as s,numpy as np
from itertools import combinations_with_replacement
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta


def verify_point():
    sample=json.loads((ROOT/"reports/S61_BINARY_UNIFORM_RANK_20261101.json").read_text())["sample"];sample.update({"u4":"0","u14":"42/5","u50":"-55"});expected_quadratic_rank=23;expected_parameters=[0,1,3]
    start=time.time();p=32003;m=S61()
    r={'sample':dict(sample)}
    b={int(k[1:])-1:s.Rational(v) for k,v in r['sample'].items()};two,F=specialize(m,b,True);fixed=two.index(29);cols=[i for i in range(45) if i!=fixed];q0=s.zeros(45,1);q0[fixed]=1
    J=s.zeros(135,45)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if len(mon)==2:
       i,j=mon
       if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
       if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
      elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
    T=s.Matrix.hstack(*J.nullspace());n=T.cols;assert n==12 and J*T==s.zeros(135,12) and T[-1,:]==s.zeros(1,12)
    _,piv=J.to_DM().rref();free=[i for i in range(45) if i not in piv];assert T[free,:]==s.eye(n)
    z,uu,bf=primitive_beta(m);bf=bf.subs({s.Symbol(k):s.Rational(v) for k,v in r['sample'].items()});a,c,d,h,k=s.symbols('a c d h k');v1=s.Matrix([a,1,0,c]);v2=s.Matrix([-c,0,1,d]);origin={a:0,c:0,d:0,h:0,k:0}
    sec=s.Matrix([(v1.T*bf.subs(dict(zip(z,v1)))*v2)[0],(v1.T*bf.subs(dict(zip(z,v2)))*v2)[0]])
    assert sec.subs(origin)==s.zeros(2,1)
    curve_jac=sec.jacobian([a,c,d]).subs(origin);assert curve_jac.rank()==2
    direction=curve_jac.nullspace()[0]
    x,y,zz,w=X[1:];ell1=x-a*y+c*zz;ell2=w-c*y-d*zz;fam=(ell1*ell2+h*ell1**2+k*ell2**2)*(v1*v2.T-v2*v1.T)
    ds=[]
    for v in [a,c,d,h,k]:
     P=X[0]**2*J0;P[1:,1:]+=fam.diff(v).subs(origin);qs=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);qs-=q0*qs[fixed];ds.append(qs[cols,:].col_join(s.zeros(1,1)))
    family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4]);coordinates=family[free,:];assert T*coordinates==family and J*family==s.zeros(135,3)
    _,params=coordinates.T.to_DM().rref();params=list(params);normal=[i for i in range(n) if i not in params];assert params==expected_parameters and len(normal)==9
    def mod(v):v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,p)%p
    def arr(M):return np.array([[mod(v) for v in row] for row in M.tolist()],dtype=np.int64)
    def rref(A):
     A=A.copy()%p;row=0;piv=[]
     for col in range(A.shape[1]):
      choices=np.flatnonzero(A[row:,col])
      if len(choices)==0:continue
      j=row+int(choices[0]);A[[row,j]]=A[[j,row]];A[row]=A[row]*pow(int(A[row,col]),-1,p)%p
      mask=np.flatnonzero(A[:,col]);mask=mask[mask!=row]
      if len(mask):A[mask]=(A[mask]-A[mask,col,None]*A[row,None,:])%p
      piv.append(col);row+=1
      if row==A.shape[0]:break
     return A,piv
    j=arr(J);t=arr(T);_,C=rref(j);_,R=rref(j[:,C].T);rank=len(C);rr,_=rref(np.concatenate([j[np.ix_(R,C)],np.eye(rank,dtype=np.int64)],axis=1));U=np.zeros((45,rank),dtype=np.int64);U[C,:]=rr[:,rank:];JU=j@U%p
    assert len(C)==len(R)==33 and np.array_equal(j[R,:]@U%p,np.eye(33,dtype=np.int64))
    H=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      i,l=[cols.index(v) for v in mon]+[44]*(2-len(mon));v=mod(c)
      if i==l:H[row,i,l]+=v
      else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p

    # Exact QQ quadratic rank in the nine transverse variables.
    Tn=T[:,normal];exact_mons=list(combinations_with_replacement(range(9),2));em={mon:i for i,mon in enumerate(exact_mons)};ecache={}
    def eprod(i,j):
     key=tuple(sorted((i,j)))
     if key not in ecache:
      out=[s.S.Zero]*45
      for a,v in enumerate(Tn[i,:]):
       if not v:continue
       for bb,w in enumerate(Tn[j,:]):
        if w:out[em[tuple(sorted((a,bb)))]]+=v*w
      ecache[key]=out
     return ecache[key]
    Qe=[]
    for f in F:
     out=[s.S.Zero]*45
     for mon,c in f.items():
      if fixed in mon:continue
      i,jj=[cols.index(v) for v in mon]+[44]*(2-len(mon))
      for k,v in enumerate(eprod(i,jj)):out[k]+=c*v
     Qe.append(out)
    coe=s.Matrix.hstack(*J.T.nullspace()).T;Ae=(coe.to_DM()*s.Matrix(Qe).to_DM()).to_Matrix();assert len(Ae.to_DM().rref()[1])==expected_quadratic_rank

    H%=p;t=t[:,normal];n=len(normal)
    m2=list(combinations_with_replacement(range(n),2));m3=list(combinations_with_replacement(range(n),3));m4=list(combinations_with_replacement(range(n),4));i3={mon:i for i,mon in enumerate(m3)};i4={mon:i for i,mon in enumerate(m4)}
    Q=np.zeros((135,len(m2)),dtype=np.int64)
    for row in range(135):
     q=t.T@(H[row]@t%p)%p
     for k,(a,c) in enumerate(m2):Q[row,k]=q[a,c]*(1 if a==c else 2)%p
    w2=-U@Q[R,:]%p;A2=(Q-JU@Q[R,:])%p
    C3=np.zeros((135,len(m3)),dtype=np.int64)
    for row in range(135):
     cross=2*t.T@(H[row]@w2%p)%p
     for i in range(n):
      for k,mon in enumerate(m2):C3[row,i3[tuple(sorted((i,)+mon))]]+=cross[i,k]
    C3%=p;w3=-U@C3[R,:]%p;A3=(C3-JU@C3[R,:])%p
    C4=np.zeros((135,len(m4)),dtype=np.int64)
    for row in range(135):
     cross=2*t.T@(H[row]@w3%p)%p;square=w2.T@(H[row]@w2%p)%p
     for i in range(n):
      for k,mon in enumerate(m3):C4[row,i4[tuple(sorted((i,)+mon))]]+=cross[i,k]
     for i,mi in enumerate(m2):
      for k,mk in enumerate(m2):C4[row,i4[tuple(sorted(mi+mk))]]+=square[i,k]
    C4%=p;w4=-U@C4[R,:]%p;A4=(C4-JU@C4[R,:])%p
    base=np.concatenate([A2,A3,A4],axis=1);ncols=base.shape[1];rows=list(base)
    for row2,row3 in zip(A2,A3):
     for a in range(n):
      out=np.zeros(ncols,dtype=np.int64)
      for mon,v in zip(m2,row2):out[len(m2)+i3[tuple(sorted((a,)+mon))]]+=v
      for mon,v in zip(m3,row3):out[len(m2)+len(m3)+i4[tuple(sorted((a,)+mon))]]+=v
      if np.any(out):rows.append(out%p)
     for ma in m2:
      out=np.zeros(ncols,dtype=np.int64)
      for mon,v in zip(m2,row2):out[len(m2)+len(m3)+i4[tuple(sorted(ma+mon))]]+=v
      if np.any(out):rows.append(out%p)
    rels=np.array(rows);rr,piv=rref(rels);rk=len(piv);length=1+n+ncols-rk
    trow=np.concatenate([w2[-1,:],w3[-1,:],w4[-1,:]])%p
    rankt=len(rref(np.concatenate([rr[:rk,:],trow[None,:]],axis=0))[1])
    tmult=[trow]
    for a in range(n):
     out=np.zeros(ncols,dtype=np.int64)
     for mon,v in zip(m2,w2[-1,:]):out[len(m2)+i3[tuple(sorted((a,)+mon))]]+=v
     for mon,v in zip(m3,w3[-1,:]):out[len(m2)+len(m3)+i4[tuple(sorted((a,)+mon))]]+=v
     tmult.append(out%p)
    for ma in m2:
     out=np.zeros(ncols,dtype=np.int64)
     for mon,v in zip(m2,w2[-1,:]):out[len(m2)+len(m3)+i4[tuple(sorted(ma+mon))]]+=v
     tmult.append(out%p)
    rankideal=len(rref(np.concatenate([rr[:rk,:],np.array(tmult)],axis=0))[1])
    length3=1+n+len(m2)-len(rref(rels[:,:len(m2)])[1]);length4=1+n+len(m2)+len(m3)-len(rref(rels[:,:len(m2)+len(m3)])[1])

    assert length==67 and rk==648 and rankt==649 and rankideal==649
    expected_hilbert=[1,9,22,35];assert [1,n,length3-1-n,length4-length3,length-length4]==expected_hilbert+[0]
    # Save nonzero minors of the actual nonhomogeneous relation matrix, and of
    # that matrix with t appended. Their characteristic-zero use is in the note.
    _,selected_rows=rref(rels.T);assert len(selected_rows)==648
    selected_columns=piv
    minor648=rels[np.ix_(selected_rows,selected_columns)].copy()
    def determinant_mod(M):
     M=M.copy();det=1
     for col in range(len(M)):
      ii=col+int(np.flatnonzero(M[col:,col])[0])
      if ii!=col:M[[col,ii]]=M[[ii,col]];det=-det
      value=int(M[col,col]);det=det*value%p;M[col]=M[col]*pow(value,-1,p)%p
      mask=np.flatnonzero(M[col+1:,col])+col+1
      if len(mask):M[mask]=(M[mask]-M[mask,col,None]*M[col,None,:])%p
     return det%p
    det648=determinant_mod(minor648);assert det648

    # Initial cubics lift to QQ because the exact and modular quadratic ranks agree.
    initial,ipiv=rref(np.concatenate([A2,A3],axis=1));assert sum(i<45 for i in ipiv)==expected_quadratic_rank
    initialq=initial[:expected_quadratic_rank,:45];initialc=initial[expected_quadratic_rank:len(ipiv),45:]
    ms={d:list(combinations_with_replacement(range(9),d)) for d in range(6)};idx5={mon:i for i,mon in enumerate(ms[5])}
    mat5=[]
    for d,gens in [(2,initialq),(3,initialc)]:
     for row in gens:
      for mul in ms[5-d]:
       out=np.zeros(1287,dtype=np.int64)
       for mon,c in zip(ms[d],row):out[idx5[tuple(sorted(mon+mul))]]+=c
       mat5.append(out%p)
    mat5=np.array(mat5);_,piv5=rref(mat5);_,rows5=rref(mat5.T);assert len(piv5)==len(rows5)==1287
    quintic_minor=determinant_mod(mat5[rows5,:]);assert quintic_minor

    aug=np.concatenate([rels[selected_rows,:],trow[None,:]],axis=0)
    _,augcols=rref(aug);assert len(augcols)==649
    det649=determinant_mod(aug[:,augcols]);assert det649
    # A regular binary point on the SAME fixed-beta threefold: Q1=xw+x^2+w^2.
    P=X[0]**2*J0;xx,yy,zz,ww=X[1:];P[2,3]+=xx*ww+xx*xx+ww*ww;P[3,2]-=xx*ww+xx*xx+ww*ww
    q1=s.Matrix([m.normalize(P)['parameters'][i] for i in two]);assert q1[fixed]==1
    J1=s.zeros(135,45)
    for row,f in enumerate(F):
     value=0
     for mon,c in f.items():
      if len(mon)==2:
       i,jj=mon;value+=c*q1[i]*q1[jj]
       if i!=fixed:J1[row,cols.index(i)]+=c*q1[jj]
       if jj!=fixed:J1[row,cols.index(jj)]+=c*q1[i]
      elif len(mon)==1:J1[row,44]+=c*q1[mon[0]]
     assert value==0
    T1=s.Matrix.hstack(*J1.nullspace());assert T1.shape==(45,10) and T1[-1,:]==s.zeros(1,10)
    mons1=list(combinations_with_replacement(range(10),2));idx1={mon:i for i,mon in enumerate(mons1)};cache={}
    def product1(i,j):
     key=tuple(sorted((i,j)))
     if key not in cache:
      vals=[s.S.Zero]*55
      for a,v in enumerate(T1[i,:]):
       if not v:continue
       for bb,w in enumerate(T1[j,:]):
        if w:vals[idx1[tuple(sorted((a,bb)))]]+=v*w
      cache[key]=vals
     return cache[key]
    Q1=[]
    for f in F:
     row=[s.S.Zero]*55
     for mon,c in f.items():
      if fixed in mon:continue
      i,jj=[cols.index(v) for v in mon]+[44]*(2-len(mon))
      for k,v in enumerate(product1(i,jj)):row[k]+=c*v
     Q1.append(row)
    co1=s.Matrix.hstack(*J1.T.nullspace()).T;A1=(co1.to_DM()*s.Matrix(Q1).to_DM()).to_Matrix();_,goodrows=A1.T.to_DM().rref();assert len(goodrows)==5
    hess=[]
    for row in A1[list(goodrows),:].tolist():
     M=s.zeros(10)
     for c,(i,jj) in zip(row,mons1):M[i,jj]+=c;M[jj,i]+=c
     hess+=M.tolist()
    assert len(s.Matrix(hess).to_DM().rref()[1])==7
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    return dict(sample=r['sample'],linear_parameters=params,normal_variables=normal,
     actual_binary_family_tangent=enc(coordinates),parameter_minor=str(coordinates[params,:].det()),curve_derivative=enc(curve_jac),
     jacobian_pivot_rows=R,jacobian_pivot_columns=C,exact_transverse_quadratic_rank=expected_quadratic_rank,
     quintic_initial_minor_rows=rows5,quintic_initial_minor_mod_prime=quintic_minor,
     residual_quadratic_mod_prime=A2.tolist(),residual_cubic_mod_prime=A3.tolist(),residual_quartic_mod_prime=A4.tolist(),
     homogenizer_quadratic_mod_prime=w2[-1,:].tolist(),homogenizer_cubic_mod_prime=w3[-1,:].tolist(),homogenizer_quartic_mod_prime=w4[-1,:].tolist(),
     relation_matrix_shape=list(rels.shape),relation_minor_rows=selected_rows,relation_minor_columns=selected_columns,relation_minor_mod_prime=det648,
     augmented_homogenizer_minor_columns=augcols,augmented_homogenizer_minor_mod_prime=det649,
     same_beta_regular_point_q=list(map(str,q1)),same_beta_regular_full_rank=35,same_beta_regular_quadratic_rank=5,same_beta_regular_hessian_radical=3,
     characteristic_zero_transverse_length_upper_bound=67,modular_transverse_length=67,modular_transverse_hilbert_function=expected_hilbert,
     characteristic_zero_transverse_fifth_power_zero=True,modular_homogenizer_nonzero=True,modular_homogenizer_ideal_dimension=1,
     full_completed_ring_nilpotence_proved=False,
     seconds=time.time()-start)


def verify():
    start=time.time();point=verify_point();m=S61();u=s.symbols('u1:62')
    one=[i for i,v in enumerate(m.reduced['parameter_weights']) if v==1]
    b={i:u[i] for i in one}
    b.update({3:0,11:0,12:0,13:s.Rational(14,25)*u[48],49:-5*u[34]})
    two,F=specialize(m,b,True);x,y,z,w=X[1:];a=s.Symbol('a')
    ell=z-25*u[4]*w/(3*u[48])-5*(5*u[25]-u[50])*x/(6*u[48])
    G=x*w*w/2+a*ell**3
    P=X[0]**2*J0
    P[2,3]+=s.diff(G,w);P[3,2]-=s.diff(G,w)
    P[2,4]-=s.diff(G,z);P[4,2]+=s.diff(G,z)
    q=s.Matrix([m.normalize(P)['parameters'][i] for i in two])
    for f in F:
        value=sum(c*q[mon[0]]*q[mon[1]] for mon,c in f.items() if len(mon)==2)
        assert s.cancel(s.expand(value))==0
    # The whole eleven-dimensional source is stable under the fixed Q stabilizer.
    family_path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    rank_path=ROOT/'reports/S61_BINARY_UNIFORM_RANK_20261101.json'
    family=json.loads(family_path.read_text());rank=json.loads(rank_path.read_text())
    uu=list(map(s.Symbol,family['sample'].keys()))
    action=s.Matrix(family['primitive_action']).applyfunc(s.sympify)
    stab=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in rank['primitive_actions']])
    cons=s.Matrix([u[11],u[12],u[3],5*u[34]+u[49],25*u[13]-14*u[48]])
    gb=s.groebner(list(cons),*uu,domain=s.QQ)
    assert all(gb.reduce(s.expand(v))[1]==0 for v in cons.jacobian(uu)*stab)
    sample={s.Symbol(k):s.Rational(v) for k,v in point['sample'].items()}
    assert cons.subs(sample)==s.zeros(5,1)
    normal=(cons.jacobian(uu)*action).subs(sample)
    columns=list(normal.to_DM().rref()[1]);assert len(columns)==5
    determinant=normal[:,columns].det();assert determinant!=0
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    paths=[Path(__file__),family_path,rank_path,
        ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
        ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
        ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
        ROOT/'reports/S61_BINARY_UNIFORM_CUBIC_20261101.json',
        ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_BINARY_MIXED_CONNECTION_20261101',all_passed=True,prime=32003,
        point=point,source_equations=list(map(str,cons)),source_dimension=11,
        joint_image_dimension=16,beta_dominance_minor=str(determinant),
        beta_normal_action=enc(normal),beta_minor_columns=columns,
        universal_tail_potential=str(G),tail_linear_form=str(ell),
        all_135_leading_equations_identically_zero=True,
        proof_note='research/notes/S61_BINARY_MIXED_CONNECTION_20261101.md',
        limitations=['The full completed ring at the connection has not been identified.',
            'The transverse upper bound67 does not imply a free rank65 ring over the binary threefold.',
            'The mixed tail is an actual additional reduced branch; its global saturation is separate.',
            'The modular Hilbert function and nonzero t are not asserted over QQ without a lower bound.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path);args=parser.parse_args();report=verify()
    with args.output.open('x') as out:json.dump(report,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,beta_minor=report['beta_dominance_minor'],
        point={k:report['point'][k] for k in ['exact_transverse_quadratic_rank',
            'quintic_initial_minor_mod_prime','relation_minor_mod_prime',
            'augmented_homogenizer_minor_mod_prime','modular_transverse_hilbert_function']},
        seconds=report['seconds'])))
