#!/usr/bin/env python3
"""Remove the semisimple minimal-tail binary endpoint A using good reduction."""
import argparse,hashlib,json,os,shlex,shutil,subprocess,tempfile,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
import numpy as np
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_minimal_baseline_thickenings_20261101 import coords,enc
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_infinity_local_data_20261101 import calculate
from weight1_structure import sp_basis
from truncated_local_ideal_io_20261101 import write_input,read_initial
ROOT=Path(__file__).resolve().parents[1]
PRIME=32003
def point_jets():
    sample=dict(zip(('u4','u5','u12','u13','u14','u25','u26','u27','u34','u35','u36','u37','u48','u49','u50','u51'),map(str,(0,0,0,0,-1,2,2,3,5,3,7,11,13,10,-15,5))));expected_parameters=[0,1,3]
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
    T=s.Matrix.hstack(*J.nullspace());n=T.cols;assert n==11 and J*T==s.zeros(135,11) and T[-1,:]==s.zeros(1,11)
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
    _,params=coordinates.T.to_DM().rref();params=list(params);normal=[i for i in range(n) if i not in params];assert params==expected_parameters and len(normal)==8
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
    assert len(C)==len(R)==34 and np.array_equal(j[R,:]@U%p,np.eye(34,dtype=np.int64))
    H=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      i,l=[cols.index(v) for v in mon]+[44]*(2-len(mon));v=mod(c)
      if i==l:H[row,i,l]+=v
      else:H[row,i,l]+=v*pow(2,-1,p)%p;H[row,l,i]+=v*pow(2,-1,p)%p

    normal=list(range(T.cols))
    H%=p;n=len(normal);W={1:t[:,normal]};A={};ms={d:list(combinations_with_replacement(range(n),d)) for d in range(1,6)};idx={d:{m:i for i,m in enumerate(ms[d])} for d in range(1,6)}
    for degree in range(2,6):
     Qd=np.zeros((135,len(ms[degree])),dtype=np.int64)
     for da in range(1,degree//2+1):
      db=degree-da;fac=1 if da==db else 2
      indices=np.array([idx[degree][tuple(sorted(ma+mb))] for ma in ms[da] for mb in ms[db]])
      for row in range(135):
       cross=(W[da].T@(H[row]@W[db]%p))%p
       acc=np.bincount(indices,weights=cross.ravel(),minlength=len(ms[degree])).astype(np.int64)
       Qd[row]=(Qd[row]+fac*acc)%p
     W[degree]=-U@Qd[R,:]%p;A[degree]=(Qd-JU@Qd[R,:])%p
     print('degree',degree,'computed',time.time()-start,flush=True)
    return dict(sample=sample,J=J,T=T,F=F,two=two,fixed=fixed,chart=cols,free=free,rows=R,columns=C,parameters=params,family=coordinates,curve_jacobian=curve_jac,A=A,tseries={d:W[d][-1,:].tolist() for d in A},monomials=ms)

def source_geometry():
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2]
    b=s.symbols('b0:10');e=b[0];bn=list(b[1:7]);ut,vt,ww=b[7:]
    beta=s.Matrix([0,0,0,0,-e/10,bn[0],ut,bn[1],bn[2],vt,bn[3],bn[4],bn[5],e,-5*vt,ww])
    x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=x*w;Q-=Q.T;q0=s.Matrix(coords(m,Q))[two,:];fi=two.index(29)
    delta0=s.zeros(45,1);delta0[0]=1;delta0[19]=-4*(5*ut+2*ww)/(15*e);delta0[25]=-4*vt/(5*e);delta0[28]=-s.Rational(4,5)
    delta1=s.zeros(45,1);delta1[0]=s.Rational(1,2);delta1[19]=(5*ut+7*ww)/(15*e);delta1[25]=8*vt/(15*e);delta1[28]=-s.Rational(3,10);delta1[35]=s.Rational(1,2)
    assert s.Matrix.hstack(delta0,delta1)[[0,35],:].det()==s.Rational(1,2)
    _,F=specialize(m,dict(zip(one,beta)),True)
    for f in F:
        for a,bb in [(q0,q0),(q0,delta0),(delta0,delta0),(q0,delta1),(delta1,delta1)]:
            assert s.cancel(sum(c*(a[mon[0]]*bb[mon[1]]+bb[mon[0]]*a[mon[1]]) for mon,c in f.items() if len(mon)==2))==0
    orbit=[]
    for A in sp_basis():
        dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q0*dq[fi];orbit.append(dq)
    O=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*O.nullspace());assert O.rank()==5 and stab.shape==(10,5)
    uu,action=beta_actions(m);H=beta.jacobian(b);act=action.subs(dict(zip(uu,beta)));ann=s.Matrix.hstack(*H.T.nullspace()).T;normal=(ann*act*stab).applyfunc(s.factor)
    assert normal.to_DM().rank()==1
    sample=dict(zip(b,[10,2,3,5,7,11,13,2,3,5]));D=H.row_join(act);Ds=D.subs(sample);pc=list(Ds.to_DM().rref()[1]);assert len(pc)==16;det=Ds[:,pc].det();assert det
    joint=D.col_join(s.zeros(45,10).row_join(O));assert joint.subs(sample).to_DM().rank()==16
    # Both k=-11/12 tail endpoints, before any specialization of nuisance beta.
    old=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text())
    K,L,U,Vv=s.symbols('K L U V');ns=s.symbols('n0:6')
    src=s.Matrix([0,0,0,0,-s.Rational(11,12),ns[0],U,ns[1],ns[2],Vv,ns[3],ns[4],ns[5],-s.Rational(5,6),K-5*Vv,5*U-L])
    T=s.Matrix([[1,0,0,0],[0,1,-4*L/5,0],[4*L/5,0,0,-1],[0,0,1,0]])
    assert T*J0[1:,1:]*T.T==J0[1:,1:]
    zz,ub,B=primitive_beta(m);from weight1_structure import COORDS
    vec=s.Matrix([B[i,j].diff(zz[h]) for (i,j),h in COORDS]);basis=vec.jacobian(ub);rows=list(basis.T.to_DM().rref()[1]);left=basis[rows,:].inv()
    Bt=B.subs(dict(zip(ub,src)));Bp=(T.T*Bt.subs(dict(zip(zz,T*s.Matrix(zz))),simultaneous=True)*T).applyfunc(s.expand)
    nv=(left*s.Matrix([Bp[i,j].diff(zz[h]) for (i,j),h in COORDS])[rows,:]).applyfunc(s.factor)
    assert all(nv[i]==0 for i in [0,1,2,3])
    assert s.expand(10*nv[4]+nv[13])==0 and s.expand(5*nv[9]+nv[14])==0
    scale=-6*e/5
    sol=s.solve(list((nv-beta/scale).subs(U,0)),[K,L,Vv]+list(ns),dict=True);assert len(sol)==1;sol={U:s.S.Zero,**sol[0]}
    assert (nv.subs(sol)*scale-beta).applyfunc(s.cancel)==s.zeros(16,1)
    for kind,expected in [('M0',delta0),('M1',delta1)]:
        fam=next(v for v in old['families'] if v['label']==kind);tail=next(v for v in fam['rational_tails'] if v['axis']==2 and v['k']=='-11/12')
        g=s.symbols('g0:9');h=g[2]
        G=sum(s.sympify(tail['coordinate_formulas'].get(str(gg),str(gg)))*mon for gg,mon in zip(g,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]))
        Qtail=s.zeros(4);Qtail[0,1]=x*x if kind=='M0' else x*x/2;Qtail[2,3]=0 if kind=='M0' else x*x/2
        Qtail[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Qtail[1,3]=s.diff(G,z);Qtail-=Qtail.T
        Qt=(T.inv()*Qtail.subs(dict(zip(V,T*V)),simultaneous=True)*T.inv().T).applyfunc(s.expand)
        qs=s.Matrix(coords(m,Qt))[two,:];assert qs.diff(h)==-2*q0
        qc=qs.subs(h,0);dd=(qc-q0*qc[fi]).subs(sol).applyfunc(s.factor);assert (dd-expected).applyfunc(s.cancel)==s.zeros(45,1)
    return dict(beta=list(map(str,beta)),parameters=list(map(str,b)),open_condition='b0 != 0',actual_line_deltas=[list(map(str,delta0)),list(map(str,delta1))],all135_actual_line_identities=True,line_embedding_minor='1/2',projective_orbit_rank=5,stabilizer=enc(stab),source_annihilator=enc(ann),normal_stabilizer_matrix=enc(normal),normal_stabilizer_rank=1,beta_sample={str(k):str(v) for k,v in sample.items()},beta_minor_columns=pc,beta_minor=str(det),joint_rank=16,image_dimension=16,source_tail_normalization=enc(T),normalized_source_beta=list(map(str,nv)),source_recovery={str(k):str(v) for k,v in sol.items()},source_beta_scale=str(scale),covered_tails=['M0/M1 g2 k=-11/12','their Weyl partners M0/M1 g4 k=-1/12'])

PRIMARY_PROGRAM='''
assert(dim(R/J)==3 and degree(R/J)==65);
B=saturate(J,ideal(v_1));D=saturate(J,ideal(v_4));E=intersect(B,D);
D0=ideal(v_1,v_2,v_3,v_6,v_7,v_8,v_9,v_10,v_0-(2/3)*v_5,v_4-(3/10)*v_5);
D1=ideal(v_1,v_2,v_3,v_6,v_8,v_9,v_10,v_0-(3/5)*v_7,v_4-(8/25)*v_7,v_5+(3/5)*v_7);
assert(D==intersect(D0,D1) and degree(R/D)==2 and degree(R/B)==65);
assert(dim(E/J)==0 and degree(E/J)==1);assert((gens((ideal vars R)*E))%J==0);
assert(hilbertFunction(5,E/J)==1);
C=ideal(v_4,v_6,v_7,v_8,v_10,v_2-(3/2)*v_1,v_5-(15/8)*v_1,v_9-v_1);
assert((gens B)%C==0);assert all(flatten entries gens C,g->g^4%B==0);
assert(dim(R/(J+ideal(v_0,v_1,v_3)))==0 and degree(R/(J+ideal(v_0,v_1,v_3)))==67);
print("initial primary part65, two actual lines, single degree5 socle verified");stdio<<flush;
S=ZZ/32003[w_2,w_4,w_5,w_6,w_7,w_8,w_9,w_10,w_0,w_3,MonomialOrder=>{GRevLex=>8,GRevLex=>2}];
ff=map(S,R,{w_0,1_S,w_2,w_3,w_4,w_5,w_6,w_7,w_8,w_9,w_10});B1=ff B;GG=gens gb B1;cc={};
scan(flatten entries GG,g->(ee:=take((exponents leadTerm g)#0,8);mm:=product(8,j->S_j^(ee#j));cf:=(sum(select(terms g,tt->take((exponents tt)#0,8)==ee)))//mm;if first degree cf>0 then cc=append(cc,cf/(leadCoefficient cf))));
assert(cc=={});normalInitial=ideal apply(flatten entries GG,g->product(8,j->S_j^((exponents leadTerm g)#0#j)));
assert(dim(S/normalInitial)==2 and degree(S/normalInitial)==65);
print("C-primary constant-leading-coefficient normal length65 verified");print("CERTIFICATE PASS");stdio<<flush;
'''

def mod(v):
    q=s.Rational(v);assert int(q.q)%PRIME;return int(q.p)*pow(int(q.q),-1,PRIME)%PRIME

def rank(M):
    A=np.array([[mod(v) for v in row] for row in s.Matrix(M).tolist()],dtype=np.int64);j=0
    for c in range(A.shape[1]):
        ix=np.flatnonzero(A[j:,c])
        if not len(ix):continue
        a=j+int(ix[0]);A[[j,a]]=A[[a,j]];A[j]=A[j]*pow(int(A[j,c]),-1,PRIME)%PRIME
        inds=np.flatnonzero(A[:,c]);inds=inds[inds!=j];A[inds]=(A[inds]-A[inds,c,None]*A[j,None,:])%PRIME;j+=1
        if j==len(A):break
    return j

def run(m2_command=None):
    start=time.time();r=point_jets();geometry=source_geometry();m=S61();pt=[s.S.Zero]*61
    for k,v in r['sample'].items():pt[int(k[1:])-1]=s.Rational(v)
    pt[29]=1;local=calculate(m,True,s.Matrix(pt))
    assert(local['leading_jacobian_rank'],local['homogenized_jacobian_rank'],local['quadratic_span_rank'],local['common_hessian_rank'])==(33,34,14,9)
    assert rank(r['J'][:,:44])==33 and rank(r['J'])==34 and rank(r['curve_jacobian'])==2 and mod(r['family'][r['parameters'],:].det())
    # Same smooth binary threefold, moved only in its nonsquare quadratic coefficient.
    x,y,z,w=X[1:];Q=s.zeros(4);Q[1,2]=x*w+x*x+w*w;Q-=Q.T;p1=s.Matrix(coords(m,Q))
    for k,v in r['sample'].items():p1[int(k[1:])-1]=s.Rational(v)
    same=calculate(m,True,p1)
    assert(same['leading_jacobian_rank'],same['homogenized_jacobian_rank'],same['quadratic_span_rank'],same['common_hessian_rank'])==(34,35,5,7)
    assert rank(s.Matrix(same['jacobian'])[:,:44])==34 and rank(same['jacobian'])==35 and rank(same['quadratic_obstruction_basis'])==5
    hs=[]
    for row in s.Matrix(same['quadratic_obstruction_basis']).tolist():
        H=s.zeros(10)
        for c,(i,j) in zip(row,same['quadratic_monomials']):H[i,j]+=s.Rational(c);H[j,i]+=s.Rational(c)
        hs+=H.tolist()
    assert rank(s.Matrix(hs))==7
    dcs=[]
    for dd in geometry['actual_line_deltas']:
        delta=s.Matrix(dd).applyfunc(s.sympify).subs({s.Symbol(k):s.Rational(v) for k,v in geometry['beta_sample'].items()})
        dt=delta[r['chart'],:].col_join(s.zeros(1,1));dc=dt[r['free'],:];assert r['T']*dc==dt;dcs.append(list(map(str,dc)))
    assert dcs==[['-8/15','0','0','0','-6/25','-4/5','0','0','0','0','0'],['3/10','0','0','0','4/25','-3/10','0','1/2','0','0','0']]
    with tempfile.TemporaryDirectory(prefix='poisson-semi-A-') as tmp:
        tmp=Path(tmp);inp=tmp/'input.txt';basis=tmp/'basis.txt';binary=tmp/'closure'
        info=write_input(inp,PRIME,11,{d:A.tolist() for d,A in r['A'].items()},r['tseries'],maximum_degree=5)
        cc=subprocess.run(['c++','-O3','-std=c++17','-Wno-unqualified-std-cast-call',str(ROOT/'code/truncated_local_ideal_20261101.cpp'),'-o',str(binary)],capture_output=True,text=True,timeout=30);assert cc.returncode==0,cc.stderr
        cr=subprocess.run([str(binary),str(inp),str(basis)],capture_output=True,text=True,timeout=60);assert cr.returncode==0 and 'rank 3247 length 1121' in cr.stdout,cr.stdout+cr.stderr
        allinitial=read_initial(basis,11,5);assert len(allinitial)==3247
        selected=list(range(16))+[20,21,23,24,25,26,28,813,821,822,823,824]
        chosen=[allinitial[i] for i in selected];assert len(chosen)==28
        assert [sum(len(next(iter(f)))==d for f in chosen) for d in [2,3,4,5]]==[14,9,0,5]
        expr=lambda f:'+'.join(str(c)+'*'+'*'.join('v_'+str(i) for i in mon) for mon,c in f.items())
        program='R=ZZ/32003[v_0..v_10];J=ideal('+','.join(expr(f) for f in chosen)+');\n'+PRIMARY_PROGRAM
        source=tmp/'certificate.m2';source.write_text(program);env=os.environ.copy();local_binary=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary')
        if m2_command:cmd=shlex.split(m2_command)
        elif local_binary.exists():cmd=[str(local_binary)];env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        else:cmd=[shutil.which('M2') or 'M2']
        mr=subprocess.run(cmd+['--script',str(source)],env=env,capture_output=True,text=True,timeout=60);output=mr.stdout+mr.stderr;assert mr.returncode==0 and 'CERTIFICATE PASS' in output and 'error:' not in output,output
    prime_path=ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json';prime=json.loads(prime_path.read_text());assert prime['all_passed'] and prime['prime']==PRIME
    dependencies=['code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','code/s61_minimal_baseline_thickenings_20261101.py','code/s61_minimal_rankone_A1_removal_20261101.py','code/s61_minimal_infinity_local_data_20261101.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py','code/truncated_local_ideal_20261101.cpp','code/truncated_local_ideal_io_20261101.py','reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json','reports/S61_BINARY_UNIFORM_CUBIC_20261101.json','reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json','reports/S61_BINARY_UNIFORM_RANK_20261101.json','reports/S61_BINARY_PRIME_SPREADING_20261101.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json','research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md','research/notes/FINITE_BASE_TORSION_SPECIALIZATION_20261101.md']
    paths=[Path(__file__)]+[ROOT/f for f in dependencies]
    return dict(id='S61_MINIMAL_SEMI_A_ENDPOINT_20261101',all_passed=True,prime=PRIME,sample=r['sample'],source_geometry=geometry,point_local_data=local,same_binary_open_point=same,same_binary_modular_ranks=[34,35,5,7],curve_jacobian=enc(r['curve_jacobian']),actual_binary_family_tangent=enc(r['family']),binary_parameters=r['parameters'],actual_line_tangent_coordinates=dcs,jacobian_pivot_rows=r['rows'],jacobian_pivot_columns=r['columns'],truncated_input=info,truncated_closure_output=cr.stdout.splitlines()[:7],initial_selected_indices=selected,initial_generators=[[[list(mon),int(c)] for mon,c in f.items()] for f in chosen],initial_generator_degrees=[len(next(iter(f))) for f in chosen],macaulay2_program=program,macaulay2_output=output,initial_model='J subset B65 intersect two actual reduced lines; quotient k(-5)',normal_component_primary=True,transverse_initial_length=67,conclusion='The whole sufficiently general-beta fiber of the dominant16D semisimple endpoint A is removed via finite-base torsion.',proof_note='research/notes/S61_MINIMAL_SEMI_A_ENDPOINT_20261101.md',limitations=['The good-reduction/DVR argument proves finite-length torsion after characteristic-zero descent; it does not identify the full characteristic-zero tangent cone.','No relative ambient t^2=0 or freeness assertion.','Other binary/flag strata, full affine finiteness and five-point exhaustion remain unresolved.'],input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},seconds=time.time()-start)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--m2-command');args=parser.parse_args();report=run(args.m2_command)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
