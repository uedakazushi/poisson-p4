#!/usr/bin/env python3
"""Remove the pure-cubic minimal-tail binary endpoint using good reduction."""
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
    start=time.time();p=PRIME;m=S61();x,y,zz,w=X[1:]
    Q=s.zeros(4);Q[1,3]=zz**2+105*x*x;Q-=Q.T;pt=s.Matrix(coords(m,Q))
    one=[i for i,wg in enumerate(m.reduced['parameter_weights']) if wg==1]
    vals=[0,0,0,0,0,2,119,3,5,3,7,11,13,6,5,-35]
    sample={'u'+str(i+1):str(v) for i,v in zip(one,vals)}
    for i,v in zip(one,vals):pt[i]=v
    local=calculate(m,True,pt);assert(local['leading_jacobian_rank'],local['homogenized_jacobian_rank'],local['quadratic_span_rank'],local['common_hessian_rank'])==(34,35,5,7)
    two,F=specialize(m,dict(zip(one,vals)),True);fixed=two.index(local['chart_fixed_u']-1);cols=local['chart_columns'];q0=pt[two,:]
    J=s.Matrix(local['jacobian']);T=s.Matrix(local['tangent_basis']);n=T.cols;assert n==10 and T[-1,:]==s.zeros(1,10)
    _,piv=J.to_DM().rref();free=[i for i in range(45) if i not in piv];assert T[free,:]==s.eye(10)
    z,u,B=primitive_beta(m);B=B.subs(dict(zip(u,vals)))
    a,b,d,h,k=s.symbols('a b d h k');v1=s.Matrix([a,1,b,0]);v2=s.Matrix([b,0,d,1]);origin={a:0,b:0,d:0,h:0,k:0}
    sec=s.Matrix([(v1.T*B.subs(dict(zip(z,v1)))*v2)[0],(v1.T*B.subs(dict(zip(z,v2)))*v2)[0]])
    assert sec.subs(origin)==s.zeros(2,1)
    curve_jac=sec.jacobian([a,b,d]).subs(origin);assert curve_jac.rank()==2;direction=curve_jac.nullspace()[0]
    e1=x-a*y-b*w;e2=zz-b*y-d*w;fam=(e2**2+105*e1**2+h*e1*e2+k*e1**2)*(v1*v2.T-v2*v1.T);ds=[]
    for vv in [a,b,d,h,k]:
        dq=s.Matrix(coords(m,fam.diff(vv).subs(origin)))[two,:];dq-=q0*dq[fixed]/q0[fixed];ds.append(dq[cols,:].col_join(s.zeros(1,1)))
    family=s.Matrix.hstack(s.Matrix.hstack(*ds[:3])*direction,ds[3],ds[4]);coordinates=family[free,:]
    assert T*coordinates==family and J*family==s.zeros(135,3)
    params=[2,5,9];assert coordinates[params,:].det()!=0
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
    assert len(C)==len(R)==35 and np.array_equal(j[R,:]@U%p,np.eye(35,dtype=np.int64))
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
    return dict(sample=sample,J=J,T=T,F=F,two=two,fixed=fixed,chart=cols,free=free,rows=R,columns=C,parameters=params,family=coordinates,curve_jacobian=curve_jac,local=local,q0=q0,A=A,tseries={d:W[d][-1,:].tolist() for d in A},monomials=ms)

def source_geometry():
    m=S61();one=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==1];two=[i for i,wg in enumerate(m.reduced['parameter_weights'])if wg==2]
    b=s.symbols('b0:9');e=b[0];beta=s.Matrix([0,0,0,0,0,b[1],119*e/6,b[2],b[3],b[7],b[4],b[5],b[6],e,b[8],-35*e/6])
    x,y,z,w=X[1:];vv=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,3]=z*z+105*x*x;Q-=Q.T;q0=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q0)if v)
    G=(-b[8]+7*b[7])/(21*e)*x*x*z-s.Rational(14,3)*x*x*w-s.Rational(1,3)*x*z*w
    P=s.zeros(4);P[0,1]=x*x/2;P[2,3]=x*x/2;P[1,2]=-s.diff(G,w)-x*z;P[1,3]=s.diff(G,z);P-=P.T
    delta=s.Matrix(coords(m,P))[two,:];delta-=q0*delta[fi]/q0[fi]
    assert delta[35]==s.Rational(1,2)
    _,F=specialize(m,dict(zip(one,beta)),True)
    for f in F:
        for a,bb in [(q0,q0),(q0,delta),(delta,delta)]:
            assert s.cancel(sum(c*(a[mon[0]]*bb[mon[1]]+bb[mon[0]]*a[mon[1]]) for mon,c in f.items()if len(mon)==2))==0
    orbit=[]
    for A in sp_basis():
        dq=s.Matrix(coords(m,sum(((A*vv)[i]*Q.diff(vv[i])for i in range(4)),s.zeros(4))-A*Q-Q*A.T))[two,:];dq-=q0*dq[fi]/q0[fi];orbit.append(dq)
    O=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*O.nullspace());assert O.rank()==5 and stab.shape==(10,5)
    ub,action=beta_actions(m);H=beta.jacobian(b);ann=s.Matrix.hstack(*H.T.nullspace()).T;act=action.subs(dict(zip(ub,beta)));N=(ann*act*stab).applyfunc(s.factor);assert N.to_DM().rank()==2
    sample=dict(zip(b,[6,2,3,5,7,11,13,3,5]));D=H.row_join(act);ds=D.subs(sample);pc=list(ds.to_DM().rref()[1]);assert len(pc)==16;det=ds[:,pc].det();assert det
    joint=D.col_join(s.zeros(45,9).row_join(O));assert joint.subs(sample).rank()==16
    # Entire source tail, with the square degeneration retained in the closure.
    old=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text());fam=next(f for f in old['families']if f['label']=='M1');tail=next(t for t in fam['rational_tails']if t['axis']==5 and t['k']=='0')
    K,L,U,V=s.symbols('K L U V');ns=s.symbols('n0:6');gg=s.symbols('g0:9');tau=gg[5]
    Gfull=sum(s.sympify(tail['coordinate_formulas'].get(str(g),str(g)))*mon for g,mon in zip(gg,[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]));Ginf=s.diff(Gfull,tau)
    pp=s.Poly(s.diff(Ginf,z),z);discriminant=s.factor(s.discriminant(pp.as_expr(),z)/x**2);assert discriminant==-s.Rational(4860,847)*(L-4*U)**2
    a=(17*L-90*U)/22;h=3*(L-4*U)/77
    T=s.Matrix([[1/h,0,0,0],[0,h,0,a],[-a/h,0,1,0],[0,0,0,1]]);assert (T*J0[1:,1:]*T.T-J0[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
    Qinf=s.zeros(4);Qinf[1,3]=s.diff(Ginf,z);Qinf-=Qinf.T
    pull=(T.inv()*Qinf.subs(dict(zip(vv,T*vv)),simultaneous=True)*T.inv().T).applyfunc(s.cancel);assert(pull-3*Q/h).applyfunc(s.cancel)==s.zeros(4)
    source=s.Matrix([0,0,0,0,0,ns[0],U,ns[1],ns[2],V,ns[3],ns[4],ns[5],1,K-5*V,5*U-L])
    zz,uu,B=primitive_beta(m);from weight1_structure import COORDS
    vec=s.Matrix([B[i,j].diff(zz[k])for(i,j),k in COORDS]);basis=vec.jacobian(uu);rows=list(basis.T.to_DM().rref()[1]);iv=basis[rows,:].inv()
    BB=B.subs(dict(zip(uu,source)));pullB=T.T*BB.subs(dict(zip(zz,T*s.Matrix(zz))),simultaneous=True)*T
    nv=(iv*s.Matrix([pullB[i,j].diff(zz[k])for(i,j),k in COORDS])[rows,:]).applyfunc(s.cancel)
    assert(ann*nv).applyfunc(s.cancel)==s.zeros(7,1)
    # The source section h=1,a=0 and a beta scale e realizes all H9.
    sol={U:s.Rational(119,6),L:s.Integer(105),K:(b[8]+5*b[7])/e,V:b[7]/e,**dict(zip(ns,[b[i]/e for i in range(1,7)]))}
    assert(e*nv.subs(sol)-beta).applyfunc(s.cancel)==s.zeros(16,1)
    return dict(beta=list(map(str,beta)),parameters=list(map(str,b)),open_condition='b0 != 0',actual_line_delta=list(map(str,delta)),all135_actual_line_identities=True,projective_orbit_rank=5,stabilizer=enc(stab),source_annihilator=enc(ann),normal_stabilizer_matrix=enc(N),normal_stabilizer_rank=2,beta_sample={str(k):str(v)for k,v in sample.items()},beta_minor_columns=pc,beta_minor=str(det),joint_rank=16,image_dimension=16,source_tail_normalization=enc(T),normalized_source_beta=list(map(str,nv)),source_recovery={str(k):str(v)for k,v in sol.items()},endpoint_discriminant=str(discriminant),covered_tails=['M1 g5 k=0','Weyl partner M1 g8 k=-1'],flag_boundary='L=4U lies in the proper flag locus of the same16D closure')

PRIMARY_PROGRAM='''
assert(dim(R/J)==3 and degree(R/J)==65);
B=saturate(J,ideal(v_9));D=saturate(J,ideal(v_7-105*v_9));
DD=ideal(v_1,v_2,v_4,v_6,v_8,v_9,v_0-(28/3)*v_7,v_3+(2/3)*v_7,v_5+(16/6615)*v_7);
assert(D==DD and intersect(B,D)==J and degree(R/B)==65);
C=ideal(v_8,v_7-105*v_9,v_6+(18344/357)*v_9,v_4+2*v_9,v_3+210*v_9,v_1+(180/17)*v_9,v_0-(45860/17)*v_9);
assert((gens B)%C==0);assert all(flatten entries gens C,g->g^4%B==0);
assert(dim(R/(J+ideal(v_2,v_5,v_9)))==0 and degree(R/(J+ideal(v_2,v_5,v_9)))==66);
S=ZZ/32003[w_0,w_1,w_3,w_4,w_6,w_7,w_8,w_2,w_5,MonomialOrder=>{GRevLex=>7,GRevLex=>2}];
ff=map(S,R,{w_0,w_1,w_2,w_3,w_4,w_5,w_6,w_7,w_8,1_S});GG=gens gb(ff B);cc={};
scan(flatten entries GG,g->(ee:=take((exponents leadTerm g)#0,7);mm:=product(7,j->S_j^(ee#j));cf:=(sum(select(terms g,tt->take((exponents tt)#0,7)==ee)))//mm;if first degree cf>0 then cc=append(cc,cf/(leadCoefficient cf))));
assert(cc=={});NI=ideal apply(flatten entries GG,g->product(7,j->S_j^((exponents leadTerm g)#0#j)));assert(dim(S/NI)==2 and degree(S/NI)==65);print("CERTIFICATE PASS");stdio<<flush;
'''

def run(m2_command=None):
    start=time.time();r=point_jets();geometry=source_geometry();m=S61();local=r['local'];assert(local['leading_jacobian_rank'],local['homogenized_jacobian_rank'],local['quadratic_span_rank'],local['common_hessian_rank'])==(34,35,5,7)
    def mod(v):q=s.Rational(v);assert int(q.q)%PRIME;return int(q.p)*pow(int(q.q),-1,PRIME)%PRIME
    def rank(M):
        A=np.array([[mod(v) for v in row] for row in s.Matrix(M).tolist()],dtype=np.int64);j=0
        for c in range(A.shape[1]):
            ix=np.flatnonzero(A[j:,c])
            if not len(ix):continue
            a=j+int(ix[0]);A[[j,a]]=A[[a,j]];A[j]=A[j]*pow(int(A[j,c]),-1,PRIME)%PRIME;inds=np.flatnonzero(A[:,c]);inds=inds[inds!=j];A[inds]=(A[inds]-A[inds,c,None]*A[j,None,:])%PRIME;j+=1
            if j==len(A):break
        return j
    assert rank(r['J'][:,:44])==34 and rank(r['J'])==35 and rank(r['curve_jacobian'])==2 and rank(s.Matrix(local['quadratic_obstruction_basis']))==5
    hs=[]
    for row in s.Matrix(local['quadratic_obstruction_basis']).tolist():
        H=s.zeros(10)
        for c,(i,j) in zip(row,local['quadratic_monomials']):H[i,j]+=s.Rational(c);H[j,i]+=s.Rational(c)
        hs+=H.tolist()
    assert s.Matrix(hs)*r['family']==s.zeros(50,3)
    assert rank(s.Matrix(hs))==7 and mod(r['family'][r['parameters'],:].det())
    delta=s.Matrix(geometry['actual_line_delta']).applyfunc(s.sympify).subs({s.Symbol(k):s.Rational(v) for k,v in geometry['beta_sample'].items()});dt=delta[r['chart'],:].col_join(s.zeros(1,1));dc=dt[r['free'],:];assert r['T']*dc==dt
    assert list(dc)==[s.Rational(14,3),0,0,-s.Rational(1,3),0,-s.Rational(8,6615),0,s.Rational(1,2),0,0]
    with tempfile.TemporaryDirectory(prefix='poisson-pure-binary-') as tmp:
        tmp=Path(tmp);inp=tmp/'input.txt';basis=tmp/'basis.txt';binary=tmp/'closure';info=write_input(inp,PRIME,10,{d:A.tolist() for d,A in r['A'].items()},r['tseries'],maximum_degree=5)
        cc=subprocess.run(['c++','-O3','-std=c++17','-Wno-unqualified-std-cast-call',str(ROOT/'code/truncated_local_ideal_20261101.cpp'),'-o',str(binary)],capture_output=True,text=True,timeout=30);assert cc.returncode==0,cc.stderr
        cr=subprocess.run([str(binary),str(inp),str(basis)],capture_output=True,text=True,timeout=60);assert cr.returncode==0 and 'rank 1900 length 1103' in cr.stdout,cr.stdout+cr.stderr
        allinitial=read_initial(basis,10,5);assert len(allinitial)==1900;selected=[0,1,2,3,4,5,6,8,9,11,12,13,17,18,20,21,22,23,28,33,847,848,850];chosen=[allinitial[i] for i in selected];assert len(chosen)==23
        expr=lambda f:'+'.join(str(c)+'*'+'*'.join('v_'+str(i) for i in mon) for mon,c in f.items())
        program='R=ZZ/32003[v_0..v_9];J=ideal('+','.join(expr(f) for f in chosen)+');\n'+PRIMARY_PROGRAM
        source=tmp/'certificate.m2';source.write_text(program);env=os.environ.copy();local_binary=Path('/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary')
        if m2_command:cmd=shlex.split(m2_command)
        elif local_binary.exists():cmd=[str(local_binary)];env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        else:cmd=[shutil.which('M2') or 'M2']
        mr=subprocess.run(cmd+['--script',str(source)],env=env,capture_output=True,text=True,timeout=60);output=mr.stdout+mr.stderr;assert mr.returncode==0 and 'CERTIFICATE PASS' in output and 'error:' not in output,output
    prime_path=ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json';prime=json.loads(prime_path.read_text());assert prime['all_passed'] and prime['prime']==PRIME
    dependencies=['code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','code/s61_minimal_baseline_thickenings_20261101.py','code/s61_minimal_rankone_A1_removal_20261101.py','code/s61_minimal_infinity_local_data_20261101.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py','code/truncated_local_ideal_20261101.cpp','code/truncated_local_ideal_io_20261101.py','reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json','reports/S61_BINARY_UNIFORM_CUBIC_20261101.json','reports/S61_BINARY_UNIFORM_QUADRATIC_20261101.json','reports/S61_BINARY_UNIFORM_RANK_20261101.json','reports/S61_BINARY_PRIME_SPREADING_20261101.json','reports/S61_UNIVERSAL_WEIGHT6_20261031.json','research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md','research/notes/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.md']
    paths=[Path(__file__)]+[ROOT/f for f in dependencies]
    return dict(id='S61_MINIMAL_PURE_BINARY_ENDPOINT_20261101',all_passed=True,prime=PRIME,sample=r['sample'],source_geometry=geometry,point_local_data=local,curve_jacobian=enc(r['curve_jacobian']),actual_binary_family_tangent=enc(r['family']),binary_parameters=r['parameters'],actual_line_tangent_coordinates=list(map(str,dc)),jacobian_pivot_rows=r['rows'],jacobian_pivot_columns=r['columns'],truncated_input=info,truncated_closure_output=cr.stdout.splitlines()[:7],initial_generators=[[[list(mon),int(c)] for mon,c in f.items()] for f in chosen],initial_generator_degrees=[len(next(iter(f))) for f in chosen],macaulay2_program=program,macaulay2_output=output,initial_decomposition='B65 intersect actual reduced line',normal_component_primary=True,transverse_initial_length=66,conclusion='The specified characteristic-zero ambient local ring has t^2=0; the whole sufficiently general-beta fiber of the dominant16D pure-cubic binary endpoint is removed.',proof_note='research/notes/S61_MINIMAL_PURE_BINARY_ENDPOINT_20261101.md',limitations=['Uses the good-reduction and DVR torsion-free argument in the proof note, not a finite-field sample inference.','The flag endpoints of other pure-cubic tails and the remaining binary/flag strata remain separate.','No full affine-fiber finiteness or five-point exhaustion.'],input_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},seconds=time.time()-start)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--m2-command');args=parser.parse_args();report=run(args.m2_command)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
