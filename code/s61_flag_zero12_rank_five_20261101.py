#!/usr/bin/env python3
"""A uniform actual flag support and a rigorous mixed-characteristic rank5 sample."""
import argparse,hashlib,json,time,tempfile,subprocess,os
from pathlib import Path
from itertools import combinations_with_replacement
import numpy as np
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from truncated_local_ideal_io_20261101 import write_input,read_initial
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]


def sample_jets():
    start=time.time();p=32003;Dmax=6;r=json.load(open(ROOT/'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json'));old=json.load(open(ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'));m=S61();b=s.symbols('b0:15');env={s.Symbol(k):s.sympify(v)for k,v in r['substitutions'].items()};env.update(dict(zip(map(s.Symbol,r['parameters']),[1,2,3,5,7,11,13])));H0=s.Matrix(old['flag']['beta_kernel']).applyfunc(s.Rational);beta=(H0*s.Matrix(b)).subs(env);one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2];fixed=two.index(37);cols=[i for i in range(45)if i!=fixed];q0=s.Matrix(old['flag']['point61']).applyfunc(s.Rational)[two,:];_,F=specialize(m,dict(zip(one,beta)),True);J=s.zeros(135,45)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if len(mon)==2:
       i,j=mon
       if i!=fixed:J[row,cols.index(i)]+=c*q0[j]
       if j!=fixed:J[row,cols.index(j)]+=c*q0[i]
      elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
    T=s.Matrix(r['tangent']).applyfunc(s.sympify).subs(env);assert J*T==s.zeros(135,9);R,C=r['pivot_rows'],r['pivot_columns'];U=s.zeros(45,36);iv=J[R,C].inv()
    for i,c in enumerate(C):U[c,:]=iv[i,:]
    O=s.Matrix(old['flag']['q_orbit_differential']).applyfunc(s.Rational);allowed=[i for i in range(10)if i!=5];OC=O[:,allowed];cp=list(OC.rref()[1]);OC=OC[:,cp].col_join(s.zeros(1,len(cp)));CT=OC[r['free_coordinates'],:];assert T*CT==OC;params=list(CT.T.rref()[1]);print('actual C3 tangent',CT.tolist(),'parameter rows',params,'seconds',time.time()-start,flush=True)
    def mod(v):v=s.Rational(v);return int(v.p)*pow(int(v.q),-1,p)%p
    def arr(M):return np.array([[mod(v)for v in row]for row in M.tolist()],dtype=np.int64)
    j,uu=arr(J),arr(U);JU=j@uu%p;assert np.array_equal(j[R,:]@uu%p,np.eye(36,dtype=np.int64));HH=np.zeros((135,45,45),dtype=np.int64)
    for row,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      i,l=[cols.index(k)for k in mon]+[44]*(2-len(mon));v=mod(c)
      if i==l:HH[row,i,l]+=v
      else:HH[row,i,l]+=v*pow(2,-1,p)%p;HH[row,l,i]+=v*pow(2,-1,p)%p
    HH%=p;n=9;W={1:arr(T)};A={};ms={d:list(combinations_with_replacement(range(n),d))for d in range(1,Dmax+1)};idx={d:{m:i for i,m in enumerate(ms[d])}for d in ms}
    for degree in range(2,Dmax+1):
     Q=np.zeros((135,len(ms[degree])),dtype=np.int64)
     for da in range(1,degree//2+1):
      db=degree-da;fac=1 if da==db else 2;indices=np.array([idx[degree][tuple(sorted(a+b))]for a in ms[da]for b in ms[db]])
      for row in range(135):
       cross=W[da].T@(HH[row]@W[db]%p)%p;acc=np.bincount(indices,weights=cross.ravel(),minlength=len(ms[degree])).astype(np.int64);Q[row]=(Q[row]+fac*acc)%p
     W[degree]=-uu@Q[R,:]%p;A[degree]=(Q-JU@Q[R,:])%p
     print('degree',degree,'seconds',time.time()-start,flush=True)
    cache=dict(sample={str(k):str(v)for k,v in env.items()},beta=list(map(str,beta)),J=J,T=T,C_tangent=CT,parameters=params,A=A,tseries={d:W[d][-1,:].tolist()for d in A},monomials=ms)
    return cache


def support_geometry():
    model=S61();old=json.loads((ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text())
    low=json.loads((ROOT/'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json').read_text())
    H=s.Matrix(old['flag']['beta_kernel']).applyfunc(s.Rational);b=s.symbols('b0:15')
    sub={s.Symbol(k):s.sympify(v)for k,v in low['substitutions'].items()}
    z,u,B=primitive_beta(model);vec=s.Matrix([B[i,j].diff(z[k])for(i,j),k in COORDS])
    basis=vec.jacobian(u);rows=list(basis.T.rref()[1]);iv=basis[rows,:].inv()
    action=s.Matrix.hstack(*[iv*infinitesimal_action(vec,A)[rows,:]for A in sp_basis()])
    action=action.subs(dict(zip(u,H*s.Matrix(b)))).subs(sub).applyfunc(s.expand)
    assert action[0,:]==s.Matrix([[0,0,0,0,0,-1,0,0,0,0]])
    O=s.Matrix(old['flag']['q_orbit_differential']).applyfunc(s.Rational).col_join(s.zeros(1,10))
    allowed=[i for i in range(10)if i!=5];OC=O[:,allowed];cs=list(OC.rref()[1])
    T=s.Matrix(low['tangent']).applyfunc(s.sympify);CT=OC[low['free_coordinates'],cs]
    assert (T*CT-OC[:,cs]).applyfunc(s.expand)==s.zeros(45,3)
    parameters=list(CT.T.rref()[1]);assert parameters==[0,1,3]and CT[parameters,:].det()==4
    values=-action[[2],[allowed[i]for i in cs]];L=s.zeros(1,9);lv=values*CT[parameters,:].inv()
    for i,j in enumerate(parameters):L[0,j]=lv[i]
    assert (L*OC[low['free_coordinates'],:]+action[[2],allowed]).applyfunc(s.expand)==s.zeros(1,9)
    assert list(L)==[0,-b[8],0,s.Rational(1,2),0,0,0,0,0]
    # This verifies descent through the fixed-flag stabilizer as well.
    assert (action[[2],allowed]*s.Matrix.hstack(*OC.nullspace())).applyfunc(s.expand)==s.zeros(1,6)
    return dict(beta_normal=list(map(str,action[0,:])),actual_support_dimension=3,
                tangent=[[str(v)for v in row]for row in CT.tolist()],parameter_rows=parameters,
                parameter_minor='4',b1_derivative=list(map(str,L)),same_support_generic_length_lower_bound=5,
                scope='Uniform on all seven parameters; uses the already-proved ordinary-flag free5 theorem.')


def run():
    start=time.time();geometry=support_geometry();data=sample_jets();p=32003;A=data['A'];n=9;D=6
    source=ROOT/'reports/S61_FLAG_ZERO12_QUARTIC_BOUNDARY_20261101.json'
    previous=json.loads(source.read_text());assert previous['all_passed']
    for path,h in previous['input_sha256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
    with tempfile.TemporaryDirectory(prefix='poisson-zero12-rank5-')as directory:
        directory=Path(directory);binary=directory/'closure';select=directory/'select'
        for stem,out in [('truncated_local_ideal',binary),('truncated_local_ideal_minimal_generators',select)]:
            cc=subprocess.run(['c++','-O3','-std=c++17',str(ROOT/'code'/f'{stem}_20261101.cpp'),'-o',str(out)],capture_output=True,text=True,timeout=30)
            assert cc.returncode==0,cc.stderr
        inp=directory/'full.txt';initial=directory/'initial.txt';chosen=directory/'chosen.txt'
        write_input(inp,p,n,{d:a.tolist()for d,a in A.items()},data['tseries'],maximum_degree=D)
        full=subprocess.run([str(binary),str(inp),str(initial)],capture_output=True,text=True,timeout=120)
        assert full.returncode==0 and 'rank 3858 length 1147' in full.stdout,full.stdout+full.stderr
        selected=subprocess.run([str(select),str(inp),str(initial),str(chosen)],capture_output=True,text=True,timeout=30)
        assert selected.returncode==0,selected.stderr
        indices=list(map(int,chosen.read_text().split()));assert len(indices)==64
        all_initial=read_initial(initial,n,D);polys=[all_initial[i]for i in indices]
        counts={d:sum(len(next(iter(f)))==d for f in polys)for d in range(2,7)}
        assert counts=={2:8,3:7,4:9,5:25,6:15}
        expr=lambda f:'+'.join(str(c)+'*'+'*'.join('v_'+str(j)for j in mon)for mon,c in f.items())or'0'
        program=f'R=ZZ/{p}[v_0..v_8];J=ideal('+','.join(expr(f)for f in polys)+');\n'
        program+='C=ideal(v_8,v_7,v_5,v_4,v_3+2*v_6,v_2);\n'
        program+='assert(dim(R/J)==3 and degree(R/J)==5 and radical J==C);\n'
        program+='assert(dim(R/(J+ideal(v_0,v_1,v_3)))==0);\n'
        program+='print("INITIAL DIMENSION3 DEGREE5 FINITE BASE PASS");print(hilbertSeries(R/J,Reduce=>true));\n'
        m2=directory/'verify.m2';m2.write_text(program);env=os.environ.copy()
        env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        mr=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(m2)],env=env,capture_output=True,text=True,timeout=60)
        assert mr.returncode==0 and 'INITIAL DIMENSION3 DEGREE5 FINITE BASE PASS' in mr.stdout,mr.stdout+mr.stderr
        normal=[i for i in range(9)if i not in data['parameters']];normal_jets={};normal_t={}
        for d,AA in A.items():
            keep=[i for i,mon in enumerate(data['monomials'][d])if all(j in normal for j in mon)]
            normal_jets[d]=AA[:,keep].tolist();normal_t[d]=[data['tseries'][d][i]for i in keep]
        transverse=directory/'transverse.txt';write_input(transverse,p,6,normal_jets,normal_t,maximum_degree=6)
        tr=subprocess.run([str(binary),str(transverse)],capture_output=True,text=True,timeout=30)
        assert tr.returncode==0 and 'modm 6 length 49' in tr.stdout and 'modm 7 length 49' in tr.stdout,tr.stdout+tr.stderr
        logs=dict(full_closure='\n'.join(line for line in full.stdout.splitlines()if not line.startswith('t ')),
                  selection=selected.stderr,m2=mr.stdout+mr.stderr,transverse=tr.stdout)
    paths=[Path(__file__),source,ROOT/'reports/S61_FLAG_ZERO12_LOW_JETS_20261101.json',ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',
           ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',
           ROOT/'code/truncated_local_ideal_io_20261101.py',ROOT/'code/truncated_local_ideal_20261101.cpp',
           ROOT/'code/truncated_local_ideal_minimal_generators_20261101.cpp',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
           ROOT/'code/weight1_structure.py',ROOT/'research/notes/ADIC_MODULAR_LENGTH_BOUND_20261101.md']
    return dict(id='S61_FLAG_ZERO12_RANK_FIVE_20261101',all_passed=True,uniform_support=geometry,
                sample=data['sample'],beta=data['beta'],guide_prime=p,transverse_mod_p_length=49,
                characteristic_zero_transverse_length_upper_bound=49,
                initial_relation_degrees=counts,initial_relations_mod_p=[[[list(mon),c]for mon,c in f.items()]for f in polys],
                special_initial_dimension=3,special_initial_degree=5,parameter_rows=[0,1,3],
                characteristic_zero_finite_base_rank=5,t_is_base_torsion=True,t_nilpotence_proved=False,
                full_characteristic_zero_initial_ideal_identified=False,
                limitations=['Finite R3 rank5 is proved only at the specified rational beta point.',
                             'The actual smooth C3 and its generic free5 lower bound are uniform on the entire seven-parameter chart.',
                             'No equality of QQ and finite-field initial ideals, no exact QQ transverse length49, and no t-nilpotence claim.',
                             'No entire ZERO12 or flag-boundary saturation removal follows.'],
                proof_note='research/notes/S61_FLAG_ZERO12_RANK_FIVE_20261101.md',logs=logs,
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths},seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=run()
    with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'])))
