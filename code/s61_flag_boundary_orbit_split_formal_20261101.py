#!/usr/bin/env python3
"""A fixed-beta flag-boundary germ: two group directions and nonzero square-zero t."""
import argparse,hashlib,json,os,subprocess,tempfile,time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J4full
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time()
    input_path=ROOT/'reports/S61_FLAG_BOUNDARY_RATIONAL_INITIAL_20261101.json'
    saved=json.loads(input_path.read_text());assert saved['all_passed']
    for name,digest in saved['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    model=S61();point=s.Matrix(saved['point61']).applyfunc(s.Rational)
    one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    beta_values=point[one,:]
    local=calculate(model,True,point)
    J=s.Matrix(local['jacobian']).applyfunc(s.Rational)
    T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
    q=point[two,:];C=list(J.to_DM().rref()[1]);free=[i for i in range(45) if i not in C]
    assert len(C)==36 and T[free,:]==s.eye(9)
    z,bu,beta=primitive_beta(model)
    bf=beta.subs(dict(zip(bu,beta_values)),simultaneous=True)
    a,lam=s.symbols('a lam')
    g=s.eye(4);g[0,1]=a;g[2,2]=lam;g[2,3]=lam-1/lam;g[3,3]=1/lam
    assert (g.T*J4full[1:,1:]*g-J4full[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
    transported=g.T*bf.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g
    assert (transported-bf).applyfunc(s.cancel)==s.zeros(4)
    A0=g.diff(a).subs({a:0,lam:1});A1=g.diff(lam).subs({a:0,lam:1})
    x,y,zz,w=X[1:];V=s.Matrix([x,y,zz,w]);Q=s.zeros(4);Q[1,3]=-x*x;Q-=Q.T
    def coords(form):
        P=X[0]**2*J4full;P[1:,1:]+=form
        return s.Matrix([model.normalize(P)['parameters'][i] for i in two])
    O=[]
    for A in [A0,A1]:
        dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T
        dq=coords(dQ);dq-=q*dq[fixed]/q[fixed]
        O.append(dq[cols,:].col_join(s.zeros(1,1)))
    O=s.Matrix.hstack(*O);OC=O[free,:]
    assert T*OC==O and OC[[0,1],:].det()==4
    assert OC[2:,:]==s.zeros(7,2)
    # The fixed-beta group acts transversely to the original-coordinate slice
    # v0=v1=0. Its finite formula, not merely its Lie algebra, fixes beta.
    assert saved['geometry']['b1_derivative'][0]==['0','0','-1/2','0','0','0','0','0','0']
    _,F=specialize(model,dict(zip(one,beta_values)),True)
    terms=[]
    for f in F:
        row=[]
        for mon,c in f.items():
            if fixed not in mon:
                i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon))
                row.append((i,j,c))
        terms.append(row)
    rows=list(J[:,C].T.to_DM().rref()[1]);U=s.zeros(45,36);iv=J[rows,C].inv(method='DM')
    for i,col in enumerate(C):U[col,:]=iv[i,:]
    assert J*T==s.zeros(135,9) and J[rows,:]*U==s.eye(36)
    factors=[s.ilcm(*[s.denom(v) for v in list(J[i,:])+[c for _,_,c in terms[i]]]) for i in range(135)]
    Jn=s.diag(*factors)*J;Un=U*s.diag(*[1/s.Integer(factors[i]) for i in rows])
    den=s.ilcm(*[s.denom(v) for v in Un]);Un=(den*Un).applyfunc(s.cancel)
    td=s.ilcm(*[s.denom(v) for v in T]);Tn=td*T
    # Row136 is the observable t, not an equation imposed on the scheme.
    Jn=Jn.col_join(s.zeros(1,44).row_join(s.ones(1,1)));terms.append([]);factors.append(1)
    ms={d:list(combinations_with_replacement(range(9),d)) for d in range(5)}
    ix={d:{m:i for i,m in enumerate(mm)} for d,mm in ms.items()}
    cpp=ROOT/'code/integer_implicit_jets_20261101.cpp'
    with tempfile.TemporaryDirectory(prefix='poisson-flag-orbit-split-') as temp:
        temp=Path(temp);inp=temp/'input.txt';coeff=temp/'coeff.txt';binary=temp/'jets'
        with inp.open('w') as out:
            out.write(f'136 9 4 36\n{den} {td}\n')
            out.write(' '.join(str(len(ms[d])) for d in range(5))+'\n')
            out.write(' '.join(map(str,rows))+'\n')
            for M in [Jn,Tn,Un]:
                for row in M.tolist():out.write(' '.join(map(str,row))+'\n')
            for i,ht in enumerate(terms):out.write(str(len(ht))+' '+' '.join(f'{a} {b} {c*factors[i]}' for a,b,c in ht)+'\n')
            for a in range(1,4):
                for b in range(1,5-a):out.write(' '.join(str(ix[a+b][tuple(sorted(x+y))]) for x in ms[a] for y in ms[b])+'\n')
        subprocess.run(['c++','-O3','-std=c++17','-I/opt/homebrew/include',str(cpp),'-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(binary)],check=True,capture_output=True,text=True,timeout=60)
        run=subprocess.run([str(binary),str(inp),str(coeff)],check=True,capture_output=True,text=True,timeout=60)
        assert 'EXACT INTEGER JETS PASS' in run.stderr
        dense=[line.split() for line in coeff.read_text().splitlines()]
        assert len(dense)==136 and all(len(row)==705 for row in dense)
        def expression(monco,prefix):
            return '+'.join('('+c+')*'+'*'.join(prefix+'_'+str(j) for j in mon) for mon,c in monco) or '0'
        initial=[expression([(tuple(mon),c) for mon,c in row],'v') for row in saved['initial_relations']]
        program='R=QQ[v_0..v_8];J=ideal('+','.join(initial)+');\n'
        program+='''K=J+ideal(v_0,v_1);B=saturate(K,ideal v_2);E=B/K;
assert(dim(R/K)==1 and degree(R/K)==5);
assert(dim(R/B)==1 and degree(R/B)==5 and degree(B+ideal v_2)==5);
assert(radical B==ideal(v_0,v_1,v_8,v_7,v_5,v_4,v_3+2*v_6,v_2+4*v_6));
assert(dim E==0 and degree E==111);
h=hilbertSeries(E,Reduce=>true);z=(ring numerator h)_0;
assert(numerator h==2*z+14*z^2+29*z^3+38*z^4+18*z^5+10*z^6 and value denominator h==1);
print("SLICE TORSION DEGREE6 PASS");stdio<<flush;
'''
        monomials=sum((ms[d] for d in range(2,5)),[])
        polys=[expression([(mon,c) for mon,c in zip(monomials,row) if c!='0/1' and all(j>1 for j in mon)],'w') for row in dense]
        program+='S=QQ[w_2..w_8];m=ideal vars S;I=ideal('+','.join(polys[:135])+');tjet='+polys[135]+';\n'
        program+='assert(tjet % (I+m^4)==0);assert(tjet % (I+m^5)!=0);print("NONZERO T ORDER4 PASS");\n'
        script=temp/'verify.m2';script.write_text(program)
        env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        test=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(script)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=120)
        assert test.returncode==0 and 'NONZERO T ORDER4 PASS' in test.stdout,test.stdout+test.stderr
        evidence=dict(native_log=run.stderr,m2_log=test.stdout+test.stderr,observable_jet_sha256=hashlib.sha256(coeff.read_bytes()).hexdigest())
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    paths=[Path(__file__),cpp,input_path,ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101.md',ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101',all_passed=True,field='QQ',point61=list(map(str,point)),
                fixed_beta_group=enc(g),group_parameters=['a','lam'],group_identity={'a':0,'lam':1},group_preserves_beta=True,
                orbit_tangent=enc(OC),orbit_slice_minor=4,orbit_slice_variables=[0,1],one_dimensional_base_variable=2,
                initial_slice_dimension=1,initial_slice_degree=5,primary_slice_free_rank=5,
                initial_torsion_length=111,initial_torsion_top_degree=6,homogenizer_order=4,
                homogenizer_nonzero=True,homogenizer_square_zero=True,
                proof_note='research/notes/S61_FLAG_BOUNDARY_ORBIT_SPLIT_FORMAL_20261101.md',
                limitations=['This is the fixed-beta germ at the stated sparse point.','No relative ambient t-nilpotence or entire18D boundary removal.','No full initial-ideal equality and no free49 assertion.'],
                **evidence,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
