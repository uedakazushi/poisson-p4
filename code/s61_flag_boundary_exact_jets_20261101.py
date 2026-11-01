#!/usr/bin/env python3
"""Exact flag-boundary six-jets, with modular and independent QQ cross-checks."""
import argparse,hashlib,json,os,subprocess,tarfile,tempfile,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]

def verify():
    started=time.time();archive=ROOT/'reports/S61_FLAG_BOUNDARY_SEVENTH_EXPLORATION_20261101.tar.gz'
    receipt=ROOT/'reports/S61_FLAG_BOUNDARY_SEVENTH_EXPLORATION_20261101.json';saved=json.loads(receipt.read_text())
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==saved['archive_sha256']
    def read(name):
        with tarfile.open(archive,'r:gz') as f:data=f.extractfile('work/'+name).read()
        assert hashlib.sha256(data).hexdigest()==saved['archived_sha256']['work/'+name]
        return json.loads(data)
    local=read('r0-flag-boundary-points.json')['b1']['local'];modular=read('flag-b1-jets-6.json')
    model=S61();point=s.Matrix(local['point61']).applyfunc(s.Rational)
    beta={i:point[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    two,F=specialize(model,beta,True);cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols)
    q0=point[two,:];assert all(q0[i]==0 for i in cols)
    J=s.zeros(135,45);terms=[]
    for row,f in enumerate(F):
        ht=[]
        for mon,c in f.items():
            if len(mon)==2:
                a,b=mon
                if a!=fixed:J[row,cols.index(a)]+=c*q0[b]
                if b!=fixed:J[row,cols.index(b)]+=c*q0[a]
            elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
            if fixed not in mon:
                a,b=[cols.index(i) for i in mon]+[44]*(2-len(mon));ht.append((a,b,c))
        terms.append(ht)
    assert J==s.Matrix(local['jacobian']).applyfunc(s.Rational)
    T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);rows=modular['pivot_rows'];C=modular['pivot_columns']
    U=s.zeros(45,36);inverse=J[rows,C].inv(method='DM')
    for i,j in enumerate(C):U[j,:]=inverse[i,:]
    assert J*T==s.zeros(135,9) and T[modular['free'],:]==s.eye(9) and J[rows,:]*U==s.eye(36)
    factors=[s.ilcm(*[s.denom(v) for v in list(J[i,:])+[c for a,b,c in terms[i]]]) for i in range(135)]
    Jn=s.diag(*factors)*J;Un=U*s.diag(*[1/s.Integer(factors[i]) for i in rows]);den=s.ilcm(*[s.denom(v) for v in Un]);Un=(den*Un).applyfunc(s.cancel)
    td=s.ilcm(*[s.denom(v) for v in T]);Tn=td*T;assert den%32003 and td%32003
    ms={d:list(combinations_with_replacement(range(9),d)) for d in range(7)};ix={d:{m:i for i,m in enumerate(mm)} for d,mm in ms.items()}
    cpp=ROOT/'code/integer_implicit_jets_20261101.cpp'
    with tempfile.TemporaryDirectory(prefix='poisson-exact-jets-') as temporary:
        temp=Path(temporary);source=temp/'input.txt';coeff=temp/'coefficients.txt';binary=temp/'jets'
        with source.open('w') as out:
            out.write(f'135 9 6 36\n{den} {td}\n');out.write(' '.join(str(len(ms[d])) for d in range(7))+'\n');out.write(' '.join(map(str,rows))+'\n')
            for M in [Jn,Tn,Un]:
                for row in M.tolist():out.write(' '.join(map(str,row))+'\n')
            for i,ht in enumerate(terms):out.write(str(len(ht))+' '+' '.join(f'{a} {b} {c*factors[i]}' for a,b,c in ht)+'\n')
            for a in range(1,6):
                for b in range(1,7-a):out.write(' '.join(str(ix[a+b][tuple(sorted(x+y))]) for x in ms[a] for y in ms[b])+'\n')
        subprocess.run(['c++','-O3','-std=c++17','-I/opt/homebrew/include',str(cpp),'-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(binary)],check=True,capture_output=True,text=True,timeout=60)
        run=subprocess.run([str(binary),str(source),str(coeff)],check=True,capture_output=True,text=True,timeout=60)
        assert 'EXACT INTEGER JETS PASS' in run.stderr
        native=coeff.read_text().splitlines();assert len(native)==135;checked=0;prime=32003
        for i,line in enumerate(native):
            expected=[int(v)*int(factors[i])%prime for d in range(2,7) for v in modular['A'+str(d)][i]]
            actual=line.split();assert len(actual)==len(expected)
            for v,w in zip(actual,expected):
                a,b=map(int,v.split('/'));assert a*pow(b,-1,prime)%prime==w;checked+=1
        assert checked==674325
        # A separate polynomial-ring implementation checks every QQ coefficient through degree four.
        fmt=lambda M:'matrix{'+','.join('{'+','.join(str(c) for c in row)+'}' for row in M.tolist())+'}'
        program='R=ZZ[v_0..v_8];T='+fmt(Tn)+';J='+fmt(Jn)+';U='+fmt(Un)+';jetDen='+str(den)+';rows={'+','.join(map(str,rows))+'};\n'
        cross=[];square=[]
        for row,ht in enumerate(terms):
            cross.append('{'+('+'.join(f'({c*factors[row]})*(V_({i},0)*W_({j},0)+V_({j},0)*W_({i},0))' for i,j,c in ht) or '0_R')+'}')
            square.append('{'+('+'.join(f'({c*factors[row]})*V_({i},0)*V_({j},0)' for i,j,c in ht) or '0_R')+'}')
        program+='jetB=(V,W)->matrix{'+','.join(cross)+'};jetQ=V->matrix{'+','.join(square)+'};\n'
        program+='W1=T*transpose vars R;Q2=jetQ(W1);W2=-U*Q2^rows;A2=jetDen*Q2+J*W2;Q3=jetB(W1,W2);W3=-U*Q3^rows;A3=jetDen*Q3+J*W3;Q4=jetB(W1,W3)+jetQ(W2);W4=-U*Q4^rows;A4=jetDen*Q4+J*W4;S=QQ[v_0..v_8];R=S;use R;\n'
        for d in range(2,5):
            monomials=['*'.join('v_'+str(i) for i in mon) for mon in ms[d]]
            program+=f'A{d}=sub(A{d},S)/({den**(d-1)*td**d});m{d}=matrix{{{{'+','.join(monomials)+f'}}}};C{d}=sub(transpose last coefficients(transpose A{d},Monomials=>m{d}),QQ);\n'
        exact4=temp/'m2-coefficients.txt'
        program+='C=C2|C3|C4;ff=openOut '+json.dumps(str(exact4))+';scan(entries C,rr->(scan(rr,a->(ff<<numerator a<<"/"<<denominator a<<" "));ff<<endl));close ff;print("EXACT FOURJET PASS");\n'
        m2source=temp/'fourjet.m2';m2source.write_text(program);env=os.environ.copy();env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib:/opt/homebrew/opt/givaro/lib'
        m2=subprocess.run(['/opt/homebrew/Cellar/macaulay2/1.26.05/bin/M2-binary','--script',str(m2source)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=60)
        assert m2.returncode==0 and 'EXACT FOURJET PASS' in m2.stdout,m2.stdout+m2.stderr
        comparison=exact4.read_text().splitlines();assert len(comparison)==135
        for a,b in zip(native,comparison):assert a.split()[:705]==b.split()
        evidence=dict(integer_jet_log=run.stderr,exact_coefficient_file_sha256=hashlib.sha256(coeff.read_bytes()).hexdigest(),exact_input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),independent_m2_program_sha256=hashlib.sha256(program.encode()).hexdigest())
    paths=[Path(__file__),cpp,archive,receipt,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'research/notes/S61_FLAG_BOUNDARY_EXACT_JETS_20261101.md']
    return dict(id='S61_FLAG_BOUNDARY_EXACT_JETS_20261101',all_passed=True,field='QQ',normal_formal_variables=9,highest_verified_degree=6,
        raw_equations=135,modular_coefficient_comparisons=checked,independent_exact_QQ_coefficient_comparisons=95175,
        initial_ideal_identified=False,completed_homogenizer_nilpotence_proved=False,
        limitations=['These are exact formal equation jets, not a complete initial ideal or local-ring presentation.','The large rational initial-relation replay remains incomplete.'],
        **evidence,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-started)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
