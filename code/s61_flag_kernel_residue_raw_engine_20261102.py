"""Reconstruct exact seven-variable six-jets after setting v0=v1=0."""
import subprocess
from fractions import Fraction
from itertools import combinations_with_replacement as cwr
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]

def raw_residual(initial,directory):
    point=s.Matrix(initial['point61']).applyfunc(s.Rational);model=S61()
    one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    local=calculate(model,True,point)
    J=s.Matrix(local['jacobian']).applyfunc(s.Rational)
    T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational)
    _,equations=specialize(model,dict(zip(one,[point[i]for i in one])),True)
    chart=local['chart_columns'];fixed=next(i for i in range(45)if i not in chart)
    terms=[]
    for f in equations:
        row=[]
        for mon,c in f.items():
            if fixed not in mon:
                a,b=[chart.index(i)for i in mon]+[44]*(2-len(mon));row.append((a,b,c))
        terms.append(row)
    pivot=list(J.to_DM().rref()[1]);selected=list(J[:,pivot].T.to_DM().rref()[1])
    U=s.zeros(45,36);inverse=J[selected,pivot].inv(method='DM')
    for i,col in enumerate(pivot):U[col,:]=inverse[i,:]
    assert J*T==s.zeros(135,9)and J[selected,:]*U==s.eye(36)
    scales=[s.ilcm(*[s.denom(c)for c in list(J[i,:])+[c for a,b,c in terms[i]]])for i in range(135)]
    Jn=s.diag(*scales)*J;Un=U*s.diag(*[1/s.Integer(scales[i])for i in selected])
    denominator=s.ilcm(*[s.denom(c)for c in Un]);Un=(denominator*Un).applyfunc(s.cancel)
    tangent_denominator=s.ilcm(*[s.denom(c)for c in T]);Tn=tangent_denominator*T[:,list(range(2,9))]
    ms={d:list(cwr(range(7),d))for d in range(7)}
    index={d:{m:i for i,m in enumerate(mm)}for d,mm in ms.items()}
    source=ROOT/'code/monomial_quotient_implicit_jets_20261102.cpp';binary=directory/'raw-engine'
    subprocess.run(['c++','-O3','-w','-std=c++17','-I/opt/homebrew/include',str(source),
                    '-L/opt/homebrew/lib','-lgmpxx','-lgmp','-o',str(binary)],check=True)
    path=directory/'raw-input.txt'
    with path.open('w')as out:
        out.write(f'135 7 6 36\n{denominator} {tangent_denominator}\n')
        out.write(' '.join(str(len(ms[d]))for d in range(7))+'\n'+' '.join(map(str,selected))+'\n')
        for matrix in [Jn,Tn,Un]:
            for row in matrix.tolist():out.write(' '.join(map(str,row))+'\n')
        for i,row in enumerate(terms):out.write(str(len(row))+' '+' '.join(f'{a} {b} {c*scales[i]}'for a,b,c in row)+'\n')
        for a in range(1,6):
            for b in range(1,7-a):
                out.write(' '.join(str(index[a+b][tuple(sorted(x+y))])for x in ms[a]for y in ms[b])+'\n')
    output=directory/'raw-coefficients.txt'
    run=subprocess.run([str(binary),str(path),str(output)],capture_output=True,text=True,check=True,timeout=40)
    assert 'EXACT INTEGER JETS PASS' in run.stderr
    rows=[list(map(Fraction,line.split()))for line in output.read_text().splitlines()]
    assert len(rows)==136 and all(len(row)==1708 for row in rows)
    return rows,run.stderr
