#!/usr/bin/env python3
"""Exact positive-weight finite implicit model on the nuisance-zero deeper core."""
import argparse
from collections import Counter
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import numpy as np
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61,BASIS,X
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_deeper_uniform_quadratic_20261101 import build_quadratics,elimination_rank,determinant_mod
ROOT=Path(__file__).resolve().parents[1]
NOTE=ROOT/'research/notes/S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102.md'
NW=[1,1,1,1,1,2,1,2,1]

def monomials(weight):
    def rec(pos,left,mon):
        if pos==9:
            if left==0:yield tuple(mon)
            return
        for value in range(left//NW[pos]+1):yield from rec(pos+1,left-value*NW[pos],mon+[value])
    return list(rec(0,weight,[]))

def verify():
    started=time.time();model=S61()
    chart=ROOT/'reports/S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101.json'
    quad_path=ROOT/'reports/S61_BINARY_DEEPER_UNIFORM_QUADRATIC_20261101.json'
    data=json.loads(chart.read_text())['linear'];R,Q=build_quadratics(data)
    S,e,d=ring('e,d',s.QQ);z=S.zero
    core=lambda p:S.from_dict({m[:2]:c for m,c in p.items() if not any(m[2:])})
    pol=lambda v:core(R.from_expr(s.sympify(v)))
    T=[[core(v) for v in row] for row in Q['Tn']]
    U=[[core(v) for v in row] for row in Q['Un']]
    J=[[pol(v) for v in row] for row in data['J']]
    den=core(Q['Uden']);dens=[core(v) for v in Q['Tden']]
    known=[S.from_expr(s.sympify(v)) for v in json.loads(quad_path.read_text())['irreducible_minor_factors']]
    for value in [den,*dens]:
        assert value
        for factor in known:
            while value.rem(factor)==0:value=value.exquo(factor)
        assert max(map(sum,value))==0
    sel=data['pivot_rows']
    for i,row in enumerate(J):
        for k in range(9):assert sum((v*T[j][k] for j,v in enumerate(row) if v),z)==0
    for a,i in enumerate(sel):
        for b in range(33):assert sum((v*U[j][b] for j,v in enumerate(J[i]) if v),z)==(den if a==b else z)
    chart_basis=s.Matrix(data['chart_basis']).applyfunc(s.sympify)
    assert not chart_basis.free_symbols
    chart_inverse=chart_basis.inv()
    inverse=[[S.domain.convert(v) for v in row] for row in chart_inverse.tolist()]
    Tc=[[sum((a*T[j][k] for j,a in enumerate(row) if a),z) for k in range(9)] for row in inverse]
    Uc=[[sum((a*U[j][k] for j,a in enumerate(row) if a),z) for k in range(33)] for row in inverse]
    for row in data['free']:assert all(v==0 for v in Uc[row])
    for k,idx in enumerate(data['normal']):
        assert Tc[data['free'][idx]]==[dens[k] if j==k else z for j in range(9)]
    cuts=[i for i in range(12) if i not in data['normal']]
    assert all(not any(Tc[data['free'][i]]) for i in cuts)
    tw=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    torus=[0,1,-1,-1,1];raw=[]
    for idx in model.certificate['raw_free_indices']:
        pair,mon=BASIS[idx]
        raw.append(sum(torus[i] for i in pair)-sum(a*b for a,b in zip(torus,s.Poly(mon,*X).monoms()[0])))
    assert all((raw[tw[i]]+4)%2==0 for i in data['cols'])
    cw=[(raw[tw[i]]+4)//2 for i in data['cols']]+[3]
    assert min(cw)>=0 and max(cw)==4
    beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()}
    _,F=specialize(model,beta,True);H=[];rw=[]
    for i,equation in enumerate(F):
        terms=[];weights={cw[j] for j,v in enumerate(J[i]) if v}
        for mon,c in equation.items():
            if data['fixed'] in mon:continue
            ij=[data['cols'].index(j) for j in mon]+[44]*(2-len(mon));c=pol(c)
            if c:terms.append((ij,c));weights.add(sum(cw[j] for j in ij))
        assert len(weights)==1
        rw.append(next(iter(weights)));H.append(terms)
    assert all(not T[i][j] or cw[i]==NW[j] for i in range(45) for j in range(9))
    assert all(not U[i][j] or cw[i]==rw[sel[j]] for i in range(45) for j in range(33))
    assert max(rw[i] for i in sel)<=4
    native=['code/parametric_implicit_jets_20261101.cpp','code/s61_binary_deeper_weighted_extract_20261102.cpp']
    logs=[]
    with tempfile.TemporaryDirectory(prefix='s61-weighted-finite-') as temp:
        wd=Path(temp);inp=wd/'input';rawout=wd/'raw';rowsfile=wd/'rows';meta=wd/'meta';weightsfile=wd/'weights';coeff=wd/'coeff'
        def write_poly(stream,p,normals=0):
            stream.write(str(len(p))+'\n')
            for m,c in p.items():stream.write(str(c)+' '+' '.join(map(str,m+(0,)*normals))+'\n')
        with inp.open('w') as stream:
            stream.write('135 2 9 4 33\n'+' '.join(map(str,sel))+'\n');write_poly(stream,den,9)
            for M in [J,T,U]:
                for row in M:
                    for p in row:write_poly(stream,p,9)
            for terms in H:
                stream.write(str(len(terms))+'\n')
                for (i,j),p in terms:stream.write(f'{i} {j}\n');write_poly(stream,p,9)
        weightsfile.write_text(' '.join(map(str,cw+rw+NW))+'\n')
        with coeff.open('w') as stream:
            write_poly(stream,den);stream.write(str(len(known))+'\n')
            for p in known:write_poly(stream,p)
        binaries=[]
        env=dict(os.environ);env['DYLD_LIBRARY_PATH']=str(ROOT/'work/m2-runtime/lib')+':/opt/homebrew/lib'
        for k,path in enumerate(native):
            binary=wd/f'engine{k}';binaries.append(binary)
            subprocess.run(['c++','-O3','-std=c++17','-I'+str(ROOT/'work/flint-headers/flint/3.5.0/include'),'-I'+str(ROOT/'work/flint-headers/mpfr/4.2.2/include'),'-I/opt/homebrew/include',str(ROOT/path),str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'),'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
        for command in [[binaries[0],inp,rawout],[binaries[1],rawout,coeff,rowsfile,meta,weightsfile]]:
            result=subprocess.run(list(map(str,command)),check=True,capture_output=True,text=True,env=env);logs.append(result.stderr)
        contents=meta.read_text().splitlines()[1:]
        assert all(v.split()[3]=='1' and v.split()[-2:]==['0','0'] for v in contents)
        counts=Counter(int(v.split()[1]) for v in contents);assert counts=={2:7,3:21,4:33}
        prime=32003;point=(2,3);it=iter(rowsfile.read_text().split());nr=int(next(it));rows=[]
        for _ in range(nr):
            idx,w,n=[int(next(it)) for _ in range(3)];row={}
            for _ in range(n):
                mon=tuple(int(next(it)) for _ in range(9));value=0
                for _ in range(int(next(it))):
                    c=s.Rational(next(it));a,b=int(next(it)),int(next(it));assert int(c.q)%prime
                    value=(value+int(c.p)*pow(int(c.q),-1,prime)*pow(point[0],a,prime)*pow(point[1],b,prime))%prime
                if value:row[mon]=value
            rows.append((idx,w,row))
        ranks={}
        for degree in [2,3,4]:
            cols=monomials(degree);lookup={m:i for i,m in enumerate(cols)};matrix=[];sources=[]
            for idx,w,row in rows:
                if w>degree:continue
                for shift in monomials(degree-w):
                    v=np.zeros(len(cols),dtype=np.int64)
                    for m,c in row.items():v[lookup[tuple(a+b for a,b in zip(m,shift))]]=c
                    matrix.append(v);sources.append([idx,list(shift)])
            M=np.array(matrix);rank,mr=elimination_rank(M.T,prime);_,mc=elimination_rank(M[mr],prime)
            determinant=determinant_mod(M[mr][:,mc],prime);assert determinant
            ranks[degree]=dict(rows=len(M),columns=len(cols),rank=rank,quotient_dimension=len(cols)-rank,minor_rows=mr,minor_columns=mc,minor_sources=[sources[i] for i in mr],minor_mod_prime=determinant)
        assert [ranks[i]['rank'] for i in [2,3,4]]==[7,64,269]
        fourth=monomials(4)
        assert len(fourth)==269 and all(any(all(a>=b for a,b in zip(m,q)) for q in fourth) for m in monomials(5))
        hashes={'raw_input':hashlib.sha256(inp.read_bytes()).hexdigest(),'raw_output':hashlib.sha256(rawout.read_bytes()).hexdigest(),'weighted_rows':hashlib.sha256(rowsfile.read_bytes()).hexdigest()}
    paths=[Path(__file__),NOTE,chart,quad_path,*[ROOT/p for p in native],ROOT/'code/s61_binary_deeper_quartic_field_20261101.hpp',ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py']
    return dict(id='S61_BINARY_DEEPER_WEIGHTED_FINITE_MODEL_20261102',all_passed=True,coefficient_ring='QQ[e,d,D^-1], D the thirteen previously treated factors',nuisance_parameters_zero=True,normal_weights=NW,local_coordinate_weights=cw,equation_weights=rw,selected_equations=sel,implicit_linear_identities_checked=True,all_input_and_output_monomial_weights_checked=True,implicit_solution_is_polynomial_of_ordinary_degree_at_most4=True,weighted_residual_row_counts=dict(counts),all_row_contents_supported_on_13_treated_factors=True,weight_four_columns=dict(light_four=210,light_squared_heavy=56,heavy_squared=3,total=269),all_weight_five_monomials_divisible_by_weight_four=True,sample=dict(prime=prime,point=list(point),ranks=ranks,weighted_hilbert_function=[1,7,23,34]),native_logs=logs,generated_sha256=hashes,proof_note=str(NOTE.relative_to(ROOT)),limitations=['The finite model is on nuisance=0 only.','No all-core rank7/64/269 certificate is asserted.','No new core-curve removal, ambient nilpotence, or full-fiber exhaustion.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-started)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);args=parser.parse_args()
    report=verify()
    with args.output.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print(json.dumps({k:report[k] for k in ['all_passed','seconds']}))
