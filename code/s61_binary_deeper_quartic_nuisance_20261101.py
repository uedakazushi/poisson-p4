#!/usr/bin/env python3
"""Exact quartic constraints from tracked quadratic-linear syzygies.

All native arithmetic is rational. No ignored work file is a mathematical
input. The conclusion concerns a nonempty open in the TWO core parameters;
it does not remove their exceptional divisors.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
from itertools import combinations_with_replacement

import numpy as np
import sympy as s
from sympy.polys.rings import ring

from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_binary_deeper_uniform_quadratic_20261101 import PIVOT_ROWS, PIVOT_COLS, elimination_rank, determinant_mod
from s61_binary_deeper_cubic_nuisance_20261101 import prepare as prepare_cubic, CERTIFICATE, LINEAR, QUADRATIC

ROOT = Path(__file__).resolve().parents[1]
SELECTED = [53,56,62,209,119,215,125,210,120,207,83,92,89,98,84,93,208,213,118,123,96,174,156,175,147,173,179,192]
SUPPORT = [0,1,4,7,8,9,12,18,21,25,26,27,30,32,34,35,40]
MONS = {d:list(combinations_with_replacement(range(9),d)) for d in (2,3,4)}
SOURCES = ['core_inverse','nuisance_dag','unit_relations','constant_span','products','nuisance_check']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def core_write(stream, poly):
    stream.write(str(len(poly))+' '+' '.join(str(c)+' '+' '.join(map(str,m)) for m,c in poly.items())+'\n')


def normal_exponent(indices):
    answer=[0]*9
    for i in indices: answer[i]+=1
    return answer


def write16(stream, terms):
    stream.write(str(len(terms))+'\n')
    for m,c in terms:
        assert len(m)==16
        stream.write(str(c)+' '+' '.join(map(str,m))+'\n')


def prepare_quartic(data, certificate, quad, tmp):
    R,*_ = ring('e,d,p,q,r,h,k',s.QQ)
    K=s.QQ.frac_field('e','d'); R2=K.field.ring
    R5,*_ = ring('p,q,r,h,k',K)
    basis=[[K.from_sympy(s.sympify(v)) for v in row]
           for row in certificate['quadratic_rref']['rref']]
    def outer(poly):
        grouped={}
        for m,c in poly.items(): grouped.setdefault(m[2:],{})[m[:2]]=c
        return R5.from_dict({m:K.field.raw_new(R2.from_dict(c),R2.one) for m,c in grouped.items()})
    B=[[outer(quad['Ag'][i][j]) for j in PIVOT_COLS] for i in PIVOT_ROWS]
    with (tmp/'core.txt').open('w') as f:
        f.write('24 24\n')
        for row in B:
            entries=[(j,v.get((0,)*5,K.zero)) for j,v in enumerate(row)]
            entries=[(j,v) for j,v in entries if v]
            f.write(str(len(entries))+'\n')
            for j,v in entries:
                f.write(str(j)+'\n'); core_write(f,v.numer); core_write(f,v.denom)
    with (tmp/'outer.txt').open('w') as f:
        f.write('24\n')
        for row in B:
            for v in row:
                f.write(str(len(v))+'\n')
                for m,c in v.items():
                    f.write(' '.join(map(str,m))+'\n'); core_write(f,c.numer); core_write(f,c.denom)
    ix4={m:i for i,m in enumerate(MONS[4])}
    with (tmp/'module4.txt').open('w') as f:
        f.write('1080 495\n')
        for row in basis:
            for a in MONS[2]:
                entries={ix4[tuple(sorted(a+m))]:v for m,v in zip(MONS[2],row) if v}
                f.write(str(len(entries))+'\n')
                for j,v in entries.items():
                    f.write(str(j)+'\n'); core_write(f,v.numer); core_write(f,v.denom)
    # The change of the quadratic implicit lift is in I2. Only quartics
    # induced by linear syzygies are used, so this changes them in I2*S2.
    reduced=[]
    for original in quad['W']:
        previous=[outer(v) for v in original]; row=previous[:]
        for k,pivot in enumerate(PIVOT_COLS):
            if previous[pivot]: row=[a-previous[pivot].mul_ground(b) for a,b in zip(row,basis[k])]
        assert all(not row[j] for j in PIVOT_COLS)
        reduced.append(row)
    assert not any(reduced[44]) and all(not v for v in quad['Tn'][44])
    assert all(not any(row) for i,row in enumerate(reduced) if i not in SUPPORT)
    Dw=R2.one
    for i in SUPPORT:
        for poly in reduced[i]:
            for c in poly.values(): Dw=Dw.lcm(c.denom)
    def clear(poly):
        ans={}
        for nu,c in poly.items():
            for m,v in (c.numer*Dw.exquo(c.denom)).items(): ans[m+nu]=v
        return R.from_dict(ans)
    lift=[[clear(v) for v in reduced[i]] for i in SUPPORT]
    pol=lambda x:R.from_expr(s.sympify(x)); zero=R.zero
    model=S61(); beta={int(k[1:])-1:s.sympify(v) for k,v in data['beta'].items()}
    _,equations=specialize(model,beta,True)
    T,U,Du=quad['Tn'],quad['Un'],quad['Uden']
    J=[[pol(v) for v in row] for row in data['J']]
    operators=[]
    for equation in equations:
        mat=[[zero]*len(SUPPORT) for _ in range(9)]
        for mon,c in equation.items():
            if len(mon)!=2 or data['fixed'] in mon: continue
            i,j=[data['cols'].index(v) for v in mon]; c=pol(c)
            for a in range(9):
                if i in SUPPORT and T[j][a]: mat[a][SUPPORT.index(i)]+=c*T[j][a]
                if j in SUPPORT and T[i][a]: mat[a][SUPPORT.index(j)]+=c*T[i][a]
        operators.append(mat)
    projection=[[sum((c*U[j][k] for j,c in enumerate(row) if c),zero) for k in range(33)] for row in J]
    selected=data['pivot_rows']
    with (tmp/'operators.txt').open('w') as f:
        f.write('24 17\n')
        for i in PIVOT_ROWS:
            assert all(not any(m[2:]) for m in quad['contents'][i])
            core_write(f,R2.from_expr(quad['contents'][i].as_expr()))
        for i in PIVOT_ROWS:
            for b in range(17):
                terms=[]
                for a in range(9):
                    value=Du*operators[i][a][b]-sum((c*operators[selected[k]][a][b] for k,c in enumerate(projection[i]) if c),zero)
                    terms.extend((list(m)+normal_exponent((a,)),c) for m,c in value.items())
                write16(f,terms)
        for row in lift:
            write16(f,[(list(m)+normal_exponent(mon),c) for mon,poly in zip(MONS[2],row) for m,c in poly.items()])
        core_write(f,Dw*R2.from_expr(Du.as_expr())**2)


def point_check(data, quad, relation_path, quartic_path):
    prime=32003; point=[2,3,5,7,11,13,17]
    values=dict(zip(s.symbols('e d p q r h k'),point))
    def mod(v):
        v=s.Rational(v); return int(v.p)%prime*pow(int(v.q),-1,prime)%prime
    def evaluate(poly):
        return sum(mod(c)*np.prod([pow(a,e,prime) for a,e in zip(point,m)],dtype=object) for m,c in poly.items())%prime
    T=np.array([[evaluate(v) for v in row] for row in quad['Tn']],dtype=np.int64)
    Du=int(evaluate(quad['Uden']))
    U=np.array([[evaluate(v) for v in row] for row in quad['Un']],dtype=np.int64)*pow(Du,-1,prime)%prime
    J=np.array([[mod(s.sympify(v).subs(values)) for v in row] for row in data['J']],dtype=np.int64)
    beta={int(k[1:])-1:s.sympify(v).subs(values) for k,v in data['beta'].items()}
    _,equations=specialize(S61(),beta,True); H=np.zeros((135,45,45),dtype=np.int64)
    for i,equation in enumerate(equations):
        for mon,c in equation.items():
            if data['fixed'] in mon: continue
            a,b=[data['cols'].index(v) for v in mon]+[44]*(2-len(mon)); c=mod(c)
            if a==b: H[i,a,b]+=c
            else: H[i,a,b]+=c*pow(2,-1,prime)%prime; H[i,b,a]+=c*pow(2,-1,prime)%prime
    H%=prime; Q=np.zeros((135,45),dtype=np.int64)
    for i in range(135):
        matrix=T.T@(H[i]@T%prime)%prime
        for j,(a,b) in enumerate(MONS[2]): Q[i,j]=matrix[a,b]*(1 if a==b else 2)%prime
    selected=data['pivot_rows']; P=J@U%prime; W=-U@Q[selected,:]%prime
    A2=(Q-P@Q[selected,:])%prime; C=np.zeros((135,165),dtype=np.int64)
    indices={d:{m:i for i,m in enumerate(MONS[d])} for d in (3,4)}
    for i in range(135):
        matrix=2*T.T@(H[i]@W%prime)%prime
        for a in range(9):
            for j,mon in enumerate(MONS[2]): C[i,indices[3][tuple(sorted((a,)+mon))]]+=matrix[a,j]
    C%=prime; A3=(C-P@C[selected,:])%prime
    def rref(A):
        A=A.copy()%prime; rank=0; piv=[]
        for col in range(A.shape[1]):
            ids=np.flatnonzero(A[rank:,col])
            if not len(ids): continue
            j=rank+int(ids[0]); A[[j,rank]]=A[[rank,j]]; A[rank]=A[rank]*pow(int(A[rank,col]),-1,prime)%prime
            ids=np.flatnonzero(A[:,col]); ids=ids[ids!=rank]
            if len(ids): A[ids]=(A[ids]-A[ids,col,None]*A[rank,None,:])%prime
            piv.append(col); rank+=1
            if rank==len(A): break
        return A[:rank],piv
    B=A2[np.ix_(PIVOT_ROWS,PIVOT_COLS)]
    aug,piv=rref(np.concatenate([B,np.eye(24,dtype=np.int64)],axis=1)); assert piv==list(range(24)); inverse=aug[:,24:]
    basis=inverse@A2[PIVOT_ROWS,:]%prime; F3=inverse@A3[PIVOT_ROWS,:]%prime
    trueC=(A3-A2[:,PIVOT_COLS]@F3)%prime
    def multiples(A,degree):
        out=np.zeros((len(A)*9,len(MONS[degree+1])),dtype=np.int64)
        for i,row in enumerate(A):
            for a in range(9):
                for j,m in enumerate(MONS[degree]): out[i*9+a,indices[degree+1][tuple(sorted((a,)+m))]]=row[j]
        return out
    M3=multiples(basis,2); M4=multiples(M3,3); R4,piv4=rref(M4); assert len(R4)==453
    lower=np.concatenate([M4,multiples(trueC[[92,134]],3)])
    assert elimination_rank(M3,prime)[0]==129 and elimination_rank(lower,prime)[0]==467
    stream=iter(relation_path.read_text().split()); nr=int(next(stream)); nc=int(next(stream)); assert nr==86 and nc==216
    def core_value():
        result=0
        for _ in range(int(next(stream))):
            c=mod(next(stream)); e=int(next(stream)); d=int(next(stream)); result=(result+c*pow(2,e,prime)*pow(3,d,prime))%prime
        return result
    relations={}
    for _ in range(nr):
        i=int(next(stream)); row=np.zeros(216,dtype=np.int64)
        for _ in range(int(next(stream))):
            j=int(next(stream)); num=core_value(); den=core_value(); row[j]=num*pow(den,-1,prime)%prime
        assert not np.any(row@M3%prime); relations[i]=row
    tails=multiples(F3,3)*pow(Du,-1,prime)%prime; quartics=[]
    with quartic_path.open() as stream:
        assert int(stream.readline())==28
        for expected in SELECTED:
            i=int(stream.readline()); assert i==expected; row=np.zeros(495,dtype=np.int64)
            for _ in range(int(stream.readline())):
                raw=stream.readline().split(); m=list(map(int,raw[1:])); value=mod(raw[0])
                for a,e in zip(point,m[:7]): value=value*pow(a,e,prime)%prime
                mon=tuple(j for j,e in enumerate(m[7:]) for _ in range(e)); row[indices[4][mon]]=(row[indices[4][mon]]+value)%prime
            den=0
            for _ in range(int(stream.readline())):
                raw=stream.readline().split(); m=list(map(int,raw[1:])); value=mod(raw[0]); assert not any(m[7:])
                for a,e in zip(point,m[:7]): value=value*pow(a,e,prime)%prime
                den=(den+value)%prime
            row=row*pow(den,-1,prime)%prime
            delta=(row-relations[i]@tails)%prime; delta=(delta-delta[piv4]@R4)%prime
            assert not np.any(delta),i
            quartics.append(row)
        assert not stream.read().strip()
    full=np.concatenate([lower,np.array(quartics,dtype=np.int64)])
    rank,rows=elimination_rank(full.T,prime); assert rank==495
    determinant=determinant_mod(full[rows,:],prime); assert determinant
    return dict(prime=prime,point=point,quadratic_rank=24,quadratic_linear_rank=129,quadratic_square_rank=453,
                known_lower_quartic_rank=467,all_28_direct_jet_comparisons_pass=True,
                full_quartic_rank=495,minor_rows=rows,minor_mod_prime=determinant)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',required=True)
    parser.add_argument('--flint-include',default=str(ROOT/'work/flint-headers/flint/3.5.0/include'))
    parser.add_argument('--mpfr-include',default=str(ROOT/'work/flint-headers/mpfr/4.2.2/include'))
    parser.add_argument('--flint-library',default=str(ROOT/'work/m2-runtime/lib/libflint.24.0.dylib'))
    args=parser.parse_args(); output=Path(args.output); assert not output.exists()
    started=time.time(); data=json.loads(LINEAR.read_text())['linear']; certificate=json.loads(CERTIFICATE.read_text())
    env=os.environ.copy(); env['DYLD_LIBRARY_PATH']=str(Path(args.flint_library).parent)+':/opt/homebrew/lib'
    dependencies=[Path(__file__),CERTIFICATE,LINEAR,QUADRATIC,
        ROOT/'reports/S61_BINARY_DEEPER_CUBIC_NUISANCE_20261101.json',
        ROOT/'code/s61_binary_deeper_quartic_field_20261101.hpp',
        ROOT/'code/s61_binary_deeper_cubic_nuisance_20261101.py',
        ROOT/'code/s61_binary_deeper_cubic_products_20261101.cpp',
        ROOT/'code/s61_binary_deeper_uniform_quadratic_20261101.py',ROOT/'code/s61_model.py',
        ROOT/'code/s61_flag_transverse_obstruction_20261101.py']
    with tempfile.TemporaryDirectory(prefix='poisson_deeper_quartic_',dir='/tmp') as directory:
        tmp=Path(directory)
        binaries={}
        for name in SOURCES+['cubic_products']:
            source=ROOT/'code'/('s61_binary_deeper_cubic_products_20261101.cpp' if name=='cubic_products' else f's61_binary_deeper_quartic_{name}_20261101.cpp')
            dependencies.append(source); binary=tmp/name
            subprocess.run(['clang++','-O3','-std=c++17','-I'+args.flint_include,'-I'+args.mpfr_include,
                            '-I/opt/homebrew/include',str(source),args.flint_library,'-L/opt/homebrew/lib','-lgmp','-o',str(binary)],check=True)
            binaries[name]=binary
        def run(name,*paths):
            log=tmp/(name+'.log')
            with log.open('w') as stream: subprocess.run([str(binaries[name]),*map(str,paths)],env=env,stdout=stream,stderr=subprocess.STDOUT,check=True)
            print(name,'PASS',round(time.time()-started,3),flush=True)
        quad,_=prepare_cubic(data,certificate,tmp/'cubic_input.txt',tmp/'module3.txt')
        run('cubic_products',tmp/'cubic_input.txt',tmp/'cubics.txt')
        prepare_quartic(data,certificate,quad,tmp)
        run('core_inverse',tmp/'core.txt',tmp/'core_inverse.txt')
        run('nuisance_dag',tmp/'outer.txt',tmp/'core_inverse.txt',tmp/'inverse.txt')
        run('unit_relations',tmp/'module3.txt',tmp/'relations.txt')
        run('constant_span',tmp/'module4.txt',tmp/'quartic_basis.txt')
        run('products',tmp/'operators.txt',tmp/'inverse.txt',tmp/'relations.txt',tmp/'quartics.txt')
        run('nuisance_check',tmp/'quartic_basis.txt',tmp/'quartics.txt',tmp/'nuisance.txt',28)
        witness=point_check(data,quad,tmp/'relations.txt',tmp/'quartics.txt')
        note=ROOT/'research/notes/S61_BINARY_DEEPER_QUARTIC_NUISANCE_20261101.md'
        if note.exists(): dependencies.append(note)
        report=dict(id='S61_BINARY_DEEPER_QUARTIC_NUISANCE_20261101',all_passed=True,
            coefficient_field='QQ(e,d)',nuisance_variables=['p','q','r','h','k'],
            quadratic_inverse_DAG_edges=144,quadratic_inverse_identity_checked=True,
            constant_cubic_pivots=114,quadratic_linear_sources_in_span=200,tracked_linear_syzygies=86,
            selected_quartic_syzygies=SELECTED,quartic_constant_pivots=389,
            all_28_nuisance_coefficients_zero_mod_quadratic_products=True,
            no_raw_fourth_jet_needed=True,point_witness=witness,
            generated_subideal_hilbert_function=[1,9,21,34],generated_subideal_length=65,
            scope='There is a nonempty core open, independent of all five nuisance variables, on which the actual transverse colength is at most65.',
            limitations=['The core exceptional locus has not been exhausted.','No entire deeper-stratum saturation removal or global component exhaustion is asserted.'],
            proof_note=str(note.relative_to(ROOT)),input_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(dependencies))},seconds=time.time()-started)
        with output.open('x') as stream: json.dump(report,stream,ensure_ascii=False,indent=2); stream.write('\n')
    print('PASS',round(time.time()-started,3),flush=True)


if __name__=='__main__': main()
