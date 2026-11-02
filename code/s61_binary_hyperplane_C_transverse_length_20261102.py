#!/usr/bin/env python3
"""Prove exact transverse length65 at the rational C-root sample."""
from fractions import Fraction
from pathlib import Path
import json
import math
import subprocess
import tempfile
import time

import numpy as np
import sympy as sp

from s61_binary_deeper_uniform_quadratic_20261101 import determinant_mod
from s61_binary_hyperplane_cubic_matrix_20261102 import compile_native
from s61_binary_hyperplane_root_model_20261102 import (
    ROOT,CODE,P,S61,actual_binary_tangent,digest,hashes,old_cut_rows,original_point,
    output_argument,residue,runtime_environment,source_paths,write_new,
)

CERT=CODE/'certificates/s61_binary_hyperplane_C_transverse_pivots_20261102.json'
NATIVE=CODE/'s61_binary_hyperplane_C_integer_rank_20261102.cpp'
NOTE=ROOT/'research/notes/S61_BINARY_HYPERPLANE_C_TRANSVERSE_LENGTH_20261102.md'


def add(poly, monomial, value):
    result=poly.get(monomial,0)+value
    if result:
        poly[monomial]=result
    else:
        poly.pop(monomial,None)


def multiply(left,right):
    product={}
    for m,c in left.items():
        for n,v in right.items():
            add(product,tuple(sorted(m+n)),c*v)
    return product


def substitute(poly,variable,replacement):
    degree=max((m.count(variable) for m in poly),default=0)
    if not degree:
        return poly.copy()
    powers=[{():Fraction(1)},replacement]
    for _ in range(2,degree+1):
        powers.append(multiply(powers[-1],replacement))
    result={}
    for m,c in poly.items():
        tail=tuple(j for j in m if j!=variable)
        for n,v in powers[m.count(variable)].items():
            add(result,tuple(sorted(tail+n)),c*v)
    return result


def monomials(variables,weights,weight,first=0):
    if weight==0:
        yield ()
        return
    for k in range(first,len(variables)):
        j=variables[k]
        if weights[j]<=weight:
            for tail in monomials(variables,weights,weight-weights[j],k):
                yield (j,)+tail


def verify():
    start=time.time()
    cert=json.loads(CERT.read_text())
    assert cert['point']==['-13/25','2','3'] and cert['cuts']==[20,19,23]
    model=S61()
    beta,two,old,weights,row_weights,full=original_point(model,Fraction(-13,25))
    jacobian,direction,tangent=actual_binary_tangent(model,beta,two,old)
    linear=sp.Matrix([[q.get((j,),0) for j in range(45)] for q in full])
    assert linear*tangent==sp.zeros(135,3)
    assert tangent[old_cut_rows(old),:].det()==0
    cuts=cert['cuts']
    transverse=tangent[cuts,:].det()
    assert transverse==sp.Rational(79,25) and residue(transverse)
    equations=[{m:c for m,c in q.items() if not set(m).intersection(cuts)}
               if row_weights[i]<=4 else {} for i,q in enumerate(full)]
    variables=[j for j in range(45) if j not in cuts]
    pivots=[]
    for i,j in cert['pivot_steps']:
        q=equations[i]
        pivot=q[(j,)]
        assert residue(pivot) and all(j not in m for m in q if m!=(j,))
        replacement={m:-c/pivot for m,c in q.items() if m!=(j,)}
        assert all(c.denominator%P for c in replacement.values())
        equations=[{} if k==i else substitute(f,j,replacement) for k,f in enumerate(equations)]
        variables.remove(j)
        pivots.append(str(pivot))
    assert len(pivots)==34 and sorted(weights[j] for j in variables)==[1]*7+[2]
    for i,q in enumerate(equations):
        assert all(len(m)>=2 and set(m).issubset(variables)
                   and sum(weights[j] for j in m)==row_weights[i] for m in q)
    mons={w:list(monomials(variables,weights,w)) for w in range(5)}
    piece={str(w):{'columns':len(mons[w]),'rank_QQ':0,'dimension':len(mons[w])} for w in [0,1]}
    with tempfile.TemporaryDirectory(prefix='poisson-C-length65-',dir='/tmp') as tmp:
        wd=Path(tmp)
        matrix_path=wd/'matrices'
        output=wd/'ranks'
        with matrix_path.open('w') as stream:
            for w in [2,3,4]:
                index={m:i for i,m in enumerate(mons[w])}
                rows=[]
                for i,q in enumerate(equations):
                    if not q or row_weights[i]>w:
                        continue
                    for shift in mons[w-row_weights[i]]:
                        rows.append({index[tuple(sorted(m+shift))]:c for m,c in q.items()})
                stream.write(f'{w} {len(rows)} {len(index)}\n')
                A=np.zeros((len(rows),len(index)),dtype=np.int64)
                for i,row in enumerate(rows):
                    scale=math.lcm(*(c.denominator for c in row.values()))
                    assert scale%P
                    stream.write(str(len(row)))
                    for j,c in sorted(row.items()):
                        stream.write(f' {j} {int(scale*c)}')
                        A[i,j]=residue(c)
                    stream.write('\n')
                minor=cert['minors'][str(w)]
                value=int(determinant_mod(A[np.ix_(minor['minor_rows'],minor['minor_columns'])],P))
                assert value==minor['minor_modp'] and value
                piece[str(w)]={'rows':len(rows),'columns':len(index),'minor_modp':value}
        exe=wd/'rank'
        compile_native(NATIVE,exe)
        subprocess.run([str(exe),str(matrix_path),str(output)],env=runtime_environment(),
                       check=True,capture_output=True,text=True,timeout=50)
        for line in output.read_text().splitlines():
            w,rows,cols,rank=map(int,line.split())
            data=piece[str(w)]
            assert (rows,cols)==(data['rows'],data['columns'])
            data.update(rank_QQ=rank,dimension=cols-rank)
        generated={'integer_matrices':digest(matrix_path),'rank_output':digest(output)}
    hf=[piece[str(w)]['dimension'] for w in range(5)]
    assert hf==[1,7,23,34,0]
    static=source_paths()+[Path(__file__),CERT,NATIVE,NOTE,
                           CODE/'s61_binary_hyperplane_cubic_matrix_20261102.py',
                           CODE/'s61_binary_deeper_uniform_quadratic_20261101.py']
    return {'all_passed':True,'id':'S61_BINARY_HYPERPLANE_C_TRANSVERSE_LENGTH_20261102',
            'scope':'Exact transverse coordinate algebra at a=-13/25,e=2,d=3 in the three-core slice.',
            'point':cert['point'],'cuts':cuts,'old_cut_determinant':'0',
            'new_cut_determinant_QQ':str(transverse),'new_cut_determinant_modp':residue(transverse),
            'curve_jacobian_QQ':[[str(x) for x in row] for row in jacobian.tolist()],
            'plane_tangent_QQ':list(map(str,direction)),
            'actual_binary_tangent_45x3_QQ':[[str(x) for x in row] for row in tangent.tolist()],
            'all135_linear_rows_annihilate_tangent':True,'pivot_values_QQ':pivots,
            'remaining_variables':variables,'remaining_weights':[weights[j] for j in variables],
            'weighted_pieces':piece,'weighted_Hilbert_function_QQ':hf,'length_QQ':sum(hf),
            'all_weights_at_least_four_zero':True,'input_sha256':hashes(static),
            'generated_sha256':generated,'seconds':time.time()-start}


if __name__=='__main__':
    path=output_argument(__doc__)
    report=verify()
    write_new(path,report)
    print('PASS',report['weighted_Hilbert_function_QQ'],report['length_QQ'],report['seconds'])
