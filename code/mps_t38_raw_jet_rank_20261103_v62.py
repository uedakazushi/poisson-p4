#!/usr/bin/env python3
"""T38 first-jet image rank before rational value/perspective normalization.

The full constant/linear jet rank equals primitive first-jet rank+6+4 on
the nondegenerate-value locus. All polynomial entries have degree at most3.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s
from s61_model import X,PAIRS
from mps_s61_dominance_20261103_v62 import rank_minor
from known_rank4_firstjet_native_rank_20261103_v62 import write_matrix
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);parser.add_argument('--binary',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    source=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json';body=json.loads(source.read_text());P=s.zeros(5)
    for pair,expression in body['homogeneous_bivector'].items():
        i,j=map(int,pair.split(','));P[i,j]=s.sympify(expression);P[j,i]=-P[i,j]
    variables=sorted(set().union(*(v.free_symbols for v in P))-set(X),key=str);sample={v:i+2 for i,v in enumerate(variables)};P=P.subs(sample)
    y=s.Matrix(s.symbols('b0:4'));substitution={X[0]:1}|dict(zip(X[1:],y))
    affine=s.Matrix(4,4,lambda i,j:s.expand((P[i+1,j+1]-X[i+1]*P[0,j+1]+X[j+1]*P[0,i+1]).subs(substitution)))
    C=affine;L=[affine.diff(v) for v in y];pairs=list(itertools.combinations(range(4),2))
    def vector(c,l):return s.Matrix([c[i,j] for i,j in pairs]+[lk[i,j] for lk in l for i,j in pairs])
    columns=[];labels=[]
    for i in range(4):
        for j in range(4):
            g=s.zeros(4);g[i,j]=1
            columns.append(vector(g*C+C*g.T,[g*l+l*g.T-sum((g[h,k]*L[h] for h in range(4)),s.zeros(4)) for k,l in enumerate(L)]));labels.append(f'gl_{i}_{j}')
    for j in range(4):
        a=s.eye(4)[:,j];ca=C*a
        columns.append(vector(s.zeros(4),[-2*(1 if j==k else 0)*C+s.eye(4)[:,k]*ca.T-ca*s.eye(4)[:,k].T for k in range(4)]));labels.append(f'perspective_{j}')
    for j in range(4):columns.append(vector(C.diff(y[j]),[l.diff(y[j]) for l in L]));labels.append(f'point_{j}')
    columns.append(vector(C,L));labels.append('scalar')
    matrix=s.Matrix.hstack(*columns).applyfunc(s.expand);minor=rank_minor(matrix.subs(dict(zip(y,[2,3,4,5]))))
    print('T38 RAW JET MODULAR LOWER',minor['rank'],flush=True)
    K=s.QQ.frac_field(*map(str,y));data=[[K.from_sympy(v) for v in row] for row in matrix.tolist()]
    input_path=args.output/'T38.in';write_matrix(input_path,K,data)
    with (args.output/'T38-native.log').open('w') as log:
        subprocess.run([str(args.binary.resolve()),str(input_path.resolve())],stdout=log,stderr=subprocess.STDOUT,timeout=300,check=True)
    log=(args.output/'T38-native.log').read_text();rank=int(next(line for line in log.splitlines() if line.startswith('EXACT_FUNCTION_FIELD_RANK ')).split()[1])
    assert rank==minor['rank'];encoded=input_path.read_bytes();archive=args.output/'T38.in.gz'
    with archive.open('wb') as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0) as out:out.write(encoded)
    input_path.unlink();inputs=[Path(__file__),source,ROOT/'code/rational_function_matrix_rank_20261103_v62.cpp',ROOT/'code/known_rank4_firstjet_native_rank_20261103_v62.py']
    report=dict(id='MPS_T38_RAW_JET_RANK_20261103_v62',all_passed=True,constant_linear_jet_image_dimension=rank,
        exact_primitive_firstjet_image_dimension=rank-10,s61_component_dimension=15,generic_firstjet_fiber_dimension=25-rank,
        normalization_dimension_removed=dict(value=6,perspective=4),column_order=labels,
        full_polynomial_matrix=[[str(v) for v in row] for row in matrix.tolist()],function_variables=list(map(str,y)),
        modular_lower_witness=minor,fixed_parameter_values={str(k):str(v) for k,v in sample.items()},
        proof_dependencies=['The same fixed source tensor has dense projective PGL orbit by the separate rank25 certificate.',
            'GL4 normalizes all6 nondegenerate values with residualSp4; perspective removes4 free trace directions.',
            'Scalar parameter is retained through its Lie column at1; uniform scalar action preserves all ranks.'],
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        native_input_sha256=hashlib.sha256(encoded).hexdigest(),native_archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
        native_log_sha256=hashlib.sha256((args.output/'T38-native.log').read_bytes()).hexdigest(),seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('T38 EXACT RAW JET RANK PASS',rank,'PRIMITIVE',rank-10,report['seconds'],flush=True)
