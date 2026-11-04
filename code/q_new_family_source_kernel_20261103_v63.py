#!/usr/bin/env python3
"""One bounded exact raw-jet 30x26 tangent kernel over five rational parameters."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import shlex
import signal
import subprocess
import time

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63'


def read(p):
    with p.open() as stream:return json.load(stream)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def savegzip(path,payload):
    with gzip.open(path,'wt') as stream:json.dump(payload,stream,indent=2);stream.write('\n')


def cap():
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))


def main(out):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    sourcepath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
    source=read(sourcepath)
    u,v,X,W,B=s.symbols('u v X W B');parameters=[u,v,X,W,B]
    q=s.Matrix(s.symbols('q0:4'));C=s.Matrix(source['raw_value4']);L=s.Matrix(source['raw_linear4'])
    pairs=list(itertools.combinations(range(4),2))
    def jet(A,Z):return s.Matrix([A[i,j] for i,j in pairs]+[s.expand(Z[i,j]).coeff(z) for i,j in pairs for z in q])
    columns=[jet(C.diff(z),L.diff(z)) for z in parameters]
    labels=list(map(str,parameters))
    columns.append(jet(C,L));labels.append('c_d_dc')
    for i,j in itertools.product(range(4),repeat=2):
        A=s.zeros(4);A[i,j]=1
        columns.append(jet(A*C+C*A.T,A*L+L*A.T-sum(((A*q)[k]*L.diff(q[k]) for k in range(4)),s.zeros(4))))
        labels.append(f'GL4_{i}_{j}')
    for i in range(4):
        alpha=s.eye(4)[:,i];z=C*alpha
        columns.append(jet(s.zeros(4),-2*alpha.dot(q)*C+q*z.T-z*q.T));labels.append(f'perspective_{i}')
    M=s.Matrix.hstack(*columns);assert M.shape==(30,26)
    RR,*rv=ring('u,v,X,W,B',s.QQ)
    original=[[RR.from_expr(s.expand(z)) for z in row] for row in M.tolist()]
    matrix=[row.copy() for row in original];rowids=list(range(30));colids=list(range(26));pivots=[];saved=[]
    while matrix and matrix[0]:
        candidates=[]
        for i,row in enumerate(matrix):
            for j,f in enumerate(row):
                if len(f)==1 and (0,)*5 in f:
                    cost=sum(len(z) for z in row)*sum(len(r[j]) for r in matrix)
                    candidates.append((cost,i,j,f[(0,)*5]))
        if not candidates:break
        cost,i,j,p=min(candidates);row=matrix[i]
        saved.append((colids[j],colids.copy(),row.copy(),p))
        pivots.append(dict(original_row=rowids[i],original_column=colids[j],pivot=str(p),
                           current_columns=colids.copy(),row_polynomials=list(map(str,row))))
        matrix=[[f-r[j]/p*row[k] for k,f in enumerate(r) if k!=j] for ii,r in enumerate(matrix) if ii!=i]
        del rowids[i];del colids[j]
    assert len(pivots)==20 and (len(matrix),len(matrix[0]))==(10,6)
    before=dict(parameters=list(map(str,parameters)),columns=labels,source_column_count=6,
                row_convention='six value entries then24 linear coefficients(pair first,variable second)',
                full30x26_polynomial_matrix=[[str(z) for z in row] for row in original],
                all_twenty_QQ_unit_pivots=pivots,remaining_rows=rowids,remaining_columns=colids,
                residual10x6_polynomial_matrix=[[str(z) for z in row] for row in matrix],
                original_scalar_c_handling='Use c*d/dc source column and divide all30 rows by c; valid c!=0; leaves five polynomial parameters.')
    savegzip(out/'matrix-certificate.json.gz',before)
    inp=out/'native-input.txt'
    with inp.open('w') as stream:
        stream.write('5 10 6\n')
        for row in matrix:
            for f in row:
                den=s.ilcm(1,1,*(s.denom(s.Rational(z)) for z in f.values()))
                stream.write(str(len(f))+'\n')
                for mon,z in f.items():stream.write(str(int(s.Rational(z)*den))+' '+' '.join(map(str,mon))+'\n')
                stream.write('1\n'+str(den)+' 0 0 0 0 0\n')
    native=ROOT/'code/rational_function_one_kernel_20261103_v63.cpp'
    binary=ROOT/'work/new25-one-kernel-native-v63'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    subprocess.run(['c++','-std=c++17','-O2',str(native),*flags,'-lgmp','-o',str(binary)],check=True)
    answer=out/'native-kernel.json'
    try:
        with (out/'native.log').open('w') as stream:
            subprocess.run([str(binary),str(inp),str(answer)],stdout=stream,stderr=subprocess.STDOUT,
                           check=True,timeout=62,preexec_fn=cap)
    except Exception as error:
        (out/'blocked.json').write_text(json.dumps(dict(error=repr(error),
            prepared_matrix_sha256=sha(out/'matrix-certificate.json.gz'),native_input_sha256=sha(inp),
            native_run_count=1,no_unchecked_kernel_claimed=True),indent=2)+'\n')
        raise
    result=read(answer);assert result['all_passed'] and result['rank']==5
    loc={str(x):x for x in parameters}
    kernel=[RR.from_expr(s.sympify(z.replace('^','**'),locals=loc)) for z in result['primitive_polynomial_kernel']]
    assert len(kernel)==6 and any(kernel)
    full=dict(zip(colids,kernel))
    for col,ids,row,p in reversed(saved):
        full[col]=-sum((f*full[k] for k,f in zip(ids,row) if k!=col),RR.zero)/p
    vec=[full[k] for k in range(26)]
    den=s.ilcm(1,*(s.denom(s.Rational(c)) for f in vec for c in f.values()))
    vec=[f*RR.domain.convert(den) for f in vec]
    content=s.igcd(*(int(c) for f in vec for c in f.values()))
    vec=[f/RR.domain.convert(content) for f in vec]
    assert all(sum((f*k for f,k in zip(row,vec)),RR.zero)==RR.zero for row in original)
    assert any(vec[:6])
    source_vec=dict(zip(labels[:6],map(str,vec[:6])))
    after=dict(full26_polynomial_kernel=list(map(str,vec)),source_six_direction=source_vec,
               interpretation='du,dv,dX,dW,dB, dlog(c); actual dc=c*dlog(c)',
               compensating_GL4_16=list(map(str,vec[6:22])),compensating_perspective4=list(map(str,vec[22:])),
               all30_original_polynomial_dot_products_exactly_zero=True,
               native_primitive_polynomial_kernel=result['primitive_polynomial_kernel'])
    savegzip(out/'kernel-certificate.json.gz',after)
    files=[Path(__file__),native,sourcepath,ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/report.json']
    report=dict(id=STEM,all_passed=True,field='QQ(u,v,X,W,B)',rawjet_matrix_shape=[30,26],
                rational_unit_pivots=20,residual_shape=[10,6],residual_exact_function_field_rank=5,
                full_exact_function_field_rank=25,one_dimensional_generic_kernel=True,
                source_part_nonzero=True,full26_polynomial_kernel_and_all30_identities_verified=True,
                no_point_evaluation_no_high_degree_fit_no_GB=True,
                five_rational_first_integrals_or_quotient_or_integrability_claimed=False,
                native_run_count=1,source_kernel_direction=source_vec,
                input_sha256={str(p.relative_to(ROOT)):sha(p) for p in files},
                evidence_sha256={p.name:sha(p) for p in [inp,answer,out/'native.log',out/'matrix-certificate.json.gz',out/'kernel-certificate.json.gz']},
                seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','full_exact_function_field_rank','source_kernel_direction','seconds']},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output.resolve())
