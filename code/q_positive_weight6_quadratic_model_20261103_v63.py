#!/usr/bin/env python3
"""Full physical350 equation, orbit23 slice and exact eight-normal quadratic data."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(5),2))
MON2=list(itertools.combinations_with_replacement(range(5),2))
MON3=list(itertools.combinations_with_replacement(range(5),3))
TRIPLES=list(itertools.combinations(range(5),3))
BASIS=list(itertools.product(PAIRS,MON2))


def read(path):
    data=path.read_bytes()
    return json.loads(gzip.decompress(data)if path.suffix=='.gz'else data)


def pivot_columns(M):return list(DomainMatrix.from_Matrix(M).rref()[1])
def matrix(M):return [[str(v)for v in row]for row in M.tolist()]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_linearization(q):
    point={}
    for (pair,mon),v in zip(BASIS,q):
        if v:point.setdefault(pair,{})[mon]=v
    def get(i,j):
        return point.get((i,j),{})if i<j else({m:-v for m,v in point.get((j,i),{}).items()}if i>j else{})
    positions={m:i for i,m in enumerate(MON3)}
    J=s.zeros(350,150)
    for col,((u,v),mon)in enumerate(BASIS):
        def sign(i,j):return int((i,j)==(u,v))-int((j,i)==(u,v))
        for r,(i,j,k)in enumerate(TRIPLES):
            for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
                for ell in range(5):
                    sg=sign(c,d)
                    if sg and ell in mon:
                        rem=list(mon);rem.remove(ell)
                        for other,coefficient in get(b,ell).items():
                            J[35*r+positions[tuple(sorted(other+tuple(rem)))],col]+=sg*mon.count(ell)*coefficient
                    sg=sign(b,ell)
                    if sg:
                        for other,coefficient in get(c,d).items():
                            if ell in other:
                                rem=list(other);rem.remove(ell)
                                J[35*r+positions[tuple(sorted(mon+tuple(rem)))],col]+=sg*other.count(ell)*coefficient
    return J


def physical_quadratic(T):
    point={}
    for row,(pair,mon)in enumerate(BASIS):
        coefficients={j:T[row,j]for j in range(T.cols)if T[row,j]}
        if coefficients:point.setdefault(pair,{})[mon]=coefficients
    def get(i,j):
        return point.get((i,j),{})if i<j else({m:{k:-v for k,v in c.items()}for m,c in point.get((j,i),{}).items()}if i>j else{})
    positions={m:i for i,m in enumerate(MON3)}
    monomials=list(itertools.combinations_with_replacement(range(T.cols),2))
    indices={m:i for i,m in enumerate(monomials)}
    Q=s.zeros(350,len(monomials))
    for r,(i,j,k)in enumerate(TRIPLES):
        for b,c,d in [(i,j,k),(j,k,i),(k,i,j)]:
            for ell in range(5):
                for mon,left in get(b,ell).items():
                    for other,right in get(c,d).items():
                        if ell not in other:continue
                        rem=list(other);rem.remove(ell)
                        row=35*r+positions[tuple(sorted(mon+tuple(rem)))]
                        multiplicity=other.count(ell)
                        for p,a in left.items():
                            for q,z in right.items():
                                Q[row,indices[tuple(sorted((p,q)))]]+=multiplicity*a*z
    return Q,monomials


def main(output,run_m2=False):
    started=time.monotonic()
    output.mkdir(parents=True,exist_ok=False)
    basis_path=ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63/certificate.json.gz'
    orbit_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63/certificate.json'
    data=read(basis_path);C=s.Matrix(data['canonical_basis']);free=data['canonical_free_raw_indices']
    assert C.shape==(150,126) and C[free,:]==s.eye(126)
    q=s.zeros(150,1)
    for i,v in {64:-s.Rational(256,3),83:-64,90:s.Rational(8,9),105:-s.Rational(2,9),116:64,121:s.Rational(4,3),137:s.Rational(2,3),140:1}.items():q[i]=v
    assert C*q[free,0]==q
    N=s.Matrix(read(orbit_path)['GL25_plus_scalar_matrix'])
    assert N[:,25]==q and C*N[free,:]==N
    orbit_columns=pivot_columns(N);assert len(orbit_columns)==23
    O=N[:,orbit_columns]
    cuts=pivot_columns(O[free,:].T);assert len(cuts)==23
    group_minor=O[[free[i]for i in cuts],:].det();assert group_minor
    remaining=[i for i in range(126)if i not in cuts]
    S=C[:,remaining]
    L=physical_linearization(q)
    assert L*q==s.zeros(350,1) and L*N==s.zeros(350,26)
    J=L*S
    columns=pivot_columns(J);assert len(columns)==95
    rows=pivot_columns(J[:,columns].T);assert len(rows)==95
    B=J[rows,columns]
    inverse=DomainMatrix.from_Matrix(B).inv().to_Matrix()
    assert B*inverse==s.eye(95)
    normal=[i for i in range(103)if i not in columns];assert len(normal)==8
    T=s.zeros(103,8)
    for j,k in enumerate(normal):T[k,j]=1
    correction=-inverse*J[rows,normal]
    for j,k in enumerate(columns):T[k,:]=correction[j,:]
    raw=S*T
    assert J*T==s.zeros(350,8) and raw[[free[i]for i in cuts],:]==s.zeros(23,8)
    print('FULL350 linear slice: orbit23, implicit95, tangent8',time.monotonic()-started,flush=True)
    Q,monomials=physical_quadratic(raw)
    implicit_quadratic=-inverse*Q[rows,:]
    obstruction=Q+J[:,columns]*implicit_quadratic
    assert obstruction[rows,:]==s.zeros(95,36)
    rrref=DomainMatrix.from_Matrix(obstruction).rref()[0].to_Matrix()
    qrank=sum(any(v for v in rrref.row(i))for i in range(rrref.rows))
    quadratic_basis=rrref[:qrank,:]
    z=s.symbols('z0:8')
    polys=[str(s.factor(sum(coefficient*z[i]*z[j]for coefficient,(i,j)in zip(row,monomials))))for row in quadratic_basis.tolist()]
    payload=dict(id='Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63',
        physical_center150=list(map(str,q)),canonical_basis_dependency=str(basis_path.relative_to(ROOT)),
        canonical_free_raw_indices=free,orbit_columns_zero_based=orbit_columns,
        orbit_cut_canonical_indices=cuts,orbit_cut_raw_indices=[free[i]for i in cuts],orbit_cut_minor=str(group_minor),
        remaining103_canonical_indices=remaining,implicit_equation_rows=rows,implicit_coordinate_columns=columns,
        implicit_inverse95=matrix(inverse),normal8_coordinate_columns=normal,
        normal8_raw_indices=[free[remaining[i]]for i in normal],tangent150_by8=matrix(raw),
        full350_by103_linear_Jacobian=matrix(J),normal_quadratic_monomials=monomials,
        full350_by36_raw_quadratic=matrix(Q),implicit95_by36_second_order_coefficients=matrix(implicit_quadratic),
        full350_by36_corrected_obstructions=matrix(obstruction),quadratic_row_basis=matrix(quadratic_basis),
        quadratic_basis_polynomials=polys)
    archive=output/'certificate.json.gz'
    with archive.open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as out:out.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    m2='R=QQ['+','.join(map(str,z))+'];\nI=ideal('+','.join(p.replace('**','^')for p in polys)+');\nG=gb I;\n'
    m2+='print("DIM="|toString dim(R/I));print("DEG="|toString degree(R/I));print("HS="|toString hilbertSeries(R/I));\n'
    m2+='K=ideal(z1,z7,z5*(z4+2*z6));assert(gens I%gb K==0);\n'
    m2+='assert(z1^3%G==0);assert(z7^3%G==0);assert((z5*(z4+2*z6))^2%G==0);\n'
    m2+='assert(dim(R/I)==5);assert(degree(R/I)==3);print("EXACT_RADICAL_TWO_FIVE_PLANES_PASS");\n'
    (output/'quadrics.m2').write_text(m2)
    if run_m2:
        executable=shutil.which('M2');assert executable
        with (output/'quadrics.log').open('w')as log:
            subprocess.run([executable,'--script',str(output/'quadrics.m2')],stdout=log,stderr=subprocess.STDOUT,timeout=30,check=True)
        assert 'EXACT_RADICAL_TWO_FIVE_PLANES_PASS'in(output/'quadrics.log').read_text()
    result=dict(id=payload['id'],all_passed=True,all350_physical_strict_Jacobi_equations_retained=True,
        physical_canonical_tangent_dimension=31,projective_tangent_dimension=30,GL_plus_scalar_orbit_dimension=23,
        orbit_cut23_unit_minor=str(group_minor),linear_implicit_rank95=True,formal_transverse_embedding_dimension=8,
        exact_quadratic_obstruction_rank=qrank,quadratic_obstruction_polynomials=polys,
        exact_quadratic_cone_GB_run=run_m2,
        quadratic_cone_radical='(z1,z7,z5*(z4+2*z6))'if run_m2 else None,
        quadratic_cone_dimension=5 if run_m2 else None,quadratic_cone_degree=3 if run_m2 else None,
        full_actual_transverse_dimension_upper_bound=5 if run_m2 else None,
        full_actual_projective_local_dimension_upper_bound=27 if run_m2 else None,
        formal8_model='Solve selected95 ORIGINAL equations formally for95 dependent coordinates, then retain the other255 original residual series.',
        finite_polynomial8_model_full_initial_ideal_dimension_or_branch_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),basis_path,orbit_path,
            ROOT/'reports/Q_FULL_FAMILY_LINEAR_SMITH_20261103_v63_REVIEW.json',ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json',
            ROOT/'research/notes/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63.md']},
        certificate_sha256=sha(archive),seconds=time.monotonic()-started)
    result['evidence_sha256']={p.name:sha(p)for p in output.iterdir()if p.is_file()and p.name!='report.json'}
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('EXACT8 NORMAL QUADRATIC OBSTRUCTION RANK',qrank,'SECONDS',result['seconds'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--m2',action='store_true',help='Run the bounded exact eight-variable quadratic-cone certificate')
    args=parser.parse_args();main(args.output,args.m2)
