#!/usr/bin/env python3
"""Independent universal Inv24 basis from full quartic Hessian and formal coefficients."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
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

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63_REVIEW'
P=32003
M2=[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)]

def read(p):
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def mod(q):
    q=F(str(q));assert q.denominator%P;return q.numerator*pow(q.denominator,-1,P)%P
def mm(A,B):return [[sum(x*y for x,y in zip(row,col))%P for col in zip(*B)]for row in A]
def det(A):
    A=[r.copy()for r in A];d=1
    for j in range(len(A)):
        i=next(i for i in range(j,len(A))if A[i][j]);
        if i!=j:A[i],A[j]=A[j],A[i];d=-d%P
        d=d*A[j][j]%P;inv=pow(A[j][j],-1,P)
        for i in range(j+1,len(A)):
            c=A[i][j]*inv%P
            if c:A[i]=[(x-c*y)%P for x,y in zip(A[i],A[j])]
    return d
def rank(columns):
    basis={}
    for row in columns:
        row=row.copy()
        for p,b in sorted(basis.items()):
            c=row[p]
            if c:row=[(x-c*y)%P for x,y in zip(row,b)]
        p=next((i for i,x in enumerate(row)if x),None)
        if p is not None:
            z=pow(row[p],-1,P);basis[p]=[x*z%P for x in row]
    return len(basis)
def pm(a,b):
    out=defaultdict(int)
    for (i,j),x in a.items():
        for (k,l),y in b.items():out[i+k,j+l]=(out[i+k,j+l]+x*y)%P
    return {m:x for m,x in out.items()if x}
def main(out,work):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63'
    ar,cert,witness=[read(folder/p)for p in ['report.json','native-certificate.json.gz','minor-witnesses.json']]
    info_path=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    spin_path=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    plane_path=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
    char_path=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    info,spin,plane,char=map(read,[info_path,spin_path,plane_path,char_path]);assert ar['all_passed']
    gates=[ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json',ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/report.json']
    assert all(read(p)['all_passed']for p in gates)
    b=plane['ambient_beta16'][:3];beta=[{(0,0):b[0][i]%P,(1,0):b[1][i]%P,(0,1):b[2][i]%P}for i in range(16)]
    Omega=s.Matrix(info['symplectic_form16']);Oi=Omega.inv();rho=[s.Matrix(r)for r in info['all_ten_representation_matrices']]
    J=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]));g=[]
    for i in range(4):
        for j in range(i,4):
            z=s.zeros(4);z[i,j]=z[j,i]=1;g.append(J*z)
    Gi=s.Matrix([[s.trace(a*z)for z in g]for a in g]).inv()
    assert all(a.T*Omega+Omega*a==s.zeros(16)for a in rho)
    products={(i,j):pm(beta[i],beta[j])for i in range(16)for j in range(i,16)}
    def quad(coefs):
        z=defaultdict(int)
        for (i,j),c in zip(spin['quadratic_monomials'],coefs):
            c=mod(c)
            for m,v in products[i,j].items():z[m]=(z[m]+c*v)%P
        return {m:v for m,v in z.items()if v}
    CD=[]
    for k in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
        pol=[quad(row)for row in spin[k]];CD.append([pol[4*i:4*i+4]for i in range(4)])
    # theta via its original moment definition; constant coefficient tensors.
    mu=[]
    Sforms=[Omega*r for r in rho]
    for A in Sforms:
        assert A==A.T
        co=[A[i,j]/(2 if i==j else 1)for i,j in spin['quadratic_monomials']]
        mu.append(quad(co))
    co=[]
    for row in Gi.tolist():
        z=defaultdict(int)
        for c,p in zip(row,mu):
            for m,v in p.items():z[m]=(z[m]+mod(c)*v)%P
        co.append(dict(z))
    theta=[[{}for _ in range(16)]for _ in range(16)]
    for i in range(16):
        for j in range(16):
            z=defaultdict(int)
            for r,p in zip(rho,co):
                for m,v in p.items():z[m]=(z[m]+mod(r[i,j])*v)%P
            theta[i][j]=dict(z)
    # H via the DIFFERENT entrance Omega^-1 Hessian(F4/2), using all117
    # independently accepted quartic coefficients, not the author's V-G-V^T.
    Hess=[[defaultdict(int)for _ in range(16)]for _ in range(16)]
    for indices,coefficient in spin['quartic_trace_moment_square']:
        counts=Counter(indices);coef=mod(F(coefficient)/2)
        for i,ei in counts.items():
            for j,ej in counts.items():
                factor=ei*(ej-(i==j))
                if not factor:continue
                rest=list(indices);rest.remove(i);rest.remove(j);k,l=sorted(rest)
                for m,v in products[k,l].items():Hess[i][j][m]=(Hess[i][j][m]+coef*factor*v)%P
    H=[[{}for _ in range(16)]for _ in range(16)]
    for i in range(16):
        for j in range(16):
            z=defaultdict(int)
            for k in range(16):
                for m,v in Hess[k][j].items():z[m]=(z[m]+mod(Oi[i,k])*v)%P
            H[i][j]=dict(z)
    # Independently check the exact Hessian identity against all6 formal
    # coefficient matrices obtained from the original low-rank moment formula.
    om=[[mod(x)for x in r]for r in Omega.tolist()];oi=[[mod(x)for x in r]for r in Oi.tolist()];gi=[[mod(x)for x in r]for r in Gi.tolist()]
    Smod=[[[mod(x)for x in row]for row in A.tolist()]for A in Sforms]
    vcoef=[]
    for q in range(3):vcoef.append([[sum(A[i][j]*b[q][j]for j in range(16))%P for A in Smod]for i in range(16)])
    for mon in M2:
        extra=[[0]*16 for _ in range(16)]
        exps=[(0,0),(1,0),(0,1)]
        for q in range(3):
            for r in range(3):
                if tuple(a+b for a,b in zip(exps[q],exps[r]))!=mon:continue
                Z=mm(mm(vcoef[q],gi),list(map(list,zip(*vcoef[r]))))
                extra=[[(a+b)%P for a,b in zip(aa,bb)]for aa,bb in zip(extra,Z)]
        extra=mm(oi,extra)
        assert all(H[i][j].get(mon,0)%P==(theta[i][j].get(mon,0)+extra[i][j])%P for i in range(16)for j in range(16))
    inp=out/'own-formal-input.txt'
    with inp.open('w')as f:
        for A in [*CD,theta,H]:
            for row in A:
                for p in row:f.write(' '.join(str(p.get(m,0)%P)for m in M2)+'\n')
    native=ROOT/'code/q_new_family_full24_basis_review_20261104_v63.cpp';binary=work/'formal-basis-review'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    with(out/'compile.log').open('w')as f:subprocess.run(['c++','-O2','-std=c++17',str(native),*flags,'-lgmp','-o',str(binary)],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=30)
    with(out/'native.log').open('w')as f:subprocess.run([str(binary),str(inp),str(out/'own-columns.json')],check=True,stdout=f,stderr=subprocess.STDOUT,timeout=40)
    own=read(out/'own-columns.json');assert own['all_passed']
    for degree in [12,16,24]:assert own[f'degree{degree}_columns']==cert[f'degree{degree}_columns']
    assert rank(own['degree16_columns'][:24])==19
    assert rank(own['degree16_columns'][:60])==19
    minors={}
    for degree in [12,16,24]:
        w=witness[str(degree)];columns=own[f'degree{degree}_columns'];A=[[columns[j][i]for j in w['columns']]for i in w['coefficient_rows']]
        minors[str(degree)]=det(A);assert minors[str(degree)]==w['minor_mod32003']
    assert minors=={'12':20149,'16':2363,'24':26663}
    assert witness['24']['columns']==[i for i in range(70)if i not in [62,68,69]]
    dims={int(r['degree']):r['invariant_dimension']for r in char['character_records_through24']}
    assert [dims[d]for d in [12,16,24]]==[9,20,67]
    paths=[Path(__file__),native,ROOT/f'research/notes/{STEM}.md',info_path,spin_path,plane_path,char_path,*gates,
        *[folder/p for p in ['report.json','native-certificate.json.gz','minor-witnesses.json','definitions.json']],
        ROOT/'code/q_new_family_degree12_14_invariant_exclusion_20261103_v63.py',
        ROOT/'research/notes/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63.md']
    result=dict(id=STEM,all_passed=True,independent_actual_constant_rep_Omega_Gram_and_same_old_formal_plane=True,
        independent_H16_from_full117_term_quartic_Hessian_matches_original_formula=True,
        separate_reversed_variable_LEX_engine_matches_all34_96_70_coefficient_columns=True,
        exact_coefficient_minors=minors,full_Inv12_9_Inv16_20_Inv24_67_QQ_bases_by_character_upper_bounds=True,
        independent_degree16_product_and_product_plus_all_CD_trace8_rank19=True,
        selected67_basis70_except=[62,68,69],CD_only_rank8_not_used_as_global_upper=True,
        no_NEW_point_source_expansion_or_NEW_equation_claim=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k]for k in ['all_passed','exact_coefficient_minors','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path);a=p.parse_args();main(a.output.resolve(),a.work.resolve())
