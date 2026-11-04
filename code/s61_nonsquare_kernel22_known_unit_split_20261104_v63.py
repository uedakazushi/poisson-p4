#!/usr/bin/env python3
"""One bounded known-unit Gaussian preprocessing on the frozen sparse K only."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import shlex
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_KERNEL22_KNOWN_UNIT_SPLIT_20261104_v63'
BASE=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    kernel=BASE/'kernel22.json';data=read(kernel);parent=read(BASE/'report.json');assert sha(kernel)==parent['evidence_sha256']['kernel22.json']
    review=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63_REVIEW/report.json';assert read(review)['all_passed']
    atlas=read(BASE/'atlas.json');assert data['shape']==[174,22]
    zero=[[],[[[0,0,0],'1']]];matrix=[[zero for _ in range(22)]for _ in range(174)]
    for j,col in enumerate(data['sparse_columns']):
        for i,num,den in col:matrix[i][j]=[num,den]
    inp=out/'native-input.txt'
    def putpoly(f,p):
        f.write(str(len(p)))
        for ex,c in p:f.write(' '+c+' '+' '.join(map(str,ex)))
        f.write('\n')
    with inp.open('w')as f:
        f.write('174 22\n')
        for row in matrix:
            for num,den in row:putpoly(f,num);putpoly(f,den)
    source=ROOT/'code/s61_nonsquare_kernel22_known_unit_split_20261104_v63.cpp';exe=out/'native'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    compile_args=['c++','-O2','-std=c++17',str(source),*flags,'-lgmp','-o',str(exe)]
    cp=subprocess.run(compile_args,capture_output=True,text=True,check=True)
    (out/'compile.log').write_text(cp.stdout+cp.stderr)
    print(json.dumps(dict(stage='native_input_bound_and_compiled',seconds=time.monotonic()-began)),flush=True)
    state=out/'native-state.txt';log=out/'native-proof.log';timedout=False
    try:cp=subprocess.run([str(exe),str(inp),str(state),str(log)],capture_output=True,text=True,timeout=55)
    except subprocess.TimeoutExpired as ex:
        timedout=True;cp=None;(out/'native-timeout.json').write_text(json.dumps(dict(timeout_seconds=55,automatic_retry=False))+'\n')
    if cp is not None:(out/'native-process.log').write_text(cp.stdout+cp.stderr)
    assert state.exists(),'No completed checkpoint; no retry authorized'
    # Read only small A,L,U prefixes; the full Schur state remains immutable.
    with state.open()as f:
        n,m,r=map(int,f.readline().split());rows=list(map(int,f.readline().split()));cols=list(map(int,f.readline().split()))
        status,verified,secs,terms=f.readline().split();verified=verified=='1';assert len(rows)==len(cols)==r
        def getpoly():
            t=f.readline().split();k=int(t[0]);assert len(t)==1+4*k
            return[[list(map(int,t[2+4*j:5+4*j])),t[1+4*j]]for j in range(k)]
        def getrat():return[getpoly(),getpoly()]
        def getmat():return[[getrat()for _ in range(r)]for _ in range(r)]
        A,L,U=getmat(),getmat(),getmat()
    assert n==174 and m==22
    full=(r==22 and verified and not timedout and cp is not None and cp.returncode==0)
    # Exact finite scalar circuit for P_selected=U^-1 L^-1 projection_rows.
    # Divisions are ONLY the saved known-unit pivots (U diagonal).
    atoms=[];aidx={};nodes=[['zero']];idx={('zero',):0}
    def atom(x):
        key=json.dumps(x,separators=(',',':'))
        if key not in aidx:aidx[key]=len(atoms);atoms.append(x)
        i=aidx[key];key2=('atom',i)
        if key2 not in idx:idx[key2]=len(nodes);nodes.append(['atom',i])
        return idx[key2]
    def mul(i,j):
        if not i or not j:return 0
        key=('mul',i,j)
        if key not in idx:idx[key]=len(nodes);nodes.append(['mul',i,j])
        return idx[key]
    def sub(i,j):
        if not j:return i
        key=('sub',i,j)
        if key not in idx:idx[key]=len(nodes);nodes.append(['sub',i,j])
        return idx[key]
    one=atom([[[[0,0,0],'1']],[[[0,0,0],'1']]])
    P=[[0]*r for _ in range(r)]
    if verified:
        for b in range(r):
            y=[]
            for i in range(r):
                z=one if i==b else 0
                for j in range(i):
                    if L[i][j][0]:z=sub(z,mul(atom(L[i][j]),y[j]))
                y.append(z)
            v=[0]*r
            for i in reversed(range(r)):
                z=y[i]
                for j in range(i+1,r):
                    if U[i][j][0]:z=sub(z,mul(atom(U[i][j]),v[j]))
                v[i]=mul(atom([U[i][i][1],U[i][i][0]]),z)
            for i in range(r):P[i][b]=v[i]
    assert len(nodes)<200000
    fullP=[P[cols.index(i)]for i in range(22)]if full else None
    deleted15=atlas['deleted15_source319_indices_in_pivot_order'];ret304=atlas['retained304_source319_indices']
    remaining=[i for i in ret304 if i not in rows]
    assert len(remaining)==304-r
    certificate=dict(source_K_sha256=sha(kernel),native_input_sha256=sha(inp),native_state_sha256=sha(state),
        allowed_unit_order=['a','d','F','B','C','J'],selected_source174_rows=rows,selected_kernel22_columns=cols,
        selected_original_matrix_A=A,unit_lower_L=L,upper_U=U,native_A_equals_LU_verified=verified,
        selected_minor_determinant='product of recorded known-unit pivots, for rows/columns in selection order',
        scalar_circuit=dict(atoms_num_den=atoms,nodes=nodes),
        P_selected_rows_in_selected_kernel_order=P if verified else None,
        P_columns_original_source174=rows,all_other_source174_columns_zero=True,
        P_full22_rows_in_original_kernel_order=fullP,
        left_inverse_proof='A=L*U and all U diagonal units; P_selected=U^-1*L^-1*pi_rows, so P_selected*K[:,selectedcols]=I_r. Reorder rows bycols for fullPK=I22.',
        original15_deleted_source319=deleted15,additional_selected_quadratic_source319=rows,
        retained_actual_source319_indices=remaining,
        resulting_actual_matrix_shape=[304-r,239],
        global282_actual_submatrix_certified=full,
        residual_source_kernel_rank=22-r,
        no_quartic_coefficient_evaluation=True)
    (out/'certificate.json').write_text(json.dumps(certificate,separators=(',',':'))+'\n')
    result=dict(id=ID,all_completed_checkpoint_checks_passed=verified,native_exit_code=None if cp is None else cp.returncode,
        timed_out=timedout,status=status,known_unit_pivots=r,remaining_kernel_columns=22-r,
        explicit_global_left_inverse_P22x174=full,
        P_is_finite_exact_triangular_solve_circuit_not_expanded_core_coefficients=True,
        actual_global_fixed282x239_same_cokernel=full,
        partial_global_split_rank=r if verified else None,
        working_open='FULLU=D(a*d*F*B*C*J)',
        no_unknown_factor_minor_or_GB_search=True,no_fullquartic_or_C122rank64_claim=True,
        arithmetic_circuit_nodes=len(nodes),scalar_atoms=len(atoms),native_seconds=float(secs),final_Schur_numden_terms=int(terms),
        total_seconds=time.monotonic()-began,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),source,ROOT/'code/s61_binary_hyperplane_cubic_qq3_20261102.hpp',kernel,BASE/'atlas.json',review]},
        evidence_sha256={name:sha(out/name)for name in ['native-input.txt','native-state.txt','native-proof.log','certificate.json']},
        frozen=True,independent_review_pending=True)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    assert(out/'certificate.json').stat().st_size<32*1024*1024
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
