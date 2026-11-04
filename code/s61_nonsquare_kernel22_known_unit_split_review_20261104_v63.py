#!/usr/bin/env python3
"""Independent fixed20 replay, source binding and exact triangular-circuit check."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import resource
import shlex
import signal
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_KERNEL22_KNOWN_UNIT_SPLIT_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(f):
    t=f.readline().split();n=int(t[0]);assert len(t)==1+4*n
    return [[list(map(int,t[2+4*i:5+4*i])),t[1+4*i]] for i in range(n)]
def rat(f):return [poly(f),poly(f)]
def canonical(p):return {tuple(e):Q(c) for e,c in p if Q(c)}
def rateq(a,b):return canonical(a[0])==canonical(b[0]) and canonical(a[1])==canonical(b[1])
def accepted(p):
    r=read(p);assert r.get('all_passed') or r.get('all_completed_checkpoint_checks_passed')
    for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/S61_NONSQUARE_KERNEL22_KNOWN_UNIT_SPLIT_20261104_v63';ap=base/'report.json';author=accepted(ap)
    cp=base/'certificate.json';c=read(cp)
    parent=ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63'
    review=parent.with_name(parent.name+'_REVIEW')/'report.json';pr=accepted(review)
    assert pr['K_basis_same_as_whole_old22_via_unit_exterior_square_and_same_last7']
    kp=parent/'kernel22.json';K=read(kp);atp=parent/'atlas.json';atlas=read(atp)
    inp=base/'native-input.txt';state=base/'native-state.txt';log=base/'native-proof.log'
    assert sha(kp)==c['source_K_sha256'] and sha(inp)==c['native_input_sha256'] and sha(state)==c['native_state_sha256']
    expected={ (i,j):[num,den] for j,col in enumerate(K['sparse_columns']) for i,num,den in col}
    zero=[[],[[[0,0,0],'1']]]
    with inp.open() as f:
        assert f.readline().split()==['174','22']
        for i in range(174):
            for j in range(22):assert rateq(rat(f),expected.get((i,j),zero))
        assert not f.read().strip()
    with state.open() as f:
        assert f.readline().split()==['174','22','20']
        rows=list(map(int,f.readline().split()));cols=list(map(int,f.readline().split()));status=f.readline().split()
        assert rows==c['selected_source174_rows'] and cols==c['selected_kernel22_columns']
        assert status[0]=='no_individual_known_unit_entry' and status[1]=='1'
        matrices=[[[rat(f) for j in range(20)] for i in range(20)] for k in range(3)]
        S=[[rat(f) for j in range(22)] for i in range(174)];assert not f.read().strip()
    A,L,U=matrices
    for actual,stored in zip(matrices,[c['selected_original_matrix_A'],c['unit_lower_L'],c['upper_U']]):
        assert actual==stored
    # Rebuild the exact serialized triangular-inverse recipe; its algebraic
    # identity follows from the separately replayed A=LU and unit diagonals.
    atoms=[];nodes=[['zero']];ai={};ni={('zero',):0}
    def atom(x):
        key=json.dumps(x,separators=(',',':'))
        if key not in ai:ai[key]=len(atoms);atoms.append(x)
        k=('atom',ai[key])
        if k not in ni:ni[k]=len(nodes);nodes.append(list(k))
        return ni[k]
    def op(tag,i,j):
        if tag=='mul' and (not i or not j):return 0
        if tag=='sub' and not j:return i
        key=(tag,i,j)
        if key not in ni:ni[key]=len(nodes);nodes.append(list(key))
        return ni[key]
    one=atom([[[[0,0,0],'1']],[[[0,0,0],'1']]])
    P=[[0]*20 for i in range(20)]
    for b in range(20):
        y=[]
        for i in range(20):
            z=one if i==b else 0
            for j in range(i):
                if L[i][j][0]:z=op('sub',z,op('mul',atom(L[i][j]),y[j]))
            y.append(z)
        x=[0]*20
        for i in range(19,-1,-1):
            z=y[i]
            for j in range(i+1,20):
                if U[i][j][0]:z=op('sub',z,op('mul',atom(U[i][j]),x[j]))
            x[i]=op('mul',atom([U[i][i][1],U[i][i][0]]),z)
        for i in range(20):P[i][b]=x[i]
    assert atoms==c['scalar_circuit']['atoms_num_den'] and nodes==c['scalar_circuit']['nodes'] and P==c['P_selected_rows_in_selected_kernel_order']
    assert len(nodes)==2044 and len(atoms)==73 and c['P_full22_rows_in_original_kernel_order'] is None
    cpp=Path(__file__).with_suffix('.cpp');exe=out/'verify-native'
    flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
    run=subprocess.run(['c++','-O2','-std=c++17',str(cpp),*flags,'-lgmp','-o',str(exe)],capture_output=True,text=True,check=True,timeout=20)
    (out/'compile.log').write_text(run.stdout+run.stderr)
    run=subprocess.run([str(exe),str(inp),str(state),str(log)],capture_output=True,text=True,check=True,timeout=max(1,55-(time.monotonic()-start)))
    nativepath=out/'fixed20-replay.log';nativepath.write_text(run.stdout+run.stderr)
    assert 'FIXED20_UNIT_LOG_ALL3828_SCHUR_ENTRIES_AND_A_L_U_PASS' in run.stdout
    remainingrows=[i for i in range(174) if i not in rows];remainingcols=[j for j in range(22) if j not in cols]
    assert remainingcols==[12,21] and len(remainingrows)==154
    frontierpath=base/'residual2-frontier.json';fr=read(frontierpath)
    assert fr['parent_certificate_sha256']==sha(cp) and fr['original_source174_rows']==remainingrows and fr['original_kernel22_columns']==remainingcols
    residual=[[S[i][j] for j in remainingcols] for i in remainingrows]
    assert residual==fr['entries_num_den']
    active=[row for row in residual if any(c[0] for c in row)];nonzero=[c for row in residual for c in row if c[0]]
    assert len(active)==106 and len(nonzero)==195 and sum(len(c[0])+len(c[1]) for c in nonzero)==7763
    retained=[i for i in atlas['retained304_source319_indices'] if i not in rows]
    assert retained==c['retained_actual_source319_indices'] and c['resulting_actual_matrix_shape']==[284,239] and len(retained)==284
    assert c['original15_deleted_source319']==atlas['deleted15_source319_indices_in_pivot_order']
    assert c['additional_selected_quadratic_source319']==rows and c['residual_source_kernel_rank']==2
    assert not c['global282_actual_submatrix_certified'] and not author['explicit_global_left_inverse_P22x174']
    evidence=dict(all_input_K_entries_bound=True,selected_rows=rows,selected_columns=cols,
        all3828_Schur_coefficients_and_A_L_U_replayed=True,scalar_DAG_nodes=2044,scalar_atoms=73,
        P20_triangular_inverse_identity_not_expanded=True,source284_indices=retained,
        residual_shape=[154,2],residual_columns=remainingcols,residual_active_rows=106,residual_nonzero_entries=195,residual_numden_terms=7763,
        no_additional_candidates_evaluated=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),cpp,ROOT/f'research/notes/{STEM}.md',ap,cp,inp,state,log,frontierpath,review,kp,atp,
        ROOT/'code/s61_binary_hyperplane_cubic_qq3_20261102.hpp',ROOT/'research/notes/S61_NONSQUARE_KERNEL22_KNOWN_UNIT_SPLIT_20261104_v63.md']
    report=dict(id=STEM,all_passed=True,independent_fixed20_replay_no_pivot_selection=True,
        arithmetic_engine='Shared accepted FLINT QQ(a,e,d) utility; independent fixed-step replay and serialized-circuit validator',
        exact_all_K_input_Schur_A_LU_and_known_unit_factor_log_verified=True,
        exact_P20_partial_left_split_2044nodes73atoms_verified=True,actual_fixed_M284_same_image_cokernel=True,
        residual154_by2_pointwise_rank2_inherited_not_new_rank_search=True,
        fullP22_global282_submatrix_or_full239_claim=False,
        no_additional_pivot_minor_Bezout_GB_point_field_or_quartic_expansion=True,
        existing_Cloud_activation_reused_after_unactivated_pkgconfig_failure_before_native_arithmetic=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep),nativepath.name:sha(nativepath)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','actual_fixed_M284_same_image_cokernel','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
