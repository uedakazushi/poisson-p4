#!/usr/bin/env python3
"""Source-bound25x32 input and fixed8 exact unit-operation replay; no search."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring
from s61_nonsquare_free281_circuit_review_20261104_v63 import read,sha
from s61_nonsquare_joint_Lell_free183_review_20261104_v63 import matrix

ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63'
ID=BASE+'_REVIEW'

def main(out):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('independent fixed8 replay60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    folder=ROOT/'reports'/BASE;pr=read(folder/'report.json');state=read(folder/'checkpoint.json')
    inp=ROOT/'reports/S61_NONSQUARE_CUBIC_ZW25_AUGMENTED_INPUT_20261104_v63';ir=read(inp/'report.json');data=read(inp/'input.json')
    assert pr['frozen'] and state['completed_pivots']==pr['completed_global_known_unit_pivots']==8
    for rp in [folder/'report.json',inp/'report.json']:
        for n,h in read(rp)['input_sha256'].items():assert sha(ROOT/n)==h,n
    for n,h in pr['evidence_sha256'].items():assert sha(folder/n)==h,n
    assert sha(inp/'input.json')==ir['input_json_sha256']==state['source_input_sha256']
    R,a,e,d=ring('a,e,d',s.QQ);K=R.to_field();zero=K.zero
    def poly(ts):return R.from_dict({tuple(m):s.QQ(c) for m,c in ts})
    def dec(z):return K(poly(z[0]))/K(poly(z[1]))
    C=100*a+18*e**2-9*e-2;F=500*a+90*e**2+3*e-24;B=(50*a+18*e**2-9*e-2)*d+60*e-40;J=5*a*d-6*e-4
    known=[a,d,F,B,C,J]
    pre=ROOT/'reports/S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63';pn=read(pre/'native.json');pc=read(pre/'certificate.json')
    Q10=matrix(pre/'quad10x26.txt');G31=matrix(pre/'genuine16x31.txt');G40=matrix(pre/'genuine16x40.txt')
    assert G31[6]==G31[1] and G40[6][5:]==G40[1][5:] and G40[6][:5]!=G40[1][:5]
    gates=[ROOT/'reports'/n/'report.json' for n in [
        'S61_NONSQUARE_ELL_LFRAME_PREFLIGHT_20261104_v63_REVIEW',
        'S61_NONSQUARE_ELL_COLON_RANK10_20261104_v63_REVIEW',
        'S61_NONSQUARE_ELL_LFRAME_BORDERS_20261104_v63_REVIEW',
        'S61_NONSQUARE_CUBIC_PROJECTED15_UNIT_PAIR_20261104_v63_REVIEW']]
    assert all(read(p)['all_passed'] for p in gates)
    br=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json')
    src=read(ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json');mons=src['normal_monomials']
    raw=ROOT/'reports/S61_NONSQUARE_ELL_SIX_BORDERS_20261104_v63/frozen57-source.txt';assert sha(raw)==br['remaining_matrix_sha256']
    hcols=[]
    for j in [6,3]:m=[0]*8;m[5]=m[j]=1;hcols.append(mons.index(m))
    H=[[R.zero,R.zero] for _ in range(16)]
    with raw.open() as f:
        assert f.readline().split()==['57','91','57']
        for i in range(57):
            for j in range(91):
                line=f.readline();assert line
                if i<41 or j not in hcols:continue
                t=line.split();assert len(t)==1+4*int(t[0])
                H[i-41][hcols.index(j)]=poly([[list(map(int,t[k+1:k+4])),t[k]] for k in range(1,len(t),4)])
        assert not f.read().strip()
    assert H==[[poly(v) for v in row] for row in data['raw16_HZ_HW_coefficients']]
    delta=-64000*C*C;z=delta*(H[6][0]-H[1][0]);w=delta*(H[6][1]-H[1][1])
    units=[-64000*C**(3-k) for k in pn['genuine_C_factors_cancelled']]
    assert units[1]==units[6]==delta
    assert units==list(map(poly,data['genuine16_source_clearing_units']))
    psi=[unit*(-w*hz+z*hw) for unit,(hz,hw) in zip(units,H)]
    assert psi==list(map(poly,data['psi16_coefficients'])) and psi[6]==psi[1]
    # Original frozen frame conventions are reused from the accepted L-frame audit.
    keep=[i for i in range(16) if i!=6]
    augmented=[row+[[]]*6 for row in Q10]+[G31[i]+[[[list(m),str(c)] for m,c in psi[i].items()]] for i in keep]
    assert [[poly(v) for v in row] for row in augmented]==[[poly(v) for v in row] for row in data['matrix25x32']]
    assert len(augmented)==25 and all(len(r)==32 for r in augmented)
    assert data['rawHZ_HW_indices_in_frozen91']==hcols
    assert data['target31_recipe']['NO_heavy5_Schur_has_been_applied']
    assert data['target31_recipe']['retained31_full40_columns']==pc['remaining26_full35_column_indices']+list(range(35,40))
    expectedlabels=[['ZW_quad10',i,pc['quad10_rows'][i]] for i in range(10)]+[['ZW_genuine',i,41+i,br['remaining_original_row_indices'][41+i]] for i in keep]
    assert expectedlabels==data['source25_labels']
    # Each step below is already selected and saved. No unit scan is repeated.
    A=[[K(poly(v)) for v in row] for row in augmented];active_rows=list(range(25));active_cols=list(range(32))
    determinant=K.one;sum_powers=[0]*6;constant=Q(1);eliminations=0
    for step in state['all_completed_unit_steps']:
        i,j=step['pivot_row'],step['pivot_column'];assert i in active_rows and j in active_cols
        pivot=dec(step['pivot']);assert A[i][j]==pivot and pivot
        unit=step['known_unit'];certified=K.from_expr(s.sympify(unit['constant']))
        assert certified
        for f,n in zip(known,unit['a_d_F_B_C_J_powers']):certified*=K(f)**n
        assert certified==pivot
        determinant*=pivot;constant*=Q(unit['constant']);sum_powers=[x+y for x,y in zip(sum_powers,unit['a_d_F_B_C_J_powers'])]
        normalized=[v/pivot for v in A[i]];assert normalized[j]==1
        saved=dict(step['eliminations']);expected_nonzero=[t for t in active_rows if t!=i and A[t][j]]
        assert list(saved)==expected_nonzero
        A[i]=normalized
        for t in expected_nonzero:
            co=dec(saved[t]);assert co==A[t][j];eliminations+=1
            for k in active_cols:
                A[t][k]=zero if k==j else A[t][k]-co*normalized[k]
            assert not A[t][j]
        active_rows.remove(i);active_cols.remove(j)
    assert len(state['all_completed_unit_steps'])==8
    assert active_rows==state['remaining_rows'] and active_cols==state['remaining_columns']
    target=[[dec(v) for v in row] for row in state['matrix25x32']]
    assert A==target
    saved_product=state['selected_minor_unit_determinant'];assert Q(saved_product['constant'])==constant and saved_product['a_d_F_B_C_J_powers']==sum_powers
    assert determinant==K(-1)/4608*K(a)**4*K(F)*K(B)**2*K(C)**2
    assert [len(active_rows),len(active_cols)]==[17,24] and not pr['global_augmented25x32_rank25_proved'] and not pr['full_cubic57_all_FULLU_proved_primary_pending_review']
    assert not pr['unfinished_step_not_in_checkpoint'] and pr['unfinished_step_coordinates'] is None
    assert state['selected_minor_rows']==[t['pivot_row'] for t in state['all_completed_unit_steps']]
    assert state['selected_minor_columns']==[t['pivot_column'] for t in state['all_completed_unit_steps']]
    out.mkdir(parents=True,exist_ok=False);ep=out/'independent-proof.json'
    evidence={'input25by32_all60380_source_coefficients_bound':True,
        'actual16_HZ_HW_and_source_clearing_units_and_psi_verified':True,
        'duplicate_G6_G1_in_pi_and_psi_only_not_full_target':True,
        'all8_fixed_pivots_and_eliminations_verified':True,'elimination_rows':eliminations,
        'all_final25by32_rational_entries_exactly_match':True,
        'selected_unit_minor_determinant':'-a^4*F*B^2*C^2/4608',
        'remaining_shape':[17,24],
        'no_remaining_unit_scan_or_new_pivot_search':True,
        'algorithm_stall_has_no_negative_geometric_inference':True,
        'no_full25_or_fullFULLU_cubic57_or_quartic239_claim':True}
    ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),folder/'report.json',folder/'checkpoint.json',inp/'report.json',inp/'input.json',
        pre/'native.json',pre/'certificate.json',pre/'quad10x26.txt',pre/'genuine16x31.txt',pre/'genuine16x40.txt',raw,*gates,
        ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json',ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json',
        ROOT/'code/s61_nonsquare_joint_Lell_free183_review_20261104_v63.py']
    result={'id':ID,'all_passed':True,'actual_augmented25x32_input_source_and_frame_bound':True,
        'all8_known_unit_fixed_operations_and_all_checkpoint_entries_verified':True,
        'global_unit8_minor_and_remaining17by24_frontier_accepted':True,
        'global_augmented25_or_fullFULLU_cubic57_or_quartic239_claim':False,
        'no_authorsolver_Gaussian_search_newpoint_prime_or_GB':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
        'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
