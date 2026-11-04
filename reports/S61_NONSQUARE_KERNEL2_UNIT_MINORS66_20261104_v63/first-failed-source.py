#!/usr/bin/env python3
"""One exact known-unit check on the approved66 residual2 row pairs only."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_KERNEL2_UNIT_MINORS66_20261104_v63'
BASE=ROOT/'reports/S61_NONSQUARE_KERNEL22_KNOWN_UNIT_SPLIT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main(out):
    began=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    frontier_path=BASE/'residual2-frontier.json';pub=read(BASE/'publish-receipt.json')
    assert sha(frontier_path)==pub['residual2_frontier_sha256']
    F=read(frontier_path);parent_path=BASE/'certificate.json';C=read(parent_path)
    assert sha(parent_path)==F['parent_certificate_sha256']
    assert C['native_A_equals_LU_verified']and len(C['selected_source174_rows'])==20
    assert C['source_K_sha256']==sha(ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json')
    pairs=F['proposed_only_first12_low_term_rows_pair_candidates'];assert len(pairs)==66
    rows=F['original_source174_rows'];lookup=dict(zip(rows,F['entries_num_den']))
    R,a,e,d=ring('a,e,d',s.QQ)
    known=[a,d,500*a+90*e**2+3*e-24,(50*a+18*e**2-9*e-2)*d+60*e-40,100*a+18*e**2-9*e-2,5*a*d-6*e-4]
    dec=lambda p:R.from_dict({tuple(ex):s.QQ(c)for ex,c in p})
    enc=lambda p:[[list(ex),str(c)]for ex,c in p.items()]
    def strip(p):
        if not p:return[0]*6,p
        powers=[]
        for k in known:
            n=0
            while not p.rem(k):p=p.exquo(k);n+=1
            powers.append(n)
        return powers,p
    def rat(c):return dec(c[0]),dec(c[1])
    outcomes=[];found=None
    for number,record in enumerate(pairs):
        i,j=record['original_source_rows'];A=[list(map(rat,lookup[x]))for x in[i,j]]
        for row in A:
            for n,q in row:
                ex,left=strip(q);assert left.is_ground and left
        (n00,d00),(n01,d01)=A[0];(n10,d10),(n11,d11)=A[1]
        # Unreduced but literal determinant numerator/denominator; denominator
        # already a known unit, so NO arbitrary multivariate gcd/factor needed.
        numerator=n00*n11*d01*d10-n01*n10*d00*d11
        denominator=d00*d11*d01*d10
        ne,nrem=strip(numerator);de,drem=strip(denominator)
        assert drem.is_ground and drem
        isunit=bool(nrem and nrem.is_ground)
        result=dict(candidate_index=number,original_source_rows=[i,j],
            matrix2_num_den=[lookup[i],lookup[j]],determinant_numerator=enc(numerator),determinant_denominator=enc(denominator),
            known_unit_factor_order=['a','d','F','B','C','J'],numerator_known_factor_powers=ne,denominator_known_factor_powers=de,
            remaining_numerator_after_known_factor_division=enc(nrem),remaining_denominator=enc(drem),
            determinant_known_unit=isunit,
            numerator_degree=numerator.total_degree()if numerator else-1,
            numerator_nonzero_terms=len(numerator),denominator_nonzero_terms=len(denominator))
        outcomes.append(result)
        (out/'progress.json').write_text(json.dumps(dict(tested=len(outcomes),unit_found=isunit,seconds=time.monotonic()-began))+'\n')
        if isunit:found=result;break
    packet=dict(parent_certificate_sha256=sha(parent_path),residual2_frontier_sha256=sha(frontier_path),
        approved_candidate_count=66,tested_count=len(outcomes),outcomes=outcomes,
        first_known_unit_minor=found,all66_failed_known_unit_test=(found is None and len(outcomes)==66),
        no_claim_of_global_nonexistence_or_no_Bezout_split=True)
    (out/'outcomes.json').write_text(json.dumps(packet,separators=(',',':'))+'\n')
    full=False;node_count=0
    if found:
        # Compose the accepted partial20 split with the2-row Schur inverse.
        # All arithmetic stays a finite scalar circuit; no P22 core expansion.
        old=C['scalar_circuit'];atoms=list(old['atoms_num_den']);nodes=list(old['nodes'])
        ai={json.dumps(v,separators=(',',':')):i for i,v in enumerate(atoms)}
        ni={tuple(n):i for i,n in enumerate(nodes)}
        def at(v):
            key=json.dumps(v,separators=(',',':'))
            if key not in ai:ai[key]=len(atoms);atoms.append(v)
            k=('atom',ai[key])
            if k not in ni:ni[k]=len(nodes);nodes.append(list(k))
            return ni[k]
        def op(t,i,j):
            if t in['add','sub']and not j:return i
            if t=='add'and not i:return j
            if t=='mul'and(not i or not j):return 0
            k=(t,i,j)
            if k not in ni:ni[k]=len(nodes);nodes.append(list(k))
            return ni[k]
        one=at([[[[0,0,0],'1']],[[[0,0,0],'1']]])
        p20=C['P_selected_rows_in_selected_kernel_order'];sigma=C['selected_source174_rows'];tau=C['selected_kernel22_columns'];J=[12,21]
        extra=found['original_source_rows'];sigma22=sigma+extra;assert len(set(sigma22))==22
        kd=read(ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json')
        K={(i,j):[num,den]for j,col in enumerate(kd['sparse_columns'])for i,num,den in col}
        def kn(i,j):return at(K[i,j])if(i,j)in K else 0
        P20=[row+[0,0]for row in p20]
        c20=[]
        for i in range(20):
            row=[]
            for j in J:
                z=0
                for b in range(20):z=op('add',z,op('mul',p20[i][b],kn(sigma[b],j)))
                row.append(z)
            c20.append(row)
        Q=[]
        for t,rowno in enumerate(extra):
            row=[]
            for b in range(22):
                z=one if b==20+t else 0
                for k in range(20):z=op('sub',z,op('mul',kn(rowno,tau[k]),P20[k][b]))
                row.append(z)
            Q.append(row)
        S=found['matrix2_num_den'];invdet=at([found['determinant_denominator'],found['determinant_numerator']])
        adj=[[at(S[1][1]),op('sub',0,at(S[0][1]))],[op('sub',0,at(S[1][0])),at(S[0][0])]]
        G=[]
        for i in range(2):
            row=[]
            for b in range(22):
                z=0
                for t in range(2):z=op('add',z,op('mul',adj[i][t],Q[t][b]))
                row.append(op('mul',invdet,z))
            G.append(row)
        ordered=[]
        for i in range(20):
            row=[]
            for b in range(22):
                z=P20[i][b]
                for j in range(2):z=op('sub',z,op('mul',c20[i][j],G[j][b]))
                row.append(z)
            ordered.append(row)
        ordered+=G;colorder=tau+J
        P22=[ordered[colorder.index(i)]for i in range(22)]
        retained=[i for i in C['retained_actual_source319_indices']if i not in extra];assert len(retained)==282
        proof=dict(scalar_circuit=dict(atoms_num_den=atoms,nodes=nodes),P22_roots_in_original_kernel_column_order=P22,
            P22_nonzero_source174_columns=sigma22,all_other_columns_zero=True,
            block_formula='Q=pi_extra*(I-K_I*P20); G=S^-1*Q; Pfull_orderedI_J=(P20-(P20*K_J)*G,G)',
            exact_proof='P20*K_I=I20; Q*K_I=0; Q*K_J=S from exact Schur state. detS is known unit. Thus G*K_I=0,G*K_J=I2 and Pfull*K=I22.',
            selected22_minor_determinant='accepted20-unit determinant times the certified residual2 determinant (column order tau+J)',
            retained_actual_source319_indices=retained,actual_global_matrix_shape=[282,239],
            no_quartic_coefficients_evaluated=True)
        assert len(nodes)<200000;(out/'global-P22.json').write_text(json.dumps(proof,separators=(',',':'))+'\n')
        node_count=len(nodes);full=True
    result=dict(id=ID,all_exact_determinant_and_known_factor_checks_passed=True,
        approved66pairs_only=True,tested_pairs=len(outcomes),known_unit_pair_found=found is not None,
        global_P22_finite_circuit_constructed=full,actual_global_fixed282_same_cokernel=full,
        no_automatic_pair_list_enlargement_Bezout_GB_field_point_or_unknown_inversion=True,
        no_full239_C122rank64_or_new_kernel_direction_claim=True,
        arithmetic_circuit_nodes=node_count,seconds=time.monotonic()-began,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),parent_path,BASE/'report.json',frontier_path,BASE/'publish-receipt.json',ROOT/'reports/S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63/kernel22.json']},
        evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()and p.name!='report.json'},
        frozen=True,independent_review_pending=True)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    for p in out.iterdir():assert p.stat().st_size<32*1024*1024
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
