#!/usr/bin/env python3
"""One bounded exact probe of448 residual cubic cofactor coefficients only."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from s61_nonsquare_lazy_quartic_oracle_20261103_v63 import Oracle

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63'
BASE=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63'

def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    review_path=ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW/report.json'
    review=read(review_path);assert review['all_passed']and review['all15_exact354_vectors_and35unit_deletions_match_all319_coefficients']
    parent=read(BASE/'report.json');assert parent['all_passed']and parent['known_unit_pivot_rank']==11
    for name,h in parent['evidence_sha256'].items():assert sha(BASE/name)==h
    o,C=Oracle.load(BASE/'certificate.json');K=o.K;initial_nodes=len(o.nodes)
    P=read(BASE/'block-pivots.json');assert P['rank_by_global_known_unit_pivots']==11
    frac=lambda c:K.from_expr(s.sympify(c))
    U=[[frac(c)for c in row]for row in P['elementary_transform']]
    rows=C['actual319_source_definitions'];assert len(rows)==319
    cubic=[i for i,row in enumerate(rows)if row[0]=='cubic_normalized'];quartic=[i for i,row in enumerate(rows)if row[0]=='quartic']
    assert len(cubic)==112 and len(quartic)==33
    sources=[r['after35_deletions_coefficients319']for r in C['raw15_source_lifts']]
    assert len(sources)==15 and all(len(v)==319 for v in sources)
    residualU=U[11:];assert len(residualU)==4
    basis_sources=[]
    for row in residualU:
        indices=[j for j,c in enumerate(row)if c]
        assert len(indices)==1 and row[indices[0]]==K.one
        basis_sources.append(indices[0])
    assert basis_sources==[0,4,13,14]
    # Shared exact cache expands only the specifically selected cofactor roots.
    # No pure quadratic/light cofactor coefficient is evaluated here.
    cache={0:K.zero}
    def value(i):
        if i in cache:return cache[i]
        assert o.weights[i]==0
        nd=o.nodes[i]
        if nd[0]=='c':z=o.atoms[nd[1]]
        elif nd[0]=='m':z=value(nd[1])*value(nd[2])
        else:
            assert nd[0]=='l';z=sum((o.atoms[c]*value(n)for c,n in nd[1]),K.zero)
        cache[i]=z;return z
    def combined(row,k):return o.lin([(c,sources[j][k])for j,c in enumerate(row)if c])
    quartic_nodes=[];cubic_nodes=[];matrix=[]
    for row in residualU:
        qn=[combined(row,k)for k in quartic]
        assert all(not value(n)for n in qn)
        cn=[combined(row,k)for k in cubic]
        mat=[value(n)for n in cn]
        quartic_nodes.append(qn);cubic_nodes.append(cn);matrix.append(mat)
    assert len(matrix)*len(matrix[0])==448
    (out/'progress.json').write_text(json.dumps(dict(stage='all132quarticzero_and448cubic_coefficients_expanded',seconds=time.monotonic()-start,
        nonzero_cubic_coefficients=sum(bool(c)for row in matrix for c in row)),indent=2)+'\n')
    # Only numerators and denominators supported on a,F,B,C may be inverted.
    def unit(c):
        c=K(c)
        if not c:return None
        powers=[];ground=[]
        for p in[c.numer,c.denom]:
            row=[]
            for q in o.known:
                n=0
                while not p.rem(q):p=p.exquo(q);n+=1
                row.append(n)
            if not(p.is_ground and p):return None
            powers.append(row);ground.append(p)
        return dict(constant=str(K(ground[0])/K(ground[1])),a_F_B_C_powers=[a-b for a,b in zip(*powers)])
    B=[row[:]for row in matrix];T=[[K.one if i==j else K.zero for j in range(4)]for i in range(4)]
    pivots=[];used=[];rank=0
    while rank<4:
        candidates=[]
        for i in range(rank,4):
            for j,c in enumerate(B[i]):
                if j in used:continue
                proof=unit(c)
                if proof is not None:candidates.append((len(c.numer)+len(c.denom),i,j,proof))
        if not candidates:break
        _,i,j,proof=min(candidates,key=lambda x:x[:3]);p=B[i][j]
        if i!=rank:B[i],B[rank]=B[rank],B[i];T[i],T[rank]=T[rank],T[i]
        B[rank]=[c/p for c in B[rank]];T[rank]=[c/p for c in T[rank]]
        operations=[]
        for k in range(4):
            if k==rank or not B[k][j]:continue
            c=B[k][j];operations.append([k,str(c)])
            B[k]=[a-c*b for a,b in zip(B[k],B[rank])];T[k]=[a-c*b for a,b in zip(T[k],T[rank])]
        pivots.append(dict(row=rank,swap_row=i,cubic_column=j,source319_column=cubic[j],pivot=str(p),unit=proof,eliminations=operations))
        used.append(j);rank+=1
    assert [[sum((T[i][k]*matrix[k][j]for k in range(4)),K.zero)for j in range(112)]for i in range(4)]==B
    allzero=all(not c for row in matrix for c in row)
    assert len(o.nodes)<200000
    record=dict(parent_certificate_sha256=sha(BASE/'certificate.json'),residual_transform_rows15to4=[[str(c)for c in row]for row in residualU],
        residual_original_block_rows=basis_sources,residual_parent_piece_indices=[C['raw15_source_lifts'][i]['source_piece_index']for i in basis_sources],
        cubic_source319_indices=cubic,quartic_source319_indices=quartic,
        cubic_coefficient_nodes=cubic_nodes,quartic_coefficient_nodes=quartic_nodes,
        all448_exact_cubic_values=[[str(c)for c in row]for row in matrix],
        unit_pivot_steps=pivots,pivot_row_transform=[[str(c)for c in row]for row in T],
        post_pivot_cubic_matrix=[[str(c)for c in row]for row in B],
        all_core_values_evaluated_at_actual_selected_roots=True,
        pure_quadratic_cofactors_not_expanded=True,
        parent_DAG_node_count=initial_nodes,additional_DAG_nodes=o.nodes[initial_nodes:])
    (out/'certificate.json').write_text(json.dumps(record,indent=2)+'\n')
    report=dict(id=ID,all_passed=True,base_open_new_arithmetic='a*F*B*C!=0',
        exact_expanded_cubic_coefficient_count=448,exact_quartic_zero_checks=132,
        all448_cubic_coefficients_literal_zero=allzero,nonzero_cubic_coefficient_count=sum(bool(c)for row in matrix for c in row),
        new_known_unit_directions=rank,remaining_block_exact_zero=all(not c for row in B[rank:]for c in row),
        direct_operator_exhaustion_on_FULLU_if_cubic_zero=allzero,
        exhaustion_scope='Only the15 DIRECT weight4 pieces; not all weighted Bianchi multiples, all syzygies or the80gate',
        theoretical_Hq_coefficient_zero_if_cubic_zero='Coefficient of H is a constant linear relation among the six pointwise-independent quadrics, so all six coefficients vanish on FULLU',
        theoretical_remaining_pure_quadratic_membership='The accepted whole22 quadratic-cofactor kernel module on FULLU contains every residual pure-quadric relation',
        old22_plus_direct15_rank_on_FULLU=33+rank if allzero else None,
        no_pure_quadratic_H_coefficient_or239_column_expansion=True,
        no_full239_C122rank64_80kernel_or_affine_exhaustion_claim=True,
        no_new_point_field_GB_or_automatic_enlargement=True,
        evaluated_core_DAG_nodes=len(cache),initial_DAG_nodes=initial_nodes,final_DAG_nodes=len(o.nodes),
        seconds=time.monotonic()-start,frozen=True,independent_review_pending=True)
    paths=[Path(__file__),BASE/'report.json',BASE/'certificate.json',BASE/'block-pivots.json',review_path,
        ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/report.json',
        ROOT/'research/notes/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63.md',
        ROOT/'research/notes/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63.md']
    report['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths}
    report['certificate_sha256']=sha(out/'certificate.json')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
