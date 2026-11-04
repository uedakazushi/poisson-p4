#!/usr/bin/env python3
"""Pure module-logic review and frozen support binding, no coefficient/rank job."""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_BIANCHI_WEIGHT4_CONTAINMENT_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
    r=read(p);assert r['all_passed']
    for n,h in r.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r

def main(out):
    start=time.monotonic();assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/S61_NONSQUARE_BIANCHI_WEIGHT4_CONTAINMENT_20261104_v63'
    arpath=base/'report.json';author=accepted(arpath);acpath=base/'certificate.json';ac=read(acpath)
    assert author['certificate_sha256']==sha(acpath)
    reviewnames=['S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW',
        'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW',
        'S61_NONSQUARE_CUBIC_QUARTIC_STRUCTURAL_REDUCTIONS_20261103_v63_REVIEW',
        'S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW',
        'S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW',
        'S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63_REVIEW']
    reviews=[ROOT/f'reports/{n}/report.json' for n in reviewnames]
    op,cb,st,q22,b11,b4=[accepted(p) for p in reviews]
    assert cb['full_localized_row_kernel_generated_by_six_relations']
    assert st['quartic354_to319_row_span_equal_on_aFBC'] and st['all35_lifted_rational_row_deletions_independently_checked']
    assert q22['entire_kernel_module_basis_accepted']
    assert not q22['ambient_specialized_rank146_or_cokernel_flatness_proved']
    assert b11['all15_exact354_vectors_and35unit_deletions_match_all319_coefficients']
    assert b4['new_direct15_independent_everywhere_on_aFBC'] and b4['combined_with_old22_only_on_FULLU']==37
    dp=ROOT/'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json';data=read(dp);rels=data['relations']
    deleted=[r['deleted_original_row'] for r in rels];assert deleted==[48,55,59,57,61,14]
    previous=set();piv=[]
    for r in rels:
        cofs={v['original_row']:v['terms'] for v in r['coefficients']}
        assert not (previous & cofs.keys())
        p=cofs[r['deleted_original_row']];assert len(p)==1
        term=p[0];assert term[2:]==[0,0];piv.append((Fraction(term[0]),term[1]));previous.add(r['deleted_original_row'])
    assert piv==[(Fraction(1),0)]*3+[(Fraction(1),1),(Fraction(-1),1),(Fraction(12,5),1)]
    assert sum(n for c,n in piv)==3
    det=Fraction(1)
    for c,n in piv:det*=c
    assert det==Fraction(-12,5)
    assert all(v['original_row']<42 for v in rels[5]['coefficients'])
    sourcepath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json';source=read(sourcepath)
    assert all(v[0]=='product' for v in source['labels'][:42]) and all(v[0]=='cubic' for v in source['labels'][42:])
    oppath=ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63/operator.json';O=read(oppath)
    oraclepath=ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json';A=read(oraclepath)
    framepath=ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json';frame=read(framepath)
    selected=set(A['selected34_rows']);cw=O['local_coordinate_weights'];rw=O['equation_weights']
    assert cw==A['local_coordinate_weights'] and rw==A['equation_weights']
    for j in [19,23]:
        assert all(A[k][j]==0 for k in ['linear_graph_roots','local_Z2_roots','local_Z3_roots','local_Z4_roots'])
    qrows=set(frame['quadric_original_rows']);genuine={v[1] for v in source['labels'][42:]}
    low=[];direct=[]
    for i,p in enumerate(O['local_weight_pieces']):
        assert p['weight']>=3
        if p['weight']==4:direct.append(i)
        if p['weight']!=3:continue
        assert p['spatial_degree']==5
        qsupport=set();csupport=set()
        for g,j,c in p['coefficients']:
            if g in selected or j in [19,23]:continue
            if rw[g]==2:
                assert g in qrows and cw[j]==1;qsupport.add(g)
            else:
                assert rw[g]==3 and g in genuine and j==-1;csupport.add(g)
        low.append(dict(piece_index=i,operator_row=p['operator_row'],spatial_monomial=p['spatial_monomial'],
                        original_quadric_support=sorted(qsupport),original_genuine_cubic_support=sorted(csupport)))
    assert [v['piece_index'] for v in low]==[35,56,71,81,87,90]
    assert low==ac['six_weight3_pieces'] and direct==ac['direct15_weight4_piece_indices']
    assert len(direct)==15 and direct==O['weight4_quartic_constant_block_row_piece_indices']
    rowsfile=ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json';rows=read(rowsfile)
    light=[0,1,2,3,4,6,7];mults=[]
    for r in rows['lifted_relations']:
        dr=rows['source_rows'][r['deleted_row']]
        assert dr[0]=='cubic_normalized';mults.append((dr[1],dr[3]))
        for j,c in r['coefficients']:
            row=rows['source_rows'][j]
            if row[0]=='quadric':assert row[2][5]==0
            else:assert row[0]=='cubic_normalized'
    assert mults==[(g,j) for g in deleted[:5] for j in light]
    evidence=dict(six_pivot_monomials=[[str(c),n] for c,n in piv],triangular_minor_coefficient=str(det),triangular_minor_a_power=3,
        operator_piece_weight_histogram=dict(Counter(p['weight'] for p in O['local_weight_pieces'])),
        exact6_weight3_pieces=low,direct15_piece_indices=direct,exact35_deleted_lift_labels=mults,
        no_piece_below3=True,whole_cubic_kernel_gate_reused=True,whole_purequad_kernel22_gate_reused=True,
        no_pointwise_57_or146_or_cokernel_flatness_inference=True)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',arpath,acpath,*reviews,dp,sourcepath,oppath,oraclepath,framepath,rowsfile]
    result=dict(id=STEM,all_passed=True,review_kind='Independent pure module proof with frozen source/support/hash binding',
        whole_cubic_kernel6_by_unit_quotient_generic_rank_torsionfree=True,
        six_weight3_pieces_in_same_actual63_source_kernel=True,
        first5_times7_killed_by_35deletions_and_sixth_times_light_in_whole22=True,
        direct15_subset_B4_subset_old22_plus_direct15_on_FULLU=True,
        old22_plus_B4_equals_direct_sum_old22_direct15_free_rank37_on_FULLU=True,
        Bianchi_B4_alone_rank37_claimed=False,missing_for80_gate=43,
        higher_weight_division_saturation_other_syzygy_or_full80_claimed=False,
        no_new_field_point_coefficient_expansion_pivot_rank_GB_or_239_job=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','old22_plus_B4_equals_direct_sum_old22_direct15_free_rank37_on_FULLU','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
