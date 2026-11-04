#!/usr/bin/env python3
"""Bind a pure module argument; no coefficient search or new rank calculation."""
import argparse
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
ID='S61_NONSQUARE_BIANCHI_WEIGHT4_CONTAINMENT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    paths={
      'operator':ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63/operator.json',
      'operator_review':ROOT/'reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW/report.json',
      'cubic_basis':ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json',
      'cubic_review':ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW/report.json',
      'cubic_relations':ROOT/'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json',
      'source63':ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/source-report.json',
      'oracle':ROOT/'reports/S61_NONSQUARE_LAZY_QUARTIC_ORACLE_20261103_v63/oracle.json',
      'frame':ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63/adapter.json',
      'adapter_review':ROOT/'reports/S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63_REVIEW/report.json',
      'source319':ROOT/'reports/S61_NONSQUARE_QUARTIC_319_ROWS_20261103_v63/certificate.json',
      'source319_review':ROOT/'reports/S61_NONSQUARE_CUBIC_QUARTIC_STRUCTURAL_REDUCTIONS_20261103_v63_REVIEW/report.json',
      'quartic11':ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/report.json',
      'quartic11_certificate':ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63/certificate.json',
      'quartic11_review':ROOT/'reports/S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63_REVIEW/report.json',
      'cubic4':ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/report.json',
      'cubic4_scope':ROOT/'reports/S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63/scope-receipt.json',
      'quad22':ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63/report.json',
      'quad22_review':ROOT/'reports/S61_NONSQUARE_UNIVERSAL_QUAD22_20261103_v63_REVIEW/report.json'}
    D={k:read(p)for k,p in paths.items()}
    assert all(D[k]['all_passed']for k in ['operator_review','cubic_basis','cubic_review','adapter_review','source319_review','quartic11','quartic11_review','cubic4'])
    assert D['quad22']['exact_kernel_module_basis_uses_reviewed_global_upper_and_generic_rank146']
    assert D['cubic_review']['full_localized_row_kernel_generated_by_six_relations']
    assert D['cubic_basis']['global_free_row_module_rank']==57
    assert D['cubic_basis']['generic_rank_independently_verified']==57
    assert not D['cubic_basis']['all_point_rank_lower57_claimed']
    rels=D['cubic_relations']['relations'];assert len(rels)==6
    deleted=[r['deleted_original_row']for r in rels];assert deleted==[48,55,59,57,61,14]
    pivots=[];previous=set()
    for r in rels:
        coeffs={c['original_row']:c['terms']for c in r['coefficients']}
        assert not previous.intersection(coeffs)
        pivots.append(coeffs[r['deleted_original_row']]);previous.add(r['deleted_original_row'])
    assert pivots==[[['1',0,0,0]],[['1',0,0,0]],[['1',0,0,0]],
                    [['1',1,0,0]],[['-1',1,0,0]],[['12/5',1,0,0]]]
    assert D['cubic_basis']['primary_triangular_minor']=='-12*a^3/5'
    assert all(c['original_row']<42 for c in rels[5]['coefficients'])
    assert all(D['source63']['labels'][i][0]=='product'for i in range(42))
    assert all(D['source63']['labels'][i][0]=='cubic'for i in range(42,63))
    assert all(D['source63']['labels'][c['original_row']][0]=='product'for c in rels[5]['coefficients'])
    O=D['operator'];A=D['oracle'];selected=set(A['selected34_rows']);cw=O['local_coordinate_weights'];rw=O['equation_weights']
    assert cw==A['local_coordinate_weights']and rw==A['equation_weights']
    qrows=set(D['frame']['quadric_original_rows']);crows={D['source63']['labels'][i][1]for i in range(42,63)}
    low=[];direct=[]
    for i,p in enumerate(O['local_weight_pieces']):
        if p['weight']<3:raise AssertionError('unexpected lower-weight piece')
        if p['weight']==4:direct.append(i)
        if p['weight']!=3:continue
        assert p['spatial_degree']==5
        support={2:set(),3:set()}
        for g,j,c in p['coefficients']:
            if g in selected:continue
            if j in[19,23]:continue
            if rw[g]==2:
                assert g in qrows and cw[j]==1
            else:assert rw[g]==3 and g in crows and j==-1
            support[rw[g]].add(g)
        low.append(dict(piece_index=i,operator_row=p['operator_row'],spatial_monomial=p['spatial_monomial'],
            original_quadric_support=sorted(support[2]),original_genuine_cubic_support=sorted(support[3])))
    assert len(low)==6 and len(direct)==15
    assert direct==O['weight4_quartic_constant_block_row_piece_indices']
    C=D['source319'];assert len(C['lifted_relations'])==35
    light=[0,1,2,3,4,6,7]
    assert[C['source_rows'][r['deleted_row']][1]for r in C['lifted_relations']]==[g for g in deleted[:5]for _ in light]
    assert[C['source_rows'][r['deleted_row']][3]for r in C['lifted_relations']]==light*5
    for rel in C['lifted_relations']:
        for i,c in rel['coefficients']:
            row=C['source_rows'][i]
            if row[0]=='quadric':assert row[2][5]==0
            else:assert row[0]=='cubic_normalized'
    assert D['cubic4']['new_known_unit_directions']==4
    assert D['cubic4_scope']['combined_uniform_kernel_directions_on_FULLU']==37
    # Proof: the accepted unit63->57 quotient is free and generically injective,
    # hence every R-linear kernel relation is a combination of these six.
    # Its first five light multiples are EXACTLY the35 deleted source vectors;
    # the sixth has only product rows, so its light multiples are pure-quadric
    # quadratic-cofactor relations, inside the accepted whole22 kernel.
    evidence=dict(six_cubic_deleted_rows=deleted,six_unit_pivot_coefficients=pivots,
        unit_triangular_minor='-12*a^3/5',sixth_relation_product_source_indices=[c['original_row']for c in rels[5]['coefficients']],
        six_weight3_pieces=low,direct15_weight4_piece_indices=direct,
        lower_piece_normal_graph_weight1='T*D*v only; no higher positive-weight correction',
        lift5x7_equals_deleted35=True,sixth_light_lifts_pure_quadric_no_H=True,
        source_module_kernel_argument='R63/(six unit relations)=R57; generic rank57; torsionfree free57 has zero kernel intoR91',
        all_possible_positive_weight_multipliers_to4='constants for weight4, seven light-linear multipliers for weight3',
        containing_kernel_module='old22 + direct15',containing_module_exact_free_rank_on_FULLU=37,
        Bianchi_alone_rank_not_asserted=True,
        excluded='higher-weight, divided or saturated operator relations; other syzygies; full80kernel; actual239/rank64')
    (out/'certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/f'research/notes/{ID}.md'
    result=dict(id=ID,all_passed=True,kind='Pure module proof and frozen-data binding; no coefficient search',
        cubic_whole_module_kernel_basis6_on_aFBC=True,pointwise_cubic_rank57_not_assumed=True,
        weight3_piece_count=6,direct_weight4_piece_count=15,
        degree4_operator_generated_relations_contained_in_old22_plus_direct15=True,
        old22_plus_degree4_operator_module_free_rank_on_FULLU=37,missing_for_80=43,
        operator_alone_rank37_not_claimed=True,
        no_higher_divided_saturated_full_syzygy_exhaustion=True,
        no_new_field_point_GB_coefficient_expansion_or_rank_computation=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,*paths.values()]},
        certificate_sha256=sha(out/'certificate.json'),seconds=time.monotonic()-start,
        frozen=True,independent_review_pending=True)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
