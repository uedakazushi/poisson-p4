#!/usr/bin/env python3
"""Exact Q1/Q2 column-subset honest lifts/pivot/first spans, no new jets."""
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import time
from ambient_B18_top_graph_telescoping_20261104_v63 import independent_rows
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T1_T2_FIRSTJET_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    full,honest,unit,jet=map(read,[fp,hp,up,jp]);Q3=full['actual_Qn_R0_squared_generator_matrices']['3']
    lookup={tuple(pair):i for i,pair in enumerate(Q3)}
    active=unit['remaining_D3_relation_original_indices'];activepos={i:j for j,i in enumerate(active)}
    step=unit['all_constant_unit_pivot_steps'][0];pivot=step['eliminated_relation_original_column'];assert pivot==1389
    assert step['eliminated_free_generator_original_index']==5 and F(step['pivot_origin'])==F(10,3)
    pivotQ=pivot-36;pivotpair=Q3[pivotQ];assert honest['all2034_Q3_columns_lifted_R0_10'][pivotQ]==step['relation_before_projection']
    results={}
    for n,size in [(1,1513),(2,1926)]:
        Q=full['actual_Qn_R0_squared_generator_matrices'][str(n)];assert len(Q)==size
        q3ids=[lookup[tuple(p)] for p in Q]
        assert all(Q3[i]==p for i,p in zip(q3ids,Q))
        at=next(i for i,p in enumerate(Q) if p==pivotpair)
        assert q3ids[at]==pivotQ
        alloriginal=list(range(36))+[36+i for i in q3ids]
        chosen=[activepos[i] for i in alloriginal if i!=pivot]
        cols=[unit['all_reduced_D3_relation_columns'][i] for i in chosen]
        jets=[list(map(F,jet['all2069_known_firstjet90_columns'][i])) for i in chosen]
        assert len(cols)==36+size-1
        assert all(not any(F(v) for r in jet['all2069_four_g_sensitivities'][i] for v in r) for i in chosen)
        assert all(all(v==0 for v in row[:10]) for row in jets)
        piv=independent_rows(jets)
        results[str(n)]={'raw_Q_column_count':size,'all_raw_pair_to_Q3_column_indices':q3ids,
                         'unit_graph52_column_index_in_raw_Q':at,'same_unit_pivot10_by_actual_membership':True,
                         'honest10_relation_count':36+size,'unit_eliminated_original_generator':5,
                         'honest9_relation_columns':cols,'honest9_column_source_positions_in_Q3_reduced':chosen,
                         'honest9_relation_count':len(cols),'all_exact_firstjet90_columns':[[str(v) for v in row] for row in jets],
                         'independent_firstjet_column_indices':piv,'exact_first_span_rank':len(piv),
                         'all_four_g_sensitivities_zero':True,'all_E0_firstjets_zero':True}
    (out/'Q1-Q2-input-and-firstjets.json').write_text(json.dumps({'same_unit_projection9_by10':unit['exact_projection_matrix_remaining_by_original10'],
        'retained_original_free_generators':unit['remaining_original_free_generator_indices'],'families':results,
        'no_actual_membership_or_T1_T2_zero_from_first_span_alone':True},separators=(',',':'))+'\n')
    rep={'id':STEM,'all_passed':True,'status':'actual_honest_Q1_Q2_inputs_and_exact_firstspans',
         'Q1_raw_columns':1513,'Q2_raw_columns':1926,'Q1_Q2_same_unit_graph52_pivot_present':True,
         'Q1_honest9_relations':results['1']['honest9_relation_count'],'Q2_honest9_relations':results['2']['honest9_relation_count'],
         'Q1_first_span_rank':results['1']['exact_first_span_rank'],'Q2_first_span_rank':results['2']['exact_first_span_rank'],
         'T1_T2_zero_cyclic_or_depth_not_claimed_without_quadratic_gate':True,
         'no_new_prefix_IFT_source_jets_GB_or_points':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),fp,hp,up,jp]},
         'evidence_sha256':{'Q1-Q2-input-and-firstjets.json':sha(out/'Q1-Q2-input-and-firstjets.json')},'seconds':time.monotonic()-t}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':main()
