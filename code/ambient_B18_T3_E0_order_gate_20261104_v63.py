#!/usr/bin/env python3
"""Exact pointer/product gate for E0 order>=3; no secondjet calculation."""
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import time
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T3_E0_ORDER_GATE_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    t=time.monotonic();out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    up=ROOT/'reports/AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63/unit-reduced-presentation-DAG.json'
    op=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    full,honest,unit,origin,src,jet=map(read,[fp,hp,up,op,sp,jp])
    steps=unit['all_constant_unit_pivot_steps'];assert len(steps)==1 and steps[0]['eliminated_free_generator_original_index']==5
    assert F(steps[0]['pivot_origin'])==F(10,3)
    n=full['DAG_scalar_nodes'];values=list(map(F,origin['all339669_scalar_origin_values_evaluated']))
    ca=src['families']['CA'];nu=ca['actual_tracefree_nu'];js=[i for row,i in ca['actual_positive_odd_J0_generators']]
    cubic=[q for q in full['all99_CB_kernel_and87_graph_generators_evaluated_on_CA'] if q['kind']=='graph_CB_minus_CA' and q['weight']==3]
    assert len(cubic)==6
    records=[]
    for alpha,q in enumerate(cubic):
        b=q['normalized_CA_trace_pair'][1]
        for e,j in enumerate(js):
            k=2*alpha+e;raw=full['actual12_top_T3_additional_columns'][k][0]
            assert honest['all12_top_columns_lifted_R0_10'][k][5]==0
            assert unit['all12_transformed_top_columns'][k][0]==raw
            atoms={nu,j,b}
            def flat(i):
                if i in atoms:return F(1),[i]
                x=n[i]
                if x[0]=='c':return F(x[1],x[2]),[]
                if x[0]=='n':c,l=flat(x[1]);return-c,l
                assert x[0]=='m',(i,x)
                c,l=flat(x[1]);d,r=flat(x[2]);return c*d,l+r
            if raw!=0:
                c,l=flat(raw);assert c==1 and sorted(l)==sorted([nu,j,b])
            assert values[nu]==values[j]==values[b]==0
            records.append({'top':k,'physical_cubic_row':q['raw_row'],'old_component5_node':0,
                            'raw_E0_and_reduced_E0_same_node':raw,'three_factor_nodes':[nu,j,b],
                            'three_factor_origin_values':['0','0','0'],'E0_order_lower_bound':3})
    assert all(all(F(v)==0 for v in r['known_firstjet90']) for r in jet['all12_top_results'])
    assert all(not any(map(F,r)) for q in jet['all12_top_results'] for r in q['four_g_sensitivity9_by4'])
    assert all(all(F(v)==0 for v in row[:10]) for row in jet['all2069_known_firstjet90_columns'])
    ev={'all12_exact_pointer_and_product_identities':records,
        'unitSchur_E0_formula':'rawE0-(pivotE0/unit)*oldtop_component5 = rawE0 EXACTLY',
        'all_tops_in_maximal_ideal_squared_free9_by_full9_firstjet':True,
        'all_tops_E0_in_maximal_ideal_cubed':True,
        'all2069_relation_E0_firstjets_zero':True,
        'no_source3_Honest2_or_numeric_zero_membership_inference':True}
    (out/'certificate.json').write_text(json.dumps(ev,indent=2)+'\n')
    rep={'id':STEM,'all_passed':True,'status':'typed_order_gate_for_pending_minimal_syzygy_argument',
         'all12_old_component5_literal_zero':True,'all12_reduced_E0_equals_raw_E0_pointer':True,
         'E0_order_atleast3_from_three_zero_origin_factors':True,
         'all_top_firstjets_and_all_relation_E0_firstjets_zero':True,
         'T3_zero_or_B18_depth_not_claimed_by_this_gate_alone':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),fp,hp,up,op,sp,jp]},
         'evidence_sha256':{'certificate.json':sha(out/'certificate.json')},'seconds':time.monotonic()-t}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':main()
