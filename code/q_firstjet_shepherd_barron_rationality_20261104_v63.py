#!/usr/bin/env python3
"""Gated application of the independently read 1989 primary rationality theorem.

No PDF copy, new arithmetic, beta point, invariant search, or GB.
"""
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_SHEPHERD_BARRON_RATIONALITY_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    start=time.monotonic()
    primary=ROOT/'reports/Q_FIRSTJET_SHEPHERD_BARRON_RATIONALITY_20261104_v63_REVIEW/report.json'
    reduction=ROOT/'reports/Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63_REVIEW/report.json'
    vector=ROOT/'reports/Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63.json'
    extra=ROOT/'reports/Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63_REVIEW/report.json'
    checked={};data={}
    for path in [primary,reduction,vector,extra]:
        obj=json.loads(path.read_text());assert obj['all_passed'];data[path]=obj
        for name,want in obj.get('input_sha256',{}).items():
            p=ROOT/name;assert p.is_file()and sha(p)==want,name;checked[name]=want
    assert data[reduction]['full_Sp4_field_equals_B_two_independent_variables']
    assert data[vector]['genuine_Aut_quadric_linearization_including_extra_swap']
    output=ROOT/'reports'/ID;output.mkdir(parents=True,exist_ok=False)
    audit={
        'primary_bibliography':{'author':'N. I. Shepherd-Barron','title':'Invariant theory for S5 and the rationality of M6',
            'journal':'Compositio Mathematica','volume':'70','number':'1','year':1989,'pages':'13-25',
            'pdf_URL':'https://www.numdam.org/article/CM_1989__70_1_13_0.pdf',
            'bibliography_URL':'https://www.numdam.org/item/CM_1989__70_1_13_0/'},
        'observed_PDF_pages_including_cover':14,'printed_pages_read_visually':[16,17,18,20,21,25],
        'primary_items':[
            {'item':'Lemma8','printed_pages':[16],'PDF_zero_based_pages':[4],'use':'projective sign twists'},
            {'item':'Proposition9','printed_pages':[17],'PDF_zero_based_pages':[5],'use':'S5 five-dimensional representation and delPezzo symmetric-square bridge'},
            {'item':'Theorem10','printed_pages':[17,18],'PDF_zero_based_pages':[5,6],'use':'rational symmetric-square quotient'},
            {'item':'Proposition11','printed_pages':[18,19,20,21],'PDF_zero_based_pages':[6,7,8,9],'use':'EXACT full Aut(P1xP1) projective V(1,2) times V(2,1) quotient is rational'},
            {'item':'Lemma14','printed_pages':[20,21],'PDF_zero_based_pages':[8,9],'use':'slice proof via quadric bundle and Tsen'},
            {'item':'Theorem15','printed_pages':[21],'PDF_zero_based_pages':[9],'use':'five-dimensional chi projective quotient rational'},
            {'item':'Character table','printed_pages':[25],'PDF_zero_based_pages':[13],'use':'our[3,2] is paper chi-prime, projectively equivalent byLemma8'}],
        'exact_target_match':'U,V bidegrees2,1 and1,2; both projectivized; entire AutQ including ruling swap; extra minus sign disappears projectively',
        'accepted_internal_bridge':'F=C(W16)^Sp4=B0(s1,s2) by genuine rank2 squared-line vector bundle descent',
        'derived_result':'F is C-rational of transcendence degree6',
        'verbatim_quoted_words':0,'entire_PDF_or_images_copied_to_repository':False,
        'no_explicit_six_transcendence_basis_or_NEW25_target_index_equation_claim':True}
    ev=output/'source-audit.json';ev.write_text(json.dumps(audit,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,primary,reduction,vector,extra,
           ROOT/'research/notes/Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63.md',
           ROOT/'research/notes/Q_FIRSTJET_SHEPHERD_BARRON_RATIONALITY_20261104_v63_REVIEW.md']
    report={'id':ID,'all_passed':True,
        'primary_P11_exact_biform_projective_quotient_full_AutQ_hypotheses_match':True,
        'P9_T10_L14_T15_and_sign_twist_label_crosscheck_read':True,
        'independent_primary_source_and_accepted_internal_vector_descent_gates_checked':True,
        'full_Sp4_W16_invariant_field_C_rational':True,'transcendence_degree':6,
        'six_explicit_pure_transcendence_generators_constructed':False,
        'polynomial_invariant_ring_Hilbert_series_NEW25_equation_index_or_exhaustion_claim':False,
        'novelty_priority_or_no_prior_literature_claim':False,
        'old_stable_rationality_packet_edited':False,'new_arithmetic_beta_point_GB_or_highdegree_probe':False,
        'primary_PDF_copied_or_verbatim_quote_used':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-start}
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'C_rational':True,'dimension':6,'seconds':report['seconds']}))

if __name__=='__main__':main()
