#!/usr/bin/env python3
"""Read-only source/live bound and exact existing raw33 multiplication table.

Does not evaluate any physical source unknown or physical difference action.
"""
import argparse
import collections
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,mv,product_columns,eye,srec
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def parse(v):return {i:F(c) for i,c in v}
def main(out):
    out=out.resolve();out.mkdir(parents=True,exist_ok=False);hashes={}
    def read(path):
        p=ROOT/path;hashes[path]=sha(p);return json.loads(p.read_text())
    rawreport=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/report.json')
    peer=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json')
    assert rawreport['all_passed'] and peer['all_passed']
    raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json')
    source=read('reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    joint=read('reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json')
    old=read('reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json')
    artpeer=read('reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63_REVIEW/report.json')
    assert artpeer['all_passed']
    basis=raw['actual33_basis'];actions=[[parse(c) for c in A] for A in raw['coordinate_actions9_sparse_columns']]
    degrees=[a+l+(t!=0) for t,a,l in basis]
    assert max(degrees)==4 and collections.Counter(degrees)=={0:1,1:9,2:16,3:5,4:2}
    assert basis[0]==[0,0,0]
    for A in actions:
        for j,v in enumerate(A):assert all(degrees[k]>=degrees[j]+1 for k in v)
    table=[];table_terms=0
    for t,a,l in basis:
        M=eye(33)
        for _ in range(a):M=product_columns(actions[7],M)
        for _ in range(l):M=product_columns(actions[8],M)
        if t:M=product_columns(actions[t-1],M)
        table.append(M);table_terms+=sum(map(len,M))
    assert all(M[0]=={i:F(1)} for i,M in enumerate(table))
    assert all(table[i][j]==table[j][i] for i in range(33) for j in range(33))
    # Associativity on every source basis triple, no source evaluation.
    for i in range(33):
        for j in range(33):
            for k in range(33):
                lhs={};rhs={}
                for a,c in table[i][j].items():lhs=add(lhs,table[a][k],c)
                for a,c in table[j][k].items():rhs=add(rhs,table[i][a],c)
                assert lhs==rhs
    nodes=source['DAG_scalar_nodes'];order=old['needed_closure_order'];roots=old['requested_root_nodes']
    ids=joint['all105_unknown_original_node_order'];eqs=joint['all105_equation_original_node_order']
    assert len(ids)==len(eqs)==105 and len(order)==18584 and len(roots)==164
    assert old['unknown105_order']==ids and old['equation105_order']==eqs
    uses=collections.Counter(roots)
    for i in order:
        if nodes[i][0] in ['a','m','n','i']:uses.update(nodes[i][1:])
    live=set();peak=0
    for i in order:
        live.add(i);peak=max(peak,len(live))
        if nodes[i][0] in ['a','m','n','i']:
            for j in nodes[i][1:]:
                uses[j]-=1
                if not uses[j]:live.remove(j)
    assert peak==414 and live==set(roots)
    jterms=sum(bool(F(c)) for row in joint['fixed105_Jacobian_inverse'] for c in row);assert jterms==866
    # Conservative bound includes streamed artifacts plus simultaneous work and
    # every source-DAG metadata record, even though many do not contain a QQ.
    bounds={'raw33_multiplication_table':table_terms,'source_DAG_metadata':len(order),
      'live_source_nodes_dense':33*peak,'current_unknowns':33*105,
      'residual_and_correction':2*33*105,'all4_saved_residual_and_correction_arrays':4*2*33*105,
      'final164_source_roots':33*164,'final6_inverse_products':6*2*33,
      'all12_physical45actions_dense':12*45*45,'quotient_echelon_dense':45*45,
      'accepted_raw33_metadata_QQ_bound':602,'J105_inverse_nonzeros':jterms}
    total=sum(bounds.values());assert total<150000
    # Save the table for reuse: the source job need not rebuild any raw base.
    tablefile=out/'exact-existing-raw33-multiplication-table.json'
    tablefile.write_text(json.dumps({'basis33':basis,'basis_degrees':degrees,
      'table_columns_for_each_left_basis':[[srec(v) for v in M] for M in table],
      'table_nonzeros':table_terms,'all_associativity333':True,'m5_zero':True},separators=(',',':'))+'\n')
    gate={
      'job_kind':'ONE NEW finite-algebra source105 evaluation; no full formal jet or gauge solve',
      'source_unknown105':ids,'source_equation105':eqs,'source_closure_order':order,
      'source_requested164_roots':roots,'unit_inverse_nodes':[i for i in order if nodes[i][0]=='i'],
      'target_CA_F9':'zA0,...,zA6,rA,rA+ell/5',
      'source_uniqueness':'same original unit J105; in actual finite local raw33, four constant-J corrections solve all105 exactly because m5=0',
      'stage_protocol':'save residual and correction only for4 stages; reconstruct successive unknowns by their sum; final105unknowns stored once',
      'raw45_action_formula':'[[delta,0,-nuB*bB],[0,delta,nuA*bA],[-bB,bA,delta]] with actual6quotients/common6lifts',
      'required_checks':['all105equations exact','all6 unit inverses exact','allphysical9readbacks same source',
        'all64 actual conductor source products zero','piA(bA)=0 and piB(bB)=0 before branch lifts',
        'IA*bB andIB*bA annihilation','all12×9 coordinate intertwiners','all66 physical commutators',
        'all540 source relation columns and quotient receipts','ell0 source/action quotient agrees with accepted raw26/Cbar8'],
      'right_exact_target_only':'actual C/ell²C from raw45 quotient; NO kernel specialization',
      'QQ_slots_and_metadata_bounds':bounds,'total_bound':total,
      'limits':{'seconds':60,'address_bytes':2_000_000_000,'aggregate_artifact_bytes':32*1024*1024,
        'expanded_live_total_slots':150000,'coefficient_numerator_or_denominator_bits':384},
      'all_numerical_physical_source_evaluation_not_run':True,
      'T2_comparison_H0_length_and_torsion_exponents_not_computed':True}
    gp=out/'whole-C2-source-and-budget-gate.json';gp.write_text(json.dumps(gate,separators=(',',':'))+'\n')
    report={'id':STEM,'all_passed':True,'scope':'source/budget preflight; C2 worker not executed',
      'source_nodes':len(order),'live_nodes':peak,'raw33_table_nonzeros':table_terms,'whole_slot_bound':total,
      'new_physical_source_values_evaluated':False,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in [tablefile,gp]}}
    (out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
    print(json.dumps({'all_passed':True,'table_nonzeros':table_terms,'whole_slot_bound':total,'source_values_run':False}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
