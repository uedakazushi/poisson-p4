#!/usr/bin/env python3
"""ONE exact-origin preflight and literal unit Schur reduction of honest D3."""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_full_relative_input_20261104_v63 import RegularDAG
from ambient_B18_top_graph_telescoping_20261104_v63 import independent_rows
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_HONEST10_UNIT_REDUCTION_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    hp=ROOT/'reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json'
    hr=hp.parent/'report.json'
    oldp=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/exact-residues-and-two-ideals.json'
    d=read(hp);assert read(hr)['all_passed'];old=read(oldp)
    nodes=[tuple(n) for n in d['DAG_scalar_nodes']];vals=list(map(F,old['all339669_scalar_origin_values_evaluated']))
    rvalues={q['node']:F(q['origin']) for t in d['all_pure_source_Hadamard_membership_telescopes'] for q in t['all9_regular_coordinate_divisions']}
    assert len(rvalues)==792
    new_inverse_units=[]
    for i in range(len(vals),len(nodes)):
        n=nodes[i];op=n[0]
        if op=='c':v=F(n[1],n[2])
        elif op=='v':v=F(0)
        elif op=='a':v=vals[n[1]]+vals[n[2]]
        elif op=='m':v=vals[n[1]]*vals[n[2]]
        elif op=='n':v=-vals[n[1]]
        elif op=='i':assert vals[n[1]];v=1/vals[n[1]];new_inverse_units.append([i,n[1],str(vals[n[1]])])
        elif op=='r':v=rvalues[i]
        else:raise ValueError((i,n))
        vals.append(v)
    D=RegularDAG();D.nodes=nodes;D.values=vals;D.memo={n:i for i,n in enumerate(nodes)}
    A=[list(c) for c in d['honest_D3_presentation_relation_columns_R0_10']]
    T=[list(c) for c in d['all12_top_columns_lifted_R0_10']]
    J=[list(c) for c in d['all2046_Jinf_columns_lifted_R0_10']]
    const=lambda C:[[D.values[i] for i in c] for c in C]
    rankA=len(independent_rows(const(A)));rankT=len(independent_rows(const(A+T)))
    rankJ=len(independent_rows(const(d['all36_Koszul_kernel_columns_R0_10']+J)))
    assert rankT==rankJ
    originalA=[c[:] for c in A];originalT=[c[:] for c in T]
    rownames=list(range(10));colnames=list(range(len(A)))
    P=[[D.one if i==j else D.zero for j in range(10)] for i in range(10)]
    steps=[]
    while True:
        choice=next(((j,i) for j,c in enumerate(A) for i,node in enumerate(c) if D.values[node]),None)
        if choice is None:break
        j,i=choice;pivot=A[j][:];u=pivot[i];ui=D.inverse(u)
        keep=[k for k in range(len(rownames)) if k!=i]
        def project(v):return [D.add(v[k],D.neg(D.mul(D.mul(pivot[k],ui),v[i]))) for k in keep]
        previousP=P;P=[[D.add(previousP[k][ell],D.neg(D.mul(D.mul(pivot[k],ui),previousP[i][ell]))) for ell in range(10)] for k in keep]
        A=[project(c) for k,c in enumerate(A) if k!=j];T=[project(c) for c in T]
        steps.append({'eliminated_free_generator_original_index':rownames[i],
                      'eliminated_relation_original_column':colnames[j],
                      'relation_before_projection':pivot,'pivot_node':u,'inverse_unit_node':ui,
                      'pivot_origin':str(D.values[u]),'remaining_generator_original_indices':[rownames[k] for k in keep]})
        rownames=[rownames[k] for k in keep];colnames.pop(j)
    assert len(steps)==rankA and len(rownames)==10-rankA
    assert all(D.values[i]==0 for c in A for i in c)
    # Projection is identity on each retained original generator.
    assert all(P[i][j]==(D.one if rownames[i]==j else D.zero) for i in range(len(rownames)) for j in rownames)
    top_rank=len(independent_rows(const(T)))
    assert top_rank==rankT-rankA
    exactzeroA=[j for j,c in enumerate(A) if all(i==D.zero for i in c)]
    # Dropping only literally zero columns is exact, not a jet-based decision.
    nonzeroA=[c for c in A if any(i!=D.zero for i in c)]
    ev={'DAG_scalar_nodes':D.nodes,'all792_coordinate_quotient_origin_values':[[i,str(v)] for i,v in sorted(rvalues.items())],
        'honest10_origin_relation_matrix':[[str(D.values[i]) for i in c] for c in originalA],
        'honest10_origin_top12_matrix':[[str(D.values[i]) for i in c] for c in originalT],
        'origin_rank_D3_relation_matrix':rankA,'origin_rank_D4_relation_matrix':rankT,
        'all_constant_unit_pivot_steps':steps,'remaining_original_free_generator_indices':rownames,
        'exact_projection_matrix_remaining_by_original10':P,
        'all_reduced_D3_relation_columns':A,'remaining_D3_relation_original_indices':colnames,
        'literally_zero_reduced_D3_column_positions':exactzeroA,
        'reduced_D3_relations_after_dropping_only_literal_zero_columns':nonzeroA,
        'all12_transformed_top_columns':T,
        'honest_reduced_D4_relations':nonzeroA+T,
        'reduced_top_origin_rank':top_rank,
        'all_old_and_new_nonunit_regular_functions_remain_in_same_original_base':True,
        'actual_T3':'imageoftransformed12top incokerreducedD3; relationkernelnotcomputed',
        'no_T3_zero_nonzero_colon_or_depth_inferred_from_origin_ranks':True}
    encoded=json.dumps(ev,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'unit-reduced-presentation-DAG.json').write_text(encoded)
    report={'id':STEM,'all_passed':True,'status':'exact_unit_reduced_honest_presentation_not_T3_solution',
            'D3_relation_origin_rank':rankA,'D4_relation_origin_rank':rankT,
            'constant_unit_pivots_eliminated':len(steps),'remaining_free_module_rank':len(rownames),
            'remaining_D3_relation_columns_before_literal_zero_removal':len(A),
            'remaining_D3_relation_columns_after_literal_zero_removal':len(nonzeroA),
            'remaining_top_columns':len(T),'reduced_top_origin_rank':top_rank,
            'exact_Schur_keeps_all_relations_up_to_recorded_unit_equivalence':True,
            'T3_zero_nonzero_or_depth_from_origin_claimed':False,
            'no_GB_jets_new_points_or_nonunit_inverses':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),hp,hr,oldp]},
            'evidence_sha256':{'unit-reduced-presentation-DAG.json':sha(out/'unit-reduced-presentation-DAG.json')},
            'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','D3_relation_origin_rank','D4_relation_origin_rank','remaining_free_module_rank','remaining_D3_relation_columns_after_literal_zero_removal','reduced_top_origin_rank','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
