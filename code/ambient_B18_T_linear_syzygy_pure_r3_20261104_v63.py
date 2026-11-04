#!/usr/bin/env python3
"""One existing-second-coefficient channel on linear syzygies of37 relations."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
from ambient_B18_comparison_T_next3_source_preflight_20261104_v63 import receipt_combos,compose
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'quad':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
      'mu':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'cache':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'cache_peer':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json',
      'Xp':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
      'known_closure':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/actual-known-quadratic-cubic-closure-certificate.json',
      'complete5':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_arrays_bound_before_kernel_replay','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};L,L1,Q,M,H,A,AP,X,CUB,CP=[D[k] for k in paths]
    assert AP['all_passed'] and CP['all_passed']
    assert sha(ROOT/paths['cache_peer'])=='143af8f9218d2e43626f58db9abd34827dee9e419f6ec4a9626b5ca5d36caae1'
    assert L1['five_remaining_free_source_rows'][4]==55
    # These literal gates remove ALL first27/84 unknown Schur corrections
    # to the source55*r² coefficient for each657 minimal cycle, not just
    # the620 zero-linear combinations considered in the previous packet.
    assert all(not rec(v).get(55,F(0)) for v in X['actual_pivot_X0_columns27_sparse56'])
    assert all(not rec(v).get(8*55+7,F(0)) for v in X['actual_pivot_X1_by27_sparse56_times8'])
    assert all(not rec(v).get(8*4+7,F(0)) for v in L1['normalized84_unitcycle_Q1_columns_sparse'])
    selected=L['complete_linear_rank_certificate']['selected_original_columns'];assert len(selected)==37
    qcols=list(map(rec,Q['all464_certified_quadratic_columns_sparse'][:296]));assert len(qcols)==296
    assert all(Q['source_columns464'][8*j+k]['kind']=='variable_times_actual_linear_relation' and
      Q['source_columns464'][8*j+k]['variable']==k and
      Q['source_columns464'][8*j+k]['source_relation']==L1['actual657_minimal_relation_source_labels'][selected[j]]
      for j in range(37) for k in range(8))
    cert=L['quadratic_certificate'];ops=cert['fixed_elimination_receipts'][:296]
    rank=sum('pivot_row' in r for r in ops);assert rank==177 and sum(bool(r.get('zero')) for r in ops)==119
    bound=296*177+119*296+sum(map(len,qcols))+37*768+5000
    assert bound<150000
    dump(out,'preflight.json',{'constraint_map_shape':[180,296],'rank_from_frozen_receipts':177,
      'kernel_dimension':119,'input_QQ_nonzeros':sum(map(len,qcols)),'TOTAL_support_bound':bound,'cap':150000,
      'no_new_source_function_or_jet':True})
    normalized=[];kernel=[];source_labels=[]
    for j,op in enumerate(ops):
        c={j:F(1)}
        for p,a in op['subtract']:addscaled(c,normalized[p],-F(a))
        if 'pivot_row' in op:normalized.append({i:a/F(op['divide']) for i,a in c.items()})
        else:
            assert not compose(c,qcols);kernel.append(c);source_labels.append(j)
    assert len(kernel)==119
    wn,wzero,wlabels=receipt_combos(M['complete_rank_certificate'],741)
    base768=list(map(rec,M['actual_constant_source_combinations']))
    source37=[compose(wzero[i],base768) for i in selected]
    values={i:(rec(p),rec(g)) for i,p,g in A['all108_composite_values_and_opaque_second_gradients']}
    values[0]=({},{});assert H['DAG_scalar_nodes'][0]==['c',0,1]
    h=[];h_source=[]
    for i,c in zip(selected,source37):
        rootcoeff={}
        for j,a in c.items():
            if j>=640:continue # native bottom eight is identically zero
            root=H['honest_D0_relation_columns'][j][8]
            v=rootcoeff.get(root,F(0))+a
            if v:rootcoeff[root]=v
            else:rootcoeff.pop(root,None)
        value=F(0)
        for root,a in rootcoeff.items():
            assert root in values,('MISSING_EXISTING_ROOT',root)
            val,opaque=values[root];assert not opaque,('UNKNOWN_SECOND_ROOT',root)
            value+=a*val.get(2,F(0)) # B11 index2 is r²
        h.append(value);h_source.append({'minimal_cycle_index657':i,'actual768_constant_source_combination':encode(c),
          'literal_rawP8_last_root_coefficients':encode(rootcoeff),'actual_e5_r2_coefficient':str(value)})
    # Every known true quadratic initial form has zero literal e5*r²
    # coefficient; consequently e5*r³ annihilates its variable multiples.
    assert all(not rec(q).get(179,F(0)) for q in CUB['quadratic_basis178'])
    functionals=[];nonzero=[]
    for i,c in enumerate(kernel):
        z=sum((a*h[j//8] for j,a in c.items() if j%8==7),F(0));functionals.append(z)
        if z:nonzero.append(i)
    evidence={'selected_actual_minimal_cycle_indices37':selected,
      'actual37_e5_r2_source_coefficients':h_source,
      'linear_syzygy_coeff_order':'8*selected_relation_index+variable (z0..z6,r)',
      'all119_linear_coefficient_syzygies':[
        {'fixed_zero_receipt_source296':source_labels[i],'coefficients296_sparse':encode(c),
         'quadratic_constraint_exact_zero':True,'actual_e5_r3_cubic_relation_value':str(functionals[i])}
        for i,c in enumerate(kernel)],
      'nonzero_pure_r3_channel_syzygy_indices':nonzero,
      'all178_known_quadratics_have_e5_r2_coefficient_zero':True,
      'first27_source55_r_linear_coefficients_all_zero':True,
      'all84_unitcycle_source55_r_linear_coefficients_all_zero':True,
      'all37_source_values_from_PREEXISTING_108_root_cache':True,
      'additional_needed_root_only_literal_zero':True,
      'no_X3_new_source_function_or_jet':True,
      'nonzero_cubic_relation_would_kill_a_candidate_not_prove_lower_length':True,
      'zero_channel_does_not_exclude_other_cubic_relations':True,
      'actual_degree3_HF_or_total_T_length_not_computed':True}
    dump(out,'actual-linear-syzygies-and-pure-r3-source-values.json',evidence)
    support=sum(map(len,kernel))+sum(map(len,source37))+sum(len(r['literal_rawP8_last_root_coefficients']) for r in h_source)+sum(bool(a) for a in h)+sum(bool(a) for a in functionals)
    assert support<150000
    dump(out,'phase.json',{'phase':'FINAL','kernel_dimension':119,'nonzero_r3_functional_count':len(nonzero)})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_linear_syzygy_kernel_dimension':119,
      'selected_actual_R_relations':37,'nonzero_pure_r3_relation_functional_count':len(nonzero),
      'all_source_r2_values_from_existing_cache':True,'expanded_QQ_support':support,
      'actual_degree3_HF_total_T_length_or_nilcube_claimed':False,
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual_linear_syzygy_kernel_dimension','nonzero_pure_r3_relation_functional_count','expanded_QQ_support','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
