#!/usr/bin/env python3
"""Two actual cubic relation channels from preexisting B11 coefficient caches."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,sparse_rank
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_VALUES_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'functional':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase2-linear-syzygy-source-root-contractions.json',
      'adjoint_report':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/report.json',
      'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'composite':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'native_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
      'cache_peer':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json',
      'theorem':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/report.json',
      'pure':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63/actual-linear-syzygies-and-pure-r3-source-values.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'source_cache_hashes_before_array_evaluation','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};A,AR,B,C,N,NP,CP,TH,P=[D[k] for k in paths]
    assert all(d['all_passed'] for d in [AR,NP,CP,TH])
    assert NP['all_passed'] and sha(ROOT/paths['native_peer'])=='bb6e7e7c5327528f2be31b1d29db5891461d12adb319daa184b6c6addd63d693'
    assert sha(ROOT/paths['cache_peer'])=='143af8f9218d2e43626f58db9abd34827dee9e419f6ec4a9626b5ca5d36caae1'
    pairs=[tuple(p) for p in A['pairs36_order']];allowed=[(0,7),(1,7),(4,7),(5,7),(7,7)]
    pair_to_B=dict(zip(allowed,[7,8,9,10,2]));pair_to_AD={p:j for j,p in enumerate(allowed)}
    assert all(tuple(p) in pair_to_B for p in A['actual_requested_X2_pairs'])
    td={i:(F(v),rec(g),rec(f),rec(q)) for i,v,g,f,q in B['node_data']}
    cv={i:(rec(v),rec(g)) for i,v,g in C['all108_composite_values_and_opaque_second_gradients']}
    nv=[{row:rec(v) for row,v in c} for c in N['source_columns128_rows48_sparse_B11_effective_second']]
    assert N['valid_for_actual_zero_origin_source_combinations_only_by_L2_cancellation']
    all_labels={tuple(label) for row in A['all119_two_cubic_source_functionals'] for f in row['two_scalar_root_X2_functionals'] for label,c in f}
    typedroots={q[1] for q in all_labels if q[0]==0};directroots={q[1] for q in all_labels if q[0]==1}
    assert typedroots<=td.keys() and directroots<=cv.keys()
    atoms={i for root in typedroots for i in td[root][1]};assert atoms<=cv.keys()
    assert all(not cv[i][1] for i in atoms|directroots),'opaque quotient second sensitivity not eliminated'
    coefficient_bound=sum(len(v) for col in nv for v in col.values())+len(all_labels)+2*119+sum(len(td[i][1]) for i in typedroots)+5000
    assert coefficient_bound<150000
    dump(out,'preflight.json',{'all_requested_pairs_in_B11':True,'typed_roots46_or_fewer':len(typedroots),
      'direct_C0_roots':len(directroots),'all_required_typed_gradient_atoms':len(atoms),
      'all_atom_and_direct_opaque_second_gradients_zero':True,'TOTAL_support_bound':coefficient_bound,'cap':150000,
      'no_new_source_function_value_jet_or_solve':True})
    def value(label):
        kind=label[0];p=pairs[label[-1]];b=pair_to_B[p]
        if kind==2:return nv[label[1]].get(label[2],{}).get(b,F(0))
        if kind==1:return cv[label[1]][0].get(b,F(0))
        c,g,f,q=td[label[1]]
        return q.get(pair_to_AD[p],F(0))+sum((a*cv[i][0].get(b,F(0)) for i,a in g.items()),F(0))
    values={label:value(label) for label in all_labels}
    outcols=[];evals=[]
    for j,row in enumerate(A['all119_two_cubic_source_functionals']):
        assert row['syzygy119_index']==j;vec=[]
        for t,f in enumerate(row['two_scalar_root_X2_functionals']):
            z=F(row['known_Xp1_a2_correction'][t])+sum((F(c)*values[tuple(label)] for label,c in f),F(0));vec.append(z)
        assert vec[1]==F(P['all119_linear_coefficient_syzygies'][j]['actual_e5_r3_cubic_relation_value'])==0
        outcols.append({i:v for i,v in enumerate(vec) if v});evals.append({'linear_syzygy119_index':j,'two_actual_cubic_relation_values':list(map(str,vec))})
    rank=sparse_rank(outcols,2);rho=rank['rank'];assert rho<=1
    dump(out,'actual-linear-syzygy-cubic-source-values.json',{'all_scalar_root_pair_values':[[list(k),str(v)] for k,v in sorted(values.items())],
      'all119_two_actual_cubic_relation_values':evals,'all119_pure_r3_values_match_previous_source_packet':True,
      'complete_fixed_rank_certificate':rank,'actual_new_cubic_relation_channel_rank':rho,
      'actual_remaining_gr3_T_dimension_upper_bound':2-rho,
      'actual_degree3_HF_computed':False,'total_T_length_computed':False,
      'scope':'119 linear-coefficient syzygies only; constant zero-quadratic source syzygies remain uncomputed',
      'native_effective_values_used_after_zero_origin_source_combinations':True,
      'no_new_X3_source_function_value_jet_GB_or_point':True})
    dump(out,'phase.json',{'phase':'FINAL','new_cubic_channel_rank':rho,'actual_gr3_T_dimension_upper_bound':2-rho})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_linear_syzygy_cubic_relation_channel_rank':rho,
      'actual_gr3_T_dimension_upper_bound':2-rho,'actual_degree3_HF_or_total_T_length':None,
      'no_X3_new_source_function_value_jet_or_solve':True,'independent_source_and_channel_review_status':'PENDING',
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual_linear_syzygy_cubic_relation_channel_rank','actual_gr3_T_dimension_upper_bound','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
