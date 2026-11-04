#!/usr/bin/env python3
"""Actual three-by620 quadratic image from frozen old/native values only."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,sparse_rank
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'old':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'oldreport':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/report.json',
      'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'native_report':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/report.json',
      'native_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
      'nativeL2':'reports/AMBIENT_B18_T_NATIVE_ADAPTER2_CANCELLATION_PREFLIGHT_20261104_v63/report.json',
      'complete5':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW/report.json',
      'depthprimary':'reports/AMBIENT_B18_ORIGINAL_DEPTH19_COROLLARY_20261104_v63/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()}
    dump(out,'phase.json',{'phase':'frozen_source_hashes_before_small_assembly','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};O,OR,N,NR,NP,L2,C,L,DEP=[D[k] for k in paths]
    assert all(r['all_passed'] for r in [OR,NR,NP,L2,C,L,DEP])
    assert sha(ROOT/paths['native_peer'])=='bb6e7e7c5327528f2be31b1d29db5891461d12adb319daa184b6c6addd63d693'
    assert C['actual_complete_minimal_generator_presentation']==[5,825] and L['actual_linear_relation_rank']==37
    assert N['valid_for_actual_zero_origin_source_combinations_only_by_L2_cancellation']
    assert not N['individual_actual_native_X2_claimed']
    assert OR['active_division_second_pair_count']==0
    for label in ['oldreport','native_report','native_peer','nativeL2','complete5','low']:
        d=D[label]
        for p,h in d.get('input_sha256',{}).items():assert sha(ROOT/p)==h
        for h in d.get('input_hashes',{}).values():
            if isinstance(h,dict) and 'path' in h:assert sha(ROOT/h['path'])==h['sha256']
    native=[{row:rec(v) for row,v in col} for col in N['source_columns128_rows48_sparse_B11_effective_second']]
    assert len(native)==128 and len(O['620_projected_source_functionals'])==620
    pair_to_B=[7,8,9,10,2]
    source_support=sum(len(v) for col in native for v in col.values())
    input_support=sum(sum(len(c) for c in row['three_native_X2_functionals']) for row in O['620_projected_source_functionals'])
    assert source_support+input_support+3*620*4<150000
    dump(out,'preflight.json',{'native_effective_QQ_nonzeros':source_support,'fixed_native_functional_QQ_nonzeros':input_support,
      'output_matrix_shape':[3,620],'total_bound':source_support+input_support+3*620*4,'cap':150000,
      'no_source_coefficients_or_IFT_solve':True})
    columns=[];evidence=[]
    for j,row in enumerate(O['620_projected_source_functionals']):
        assert row['index']==j and all(not r for r in row['three_unknown_division_second_functionals'])
        old=list(map(F,row['three_known_oldB_P_and_pivot_correction_values']));new=[]
        for t,f in enumerate(row['three_native_X2_functionals']):
            z=sum((F(c)*native[a].get(b,{}).get(pair_to_B[p],F(0)) for (a,b,p),c in f),F(0));new.append(z)
        values=[a+b for a,b in zip(old,new)];assert values[2]==0
        columns.append({t:a for t,a in enumerate(values) if a})
        evidence.append({'zero_linear_combination':j,'old_actual_and_pivot_correction':list(map(str,old)),
          'native_effective_second_contribution':list(map(str,new)),
          'actual_projected_quadratic_vector':list(map(str,values))})
    rank=sparse_rank(columns,3);d=rank['rank'];assert d<=2
    hf=[5,3,3-d];length=sum(hf)
    ev={'actual620_three_class_columns':list(map(encode,columns)),'all620_source_value_splits':evidence,
      'all620_actual_third_components_zero':True,'complete_fixed_rank_certificate':rank,
      'actual_projected_quadratic_rank':d,'actual_dim_gr2_T':3-d,'actual_HF_through_degree2':hf,
      'actual_length_T_mod_m3':length,'complete825_zero_linear_initial_space_argument_required':True,
      'native_effective_E_used_only_after_source_zero_origin_combinations':True,
      'actual_total_T_length':None,'higher_T_layers_computed':False,'full_nil_cube_computed':False}
    dump(out,'actual-three-by620-quadratic-image-certificate.json',ev)
    dump(out,'phase.json',{'phase':'FINAL','rank':d,'HF':hf,'length_T_mod_m3':length})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_three_by620_rank':d,'actual_HF_through_degree2':hf,
      'actual_length_T_mod_m3':length,'actual_length_T_total':None,'higher_T_layers_or_full_nil_cube_claimed':False,
      'independent_full_source_and_rank_review_status':'PENDING',
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual_three_by620_rank','actual_HF_through_degree2','actual_length_T_mod_m3','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
