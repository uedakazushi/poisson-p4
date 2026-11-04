#!/usr/bin/env python3
"""Common lower619 adjoint source inventory only; no fresh coefficient values."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_COMMON_LOWER_ADJOINT_SOURCE_PREFLIGHT_20261104_v63'
def put(v,k,a):
    c=v.get(k,F(0))+a
    if c:v[k]=c
    else:v.pop(k,None)
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    def dump(n,x):
        (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths={
      'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
      'star':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json',
      'L37':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'composite':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'P3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
      'nulls':'reports/AMBIENT_B18_T_NATIVE_L2_ADJOINT_PREFLIGHT_20261104_v63/report.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};dump('phase.json',{'phase':'source_inputs_persisted_before_common_lower_fold','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};H,S,L,B,C,AD,CV,P3,MS,NU=[D[k] for k in paths]
    assert NU['all_passed'] and NU['all84_unitstar_and37_selectedL_separate_origins_zero']
    mons2=list(itertools.combinations_with_replacement(range(8),2));mi2={m:i for i,m in enumerate(mons2)}
    mons3=list(itertools.combinations_with_replacement(range(8),3));gam=rec(H['Gamma4_56_times120'])
    ops=B['operator_scalar_nodes'];bc=list(map(dict,B['B48_by640_sparse_columns']));nodes=C['DAG_scalar_nodes']
    def key(j,row,m):
        if j>=640:return (2,j-640,row,m) if row<48 else None
        if row<48:
            root=bc[j].get(row)
            if root is None or ops[root][0]=='c':return None
            op=ops[root];return (1,op[1],m) if op[0]=='actual_C0_scalar_function' else (0,root,m)
        root=C['honest_D0_relation_columns'][j][row-48+1]
        return None if nodes[root][0]=='c' else (1,root,m)
    shared=[('star',i,rec(c)) for i,c in enumerate(S['actual_unitstar84_original768_combinations'])]+[('L',i,rec(c)) for i,c in enumerate(L['actual_selected37_original768_constant_combinations'])]
    allroots=set();p3roots=set();cases=[];support=0
    for kind,i,c in shared:
        fs=[{} for _ in range(8)];p={}
        for j,a in c.items():
            if j<640:
                root=C['honest_D0_relation_columns'][j][8]
                if nodes[root][0]!='c':put(p,root,-a)
        for index,a in gam.items():
            src,m=divmod(index,120);mon=mons3[m]
            for var in set(mon):
                pair=list(mon);pair.remove(var);q=mi2[tuple(pair)]
                for j,b in c.items():
                    label=key(j,src,q)
                    if label is not None:put(fs[var],label,a*b)
        allroots.update(k for f in fs for k in f);p3roots.update(p);support+=len(c)+len(p)+sum(map(len,fs));assert support<150000
        cases.append({'shared_source_kind':kind,'index':i,'actual_original768_constant_combination':encode(c),
          'lambda_by8variables_X2_functionals':list(map(enc,fs)),'additional_rawP3_functional_at_variable_r':encode(p)})
    td={i:(rec(g),rec(q)) for i,v,g,f,q in AD['node_data']};cv={i:(rec(v),rec(g)) for i,v,g in CV['all108_composite_values_and_opaque_second_gradients']}
    typed={q[1] for q in allroots if q[0]==0};direct={q[1] for q in allroots if q[0]==1};uncachedtyped=typed-td.keys();atoms={i for root in typed if root in td for i in td[root][0]}
    uncachedC0=(direct|atoms)-cv.keys();opaque={i:cv[i][1] for i in (direct|atoms)&cv.keys() if cv[i][1]}
    jets={i for i,p in P3['raw_source_roots_with_missing_quotient2_3_set_zero_ONLY_for_known_part']};jac={i:(rec(a),rec(b)) for i,a,b in MS['source_root_Jacobians0_and_r']}
    missingP3=p3roots-jets;opaqueP3={i:jac[i] for i in p3roots&jac.keys() if any(jac[i])};missingJac=p3roots-jac.keys()
    pairs=sorted({mons2[k[-1]] for k in allroots});allowed={(0,7),(1,7),(4,7),(5,7),(7,7)}
    # Opaque sources stay symbolic even if only a few shared λ's use them.
    pre={'shared_star_columns':84,'shared_L_columns':37,'source_QQ_support':support,'scalar_X2_labels':len(allroots),
      'needed_pairs':pairs,'pairs_outside_B11':sorted(set(pairs)-allowed),'typed_roots':len(typed),'direct_C0_roots':len(direct),
      'typed_roots_outside_frozen762_AD':sorted(uncachedtyped),'C0_roots_outside_frozen108_composite_cache':sorted(uncachedC0),
      'opaque_second_sources_retained':[[i,encode(p)] for i,p in sorted(opaque.items())],
      'rawP3_roots':len(p3roots),'rawP3_roots_outside_frozen90_cache':sorted(missingP3),
      'rawP3_roots_missing_sensitivity_reader':sorted(missingJac),
      'rawP3_opaque_sensitivities_retained':[[i,encode(a),encode(b)] for i,(a,b) in sorted(opaqueP3.items())],
      'all_required_values_already_fixed':not any([uncachedtyped,uncachedC0,opaque,missingP3,missingJac,opaqueP3,set(pairs)-allowed]),
      'all121_native_origin_parts_zero_allows_effective_S2':True,'no_shared_lambda_values_evaluated':True,
      'no_new_source_jet_coefficient_or_solve':True,'cap':150000}
    dump('all121-common-lower-lambda-source-functionals.json',{'all121':cases,'monomials2_36':mons2,'Gamma4_56_times120':encode(gam),
      'lambda_formula':'lambda_b/l(v)=Gamma4(Xshared2*v)-rawP3_shared*indicator(v=r)',
      'effective_nativeS2_valid_only_after_whole_shared_constant_column_combination':True})
    dump('preflight.json',pre);dump('phase.json',{'phase':'FINAL_SOURCE_ONLY','all_cached':pre['all_required_values_already_fixed']})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,**pre,'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump('report.json',report);print(json.dumps({k:pre[k] for k in ['source_QQ_support','scalar_X2_labels','all_required_values_already_fixed','typed_roots_outside_frozen762_AD','C0_roots_outside_frozen108_composite_cache','rawP3_roots_outside_frozen90_cache']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
