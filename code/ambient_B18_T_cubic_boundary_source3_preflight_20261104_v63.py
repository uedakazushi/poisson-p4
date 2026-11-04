#!/usr/bin/env python3
"""Known-array corrected d3 cubic channel and exact cached-source order audit."""
import argparse,collections,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
from ambient_B18_comparison_T_next3_source_preflight_20261104_v63 import receipt_combos,compose
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'quad':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'pure':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63/actual-linear-syzygies-and-pure-r3-source-values.json',
      'resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
      'mu':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'origin':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'cache':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'source3':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
      'source3peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
      'source4':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
      'source4peer':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json'}
    # Older source3 review was named JOINT_CURVE rather than RESTRICTIONS.
    if not(ROOT/paths['source3peer']).exists():paths['source3peer']='reports/AMBIENT_B18_T3_JOINT_CURVE_20261104_v63_REVIEW/report.json'
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'all_frozen_input_hashes_before_known_array_check','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};L,Q,L1,P,R,M,X0,H,C,J,JP,K,KP=[D[k] for k in paths]
    assert JP['all_passed'] and KP['all_passed']
    qcols=list(map(rec,Q['all464_certified_quadratic_columns_sparse']))
    ops=L['quadratic_certificate']['fixed_elimination_receipts'][:296]
    h=list(map(F,[r['actual_e5_r2_coefficient'] for r in P['actual37_e5_r2_source_coefficients']]))
    normalized=[];theta_basis=[];kernel_checks=0
    for j,op in enumerate(ops):
        c={j:F(1)}
        for p,a in op['subtract']:addscaled(c,normalized[p],-F(a))
        value=sum((a*h[i//8] for i,a in c.items() if i%8==7),F(0))
        if 'pivot_row' in op:
            c={i:a/F(op['divide']) for i,a in c.items()};normalized.append(c);theta_basis.append(value/F(op['divide']))
        else:assert not value;kernel_checks+=1
    assert len(theta_basis)==177 and kernel_checks==119
    # The complete actual d3 has no constant terms; its linear coefficients
    # involve only z0..z6, and the source55*r part of every Ystar1 is zero.
    assert all(not rec(v).get(39,F(0)) for v in L1['normalized84_unitcycle_Q1_columns_sparse'])
    d3=R['differentials8_column_sparse'][2];assert len(d3)==168
    for col in d3:
        for row,poly in col:
            for mon,c in poly:
                assert sum(mon)>0
                assert not(sum(mon)==1 and mon[7]==1 and F(c))
    factors=list(map(rec,L['quadratic_certificate']['all_column_factorization_in_normalized_basis']))
    boundary=[]
    for j in range(168):
        f=factors[296+j];value=sum((c*theta_basis[p] for p,c in f.items()),F(0))
        boundary.append({'literal_d3_source_column':j,'quadratic_lift_coefficients296':encode(compose(f,normalized)),
          'actual_minimal_d3_pure_r3':str(F(0)),'linear_lift_pure_r3':str(value),
          'corrected_boundary_pure_r3':str(-value)})
    nonzero=[b['literal_d3_source_column'] for b in boundary if F(b['corrected_boundary_pure_r3'])]
    dump(out,'phase1-corrected-d3-pure-r3-channel.json',{'theta_on_normalized177_quadratic_basis':list(map(str,theta_basis)),
      'all119_linear_syzygies_annihilated_so_theta_is_well_defined':True,
      'actual_d3_linear_r_coefficient_identically_zero':True,'all84_Ystar1_source55_r_zero':True,
      'all168_corrected_boundaries':boundary,'nonzero_boundary_channel_indices':nonzero,
      'no_actual_source3_coefficients_evaluated':True})
    # Existing second gates for later rawP3 source dependence. Native bottom
    # eight rows are zero at all orders. Keep missing roots explicit.
    vals={i:(rec(v),rec(g)) for i,v,g in C['all108_composite_values_and_opaque_second_gradients']};vals[0]=({}, {})
    base=list(map(rec,M['actual_constant_source_combinations']));wn,wz,wl=receipt_combos(M['complete_rank_certificate'],741)
    originalpivots=[r for r in X0['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
    xp=[rec(r['constant_original_source_combination']) for r in originalpivots];star=[compose(c,base) for c in wn]
    def rawp_rr(c):
        roots={}
        for j,a in c.items():
            if j>=640:continue
            root=H['honest_D0_relation_columns'][j][8]
            z=roots.get(root,F(0))+a
            if z:roots[root]=z
            else:roots.pop(root,None)
        missing=[i for i in roots if i not in vals or vals[i][1]]
        value=None if missing else sum((a*vals[i][0].get(2,F(0)) for i,a in roots.items()),F(0))
        return {'rawP8_root_coefficients':encode(roots),'missing_or_opaque_second_roots':missing,'known_rr':None if value is None else str(value)}
    gates={'normalized_Xp27_source55_rr':[rawp_rr(c) for c in xp],
      'normalized_Ystar84_source55_rr':[rawp_rr(c) for c in star],
      'normalization_uses_only_constant_changes_for_rr':'all source55*r firsts zero, source55 constant zero',
      'no_new_source_values_computed':True}
    dump(out,'phase2-known-Schur-second-r-gates.json',gates)
    # Inventory exact source orders already saved. The JOINT105 is in actual
    # common YA; own-plane57 caches are in each family's own last2 F9 coords.
    assert J['all105_equations_zero_through_degree3']
    mon3=J['monomials'];assert len(mon3)==10 and max(map(sum,mon3))==3
    common=[]
    for i,coeff in J['all105_unknown_jets']:
        assert len(coeff)==10;rline=[F(0)]*4
        for m,c in zip(mon3,coeff):rline[sum(m)]+=F(c)
        common.append([i,list(map(str,rline))])
    assert len(common)==105
    own=[]
    for family,f in K['families'].items():
        assert len(f['all57_source4_jets'])==57 and f['all57_equations_all15_coefficients_verified']
        assert all(len(p)==15 for i,p in f['all57_source4_jets'])
        own.append({'family':family,'source_unknowns57':[i for i,p in f['all57_source4_jets']],
          'coordinate_scope':'own F9 first7 zero, own last2 arbitrary','cached_maximum_total_degree':4,
          'all15_monomials_and_all57_equations_verified':True})
    profile={'actual_common105_source3_rline':common,'actual_common105_scope':'YA first7 zero, rA=sA=r (ell=0)',
      'own_CA_CB_source4_inventory':own,
      'copied_FULLINV15_third_recipe':'compose ownCB57 source3 (contained in source4) with actual commonYB7/8 source3',
      'actual_R0_function_r3_values_are_existing_cache_evaluations_not_new_source_jets':True,
      'generic_missing_Hadamard_h2_r2_requires':'mixed zi*r² numerator3 (not in cached own-plane source3)',
      'generic_missing_Hadamard_h3_r3_requires':'mixed zi*r³ numerator4 (not in cached own-plane source4)',
      'generic_missing_g2_r2_requires':'own-plane numerator5 after cubic denominator',
      'generic_missing_g3_r3_requires':'own-plane numerator6 after cubic denominator',
      'these_missing_coefficients_must_be_retained_until_actual_adjoint_cancellation':True,
      'new_source_jet_solve_authorized_or_run':False}
    dump(out,'phase3-exact-existing-source-order-inventory.json',profile)
    note=ROOT/'research/notes'/f'{STEM}.md'
    dump(out,'phase.json',{'phase':'FINAL','nonzero_corrected_boundary_channels':nonzero,'new_source3_evaluation':False})
    report={'id':STEM,'all_passed':True,'corrected_d3_pure_r3_nonzero_count':len(nonzero),
      'normalized_Xp27_rr_known_count':sum(r['known_rr'] is not None for r in gates['normalized_Xp27_source55_rr']),
      'normalized_Ystar84_rr_known_count':sum(r['known_rr'] is not None for r in gates['normalized_Ystar84_source55_rr']),
      'existing_common_source3_unknowns':105,'existing_own_plane_source4_unknowns_each':57,
      'no_new_sourcejet_coefficient_solve_or_true_HF3_claim':True,
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','corrected_d3_pure_r3_nonzero_count','normalized_Xp27_rr_known_count','normalized_Ystar84_rr_known_count','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
