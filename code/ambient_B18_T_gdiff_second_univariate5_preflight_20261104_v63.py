#!/usr/bin/env python3
"""Input/cost only for the genuinely needed four g-difference second values."""
import argparse,collections,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def bits(q):return max(q.numerator.bit_length(),q.denominator.bit_length())
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'fourth':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
      'fourth_report':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/report.json',
      'fourth_peer':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json',
      'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
      'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'exact_source5_inputs_bound_no_new_coefficient','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};S,K,KR,KP,C,M,H=[D[k] for k in paths]
    assert KR['all_passed'] and KP['all_passed']
    assert not M['all_third_unknown_nodes_remaining']
    remain={i for i,n in M['all_second_unknown_nodes_remaining']}
    old=[330164,330289,330524,330594];new=[566400,566401,566402,566403];assert remain==set(old+new)
    y=C['actual_YB_last2_from_saved105'];yval=[rec(v) for i,v in y];assert all(not v.get(2,F(0)) for v in yval)
    slopes={'CA':(F(1),F(1)),'CB':tuple(v[1] for v in yval)};assert slopes['CB']==(F(17,12),F(7,12))
    mon=K['monomials_bivariate_through4'];nodes=S['DAG_scalar_nodes'];records=[];totalbound=0
    for family,f in K['families'].items():
        a,b=slopes[family]
        def line(coeff):
            p=[F(0)]*5
            for (i,j),c in zip(mon,coeff):p[i+j]+=F(c)*a**i*b**j
            return p
        u=f['actual_unknown57_order'];eq=f['actual_equations57'];g=f['g_original_node'];nums=f['four_numerator_original_nodes'];needed=f['source_needed_node_ids'];roots=eq+[g]+nums
        sol=[[i,list(map(str,line(p)))] for i,p in f['all57_source4_jets']]
        gp=line(f['g_all15_coefficients']);assert not any(gp[:3]) and gp[3] and not gp[4]
        np=[line(v) for v in f['four_numerator_all15_coefficients']];assert all(not any(p[:4]) for p in np)
        h1=[F(x)*a+F(z)*b for x,z in f['quotient_firsts4_by2']]
        assert all(p[4]==gp[3]*h for p,h in zip(np,h1))
        invrows=[[[j,c] for j,c in enumerate(row) if F(c)] for row in f['fixed_J57_inverse']]
        uses=collections.Counter()
        for i in needed:
            n=nodes[i]
            if n[0] in ['a','m']:uses.update(n[1:3])
            elif n[0] in ['n','i']:uses.update(n[1:])
            elif n[0]=='v':assert i in set(u)
        uses.update(roots);live=set();peak=0
        for i in needed:
            live.add(i);peak=max(peak,len(live));n=nodes[i]
            for j in n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else []:
                uses[j]-=1
                if uses[j]==0:live.discard(j)
        sourceheight=max(bits(F(c)) for i,p in sol for c in p)
        inverseheight=max(bits(F(c)) for row in invrows for j,c in row)
        # Sequential6-coefficient evaluation; no 105-dimensional new solve.
        slotbound=6*(peak+57+len(roots))+sum(map(len,invrows))+57*6+len(needed)+5000
        assert slotbound<150000;totalbound=max(totalbound,slotbound)
        records.append({'family':family,'original_unknown57_order':u,'original_equations57':eq,
          'equation_targets':'first55 zero,last2 slope*r','last2_target_slopes':list(map(str,(a,b))),
          'cached_univariate_source_through4':sol,'new_unknown_coefficients_needed':57,
          'new_order5_formula':'u5=-fixed_J57_inverse * residual_of_cached_u0..4_with_u5=0_at_degree5',
          'fixed_J57_inverse_sparse_rows':invrows,'fixed_J57_determinant':f['fixed_J57_determinant'],
          'actual_g_original_node':g,'actual_bbar4_original_nodes':nums,'cached_g_line_through4':list(map(str,gp)),
          'cached_four_bbar_line_through4':[list(map(str,p)) for p in np],
          'cached_four_quotient_first_line':list(map(str,h1)),
          'needed_g_second_line_formula':'bbar5/g3 (g4=0 and quotient0=0)',
          'source_needed_original_DAG_nodes':needed,'source_closure_node_count':len(needed),
          'peak_live_source_node_count':peak,'new_evaluation_series_slots_per_node':6,
          'cached_source4_line_max_QQ_bitheight':sourceheight,'fixed_inverse_max_QQ_bitheight':inverseheight,
          'conservative_live_QQ_slots_and_records_bound':slotbound})
    proof={'CA_slope':[str(x) for x in slopes['CA']],'CB_slope':[str(x) for x in slopes['CB']],
      'actual_common_YB7_8_quadratic_r2_zero':True,
      'quotient_origin0_from_source3':True,
      'composition_to_actual_common_rline_second':'own_slope_second + own_first·actual_YB_curvature; curvature=0 here',
      'regularity':'g has proved order3 with nonzero known g3; divide by r³, then a true unit. Never invert zero.',
      'needs_no_common105_order4_or5':True,'new_source_coefficients_total':114,
      'all114_are_ONE_new_univariate_order_only':True,'old_primary_twofamily_source4_seconds':KR['seconds'],
      'actual_fifth_source_coefficients_or_new_g2_values_computed':False,
      'strict_order5_height_bound_not_claimed':'measured cached4/inverse heights recorded; actual 32MiB/150k/live/time guards mandatory',
      'proposed_worker_caps':{'seconds':60,'address_space_bytes':2000000000,'artifact_bytes':32*1024*1024,'expanded_live_QQ':150000}}
    dump(out,'exact-two57-univariate5-input.json',{'families':records,'physical_common_rline_binding':proof})
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL_PREFLIGHT_ONLY','new_source_coefficients_needed':114,'max_sequential_slot_bound':totalbound})
    report={'id':STEM,'all_passed':True,'new_source5_unknown_coefficients_needed':114,
      'max_sequential_live_QQ_slot_bound':totalbound,'same_common_rline_proved_by_zero_curvature':True,
      'new_source5_or_g2_values_computed':False,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','new_source5_unknown_coefficients_needed','max_sequential_live_QQ_slot_bound','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
