#!/usr/bin/env python3
"""Projected head4*z5 SOURCE/cache gate; no new natural coefficient value."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_honest_nil_T_first_values_20261104_v63 import enc,nqq,mv
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_TOLD_TNIL_HEAD4_Z5_SOURCE_PREFLIGHT_20261105_v63'
PATHS={
 'natural':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63/actual-natural-Told-to-Tnil-source-map-and-origin.json',
 'natural_peer':'reports/AMBIENT_B18_TOLD_TNIL_NATURAL_MAP_20261104_v63_REVIEW/report.json',
 'oldunit':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'X1':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'newm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'first':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63/phase3-all849-minimal-first-relations-and-rank.json',
 'W':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase3-actual-W84-rank-and-mu.json',
 'C1':'reports/AMBIENT_B18_HONEST_NIL_T_PIVOT_FIRST_CACHED_20261105_v63/actual81-pivot-first-carries-and-six-H1-source-adjoints.json',
 'q3':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC_CLOSURE_PREFLIGHT_20261105_v63/literal-known-cubic-input-and-exact24-row-reduction.json',
 'cubic':'reports/AMBIENT_B18_HONEST_NIL_T_CUBIC917_SOURCE_OPERATOR_20261105_v63/complete-cubic917-typed-source-and-three-NF-adjoints.json',
 'rightgate':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63/right106-full-first-span-and-cache-footprint.json',
 'rightvalues':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63/phase1-all35-nine-convolution-projections.json',
 'rightreport':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63/report.json',
 'physics':'reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_CACHED_20261105_v63/phase1-common105-and-physical-source-base25.json',
 'gaugepeer':'reports/AMBIENT_B18_HONEST_NIL_T_GAUGE120_20261105_v63_REVIEW/report.json',
 'source2':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'sourcephysical':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'filtration':'reports/AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63/actual-F2-source-products-and-m5-zero-proof.json',
 'filtrationpeer':'reports/AMBIENT_B18_TOLD_TNIL_FILTRATION_R5_BOUND_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def sp(v):return [[i,c]for i,c in sorted(v.items())if c]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()};(out/'phase.json').write_text(json.dumps({'phase':'pinned inputs before fixed projected source checks','input_sha256':pins},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(D[k]['all_passed']for k in ['natural_peer','rightreport','gaugepeer','filtrationpeer'])
 N=D['natural']['actual_origin_data'];cold=list(map(Q,N['all_five_old_cycle_origins56'][4]));c0=list(map(Q,N['all_five_new_cycle_origins44'][4]));origin=list(map(Q,N['actual_origin_map_columns5'][4]))
 oldn=[[Q(c)for c in r]for r in N['old_n0_56_by24']];Z0=[[Q(c)for c in r]for r in N['actual_Z0_44_by56']]
 assert [i for i,c in enumerate(cold)if c]==[1,2,17,18,21,25,28,55]
 oldm0=[[Q(c)for c in r]for r in N['old_m0_24_by56']];assert [i for i,r in enumerate(oldm0)if r[55]]==[7] and oldm0[7][55]==-Q(125,24)
 # n=E(mE)^-1 implies the FULL cycle first: c1=-n0*m1*c0.
 # This retains -n1*m0*e55; it is not ZERO3's vanishing-m0 shortcut.
 X=D['X1'];raw=[{i:rec(v)for i,v in col}for col in X['all780_source_columns44_sparse_dual9']];xp=X['X0_pivot_columns']
 C=[i for i in range(44)if i not in D['newm0']['selected_source_columns24']];I=[C[i]for i in D['minimal']['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji=[[Q(c)for c in r]for r in D['minimal']['fixed_minor_data']['inverse']]
 assert all(not raw[j].get(42,{}).get(0,Q(0))for j in xp)
 a0=mv(Ji,[c0[i]for i in I]);l=[Q(i==42)for i in range(44)]
 for k,r in enumerate(I):l[r]-=sum((raw[xp[j]].get(42,{}).get(0,Q(0))*Ji[j][k]for j in range(18)),Q(0))
 assert l==[Q(i==42)for i in range(44)]
 y=[sum((l[i]*Z0[i][j]for i in range(44)),Q(0))for j in range(56)]
 oldadj=[-sum((y[j]*oldn[j][k]for j in range(56)),Q(0))for k in range(24)]
 expected={6:Q(1),8:-Q(5,9),12:-Q(125,36),13:-Q(5,9),14:-Q(4)};assert {i:c for i,c in enumerate(oldadj)if c}==expected
 # Actual first Schur carry, before projecting head4*z5 to gr2.
 W=D['W'];wi=[[Q(c)for c in r]for r in W['W_fixed_minor_inverse']];wr=W['W_pivot_rows'];c1=D['C1']['all81_actual_visible_C1_top2_by8_and_wedge3_r']
 hadj={h:sum((Q(c1[k][0][7])*wi[k][j]for k in range(81)),Q(0))for j,h in enumerate(wr)};hadj={h:c for h,c in hadj.items()if c};assert hadj=={27:-Q(8,3)}
 assert list(itertools.combinations(range(8),2))[27]==(6,7)
 # Thus H27=(5/2)c20,z6, since canonical m7 has no z derivatives.
 newm0=[[Q(c)for c in r]for r in D['newm0']['m0_24_by44']];assert newm0[7]==[Q(5,2)*Q(i==20)for i in range(44)]
 m2=list(itertools.combinations_with_replacement(range(8),2));im2={m:i for i,m in enumerate(m2)};pi=[rec(v)for v in D['q3']['quadratic_projection3_by180_sparse']]
 z5adj=[]
 for row in pi:
  z5adj.append({8*k+v:row[36*k+im2[tuple(sorted((5,v)))]]for k in range(5)for v in range(8)if 36*k+im2[tuple(sorted((5,v))) ]in row})
 factors=[Q(0),-Q(382163689,166960800),-Q(98809217,26713728)];assert z5adj==[({7:c}if c else{})for c in factors]
 # The actual ALL-ORDER R35 lift of the vanished first image is mandatory.
 F=D['first'];R35=F['rank_minor_columns'];rr=F['rank_minor_rows'];iv=[[Q(c)for c in r]for r in F['rank_minor_inverse']];wfirst={8*i+5:c for i,c in enumerate(origin)if c};lift=mv(iv,[wfirst.get(r,Q(0))for r in rr]);lin=[rec(F['all849_actual_minimal_first_columns40_sparse'][j])for j in R35]
 assert all(sum((lift[k]*lin[k].get(i,Q(0))for k in range(35)),Q(0))==wfirst.get(i,Q(0))for i in range(40))
 # r*gr2 basis maps to the accepted r^3 basis by I3.  Therefore saved
 # RIGHT106 convolution duals (f,r) are EXACT pi3 on M2(R35), not a guess.
 delta=[rec(v)for v in D['cubic']['cubic_NF_dual3_by600_sparse']];m3=list(itertools.combinations_with_replacement(range(8),3));im3={m:i for i,m in enumerate(m3)}
 assert all({36*k+j:delta[f][120*k+im3[tuple(sorted(mon+(7,)))]]for k in range(5)for j,mon in enumerate(m2)if 120*k+im3[tuple(sorted(mon+(7,))) ]in delta[f]}==pi[f]for f in range(3))
 G=D['rightgate'];res=[[Q(c)for c in r]for r in G['all35_by8_k42_residual_coefficients']];eff={(f,v):{(h,q):Q(c)for h,q,c in e}for f,v,e in G['all_actual_effective_C1_H2_callbacks']};pairs=list(itertools.combinations(range(8),2));us=[];vs=[]
 for f in range(3):
  us.append(sum((lift[p]*res[p][zi]*delta[f].get(120*(n+2)+119,Q(0))for p in range(35)for n,zi in enumerate([0,4,5])),Q(0)))
  vs.append(-sum((lift[p]*res[p][i]*eff[(f,7)].get((56+pairs.index((i,7)),7),Q(0))for p in range(35)for i in range(7)),Q(0)))
 assert us==vs==[Q(0)]*3
 mon=[tuple(m)for m in D['physics']['base25_monomials']];needed=sorted({tuple(sorted((i%8,7)))for i,c in enumerate(oldadj)if c});assert needed==[(0,7),(4,7),(5,7),(6,7)] and all(m in mon for m in needed)
 assert D['sourcephysical']['families']['CA']['actual_zero_odd_I0_generators'][-1]==[63,7790]
 assert '7790'in D['source2']['both57_source2']['CA']['pure_function2jets']
 budget={'actual_C4_physical_pairs_nu_live':8000,'fixed_origin_inverse_first_arrays':16000,'four_prefix_contractions':1000,'source_copy_f63_55':200,'saved35_projected_representatives':500,'fixed_lift_and_UV_receipts':4000,'next_value_and_F2_rank_receipts':2000,'metadata':4000};assert sum(budget.values())<150000
 data={'actual_head4_old_free56':55,'actual_old_c0_sparse':sp(dict(enumerate(cold))),'actual_new_c0_sparse':sp(dict(enumerate(c0))),'actual_minimal_origin_column5':origin,
  'old_cycle_first_exact_formula':'cold1=-n0*mold1*cold0 because n=E(mE)^-1; includes -n1*m0e55 exactly',
  'projected_top0_first_r_formula':'Jr=(Z1_r*cold0)_42 + oldadj*(mold1_r*cold0) -(Xp1_r*a0)_42 +(20/3)*(Z1_z6*cold0-Xp1_z6*a0)_20',
  'all18_first_origin_cleanup_a0':a0,'actual_top0_functional44':l,'actual_old_m1_r_adjoint24':oldadj,'actual_C1_H0_adjoint_sparse':sp(hadj),'canonical_H27_coefficient':Q(5,2),
  'old_m1_required_prefix_pairs':needed,'old_m1_required_sector_and_zero_index':[[i//8,i%8,c]for i,c in enumerate(oldadj)if c],
  'all_old_last_r_division_rows_zero':all(not oldadj[i]for i in[7,15,23]),'old_gA4_n1_required_unprojected_but_drops_exactly_here':True,
  'source2_binding':'Vold_Q48=(delta+bA*tauA-bB*tauB)*basis_mask, rows18,19,28,29,60,67,72,63,64,68,73,74; negative_U8 last=-gA*tauA, gA pure source7790. Same R0 and trace monics from accepted105.',
  'sectorB_prefix_frame':'zBcanonical=zA+c*r^2; substitution changes only pure-r second coefficient, hence every requested zi*r coefficient is unchanged. No linear coordinate/Gm substitution.',
  'quadratic_projection_z5_times_first_columns':list(map(sp,z5adj)),'quadratic_projection_factors_times_Jr':factors,
  'actual_R35_source_labels':R35,'literal_first_image40':sp(wfirst),'exact_R35_constant_lift35':lift,'all40_first_image_lift_equalities':True,
  'right_r_convolution_equals_actual_quad_projection':True,'lifted_R35_U_sensitivity3':us,'lifted_R35_V_sensitivity3':vs,'all_missing_third_contractions_cancel_after_actual_linear_lift':True,
  'actual_w_image_formula':'factors*Jr - sum_p lift[p]*saved_M2R35_lower_projection(f,r). These individual lower representatives become actual ONLY after this full U,V cancellation.',
  'next_F2_recipe':'Use independently saved i(e3*r),i(e4*r^2); their2x2 minor is nonzero. Solve actual w image=a*i(e3*r)+b*i(e4*r^2); verify every third coordinate. Then F3=oldm3 + span(qC,w-a*v-b*q2), with oldm2 basis(q2,qC,g).',
  'completion_criterion':'Here qC=e4*z5*r, q2=e4*r^2, g=e4*r^3: r*qC=0 and r*q2=g, while r*v=-(11767/3652)*qC+alpha*g. Thus r*(w-a*v-b*q2)=(1+a*(11767/3652))*qC-(a*alpha+b)*g; if its qC coefficient is nonzero then rF3 nonzero irrespective alpha, forcing exact HF55321. If zero retain only actual alpha source profile, no guessed value.',
  'whole_worker_QQ_upper':sum(budget.values()),'budget':budget,'all_inputs_already_cached':True,'new_physical_source_gauge_jet_or_higher_callback':False,'actual_Jr_w_or_tail_values_evaluated':False}
 assert nqq(data)<150000
 p=out/'head4-z5-projected-source-and-cache-gate.json';blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'SOURCE_GATE_CLOSED_projected_head4_z5_cache_only','retained_old_n1_term_exact':True,'old_last_division_and_new_m2_callbacks_project_out':True,'all_R35_lift_UV_sensitivities_zero':True,'whole_worker_QQ_upper':sum(budget.values()),'actual_natural_w_or_HF_values_computed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
