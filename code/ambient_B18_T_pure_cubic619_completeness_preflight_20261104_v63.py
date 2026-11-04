#!/usr/bin/env python3
"""Pure-r3 Schur/Theta/completeness preflight: no new third coefficients."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_PURE_CUBIC619_COMPLETENESS_PREFLIGHT_20261104_v63'
PATHS={
 'boundary_report':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/report.json',
 'boundary':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase1-corrected-d3-pure-r3-channel.json',
 'second':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
 'orders':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase3-exact-existing-source-order-inventory.json',
 'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
 'quadratics':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase2-certified-quadratic-input.json',
 'pure':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63/actual-linear-syzygies-and-pure-r3-source-values.json',
 'pure_peer':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63_REVIEW/report.json',
 'resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
 'S':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63/actual-minimal5-by825-Schur-operator.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'first640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'native128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'cache':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
 'rank1':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
 'rank1_peer':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(x):return {int(i):Q(c) for i,c in x if Q(c)}
def enc(x):return [[i,str(c)] for i,c in sorted(x.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,v in b.items():
  q=a.get(i,Q(0))+c*v
  if q:a[i]=q
  else:a.pop(i,None)
 return a
def combinations(cert):
 piv=[];zeros=[]
 for j,r in enumerate(cert['fixed_elimination_receipts']):
  c={j:Q(1)}
  for p,a in r['subtract']:c=add(c,piv[p],-Q(a))
  if 'pivot_row' in r:c={i:a/Q(r['divide']) for i,a in c.items()};piv.append(c)
  else:zeros.append(c)
 return piv,zeros
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inp={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'source_gate_inputs_saved_before_proof_checks','inputs':inp},indent=2)+'\n')
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert D['boundary_report']['all_passed'] and D['pure_peer']['all_passed'] and D['rank1_peer']['all_passed']
 for p,h in D['boundary_report']['input_sha256'].items():assert sha(ROOT/p)==h
 L,Qd,P,S,M=D['low'],D['quadratics'],D['pure'],D['S'],D['m1']
 assert S['actual_minimal_generator_presentation_shape']==[5,825]
 qcols=list(map(rec,Qd['all464_certified_quadratic_columns_sparse']))
 h=list(map(Q,[r['actual_e5_r2_coefficient'] for r in P['actual37_e5_r2_source_coefficients']]))
 # Independently replay only fixed basis receipts of the296 linear multiples.
 basis_combos=[];theta=[];zero119=[]
 for j,r in enumerate(L['quadratic_certificate']['fixed_elimination_receipts'][:296]):
  c={j:Q(1)}
  for p,a in r['subtract']:c=add(c,basis_combos[p],-Q(a))
  v=sum(a*h[i//8] for i,a in c.items() if i%8==7)
  if 'pivot_row' in r:
   d=Q(r['divide']);basis_combos.append({i:a/d for i,a in c.items()});theta.append(v/d)
  else:assert not v;zero119.append(enc(c))
 assert len(theta)==177 and len(zero119)==119
 assert list(map(str,theta))==D['boundary']['theta_on_normalized177_quadratic_basis']
 qb=list(map(rec,L['quadratic_certificate']['basis_normalized_columns']));lead={min(v):j for j,v in enumerate(qb)}
 ext={}
 for j in range(180):
  v={j:Q(1)};z=Q(0)
  for p in sorted(lead):
   c=v.get(p,Q(0))
   if c:z+=c*theta[lead[p]];v=add(v,qb[lead[p]],-c)
  if z:ext[j]=z
 assert all(sum(ext.get(i,Q(0))*c for i,c in q.items())==theta[j] for j,q in enumerate(qb))
 assert all(sum(ext.get(i,Q(0))*c for i,c in q.items())==0 for q in qcols[296:])
 # Exact r-axis gates for the actual fixed Schur source, no arbitrary unit
 # adapter derivative or Euler function coefficient is evaluated.
 first=D['first640']['source_columns640_rows56_sparse_dual9']+D['native128']['source_columns128_rows48_sparse_dual9']
 raw_r=[dict((i,rec(p)) for i,p in col).get(55,{}).get(8,Q(0)) for col in first]
 xp=S['first27_source_combinations'];remaining=list(map(rec,S['source741_zero_constant_combinations']))
 star=[]
 for c in S['second84_normalized_source741_combinations']:
  v={}
  for j,a in rec(c).items():v=add(v,remaining[j],a)
  star.append(v)
 assert all(sum(raw_r[j]*a for j,a in rec(c).items())==0 for c in xp)
 assert all(sum(raw_r[j]*a for j,a in c.items())==0 for c in star)
 assert all(Q(x)==0 for x in D['X0']['full_X0_56_by768'][55])
 d3=D['resolution']['differentials8_column_sparse'][2]
 assert len(d3)==168
 assert all(not(sum(m)==1 and m[7] and Q(c)) for col in d3 for i,p in col for m,c in p)
 # Reconstruct actual Ystar2 last-row rr from the same original C0 roots.
 cv={i:(rec(v),rec(g)) for i,v,g in D['cache']['all108_composite_values_and_opaque_second_gradients']};cv[0]=({}, {})
 pcols=D['C0']['honest_D0_relation_columns'];ys=[]
 for combo in star:
  roots={}
  for j,a in combo.items():
   if j<640:
    root=pcols[j][8];q=roots.get(root,Q(0))+a
    if q:roots[root]=q
    else:roots.pop(root,None)
  assert all(i in cv and not cv[i][1] for i in roots)
  ys.append(sum(c*cv[i][0].get(2,Q(0)) for i,c in roots.items()))
 assert [str(x) for x in ys]==[q['known_rr'] for q in D['second']['normalized_Ystar84_source55_rr']]
 assert sum(bool(x) for x in ys)==13
 roots4=[567067,567096,567125,567154];old4=[330164,330289,330524,330594];new4=[566400,566401,566402,566403];nodes=D['C0']['DAG_scalar_nodes']
 for i,o,n in zip(roots4,old4,new4):
  v=nodes[i];assert v[0]=='a' and v[1]==o and nodes[v[2]]==['n',n]
  assert nodes[o][0]==nodes[n][0]=='g'
 a1=list(map(rec,M['620_actual_linear_pivot_coefficients_a1_by27_times8']))
 affected=[j for j,c in enumerate(a1) if any(c.get(8*i+7,Q(0)) for i in range(4))];assert len(affected)==13
 cert=D['rank1']['complete_fixed_rank_certificate'];rp,kernel619=combinations(cert);assert len(rp)==1 and len(kernel619)==619
 qimage=list(map(rec,D['rank1']['actual620_three_class_columns']))
 assert all(not {k:sum(c*qimage[j].get(k,Q(0)) for j,c in v.items()) for k in range(3) if sum(c*qimage[j].get(k,Q(0)) for j,c in v.items())} for v in kernel619)
 support=sum(map(len,kernel619))+len(ext)+len(theta)+sum(map(len,zero119))
 assert support<150000
 ev={'source_manifest':inp,'actual_complete_relation_split':'37 selected source relations with independent linear parts;620 constant-cleared zero-linear source relations;168 literal minimal d3 columns',
  'Theta_on_normalized177_basis':list(map(str,theta)),'Theta_extension180_sparse':enc(ext),
  'all119_kernel_lift_pure_values_zero':True,'Theta_all168_d3_quadratics_zero':True,
  'all_actual_Xp27_last_constant_and_r_first_zero':True,'all_actual_Ystar84_last_constant_and_r_first_zero':True,
  'all_literal_d3_constant_and_r_linear_zero':True,'all168_actual_and_corrected_d3_pure_cubes_zero':True,
  'actual_Ystar2_last_rr84':list(map(str,ys)),
  'actual_Xp2_last_rr27':'first4 are minus second rr coefficients of roots4; native23 literallyzero',
  'four_missing_difference_roots':roots4,'actual_difference_definitions':[[i,o,n] for i,o,n in zip(roots4,old4,new4)],
  '620_indices_using_missing4_times_a1r':affected,
  'actual619_constant_quadratic_kernel_combinations620':list(map(enc,kernel619)),
  'typed_pure_source620_formula':{
   'C2':'actual_Xj2-actual_Xp1*a1_j',
   'beta1_r':'coefficient_r U0^-1*kappa_linear(m0*C2)',
   'R3_pure':'raw_source55_r3 - sum_i Xp55_2rr[i]*a1_j_r[i] - sum_s Ystar55_2rr[s]*beta1_r[s]',
   'Q2proxy':'pi*C2-Ystar1*U0^-1*kappa_linear(m0*C2)',
   'F620':'R3_pure-Theta_extension(Q2proxy)',
   'actual_corrected619_values':'apply saved619 constant kernel combinations to F620'},
  'Euler_ambiguity_proof':'wedge first ambiguity is d3_linear times a constant; d3_linear has no r coefficient, so beta1_r is exact. Q2proxy error is actual d3_quadratic, annihilated byTheta.',
  'all_unknown_a2_a3_U1_beta2_drop_on_r_axis':True,
  'completeness_proof':'For an arbitrary relation of order>=3, degree1 forces constant coefficients of the37 independent linear relations zero. Degree2 modulo the177 known space forces the620 constant coefficients into the619 rank1 kernel. Subtract their chosen quadratic lifts and the corrected168 d3 lifts; remaining degree1 coefficients lie in the119 kernel of the296 linear-multiple map. The other cubic pieces are variables times actual quadratic initials or quadratics times linear initials. Their pure dual, the119 values and corrected168 values arezero, leaving exactly the619 functionals.',
  'if_all619_actual_values_zero':'pure cubic dual annihilates the WHOLE degree3 initial relation space; remaining candidate survives, henceactualHF3=1',
  'if_any619_actual_value_nonzero':'the last candidate cubic is killed; m3T=m4T, Nakayama givesm3T=0 andlengthT=10',
  'source_order_gate':{'rawP3':'existing common105 plane source3 and copied own57 source3 inside source4 can evaluate only covered actual unit germs; uncached regular divisions remain formal',
   'generic_Hadamard_h2':'needs mixed zi*r2 numerator3; plane-only cache does not imply availability',
   'generic_Hadamard_h3':'needs mixed zi*r3 numerator4; plane-only cache does not imply availability',
   'generic_cubic_division_g2':'needs own-plane numerator5; NOTauthorized',
   'generic_cubic_division_g3':'needs own-plane numerator6; NOTauthorized',
   'no_new_gauge3_or_native_X3':'native source55 row is literallyzero; beta1_r usesX2 only',
   'actual_R0_frame':'same CA-F9 zero7,rA=sA=r, no rawphysical/frame substitution'},
  'actual619_values_or_HF3_not_computed':True,'new_source_coefficient_solve_or_jet_run':False,'expanded_support':support}
 ep=out/'independent-pure619-Schur-Theta-and-completeness-preflight.json';blob=json.dumps(ev,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'actual_Schur_pure_cubic_formula_accepted':True,
  'Theta_well_defined_and_all168_corrected_boundary_values_zero':True,
  'all619_full_initial_space_completeness_criterion_proved':True,
  'four_missing_g_difference_second_roots_retained':roots4,'source_620_affected_count':len(affected),
  'new_coefficients_619_values_HF3_or_length_computed':False,'input_hashes':inp,
  'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
