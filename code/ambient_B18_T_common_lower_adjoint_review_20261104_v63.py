#!/usr/bin/env python3
"""Independent common lambda values and fixed27/177/Euler D pullback."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_cached_source3_B11_review_helpers_20261104_v63 import ordinary3,typed_B11,coefficient
from ambient_B18_comparison_X_origin_review_20261104_v63 import rational_inverse
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_COMMON_LOWER_ADJOINT_20261104_v63_REVIEW'
A='AMBIENT_B18_T_COMMON_LOWER_ADJOINT_20261104_v63'
PATHS={
 'report':f'reports/{A}/report.json','lambda_values':f'reports/{A}/phase1-actual-common-lambda-values.json',
 'D':f'reports/{A}/phase2-common-D56-by36.json','fold':f'reports/{A}/phase3-all619-common-D-contractions.json',
 'lambda_source':'reports/AMBIENT_B18_T_COMMON_LOWER_ADJOINT_SOURCE_PREFLIGHT_20261104_v63/all121-common-lower-lambda-source-functionals.json',
 'Gamma':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63_REVIEW/independent619-high-order-source-folds.json',
 'actual619':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/actual619-quartic-source-order-and-unknown-gates.json',
 'A119':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
 'star':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json',
 'm1':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'V':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
 'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'w0':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
 'gates':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'cache':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
 'cache_peer':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63_REVIEW/report.json',
 'atom1':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'copy1':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
 'copy2':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
 'fourth':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
 'native_peer':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63_REVIEW/report.json',
 'g2_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
 'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
 'helper':'code/ambient_B18_cached_source3_B11_review_helpers_20261104_v63.py'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rec(p):return {int(i):Q(c) for i,c in p if Q(c)}
def enc(p):return [[list(k) if isinstance(k,tuple) else k,str(c)] for k,c in sorted(p.items()) if c]
def put(p,k,c):
 v=p.get(k,Q(0))+c
 if v:p[k]=v
 else:p.pop(k,None)
def add(p,q,c=Q(1)):
 p=dict(p)
 for k,v in q.items():put(p,k,c*v)
 return p
def dump(out,name,p):(out/name).write_text(json.dumps(p,separators=(',',':'))+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 dump(out,'phase.json',{'phase':'common_lower_independent_inputs_saved','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items() if not p.endswith('.py')}
 assert all(D[k]['all_passed'] for k in ['report','cache_peer','native_peer','g2_peer'])
 for p,h in D['report']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['report']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 mons=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(mons)};m3=list(itertools.combinations_with_replacement(range(8),3));gam=rec(D['Gamma']['own_Gamma4'])
 ops=D['chain']['operator_scalar_nodes'];bc=list(map(dict,D['chain']['B48_by640_sparse_columns']));nodes=D['C0']['DAG_scalar_nodes'];pcols=D['C0']['honest_D0_relation_columns']
 def key(j,row,q):
  if j>=640:return (2,j-640,row,q) if row<48 else None
  if row>=48:
   i=pcols[j][row-48+1];return None if nodes[i][0]=='c' else (1,i,q)
  i=bc[j].get(row)
  if i is None or ops[i][0]=='c':return None
  return (1,ops[i][1],q) if ops[i][0]=='actual_C0_scalar_function' else (0,i,q)
 shared=[('star',i,rec(c)) for i,c in enumerate(D['star']['actual_unitstar84_original768_combinations'])]+[('L',i,rec(c)) for i,c in enumerate(D['A119']['actual_selected37_original768_constant_combinations'])]
 functionals=[];allkeys=set();p3roots=set()
 for ix,(kind,i,c) in enumerate(shared):
  fs=[{} for _ in range(8)];p3={}
  for j,a in c.items():
   if j<640 and nodes[pcols[j][8]][0]!='c':put(p3,pcols[j][8],-a)
  for index,a in gam.items():
   row,q=divmod(index,120);mon=m3[q]
   for v in set(mon):
    rest=list(mon);rest.remove(v);pair=mi[tuple(rest)]
    for j,b in c.items():
     k=key(j,row,pair)
     if k is not None:put(fs[v],k,a*b)
  expected={'shared_source_kind':kind,'index':i,'actual_original768_constant_combination':enc(c),
   'lambda_by8variables_X2_functionals':list(map(enc,fs)),'additional_rawP3_functional_at_variable_r':enc(p3)}
  assert expected==D['lambda_source']['all121'][ix]
  functionals.append((kind,i,fs,p3));allkeys.update(k for f in fs for k in f);p3roots.update(p3)
 xproots={i for r in D['gates']['normalized_Xp27_source55_rr'] for i,c in r['rawP8_root_coefficients'] if nodes[i][0]!='c'}
 raw,sourceproof=ordinary3(D,sorted(p3roots|xproots));bv,cv=typed_B11(D,{k[1] for k in allkeys if k[0]==0})
 native=[{i:rec(p) for i,p in col} for col in D['native']['source_columns128_rows48_sparse_B11_effective_second']]
 bp={(0,7):7,(1,7):8,(4,7):9,(5,7):10,(7,7):2};g2={i:Q(c) for i,c in D['g2_peer']['actual_eight_regular_g2_verified']}
 lb={};ll={};sourcevalues=[]
 for kind,i,fs,p3 in functionals:
  for var,f in enumerate(fs):
   z=Q(0);u={}
   for k,c in f.items():
    b=bp[mons[k[-1]]]
    if k[0]==2:a=native[k[1]].get(k[2],{}).get(b,Q(0));g={}
    elif k[0]==1:a=cv[k[1]].get(b,Q(0));g={}
    else:a,g=coefficient(bv[k[1]],b)
    z+=c*a;u=add(u,g,c)
   if var==7:
    for root,c in p3.items():a,g=coefficient(raw[root],3);z+=c*a;u=add(u,g,c)
   assert not u,('uncanceled_shared_source_sensitivity',kind,i,var,u)
   if z:(lb if kind=='star' else ll)[8*i+var]=z
   sourcevalues.append([kind,i,var,str(z)])
 xp0=list(map(rec,D['m1']['actual_pivot_X0_columns27_sparse56']));xp1=list(map(rec,D['m1']['actual_pivot_X1_by27_sparse56_times8']));la={}
 for index,a in gam.items():
  row,q=divmod(index,120)
  for v in set(m3[q]):
   rest=list(m3[q]);rest.remove(v);p=mi[tuple(rest)]
   for j,first in enumerate(xp1):put(la,36*j+p,a*first.get(8*row+v,Q(0)))
 for j,r in enumerate(D['gates']['normalized_Xp27_source55_rr']):
  z=Q(0);u={}
  for root,c in r['rawP8_root_coefficients']:
   if nodes[root][0]=='c':continue
   a,g=coefficient(raw[root],2);z+=Q(c)*a;u=add(u,g,Q(c))
  assert all(tag==(7,7) and root in g2 for root,tag in u)
  z+=sum(c*g2[root] for (root,tag),c in u.items());put(la,36*j+35,-z)
 LV=D['lambda_values'];assert enc(la)==LV['lambda_a2_27_times36'] and enc(lb)==LV['lambda_b1_84_times8'] and enc(ll)==LV['lambda_l1_37_times8']
 # Independent unit27 inverse, then descending triangular dual177 solve.
 free=D['X0']['free_kernel_rows32'];receipts=[r for r in D['X0']['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
 selected=[free[r['pivot']] for r in receipts];Ai=rational_inverse([[xp0[j].get(i,Q(0)) for j in range(27)] for i in selected]);zeta=[{selected[k]:a for k,a in enumerate(row) if a} for row in Ai]
 assert list(map(enc,zeta))==LV['zeta27_by56_rows']
 delta={}
 for index,a in la.items():
  j,p=divmod(index,36)
  for row,b in zeta[j].items():put(delta,36*row+p,a*b)
 def subtract(lam,cols,dual):
  lam=dict(lam)
  for j,col in enumerate(cols):
   for v in range(8):
    a=sum(c*dual.get(36*(i//8)+mi[tuple(sorted((v,i%8)))],Q(0)) for i,c in col.items());put(lam,8*j+v,-a)
  return lam
 star1=list(map(rec,D['star']['actual_unitstar84_first_sparse56_times8']));lin1=list(map(rec,D['A119']['actual_selected37_X1_sparse56_by8']))
 lb1=subtract(lb,star1,delta);ll1=subtract(ll,lin1,delta)
 assert enc(delta)==LV['delta_from_a2_source56_times36'] and enc(lb1)==LV['lambda_b1_after_a2'] and enc(ll1)==LV['lambda_l1_after_a2']
 checks=[sum(c*ll1.get(i,Q(0)) for i,c in rec(r['coefficients296_sparse']).items()) for r in D['A119']['all119_linear_syzygies_source296']];assert not any(checks)
 norm=[];cert=D['low']['quadratic_certificate']
 for j,op in enumerate(cert['fixed_elimination_receipts'][:296]):
  c={j:Q(1)}
  for p,a in op['subtract']:c=add(c,norm[p],-Q(a))
  if 'pivot_row' in op:norm.append({i:a/Q(op['divide']) for i,a in c.items()})
 basis=list(map(rec,cert['basis_normalized_columns']));alpha=[sum(c*ll1.get(i,Q(0)) for i,c in row.items()) for row in norm];mu={}
 for p,j in sorted([(min(col),j) for j,col in enumerate(basis)],reverse=True):
  a=alpha[j]-sum(c*mu.get(i,Q(0)) for i,c in basis[j].items() if i!=p)
  if a:mu[p]=a
 assert all(sum(c*mu.get(i,Q(0)) for i,c in b.items())==a for b,a in zip(basis,alpha))
 qnorm=list(map(rec,D['linear']['normalized84_unitcycle_Q1_columns_sparse']));lb2=subtract(lb1,qnorm,mu)
 U=D['m1']['free_kernel_rows32'];pi=[{U[j]:c for j,c in rec(p).items()} for p in D['m1']['pi_free32_rows']];cov=dict(delta)
 for index,a in mu.items():
  t,p=divmod(index,36)
  for row,b in pi[t].items():put(cov,36*row+p,a*b)
 wb=list(map(rec,D['w0']['complete_rank_certificate']['basis_normalized_columns']));wp={min(c):j for j,c in enumerate(wb)};kap={}
 for w,(s,pair) in enumerate(D['w0']['literal_wedge_basis84']):
  rem={w:Q(1)};coeff={}
  while rem:
   k=min(rem);j=wp[k];a=rem[k];put(coeff,j,a);rem=add(rem,wb[j],-a)
  i,j=pair;terms=[(8*s+7,i,1,1)] if j==7 else [(8*s+j,i,1,2),(8*s+i,j,-1,2)]
  for row,k,sgn,shift in terms:
   for p,mon in enumerate(mons):
    if k not in mon:continue
    rest=list(mon);rest.remove(k);v=rest[0];weight=sum(c*lb2.get(8*t+v,Q(0)) for t,c in coeff.items())
    put(kap,36*row+p,weight*Q(sgn*mon.count(k),int(v<7)+shift))
 m0=[list(map(Q,row)) for row in D['V']['all24_by56_m0']]
 for index,a in kap.items():
  row,p=divmod(index,36)
  for src,b in enumerate(m0[row]):
   if b:put(cov,36*src+p,a*b)
 assert enc(mu)==D['D']['mu177_extension180'] and enc(lb2)==D['D']['lambda_b1_after_a2_and_l1']
 assert enc(kap)==D['D']['m0_kappa2_adjoint24_times36'] and enc(cov)==D['D']['D_source56_times36'] and len(cov)==106
 # Right source coefficients only: unbound zz values stay uncomputed.
 folded=[];labels=set()
 for case in D['actual619']['all619']:
  c=rec(case['actual_original768_constant_combination_c0']);a1=rec(case['actual_Xp27_linear_a1_sparse27_times8']);f={};k0=Q(0)
  for index,a in cov.items():
   row,p=divmod(index,36);u,v=mons[p]
   for j,b in c.items():
    q=key(j,row,p)
    if q is not None:put(f,q,a*b)
   for j,first in enumerate(xp1):
    k0-=a*a1.get(8*j+u,Q(0))*first.get(8*row+v,Q(0))
    if u!=v:k0-=a*a1.get(8*j+v,Q(0))*first.get(8*row+u,Q(0))
  folded.append({'kernel619_index':case['constant619_index'],'X2_source_functional':enc(f),'known_Xp1_a1':str(k0)});labels.update(f)
 assert folded==D['fold']['all619_source_functionals'] and not D['fold']['all_required_values_fixed']
 missing=sorted({mons[i%36] for i in cov}-set(bp));assert list(map(list,missing))==D['fold']['missing_pairs']
 ev={'source_manifest':inputs,'all121_source_functionals_and_symbolic_sensitivities_checked':True,'all121_lambda_values':sourcevalues,
  'lambda_a':enc(la),'lambda_b':enc(lb),'lambda_l':enc(ll),'fixed27_inverse_method':'inverse of actual27 selected row minor, rather than primary column-by-column projector replay',
  'all119_lambda_l_kernel_values':list(map(str,checks)),'independent_D_source56_times36':enc(cov),
  'fixed177_dual_method':'descending triangular coefficient solve with free three dual coordinates zero',
  'all619_D_source_folds_and_first_correction_signs_checked':True,'new_unbound_zz_pairs':missing,
  'source3_reader_checks':sourceproof,'all_actual619_lower_values_not_computed':True,
  'no_new_source_coefficients_or_rawP4_or_X3_evaluation':True}
 ep=out/'independent-common-lambda-D-and619-source-folds.json';dump(out,ep.name,ev)
 note=ROOT/'research/notes'/f'{STEM}.md'
 report={'id':STEM,'all_passed':True,'actual_common_lambda_values_and_D106_independently_accepted':True,
  'all119_kernel_and_all619_source_fold_checks_passed':True,'actual619_lower_values_or_completequartics_or_HF4_claimed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in out.iterdir())<32<<20
 dump(out,'report.json',report);dump(out,'phase.json',{'phase':'FINAL','report_sha256':sha(out/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
