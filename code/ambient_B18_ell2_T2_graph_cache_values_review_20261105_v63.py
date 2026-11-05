#!/usr/bin/env python3
"""Independent primitive raw33 dots and all46 source A23/normalizedgraph values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,resource,signal
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','primitive':f'reports/{A}/phase1-actual-two-invariants-from-raw33-primitive.json','values':f'reports/{A}/actual-two-invariants-A23-and-normalized-linear-graph.json',
'gate':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63/actual-combined-invariant-and-A23-cache-source-input.json',
'gate_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63_REVIEW/report.json',
'height_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_HEIGHT_BOUND_20261105_v63_REVIEW/report.json',
'factor_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW/report.json',
'fold':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json',
'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
'atom':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
'nativeN':'reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json',
'S':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
'graph':'reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json',
'gauge':'reports/AMBIENT_B18_ELL2_T2_FIRST_GRAPH_COMPARISON_GAUGE_PREFLIGHT_20261105_v63/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def dot(a,b):return sum((c*b.get(i,Q(0)) for i,c in a.items()),Q(0))
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'INDEPENDENT_EXISTING_CACHE_SOURCE_VALUES_ONLY','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k not in {'C0','S'}};assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in D if k.endswith('peer'))
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E={};refs=D['gate']['existing_cache_paths']
 for k in ['raw33','physical12','source105','actual_old_C0','typed_AD','native','Mugate','old640','old128','old640ell','mini']:
  p=refs[k];inputs[k]={'path':p,'sha256':sha(ROOT/p)};E[k]=json.loads((ROOT/p).read_text())
 R=E['raw33'];actions=[list(map(sp,M)) for M in R['coordinate_actions9_sparse_columns']];table=[]
 for s,r,e in R['actual33_basis']:
  cols=[]
  for j in range(33):
   q={j:Q(1)}
   for k in [8]*e+[7]*r+([s-1] if s else []):q=mv(actions[k],q)
   cols.append(q)
  table.append(cols)
 def product(a,b):return add(*(scale(table[i][j],c*d) for i,c in a.items() for j,d in b.items()))
 source=E['source105'];roots={i:sp(v) for i,v in source['all_requested164_root_values']}
 for i,p in zip(source['joint105_unknown_order'],source['all105_unknown_final_values']):
  if i in roots:assert roots[i]==sp(p)
  roots[i]=sp(p)
 nu=[add(roots[t],scale(product(roots[s],roots[s]),Q(1,4))) for s,t in D['gate']['raw33_primitive_source_formula']['actual_monic_trace_s_t_nodes']]
 p67,p63=E['physical12']['physical12_coefficients'][5],E['physical12']['physical12_coefficients'][7];assert p67['physical_row']==67 and p63['physical_row']==63
 phi=[add(scale(product(nu[1],sp(p67['bB'])),Q(3375,11)),scale(sp(p63['delta']),-50)),add(scale(product(nu[0],sp(p67['bA'])),Q(-3375,11)),scale(sp(p63['delta']),Q(25,2))),add(scale(sp(p67['delta']),Q(-3375,11)),scale(sp(p63['bB']),50),scale(sp(p63['bA']),Q(25,2)))]
 assert list(map(enc,nu))==D['primitive']['nu_A_B_full33'] and list(map(enc,phi))==D['primitive']['primitive_V_Wv_full33']
 bars=[[sp(p) for p in q] for q in D['gate']['actual_two_three_component_bar_rho']];invariants=[sum((dot(p,q) for p,q in zip(row,phi)),Q(0)) for row in bars]
 assert list(map(str,invariants))==D['primitive']['actual_two_invariants_alpha_beta']==D['values']['actual_two_invariants_alpha_beta']
 assert invariants==[Q(-198005215625,261453984),Q(-738247661671905875,6408579914013024)]
 dump(o,'phase.json',{'phase':'PRIMITIVE_NU_PHI_AND_TWO_DOTS_VERIFIED'})
 mu=E['Mugate'];native=E['native'];L0=list(map(sp,mu['L0_Res4_Ntranspose_columns']));S0=list(map(sp,mu['native_S0_columns']));Se=list(map(sp,mu['native_S1ell_columns']));S2=list(map(sp,native['intrinsic_S2rell_columns']));Lr=list(map(sp,native['actual_L1r_columns']));Le=list(map(sp,native['actual_L1ell_columns']));Sr=[{} for _ in range(128)]
 Sop=fields(ROOT/P['S'],['native_S12_by32'])['native_S12_by32'];gs=D['nativeN']['actualG8']
 for row,col,g,sign in Sop['koszul28_in_first8_rows']:
  for inp in range(4):
   for mask,c in enumerate(map(Q,gs[g][4:8])):
    if c and not mask&inp:Sr[4*col+inp]=add(Sr[4*col+inp],{4*row+(mask|inp):Q(sign)*c})
 nativeR=[add(mv(L0,p),mv(Lr,q)) for p,q in zip(Sr,S0)];nativeE=[add(mv(L0,p),mv(Le,q)) for p,q in zip(Se,S0)]
 allfirst=E['old640']['source_columns640_rows56_sparse_dual9']+E['old128']['source_columns128_rows48_sparse_dual9'];x0=[{i:sp(p)[0] for i,p in col if sp(p).get(0,Q(0))} for col in allfirst];xr=[{i:sp(p)[8] for i,p in col if sp(p).get(8,Q(0))} for col in allfirst];assert xr[640:]==nativeR
 xe=list(map(sp,E['old640ell']['actual_first640_ell_columns']))+nativeE;mini=E['mini'];piv=mini['original27_X_source_columns'];Xp=[x0[j] for j in piv];inv=list(map(sp,mini['fixed_bottom30_inverse_columns']));C=mini['free30_V_indices'];M=list(map(sp,mu['Mnew0_columns']));W=list(map(sp,D['graph']['H0_free3_columns']));gammas=list(map(sp,D['fold']['two_Gamma_on56_actual_X2']))
 chain=D['chain'];ops=chain['operator_scalar_nodes'];B=list(map(dict,chain['B48_by640_sparse_columns']));C0=fields(ROOT/P['C0'],['honest_D0_relation_columns','actual_C0_basis8_scalar_nodes']);Pcols=C0['honest_D0_relation_columns'];CV={i:list(map(Q,v)) for i,v,g in E['actual_old_C0']['all104_original_C0_values_and_opaque_gradients']};assert all(not g for i,v,g in E['actual_old_C0']['all104_original_C0_values_and_opaque_gradients'])
 for i in C0['actual_C0_basis8_scalar_nodes']:CV.setdefault(i,[Q(0)]*4)
 qknown={i:Q(v[3]) for i,v in E['typed_AD']['typed_node_known4_coefficients']};grad={i:sp(p) for i,p in E['typed_AD']['typed_node_original_second_gradients']};ar0={i:sp(v).get(8,Q(0)) for i,v in D['atom']['all79_atom_values']};memo={r['original_node']:sp(r['origin_and_ds']) for r in mu['all_exact_cached_leaf_bindings_origin_and_ds']}
 for i,n in mu['old_exact84_original_arithmetic_nodes']:
  if n[0]=='c':v={0:Q(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=add(memo[n[1]],memo[n[2]])
  elif n[0]=='n':v=scale(memo[n[1]],-1)
  else:
   a,b=memo[n[1]],memo[n[2]];v={0:a.get(0,Q(0))*b.get(0,Q(0)),1:a.get(0,Q(0))*b.get(1,Q(0))+a.get(1,Q(0))*b.get(0,Q(0))};v={k:c for k,c in v.items() if c}
  memo[i]=v
 for root in [516,785,846]:
  x,y=[ops[j][1] for j in ops[root][3]];assert memo[x].get(0,Q(0))==memo[y].get(0,Q(0))==0;qknown[root]=-(ar0[x]*memo[y].get(1,Q(0))+memo[x].get(1,Q(0))*ar0[y])/5;grad[root]={}
 values=[];decompositions=[]
 for number,rec in enumerate(D['gate']['all23_old_source_sensitivity_proofs']):
  c=sp(rec['original768']);cn={j-640:a for j,a in c.items() if j>=640};assert not mv(S0,cn);firsts=[mv(xr,c),mv(xe,c)];ab=[]
  for vector in firsts:
   assert not mv(M,vector);q=mv(inv,{i:vector[row] for i,row in enumerate(C) if vector.get(row,Q(0))});a={i:v for i,v in q.items() if i<27};b={i-27:v for i,v in q.items() if i>=27};assert add(mv(Xp,a),mv(W,b))==vector;ab.append((a,b))
  ar,ae=ab[0][0],ab[1][0];author=D['values']['all23_actual_source_pivot_first_proofs'][number];assert enc(ar)==author['a1_r'] and enc(ae)==author['a1_ell']
  bknown=sp(D['graph']['all23_linear_graph_source_records'][number]['known_b1_3_times9']);assert ab[0][1]=={i:bknown[i*9+7] for i in range(3) if i*9+7 in bknown} and ab[1][1]=={i:bknown[i*9+8] for i in range(3) if i*9+8 in bknown}
  vn=add(mv(L0,mv(S2,cn)),mv(Lr,mv(Se,cn)),mv(Le,mv(Sr,cn)));row=[];parts=[]
  for g in gammas:
   f={};old=Q(0)
   for j,a in c.items():
    if j>=640:continue
    for i,b in g.items():
     coefficient=a*b
     if i>=48:f=add(f,{Pcols[j][i-47]:coefficient});continue
     root=B[j].get(i)
     if root is None:continue
     n=ops[root]
     if n[0]=='c':continue
     if n[0]=='actual_C0_scalar_function':f=add(f,{n[1]:coefficient})
     else:old+=coefficient*qknown[root];f=add(f,scale(grad[root],coefficient))
   assert all(i in CV for i in f);old+=sum((a*CV[i][3] for i,a in f.items()),Q(0));nv=dot(g,vn);correction=sum((dot(g,xr[j])*ae.get(k,Q(0))+dot(g,xe[j])*ar.get(k,Q(0)) for k,j in enumerate(piv)),Q(0));row.append(old+nv-correction);parts.append(list(map(str,[old,nv,correction])))
  assert parts==author['two_old_native_pivot_pieces'];values.append(row);decompositions.append({'source642':rec['source642'],'actual_a1_r':enc(ar),'actual_a1_ell':enc(ae),'all56_first_decompositions_verified':True})
 assert list(map(lambda p:list(map(str,p)),values))==D['values']['actual23_two_A_values']
 normal=[]
 for record,row in zip(D['graph']['all23_linear_graph_source_records'],values):
  ce=Q(record['four_unknown_mu_coefficients_per_channel']['M_r_H0v']);cr=Q(record['four_unknown_mu_coefficients_per_channel']['M_ell_H0v']);normal.append([row[0]+ce*invariants[0]+cr*invariants[1],row[1]])
 assert [list(map(str,p)) for p in normal]==D['values']['actual23_normalized_graph_values'];nz=[j for j,p in enumerate(normal) if any(p)];assert nz==[19,21]
 ep=o/'independent-primitive-two-dots-all46-A-and-normalized-graph.json';dump(o,ep.name,{'inputs':inputs,'actual_two_invariants_alpha_beta':list(map(str,invariants)),'all33_monic_and_primitive_Vw_coefficients_verified':True,'all23_nativeS0_and_full56_first_pivot_decompositions':decompositions,'all46_old_native_pivot_source_pieces_verified':True,'all46_A23_values_verified':True,'all23_normalized_graph_values_verified':True,'nonzero_normalized_columns':nz,
'normalization_scope':'TheactualnormalizedlineargraphexistsunderacoherentconstantunitfromMprime=M+G2A0. Itsrecipeusescanonicalmu1r/mu1ell,whichremainwell-definedtypedscalarsnotnumericallycomputed. No canonicalfourmu wereassigned0andno allorder numericalpresentationisclaimed.',
'all_honest_second_sensitivities_canceled_before_any_unavailable_leaf_value_assignment':True,'no_new_source_implicit_gauge_or_jet_coefficients':True,'actual_T2_total_H0_or_K_exponents_not_computed':True})
 report={'id':STEM,'all_passed':True,'actual_two_invariants_alpha_beta':list(map(str,invariants)),'all46_A23_and_existing_frame_normalized_graph_independently_accepted':True,'canonical_allfour_mu_or_T2_total_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'alpha_beta':list(map(str,invariants)),'nonzero_normalized_columns':nz,'seconds':report['seconds']}))
if __name__=='__main__':main()
