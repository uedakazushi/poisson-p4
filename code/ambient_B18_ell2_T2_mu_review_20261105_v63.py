#!/usr/bin/env python3
"""Independent full9 source/X1 and supplied-Y0 rank postcheck; no solves."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,resource,signal,time
from ambient_B18_ell2_full9_X1_review_helpers_20261105_v63 import first640_full10,add,scale,encode
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,mv,mp,eye
from ambient_B18_ell2_T2_origin_bridge_review_20261105_v63 import mu_transpose
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_MU_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63'
P={
 'primary':f'reports/{A}/report.json','values':f'reports/{A}/actual741-Y0-fixed-rank-and-T2-minimal-generators.json',
 'ell640':f'reports/{A}/phase1-actual-first640-ell-from-existing-full9-cache.json',
 'gate_report':f'reports/{G}/report.json','gate':f'reports/{G}/actual-T2-mu-input-source-and-cache-gate.json',
 'resource':'reports/AMBIENT_B18_ELL2_T2_MU_RESOURCE_ADDENDUM_20261105_v63/report.json',
 'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
 'quartic':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
 'sourcevalues':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json',
 'old640':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json',
 'old128':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json',
 'origin':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'origin_bridge':'reports/AMBIENT_B18_ELL2_T2_ORIGIN_BRIDGE_20261105_v63_REVIEW/independent27units-741source-combinations-and-nativeS0zero.json',
 'origin_peer':'reports/AMBIENT_B18_ELL2_T2_ORIGIN_BRIDGE_20261105_v63_REVIEW/report.json',
 'frame':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
 'frame_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
 'active':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
 'res':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
 'G2':'reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63/actual-G2-linear-fixed-source-profile.json',
 'G2_peer':'reports/AMBIENT_B18_ELL2_C2_G2_LINEAR_SOURCE_PROFILE_20261105_v63_REVIEW/report.json',
 'adapter':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
 'S':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'native_ell':'reports/AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63/actual30-firstell-and-native-S1ell.json',
 'native_ell_peer':'reports/AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63_REVIEW/report.json',
 'DVR_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_DVR_TYPE_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):
 (o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 assert sum(p.stat().st_size for p in o.iterdir() if p.is_file())<32*1024*1024
def fields(path,wanted):
 """Parse selected TOP-level fields without retaining unused dense arrays."""
 text=Path(path).read_text();decoder=json.JSONDecoder();i=text.index('{')+1;out={}
 def ws(i):
  while text[i].isspace():i+=1
  return i
 while True:
  i=ws(i)
  if text[i]=='}':break
  if text[i]==',':i=ws(i+1)
  key,j=decoder.raw_decode(text,i);assert isinstance(key,str);i=ws(j);assert text[i]==':';i=ws(i+1)
  if key in wanted:
   out[key],i=decoder.raw_decode(text,i)
  else:
   quoted=False;escape=False;depth=0
   while i<len(text):
    c=text[i]
    if quoted:
     if escape:escape=False
     elif c=='\\':escape=True
     elif c=='"':quoted=False
    elif c=='"':quoted=True
    elif c in '[{':depth+=1
    elif c in ']}':
     if not depth:break
     depth-=1
    elif c==',' and not depth:break
    i+=1
 assert set(out)==set(wanted),(path,set(out),wanted)
 return out
def full_G2(frame,pref,res,g2):
 active=[0,1,4,5,7];zero=[2,3,6];cols=[];varmap=active+[8]
 for p in map(sp,g2['actual_linear_G2_active102_by42_columns']):cols.append({(i//6)*9+varmap[i%6]:c for i,c in p.items()})
 d1=list(map(sp,pref['active_d1_columns']));a0=list(map(sp,frame['active_F1_fixed67_columns'][:11]));eta=res['physical_target3_basis_positions'];top={i:j for j,i in enumerate(eta)};nonphys=[i for i in range(14) if i not in top]
 linear=[]
 for j in frame['active_F1_free17']:
  v={j:Q(1)}
  for k,i in enumerate(nonphys):v=add(v,scale(a0[k],-d1[j].get(i,Q(0))))
  assert not mv(d1,v);lin={}
  for row,c in v.items():
   block,i=divmod(row,14)
   if i in top:lin=add(lin,{top[i]*9+(8 if block==0 else active[block-1]):c})
  linear.append(lin)
 for iz,z in enumerate(zero):
  for j in range(17):
   v={j*9+z:Q(-1)}
   for k,c in linear[j].items():
    i,var=divmod(k,9);v=add(v,{(17+3*iz+i)*9+var:c})
   cols.append(v)
 for a,b in itertools.combinations(range(3),2):
  for j in range(3):cols.append({(17+3*b+j)*9+zero[a]:Q(1),(17+3*a+j)*9+zero[b]:Q(-1)})
 assert len(cols)==102 and sum(map(len,cols))==226
 return cols
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'INDEPENDENT_FULL9_SOURCE_THEN_SUPPLIED_Y0_RANK','inputs':inputs})
 small=set(P)-{'C0','FJ','quartic','origin'};D={k:json.loads((ROOT/P[k]).read_text()) for k in small}
 assert D['primary']['all_passed'] and D['gate_report']['all_passed'] and all(D[k]['all_passed'] for k in small if k.endswith('_peer'))
 for name,folder in [('primary',A),('gate_report',G)]:
  for p,h in D[name]['input_sha256'].items():assert sha(ROOT/p)==h
  for p,h in D[name]['evidence_sha256'].items():assert sha(ROOT/'reports'/folder/p)==h
 R=D['resource'];assert R['all_passed'] and R['primary_sha256']==inputs['primary']['sha256']
 assert R['corrected_conservative_peak_QQ_source_records']==R['original_reported_peak']+56*768==104911<136665<150000
 O=fields(ROOT/P['origin'],['full_X0_56_by768','all768_fixed_elimination_receipts_and_constant_combinations','selected_independent_columns'])
 X0=[{i:Q(O['full_X0_56_by768'][i][j]) for i in range(56) if Q(O['full_X0_56_by768'][i][j])} for j in range(768)]
 old_combos=[(r['source_column'],sp(r['constant_original_source_combination'])) for r in O['all768_fixed_elimination_receipts_and_constant_combinations'] if r.get('zero')]
 assert len(O['selected_independent_columns'])==27 and len(old_combos)==741;del O
 gate=D['gate'];cases=gate['all741_zero_origin_source_combinations'];assert [(r['source_column'],sp(r['original768'])) for r in cases]==old_combos
 assert len(cases)==741 and sum(len(c) for j,c in old_combos)==759
 assert gate['whole_slot_bound']==sum(gate['slot_breakdown'].values())==136665 and gate['limits']['QQ_coefficients']==150000
 # Selected-field input keeps the unused dense FULL9 arrays out of memory.
 C=json.loads((ROOT/P['C0']).read_text())
 J=fields(ROOT/P['FJ'],['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA','source2_monomials','both57_source2','YB_to_YA_linear'])
 Q4=fields(ROOT/P['quartic'],['monomials_bivariate_through4','families','all8_regular_function_first_gradients_in_FULL9_CA_F9','original_regular_function_nodes_unmodified','all640_exact_D0_firstjet81_columns'])
 for i,node in gate['old_exact84_original_arithmetic_nodes']:assert C['DAG_scalar_nodes'][i]==node
 assert len(gate['old_exact84_original_arithmetic_nodes'])==84
 bindings=gate['all_exact_cached_leaf_bindings_origin_and_ds'];boundary=[r['original_node'] for r in bindings]
 first,proof=first640_full10(C,D['chain'],J,Q4,D['sourcevalues'],D['old640'],{'columns_sparse':X0},boundary)
 bvals={i:sp(p) for i,p in proof['source_boundary_values_full10']}
 for r in bindings:
  p=bvals[r['original_node']];pair={}
  if p.get(0):pair[0]=p[0]
  if p.get(9):pair[1]=5*p[9]
  assert pair==sp(r['origin_and_ds'])
 ell640=[{row:p[9] for row,p in col.items() if p.get(9)} for col in first]
 assert list(map(enc,ell640))==D['ell640']['actual_first640_ell_columns']
 dump(o,'phase1-independent-full9-first640-source-proof.json',{'all_source_bindings_and84_original_nodes_exact':True,'first640':proof,'all640_actual_ell_columns_match':True,'all640_origins_and_old8_projections_checked':True})
 del C,J,Q4,bvals,first # Keep only sparse first640 derivative columns below.
 # Rebuild full9 old source derivatives from the independently checked
 # existing old8 values and the newly source-checked ell entries.
 derivatives=[]
 for j,col in enumerate(D['old640']['source_columns640_rows56_sparse_dual9']):
  p={}
  for row,v in col:
   for k,c in sp(v).items():
    if k:p[row*9+k-1]=c
  p=add(p,{row*9+8:c for row,c in ell640[j].items()});derivatives.append(p)
 L=mu_transpose(D['adapter']['N']);Li=mu_transpose(D['adapter']['N_inverse']);assert mp(L,Li)==mp(Li,L)==eye(48)
 assert list(map(enc,L))==gate['L0_Res4_Ntranspose_columns']
 S0=list(map(sp,gate['native_S0_columns']));S1=list(map(sp,D['native_ell']['actual_native_S1_ell_48_by128_sparse_columns']))
 assert gate['native_S1ell_columns']==D['native_ell']['actual_native_S1_ell_48_by128_sparse_columns']
 assert all(mv(L,p)==X0[640+j] for j,p in enumerate(S0))
 effective=[mv(L,p) for p in S1]
 for j,col in enumerate(D['old128']['source_columns128_rows48_sparse_dual9']):
  p={}
  for row,v in col:
   for k,c in sp(v).items():
    if k:p[row*9+k-1]=c
    else:assert X0[j+640].get(row,Q(0))==c
  p=add(p,{row*9+8:c for row,c in effective[j].items()});derivatives.append(p)
 assert len(derivatives)==768
 M=list(map(sp,D['frame']['actual_Mnew0_columns']));assert gate['Mnew0_columns']==D['frame']['actual_Mnew0_columns'] and all(not mv(M,p) for p in X0)
 G2=full_G2(D['frame'],D['active'],D['res'],D['G2']);assert list(map(enc,G2))==gate['full_G2_linear234_by102_columns']
 E=D['values'];selected=E['actual_selected_Y0_source_columns'];B=list(map(sp,E['actual_selected_Y0_columns102']));factors=list(map(sp,E['all741_factorizations_in_original_selected_Y0_basis']))
 assert len(B)==len(selected)==99 and len(factors)==741 and len(set(selected))==99
 selected_index={k:j for j,k in enumerate(selected)};fingerprint=hashlib.sha256();checks=0
 for k,(source,c) in enumerate(old_combos):
  assert not mv(X0,c)
  native={j-640:a for j,a in c.items() if j>=640};old={j:a for j,a in c.items() if j<640}
  assert not mv(S0,native) and not mv(X0,old)
  x1=mv(derivatives,c);rhs={}
  for i,a in x1.items():
   row,var=divmod(i,9)
   for target,b in M[row].items():rhs=add(rhs,{target*9+var:a*b})
  y=mv(B,factors[k]);assert mv(G2,y)==rhs,('ACTUAL_ALL9_SOURCE_EQUATIONS',k)
  if k in selected_index:assert factors[k]=={selected_index[k]:Q(1)} and y==B[selected_index[k]]
  fingerprint.update(json.dumps([k,enc(rhs),enc(y)],separators=(',',':')).encode());checks+=234
  if (k+1)%128==0:dump(o,'phase.json',{'phase':'SOURCE_ALL9_G2_EQUATIONS','completed741_columns':k+1})
 assert checks==741*234==E['full_G2_equations_verified'] and fingerprint.hexdigest()==E['all741_unique_Y0_source_equation_fingerprint']
 # No new rank elimination: replay exactly the supplied99 original-column
 # receipts and verify all741 factors plus the sparse original minor.
 basis={};reprs={};seq=[]
 for j,r in enumerate(E['fixed_rank_receipts']):
  assert r['basis_column']==j and r['kernel_column']==selected[j] and r['original_source_column']==cases[selected[j]]['source_column']
  v=dict(B[j]);w={j:Q(1)}
  for p,c in r['subtract']:
   assert v.get(p)==Q(c);v=add(v,scale(basis[p],-Q(c)));w=add(w,scale(reprs[p],-Q(c)))
  p=r['pivot_row'];c=Q(r['divide_by']);assert c and min(v)==p and v[p]==c and p not in basis
  basis[p]=scale(v,1/c);reprs[p]=scale(w,1/c);seq.append(p);assert mv(B,reprs[p])==basis[p]
 assert len(basis)==99
 assert [[p,enc(v)] for p,v in sorted(basis.items())]==E['normalized_basis_rows']
 assert [[p,enc(v)] for p,v in sorted(reprs.items())]==E['normalized_basis_original_selected_source_representations']
 rows=sorted(seq);ix={p:i for i,p in enumerate(rows)}
 assert rows==E['fixed_original_minor_rows'] and selected==E['fixed_original_minor_columns']
 assert [enc({ix[p]:c for p,c in v.items() if p in ix}) for v in B]==E['fixed_original_minor_sparse_columns']
 assert all(basis[p][p]==1 and not any(i<p for i in basis[p]) for p in rows)
 sign=(-1)**sum(seq[i]>seq[j] for i in range(99) for j in range(i+1,99))
 assert E['fixed_original_minor_determinant_nonzero_product']=={'sign':sign,'QQ_factors':[r['divide_by'] for r in E['fixed_rank_receipts']]}
 rank=99;mu=132-27-rank;assert 94<=rank<=102 and mu==6<=11 and E['actual_mu_B_T2']==mu
 ep=o/'independent-all9-741-source-equations-and-rank99-mu6-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all27_units_and741_constant_zero_source_combinations_bound':True,'all79_atoms_84restoration_and_full9_first640_source_reconstructed':True,'old_f63_restoration_and_new_gA_plane_kept_distinct':True,'all741_nativeS0_separately_zero':True,'all741_full9_G2_equations_checked_without_solving_Y0':checks,'all99_original_selected_Y0_columns_and741_factors_source_bound':True,'fixed99_minor_nonzero_product_verified':True,'actual_Y0_rank':rank,'actual_mu_B_T2':mu,'DVR_mu_eleven_consistency_rank_atleast94':True,'corrected_primary_source_peak':104911,'original_resource_report_addendum_pinned':True,'no_new_physical_or_implicit_solve_or_N1_or_rank_search':True,'actual_T2_length_H0_exponents_or_nilcube_not_computed':True})
 report={'id':STEM,'all_passed':True,'actual_full9_source_Y0_rank99_mu_B_T2_six_accepted':True,'all741_equations_source_combinations_minor_and_factors_accepted':True,'actual_T2_length_or_H0_or_exponents_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'helper_sha256':sha(ROOT/'code/ambient_B18_ell2_full9_X1_review_helpers_20261105_v63.py'),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep),'phase1-independent-full9-first640-source-proof.json':sha(o/'phase1-independent-full9-first640-source-proof.json')},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
