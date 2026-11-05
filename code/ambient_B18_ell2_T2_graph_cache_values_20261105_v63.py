#!/usr/bin/env python3
"""One final cache-only job: two invariants,A23,normalized linear graph."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63'
def main(out,gate):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'conditional_source_and_height_pass_cache_only_worker'})
 hs={};maxbits=0
 def check(v):
  nonlocal maxbits
  vals=v.values() if isinstance(v,dict) else v
  for c in vals:
   if isinstance(c,F):maxbits=max(maxbits,abs(c.numerator).bit_length(),c.denominator.bit_length())
  assert maxbits<=538,'source-derived538bitbound';return v
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 gr=read(gate);assert gr['all_passed'] and gr['whole_cache_slot_bound']==47182;ip=str(Path(gate).parent/'actual-combined-invariant-and-A23-cache-source-input.json');D=read(ip);assert hs[ip]==gr['evidence_sha256'][Path(ip).name]
 sp=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_INPUT_GATE_20261105_v63_REVIEW/report.json');hp=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_HEIGHT_BOUND_20261105_v63_REVIEW/report.json');assert sp['all_passed'] and hp['all_passed']
 height=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_HEIGHT_BOUND_20261105_v63/report.json');assert height['all_passed'] and height['maximum_derived_intermediate_or_output_height_bits']==538
 factorpeer=read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW/report.json');assert factorpeer['all_passed']
 for p,h in gr['input_sha256'].items():assert H.sha(ROOT/p)==h;hs[p]=h
 refs=D['existing_cache_paths'];raw=read(refs['raw33']);A=[list(map(H.parse,c)) for c in raw['coordinate_actions9_sparse_columns']];basis=raw['actual33_basis'];table={}
 for i,(s,r,e) in enumerate(basis):
  for j in range(33):
   v={j:F(1)}
   for k in ([s-1] if s else [])+[7]*r+[8]*e:v=mv(A[k],v)
   if v:table[i,j]=v
 def product(a,b):
  v={}
  for i,c in a.items():
   for j,d in b.items():
    if (i,j) in table:v=check(add(v,table[i,j],c*d))
  return v
 source=read(refs['source105']);root=dict(source['all_requested164_root_values']);root.update(zip(source['joint105_unknown_order'],source['all105_unknown_final_values']));traces=D['raw33_primitive_source_formula']['actual_monic_trace_s_t_nodes'];nu=[]
 for s,t in traces:nu.append(check(add(H.parse(root[t]),product(H.parse(root[s]),H.parse(root[s])),F(1,4))))
 phys=read(refs['physical12'])['physical12_coefficients'];p5=phys[5];p7=phys[7];assert p5['physical_row']==67 and p7['physical_row']==63
 phi=[check(add(scale(product(nu[1],H.parse(p5['bB'])),F(3375,11)),H.parse(p7['delta']),F(-50))),
  check(add(scale(product(nu[0],H.parse(p5['bA'])),F(-3375,11)),H.parse(p7['delta']),F(25,2))),
  check(add(add(scale(H.parse(p5['delta']),F(-3375,11)),H.parse(p7['bB']),F(50)),H.parse(p7['bA']),F(25,2)))]
 bars=[[H.parse(v) for v in row] for row in D['actual_two_three_component_bar_rho']];invariants=[sum((c*phi[t].get(i,F(0)) for t,bar in enumerate(row) for i,c in bar.items()),F(0)) for row in bars];check(invariants)
 H.dump(out/'phase1-actual-two-invariants-from-raw33-primitive.json',{'nu_A_B_full33':[srec(v) for v in nu],
  'primitive_V_Wv_full33':[srec(v) for v in phi],'actual_two_invariants_alpha_beta':list(map(str,invariants)),
  'all_highest_source4_handled_by_proven_functional_descent':True,'individual_new_Phi4_coefficients_not_queried':True})
 # Exact oldfirst data andnative individualfirsts, nevereffectiveXp shortcuts.
 mg=read(refs['Mugate']);L0=list(map(H.parse,mg['L0_Res4_Ntranspose_columns']));S0=list(map(H.parse,mg['native_S0_columns']));Se=list(map(H.parse,mg['native_S1ell_columns']));native=read(refs['native']);S2=list(map(H.parse,native['intrinsic_S2rell_columns']));Lr=list(map(H.parse,native['actual_L1r_columns']));Le=list(map(H.parse,native['actual_L1ell_columns']))
 g8=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json')['actualG8'];Sop=H.selected(H.PATHS['S'],'native_S12_by32');Sr=[{} for _ in range(128)]
 for row,col,k,sign in Sop['koszul28_in_first8_rows']:
  for inp in range(4):
   for mask,c in enumerate(map(F,g8[k][4:8])):
    if c and not(mask&inp):Sr[4*col+inp]=add(Sr[4*col+inp],{4*row+(mask|inp):F(sign)*c})
 nativeE=[check(add(mv(L0,Se[j]),mv(Le,S0[j]))) for j in range(128)];nativeR=[check(add(mv(L0,Sr[j]),mv(Lr,S0[j]))) for j in range(128)]
 old640=read(refs['old640'])['source_columns640_rows56_sparse_dual9'];old128=read(refs['old128'])['source_columns128_rows48_sparse_dual9'];x0=[];xr=[]
 for col in old640+old128:x0.append({row:H.parse(v)[0] for row,v in col if H.parse(v).get(0,F(0))});xr.append({row:H.parse(v)[8] for row,v in col if H.parse(v).get(8,F(0))})
 assert xr[640:]==nativeR;xe=list(map(H.parse,read(refs['old640ell'])['actual_first640_ell_columns']))+nativeE
 mini=read(refs['mini']);piv=mini['original27_X_source_columns'];Xp=[x0[j] for j in piv];B30inv=list(map(H.parse,mini['fixed_bottom30_inverse_columns']));C=mini['free30_V_indices'];M=list(map(H.parse,mg['Mnew0_columns']))
 graphsource=read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json');W=list(map(H.parse,graphsource['H0_free3_columns']))
 fold=read('reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63/actual-two-rell-source-folds-and-retained-native-first-operators.json');Gamma=list(map(H.parse,fold['two_Gamma_on56_actual_X2']));chain=read(H.PATHS['chain']);ops=chain['operator_scalar_nodes'];B=list(map(dict,chain['B48_by640_sparse_columns']));P=H.selected(H.PATHS['C0'],'honest_D0_relation_columns')
 cv=read(refs['actual_old_C0']);CV={i:list(map(F,v)) for i,v,g in cv['all104_original_C0_values_and_opaque_gradients']};assert all(not g for i,v,g in cv['all104_original_C0_values_and_opaque_gradients']);Cbasis=H.selected(H.PATHS['C0'],'actual_C0_basis8_scalar_nodes')
 for i in Cbasis:CV.setdefault(i,[F(0)]*4)
 ad=read(refs['typed_AD']);qknown={i:F(v[3]) for i,v in ad['typed_node_known4_coefficients']};grad={i:H.parse(g) for i,g in ad['typed_node_original_second_gradients']}
 af=read('reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json');ar0={i:H.parse(v).get(8,F(0)) for i,v in af['all79_atom_values']};memo={v['original_node']:H.parse(v['origin_and_ds']) for v in mg['all_exact_cached_leaf_bindings_origin_and_ds']}
 for i,n in mg['old_exact84_original_arithmetic_nodes']:
  if n[0]=='c':v={0:F(n[1],n[2])} if n[1] else {}
  elif n[0]=='a':v=add(memo[n[1]],memo[n[2]])
  elif n[0]=='n':v=scale(memo[n[1]],-1)
  else:v=H.pairmul(memo[n[1]],memo[n[2]])
  memo[i]=v
 for i in [516,785,846]:
  x,y=[ops[j][1] for j in ops[i][3]];assert not memo[x].get(0,F(0)) and not memo[y].get(0,F(0));qknown[i]=-(ar0[x]*memo[y].get(1,F(0))/5+memo[x].get(1,F(0))/5*ar0[y]);grad[i]={}
 values=[];proof=[]
 for record in D['all23_old_source_sensitivity_proofs']:
  c=H.parse(record['original768']);cn={j-640:a for j,a in c.items() if j>=640};assert not mv(S0,cn);ur=mv(xr,c);ue=mv(xe,c);ars=[]
  for u in [ur,ue]:
   assert not mv(M,u);coordinates=mv(B30inv,{i:u.get(row,F(0)) for i,row in enumerate(C) if u.get(row,F(0))});a={i:v for i,v in coordinates.items() if i<27};b={i-27:v for i,v in coordinates.items() if i>=27};assert add(mv(Xp,a),mv(W,b))==u;ars.append(a)
  ar,ae=ars;vnat=check(add(add(mv(L0,mv(S2,cn)),mv(Lr,mv(Se,cn))),mv(Le,mv(Sr,cn))))
  rowvalues=[];parts=[]
  for gamma in Gamma:
   f={};old=F(0)
   for j,a in c.items():
    if j>=640:continue
    for row,b in gamma.items():
     coefficient=a*b
     if row<48:
      i=B[j].get(row)
      if i is None:continue
      n=ops[i]
      if n[0]=='c':continue
      if n[0]=='actual_C0_scalar_function':f=add(f,{n[1]:coefficient})
      else:old+=coefficient*qknown[i];f=add(f,grad[i],coefficient)
     else:f=add(f,{P[j][row-48+1]:coefficient})
   assert all(i in CV for i in f);old+=sum((a*CV[i][3] for i,a in f.items()),F(0));nv=sum((a*vnat.get(i,F(0)) for i,a in gamma.items()),F(0));correction=F(0)
   for k,j in enumerate(piv):correction+=sum((a*xr[j].get(i,F(0)) for i,a in gamma.items()),F(0))*ae.get(k,F(0))+sum((a*xe[j].get(i,F(0)) for i,a in gamma.items()),F(0))*ar.get(k,F(0))
   z=old+nv-correction;check([old,nv,correction,z]);rowvalues.append(z);parts.append(list(map(str,[old,nv,correction])))
  values.append(rowvalues);proof.append({'basis_column':record['basis_column'],'source642':record['source642'],'a1_r':srec(ar),'a1_ell':srec(ae),'two_old_native_pivot_pieces':parts})
 normalized=[]
 for row,Aval in zip(graphsource['all23_linear_graph_source_records'],values):
  ce=F(row['four_unknown_mu_coefficients_per_channel']['M_r_H0v']);cr=F(row['four_unknown_mu_coefficients_per_channel']['M_ell_H0v']);normalized.append(check([Aval[0]+ce*invariants[0]+cr*invariants[1],Aval[1]]))
 ep=out/'actual-two-invariants-A23-and-normalized-linear-graph.json';H.dump(ep,{'actual_two_invariants_alpha_beta':list(map(str,invariants)),
  'actual23_two_A_values':[list(map(str,v)) for v in values],'actual23_normalized_graph_values':[list(map(str,v)) for v in normalized],
  'all23_actual_source_pivot_first_proofs':proof,'all23_nativeS0_zero_verified':True,
  'all23_honest_second_sensitivity_zero_source_proof_pinned':True,
  'canonical_other_two_mu_values_not_computed_or_set_zero':True,
  'normalizing_unit':'typed coherentMprime=M+G2*A0 fromGRAPH4 source; A0 definedbycanonicalmu1r,mu1ell andnonzerovcomponent. Existenceunitisproved; notnumericallyserializedcanonicalallorderpresentation.',
  'normalized_graph_is_actual_initial_space_in_an_existing_constant_unit_frame':True,
  'actual_T2_total_length_or_H0_or_K_exponents_not_claimed':True})
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'actual_two_invariants_alpha_beta':list(map(str,invariants)),
  'all46_A23_values_and_normalized_graph_source_bound':True,'canonical_four_M1_values_computed':False,
  'new_physical_implicit_gauge_or_jet_coefficients':False,'T2_total_length_H0_or_K_exponents_computed':False,
  'source_proven_height_bound_bits':538,'declared_height_guard_bits':2048,'actual_computed_max_height_bits':maxbits,
  'owned_QQ_source_bound':47182,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),
  'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{p.name:H.sha(p) for p in out.iterdir() if p.name!='phase.json'},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','actual_two_invariants_alpha_beta','actual_computed_max_height_bits','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--gate',required=True);a=p.parse_args();main(a.output.resolve(),a.gate)
