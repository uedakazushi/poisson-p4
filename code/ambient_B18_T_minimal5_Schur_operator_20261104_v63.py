#!/usr/bin/env python3
"""Complete five-target T presentation by two actual unit Schur complements.

This is a typed formal-function operator builder.  It evaluates no new source
coefficient and performs no rank search.  The fixed27/84 pivots come from the
immutable actual X0 and first-order wedge certificates.
"""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63'
PATHS={
 'T':'reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63/actual-honest-T116-936-Euler-operator.json',
 'T_peer':'reports/AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63_REVIEW/report.json',
 'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'X0_peer':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63_REVIEW/report.json',
 'w0':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
 'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
 'low_report':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/report.json'}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def sparse(x):return {int(i):Q(c) for i,c in x if Q(c)}
def enc(x):return [[i,str(c)] for i,c in sorted(x.items()) if c]
def add(a,b,c=Q(1)):
 a=dict(a)
 for i,x in b.items():
  v=a.get(i,Q(0))+c*x
  if v:a[i]=v
  else:a.pop(i,None)
 return a
def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def triangular_unit(cols,n):
 assert len(cols)==n
 leads=[]
 for col in cols:
  assert col
  p=min(col);assert col[p]==1 and all(i>=p for i in col)
  leads.append(p)
 assert sorted(leads)==list(range(n))
 return sign(leads)

def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=False)
 source={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}
 (out/'phase.json').write_text(json.dumps({'phase':'source_manifest_saved_before_operator_construction','inputs':source},indent=2)+'\n')
 D={k:read(ROOT/p) for k,p in PATHS.items()}
 assert D['T_peer']['all_passed'] and D['X0_peer']['all_passed']
 T,X,W,L=D['T'],D['X0'],D['w0'],D['linear']
 assert len(T['honest_T116_by936_sparse_columns'])==936
 assert len(T['scalar_and_vector_operator_nodes'])==74176
 assert X['rank_X0']==27 and W['rank_w0']==84 and W['mu_Rbar_T_exact']==5
 rows=X['free_kernel_rows32'];assert rows==T['actual_source32_complement_rows']
 receipts=X['all768_fixed_elimination_receipts_and_constant_combinations']
 piv=[r for r in receipts if 'pivot' in r];rest=[r for r in receipts if r.get('zero')]
 assert len(piv)==27 and len(rest)==741
 pb={r['pivot']:sparse(r['normalized']) for r in piv};porder=sorted(pb)
 tail=[i for i in range(32) if i not in pb];assert len(tail)==5
 # Q is an actual constant row unit on the free32 source coordinates.
 # Its lower five rows are the quotient projection used in the small test.
 projection=[]
 for t in tail:
  row={}
  for j in range(32):
   v={j:Q(1)}
   for p in porder:v=add(v,pb[p],-v.get(p,Q(0)))
   if v.get(t,Q(0)):row[j]=v[t]
  projection.append(row)
 assert [[p.get(j,Q(0)) for j in tail] for p in projection]==[[Q(i==j) for j in range(5)] for i in range(5)]
 assert all(not sum(row.get(i,Q(0))*c for i,c in col.items()) for row in projection for col in pb.values())
 A0cols=[{i:c for i,p in enumerate(porder) if (c:=sparse(r['normalized']).get(p,Q(0)))} for r in piv]
 # Ordered by the pivot labels, this is lower triangular with diagonal one.
 colperm=sorted(range(27),key=lambda j:piv[j]['pivot'])
 assert triangular_unit([A0cols[j] for j in colperm],27)==1
 A0det=sign(colperm)
 wc=W['complete_rank_certificate'];wops=wc['fixed_elimination_receipts']
 assert len(wops)==741 and wc['rank']==84
 W0cols=list(map(sparse,wc['basis_normalized_columns']));U0det=triangular_unit(W0cols,84)
 assert [rows[i] for i in tail]==L['five_remaining_free_source_rows']
 assert porder==L['constant_projection_pivot_rows27'] and tail==L['constant_projection_tail5']
 # Replaying only the fixed column operations checks the literal source labels
 # and all first-coefficient interfaces, without a new rank calculation.
 q1=list(map(sparse,L['all741_Q1_columns_sparse']));qn=[];qr=[];restlabels=[]
 normalized_combos=[];zero_combos=[]
 for j,op in enumerate(wops):
  assert op['column']==j
  v=dict(q1[j]);combo={j:Q(1)}
  for p,c in op['subtract']:
   v=add(v,qn[p],-Q(c));combo=add(combo,normalized_combos[p],-Q(c))
  if 'pivot_row' in op:
   c=Q(op['divide']);qn.append({i:a/c for i,a in v.items()});normalized_combos.append({i:a/c for i,a in combo.items()})
  else:qr.append(v);zero_combos.append(combo);restlabels.append(rest[j]['source_column'])
 assert list(map(enc,qn))==L['normalized84_unitcycle_Q1_columns_sparse']
 assert list(map(enc,qr))==L['actual40_by657_linear_relation_columns']
 assert restlabels==L['actual657_minimal_relation_source_labels']
 assert len(qn)==84 and len(qr)==657
 (out/'phase.json').write_text(json.dumps({'phase':'fixed_constant_interfaces_verified','A27_origin_determinant':str(A0det),'U84_origin_determinant':str(U0det),'no_rank_search':True})+'\n')
 nodes=[]
 def put(kind,shape,**kw):
  n=len(nodes);nodes.append({'type':kind,'shape':shape,**kw});return n
 original=put('immutable_actual_formal_function_matrix',[116,936],artifact='T',key='honest_T116_by936_sparse_columns',coefficient_node_key='scalar_and_vector_operator_nodes')
 d3=put('row_and_column_slice',[84,168],matrix=original,rows=list(range(84)),columns=list(range(168)))
 sourcecols=put('column_slice',[116,768],matrix=original,columns=list(range(168,936)))
 pcols=put('constant_column_combinations',[116,27],matrix=sourcecols,combinations=[r['constant_original_source_combination'] for r in piv])
 ccols=put('constant_column_combinations',[116,741],matrix=sourcecols,combinations=[r['constant_original_source_combination'] for r in rest])
 def top(m,n):return put('row_slice',[84,n],matrix=m,rows=list(range(84)))
 def free(m,n):return put('row_slice',[32,n],matrix=m,rows=list(range(84,116)))
 Pfree,Cfree=free(pcols,27),free(ccols,741)
 A=put('row_slice',[27,27],matrix=Pfree,rows=porder)
 C=put('row_slice',[27,741],matrix=Cfree,rows=porder)
 B=put('constant_left_multiply',[5,27],matrix=Pfree,left_rows=list(map(enc,projection)),inner_dimension=32)
 Y=put('constant_left_multiply',[5,741],matrix=Cfree,left_rows=list(map(enc,projection)),inner_dimension=32)
 Wp,Wr=top(pcols,27),top(ccols,741)
 Ai=put('formal_unit_matrix_inverse',[27,27],matrix=A,origin_determinant=str(A0det),origin_columns=[enc(c) for c in A0cols])
 AC=put('matrix_product',[27,741],left=Ai,right=C)
 Wh=put('matrix_difference',[84,741],left=Wr,right=put('matrix_product',[84,741],left=Wp,right=AC))
 Yh=put('matrix_difference',[5,741],left=Y,right=put('matrix_product',[5,741],left=B,right=AC))
 Us=put('constant_column_combinations',[84,84],matrix=Wh,combinations=list(map(enc,normalized_combos)))
 Ys=put('constant_column_combinations',[5,84],matrix=Yh,combinations=list(map(enc,normalized_combos)))
 Ur=put('constant_column_combinations',[84,657],matrix=Wh,combinations=list(map(enc,zero_combos)))
 Yr=put('constant_column_combinations',[5,657],matrix=Yh,combinations=list(map(enc,zero_combos)))
 Ui=put('formal_unit_matrix_inverse',[84,84],matrix=Us,origin_determinant=str(U0det),origin_columns=list(map(enc,W0cols)))
 P=put('matrix_product',[5,84],left=Ys,right=Ui)
 Rd=put('matrix_negative',[5,168],matrix=put('matrix_product',[5,168],left=P,right=d3))
 Rr=put('matrix_difference',[5,657],left=Yr,right=put('matrix_product',[5,657],left=P,right=Ur))
 Rmin=put('horizontal_concatenation',[5,825],matrices=[Rd,Rr])
 roots=[[put('scalar_matrix_entry',[],matrix=Rmin,row=i,column=j) for j in range(825)] for i in range(5)]
 coefficient_slots=sum(len(c) for c in normalized_combos+zero_combos)+sum(len(x) for x in projection)+len(nodes)+74176
 assert coefficient_slots<150000
 evidence={
  'actual_base_ring':'QQ[[zA0,...,zA6,rA]]=R0/(ell), ell=5(sA-rA)',
  'source_manifest':source,'typed_matrix_and_scalar_operator_nodes':nodes,
  'actual_minimal_generator_presentation_shape':[5,825],
  'actual_five_by825_scalar_operator_root_ids':roots,
  'actual_relation_matrix_node':Rmin,'actual_relation_order':['168 literal d3 images','657 remaining actual source cycles'],
  'actual_remaining_source_labels657':restlabels,'actual_free_five_source_rows':[rows[i] for i in tail],
  'constant_free32_quotient_projection5_rows':list(map(enc,projection)),
  'constant_free32_pivot_rows27':porder,'first27_source_combinations':[r['constant_original_source_combination'] for r in piv],
  'source741_zero_constant_combinations':[r['constant_original_source_combination'] for r in rest],
  'second84_normalized_source741_combinations':list(map(enc,normalized_combos)),
  'remaining657_source741_combinations':list(map(enc,zero_combos)),
  'unit27_origin_determinant':str(A0det),'unit84_origin_determinant':str(U0det),
  'actual_second_stage_U0_columns84':list(map(enc,W0cols)),
  'first_schur_formula':{'W_hat':'Wr-Wp*A^-1*C','Y_hat':'Y-B*A^-1*C','orders':['B in m','C in m','Y in m']},
  'second_schur_formula':{'unit':'U=W_hat*normalized84','d3':'-Ystar*U^-1*d3','source657':'Yrest-Ystar*U^-1*Wrest','orders':['Ystar in m','Yrest in m','Wrest in m']},
  'known_low_coefficients':{
   'source657_linear':'Yrest_1 (Wrest_0=0 after the exact frozen constant column operations)',
   'd3_constant_and_linear':'zero',
   'd3_quadratic':'-Ystar_1*U0^-1*d3_1',
   'Yhat_1':'Y_1 because B*A^-1*C is in m^2',
   'What_0':'Wr_0 because C0=0',
   'unknown_m1_X2_or_wedge_unit_firsts_not_needed_for_these_coefficients':True,
   'all_remaining_source_quadratics_known':False},
  'full_minimal_generator_count':5,'relation_count_not_claimed_minimal':825,
  'exact_completeness':'Two invertible source row/column operations and two full formal unit Schur complements preserve the actual cokernel of the complete116/936 source matrix.',
  'no_new_scalar_coefficient_evaluation_source_jet_solve_rank_search_or_348_matrix':True,
  'combined_source_and_new_operator_support_bound':coefficient_slots,
  'frozen_A0_linear_interfaces_replayed_without_rank_search':True,
  'total_T_length_original_depth_and_nil_cube_not_decided':True}
 p=out/'actual-minimal5-by825-Schur-operator.json';blob=json.dumps(evidence,separators=(',',':'))+'\n';assert len(blob.encode())<32<<20;p.write_text(blob)
 report={'id':STEM,'all_passed':True,'actual_complete_presentation':[5,825],'cancelled_unit_pairs':111,
  'actual_original_complete_presentation':[116,936],'unit_blocks':[27,84],'five_generators_minimal':True,
  'actual_full_cokernel_preserved':True,'known_source_linear_and_d3_quadratic_formulas_proved':True,
  'new_operator_nodes':len(nodes),'combined_support_bound':coefficient_slots,'bytes':len(blob.encode()),
  'no_new_source_or_scalar_coefficients_or_rank_search':True,'T_length_or_depth_decided':False,
  'input_hashes':source,'evidence_sha256':{p.name:sha(p)},'code_sha256':sha(Path(__file__)),
  'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'phase.json').write_text(json.dumps({'phase':'FINAL','report_sha256':sha(out/'report.json'),'actual_shape':[5,825]},indent=2)+'\n')
 print(json.dumps(report),flush=True)

if __name__=='__main__':main()
