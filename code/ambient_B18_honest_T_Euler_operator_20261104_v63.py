#!/usr/bin/env python3
"""Actual768 cycle columns and complete formal T presentation, without Taylor evaluation."""
from fractions import Fraction as Q
import hashlib,itertools,json,resource,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_HONEST_T_EULER_OPERATOR_20261104_v63'
PATHS={
 'm':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
 'm_peer':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63_REVIEW/report.json',
 'B':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
 'B_peer':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_20261104_v63_REVIEW/report.json',
 'P':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
 'SyzQ':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
 'SyzQ_peer':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63_REVIEW/report.json',
 'resolution':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
 'resolution_peer':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def add(a,b,c=Q(1)):
 a=dict(a)
 for m,x in b.items():
  z=a.get(m,Q(0))+c*x
  if z:a[m]=z
  else:a.pop(m,None)
 return a
def mul(a,b):
 z={}
 for m,x in a.items():
  for n,y in b.items():
   e=tuple(p+q for p,q in zip(m,n));v=z.get(e,Q(0))+x*y
   if v:z[e]=v
   else:z.pop(e,None)
 return z
def var(i,n=1):return {tuple(n if j==i else 0 for j in range(8)):Q(1)}
def deriv(p,i):
 z={}
 for m,c in p.items():
  if m[i]:
   n=list(m);n[i]-=1;z[tuple(n)]=c*m[i]
 return z
def euler_inverse(p,k):return {m:c/Q(sum(m[:7])+k) for m,c in p.items()}
PAIRS=list(itertools.combinations(range(8),2))
def lift(v):
 out=[]
 for i,j in PAIRS:
  out.append(euler_inverse(deriv(v[7],i),1) if j==7 else euler_inverse(add(deriv(v[j],i),deriv(v[i],j),-1),2))
 return out
def boundary(w,n):
 g=[var(i) for i in range(7)]+[var(7,n)];v=[{} for _ in range(8)]
 for (i,j),c in zip(PAIRS,w):v[i]=add(v[i],mul(g[j],c),-1);v[j]=add(v[j],mul(g[i],c))
 return v
def main():
 t=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/STEM;out.mkdir(parents=True,exist_ok=True)
 D={k:json.loads((ROOT/p).read_text()) for k,p in PATHS.items()}
 assert all(D[k]['all_passed'] for k in ['m_peer','B_peer','SyzQ_peer','resolution_peer'])
 M,B,P,S,R=(D[k] for k in ['m','B','P','SyzQ','resolution'])
 assert M['h_is_zero'] and B['all640_UP_equals_QB_verified']
 assert S['actual_d2']['shape']==[48,128]
 assert R['ranks']==[3,24,84,168,210,168,84,24,3]
 # Sign/order tests supplement the all-degree Euler identity proof, not replace it.
 fixtures=[{(0,)*8:Q(1)}]+[var(i) for i in range(8)]+[var(7,2)]
 fixtures += [mul(var(i),var(j)) for i in range(7) for j in range(i,7)]
 tested=0
 for n in [2,3]:
  for k in range(28):
   for p in fixtures:
    w=[{} for _ in range(28)];w[k]=p;v=boundary(w,n)
    assert boundary(lift(v),n)==v;tested+=1
 # All source coefficients remain actual formal functions before any quotient.
 ops=[];cache={}
 def put(n):
  k=json.dumps(n,separators=(',',':'))
  if k in cache:return cache[k]
  i=len(ops);ops.append(n);cache[k]=i;return i
 zero=put(['zero_function']);X=[]
 for j,col in enumerate(B['B48_by640_sparse_columns']):
  x={i:put(['actual_chain_scalar_function_restrict_ell0',a]) for i,a in col}
  for i,a in enumerate(P['honest_D0_relation_columns'][j][1:]):
   if a:x[48+i]=put(['actual_C0_scalar_function_restrict_ell0',a])
  X.append(x)
 for j in range(128):
  x={i:put(['actual_Res4_SyzQ_scalar_entry_restrict_ell0',i,j]) for i in range(48)}
  X.append(x)
 assert len(X)==768
 source_manifest={
  'actual_X_shape':[56,768],'columns640':'(B48,P8), P0 removed by actual UP=QB constant-row identity',
  'columns128':'(actual_SyzQ48_by128,zero8)',
  'actual_SyzQ_entry_definition':{
   'native_value':'sum_k N[k,i_native]*S12_by32[k,j_native]',
   'R0_scalar_entry':'mu(native_value)[i_mod4,j_mod4]',
   'i_native':'floor(row/4)','j_native':'floor(column/4)',
   'source':'the immutable actual_d2 operator, including its source120 unit coefficient lowering',
   'source_signs':'exact native_S12_by32 signs retained; not replaced by another convention'},
  'source_X_sparse_columns':[[[i,a] for i,a in sorted(c.items())] for c in X],
  'actual_VX_zero':'UP=QB plus Q*SyzQ=0; target constant1 removal preserves all columns',
  'old_actual_regular_functions_restricted_before_any_quotient':True}
 phase={'operator_nodes':ops,'source_manifest':source_manifest,
        'source_hashes':{k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()}}
 (out/'phase1-actual-X768-source.json').write_text(json.dumps(phase,separators=(',',':'))+'\n')
 (out/'progress.json').write_text(json.dumps({'phase':'actual_X768_source_bound','nodes':len(ops)})+'\n')
 selected=M['fixed_selected24_columns'];complement=[j for j in range(56) if j not in selected]
 assert len(complement)==32
 templates=[]
 for sector,n in enumerate([3,3,2]):
  outputs=[]
  for i,j in PAIRS:
   outputs.append({'pair':[i,j],
    'formula':{'Euler_output_shift':1,'partial':i,'input_component':7} if j==7 else
              {'Euler_output_shift':2,'difference_of_partials':[[i,j],[j,i]]}})
  templates.append({'sector':sector,'r_power':n,
    'to_coordinates':'z=u-c*r2' if sector==1 else 'z=u',
    'from_coordinates':'u=z+c*r2' if sector==1 else 'u=z',
    'derivatives':'partial in u at fixed r','Euler':'sum_i u_i*partial_u_i',
    'outputs28':outputs})
 template=put(['formal_Euler_Koszul_template24_to84',templates])
 W=[];Tcols=[]
 # 168 actual d3 columns enter in the first84 rows of the T generators.
 for col in R['differentials8_column_sparse'][2]:
  Tcols.append({i:put(['actual_literal_Rbar_polynomial',p]) for i,p in col})
 for j,x in enumerate(X):
  cycle=put(['actual_matrix_product_column',{'matrix_artifact':'m','matrix_key':'typed_m24_by56'},j])
  lifted=put(['apply_formal_Euler_Koszul_template',template,cycle])
  w=[put(['Euler_lift_component',lifted,k]) for k in range(84)];W.append(w)
  col={i:a for i,a in enumerate(w)}
  for k,old in enumerate(complement):
   if old in x:col[84+k]=x[old]
  Tcols.append(col)
 assert len(Tcols)==936 and len(ops)<150000
 nnz=sum(map(len,Tcols));assert nnz<150000
 evidence={
  'scalar_and_vector_operator_nodes':ops,'actual_X_source':source_manifest,
  'source_Cbar_minimal_target_order':['tauA','tauB','tauA*tauB'],
  'actual_m_source':'typed_m24_by56 from immutable comparison unit packet',
  'actual_source24_selected_columns':selected,'actual_source32_complement_rows':complement,
  'exact_source_column_change':'m=[M,N], old=(M^-1*(new24-N*new32),new32); transformed X=(mX,X_complement)',
  'formal_Koszul_Euler_templates':templates,'all768_source_cycles_mX':True,
  'actual_CB_shift_coefficients7':R['CB_shift_coefficients7'],
  'w84_by768_component_operator_nodes':W,
  'honest_T_generator_order':'84 G2 sector/wedge generators, then32 actual unselected source coordinates',
  'honest_T_relation_order':'168 exact literal d3 columns, then768 actual source cycle columns',
  'honest_T116_by936_sparse_columns':[[[i,a] for i,a in sorted(c.items())] for c in Tcols],
  'all_degree_contraction_proof':{
   'cycle':'sum u_i*v_i+r^n*v8=0',
   'w_i8':'(E+1)^-1 partial_i(v8)',
   'w_ij':'(E+2)^-1(partial_i(vj)-partial_j(vi))',
   'first7':'sum_i u_i*w_ij=(E+1)^-1((E+1)vj+r^n*partial_j(v8))=vj+r^n*w_j8; d2 subtracts the final term',
   'last':'sum u_i*w_i8=v8-v8|u0=v8; r^n*v8|u0=0 in the domainQQ[[r]]',
   'output_order':'Euler acts AFTER derivatives; output degree d uses1/(d+1) or1/(d+2)',
   'integer_denominators_only':True,
   'not_Rbar_linear_contraction':'choose each of768 coefficient columns first, then extend that matrixRbar-linearly',
   'must_apply_before_I_squared_reduction':True,
   'formal_series_not_algebraic_function_claim':True},
  'fixture_count':tested,'all_fixture_boundaries_recovered':True,
  'completeness_proof':'unit column change identifies kerV with im(d2) direct_sum Rbar32; d3 is complete kernel of d2; Euler selects a true preimage for everyactual mX, so the displayed936 columns give the actual quotientT',
  'operator_coefficient_evaluation_or_T_length_computed':False,
  'X0_X1_or_minimalization_computed':False,
  'no_new_source_solve_jet_or_348_scalar_matrix':True,
  'input_hashes':phase['source_hashes']}
 text=json.dumps(evidence,separators=(',',':'))+'\n';assert len(text.encode())<32<<20
 ep=out/'actual-honest-T116-936-Euler-operator.json';ep.write_text(text)
 report={'all_passed':True,'actual_complete_T_typed_presentation':[116,936],
  'all768_source_cycles_and_lifts_bound':True,'correct_output_Euler_shifts':[1,2],
  'operator_nodes':len(ops),'T_operator_entries':nnz,'bytes':len(text.encode()),
  'fixtures':tested,'actual_T_length_kernel_torsion_or_depth_computed':False,
  'input_hashes':phase['source_hashes'],'evidence':{'path':str(ep.relative_to(ROOT)),'sha256':sha(ep)},
  'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-t,
  'rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
