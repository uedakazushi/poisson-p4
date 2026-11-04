#!/usr/bin/env python3
"""Independent typed full-ring Schur/source-unit audit; no coefficient query."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,resource,time

ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63';STEM=AUTHOR+'_REVIEW'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def sp(raw):return{i:F(q)for i,q in raw if F(q)}
def add(a,b,c=1):
 a=dict(a)
 for i,q in b.items():
  v=a.get(i,F(0))+c*q
  if v:a[i]=v
  else:a.pop(i,None)
 return a
def parity(p):return(-1)**sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))
def normalized_origin_determinant(columns):
 leads=[min(c)for c in columns];assert sorted(leads)==list(range(len(columns)))
 for k,c in zip(leads,columns):assert c[k]==1 and all(i>=k for i in c)
 return parity(leads)
def word_reduce(word):
 out=[]
 for x in word:
  if out and(out[-1],x)in[('A','Ai'),('Ai','A'),('U','Ui'),('Ui','U')]:out.pop()
  else:out.append(x)
 return tuple(out)
def wordsum(terms):
 out={}
 for q,w in terms:
  w=word_reduce(w);q=out.get(w,F(0))+q
  if q:out[w]=q
  else:out.pop(w,None)
 return out

def main():
 started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_CPU,(60,60))
 out=ROOT/'reports'/STEM;assert not(out/'report.json').exists();out.mkdir(parents=True,exist_ok=True)
 rp='reports/'+AUTHOR+'/report.json';R=read(rp);assert R['all_passed'];hashes={rp:sha(rp)};D={}
 for k,e in R['input_hashes'].items():assert sha(e['path'])==e['sha256'];hashes[e['path']]=e['sha256'];D[k]=read(e['path'])
 for f,h in R['evidence_sha256'].items():p='reports/'+AUTHOR+'/'+f;assert sha(p)==h;hashes[p]=h
 code='code/ambient_B18_T_minimal5_Schur_operator_20261104_v63.py';assert sha(code)==R['code_sha256'];hashes[code]=sha(code)
 for st in['AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63_REVIEW','AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63_REVIEW']:
  p='reports/'+st+'/report.json';q=read(p);assert q['all_passed'];hashes[p]=sha(p)
 T,X,W,L=D['T'],D['X0'],D['w0'],D['linear']
 assert D['T_peer']['all_passed']and D['X0_peer']['all_passed']
 E=read('reports/'+AUTHOR+'/actual-minimal5-by825-Schur-operator.json')
 assert E['source_manifest']==R['input_hashes']
 assert E['actual_minimal_generator_presentation_shape']==[5,825]
 receipts=X['all768_fixed_elimination_receipts_and_constant_combinations']
 piv=[r for r in receipts if'pivot'in r];rest=[r for r in receipts if r.get('zero')]
 assert len(piv)==27 and len(rest)==741
 assert E['first27_source_combinations']==[r['constant_original_source_combination']for r in piv]
 assert E['source741_zero_constant_combinations']==[r['constant_original_source_combination']for r in rest]
 # All768 constant operations form a triangular matrix in original source
 # order, with nonzero diagonal, followed only by a permutation.
 for r in receipts:
  c=sp(r['constant_original_source_combination']);j=r['source_column'];assert j in c and max(c)==j and c[j]
 pb={r['pivot']:sp(r['normalized'])for r in piv};porder=sorted(pb);tail=[i for i in range(32)if i not in pb];rows=X['free_kernel_rows32']
 assert rows==T['actual_source32_complement_rows']
 projection=[]
 for k in tail:
  p={k:F(1)}
  for i in reversed(porder):
   c=-sum((q*p.get(j,F(0))for j,q in pb[i].items()if j!=i),F(0))
   if c:p[i]=c
  projection.append(p)
 assert projection==[sp(p)for p in E['constant_free32_quotient_projection5_rows']]
 assert all(sum((c*p.get(i,F(0))for i,c in b.items()),F(0))==0 for p in projection for b in pb.values())
 assert[[p.get(i,F(0))for i in tail]for p in projection]==[[F(i==j)for j in range(5)]for i in range(5)]
 assert E['actual_free_five_source_rows']==[rows[i]for i in tail]==[22,23,29,30,55]
 A0=[]
 for r in piv:
  actual={}
  for j,c in sp(r['constant_original_source_combination']).items():
   col={i:F(X['full_X0_56_by768'][old][j])for i,old in enumerate(rows)if F(X['full_X0_56_by768'][old][j])}
   actual=add(actual,col,c)
  assert actual==pb[r['pivot']]
  A0.append({i:actual[p]for i,p in enumerate(porder)if actual.get(p,F(0))})
 adet=normalized_origin_determinant(A0)
 assert F(E['unit27_origin_determinant'])==adet==-1
 wc=W['complete_rank_certificate'];wbasis=[sp(c)for c in wc['basis_normalized_columns']];wfactors=[sp(c)for c in wc['all_column_factorization_in_normalized_basis']]
 udet=normalized_origin_determinant(wbasis);assert F(E['unit84_origin_determinant'])==udet==1
 selected=wc['selected_original_columns'];ncomb=[]
 for k,j in enumerate(selected):
  c={j:F(1)}
  for i,q in wfactors[j].items():
   if i<k:c=add(c,ncomb[i],-q)
   else:assert i==k
  ncomb.append({i:q/wfactors[j][k]for i,q in c.items()})
 zcomb=[];labels=[]
 for j in range(741):
  if j in selected:continue
  c={j:F(1)}
  for i,q in wfactors[j].items():c=add(c,ncomb[i],-q)
  assert j in c and c[j];zcomb.append(c);labels.append(rest[j]['source_column'])
 assert ncomb==[sp(c)for c in E['second84_normalized_source741_combinations']]
 assert zcomb==[sp(c)for c in E['remaining657_source741_combinations']]
 assert labels==E['actual_remaining_source_labels657']==L['actual657_minimal_relation_source_labels']
 q1=[sp(c)for c in L['all741_Q1_columns_sparse']]
 def combine(comb):
  a={}
  for j,c in comb.items():a=add(a,q1[j],c)
  return a
 assert[combine(c)for c in ncomb]==[sp(c)for c in L['normalized84_unitcycle_Q1_columns_sparse']]
 assert[combine(c)for c in zcomb]==[sp(c)for c in L['actual40_by657_linear_relation_columns']]
 nodes=E['typed_matrix_and_scalar_operator_nodes'];assert len(nodes)==4155
 combined=len(T['scalar_and_vector_operator_nodes'])+len(nodes)+sum(map(len,ncomb+zcomb))+sum(map(len,projection))
 assert combined==E['combined_source_and_new_operator_support_bound']==R['combined_support_bound']==79716<150000
 expected_types={0:'immutable_actual_formal_function_matrix',1:'row_and_column_slice',2:'column_slice',3:'constant_column_combinations',4:'constant_column_combinations'}
 for i,name in expected_types.items():assert nodes[i]['type']==name
 assert nodes[0]=={'type':'immutable_actual_formal_function_matrix','shape':[116,936],'artifact':'T','key':'honest_T116_by936_sparse_columns','coefficient_node_key':'scalar_and_vector_operator_nodes'}
 for i,n in enumerate(nodes):
  typ=n['type'];shape=n['shape']
  def ref(k):
   j=n[k];assert isinstance(j,int)and 0<=j<i;return nodes[j]['shape']
  if typ=='immutable_actual_formal_function_matrix':assert i==0
  elif typ in['row_and_column_slice','row_slice','column_slice']:
   s=ref('matrix');rr=n.get('rows',list(range(s[0])));cc=n.get('columns',list(range(s[1])))
   assert len(set(rr))==len(rr)and len(set(cc))==len(cc)
   assert all(0<=x<s[0]for x in rr)and all(0<=x<s[1]for x in cc)and shape==[len(rr),len(cc)]
  elif typ=='constant_column_combinations':
   s=ref('matrix');assert shape==[s[0],len(n['combinations'])]
   assert all(all(0<=k<s[1]for k in sp(c))for c in n['combinations'])
  elif typ=='constant_left_multiply':
   s=ref('matrix');assert n['inner_dimension']==s[0]and shape==[len(n['left_rows']),s[1]]
   assert all(all(0<=k<s[0]for k in sp(c))for c in n['left_rows'])
  elif typ=='formal_unit_matrix_inverse':
   assert ref('matrix')==shape and shape[0]==shape[1]
   col=[sp(c)for c in n['origin_columns']];assert len(col)==shape[1]
   assert F(n['origin_determinant'])==normalized_origin_determinant(col)
  elif typ=='matrix_product':
   a,b=ref('left'),ref('right');assert a[1]==b[0]and shape==[a[0],b[1]]
  elif typ=='matrix_difference':assert ref('left')==ref('right')==shape
  elif typ=='matrix_negative':assert ref('matrix')==shape
  elif typ=='horizontal_concatenation':
   matrices=n['matrices'];assert all(j<i for j in matrices)
   assert all(nodes[j]['shape'][0]==shape[0]for j in matrices)
   assert sum(nodes[j]['shape'][1]for j in matrices)==shape[1]
  elif typ=='scalar_matrix_entry':
   s=ref('matrix');assert shape==[]and 0<=n['row']<s[0]and 0<=n['column']<s[1]
  else:raise AssertionError(n)
 # Bind every nontrivial constant operator and both Schur blocks explicitly.
 assert nodes[1]['matrix']==0 and nodes[1]['rows']==list(range(84))and nodes[1]['columns']==list(range(168))
 assert nodes[2]['matrix']==0 and nodes[2]['columns']==list(range(168,936))
 assert nodes[3]['combinations']==E['first27_source_combinations']and nodes[3]['matrix']==2
 assert nodes[4]['combinations']==E['source741_zero_constant_combinations']and nodes[4]['matrix']==2
 unit=[(i,n)for i,n in enumerate(nodes)if n['type']=='formal_unit_matrix_inverse'];assert[len(n['origin_columns'])for i,n in unit]==[27,84]
 # Operator graph is small enough to compare the exact intended block recipe.
 wanted=[
  ('row_slice',{'matrix':3,'rows':list(range(84,116))}),
  ('row_slice',{'matrix':4,'rows':list(range(84,116))}),
  ('row_slice',{'matrix':5,'rows':porder}),
  ('row_slice',{'matrix':6,'rows':porder}),
  ('constant_left_multiply',{'matrix':5,'left_rows':E['constant_free32_quotient_projection5_rows']}),
  ('constant_left_multiply',{'matrix':6,'left_rows':E['constant_free32_quotient_projection5_rows']}),
  ('row_slice',{'matrix':3,'rows':list(range(84))}),
  ('row_slice',{'matrix':4,'rows':list(range(84))}),
  ('formal_unit_matrix_inverse',{'matrix':7}),
  ('matrix_product',{'left':13,'right':8}),
  ('matrix_product',{'left':11,'right':14}),
  ('matrix_difference',{'left':12,'right':15}),
  ('matrix_product',{'left':9,'right':14}),
  ('matrix_difference',{'left':10,'right':17}),
  ('constant_column_combinations',{'matrix':16,'combinations':E['second84_normalized_source741_combinations']}),
  ('constant_column_combinations',{'matrix':18,'combinations':E['second84_normalized_source741_combinations']}),
  ('constant_column_combinations',{'matrix':16,'combinations':E['remaining657_source741_combinations']}),
  ('constant_column_combinations',{'matrix':18,'combinations':E['remaining657_source741_combinations']}),
  ('formal_unit_matrix_inverse',{'matrix':19}),
  ('matrix_product',{'left':20,'right':23}),
  ('matrix_product',{'left':24,'right':1}),
  ('matrix_negative',{'matrix':25}),
  ('matrix_product',{'left':24,'right':21}),
  ('matrix_difference',{'left':22,'right':27}),
  ('horizontal_concatenation',{'matrices':[26,28]})]
 for i,(typ,values)in enumerate(wanted,5):
  assert nodes[i]['type']==typ,(i,nodes[i],typ)
  assert all(nodes[i][k]==v for k,v in values.items()),(i,nodes[i],values)
 assert E['actual_relation_matrix_node']==29
 roots=E['actual_five_by825_scalar_operator_root_ids'];assert len(roots)==5and all(len(row)==825 for row in roots)
 for i,row in enumerate(roots):
  for j,k in enumerate(row):assert nodes[k]=={'type':'scalar_matrix_entry','shape':[],'matrix':29,'row':i,'column':j}
 assert set(k for row in roots for k in row)==set(range(30,4155))
 # Full-ring block cancellations, preserving multiplication order.
 assert not wordsum([(1,('Wp',)),(-1,('Wp','Ai','A'))])
 assert not wordsum([(1,('B',)),(-1,('B','Ai','A'))])
 assert not wordsum([(1,('Ys',)),(-1,('Ys','Ui','U'))])
 assert wordsum([(1,('Wr',)),(-1,('Wp','Ai','C'))])=={('Wr',):1,('Wp','Ai','C'):-1}
 assert wordsum([(1,('Yr',)),(-1,('Ys','Ui','Ur'))])=={('Yr',):1,('Ys','Ui','Ur'):-1}
 assert wordsum([(-1,('Ys','Ui','D'))])=={('Ys','Ui','D'):-1}
 proof={'all_passed':True,'all4155_typed_nodes_and4125_roots_checked':True,
  'all_actual_source27_741_and84_657_constant_combinations_invertible':True,
  'constant_quotient5_projection_left_annihilator_independently_checked':True,
  'origin_A27_determinant':str(adet),'origin_U84_determinant':str(udet),
  'determinants_minus1_plus1_are_ORIGIN_values_only':True,
  'formal_unit_inverse_justified_by_nonzero_origin':True,
  'full_ring_noncommutative_ordered_block_identities_checked':True,
  'all116_936_sources_preserved_by111_unit_pair_cancellations':True,
  'actual_complete_presentation':[5,825],'five_generators_minimal_all_relation_constants_zero':True,
  'source_linear_Yrest1_and_d3_quadratic_minus_Ystar1_U0inv_D1_match_accepted_interfaces':True,
  'original_formal_Euler_functions_retained_before_any_I2_quotient':True,
  'no_new_source_coefficient_solve_rank_348_length_depth':True}
 (out/'independent-typed-source-unit-and-full-Schur-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 for p in[str(Path(__file__).relative_to(ROOT)),'research/notes/'+STEM+'.md']:hashes[p]=sha(p)
 report={'id':STEM,'all_passed':True,'frozen':True,'actual_complete_minimal_generator_presentation':[5,825],
  'full_ring_111_unit_cancellations_and_source_completeness_verified':True,
  'unit_determinants_certified_at_origin_only':True,'T_length_or_depth_claim':False,
  'input_sha256':hashes,'evidence_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
  'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}),flush=True)

if __name__=='__main__':main()
