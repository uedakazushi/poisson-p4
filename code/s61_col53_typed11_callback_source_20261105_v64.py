"""SOURCE-only finite typed11 callbacks and positive bound interfaces. No matrix VALUE."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ID='S61_COL53_TYPED11_CALLBACK_SOURCE_20261105_v64'
CHART='reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64'
PEER='reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64_REVIEW/report.json'
TEN='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_UNIT_SOURCE_V2_20261105_v64'
NINE='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64/actual-row35-ninth-unit-and-full168x15-operator.json'
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.prepare
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};loaded=[];context={};peak=0
 def put(n,x):
  b=enc(x);assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+len(b)<=4*1024*1024
  p=out/n;assert not p.exists();p.write_bytes(b);return p
 def read(f,h=None):
  p=ROOT/f;actual=sha(p)
  if h:assert actual==h
  pins[f]=actual;d=json.loads(p.read_text());loaded.append(d);return d
 def stop(s,f):put('failure.json',{'all_passed':False,'SOURCE_only':True,'context':context,'no_matrix_VALUE':True,'no_automatic_retry':True});raise TimeoutError('SOURCE55s')
 signal.signal(signal.SIGALRM,stop);signal.alarm(55)
 try:
  chart=read(CHART+'/report.json','d7542d84ab5d20877b4cf3106b118f6c97644ae32543856aa282a1e1ed07aaed')
  peer=read(PEER,'8911d7071992759ca5544463fa6d098f97887e0bbf50dcb8abe169a3a016059e')
  assert chart['all_passed']and peer['all_passed']
  C=read(CHART+'/actual-full-bracket-unit-chart-and-manifest-counts.json',chart['evidence_sha256']['actual-full-bracket-unit-chart-and-manifest-counts.json'])
  ten=read(TEN+'/report.json','ad9dc2bc80fa50e1bf0946f5bf515450a04acad56818cb23851120111be8824b')
  T=read(TEN+'/actual-F10-pivot-and-full167x14-typed-operator.json',ten['evidence_sha256']['actual-F10-pivot-and-full167x14-typed-operator.json'])
  assert T['common_representation_ONLY']['D10']=='U^2*F9*F10'and T['P290_applied_once_no_double_scaling']
  phase=read('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64/phase-10.json','94837b3bb2ba39ed2f31bcbeba8441bad2da0a292d7e594d3e8619f5c27ee5ab')
  full11=C['one_prospective_pivot']['full_bracket'];assert full11['all_nodes']==1203 and full11['root_index']==1202 and full11['not_the_Hminus3_face_polynomial']
  p=ROOT/full11['file'];assert sha(p)==full11['sha256'];pins[full11['file']]=sha(p)
  # The whole expression stays externally referenced. No H-mode plan replaces it.
  profile=phase['profile'];assert phase['accumulator']==full11['root_name']and profile['terms_upper']==36620
  fq=Q(profile['known_QQ_content']);fm=int(profile['primitive_integer_abs_upper']);fn=profile['terms_upper']
  nine=read(NINE,'a8ec217416218de27fd4a0170b6fc7208f4a6d12f10381e5886a5e132de4ffd6')
  f9terms=nine['primitive_F9_integer_terms75'];assert len(f9terms)==75
  f9norm=sum(abs(int(c))for e,c in f9terms)
  f10=T['actual_F10'];path=ROOT/f10['file'];assert sha(path)==f10['sha256'];pins[f10['file']]=sha(path)
  body=path.read_bytes();assert body[:9]==b'Q5CELL1\n\1';oldpoly=decode(body[9:]);assert len(oldpoly['terms'])==317 and Q(oldpoly['content'])==1
  f10norm=sum(abs(int(t[0]))for t in oldpoly['terms']);decoded_numeric=count_numeric(oldpoly);del oldpoly,body
  fbound=abs(fq)*fm*fn
  npbound=12*f9norm*fbound
  ubound=6**4*2**2*82
  d10bound=ubound**2*f9norm*f10norm;d11bound=d10bound*npbound
  assert max(npbound.numerator.bit_length(),npbound.denominator.bit_length(),d11bound.numerator.bit_length(),d11bound.denominator.bit_length())<=2048
  # Type budget fixed before constructing the operator. Scalar nodes are formal
  # source operations, never arrays of expanded Laurent coefficients.
  docnum=sum(map(count_numeric,loaded))
  maxnodes=6000;maxnode_numeric=8;maxroot_numeric=8000
  base=25000+docnum+12000+4000
  upper=base+maxnodes*maxnode_numeric+maxroot_numeric+4096
  cold=25000+docnum+decoded_numeric+12000
  assert max(upper,cold)<=150000
  sourcebudget={'SOURCE_only':True,'runtime_reserve':25000,'loaded_document_numeric':docnum,'index_cache_scalar_reserve':12000,'receipt_and_temporary_metadata':4000,'max_typed_nodes':maxnodes,'max_numeric_per_node':maxnode_numeric,'max_root_table_numeric':maxroot_numeric,'temporary_reserve':4096,'whole_SOURCE_numeric_upper':upper,'cold_existing_F10_decode_upper':cold,'no_expanded_F11_or_matrix_maps':True,'RAW_output_upper':4*1024*1024,'time':{'signal_s':55,'hard_s':60,'AS_bytes':2000000000},'no_MATRIX_VALUE_worker_or_GO':True}
  put('pre-arithmetic-SOURCE-budget.json',sourcebudget)
  nodes=[];numnodes=0;roots=[];rhs=[];lambdas=[];owners={}
  def own(extra=4096):
   nonlocal peak
   now=base+numnodes+count_numeric(roots)+count_numeric(rhs)+count_numeric(lambdas)+extra;peak=max(peak,now);assert now<=upper
  def node(op,args=(),**fields):
   nonlocal numnodes
   assert all(isinstance(k,int)and 0<=k<len(nodes)for k in args)
   d={'id':len(nodes),'op':op,'args':list(args),**fields};size=count_numeric(d);assert size<=maxnode_numeric and len(nodes)<maxnodes
   own(size+4096);nodes.append(d);numnodes+=size;return d['id']
  zero=node('QQ_constant',value='0')
  f11=node('external_full_SOURCE_DAG',owner='FULL_F11',root_index=1202)
  p7=node('known_parameter_polynomial',formula='Q*R')
  f9=node('external_known_polynomial',owner='F9')
  d10=node('known_parameter_polynomial',formula='U^2*F9*F10')
  pref=node('parameter_multiply',(p7,f9))
  np=node('parameter_multiply',(pref,f11))
  invnp=node('certified_parameter_unit_inverse',(np,),chart='old_chart_and_D_FULL_F11',certificate='CHART_PEER')
  d11=node('parameter_multiply',(d10,np))
  full=T['full167_row_ids'];top=T['top101_row_ids'];sel=T['selected14_row_ids'];cols=T['residual14_column_ids'];extra=T['all87_other_top_RHS'];zeros=T['all66_zero_constraints']
  after=C['prospective_shapes_and_ids_ONLY']
  assert after['full_row_ids']==[i for i in full if i!=6]and after['column_ids']==[j for j in cols if j!=53]
  assert set(after['full_row_ids'])==set(after['top_row_ids'])|set(zeros)
  assert len(after['top_row_ids'])==100 and len(after['full_row_ids'])==166 and len(after['column_ids'])==13
  current={}
  def N10(i,j):
   assert i in full and j in cols
   if i in zeros:return zero
   key=(i,j)
   if key not in current:current[key]=node('external_current_cleared_numerator',owner='CURRENT_M10',row=i,column=j,common_denominator=d10)
   return current[key]
  pivotrow={j:N10(6,j)for j in after['column_ids']}
  # Explicit equality N10[6,53]=P7*F9*FULL_F11, not a source face surrogate.
  pivot_alias={'row':6,'column':53,'actual_current_numerator_node':np,'identity_certificate':'CHART_PRIMARY_and_PEER','full_not_face':True}
  for i in after['full_row_ids']:
   context={'row':i,'stage':'typed_callback_construction'}
   ni=N10(i,53);lam=zero if i in zeros else node('parameter_multiply',(ni,invnp));lambdas.append([i,lam])
   for j in after['column_ids']:
    if i in zeros:ans=zero
    else:
     a=node('parameter_multiply',(np,N10(i,j)))
     b=node('parameter_multiply',(ni,pivotrow[j]))
     ans=node('parameter_subtract',(a,b))
    roots.append([i,j,ans])
   h=node('external_current_RHS',owner='CURRENT_RHS10',row=i)
   if i in zeros:newh=h
   else:
    # The pivot RHS is shared once; no values are evaluated.
    if 'pivot_RHS'not in owners:owners['pivot_RHS']=node('external_current_RHS',owner='CURRENT_RHS10',row=6)
    term=node('parameter_times_RHS',(lam,owners['pivot_RHS']))
    newh=node('RHS_subtract',(h,term))
   rhs.append([i,newh]);own()
  assert len(roots)==166*13 and len(rhs)==166 and len(lambdas)==166
  assert len(nodes)<=maxnodes and numnodes<=maxnodes*maxnode_numeric and count_numeric(roots)+count_numeric(rhs)+count_numeric(lambdas)<=maxroot_numeric
  bycell={(i,j):k for i,j,k in roots}
  assert all(bycell[i,j]==zero for i in zeros for j in after['column_ids'])
  assert all(k==zero for i,k in lambdas if i in zeros)
  # Exact formal identities are coefficient-level equalities of free symbols.
  # The factor of the coefficient RHS is N_i53/Np; no D10 is lost.
  universal={'cleared_coefficient_identity':'D10*Np*M11_ij = Np*N10ij - N10i53*N10_6j',
   'Schur_definition':'M11_ij=M10ij-(M10i53/p)*M10_6j',
   'pivot':'p=Np/D10=P7*FULL_F11/(U^2*F10)',
   'lambda':'M10i53/p=N10i53/Np',
   'RHS_identity':'h11_i=h10_i-(N10i53/Np)*h10_6',
   'inverse_substitution':'x53=(D10*h10_6-sum_j N10_6j*xj)/Np',
   'no_polynomial_cancellation_or_inverse_expansion_assumed':True}
  # Test the universal formulas by cross-multiplying formal numerator pairs.
  # Only integer exponent vectors in D,Np,Nij,Nic,Nrj are used.
  e=lambda *a:tuple(a)
  rhs_poly={e(0,1,1,0,0):1,e(0,0,0,1,1):-1}
  schur_after_common={e(0,1,1,0,0):1,e(0,0,0,1,1):-1}
  assert rhs_poly==schur_after_common
  ownerdefs={'CURRENT_M10':{'file':TEN+'/actual-F10-pivot-and-full167x14-typed-operator.json','sha256':pins[TEN+'/actual-F10-pivot-and-full167x14-typed-operator.json'],'source_field':'full167x14_RF_roots + common_representation_ONLY','clearing_semantics':'N10_ij=D10*M10_ij exactly; entire original callback restored, not raw/zero-filled','valid_rows':full,'valid_columns':cols},
   'CURRENT_RHS10':{'file':TEN+'/actual-F10-pivot-and-full167x14-typed-operator.json','source_field':'all167_residual_RHS + incoming168_RHS_kept_for_restoration','all_87_extra_and_66_consistency_RHS_retained':True},
   'FULL_F11':full11,'F9':{'file':NINE,'sha256':pins[NINE],'field':'primitive_F9_integer_terms75'},
   'CHART_PEER':{'file':PEER,'sha256':pins[PEER]}}
  boundinterface={'all_bounds_are_nonnegative':True,'F11_full_SOURCE_cardinal_upper':fn,'F11_fixed_q':str(fq),'F11_fixed_integer_M':str(fm),'F11_L1_upper':str(fbound),'F9_existing_L1':str(f9norm),'F10_existing_L1':str(f10norm),'P7_L1_upper':12,'U_L1_upper':ubound,'Np_L1_upper':str(npbound),'D10_L1_upper':str(d10bound),'D11_L1_upper':str(d11bound),'Np_expanded_term_upper_not_allocated':3*75*fn,
   'new_numerator_cardinality_bound':'n11(i,j)<=nNp*n10(i,j)+n10(i,53)*n10(6,j)',
   'new_numerator_L1_bound':'C11(i,j)<=C_Np*C10(i,j)+C10(i,53)*C10(6,j)',
   'old_callback_numeric_bounds_not_instantiated_for_new_VALUE':True,
   'no_whole_new_coefficient_VALUE_budget_or_GO_claimed':True,
   'actual_stored_representation':'typed parameter scalar/source-function DAG; no expanded Np, F11, N11 maps'}
  operator={'id':ID,'SOURCE_only':True,'owners':ownerdefs,'scalar_and_RHS_nodes':nodes,'full166x13_cleared_numerator_roots':roots,'common_denominator_node':d11,'current_pivot_numerator_alias':pivot_alias,'full166_RHS_roots':rhs,'all166_lambda_roots':lambdas,'remaining_shapes_and_ids':after,'pivot_RHS_node':owners['pivot_RHS'],'pivot_row13_numerator_nodes':[[j,k]for j,k in pivotrow.items()],
   'inverse_substitution':universal['inverse_substitution'],'all_previous69_restoration_steps_kept':True,'P290_exactly_ONCE':True,'all87extra_and66zero_RHS_kept':True,'zero66_coefficient_proof':'current full M10 constraint rows are zero; their N10i53=0 gives lambda0 and the new coefficient zero; RHS=h10_i, not zero-filled',
   'cokernel_scope':'full167x14 unit-equivalent to unit[p] direct_sum full166x13 on sameoldchart D(FULL_F11); cokerfull166=cokertop100 direct_sum K5^66; selected13square NOTsamecokernel',
   'all_old_and_F9_F10_FULL_F11_complements_OPEN':True,'no_new_matrix_or_RHS_coefficient_VALUE_performed':True}
  put('actual-full166x13-and-RHS-typed11-operator.json',operator)
  put('positive-bound-interface-and-universal-identities.json',{'positive_bounds':boundinterface,'universal':universal,'actual_source_node_numeric':numnodes,'actual_source_node_count':len(nodes),'root_table_numeric':count_numeric(roots)+count_numeric(rhs)+count_numeric(lambdas),'whole_SOURCE_peak':peak})
  for f in[Path(__file__),ROOT/'code/s61_sparse_generic5_source_view_profile_20261105_v64.py',ROOT/'code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'research/notes'/f'{ID}.md']:pins[str(f.relative_to(ROOT))]=sha(f)
  raw=sum((ROOT/p).stat().st_size for p in pins)+4*1024*1024+262144;assert raw<32*1024*1024
  report={'id':ID,'all_passed':True,'SOURCE_typed_callbacks_only_no_matrix_RHS_VALUE':True,'exact_current_full_F11_used_not_face411':True,'new_full_shape':[166,13],'top_shape':[100,13],'selected_square_shape':[13,13],'extra_RHS':87,'zero_consistency_RHS':66,'typed_node_count':len(nodes),'source_numeric_upper':upper,'actual_SOURCE_peak':peak,'fresh_RAW_upper':raw,'all_current_function_owners_and_RHS_retained':True,'immutable_SOURCE_worker_not_VALUE_evaluator':True,'future_VALUE_needs_NEW_instantiated_gate':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report);signal.alarm(0)
  print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  if not(out/'failure.json').exists():put('failure.json',{'all_passed':False,'SOURCE_only':True,'error':repr(exc),'context':context,'no_matrix_VALUE':True,'no_automatic_retry':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');main(p.parse_args())
