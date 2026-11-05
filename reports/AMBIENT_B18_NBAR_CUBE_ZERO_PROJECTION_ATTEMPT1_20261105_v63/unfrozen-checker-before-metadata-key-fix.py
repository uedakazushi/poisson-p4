#!/usr/bin/env python3
"""Exactthree-SchurΠ0 row plusactualquarticdual impliesallbarcubeszero."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63'
PATHS={
 'minimal':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json',
 'minimalpeer':'reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW/report.json',
 'Xfirst':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json',
 'Xpeer':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63_REVIEW/report.json',
 'm0':'reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json',
 'archive':'reports/AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63/actual512-origin-values-and-Pi0.json',
 'quartic':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63/actual-eight-relations-and-complete-quartic-dual.json',
 'quarticpeer':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63_REVIEW/report.json',
 'deep':'reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_REDUCTION_20261105_v63/report.json',
 'deeppeer':'reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_20261105_v63_REVIEW/report.json',
 'Nbarpeer':'reports/AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'finiteN':'reports/AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW/report.json',
 'finitecube':'reports/AMBIENT_B18_NIL_CUBE_FINITE_ORIGIN_20261104_v63_REVIEW/report.json',
 'conditional':'reports/AMBIENT_B18_NIL_CUBE_280_AND_ELL_DETECTION_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def field(text,key):
 token=json.dumps(key);assert text.count(token)==1;i=text.index(':',text.index(token))+1
 while text[i].isspace():i+=1
 value,j=json.JSONDecoder().raw_decode(text,i);return value,[i,j]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={p:sha(ROOT/p)for p in PATHS.values()};peak=0
 def record(phase,live):
  nonlocal peak
  peak=max(peak,count(live));assert peak<=60000,peak;(out/'phase.json').write_text(json.dumps({'phase':phase,'peak_numeric_slots':peak,'seconds':time.monotonic()-t})+'\n')
 reports={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k not in['minimal','Xfirst','m0','archive','quartic']};assert all(d['all_passed']for d in reports.values())
 text=(ROOT/PATHS['minimal']).read_text();nodes,nint=field(text,'operator_nodes');record('decoded_minimal_nodes_before_pruning',[reports,nodes]);unitmeta=[]
 for i,n in enumerate(nodes):
  if n[0]=='origin_unit_matrix_inverse':
   assert Q(n[2]['determinant']);unitmeta.append([i,n[1],n[2]['determinant']]);n[2]={'determinant':n[2]['determinant'],'origin_matrix_source_retained_by_input_hash':True}
 steps,sint=field(text,'three_exact_unit_Schur_steps');record('decoded_steps_after_node_pruning',[reports,nodes,steps]);slim=[]
 for s in steps:
  assert Q(s['unit_origin_determinant']);slim.append({k:s[k]for k in ['input','unit_rows','unit_columns','retained_rows','target_projection','unit_origin_determinant']})
 del steps,text;record('pruned_minimal_metadata',[reports,nodes,slim])
 m0=json.loads((ROOT/PATHS['m0']).read_text());E=m0['selected_source_columns24'];assert 42 not in E
 assert nodes[5]==['QQ_matrix_rows',[[str(int(i==j))for j in E]for i in range(44)]]and not any(Q(c)for c in nodes[5][1][42])
 assert nodes[8]==['matrix_product',[5,7]]and nodes[7][0:2]==['origin_unit_matrix_inverse',6]and nodes[6]==['matrix_product',[0,5]]
 assert nodes[10]==['horizontal_blocks',[1,8,9]]and nodes[9]==['zero']and nodes[13]==['vertical_blocks',[10,12]]
 # Streamoneactualsourcecolumnatatime; neverdecode780-columnwholearray.
 raw=(ROOT/PATHS['Xfirst']).read_text();token=json.dumps('all780_source_columns44_sparse_dual9');assert raw.count(token)==1;i=raw.index(':',raw.index(token))+1
 while raw[i].isspace():i+=1
 assert raw[i]=='[';i+=1;start=i;dec=json.JSONDecoder();digest=hashlib.sha256();zero_count=0;max_col=0
 for j in range(780):
  while raw[i].isspace()or raw[i]==',':i+=1
  col,i=dec.raw_decode(raw,i);max_col=max(max_col,count(col));record('stream_X0_column',[reports,nodes,slim,m0,col]);matches=[v for row,v in col if row==42];assert len(matches)<=1
  coeff=Q(0)if not matches else Q(dict(matches[0]).get(0,'0'));assert coeff==0,(j,coeff)
  digest.update(json.dumps([j,matches],separators=(',',':')).encode()+b'\n');zero_count+=1
 while raw[i].isspace():i+=1
 assert raw[i]==']';del raw,col,m0
 labels=list(range(128));route=[]
 for stage,s in enumerate(slim):
  at=labels.index(42);assert at not in s['unit_rows']and at in s['retained_rows'];new=s['retained_rows'].index(at);p=s['target_projection'];pn=nodes[p];assert pn[0]=='matrix_sum';sel,neg=pn[1:];sn=nodes[sel];assert sn[0]=='row_column_selection'and sn[2]==s['retained_rows']and nodes[sn[1]]==['identity']
  assert nodes[neg][0]=='matrix_neg';product=nodes[nodes[neg][1]];assert product[0]=='matrix_product'and len(product[1])==3;C,U,I=product[1]
  assert nodes[C]==['row_column_selection',s['input'],s['retained_rows'],s['unit_columns']]
  assert nodes[U][0]=='origin_unit_matrix_inverse';unit=nodes[nodes[U][1]];assert unit==['row_column_selection',s['input'],s['unit_rows'],s['unit_columns']]
  assert nodes[I][0]=='row_column_selection'and nodes[I][2]==s['unit_rows']and nodes[nodes[I][1]]==['identity']
  route.append({'stage':stage+1,'raw_label':42,'incoming_position':at,'retained_position':new,'unit_size':len(s['unit_rows']),'cross_C_origin_row_zero':'InitialA13row42origin0; inductionthroughpreviouszerocrossSchurrow','projection_row_is_literal_retained_selector_at_origin':True,'unit_origin_determinant':s['unit_origin_determinant']});labels=[labels[j]for j in s['retained_rows']]
 assert labels==[42,43,50,68,70]and nodes[59]==['matrix_product',[58,43,28]]and[s['target_projection']for s in slim]==[28,43,58]
 archive,ainter=field((ROOT/PATHS['archive']).read_text(),'constant_minimal_projection5_by128_sparse');assert archive[0]==[[42,'1']]
 quartic=json.loads((ROOT/PATHS['quartic']).read_text());assert quartic['actual_complete_quartic_dual_proved']and quartic['actual_complete_quartic_dual1650_sparse']==[[329,'1']]
 data={'literal_direct_n_row42_zero_allorders':True,'initial_complete_A128_by972_row42_origin_zero':True,'all780_X0_row42_coefficients_exact_zero':zero_count,'all780_row42_jet_records_sha256':digest.hexdigest(),
  'all3_exact_Schur_row_routes':route,'final_retained_original_labels5':labels,'full_minimal_Pi0_first_row128_sparse':[[42,'1']],'archive_first_row_matches_independent_source_route':True,
  'actual_complete_quartic_functional1650_sparse':[[329,'1']],'quartic_functional_meaning':'coefficientrA^4 ofFIRSTminimalgenerator, allotherdegree4coefficientsannihilated',
  'all_cycle_cube_vanishing_proof':'DEEPsourcegivesactualhonest128representative(c,Hmc)∈m4, withrawfree44coord42=c_zero7_last∈m5. Regularminimalprojectiondegree4isPi0 timesinputdegree4, andPi0firstrow=e42. Henceitsfirstminimalr4coefficientzero. Actualcompletequarticdualistheonecoordinate, som4classzero. m5Tnil=0 eliminatesallhigherterms. ThisquantifiesALLtripleproducts ofgenuineNbar, notonly35selectedinputs.',
  'proved_actual_Nbar_cube_zero':True,'no35_or560_cube_value_job_needed':True,
  'original_nilcube_consequence':'N³⊂ellN andknownN³⊂H0(D0). ForfiniteDmoduleN, H0(D0)=D-torsion(N) andH0∩ellN=ellH0: ifell*yistorsionthenyistorsion. ThusN³⊂ellH0(D0).',
  'conditional_original_exact_index3':'IfellH0(D0)=0 (acceptedcriterionT2length17), thenoriginalN³=0 andknownQbN²≠0givesindex3. ThisconditionisNOTprovedhere.',
  'unconditional_original_index_unchanged':[3,7],'original_nilcube_zero_nonzero_not_decided':True,'no_new_source_gauge_jet_or_cube_scalar_values':True,
  'source_custody':{'minimal_operator_nodes_bytes':nint,'minimal_Schur_steps_bytes':sint,'streamed_X_columns_bytes_start':start,'all780_decoded_separately':True,'maximum_single_column_numeric_slots':max_col}}
 ev=out/'all-three-Schur-row42-and-structural-bar-cube-zero.json';blob=json.dumps(data,indent=2)+'\n';assert sum((ROOT/p).stat().st_size for p in pins)+len(blob.encode())<32<<20;ev.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_structural_all_actual_bar_nil_cubes_zero','genuine_Nbar_cube_zero':True,'original_nilcube_subset_ell_H0':True,'original_nilcube_vanishing_decided':False,'unconditional_original_nil_index_interval':[3,7],
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'peak_numeric_slots':peak,'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':main()
