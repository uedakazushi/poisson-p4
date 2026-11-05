"""SOURCE preflight of ONE prospective minor coefficient. No new signed values."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,re,resource,time
ROOT=Path(__file__).resolve().parents[1]
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ID='S61_COL53_NEXT17_79_FUNCTIONAL_SOURCE_20261105_v64'
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
BASE='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_STREAM49_ACTIVE_FRAGMENT_MIXED_BINARY_COMPACTION_20261105_v64'
SIDE='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.prepare
 resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));start=time.monotonic()
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};documents={};numeric=0
 def put(n,d):
  b=enc(d);assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+len(b)<4*1024*1024;(out/n).write_bytes(b)
 def read(f):
  nonlocal numeric
  if f not in documents:
   p=ROOT/f;pins[f]=sha(p);d=json.loads(p.read_text());numeric+=count_numeric(d);assert 25000+numeric+40000<=150000;documents[f]=d
  return documents[f]
 empty=lambda:{'corner':None,'n':0,'norm':Q(0)}
 def qm(a,b):
  assert max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length())<=2048
  return a*b
 def qa(a,b):
  assert max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1<=2048
  assert a.denominator.bit_length()+b.denominator.bit_length()<=2048
  return a+b
 def leaf(p,scale=Q(1),shift=(0,0,0,0,0)):
  q=qm(Q(p['content']),scale)
  if not p['terms']or not q:return empty()
  ex=[tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5))for t in p['terms']]
  return{'corner':min(ex),'n':len(ex),'norm':qm(abs(q),Q(sum(abs(int(t[0]))for t in p['terms'])))}
 def cp(A,scale=Q(1),shift=(0,0,0,0,0)):
  if not A['n']or not scale:return empty()
  return{'corner':tuple(A['corner'][j]+shift[j]for j in range(5)),'n':A['n'],'norm':qm(A['norm'],abs(scale))}
 def mul(A,B):
  if not A['n']or not B['n']:return empty()
  return{'corner':tuple(A['corner'][j]+B['corner'][j]for j in range(5)),'n':A['n']*B['n'],'norm':qm(A['norm'],B['norm'])}
 def add(A,B):
  if not A['n']:return B
  if not B['n']:return A
  return{'corner':min(A['corner'],B['corner']),'n':A['n']+B['n'],'norm':qa(A['norm'],B['norm'])}
 def literal(scope,label):
  d=read(OLD+'/'+scope);p=(d['common_operands']if'common_operands'in d else d['row_operands'])[label]
  return leaf(p['original_literal'],Q(p['QQ_transform_scalar']),p['transform_monomial_exponent5'])
 try:
  chart=read('reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64/report.json')
  assert pins['reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64/report.json']=='d7542d84ab5d20877b4cf3106b118f6c97644ae32543856aa282a1e1ed07aaed'and chart['all_passed']
  E=read('reports/S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64/actual-full-bracket-unit-chart-and-manifest-counts.json')
  assert 17 in E['prospective_shapes_and_ids_ONLY']['selected_row_ids']and 79 in E['prospective_shapes_and_ids_ONLY']['column_ids']
  cost=read('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64/remaining15-literal-source-cost-comparison.json')
  remaining=[x for x in cost['all15_fixed_source_costs']if x['column']not in(56,53)]
  assert min(remaining,key=lambda x:x['next15_literal_source_tuple_majorant'])['column']==79
  api=read('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64/actual15-col56-and-once-shared-row35.json')
  I={x['row']:x for x in api['actual15_ninth_col56_numerators']}
  rrows=E['prospective_shapes_and_ids_ONLY']['selected_row_ids']
  assert min(rrows,key=lambda i:I[i]['terms'])==17 and I[17]['terms']==1319
  # Index is streamed one literal descriptor at a time; no complete index tree.
  base_report=read(BASE+'/report.json');assert pins[BASE+'/report.json']=='62a95473e52e2a85054586bd63e28b79ae58850bebf49d6e5a0a4daeb69e7970'and base_report['all_passed']
  ip=ROOT/BASE/'literal-file-index.json';pins[str(ip.relative_to(ROOT))]=sha(ip);assert pins[str(ip.relative_to(ROOT))]==base_report['current_index_sha256'];text=ip.read_text()
  marker='"all_literal_cells":';at=text.index(marker)+len(marker);assert text[at]=='[';at+=1
  dec=json.JSONDecoder();required={(r,c)for r in[6,17,58,31,83,14,35,104,51,11,44,52]for c in[49,51,53,79]};thin={}
  while True:
   while text[at]in' \n\r\t,':at+=1
   if text[at]==']':break
   d,end=dec.raw_decode(text,at);at=end
   key=(d['row'],d['column'])
   if key in required:assert key not in thin;thin[key]=d
  rowids=dec.raw_decode(text,text.index('"all_original_row_ids":')+len('"all_original_row_ids":'))[0]
  colids=dec.raw_decode(text,text.index('"all_original_column_ids":')+len('"all_original_column_ids":'))[0]
  assert all(r in rowids and c in colids for r,c in required);del text
  base=ROOT/BASE/'current-matrix.bin';pins[str(base.relative_to(ROOT))]=sha(base)
  assert pins[str(base.relative_to(ROOT))]=='a75785cd3bd3509e1ce56fe0a1c720b964ca4b81f36a316272335b4d8970e951'
  raw={};descriptor=[];pack=out/'selected48-existing-cell-payloads.bin';offset=0
  with base.open('rb')as f,pack.open('wb')as w:
   for r,c in sorted(required):
    if(r,c)not in thin:raw[r,c]=empty();descriptor.append({'row':r,'column':c,'known_existing_zero_by_complete_index':True});continue
    d=thin[r,c];f.seek(d['payload_offset']);body=f.read(d['payload_bytes']);assert hashlib.sha256(body).hexdigest()==d['payload_sha256']
    p=decode(body);assert len(p['terms'])==d['terms'];assert 25000+numeric+40000+count_numeric(p)+6*len(p['terms'])<150000
    raw[r,c]=leaf(p);w.write(body);descriptor.append({'row':r,'column':c,'old':d,'new_offset':offset,'bytes':len(body),'payload_sha256':d['payload_sha256']});offset+=len(body)
  prog=read(OLD+'/staged-literal-col56-program-and-positive-profiles.json')
  ops={x['output']:x for x in prog['literal_ops']};cache={};needed=set()
  def visit(name):
   if name in needed:return
   needed.add(name);op=ops[name]
   for k in['A','B']:
    if k in op:visit(op[k])
  for row in[6,17,58]:visit('row'+str(row)+'_N8_sum_m')
  visit(prog['literal_shared_N8_35_56_node'])
  for op in prog['literal_ops']:
   name=op['output']
   if name not in needed:continue
   kind=op['op']
   if kind=='input':p=literal(op['binding']['scope'],op['binding']['label'])
   elif kind=='copy':p=cp(cache[op['A']],Q(op.get('scalar','1')),op.get('shift5',[0]*5))
   elif kind=='multiply':p=cp(mul(cache[op['A']],cache[op['B']]),Q(op.get('scalar','1')),op.get('shift5',[0]*5))
   else:assert kind=='add';p=add(cache[op['A']],cp(cache[op['B']],Q(op.get('scalarB','1'))))
   old=prog['positive_support_and_height_profiles'][name]
   if not Q(old['known_QQ_content']):p=empty()
   elif p['n']:
    p['n']=old['terms_upper'];p['norm']=qm(qm(abs(Q(old['known_QQ_content'])),Q(int(old['primitive_integer_abs_upper']))),Q(old['terms_upper']))
   cache[name]=p
  c8=read('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64/actual-cleared-col82-and-commonQ3RE-representation.json')
  B={x['row']:x for x in c8['actual16_col82_numerator_values']}
  def binary(stem,rec):
   p=ROOT/stem/rec['file'];pins[str(p.relative_to(ROOT))]=sha(p);assert pins[str(p.relative_to(ROOT))]==rec['sha256'];b=p.read_bytes();assert b[:9]==b'Q5CELL1\n\1';return leaf(decode(b[9:]))
  common={x:literal('common-source-operands.json',x)for x in['Q','P7','E']}
  for j,name in enumerate(['H','L','s','q','t']):common[name]={'corner':tuple(int(i==j)for i in range(5)),'n':1,'norm':Q(1)}
  p58=binary('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64',I[58])
  C10=Q(353894400000);C9C10=Q(62208000000)*C10
  common['UF10']=cp(p58,C10,[1,0,0,0,1])
  bb58=binary('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64',B[58])
  side=read(SIDE+'/exact-11-row-coefficient-operator-and-callback-owners.json')
  variants={}
  for row in[6,17]:
   V=dict(common)
   for tag in['35','58','i']:
    rr=row if tag=='i'else int(tag)
    for name in['U','V','C','F','N','M']:V[name+tag]=literal('rows/row-'+str(rr)+'-source-operands.json',name+str(rr))
   for rr in['i','58','31','83','14','35','104']:
    actual=row if rr=='i'else int(rr)
    for col in[49,51]:V['A'+rr+'_'+str(col)]=raw[actual,col]
   n9=binary('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64',I[row])
   bi=binary('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64',B[row])
   V['C10_Ni56']=cp(n9,C10,[1,0,0,0,1])
   V['C9_C10_alpha_i']=cp(add(mul(bb58,cache['row'+str(row)+'_N8_sum_m']),mul(bi,cache['row58_N8_sum_m'])),C9C10,[1,0,1,1,3])
   for col in[53,79]:
    total=empty()
    for rawrow,poly in zip(side['raw_row_ids11'],side['raw_coefficient_side11']):
     rr=row if rawrow=='i'else int(rawrow);a=raw[rr,col]
     for mon,coef in poly:
      term=cp(a,Q(coef))
      for name,power in mon:
       if power<0:
        assert name in['H','s','q','t'];term=cp(term,Q(1),[power*x for x in V[name]['corner']])
       else:
        for _ in range(power):term=mul(term,V[name])
      total=add(total,term)
    variants[str(row)+','+str(col)]=total
  D=variants['17,79'];assert D['corner']is not None
  witness=[-3,0,6,-1,10];target=[witness[j]+D['corner'][j]for j in range(5)]
  A_source_corner=[-3,0,4,-1,10]
  diagcorner=tuple(A_source_corner[j]+D['corner'][j]for j in range(5))
  Bp,Cp=variants['17,53'],variants['6,79'];crosscorner=tuple(Bp['corner'][j]+Cp['corner'][j]for j in range(5));minorcorner=min(diagcorner,crosscorner)
  delta=[target[j]-minorcorner[j]for j in range(5)]
  serial=lambda p:{'positive_lex_corner5':p['corner'],'term_upper_not_allocated':p['n'],'L1_upper':str(p['norm'])}
  sourceowned=max(25000+numeric+40000,25000+numeric+16*len(cache)+len(thin)*12+15000)
  assert sourceowned<=150000 and time.monotonic()-start<50
  ev={'SOURCE_ONLY':True,'candidate':{'row':17,'column':79,'current_full_shape':[166,13],'prospective_pivot_NOT_proved_nonzero':True,'selection_proxy':'remaining recipe col79minimum; remainingoldNi56termcount row17minimum1319; NOTnewrank/minimumactualsupport'},'actual_minor_identity':'C17,79=(P7*F9)^2*B; B=FULLF11*F17,79-F17,53*F6,79','FULL_F11_preserved_1203DAG_NOT_face411':True,'new_functional_request':{'kind':'ONE coefficient functional of B, no VALUE executed','exponent5':target,'definition':'knownAactualwitness(-3,0,6,-1,10)+positiveSourceLexCorner(F17,79)','all5coords_retained':True,'source_positive_minor_corner5':minorcorner,'target_minus_positive_corner5':delta,'zero_vertex_OR_face_pruning_FORBIDDEN':True},'positive_source_profiles':{k:serial(v)for k,v in variants.items()},'old700_callback_closure_source_nodes':len(needed),'selected48_literal_payloads':descriptor,'old_functions_unexpanded_not_cached_value_guess':True,'future_projection_support_lifetime_VALUE_budget_NOT_completed':True,'no_face_or_matrix_coefficient_VALUE_GO':True,'source_whole_numeric_upper':sourceowned}
  put('chosen-functional-and-thin-source-inputs.json',ev)
  for f in[Path(__file__),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/'code/s61_sparse_generic5_source_view_profile_20261105_v64.py',ROOT/'code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py',ROOT/'research/notes'/f'{ID}.md']:pins[str(f.relative_to(ROOT))]=sha(f)
  rawbytes=sum((ROOT/p).stat().st_size for p in pins)+4*1024*1024+262144;assert rawbytes<32*1024*1024
  report={'id':ID,'all_passed':True,'SOURCE_PREFLIGHT_only_no_signed_coefficients':True,'chosen_candidate':[17,79],'fixed_functional_exponent5':target,'source_whole_numeric_upper':sourceowned,'fresh_RAW_upper':rawbytes,'old_callback_closure_nodes':len(needed),'future_projection_SOURCE_and_VALUE_gate_pending':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report)
  print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  if not(out/'failure.json').exists():put('failure.json',{'all_passed':False,'SOURCE_only':True,'error':repr(exc),'no_signed_values':True,'no_automatic_retry':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');main(p.parse_args())
