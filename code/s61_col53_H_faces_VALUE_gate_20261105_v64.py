"""SOURCE/code gate only for the two complete fixed H faces. No signed face VALUE."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,ast,hashlib,json,resource,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
ID='S61_COL53_H_FACES_VALUE_GATE_20261105_v64'
SOURCE='S61_COL53_H_FACES_SOURCE_20261105_v64'
WORKER='code/s61_col53_H_faces_values_20261105_v64.py'
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.prepare
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={}
 def pin(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return p
 def read(p):return json.loads(pin(p).read_text())
 def put(n,x):
  p=out/n;assert not p.exists();p.write_bytes(enc(x));return p
 try:
  sb=ROOT/'reports'/SOURCE;S=read(sb/'report.json')
  assert inputs[str((sb/'report.json').relative_to(ROOT))]=='5e66838fb845def468a9a1c883fd08a5836f21fdfd3a88de8f4d62109013d718'and S['all_passed']
  docs=[]
  for name in['compact-face-graph.json','positive-face-profiles-and-VALUE-lifetimes.json','two-complete-positive-H-faces.json']:
   path=sb/name;assert sha(path)==S['evidence_sha256'][name];docs.append(read(path))
  D,P,F=docs;ns=D['nodes'];N=len(ns);root=D['root'];assert N==1203 and root==N-1
  assert [x['terms']for x in F['faces']]==[766,214]and F['full_snapshot_equality']
  meta_parts={'runtime_and_python_code_reserve':25000,'loaded_compact_graph':count_numeric(D),'loaded_SOURCE_profiles':count_numeric(P),'loaded_both_root_faces':count_numeric(F),'indices_Mlimits_required_active_counts_uses_cache_bookkeeping':12*N,'root_allowed_5tuple_and_index_duplicate':6*766,'ROOT_document_numeric_upper':2000,'gate_budget_new_receipts_numeric_reserve':4000,'failure_numeric_reserve':8192}
  META=sum(meta_parts.values());assert META==98171
  box=[];max_width=0;max_abs_shift=0;maxq=0;maxbits=0;tempbits=0;integer_identities=0;operation_bounds=0
  for k,n in enumerate(ns):
   args=n['args'];assert all(0<=j<k for j in args)
   q=Q(n['source_profile']['known_QQ_content']);M=int(n['source_profile']['primitive_integer_abs_upper'])
   maxq=max(maxq,abs(q.numerator).bit_length(),q.denominator.bit_length());maxbits=max(maxbits,M.bit_length())
   assert n['actual_zero_by_source_q']==(not q)
   op=n['op'];children=[box[j]for j in args if box[j]is not None]
   if op=='input':b=n['literal_profile']['exponent_box5']if n['literal_profile']['terms_upper']else None
   elif not children:b=None
   elif op in['copy','immutable_SOURCE_view']:
    assert len(args)==1;b=[[box[args[0]][i][j]+n['shift'][i]for j in range(2)]for i in range(5)]
   elif op in['multiply','unmaterialized_multiply_for_SOURCE_union']:
    b=None if len(children)!=2 else[[sum(box[t][i][j]for t in args)+n['shift'][i]for j in range(2)]for i in range(5)]
   else:
    assert op in['add','SOURCE_fused_union_product']
    b=[[min(c[i][0]for c in children),max(c[i][1]for c in children)]for i in range(5)]
   box.append(b)
   if b:
    assert all(-64<=lo<=hi<=191 for lo,hi in b)
    max_width=max(max_width,max(hi-lo for lo,hi in b));max_abs_shift=max(max_abs_shift,max(abs(x)for v in b for x in v))
   if not q:continue
   if op=='input':
    if 'leaf_integer_ratio'in n:assert int(n['leaf_integer_ratio'])==1
    continue
   w=list(map(int,n['integer_weights']));qa=[Q(ns[j]['source_profile']['known_QQ_content'])for j in args];Ma=[int(ns[j]['source_profile']['primitive_integer_abs_upper'])for j in args]
   if op in['copy','immutable_SOURCE_view']:
    assert qa[0]*Q(n['scalar'])==q*w[0]and abs(w[0])*Ma[0]<=M
    temp=Ma[0]*abs(w[0])
   elif op in['multiply','unmaterialized_multiply_for_SOURCE_union']:
    assert qa[0]*qa[1]*Q(n['scalar'])==q*w[0]
    temp=Ma[0]*Ma[1]*abs(w[0])
   else:
    assert len(w)==len(args)
    for j in range(len(args)):assert qa[j]*(Q(n.get('scalarB','1'))if j else Q(1))==q*w[j]
    assert sum(abs(a)*b for a,b in zip(w,Ma))<=M
    temp=M+max([abs(a)*b for a,b in zip(w,Ma)]+[0])
   tempbits=max(tempbits,temp.bit_length());integer_identities+=len(w)
  assert maxbits<=640 and tempbits<=1921 and maxq<=92 and max_width<256
  forecasts=[];peak=META;maxterms=0;maxmaps=0;maxnode=0;all_updates=0;leafpaths=set();node_checks=0
  for mode in range(2):
   plan=D['mode_plans'][mode];required=set(plan['VALUE_needed']);active={int(k):v for k,v in plan['active_tied_args'].items()}
   counts={int(k):v for k,v in P['profiles'][mode]['face_terms_by_node'].items()};uses={k:0 for k in required};live={};rows=[];modepeak=META
   for k in required:
    if not ns[k]['actual_zero_by_source_q']:
     for c in active[k]:assert c in required;uses[c]+=1
   uses[root]+=1
   for k in sorted(required):
    n=ns[k];h=n['H_extrema'][mode];args=active[k]
    if n['actual_zero_by_source_q']:assert not Q(n['source_profile']['known_QQ_content'])
    elif n['op']=='input':assert args==[]
    elif n['op']in['copy','immutable_SOURCE_view']:assert len(args)==1 and h==ns[args[0]]['H_extrema'][mode]+n['shift'][0]
    elif n['op']in['multiply','unmaterialized_multiply_for_SOURCE_union']:
     assert len(args)==2 and h==sum(ns[t]['H_extrema'][mode]for t in args)+n['shift'][0]
     M=int(n['source_profile']['primitive_integer_abs_upper']);a,b=args
     assert min(counts[a],counts[b])*int(ns[a]['source_profile']['primitive_integer_abs_upper'])*int(ns[b]['source_profile']['primitive_integer_abs_upper'])*abs(int(n['integer_weights'][0]))<=M
     operation_bounds+=1;all_updates+=counts[a]*counts[b]
    else:assert args==[t for t in n['args']if ns[t]['H_extrema'][mode]==h]
    size=0 if n['actual_zero_by_source_q']else counts[k];old=sum(live.values())
    decoded=16152 if n['op']=='input'and size else 0
    owned=META+7*(old+size)+decoded+2048
    modepeak=max(modepeak,owned);peak=max(peak,owned);maxterms=max(maxterms,old+size);maxmaps=max(maxmaps,len(live)+1);maxnode=max(maxnode,size)
    rows.append({'node':k,'existing_live_terms':old,'new_output_before_release':size,'leaf_decode_numeric_upper':decoded,'whole_VALUE_numeric_upper':owned})
    live[k]=size;node_checks+=1
    if n['op']=='input'and not n['actual_zero_by_source_q']:
     source=n['source']
     if 'binding'in source:leafpaths.add(OLD+'/'+source['binding']['scope'])
     elif'file'in source:leafpaths.add(source['file'])
     elif source.get('QQ_monomial'):pass
     else:assert source.get('known_source_zero')
    if not n['actual_zero_by_source_q']:
     for c in args:
      uses[c]-=1
      if not uses[c]:live.pop(c,None)
    if not uses[k]:live.pop(k,None)
   assert live=={root:counts[root]}
   codec=META+10*counts[root]+4096;peak=max(peak,codec)
   forecasts.append({'mode':mode,'H':[-3,21][mode],'schedule':rows,'max_normal_owned':modepeak,'codec_owned':codec})
  assert peak==126262 and maxterms==2912 and maxmaps==26 and maxnode==766 and peak<150000
  cell=lambda n:512+128*n
  facecaps=[cell(766),cell(214)];onecell=cell(maxnode)
  failuremaps=maxmaps+1
  failraw=128*maxterms+512*failuremaps+8192*failuremaps
  generated=sum(facecaps)+failraw+262144+65536;generatedcap=2*1024*1024
  assert generated<generatedcap and 3*maxnode+4096+512<=8192
  codecbits=maxq+2*640;assert codecbits<=1921<2048
  for path in leafpaths:pin(ROOT/path)
  for path in[WORKER,'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py','code/s61_sparse_generic5_source_view_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py',str(Path(__file__).relative_to(ROOT))]:
   pin(ROOT/path);ast.parse((ROOT/path).read_text())
  for name in['S61_COL53_H_FACES_SOURCE_20261105_v64.md','S61_COL53_H_FACES_VALUE_GATE_20261105_v64.md']:pin(ROOT/'research/notes'/name)
  workerhash=sha(ROOT/WORKER)
  helpresult=subprocess.run([sys.executable,str(ROOT/WORKER),'--help'],capture_output=True,text=True,timeout=10)
  assert helpresult.returncode==0 and not(ROOT/'reports/S61_COL53_H_FACES_VALUES_20261105_v64').exists()
  put('worker-help-receipt.json',{'worker_sha256':workerhash,'argv':['--help'],'exit_code':helpresult.returncode,'stdout':helpresult.stdout,'no_VALUE_execution':True})
  import tempfile
  with tempfile.TemporaryDirectory(dir=ROOT/'work')as td:
   for i,A in enumerate([{}, {(-3,1,-7,2,10):-6,(-3,0,4,-1,7):12}]):
    path=Path(td)/('fixture'+str(i)+'.bin');q=Q(-7,13);rec=write_cell(path,q,A);assert validate_cell(path,q,A,rec)
  inputbytes=sum((ROOT/p).stat().st_size for p in inputs)
  fresh=inputbytes+512*1024+65536+generatedcap+262144;assert fresh<32*1024*1024
  b={'worker_sha256':workerhash,'VALUE_metadata_parts':meta_parts,'VALUE_metadata_base':META,'whole_VALUE_numeric_upper':peak,'VALUE_global_cap':150000,'leaf_decode_numeric_upper':16152,'failure_numeric_reserve':8192,'one_cell_RAW_upper':onecell,'face_binary_RAW_upper':facecaps,'generated_RAW_upper':generatedcap,'fresh_RAW_upper':fresh,'source_input_RAW_bytes':inputbytes,'future_gate_all_evidence_RAW_reserve':512*1024,'future_ROOT_report_RAW_cap':65536,'future_ROOT_report_numeric_cap':2000,'normal_commit_bits':640,'normal_temp_bits':1921,'stored_VALUE_bits_unchanged':2048,'actual_SOURCE_q_max_bits':maxq,'committed_SOURCE_integer_bits':maxbits,'ordinary_SOURCE_temp_bits':tempbits,'codec_naive_numerator_temp_bits':codecbits,'max_complete_cache_plus_work_terms':maxterms,'max_complete_cache_plus_work_maps':maxmaps,'max_single_face_map_terms':maxnode,'failure_cache_RAW_upper':failraw,'all5coords_VALUE7_per_term':True,'time':{'math_s':50,'signal_s':55,'hard_s':60,'AS_bytes':2000000000},'no_signed_face_VALUE_run':True}
  put('whole-integer-face-VALUE-budget.json',b)
  put('actual-code-lifetime-and-codec-source-proof.json',{'metadata':meta_parts,'reconstructed_forecasts':forecasts,'node_tie_and_lastuse_checks':node_checks,'integer_scale_identities':integer_identities,'positive_coefficient_majorant_product_checks':operation_bounds,'VALUE_pair_count_majorant':all_updates,'failure_codec_reserve_proof':{'extra_3n_plus4096_plus512':3*maxnode+4096+512,'reserve':8192,'failure_index_streamed_not_array':True,'existing_maps_never_copied':True},'codec_integer_bounds':{'sourceq_max_bits':maxq,'committed640':640,'naive_q_c_G_numerator':codecbits,'stored2048_unchanged':True},'source_box_max_width':max_width,'source_box_max_abs_coordinate':max_abs_shift,'exact_binary_RAW_per_term_bound':128,'exact_binary_fixed_header_bound':512,'all_actual_leaf_files':sorted(leafpaths),'future_VALUE_only_full_H_minus3_plus21':True,'no_reuse_of_ten_vertex_zero_for_pruning':True})
  report={'id':ID,'all_passed':True,'SOURCE_CODE_ONLY_no_signed_face_VALUE':True,'fixed_source_report_sha256':sha(sb/'report.json'),'worker_sha256':workerhash,'future_whole_VALUE_numeric_upper':peak,'future_fresh_RAW_upper':fresh,'time_60s_2GB':True,'each_complete_map_all5coords_charge7':True,'zero_q_only_pruning_not_vertexzeros':True,'failure_live_cache_preservation':True,'scope':'ONE prepared two complete H-face job only; fixed row6/col53; full167x14/87+66/P290once/F9F10/complements unchanged; no entry zero/nonzero/rank claim','input_sha256':inputs,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start}
  put('report.json',report);assert sum(p.stat().st_size for p in out.iterdir())<=512*1024
  print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  if not(out/'failure.json').exists():put('failure.json',{'all_passed':False,'SOURCE_CODE_ONLY':True,'error':repr(exc),'no_face_coefficients_evaluated':True,'no_automatic_retry':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');main(p.parse_args())
