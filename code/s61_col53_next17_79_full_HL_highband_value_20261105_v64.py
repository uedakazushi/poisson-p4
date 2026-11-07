"""Prepared incremental high-band VALUE; no execution without peer+ROOT pins."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,time
from s61_col53_projection_leaf_io_20261105_v64 import read_leaf
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
from s61_col53_next17_79_streamed_value_20261105_v64 import recover_committed
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_FULL_HL_HIGHBAND_VALUE_20261105_v64'
GATE='S61_COL53_NEXT17_79_FULL_HL_HIGHBAND_VALUE_GATE_20261105_v64'
PROPOSAL='reports/S61_COL53_NEXT17_79_FULL_HL_EXTENSION_BUDGET_SOURCE_20261105_v64'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64'
OLD='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64'
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def fixed_q_terms(body,qi,n,M):
 """Restore fixed-q integers, independently for a low or a high codec band."""
 assert body[:8]==b'Q5CELL1\n'
 if not n:
  assert body==b'Q5CELL1\n\0';return
 assert qi and body[8:9]==b'\1';p=decode(body[9:]);assert len(p['terms'])==n
 ratio=Q(p['content'])/qi;assert ratio.denominator==1 and ratio>0;G=ratio.numerator
 assert(Q(p['content'])>0)==(qi>0)and G.bit_length()<=824
 for t in p['terms']:
  e=tuple(p['shift'][j]+int(t[j+1])for j in range(5));v=int(t[0])*G
  assert v and abs(v).bit_length()<=1921 and abs(v)<=M;yield e,v
def frame_record(k,off,rec,actual,H,oo,ob,A,oh,qi,cap,lead,mathG):
 return{'node':k,'stream_offset':off,'bytes':rec['bytes'],'payload_sha256':rec['sha256'],'new_high_actual_slots':actual,'source_high_positive_bound':H,'old_low_reference':{'file':OLD+'/node-values.q5stream','offset':oo,'bytes':ob,'actual_slots':A,'sha256':oh,'unchanged_byteidentity_checked':True},'fixed_q':str(qi),'old_cap':cap,'fixed_SOURCE_lead_UNCHANGED':lead,'math_G_high':mathG,'codec':rec,'two_NEW_generations_byteidentical':True,'full_function_is_disjoint_OLD_LOW_plus_NEW_HIGH':True}
def prepare():
 out=ROOT/'reports'/GATE;out.mkdir(exist_ok=False)
 rp=ROOT/PROPOSAL/'report.json';assert sha(rp)=='f9a3bc05e0975ef20698dd6f4fde3f3c19f3c021c847edf4fcb3ec65a95f2635';r=json.loads(rp.read_text());assert r['all_passed']
 pp=rp.parent/'source-cauchy-two-band-lifetime-and-raw-proposal.json';assert sha(pp)==r['evidence_sha256'][pp.name];P=json.loads(pp.read_text())
 pins={**P['future_fixed_READ_sha256'],str(rp.relative_to(ROOT)):sha(rp),str(pp.relative_to(ROOT)):sha(pp)}
 for f,h in pins.items():assert sha(ROOT/f)==h
 gp=ROOT/FINAL/'compiled-full-SLP-and-fixed-projection-plan.json';D=json.loads(gp.read_text());assert count_numeric(D)==65701
 retired=count_numeric(D.pop('active_tied_HL_children'));assert retired==3078 and count_numeric(D)==62623
 N=1884;census={'runtime':25000,'resident_graph_after_unused_active_retirement':62623,'remaining_numeric_mode_keys':N,'order':N,'combined_index_KEY_old3_high3_TWO_hash_numerics':9*N,'current_gate_plan_codec_root_control':4096}
 assert sum(census.values())==112443 and 112443+4000+1557==118000
 plan=out/'incremental-node-plan.jsonl'
 maxline=0
 with plan.open('wb')as f:
  for row in P['all_node_records']:
   f.write(enc([row['node'],row['new_positive_slots'],row['full_actual_slots_upper'],row['old_actual_slots'],row['full_positive_slots'],row['new_codec_RAW_upper']]))
   k=row['node'];n=D['nodes'][k];q=Q(n['q']);mb=int(n['M']).bit_length()
   # 30103/100000 > log10(2); fake all9 strings are RAW length bounds,
   # not coefficients or a parameter/source point.
   digits=lambda b:max(1,(b*30103+99999)//100000)
   G=10**digits(mb)-1;content='-'+'9'*digits(abs(q.numerator).bit_length()+mb)+'/'+'9'*digits(q.denominator.bit_length())
   rr={'terms':row['new_positive_slots'],'G':G,'content':content,'shift':[-191]*5,'value_sha256':'f'*64,'sha256':'f'*64,'bytes':row['new_codec_RAW_upper']}
   worst=frame_record(k,10676690,rr,row['new_positive_slots'],row['new_positive_slots'],859628,859628,row['old_actual_slots'],'f'*64,q,D['needed_node_degree'][str(k)],n['lead5'],G)
   maxline=max(maxline,len(enc(worst)))
 assert maxline==1749 and maxline<=2048
 assert len(P['all_node_records'])==1884
 fixture='reports/S61_COL53_NEXT17_79_HIGHBAND_CODEC_FIXTURES_20261105_v64/report.json';fr=json.loads((ROOT/fixture).read_text());assert fr['all_passed']and fr['synthetic_only_no_chosen_values'];pins[fixture]=sha(ROOT/fixture)
 for f,h in fr['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
 for f in['code/s61_col53_next17_79_full_HL_highband_value_20261105_v64.py','code/s61_col53_projection_leaf_io_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py','code/s61_col53_next17_79_streamed_value_20261105_v64.py','code/s61_sparse_generic5_source_view_profile_20261105_v64.py','code/s61_sparse_generic5_source_packed5_profile_v2_20261105_v64.py']:pins[f]=sha(ROOT/f)
 anchor='reports/S61_COL53_H_FACES_VALUES_20261105_v64/H-minus3-face.bin';assert sha(ROOT/anchor)=='cc9165dc8ec941004451f7924944e839c36db3e0b4e939be17b4a8f9afd88bce';pins[anchor]=sha(ROOT/anchor)
 fixed=sum((ROOT/f).stat().st_size for f in pins)+plan.stat().st_size
 # Includes current gate/budget, two future<=64KiB reviews and stop controls.
 raw=fixed+P['ALL_generated_RAW_upper']+524288+131072+262144;assert raw<32*1024*1024
 B={'worker_sha256':sha(__file__),'metadata_base':118000,'explicit_census':census,'exact_census_total':112443,'two_review_docs_numeric_reserve':4000,'unallocated_numeric_margin':1557,'normal_upper':143570,'prefix_upper':138668,'same_global_guard':150000,'source_M_bits':824,'temporary_bits':1921,'hard_VALUE_bits':2048,'manifest_line_proved_RAW_upper':maxline,'manifest_line_fixed_RAW_cap':2048,'generated_RAW_upper':P['ALL_generated_RAW_upper'],'fresh_RAW_upper':raw,'math_s':50,'signal_s':55,'hard_s':60,'AS_bytes':2000000000,'two_NEW_highband_generations':True,'old_low_generation_never_modified':True,'unused_active_retirement_metadata_only':True}
 (out/'highband-whole-lifetime-and-storage-budget.json').write_bytes(enc(B))
 d={'id':GATE,'all_passed':True,'SOURCE_CODE_PREPARATION_ONLY_highband_VALUE_UNRUN':True,'worker_sha256':B['worker_sha256'],'normal_upper':143570,'prefix_upper':138668,'fresh_RAW_upper':raw,'all5_VALUE7_fixedq_oldlow_preserved':True,'new_peer_and_ROOT_GO_required':True,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()}}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'));print('RAW',raw)
def value(a):
 assert all([a.gate_file,a.gate_sha,a.source_peer_file,a.source_peer_sha,a.root_review_file,a.root_review_sha])
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);tmp=ROOT/'work'/(ID+'-canonical-high-cells');tmp.mkdir(exist_ok=False)
 pins={};index={};held=[];work={};memo={};transient=0;peak=0;context={};oldbyteschecked=0;stream=None;manifest=None;oldstream=None;oldmanifest=None;planstream=None;rootrecords={};B={'generated_RAW_upper':1048576};D=None
 def put(name,d):
  body=enc(d);total=sum(p.stat().st_size for p in out.iterdir()if p.is_file())+sum(p.stat().st_size for p in tmp.iterdir()if p.is_file());assert total+len(body)<=B['generated_RAW_upper'];(out/name).write_bytes(body)
 def alarm(s,f):raise TimeoutError('incremental-highVALUE55s; preserve prefix; no retry')
 signal.signal(signal.SIGALRM,alarm);signal.alarm(55)
 def own(extra=4096,new=False):
  nonlocal peak
  seen={id(v):len(v)for v in[work,*held]};n=118000+7*(sum(seen.values())+int(new))+transient+extra;peak=max(peak,n);assert n<=150000,('prospective_HIGH_VALUE_owner',n,context)
 def raw(v):assert abs(v).bit_length()<=1921;return v
 def setterm(e,v,M):
  assert len(e)==5 and all(isinstance(x,int)for x in e)and abs(v)<=M and abs(v).bit_length()<=824
  kk=context['node'];ll=D['nodes'][kk]['lead5'];assert e[:2]==tuple(ll[:2])and e[2]-ll[2]>D['needed_node_degree'][str(kk)]
  if v:
   if e not in work:own(new=True)
   work[e]=v
  else:work.pop(e,None)
 def child(k):
  nonlocal transient
  # 9N numeric inventory: key; each band's offset/bytes/count; TWO hashes.
  oo,ob,on,ho,hb,hn,oh,hh=index[k];n=D['nodes'][k];qi=Q(n['q']);M=int(n['M']);cap=D['needed_node_degree'][str(k)];lead=n['lead5'];v={};held.append(v)
  for high,nb,nc,want in[(False,ob,on,oh),(True,hb,hn,hh)]:
   transient=8*nc+2048;own()
   if high:body=(tmp/('node-'+str(k)+'.bin')).read_bytes()
   else:oldstream.seek(oo);body=oldstream.read(nb)
   assert len(body)==nb and hashlib.sha256(body).hexdigest()==want
   for e,c in fixed_q_terms(body,qi,nc,M):
    assert e[:2]==tuple(lead[:2])and((e[2]-lead[2]>cap)if high else(0<=e[2]-lead[2]<=cap))and e not in v
    own(new=True);v[e]=c
   del body;transient=0;own()
  assert len(v)==on+hn;return v
 def anchor_high(qi,lead,cap):
  nonlocal transient
  transient=8*411+2048;own();body=(ROOT/'reports/S61_COL53_H_FACES_VALUES_20261105_v64/H-minus3-face.bin').read_bytes();p=decode(body[9:]);ratio=Q(p['content'])/qi;assert ratio.denominator==1 and ratio>0;matched=0
  for t in p['terms']:
   e=tuple(p['shift'][j]+int(t[j+1])for j in range(5))
   if e[:2]==tuple(lead[:2])and e[2]-lead[2]>cap:
    c=raw(int(t[0])*ratio.numerator);assert work.get(e)==c;matched+=1
  assert matched==len(work);del p,body;transient=0;own()
 try:
  gp=ROOT/a.gate_file;assert sha(gp)==a.gate_sha;g=json.loads(gp.read_text());assert g['all_passed'];pins[a.gate_file]=a.gate_sha
  bp=gp.parent/'highband-whole-lifetime-and-storage-budget.json';assert sha(bp)==g['evidence_sha256'][bp.name];B=json.loads(bp.read_text());assert sha(__file__)==B['worker_sha256'];pins[str(bp.relative_to(ROOT))]=sha(bp)
  for filename,want in[(a.source_peer_file,a.source_peer_sha),(a.root_review_file,a.root_review_sha)]:
   p=ROOT/filename;assert sha(p)==want and p.stat().st_size<=65536;r=json.loads(p.read_text());assert r['all_passed']and count_numeric(r)<=2000
   assert r['input_sha256'][a.gate_file]==a.gate_sha and r['input_sha256']['code/s61_col53_next17_79_full_HL_highband_value_20261105_v64.py']==B['worker_sha256'];pins[filename]=want
  for f,h in g['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
  plan=gp.parent/'incremental-node-plan.jsonl';assert sha(plan)==g['evidence_sha256'][plan.name];pins[str(plan.relative_to(ROOT))]=sha(plan)
  D=json.loads((ROOT/FINAL/'compiled-full-SLP-and-fixed-projection-plan.json').read_text());assert count_numeric(D)==65701
  retired=count_numeric(D.pop('active_tied_HL_children'));assert retired==3078 and count_numeric(D)==62623
  ns=D['nodes'];caps=D['needed_node_degree'];order=sorted(map(int,caps));assert len(order)==1884;root=D['minor_B_root']
  oldfile=ROOT/OLD/'node-values.q5stream';oldindex=ROOT/OLD/'node-values-index.jsonl';oldsha=sha(oldfile);oldindexsha=sha(oldindex)
  oldstream=oldfile.open('rb');oldmanifest=oldindex.open('rb');planstream=plan.open('rb');stream=(out/'node-values.q5stream').open('w+b');manifest=(out/'node-values-index.jsonl').open('wb');tmpbytes=0
  for pos,k in enumerate(order):
   n=ns[k];qi=Q(n['q']);M=int(n['M']);lead=n['lead5'];cap=caps[str(k)];memo={};held=[];work={};context={'node':k,'order_position':pos,'phase':'fixed_old_frame_and_current_plan'};own()
   pp=json.loads(planstream.readline());assert len(pp)==6 and pp[0]==k;_,H,fullupper,A,F,rawupper=pp
   low=json.loads(oldmanifest.readline());assert low['node']==k and low['actual_nonzero_slots']==A and Q(low['fixed_SOURCE_q'])==qi and low['fixed_SOURCE_lead_UNCHANGED']==lead
   oo=low['stream_offset'];ob=low['bytes'];oh=low['payload_sha256'];oldstream.seek(oo);body=oldstream.read(ob);assert len(body)==ob and hashlib.sha256(body).hexdigest()==oh;oldbyteschecked+=ob;del body
   if H and qi:
    if n['op']=='input':
     transient=16152;own();p,shift,dec=read_leaf(n);assert dec+2048<=16152
     for ti,t in enumerate(p['terms']):
      context['input_cursor']=ti;e=tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5))
      if e[:2]==tuple(lead[:2])and e[2]-lead[2]>cap:setterm(e,int(t[0]),M)
     del p;transient=0
    else:
     args=n['args'];ch=args if n['op']in['copy','multiply']else[c for c in args if ns[c]['q']!='0'and ns[c]['lead5'][:2]==lead[:2]]
     for c in dict.fromkeys(ch):memo[c]=child(c)
     if n['op']=='copy':
      c=args[0];w=Q(ns[c]['q'])*Q(n['scalar'])/qi;assert w.denominator==1
      for ti,(e,v)in enumerate(memo[c].items()):
       context['copy_cursor']=ti;ee=tuple(e[j]+n['shift5'][j]for j in range(5))
       if ee[2]-lead[2]>cap:setterm(ee,raw(v*w.numerator),M)
     elif n['op']=='multiply':
      aa,bb=args;w=Q(ns[aa]['q'])*Q(ns[bb]['q'])*Q(n['scalar'])/qi;assert w.denominator==1
      for ai,(e,v)in enumerate(memo[aa].items()):
       for bi,(f,z)in enumerate(memo[bb].items()):
        context['product_cursor']=[ai,bi];ee=tuple(e[j]+f[j]+n['shift5'][j]for j in range(5))
        if ee[2]-lead[2]>cap:
         term=raw(raw(v*z)*w.numerator);setterm(ee,raw(work.get(ee,0)+term),M)
     else:
      assert n['op']=='add'
      for j,c in enumerate(args):
       if c not in memo:continue
       w=Q(ns[c]['q'])*(Q(n['scalar'])if j else Q(1))/qi;assert w.denominator==1
       for ti,(e,v)in enumerate(memo[c].items()):
        context['sum_cursor']=[j,ti]
        if e[2]-lead[2]>cap:setterm(e,raw(work.get(e,0)+raw(v*w.numerator)),M)
   assert len(work)<=H and A+len(work)<=fullupper
   # Exact last use: child maps are released BEFORE anchor and codec.
   memo.clear();held.clear();own()
   if k==D['FULL_F11_original_root']:anchor_high(qi,lead,cap)
   context['phase']='new_highband_codec_commit';own(extra=3*len(work)+4096)
   assert tmpbytes+rawupper+stream.seek(0,2)+rawupper+manifest.tell()+2048+1572864<=B['generated_RAW_upper']
   path=tmp/('node-'+str(k)+'.bin');rec=write_cell(path,qi,work);assert rec['bytes']<=rawupper and validate_cell(path,qi,work,rec)
   mathG=0 if not work else rec['G']
   if not work:assert rec['G']==1 and rec['content']=='0'
   else:assert Q(rec['content'])/qi==mathG>0
   tmpbytes+=rec['bytes'];off=stream.seek(0,2)
   with path.open('rb')as t:
    for body in iter(lambda:t.read(65536),b''):stream.write(body)
   stream.flush();end=stream.tell();assert end-off==rec['bytes'];stream.seek(off);body=stream.read(rec['bytes']);assert hashlib.sha256(body).hexdigest()==rec['sha256'];stream.seek(end);del body
   record=frame_record(k,off,rec,len(work),H,oo,ob,A,oh,qi,cap,lead,mathG)
   line=enc(record);assert len(line)<=2048;manifest.write(line);manifest.flush();index[k]=(oo,ob,A,off,rec['bytes'],len(work),oh,rec['sha256'])
   if k in D['four_callback_roots'].values()or k==root:rootrecords[str(k)]=record
   if k==root:
    assert A==0;witness=None
    if work:
     e=min(work);witness={'exponent5':list(e),'coefficient_QQ':str(qi*work[e])}
    put('actual-complete-leading-HL-face-summary.json',{'HL':[-5,0],'complete_positive_S_range':[7,38],'old_low_slots':A,'new_high_actual_slots':len(work),'whole_leading_HL_face_ZERO':not work,'nonzero_witness':witness,'whole_minor_ZERO_NOT_inferred':True,'full_F11_retained':True,'no_new_matrix_RHS_update':True})
   work={};assert time.monotonic()-start<50
  assert not oldmanifest.read(1)and not planstream.read(1);stream.close();manifest.close();stream=None;manifest=None;oldstream.close();oldmanifest.close();planstream.close();oldstream=None;oldmanifest=None;planstream=None
  assert sha(oldfile)==oldsha and sha(oldindex)==oldindexsha and oldbyteschecked==oldfile.stat().st_size
  cursor=recover_committed(out);assert cursor['nodes_committed']==1884 and cursor['last_node']==root
  put('old-prefix-and-new-two-generation-custody.json',{'cursor':cursor,'all1884_old_low_payloads_SHA_equal_original':True,'old_bytes_verified':oldbyteschecked,'old_global_stream_SHA_UNCHANGED':oldsha,'old_global_index_SHA_UNCHANGED':oldindexsha,'all1884_new_high_segments_and_tempfiles_preserved':True,'new_temp_directory':str(tmp.relative_to(ROOT)),'root_records':rootrecords,'same_FIXED_SOURCE_q_and_leads':True,'paired_index_numerics_per_node_including_both_hashes':9})
  generated=sum(p.stat().st_size for p in out.iterdir())+sum(p.stat().st_size for p in tmp.iterdir());assert generated<=B['generated_RAW_upper'];fresh=sum((ROOT/f).stat().st_size for f in pins)+generated+262144;assert fresh<=B['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'complete_leading_HL_face_only':True,'nodes_committed':1884,'whole_leading_HL_face_ZERO':not witness,'nonzero_witness':witness,'actual_owned_peak':peak,'same_global_guard':150000,'META':118000,'no_old_low_arithmetic_replay_or_changes':True,'only_HIGH_output_exponents_computed':True,'mixed_low_times_low_HIGH_contributions_retained':True,'new_matrix_RHS_update_NOT_performed':True,'fresh_RAW_upper_including_control_reserve':fresh,'input_sha256':pins,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};put('report.json',report);signal.alarm(0);print('report_sha256',sha(out/'report.json'));print('face_witness',witness)
 except Exception as exc:
  for f in[stream,manifest]:
   if f is not None:f.flush()
  recovered=recover_committed(out)if(out/'node-values.q5stream').exists()else{'nodes_committed':0}
  # All child functions are already immutable OLD/HIGH frame references;
  # they are not needed to preserve the current incomplete accumulator.
  memo.clear();held.clear()
  if D is not None and work:
   own(extra=3*len(work)+4096);k=context['node'];rec=write_cell(out/'incomplete-current-high-band.bin',Q(D['nodes'][k]['q']),work);validate_cell(out/'incomplete-current-high-band.bin',Q(D['nodes'][k]['q']),work,rec)
   put('incomplete-current-high-band-status.json',{'node':k,'context':context,'complete_band_NOT_claimed':True,'codec':rec})
  put('failure.json',{'all_passed':False,'error':repr(exc),'context':context,'committed_cursor':recovered,'all_old_low_immutable':True,'all_completed_new_segments_and_partial_preserved':True,'no_automatic_retry':True});raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare-gate',action='store_true');p.add_argument('--value',action='store_true');p.add_argument('--gate-file');p.add_argument('--gate-sha');p.add_argument('--source-peer-file');p.add_argument('--source-peer-sha');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');a=p.parse_args()
 if a.prepare_gate:prepare()
 else:assert a.value;value(a)
