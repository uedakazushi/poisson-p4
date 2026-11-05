#!/usr/bin/env python3
"""Independent stream code/resource proof and UNRELATED codec fixtures only."""
import ast,hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
ROOT=pathlib.Path(__file__).resolve().parents[1]
GATE='S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64';ID=GATE+'_REVIEW'
WORKER='code/s61_col53_next17_79_streamed_value_20261105_v64.py'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def uenc(n):
 z=[]
 while n>127:z.append((n%128)+128);n//=128
 return bytes(z+[n])
def ienc(n):
 a=abs(n);b=a.to_bytes((a.bit_length()+7)//8,'big');return uenc(2*len(b)+int(n<0))+b
def own_decode(data):
 assert data[:8]==b'Q5CELL1\n';at=9
 if data[8]==0:assert len(data)==9;return Q(0),[0]*5,[]
 assert data[8]==1
 def u():
  nonlocal at
  n=0;power=1;start=at
  while True:
   c=data[at];at+=1;n+=(c%128)*power
   if c<128:break
   power*=128
  assert data[start:at]==uenc(n);return n
 def i():
  nonlocal at
  tag=u();size=tag//2;b=data[at:at+size];at+=size;assert len(b)==size and(not b or b[0]);n=int.from_bytes(b,'big');assert n or not tag%2;return-n if tag%2 else n
 n,d=i(),i();q=Q(n,d);assert n and d>0 and(q.numerator,q.denominator)==(n,d)
 shift=[]
 for j in range(5):x=u();shift.append(x//2 if x%2==0 else-(x//2)-1)
 terms=[]
 for j in range(u()):
  c=i();assert c;terms.append((tuple(shift[k]+u()for k in range(5)),c))
 assert at==len(data)and terms
 assert [e for e,c in terms]==sorted({e for e,c in terms},key=lambda e:(sum(e),e),reverse=True)
 assert math.gcd(*(abs(c)for e,c in terms))==1
 return q,shift,terms
def own_encode(q,values):
 if not values:return b'Q5CELL1\n\0'
 g=math.gcd(*(abs(x)for x in values.values()));shift=[min(e[k]for e in values)for k in range(5)];c=q*g
 out=b'Q5CELL1\n\1'+ienc(c.numerator)+ienc(c.denominator)
 out+=b''.join(uenc(2*x if x>=0 else-2*x-1)for x in shift)+uenc(len(values))
 for e in sorted(values,key=lambda e:(sum(e),e),reverse=True):out+=ienc(values[e]//g)+b''.join(uenc(e[k]-shift[k])for k in range(5))
 return out
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  pp=ROOT/p;z=sha(pp)
  if h:assert z==h,(p,h,z)
  if p not in pins:raw+=pp.stat().st_size
  pins[p]=z;return json.loads(pp.read_text())if parse else None
 gp=f'reports/{GATE}/report.json';g=pin(gp,'c5823b8df2afda1ff7c4285889da3f45e06ace7d00d89fc27c72133464fa98e3');assert g['all_passed']
 for p,h in g['input_sha256'].items():pin(p,h,False)
 pin(WORKER,'a946d904a20a78bc7f9a211ee010620ac0280035bce6565eb7b729ecbaffe808',False)
 bp=f'reports/{GATE}/streamed-whole-lifetime-and-raw-budget.json';B=pin(bp,g['evidence_sha256'][pathlib.Path(bp).name]);assert B['worker_sha256']==pins[WORKER]
 custody=pin('reports/S61_COL53_NEXT17_79_STREAMED_VALUE_CUSTODY_20261105_v64/report.json','e75cf2b3576812d20022139e41f7096b701830e04703265f0255094dc9d6975f');assert custody['all_passed']
 for p,h in custody.get('input_sha256',{}).items():pin(p,h,False)
 # The lifetime companion is required independently of the unchanged gate.
 companion=pin('reports/S61_COL53_NEXT17_79_STREAMED_PHASE_COMPANION_20261105_v64/report.json','2577d4a15012a9ae99c806b93c55fe7f00c70d5c29ec25584eaa141d1aec94e9');assert companion['all_passed']
 for p,h in companion['input_sha256'].items():pin(p,h,False)
 assert(companion['normal_success_whole_upper'],companion['guard_prefix_whole_upper'],companion['corrected_child_decode_upper'],companion['anchor_prefix_upper_with_retained_children'])==(138535,138658,126982,128892)
 method=pin('reports/S61_COL53_NEXT17_79_STREAM_VALUE_METHOD_SOURCE_V2_20261105_v64/report.json','e15323af46f399a232e5ad482439b7c1e910e012ca6f85bcb6fa572cb2338c99')
 proofpath='reports/S61_COL53_NEXT17_79_STREAM_VALUE_METHOD_SOURCE_V2_20261105_v64/streamed-residency-height-codec-storage-source-proof.json'
 methodproof=pin(proofpath,method['evidence_sha256'][pathlib.Path(proofpath).name])
 final='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64/compiled-full-SLP-and-fixed-projection-plan.json';D=pin(final,g['input_sha256'][final])
 cp='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64/all-projected-SOURCE-node-counts.json';C=dict(pin(cp,g['input_sha256'][cp]));nodes=D['nodes'];order=sorted(C)
 for st,h in[('S61_COL53_NEXT17_79_FINALIZE_ROOT_SOURCE_REVIEW_20261105_v64','6174012f35862c6d12d5435c330eb3b529340f72431105b1d112def7b5200519'),('S61_COL53_NEXT17_79_PROJECTION_ROOT_SOURCE_REVIEW_20261105_v64','43b7ac645d4eb6a5184721b894d84ca9e5589e493c6add65d551cfe2927c6a32')]:assert pin(f'reports/{st}/report.json',h)['all_passed']
 root=pin('reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64/boundary-2997.json',g['input_sha256']['reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64/boundary-2997.json'])
 assert len(nodes)==2998 and len(order)==1884 and sum(C.values())==71584 and max(C.values())==129 and C[D['minor_B_root']]==129
 assert len(root)==len(set(map(tuple,root)))==129 and all(len(e)==5 for e in root)and D['target_full5']in root
 assert D['FULL_F11_original_root']==1202 and D['minor_B_root']==2997 and D['whole_F11_not_face411']
 review_prior=count([g,B,custody,companion,method,methodproof,D,C,root,pins])+25000;assert review_prior<150000 and raw+2_000_000<32*1024*1024
 save(out/'positive-before-source-code-synthetic-review.json',{'whole_review_numeric_upper':review_prior,'RAW_review_upper':raw+2000000,'no_chosen_VALUE_queries':True,'only_unrelated_synthetic_codec_tests':True})
 maxM=maxn=maxd=0;checked=[];phasepeaks={};phaseat={}
 def bound(name,k,n):
  if n>phasepeaks.get(name,0):phasepeaks[name]=n;phaseat[name]=k
 for k in order:
  n=nodes[k];args=n['args'];qi=Q(n['q']);M=int(n['M']);cap=D['needed_node_degree'][str(k)];lead=n['lead5'];assert len(lead)==5 and cap>=0
  maxM=max(maxM,M.bit_length());maxn=max(maxn,abs(qi.numerator).bit_length());maxd=max(maxd,qi.denominator.bit_length())
  if not qi or not C[k]:weights=[];partial=0
  elif n['op']=='input':weights=[];partial=M
  elif n['op']=='copy':
   a=args[0];w=Q(nodes[a]['q'])*Q(n['scalar'])/qi;assert w.denominator==1;weights=[w.numerator];partial=abs(w.numerator)*int(nodes[a]['M'])
  elif n['op']=='multiply':
   a,b=args;w=Q(nodes[a]['q'])*Q(nodes[b]['q'])*Q(n['scalar'])/qi;assert w.denominator==1;weights=[w.numerator];partial=min(C[a],C[b])*abs(w.numerator)*int(nodes[a]['M'])*int(nodes[b]['M'])
  else:
   assert n['op']=='add';weights=[]
   for j,a in enumerate(args):w=Q(nodes[a]['q'])*(Q(n['scalar'])if j else Q(1))/qi;assert w.denominator==1;weights.append(w.numerator)
   partial=sum(abs(w)*int(nodes[a]['M'])for w,a in zip(weights,args))
  assert partial<=M;checked.append({'node':k,'integer_weights':list(map(str,weights)),'positive_partial_bound_le_M':True})
  if n['op']=='add':
   children=D['active_tied_HL_children'][str(k)]
   for a in children:assert nodes[a]['lead5'][:2]==lead[:2]and D['needed_node_degree'][str(a)]>=cap+lead[2]-nodes[a]['lead5'][2]
  elif n['op']in['copy','multiply']:
   children=args;assert all(D['needed_node_degree'][str(a)]>=cap for a in children)
   expect=[sum(nodes[a]['lead5'][j]for a in args)+n['shift5'][j]for j in range(5)];assert expect==lead
  else:children=[]
  children=list(dict.fromkeys(children));assert all(a in C and a<k for a in children);ct=sum(C[a]for a in children);mc=max([C[a]for a in children]+[0]);ni=C[k]
  bound('normal_eval',k,118000+7*(ct+ni)+4096)
  bound('normal_child_decode_corrected',k,118000+7*ct+8*mc+2048+4096)
  bound('normal_codec',k,118000+7*(ct+ni)+3*ni+4096)
  if n['op']=='input'and qi:
   bound('normal_leaf_decode',k,118000+16152+7*ni+4096)
   bound('failure_leaf_then_codec',k,118000+16152+10*ni+4096)
  if k==D['FULL_F11_original_root']:
   bound('normal_anchor',k,118000+8*411+2048+7*(ct+ni)+4096)
   bound('failure_anchor_then_codec',k,118000+8*411+2048+7*ct+10*ni+4096)
 assert checked==methodproof['integer_q_normalization_and_partial_majorant_checks']and(maxM,maxn,maxd)==(824,7,164)
 assert phasepeaks['normal_leaf_decode']==138535 and phasepeaks['normal_child_decode_corrected']==126982 and phasepeaks['failure_leaf_then_codec']==138658
 assert max(phasepeaks.values())==138658<150000 and B['metadata_base']==118000 and B['whole_numeric_forecast']==138535
 # Count all numeric graph keys ourselves; the producer's65701 graph charge
 # is conservative compared with this independent51115 numeric-leaf reader.
 assert count(D)==51115<=65701 and count(C)==2*1884
 baseline=25000+65701+2*1884+1884+2*1884+4*1884
 # Both receipt documents can briefly coexist during RHS json decoding.
 control=2*2000+6*129+count([g,B])+512+512
 assert baseline==107657 and baseline+control<=118000
 assert methodproof['metadata_exact_total']==111753 and methodproof['future_fixed_metadata_base']==118000
 assert 28+2*824==1676<1921<2048 and 2*824==1648 and 28+824+164==1016
 one=512*1884+116*71584;generated=2*one+2048*1884+524288+1048576
 assert one==9268352 and generated==B['ALL_generated_RAW_upper']==23968000
 assert B['fresh_RAW_upper']==26155445 and B['math_s']==50 and B['signal_s']==55 and B['hard_s']==60 and B['AS_bytes']==2000000000
 # Inspect the actual operational code, not only the budget prose.
 text=(ROOT/WORKER).read_text();tree=ast.parse(text)
 for token in ['own(new=True)','transient=8*C[k]+2048','memo.clear();held.clear();work={}','ratio=Q(p[\'content\'])/qi','ratio.numerator>0']:assert token in text
 assert "raw(work.get(ee,0)+term)"in text and "own(extra=3*len(work)+4096)"in text
 assert "fixed_SOURCE_lead_UNCHANGED"in text and "projection_zero_only_NOT_whole_function_zero"in text
 assert "rec[2]==0"in text and "b'Q5CELL1\\n\\0'"in text
 assert "sum(p.stat().st_size for p in tmp.iterdir()"in text and "stream_prefix_bytes"in text
 # Execute ONLY recover_committed, extracted from the frozen AST, never the
 # selected signed VALUE function or its source polynomial recurrences.
 fnode=next(n for n in tree.body if isinstance(n,ast.FunctionDef)and n.name=='recover_committed');env={'Path':pathlib.Path,'json':json,'hashlib':hashlib};exec(compile(ast.Module(body=[fnode],type_ignores=[]),WORKER,'exec'),env);recover=env['recover_committed']
 fixture=out/'unrelated-codec-fixtures';fixture.mkdir();cases=[(Q(-7,30),{(-3,0,6,-1,10):12,(-3,0,7,-1,9):-18}),(Q(11,37),{(2,-2,0,4,1):-25}),(Q(-13,21),{}),(Q(0),{}),(Q(5,17),{(-64,0,2,-1,191):3,(191,-64,1,0,0):5}),(Q(-(1<<27),(1<<163)+1),{(-3,0,6,-1,10):3*(1<<400),(-3,0,7,-1,9):-5*(1<<400)})]
 receipts=[];bodies=[]
 for k,(qi,values)in enumerate(cases):
  p=fixture/f'case-{k}.bin';rec=write_cell(p,qi,values);assert validate_cell(p,qi,values,rec);data=p.read_bytes();assert data==own_encode(qi,values);stored,shift,terms=own_decode(data)
  if values:
   G=stored/qi;assert G.denominator==1 and G>0 and(stored>0)==(qi>0)
   assert{e:c*G.numerator for e,c in terms}==values and G==math.gcd(*(abs(x)for x in values.values()))
  else:assert stored==0 and rec['G']==1 and rec['content']=='0';G=Q(0)
  receipts.append({'fixture':k,'fixed_q':str(qi),'math_G':str(G),'codec_G':rec['G'],'bytes':len(data),'body_sha256':sha(p)});bodies.append(data)
 try:write_cell(fixture/'invalid-nonempty-q0.bin',Q(0),{(0,0,0,0,0):1});raise RuntimeError('q0 nonempty unexpectedly accepted')
 except AssertionError:pass
 streamdir=fixture/'prefix';streamdir.mkdir();parts=[bodies[0],bodies[2]];lines=[];offset=0
 for k,data in enumerate(parts):lines.append((json.dumps({'node':k,'stream_offset':offset,'bytes':len(data),'payload_sha256':hashlib.sha256(data).hexdigest()})+'\n').encode());offset+=len(data)
 (streamdir/'node-values.q5stream').write_bytes(b''.join(parts)+b'UNCOMMITTED');(streamdir/'node-values-index.jsonl').write_bytes(b''.join(lines)+b'{"node":2')
 cur=recover(streamdir);assert cur=={'nodes_committed':2,'last_node':1,'stream_prefix_bytes':offset,'manifest_prefix_bytes':sum(map(len,lines))}
 (streamdir/'node-values-index.jsonl').write_bytes(b''.join(lines)+(json.dumps({'node':2,'stream_offset':offset,'bytes':1,'payload_sha256':'0'*64})+'\n').encode());assert recover(streamdir)==cur
 assert b''.join((fixture/f'case-{k}.bin').read_bytes()for k in[0,2])==(streamdir/'node-values.q5stream').read_bytes()[:offset]
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 ev={'all_passed':True,'no_chosen_signed_VALUE_or_source_coeff_response_evaluated':True,'SOURCE1884_dependencies_and_all_integer_partial_majorants_rebuilt':True,'SOURCE129_root_full5keys_and_fixed_HL_sdegree2_projection_checked':True,
  'source_positive_vertices_not_changed_by_actual_cancellation':True,'fixed_q_recurrences_vs_primitive_G_restore_and_sign_proved':True,'empty_mathG0_codecG1_NO_qi_division':True,
  'independent_schema_encoder_decoder_and_synthetic_writer_fixtures':receipts,'prefix_manifest_exacthash_cursor_excludes_incomplete_or_invalid_tail':cur,
  'two_byteidentical_generations_and_all_RAW_versions_charged':True,'META_baseline_without_old4096_control':baseline,'fresh_control_bound_with_two_receipts_root129_and_misc':control,'fresh_META_sum':baseline+control,'META118000_suffices':True,
  'corrected_all1884_phase_peaks':phasepeaks,'phase_peak_nodes':phaseat,'normal_success138535_prefix138658_both_below150000':True,
  'prefix_scope':'Failure during leaf decode may retain16152 transient and add3*work codec references:138658. Child decode includes2048 omitted in original pernode formula but max126982. No capraise/code rewrite; completed manifest stream bytes always preserved. For fixed source every legitimate failure codec remains below150k even with held children.',
  'heights':'M824 source, codec qi numerator28/denominator164, G<=M: valuehash1676, restore1648,ratio1016 <=1921<2048. Allinteger node update weights and partialbounds verified before values.',
  'RAW_one_generation':one,'RAW_both_generations_and_manifests_controls':generated,'fresh_RAW_upper':B['fresh_RAW_upper'],'whole_review_numeric_upper':review_prior,
  'scope':'Source/code/custody and unrelated codec fixtures ONLY. ROOTGO separate. No chosen projection coefficient/fullentry/rank/pivot or matrix/RHS update claim.'}
 ep=out/'independent-code-lifetime-fixed-q-codec-and-unrelated-fixture-proof.json';save(ep,ev)
 rep={'id':ID,'all_passed':True,'frozen':True,'SOURCE_CODE_SYNTHETIC_ONLY':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};save(out/'report.json',rep)
 assert(out/'report.json').stat().st_size<=65536 and count(rep)<=2000;signal.alarm(0)
 print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'SOURCE_review_owned':review_prior,'META':baseline+control,'nominal138535_prefix138658':True,'seconds':rep['seconds']}))
if __name__=='__main__':main()
