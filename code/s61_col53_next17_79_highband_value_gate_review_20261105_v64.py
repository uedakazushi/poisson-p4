#!/usr/bin/env python3
"""Independent HIGH-only SOURCE/code/custody proof; no signed HIGH VALUE."""
import ast,hashlib,importlib.util,itertools,json,math,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
from s61_col53_next17_79_streamed_value_gate_review_20261105_v64 import own_decode,own_encode
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
ROOT=pathlib.Path(__file__).resolve().parents[1]
GATE='S61_COL53_NEXT17_79_FULL_HL_HIGHBAND_VALUE_GATE_20261105_v64';ID=GATE+'_REVIEW'
WORKER='code/s61_col53_next17_79_full_HL_highband_value_20261105_v64.py'
WANT='6bba032b7a22fe827539739fd9b496f0c968eec843a52fbd71b5db56ce8357fa'
FINAL='reports/S61_COL53_NEXT17_79_FINALIZED_SOURCE_20261105_v64/compiled-full-SLP-and-fixed-projection-plan.json'
FULL='reports/S61_COL53_NEXT17_79_FULL_HL_SUPPORT_20261105_v64'
OLD='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64'
LOW='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_20261105_v64/all-projected-SOURCE-node-counts.json'
PROPOSAL='reports/S61_COL53_NEXT17_79_FULL_HL_EXTENSION_BUDGET_SOURCE_20261105_v64/source-cauchy-two-band-lifetime-and-raw-proposal.json'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def numeric_legacy(x):
 if isinstance(x,bool):return 0
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):return 3 if re.fullmatch(r'[+-]?\d+(?:/\d+)?',x) else 0
 if isinstance(x,list):return sum(map(numeric_legacy,x))
 if isinstance(x,dict):return sum(map(numeric_legacy,x.values()))
 return 0
def pack(e):
 assert len(e)==5 and all(-64<=x<=191 for x in e)
 return sum((x+64)<<(8*j) for j,x in enumerate(e))
def unpack(k):return tuple(((k>>(8*j))&255)-64 for j in range(5))
def uint_size(n):return max(1,(n.bit_length()+6)//7)
def sint_size(bits):
 n=(bits+7)//8;return uint_size(2*n+1)+n
def main():
 start=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=False):
  nonlocal raw
  pp=ROOT/p;s=sha(pp)
  if h:assert h==s,(p,h,s)
  if p not in pins:raw+=pp.stat().st_size
  pins[p]=s;return json.loads(pp.read_text()) if parse else pp
 gp=f'reports/{GATE}/report.json';g=pin(gp,'29c95901ebdd16b70e3c374f37c9fdac734369d71edcedf5fa3ef147a6a2cd86',True);assert g['all_passed']
 for p,h in g['input_sha256'].items():pin(p,h)
 assert pins[WORKER]==WANT
 bp=f'reports/{GATE}/highband-whole-lifetime-and-storage-budget.json';B=pin(bp,g['evidence_sha256'][pathlib.Path(bp).name],True)
 plan=f'reports/{GATE}/incremental-node-plan.jsonl';pin(plan,g['evidence_sha256'][pathlib.Path(plan).name])
 cr=pin('reports/S61_COL53_NEXT17_79_FULL_HL_HIGHBAND_CUSTODY_20261105_v64/report.json','5ed3646ced395b844a62c5993782513a71eac7589637310fb4d0a07d5033286a',True);assert cr['all_passed']
 ce='reports/S61_COL53_NEXT17_79_FULL_HL_HIGHBAND_CUSTODY_20261105_v64/whole-actual-read-and-generation-custody.json'
 custody=pin(ce,cr['evidence_sha256'][pathlib.Path(ce).name],True)
 for p,h in custody['fixed_actual_READ_sha256'].items():pin(p,h)
 assert pin(OLD+'_LOCAL_IDENTITY_REVIEW/report.json','2b0b2212d91ad5f74481ffc5428a9c33bd1d6cfb018a78c7a3cea7130747dd88',True)['all_passed']
 fullreport=pin(FULL+'/report.json','4488f4a9425da0c0c81b93b5546411f9f03edceaaf7767be19ed4a6277350ace',True)
 for n,h in fullreport['evidence_sha256'].items():pin(FULL+'/'+n,h)
 D=json.loads((ROOT/FINAL).read_text());assert numeric_legacy(D)==65701
 # Retire only an unused in-memory table. The actual input file stays intact.
 assert numeric_legacy(D.pop('active_tied_HL_children'))==3078
 assert numeric_legacy(D)==62623
 ns=D['nodes'];caps=D['needed_node_degree'];order=sorted(map(int,caps));N=len(order);assert N==1884
 F=dict(json.loads((ROOT/FULL/'all1884-exact-complete-face-source-counts.json').read_text()))
 L=dict(json.loads((ROOT/LOW).read_text()));A={};oldsize=0
 with(ROOT/OLD/'node-values-index.jsonl').open() as f:
  for line in f:
   assert len(line)<=2048;r=json.loads(line);k=r['node'];assert k not in A
   assert Q(r['fixed_SOURCE_q'])==Q(ns[k]['q']) and r['fixed_SOURCE_lead_UNCHANGED']==ns[k]['lead5']
   assert r['stream_offset']==oldsize;oldsize+=r['bytes'];A[k]=r['actual_nonzero_slots']
 assert oldsize==(ROOT/OLD/'node-values.q5stream').stat().st_size==859628
 assert sorted(F)==sorted(L)==sorted(A)==order and all(A[k]<=L[k]<=F[k] for k in order)
 H={k:F[k]-L[k] for k in order};V={k:A[k]+H[k] for k in order}
 P=json.loads((ROOT/PROPOSAL).read_text());records=P['all_node_records'];assert len(records)==N
 prior=count((g,B,cr,custody,D,F,L,A,H,V,P,pins))+25000
 assert prior<150000 and raw+2000000<32*1024**2
 save(out/'positive-before-SOURCE-code-and-unrelated-fixtures.json',{'all_passed':True,'SOURCE_only':True,'numeric_upper':prior,'RAW_upper':raw+2000000,'no_signed_HIGH_values':True})
 phases={};ats={};prefix={};onegen=0;weights_n=0;maxbits=maxweight=maxnormalization=0;maxline=0
 spec=importlib.util.spec_from_file_location('prepared_highband_fixture_only',ROOT/WORKER);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 assert mod.__name__!='__main__'
 def bound(name,k,v):
  if v>phases.get(name,0):phases[name]=v;ats[name]=k
 def pbound(name,k,v):prefix[name]=max(prefix.get(name,0),v)
 def digits(bits):return max(1,(bits*30103+99999)//100000)
 with(ROOT/plan).open() as pf:
  for rec,k in zip(records,order):
   n=ns[k];args=n['args'];q=Q(n['q']);M=int(n['M']);lead=n['lead5'];cap=caps[str(k)]
   ch=list(dict.fromkeys(args if n['op'] in ('copy','multiply') else [c for c in args if ns[c]['q']!='0' and ns[c]['lead5'][:2]==lead[:2]]))
   assert all(c<k and c in F for c in ch) and len(lead)==5
   w=[];partial=0
   if not q:pass
   elif n['op']=='input':partial=M
   elif n['op']=='copy':
    c=args[0];ratio=Q(ns[c]['q'])*Q(n['scalar'])/q;assert ratio.denominator==1;w=[ratio.numerator]
    partial=abs(w[0])*int(ns[c]['M'])
   elif n['op']=='multiply':
    a,b=args;ratio=Q(ns[a]['q'])*Q(ns[b]['q'])*Q(n['scalar'])/q;assert ratio.denominator==1;w=[ratio.numerator]
    partial=min(F[a],F[b])*int(ns[a]['M'])*int(ns[b]['M'])*abs(w[0])
   else:
    assert n['op']=='add'
    for j,c in enumerate(args):
     if c not in ch:continue
     ratio=Q(ns[c]['q'])*(Q(n['scalar']) if j else Q(1))/q;assert ratio.denominator==1;w.append(ratio.numerator)
     partial+=abs(ratio.numerator)*int(ns[c]['M'])
   assert partial<=M and rec['integer_weights']==list(map(str,w))
   assert [rec['node'],rec['new_positive_slots'],rec['full_actual_slots_upper'],rec['old_actual_slots'],rec['full_positive_slots'],rec['new_codec_RAW_upper']]==json.loads(pf.readline())
   assert (rec['new_positive_slots'],rec['old_actual_slots'],rec['full_positive_slots'])==(H[k],A[k],F[k])
   assert rec['full_actual_slots_upper']==V[k]
   if n['op'] in ('copy','multiply'):
    assert [sum(ns[c]['lead5'][j] for c in args)+n['shift5'][j] for j in range(5)]==lead
   maxbits=max(maxbits,M.bit_length());maxweight=max(maxweight,*[abs(x).bit_length() for x in w],0);weights_n+=len(w)
   ct=sum(V[c] for c in ch);mc=max([max(A[c],H[c]) for c in ch]+[0]);hh=H[k]
   bound('arithmetic',k,118000+7*(ct+hh)+4096)
   bound('child_decode',k,118000+7*ct+8*mc+2048+4096)
   bound('codec_released_children',k,118000+10*hh+4096)
   if n['op']=='input' and q:
    bound('leaf_decode',k,118000+16152+7*hh+4096);pbound('leaf_STOP',k,118000+16152+10*hh+4096)
   if k==D['FULL_F11_original_root']:
    bound('anchor_released_children',k,118000+7*hh+8*411+2048+4096)
    pbound('anchor_STOP',k,118000+10*hh+8*411+2048+4096)
   # A child-decode failure has produced NO current output term; its
   # decoder stack unwinds before handler, which clears child functions.
   pbound('child_STOP_conservative_empty_output',k,118000+8*mc+2048+4096)
   rawcell=9 if not hh else 9+sint_size(abs(q.numerator).bit_length()+M.bit_length())+sint_size(q.denominator.bit_length())+10+uint_size(hh)+hh*(sint_size(M.bit_length())+10)
   assert rawcell==rec['new_codec_RAW_upper'];onegen+=rawcell
   gg=10**digits(M.bit_length())-1
   rr={'terms':hh,'G':gg,'content':'-'+'9'*digits(abs(q.numerator).bit_length()+M.bit_length())+'/'+'9'*digits(q.denominator.bit_length()),'shift':[-191]*5,'value_sha256':'f'*64,'sha256':'f'*64,'bytes':rawcell}
   frame=mod.frame_record(k,10676690,rr,hh,hh,859628,859628,A[k],'f'*64,q,cap,lead,gg)
   maxline=max(maxline,len(mod.enc(frame)))
  assert not pf.read(1)
 assert maxbits==824 and maxweight==31 and 2*824+31<1921<2048
 assert onegen==10676690 and maxline==1749<=2048
 assert sum(H.values())==191169 and sum(A.values())==67463 and F[D['minor_B_root']]==1012 and H[D['minor_B_root']]==883
 assert phases=={'arithmetic':140639,'child_decode':143570,'codec_released_children':130926,'leaf_decode':138542,'anchor_released_children':128629}
 assert max(phases.values())==143570 and max(prefix.values())==138668<150000
 generated=2*onegen+2048*N+524288+1048576;assert generated==B['generated_RAW_upper']==26784676
 assert custody['fixed_actual_READ_bytes']==sum((ROOT/f).stat().st_size for f in custody['fixed_actual_READ_sha256'])==3953454
 assert custody['fixed_actual_READ_bytes']+generated+2*65536+262144<=B['fresh_RAW_upper']==31647203<32*1024**2
 census={'runtime':25000,'resident_graph_after_unused_active_retirement':62623,'remaining_numeric_mode_keys':N,'order':N,'combined_index_KEY_old3_high3_TWO_hash_numerics':9*N,'current_gate_plan_codec_root_control':4096}
 assert census==B['explicit_census'] and sum(census.values())==112443 and 112443+4000+1557==118000
 del P,records
 # Independent local proof of all complete POSITIVE source supports.
 # This uses only original input support and existing SOURCE key streams;
 # no signed high coefficient is ever computed.
 ix={};offset=0
 with(ROOT/FULL/'source-node-index.jsonl').open() as f:
  for line in f:
   rec=json.loads(line);k=rec['node'];assert rec['offset']==offset and rec['positive_terms']==F[k] and rec['bytes']==5*F[k]
   ix[k]=(rec['offset'],rec['bytes'],rec['positive_terms'],rec['sha256']);offset+=rec['bytes']
 assert sorted(ix)==order and offset==(ROOT/FULL/'all-node-positive-supports.q5keys').stat().st_size
 local_prior=numeric_legacy(D)+count((g,B,cr,custody,F,L,A,H,V,ix,pins))+25000+21000+10000
 assert local_prior<150000
 def leaf(n):
  s=n['source'];ref=s.get('legacy_source',s);shift=[0]*5;scale=Q(1);terms=[]
  if 'binding' in ref:
   b=ref['binding'];p=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'/b['scope']
   doc=json.loads(p.read_text());v=(doc['common_operands'] if 'common_operands' in doc else doc['row_operands'])[b['label']]
   p=v['original_literal'];base=p['shift'];terms=[(tuple(base[j]+int(t[j+1]) for j in range(5)),int(t[0])) for t in p['terms']]
   scale=Q(v['QQ_transform_scalar']);shift=v['transform_monomial_exponent5'];content=Q(p['content'])
  elif 'thin_payload_file' in ref or ('file' in ref and 'offset' in ref):
   filename=ref.get('thin_payload_file',ref.get('file'))
   with(ROOT/filename).open('rb') as f:f.seek(ref['offset']);body=f.read(ref['bytes'])
   assert hashlib.sha256(body).hexdigest()==ref['payload_sha256'];content,base,terms=own_decode(b'Q5CELL1\n\1'+body)
  elif 'cell_file' in ref or 'file' in ref:
   filename=ref.get('cell_file',ref.get('file'));body=(ROOT/filename).read_bytes();content,base,terms=own_decode(body)
   if 'legacy_source' in s:
    prof=s['literal_profile'];low=[min(e[j] for e,c in terms) for j in range(5)];shift=[prof['exponent_box5'][j][0]-low[j] for j in range(5)]
    scale=Q(prof['known_QQ_content'])/content
    assert (scale,shift)==((Q(353894400000),[1,0,0,0,1]) if 'col56-ninth-cleared' in filename else (Q(1),[0]*5))
  elif ref.get('QQ_monomial'):
   prof=s.get('literal_profile');e=ref.get('exponent5',[b[0] for b in prof['exponent_box5']] if prof else [0]*5);terms=[(tuple(e),1)];content=Q(n['q'])
  else:assert ref.get('known_zero') or ref.get('known_source_zero');content=Q(0)
  assert content*scale==Q(n['q'])
  return {pack(tuple(e[j]+shift[j] for j in range(5))) for e,c in terms if c and tuple(e[j]+shift[j] for j in range(2))==tuple(n['lead5'][:2])}
 pair_count=0;support_hash=hashlib.sha256();sourcepeak=local_prior
 with(ROOT/FULL/'all-node-positive-supports.q5keys').open('rb') as sf:
  def load(k):
   off,nb,card,want=ix[k];sf.seek(off);body=sf.read(nb);assert len(body)==nb and hashlib.sha256(body).hexdigest()==want
   keys=[int.from_bytes(body[j:j+5],'big') for j in range(0,nb,5)]
   assert all(a<b for a,b in zip(keys,keys[1:]));assert all(unpack(z)[:2]==tuple(ns[k]['lead5'][:2]) and unpack(z)[2]>=ns[k]['lead5'][2] for z in keys)
   return set(keys)
  for k in order:
   n=ns[k];expected=load(k);assert len(expected)==F[k]
   assert sum(unpack(z)[2]-n['lead5'][2]<=caps[str(k)] for z in expected)==L[k]
   if n['q']=='0':actual=set();children=[]
   elif n['op']=='input':actual=leaf(n);children=[]
   else:
    wanted=n['args'] if n['op'] in ('copy','multiply') else [c for c in n['args'] if ns[c]['q']!='0' and ns[c]['lead5'][:2]==n['lead5'][:2]]
    memo={c:load(c) for c in dict.fromkeys(wanted)};children=list(memo.values());actual=set()
    if n['op']=='copy':
     for z in memo[n['args'][0]]:actual.add(pack(tuple(x+y for x,y in zip(unpack(z),n['shift5']))))
    elif n['op']=='multiply':
     a,b=n['args']
     for z in memo[a]:
      e=unpack(z)
      for zz in memo[b]:
       f=unpack(zz);actual.add(pack(tuple(e[j]+f[j]+n['shift5'][j] for j in range(5))));pair_count+=1
    else:
     for v in children:actual.update(v)
   assert actual==expected,(k,len(actual),len(expected))
   sourcepeak=max(sourcepeak,numeric_legacy(D)+count((g,B,cr,custody,F,L,A,H,V,ix,pins))+25000+len(actual)+len(expected)+sum(map(len,children))+4096)
   assert sourcepeak<150000
   support_hash.update((str(k)+' '+str(len(actual))+' '+hashlib.sha256(b''.join(z.to_bytes(5,'big') for z in sorted(actual))).hexdigest()+'\n').encode())
 assert pair_count==494806
 # Frozen code operational conditions, not only the gate's prose.
 text=(ROOT/WORKER).read_text();ast.parse(text)
 for token in ["retired=count_numeric(D.pop('active_tied_HL_children'))","index[k]=(oo,ob,A,off,rec['bytes'],len(work),oh,rec['sha256'])",'transient=8*nc+2048','own(new=True)','memo.clear();held.clear();own()',"mathG=0 if not work else rec['G']","assert A==0","assert sha(oldfile)==oldsha and sha(oldindex)==oldindexsha",'assert time.monotonic()-start<50',"signal.alarm(55)"]:
  assert token in text,token
 release=text.index('memo.clear();held.clear();own()');assert release<text.index("if k==D['FULL_F11_original_root']:anchor_high")<text.index("context['phase']='new_highband_codec_commit'")
 handler=text[text.index('except Exception as exc:'):];assert handler.index('memo.clear();held.clear()')<handler.index('own(extra=3*len(work)+4096)')
 assert 'mixed_low_times_low_HIGH_contributions_retained' in text and "for ai,(e,v)in enumerate(memo[aa].items())" in text
 assert 'whole_minor_ZERO_NOT_inferred' in text and 'new_matrix_RHS_update_NOT_performed' in text
 # Unrelated negative-q / different-G / empty / low*low->high fixtures.
 fix=out/'unrelated-highband-fixtures';fix.mkdir();qi=Q(-17,33)
 low={(2,-1,0,4,1):12,(2,-1,2,4,1):-18};high={(2,-1,3,4,1):35,(2,-1,5,4,1):-70};merged={};receipts=[]
 for name,vals in [('low',low),('high',high),('empty',{})]:
  p=fix/(name+'.bin');rec=write_cell(p,qi,vals);assert validate_cell(p,qi,vals,rec);body=p.read_bytes();assert body==own_encode(qi,vals)
  ownc,sh,terms=own_decode(body)
  if vals:
   G=ownc/qi;assert G>0 and G.denominator==1
   restored={e:c*G.numerator for e,c in terms};assert restored==vals
   assert dict(mod.fixed_q_terms(body,qi,len(vals),1000))==vals
   for e,c in restored.items():assert e not in merged;merged[e]=c
  else:assert list(mod.fixed_q_terms(body,Q(0),0,1000))==[] and rec['G']==1 and rec['content']=='0'
  receipts.append({'name':name,'G':rec['G'],'sha256':sha(p)})
 assert [r['G'] for r in receipts[:2]]==[6,35]
 full={};hi={};mixed=False
 for e,c in merged.items():
  for f,d in merged.items():
   z=tuple(x+y for x,y in zip(e,f));full[z]=full.get(z,0)+c*d
   if z[2]>1:
    hi[z]=hi.get(z,0)+c*d
    if e in low and f in low:mixed=True
 assert mixed and hi=={e:c for e,c in full.items() if e[2]>1}
 assert not (ROOT/'reports/S61_COL53_NEXT17_79_FULL_HL_HIGHBAND_VALUE_20261105_v64').exists()
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 ev={'all_passed':True,'SOURCE_CODE_ONLY_no_signed_HIGH_coefficients':True,'all1884_source_integer_weights_and_full_Cauchy_bounds':True,
  'all1884_complete_positive_supports_rebuilt_from_original_leaves_and_saved_source_children':True,'all_low_positive_supports_are_literal_cap_restrictions':True,
  'all1884_source_support_identity_digest':support_hash.hexdigest(),'positive_product_pairs':pair_count,'fixed_OLD_actual_slots_metadata':67463,'new_high_positive_slots':191169,'complete_root_positive_slots':1012,'root_high_positive_slots':883,
  'full_function_induction':'For each node, OLD(low<=owncap)+NEW(high>owncap) is disjoint and complete on its SOURCE leadingHL face. Products use BOTH complete child bands and retain low*low high-output terms. Sums use only SOURCE-HL tied nonzero-q branches, never prune actual vertex cancellation. All5 coordinates remain in future VALUE7 maps.',
  'fixed_q_and_codec_proof':'Each band primitive content equals qi*Gband; restore its OWN positive integer G (possibly differing between low/high), then merge. Nonempty sign(content)=sign(qi), empty mathG0/codecG1 skips qi division. Input col56 transform353894400000 and[1,0,0,0,1] independently rebound.',
  'all1884_metadata_census_9N_index':census,'census112443_plusdocs4000_margin1557_equals118000':True,
  'normal_phase_maxima':phases,'normal_peak_nodes':ats,'all_guard_prefix_bounds':prefix,'normal143570_prefix138668_below150000':True,
  'SOURCE_M_max_bits':maxbits,'integer_weight_max_bits':maxweight,'positive_temporary_bound':2*824+31,'temporary_guard1921_hard2048':True,
  'one_NEW_generation_RAW':onegen,'TWO_NEW_generations_and_manifest_controls':generated,'fixed_actual_READ_bytes':3953454,'fresh_RAW_upper':31647203,'manifest_line_upper':maxline,
  'old_low_stream_and_index_immutable_and_only_READ':True,'child_retirement_before_normal_and_STOP_codec':True,'math50_signal55_hard60_AS2e9':True,
  'unrelated_codec_and_mixed_band_fixtures':receipts,'review_SOURCE_prior':prior,'review_source_local_support_upper':sourcepeak,'review_RAW_upper':raw+2000000,
  'ROOT_GO_and_actual_values_still_pending':True,'no_whole_minor_zero_or_rank_or_matrix_update':True}
 ep=out/'independent-highband-source-code-completeness-lifetimes-and-codec-proof.json';save(ep,ev)
 report={'id':ID,'all_passed':True,'SOURCE_CODE_CUSTODY_ONLY':True,'input_sha256':pins,'code_sha256':sha(pathlib.Path(__file__)),
         'evidence_sha256':{ep.name:sha(ep)},'worker_sha256':WANT,'normal143570_prefix138668':True,'whole_SOURCE_upper':max(prior,sourcepeak),'new_signed_HIGH_values_UNCOMPUTED':True,'seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'whole_source':max(prior,sourcepeak),'seconds':report['seconds']}))
if __name__=='__main__':main()
