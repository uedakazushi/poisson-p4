#!/usr/bin/env python3
"""UNRUN until ROOT ONE GO: fixed lower-cache B63/plane10 composition."""
import argparse,hashlib,json,re,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SOURCE='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_V2_20261105_v63'
STEM='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63'
RAW=0;MAXPRE=0;MAXSTORED=0;LIMIT=2048;OWNEDMAX=0;ROOTONLY=2048
def sha(b):return hashlib.sha256(b).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(tuple,list)):return sum(map(count,x))
 if isinstance(x,dict):return count(list(x.keys()))+count(list(x.values()))
 return 0
def hb(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def pre(h):
 global MAXPRE
 MAXPRE=max(MAXPRE,h);assert h<=LIMIT and h<=1153 and h<=2048,('ROOT_ONLY_PREOP_STOP',h,LIMIT)
def keep(q):
 global MAXSTORED
 MAXSTORED=max(MAXSTORED,hb(q));assert hb(q)<=384 and hb(q)<=2048;return q
def parse(s):
 if isinstance(s,int):pre(abs(s).bit_length());return keep(Q(s))
 assert isinstance(s,str)and re.fullmatch(r'-?\d+(?:/\d+)?',s),('QQ_PARSE',s)
 a=s.split('/');n=int(a[0]);d=int(a[1])if len(a)>1 else 1
 assert d;pre(max(abs(n).bit_length(),abs(d).bit_length()));return keep(Q(n,d))
def add(a,b):pre(2*max(hb(a),hb(b))+1);return keep(a+b)
def times(a,b):pre(hb(a)+hb(b));return keep(a*b)
def divide(a,b):assert b;pre(hb(a)+hb(b));return keep(a/b)
def owned(n):
 global OWNEDMAX
 OWNEDMAX=max(OWNEDMAX,n);assert n<=100000 and n<=150000,('ROOT_ONLY_OWNED_STOP',n)
def readbytes(path):
 global RAW
 b=path.read_bytes();RAW+=len(b);assert RAW<=2000000 and RAW<=33554432,('ROOT_ONLY_RAW_STOP',RAW);return b
def load(path):return json.loads(readbytes(path))
def selective(path,fields,expected_hash=None):
 b=readbytes(path);s=b.decode();decoder=json.JSONDecoder();out={}
 if expected_hash is not None:assert sha(b)==expected_hash,('IMMUTABLE_INPUT_HASH',str(path))
 for key in fields:
  matches=list(re.finditer(r'"'+re.escape(key)+r'"\s*:',s));assert len(matches)==1,(str(path),key)
  j=matches[0].end()
  while s[j].isspace():j+=1
  out[key]=decoder.raw_decode(s,j)[0]
 return out
def write(path,data):
 global RAW
 b=(json.dumps(data,separators=(',',':'))+'\n').encode();RAW+=len(b);assert RAW<=2000000 and RAW<=33554432
 path.write_bytes(b);return sha(b)
def encode(a):return[str(q)for q in a]
def main():
 global LIMIT,ROOTONLY
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--source-peer-report',required=True)
 ap.add_argument('--head30-source-peer-report',required=True)
 args=ap.parse_args();started=time.monotonic();signal.alarm(60)
 source_dir=ROOT/'reports'/SOURCE
 sourcebytes=readbytes(source_dir/'report.json');sr=json.loads(sourcebytes);sourcehash=sha(sourcebytes);del sourcebytes
 pp=source_dir/'positive-before-actual-twohead-cache-and-whole-budget.json';pb=readbytes(pp)
 assert sha(pb)==sr['evidence_sha256'][pp.name];profile=json.loads(pb);del pb
 pb=readbytes(ROOT/args.source_peer_report);peerhash=sha(pb);peer=json.loads(pb);del pb
 pb=readbytes(ROOT/args.head30_source_peer_report);headpeerhash=sha(pb);headpeer=json.loads(pb);del pb
 headpeerpath='reports/AMBIENT_B18_C_CONNECTING_HEAD30_SOURCE_20261105_v63_REVIEW/report.json'
 assert headpeerhash==sr['input_sha256'][headpeerpath]
 assert sr['all_passed']and peer['all_passed']and headpeer['all_passed']
 assert sr['status']=='SOURCE_TEN_ROOTS_MIXED63_PLANE10_BOTH_HEADS_POSITIVE_NO_VALUES'
 assert sr['actual_root_or_edge_values_evaluated']is False
 assert profile['all_new_output_values_UNCOMPUTED']is True
 cp=source_dir/'compiled-ten-roots-mixed63-plane10-and-actual-provided-lower-inputs.json'
 raw=readbytes(cp);assert sha(raw)==sr['evidence_sha256'][cp.name];compiled=json.loads(raw);del raw
 assert profile['future_static_owned_upper']<=75000 and max(profile['SOURCE_phase_numeric_upper'])<150000
 ROOTONLY=profile['root_before_operation_bound'];LIMIT=ROOTONLY
 assert 0<LIMIT<=1153 and profile['post_before_operation_bound']<=1153
 assert profile['committed_QQ_bits']==384 and profile['temporary_QQ_bits']==1153
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 phase={'phase':'pre-arithmetic-fixed-mixed-precision-source-pinned','source_report_sha256':sourcehash,
  'source_peer_report_sha256':peerhash,'head30_source_peer_sha256':headpeerhash,
  'code_sha256':sha(readbytes(Path(__file__))),'ROOT_ONE_GO_external_to_worker_required':True,
  'source_static_upper':profile['future_static_owned_upper'],'hard_value_owned':100000,'root_preop':LIMIT,
  'post_preop':profile['post_before_operation_bound'],'hard_QQ_bits':2048,'no_new_implicit_source_gauge_or_jet':True}
 write(out/'phase.json',phase)
 try:
  mon63=[tuple(m)for m in compiled['basis63_monomials']];mon10=[tuple(m)for m in compiled['plane10_monomials']]
  mons={63:mon63,10:mon10};mi={p:{m:i for i,m in enumerate(mm)}for p,mm in mons.items()}
  i10=[mi[63][m]for m in mon10];products={}
  for p,mm in mons.items():
   products[p]=[[i,j,mi[p][tuple(sorted(a+b))]]for i,a in enumerate(mm)for j,b in enumerate(mm)if tuple(sorted(a+b))in mi[p]]
  records=compiled['original_records_precision_and_retirement'];requests=compiled['odd_trusted_root_requests']
  groots=compiled['delta_safe_root_pair'];oroots=[r['original_node']for r in requests];roots=groots+oroots
  # Convert in place so parsed decimal strings and Fraction coefficient
  # arrays do not remain as two full105x63 numeric copies.
  given={}
  for node,arr in compiled['provided105_lower_extension63']:
   for i,s in enumerate(arr):arr[i]=parse(s)
   given[node]=arr
  def restrict(v,p):return v if len(v)==p else[v[j]for j in i10]
  def conv(a,b,p):
   v=[Q(0)]*p
   for i,j,k in products[p]:v[k]=add(v[k],times(a[i],b[j]))
   return v
  cache={};outputs={};unitchecks=[];valuehash=hashlib.sha256();actual_peak=0;actual_QQ=0
  for r,p,gone in records:
   i,op=r[:2]
   if op=='input':v=restrict(given[i],p)
   elif op=='c':v=[Q(0)]*p;v[0]=divide(parse(r[2]),parse(r[3]))
   elif op=='a':v=[add(a,b)for a,b in zip(restrict(cache[r[2]],p),restrict(cache[r[3]],p))]
   elif op=='n':v=[-q for q in restrict(cache[r[2]],p)]
   elif op=='m':v=conv(restrict(cache[r[2]],p),restrict(cache[r[3]],p),p)
   elif op=='i':
    a=restrict(cache[r[2]],p);assert a[0]==1;v=[Q(0)]*p;v[0]=divide(Q(1),a[0])
    for k in range(1,p):
     s=Q(0)
     for x,y,t in products[p]:
      if t==k and x>0 and y<k:s=add(s,times(a[x],v[y]))
     v[k]=-divide(s,a[0])
    one=[Q(1)]+[Q(0)]*(p-1);assert conv(a,v,p)==one==conv(v,a,p);unitchecks.append([i,r[2],'1',p])
   else:raise AssertionError(op)
   assert len(v)==p;cache[i]=v;actual_peak=max(actual_peak,len(cache));actual_QQ=max(actual_QQ,sum(map(len,cache.values())))
   valuehash.update((json.dumps([i,p,encode(v)],separators=(',',':'))+'\n').encode())
   if i in roots:outputs[i]=restrict(v,63 if i in groots else 10).copy()
   for j in gone:del cache[j]
   owned(sum(profile['future_runtime1_numeric_components'].values()))
  assert not cache and actual_peak<=profile['lastuse_peak_nodes']and actual_QQ<=profile['lastuse_peak_QQ_slots']
  assert[[i,len(outputs[i])]for i in roots]==compiled['captured_output_precisions']
  assert[[i,j,u]for i,j,u,p in unitchecks]==compiled['unit_inverse_centers']
  delta=[add(a,-b)for a,b in zip(outputs[groots[0]],outputs[groots[1]])]
  checks=[]
  for req in requests:
   v=outputs[req['original_node']]
   for m,s in req['root_lower_checks']:
    q=parse(s);assert v[mi[10][tuple(m)]]==q,(req['physical_row'],m,'ROOT_LOWER_OVERLAP');checks.append([req['physical_row'],m,s])
  assert delta[mi[63][(7,7,7,7)]]==parse(compiled['known_delta63_plane4_r4_r3ell'][0])
  assert delta[mi[63][(7,7,7,8)]]==parse(compiled['known_delta63_plane4_r4_r3ell'][1])
  node63=next(r['original_node']for r in requests if r['physical_row']==63)
  assert outputs[node63][mi[10][(7,8,8)]]==parse(compiled['known_bB63_r_ell2'])
  assert delta[0]==0 and all(delta[mi[63][(i,)]]==0 for i in range(9))
  assert MAXSTORED<=profile['root_stored_bound']
  first={'scope':'Actual entire safe delta63 B63 difference, and trusted bB plane10 ONLY; no individual even B63 truth claim.',
   'safe_delta63_coefficients':encode(delta),'basis63_monomials':compiled['basis63_monomials'],
   'actual_bB_plane10_roots':[[r['physical_row'],r['original_node'],encode(outputs[r['original_node']])]for r in requests],
   'plane10_monomials':compiled['plane10_monomials'],'all_root_lower_overlap_checks':checks,
   'known_delta_plane4_two_and_known_bB63_lambda_checked':True,'all_unit_both_products':unitchecks,
   'all_node_chosen_lower_extension_value_hash':valuehash.hexdigest(),
   'unknown_Phi_socle_NOT_asserted_zero':'The chosen lower extension carries zeros there. Safe delta truth is by zero-gradient SOURCE theorem; bB truth only trusted plane10.',
   'no_new_physical_IFT_gauge_or_jet':True,'root_phase_max_pre':MAXPRE,'root_phase_max_stored':MAXSTORED,
   'actual_lastuse_peak_nodes':actual_peak,'actual_lastuse_QQ_slots':actual_QQ}
  p1=out/'phase1-safe-delta63-and-eight-trusted-bB-cache-readbacks.json';p1hash=write(p1,first)
  # Mandatory lifetime boundary. A post failure preserves phase1; this
  # worker is never automatically relaunched and never recomputes phase1.
  headpaths=compiled['head29_head30_paths'];saved_basis=compiled['basis63_monomials'];saved_plane=compiled['plane10_monomials']
  rowroots={r['physical_row']:r['original_node']for r in requests}
  del compiled,records,given,cache,products,unitchecks,checks,first,requests
  LIMIT=profile['post_before_operation_bound'];leaf={};headoperators=[]
  for path in headpaths:
   h=selective(ROOT/path,['known_cached_inputs','four_genuine_and_canonical_high_polynomials_UNPAIRED'],sr['input_sha256'][path]);headoperators.append(h)
   for label,s,locator in h['known_cached_inputs']:
    key=tuple(label);q=parse(s)
    if key in leaf:assert leaf[key]==q
    leaf[key]=q
  for key,mon in[(('NEW_delta63_plane4','r4'),(7,7,7,7)),(('NEW_delta63_plane4','r3ell'),(7,7,7,8)),(('NEW_delta63_plane4','r2ell2'),(7,7,8,8))]:leaf[key]=delta[mi[63][mon]]
  leaf[('NEW_CACHED_bB63_ell3',)]=outputs[node63][mi[10][(8,8,8)]]
  for m,i in mi[63].items():leaf[('SAFE_delta63',)+m]=delta[i]
  for row,node in rowroots.items():leaf[('TRUSTED_plane3','bB',row,1,2)]=outputs[node][mi[10][(7,8,8)]]
  values=[];terms=[]
  for h in headoperators:
   vv=[];tt=[]
   for poly in h['four_genuine_and_canonical_high_polynomials_UNPAIRED']:
    v=Q(0);t=[]
    for labels,weight in poly:
     q=parse(weight)
     for label in labels:q=times(q,leaf[tuple(label)])
     v=add(v,q);t.append(str(q))
    vv.append(v);tt.append(t)
   assert vv[2]==times(parse('5/3'),vv[0])and vv[3]==add(vv[1],times(parse('-121/48'),vv[0]))
   values.append(vv);terms.append(tt)
  owned(sum(profile['future_runtime2_numeric_components'].values()))
  assert MAXSTORED<=profile['post_partial_stored_bound']
  second={'scope':'Full coherent HEAD29/30 in SAME original connecting map, genuine highA/highB and canonical rg/k2. No other head highs or original H0/nil inference.',
   'head29_A_B_rg_k2':encode(values[0]),'head30_A_B_rg_k2':encode(values[1]),
   'all_eight_full_post_termwise_products':terms,'canonical_transformation_checked':True,
   'four_true_parameters_for_Kbar_source_basis':{'p2':str(values[0][2]),'b2':str(values[0][3]),'p3':str(values[1][2]),'b3':str(values[1][3])},
   'still_uncomputed_head0_head4_highs':['p0','b0','p4','b4'],
   'Kbar_basis_scope':'Existing formal source basis I4+rh3+u,s now has its four actual parameters. Full R0 nil lifts and Kbar actions/extension/H0 require separate source proof.',
   'no_new_source_response_solve':True,'phase1_sha256':p1hash}
  p2=out/'phase2-actual-full-head29-head30-and-four-Kbar-parameters.json';p2hash=write(p2,second)
  rep={'id':STEM,'all_passed':True,'status':'ACTUAL_CACHE_ONLY_FULL_HEAD29_HEAD30_FOUR_CANONICAL_HIGHS',
   'source_report_sha256':phase['source_report_sha256'],'source_peer_report_sha256':phase['source_peer_report_sha256'],
   'head30_source_peer_sha256':phase['head30_source_peer_sha256'],'code_sha256':phase['code_sha256'],
   'head29_A_B_rg_k2':encode(values[0]),'head30_A_B_rg_k2':encode(values[1]),
   'new_physical_source_IFT_gauge_or_jet':False,'full_Kbar_actions_H0_nil_inferred':False,
   'evidence_sha256':{p1.name:p1hash,p2.name:p2hash},'maximum_preoperation_bits':MAXPRE,
   'maximum_stored_QQ_bits':MAXSTORED,'whole_owned_numeric_upper':OWNEDMAX,'fresh_RAW_before_report':RAW,'seconds':time.monotonic()-started}
  write(out/'report.json',rep);print(json.dumps(rep));signal.alarm(0)
 except BaseException as e:
  failure={'all_passed':False,'exception':repr(e),'phase1_exists':(out/'phase1-safe-delta63-and-eight-trusted-bB-cache-readbacks.json').exists(),
   'preserve_phase1_no_automatic_retry':True,'maximum_preoperation_bits':MAXPRE,'maximum_stored_bits':MAXSTORED,'owned_upper':OWNEDMAX,'RAW':RAW}
  (out/'FAILURE-preserve-no-retry.json').write_text(json.dumps(failure,indent=2)+'\n');raise
if __name__=='__main__':main()
