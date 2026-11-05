#!/usr/bin/env python3
"""One gated existing-cache three-pivot ell^3 scan with the full six-coordinate HB matrix in QQ[ell]/ell^4."""
import argparse,hashlib,json,re,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63'
SOURCE='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/report.json'
SOURCE_SHA='6cbc0f67ea1b88f6d5e14d06176f38b58d876b34afb5823a45c73bc78a4023d0'
COMPILED='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/compiled-six-pivot-even-roots-and-pure-ell4-inputs.json'
PROFILE='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/full-HB-pivot03-frame-and-whole-operational-profile.json'
FORMAL='reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/full-HB6-by2-affine-source-operator-and-old27-leaves.json'
MAX_PRECHECK=MAX_VALUE=0

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def check(n):
 global MAX_PRECHECK
 MAX_PRECHECK=max(MAX_PRECHECK,n);assert n<=2048
def note(q):
 global MAX_VALUE
 MAX_VALUE=max(MAX_VALUE,abs(q.numerator).bit_length(),q.denominator.bit_length());assert MAX_VALUE<=2048
 return q
def parse(x):
 p=str(x).split('/');assert len(p)in[1,2]
 check(max(abs(int(v)).bit_length()for v in p));return note(Q(str(x)))
def add(a,b):
 if not a:return b
 if not b:return a
 check(max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a+b)
def times(a,b):
 if not a or not b:return Q(0)
 check(max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length()))
 return note(a*b)
def neg(a):return note(-a)
def conv(a,b):
 v=[Q(0)]*4
 for n in range(4):
  for k in range(n+1):v[n]=add(v[n],times(a[k],b[n-k]))
 return v
def inverse(x):
 assert x[0]==1
 v=[Q(1),Q(0),Q(0),Q(0)]
 for n in range(1,4):
  for k in range(1,n+1):v[n]=add(v[n],neg(times(x[k],v[n-k])))
 return v

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-peer-report',required=True);args=parser.parse_args()
 signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648));started=time.monotonic()
 assert sha(ROOT/SOURCE)==SOURCE_SHA
 source=json.loads((ROOT/SOURCE).read_text());peerpath=Path(args.source_peer_report)
 if not peerpath.is_absolute():peerpath=ROOT/peerpath
 peer=json.loads(peerpath.read_text());assert peer['all_passed']and SOURCE_SHA in peer['input_sha256'].values()
 for p in[COMPILED,PROFILE,FORMAL]:assert sha(ROOT/p)==source['evidence_sha256'][Path(p).name]
 c=json.loads((ROOT/COMPILED).read_text());profile=json.loads((ROOT/PROFILE).read_text());formal=json.loads((ROOT/FORMAL).read_text())
 records=c['operator_records_originalID_op_args'];meta=c['root_pairs_physical3']
 assert len(records)==7095 and[r['physical_row']for r in meta]==[18,60,67]
 roots=[r[n]for r in meta for n in['aA_node','aB_node']]
 assert roots==[20114,40892,20138,40916,20153,40931]
 owned=count(source)+count(peer)+count(c)+count(profile)+count(formal)+2*len(records)+525+4*162+16+14000
 assert owned<=75000<100000<150000
 raw=sum(p.stat().st_size for p in[ROOT/SOURCE,ROOT/COMPILED,ROOT/PROFILE,ROOT/FORMAL,peerpath]);assert raw+150000<=1000000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 (out/'pre-arithmetic-source-and-whole-budget-custody.json').write_text(json.dumps({
  'source_sha256':SOURCE_SHA,'source_peer_sha256':sha(peerpath),'code_sha256':sha(Path(__file__)),
  'owned_numeric_upper':owned,'RAW_input_bytes':raw,'hard_bits':2048,'hard_seconds':60,
  'requested_positive_root_values_not_yet_evaluated':True},indent=2)+'\n')
 phi={i:list(map(parse,v))for i,v in c['provided105_originalID_pure_ell4']};prior=dict(c['provided105_originalID_pure_ell3_prior'])
 assert len(phi)==105 and all(v[:3]==list(map(parse,prior[i]))for i,v in phi.items())
 uses=Counter()
 for r in records:
  if r[1]not in['c','input']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 cache={};peak=0;inverses=[]
 try:
  for r in records:
   i,op=r[:2]
   if op=='input':v=phi[i]
   elif op=='c':v=[parse(str(r[2])+'/'+str(r[3]))]+[Q(0)]*3
   elif op=='a':v=[add(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
   elif op=='n':v=[neg(a)for a in cache[r[2]]]
   elif op=='m':v=conv(cache[r[2]],cache[r[3]])
   elif op=='i':
    v=inverse(cache[r[2]]);one=[Q(1),Q(0),Q(0),Q(0)]
    assert conv(cache[r[2]],v)==one==conv(v,cache[r[2]]);inverses.append([i,r[2],'1'])
   else:raise AssertionError(op)
   cache[i]=v;peak=max(peak,len(cache))
   if op not in['c','input']:
    for k in r[2:]:
     uses[k]-=1
     if not uses[k]:del cache[k]
   if not uses[i]:del cache[i]
  assert peak==162 and inverses==c['unit_inverse_centers']
  expected=dict(c['expected6_roots_prior_pure_ell3'])
  assert all(cache[i][:3]==list(map(parse,expected[i]))for i in roots)
  results=[]
  for row in meta:
   a,b=cache[row['aA_node']],cache[row['aB_node']];delta=[add(x,neg(y))for x,y in zip(a,b)]
   assert delta[:3]==[Q(0)]*3
   results.append({'physical_row':row['physical_row'],'aA_node':row['aA_node'],'aB_node':row['aB_node'],
    'aA_pure_ell4':list(map(str,a)),'aB_pure_ell4':list(map(str,b)),
    'delta_pure_ell4':list(map(str,delta)),'delta_ell3':str(delta[3])})
 except BaseException as err:
  (out/'ROOT_ONLY_STOP_CHECKPOINT.json').write_text(json.dumps({'status':'STOP_no_cap_raise','error':repr(err),'completed_original_node':i,
   'max_value_bits':MAX_VALUE,'max_precheck_bits':MAX_PRECHECK,'already_available_root_values':[[j,list(map(str,cache[j]))]for j in roots if j in cache]},indent=2)+'\n');raise
 assert MAX_PRECHECK<=profile['root_phase_before_operation_bound']
 root_pre,root_bits=MAX_PRECHECK,MAX_VALUE
 phase1={'all6_pivot_roots_pure_ell4':[[i,list(map(str,cache[i]))]for i in roots],
  'all3_pivot_mu03_readbacks':results,'all315_input_overlaps':True,'all18_prior_root_overlaps':True,
  'both_inverse_products':inverses,'lastuse_peak_nodes':peak,'max_root_precheck_bits':root_pre,'max_root_value_bits':root_bits}
 rp=out/'phase1-six-cached-pivot-pure-ell4-values.json';rp.write_text(json.dumps(phase1,indent=2)+'\n')
 try:
  leaves={tuple(k):parse(v)for k,v in formal['known27_old_cached_leaf_values']}
  for r in results:leaves[('mu03',r['physical_row'],-1)]=parse(r['delta_ell3'])
  matrix=[];products=[]
  for a,row in enumerate(formal['full_HB_6_by2_affine_operators']):
   vals=[]
   for b,entry in enumerate(row):
    value=Q(0)
    for k,w in entry:
     product=times(parse(w),leaves[tuple(k)]);value=add(value,product)
     products.append([a,b,k,w,str(leaves[tuple(k)]),str(product)])
    vals.append(value)
   matrix.append(vals)
  assert matrix[0]==[parse('-3/5'),parse('6')]and matrix[3]==[parse('-121/80'),parse('121/8')]
  w=[add(row[1],times(parse('10'),row[0]))for row in matrix]
  assert w[0]==w[3]==0
  nz=[j for j in range(1,6)if w[j]]
  if nz:
   j=nz[0];block=[matrix[0],matrix[j]]
   det=add(times(block[0][0],block[1][1]),neg(times(block[0][1],block[1][0])))
   assert det==times(parse('-3/5'),w[j])and det
   # Guard BEFORE constructing the reciprocal Fraction; its height is
   # bounded by the proved generic det height4H+1.
   check(max(abs(det.numerator).bit_length(),det.denominator.bit_length()))
   invdet=note(Q(det.denominator,det.numerator))
   inv=[[times(block[1][1],invdet),times(neg(block[0][1]),invdet)],
        [times(neg(block[1][0]),invdet),times(block[0][0],invdet)]]
   for A,B in[(block,inv),(inv,block)]:
    actual=[[add(times(A[a][0],B[0][b]),times(A[a][1],B[1][b]))for b in range(2)]for a in range(2)]
    assert actual==[[1,0],[0,1]]
   rank=2;cert={'rows':[0,j],'minor':str(det),'inverse':[[str(x)for x in row]for row in inv],
    'both_inverse_products_verified':True,'all2_column_factorization':'The two source columns are the fixedbasis; coefficient matrixI2, every full6row retained.'}
  else:
   assert all(row[1]==times(parse('-10'),row[0])for row in matrix)
   rank=1;cert={'nonzero_row':0,'nonzero_entry':'-3/5','column1_equals_minus10_column0_all6rows':True,
    'factorization_full6_rows_by1_then1_by2':{'left':[str(row[0])for row in matrix],'right':['1','-10']}}
  assert MAX_PRECHECK<=profile['post_affine_rank_inverse_operational_bound']
 except BaseException as err:
  (out/'POST_ROOT_ONLY_STOP_CHECKPOINT.json').write_text(json.dumps({'status':'POST_ROOT_STOP_no_cap_raise_no_root_rerun','error':repr(err),
   'saved_phase1_sha256':sha(rp),'max_value_bits':MAX_VALUE,'max_precheck_bits':MAX_PRECHECK},indent=2)+'\n');raise
 data={'actual_full_HB6_by2_matrix':[[str(x)for x in row]for row in matrix],
  'all_termwise_affine_products':products,'col1_plus10col0_all6rows':list(map(str,w)),
  'actual_known_full_HB_rank':rank,'fixed_rank_certificate':cert,
  'actual_full_C_HF4_if_rank2':4 if rank==2 else'UNDETERMINED_4_OR5',
  'Tor_m2C_type_if_rank2':[2,2]if rank==2 else'[2,b>=2]_UNDETERMINED',
  'long_C_torsion_exponent_if_rank2':'3_OR4'if rank==2 else'e>=3_OPEN',
  'rank1_scope':'ThisisONLYtheknownTWOHBsubmodule. No completequarticrank/exactHF4/e claimifrank1.',
  'all_r_ell3_and_ell4_coordinates_retained':True,'no_literal_ell3H_cut':True,
  'new_physical_IFT_gauge_jet_delta4_mixed3_beta3_coefficients':0,
  'owned_numeric_upper':owned,'RAW_input_bytes':raw,'whole_max_value_bits':MAX_VALUE,'whole_max_precheck_bits':MAX_PRECHECK,
  'source_root_precheck_bound':profile['root_phase_before_operation_bound'],'source_post_precheck_bound':profile['post_affine_rank_inverse_operational_bound']}
 ep=out/'phase2-full-six-coordinate-HB-matrix-and-fixed-rank-certificate.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'FULL_HB_KNOWN_RANK2_ACTUAL_HF4_4'if rank==2 else'FULL_HB_KNOWN_RANK1_SUBMODULE_ONLY',
  'actual_known_full_HB_rank':rank,'actual_full_C_HF4_determined':rank==2,'exact_long_C_torsion_exponent_determined':False,
  'input_sha256':{SOURCE:SOURCE_SHA,COMPILED:sha(ROOT/COMPILED),PROFILE:sha(ROOT/PROFILE),FORMAL:sha(ROOT/FORMAL),str(peerpath.relative_to(ROOT)):sha(peerpath)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)for p in[rp,ep]},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 assert raw+sum(p.stat().st_size for p in out.iterdir())<=1000000
 print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
