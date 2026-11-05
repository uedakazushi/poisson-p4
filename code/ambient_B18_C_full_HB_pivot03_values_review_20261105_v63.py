#!/usr/bin/env python3
"""Independent fixed-cache pivot jets, reversed affine sums and full HB factor."""
import hashlib,json,pathlib,resource,signal,time
from collections import Counter
from fractions import Fraction as Q
import ambient_B18_C_full_eta4_cache_values_review_20261105_v63 as J
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='8a1c4d8e0a5a32c2ab4ff88af5af60478a4d83bf591588196aefae91a33e18c5'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None):
  nonlocal raw
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:raw+=p.stat().st_size
  pins[k]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 docs={p:pin(ROOT/p,h) for p,h in pr['input_sha256'].items()}
 peer=next(d for p,d in docs.items() if p.endswith('_REVIEW/report.json'))
 assert peer['full_source_operator_code_custody_and_whole_budget_PASS']
 assert sha(ROOT/'code/ambient_B18_C_full_HB_pivot03_values_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
 own_prior=pin(ROOT/'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63_REVIEW/report.json',
  '8dcf7359fe70c248ece84a61bca323f0820b1086e565d3687b23b59b7b70ca1b')
 pin(ROOT/'code/ambient_B18_C_full_eta4_cache_values_review_20261105_v63.py',own_prior['code_sha256'])
 P=pin(ROOT/'reports'/PRI/'phase1-six-cached-pivot-pure-ell4-values.json',
   '2bd0f488a412226a324fd2eb8fa6631995e9d075926cf1c2933ba4bc20376f7b')
 V=pin(ROOT/'reports'/PRI/'phase2-full-six-coordinate-HB-matrix-and-fixed-rank-certificate.json',
   '744bf26ef7be4f728b75476a7f2713fe41ae0d6584a41a74f61f68f97d176c2e')
 c=next(d for p,d in docs.items() if p.endswith('compiled-six-pivot-even-roots-and-pure-ell4-inputs.json'))
 F=next(d for p,d in docs.items() if p.endswith('full-HB6-by2-affine-source-operator-and-old27-leaves.json'))
 rec=c['operator_records_originalID_op_args'];meta=c['root_pairs_physical3'];roots=[r[k] for r in meta for k in ['aA_node','aB_node']]
 assert len(rec)==7095 and roots==[20114,40892,20138,40916,20153,40931]
 owned=sum(J.nc(v) for v in docs.values())+J.nc(P)+J.nc(V)+J.nc(pr)+2*len(rec)+525+4*162+16+14000
 assert owned<75000 and raw<2<<20
 save(out/'pre-arithmetic-independent-cache-and-affine-bound.json',{'owned_upper':owned,'RAW_read':raw,'hard_bits':2048,
  'hard_numeric':100000,'global_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2<<30,
  'method':'own dictionary ell4 jets closed inverse; reversed affine products; fixed rank1 source-factor check'})
 J.premax=J.valmax=0;phi={};prior=dict(c['provided105_originalID_pure_ell3_prior'])
 for i,v in c['provided105_originalID_pure_ell4']:
  a=list(map(J.parse,v));assert a[:3]==list(map(J.parse,prior[i]));phi[i]={k:q for k,q in enumerate(a) if q}
 assert len(phi)==105
 uses=Counter()
 for r in rec:
  if r[1]not in['input','c']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 cache={};peak=0;inv=[]
 for r in rec:
  i,op=r[:2]
  if op=='input':a=phi[i]
  elif op=='c':q=J.parse(str(r[2])+'/'+str(r[3]));a={0:q} if q else {}
  elif op=='a':a=J.plus(cache[r[2]],cache[r[3]])
  elif op=='n':a=J.neg(cache[r[2]])
  elif op=='m':a=J.times(cache[r[2]],cache[r[3]])
  elif op=='i':
   a=J.inv(cache[r[2]]);assert J.times(a,cache[r[2]])=={0:Q(1)}==J.times(cache[r[2]],a);inv.append([i,r[2],'1'])
  else:raise AssertionError(op)
  cache[i]=a;peak=max(peak,len(cache))
  if op not in ['c','input']:
   for k in r[2:]:
    uses[k]-=1
    if not uses[k]:del cache[k]
  if not uses[i]:del cache[i]
 assert peak==162 and inv==c['unit_inverse_centers']==P['both_inverse_products']
 assert [[i,J.dense(cache[i])] for i in roots]==P['all6_pivot_roots_pure_ell4']
 old=dict(c['expected6_roots_prior_pure_ell3'])
 assert all(J.dense(cache[i])[:3]==old[i] for i in roots)
 mu={};deltas=[]
 for r in meta:
  d=J.plus(cache[r['aA_node']],J.neg(cache[r['aB_node']]))
  assert all(not d.get(i,Q(0)) for i in range(3));mu[r['physical_row']]=d.get(3,Q(0))
  saved=next(v for v in P['all3_pivot_mu03_readbacks'] if v['physical_row']==r['physical_row'])
  assert J.dense(d)==saved['delta_pure_ell4']
  deltas.append([r['physical_row'],J.dense(d)])
 assert mu=={18:Q(22,675),60:Q(659,5400),67:Q(1,675)}
 leaves={tuple(k):J.parse(q) for k,q in F['known27_old_cached_leaf_values']}
 leaves.update({('mu03',p,-1):q for p,q in mu.items()})
 matrix=[];products={}
 for i,row in enumerate(F['full_HB_6_by2_affine_operators']):
  vals=[]
  for j,entry in enumerate(row):
   v=Q(0)
   for k,cq in reversed(entry):
    p=J.mul(J.parse(cq),leaves[tuple(k)]);v=J.add(v,p);products[(i,j,tuple(k))]=p
   vals.append(v)
  matrix.append(vals)
 savedproducts={(i,j,tuple(k)):J.parse(p) for i,j,k,cq,lv,p in V['all_termwise_affine_products']}
 assert products==savedproducts
 assert [[str(q) for q in row] for row in matrix]==V['actual_full_HB6_by2_matrix']
 assert matrix[0]==[-Q(3,5),Q(6)] and matrix[3]==[-Q(121,80),Q(121,8)]
 w=[J.add(row[1],J.mul(Q(10),row[0])) for row in matrix]
 assert not any(w) and V['col1_plus10col0_all6rows']==['0']*6
 assert all(row[1]==J.mul(-Q(10),row[0]) for row in matrix)
 assert pr['actual_known_full_HB_rank']==V['actual_known_full_HB_rank']==1
 assert not pr['actual_full_C_HF4_determined'] and not pr['exact_long_C_torsion_exponent_determined']
 cert=V['fixed_rank_certificate']['factorization_full6_rows_by1_then1_by2']
 assert cert=={'left':[str(row[0]) for row in matrix],'right':['1','-10']}
 proof={'all7095_original_cache_ops_independently_replayed':True,
  'own_dictionary_jet_closed_inverse_reverse_products':True,
  'all315_input_and18_prior_root_overlaps':True,'both_unit_inverse_products':inv,'lastuse_peak':peak,
  'all6_root4_coeffs_match':True,'all3_actual_mu03':[[p,str(mu[p])] for p in [18,60,67]],
  'all12_affine_entries_and_every_reversed_product_match':True,
  'actual_full_HB_matrix':[[str(q) for q in row] for row in matrix],
  'all6_w_exact_zero':True,'fixed_factorization_rank':1,
  'rank1_only_known_two_HB_sources':True,'actual_full_HF4_still_4_or5':True,
  'b_and_e_remain_undetermined':True,'no_ell3H_cut_or_new_source_physical_gauge_jet':True,
  'own_preoperation_max_bits':J.premax,'own_stored_max_bits':J.valmax,
  'whole_owned_upper':owned,'RAW_read':raw,'hard_bits':2048,'hard_numeric':100000,'global_numeric':150000}
 ep=out/'independent-pivot03-cache-full-affine-and-rank1-factor-certificate.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'full_HB_known_rank1_accepted':True,
  'no_actual_HF4_or_e_claim':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
