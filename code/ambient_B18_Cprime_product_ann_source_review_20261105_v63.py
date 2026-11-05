#!/usr/bin/env python3
"""Independent fixed-source annihilator/quotient check; no value/rank search."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_CPRIME_PRODUCT_ANN_SOURCE_20261105_v63_REVIEW'
OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def mat(v):return [[F(c)for c in r]for r in v]
def mul(a,b):return [[sum((x*y for x,y in zip(r,c)),F(0))for c in zip(*b)]for r in a]
def eye(n):return [[F(i==j)for j in range(n)]for i in range(n)]
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',expected=None):
  p=ROOT/'reports'/s/n;h=sha(p);pins[str(p.relative_to(ROOT))]=h
  if expected:assert h==expected
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 def checkpins(d,key='input_sha256'):
  for p,h in d.get(key,{}).items():assert sha(ROOT/p)==h;pins[p]=h
 A=load('AMBIENT_B18_CPRIME_PRODUCT_ANN_INITIAL_SUPPLEMENT_20261105_v63',expected='2007bc30234c073232db934177b55ba62ce9cb87a069225238d4da854a8f9c95');checkpins(A)
 P=load('AMBIENT_B18_CPRIME_PRODUCT_ANN_QUOTIENT_SOURCE_PREFLIGHT_20261105_v63',expected='82ded8e376d3a31d743f10cbc3a47fcd653224a8590f811dc6fa23bd3c98ae74');checkpins(P)
 oldpeer=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63_REVIEW',expected='4108f86bd3e122b92fd0488ec970fac0a38fae494194b4454617b906d3a91450');assert oldpeer['all_passed']
 hist=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_HISTORY_20261105_v63',expected='4538b10299661a210b86043ce1d702a21058ebe11f8b020d7a0a06a9909f44ef');checkpins(hist)
 hp=ROOT/hist['contentaddressed_prepared_source'];assert sha(hp)==hist['prepared_source_sha256']=='1be2cda944e1dfc0a608e0ac03035168438f4240abc279ad34739146c474b028';pins[str(hp.relative_to(ROOT))]=sha(hp)
 S=load('AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63','actual-special-fiber20-coefficients.json')
 for rec in S['input_hashes'].values():assert sha(ROOT/rec['path'])==rec['sha256'];pins[rec['path']]=rec['sha256']
 F9=load('AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63','full9-firstjet-certificate.json')
 def first(polys):
  out=[]
  for p in polys:
   d={tuple(e):F(c)for e,c in p};out.append([d.get(tuple(int(i==j)for i in range(9)),F(0))for j in range(9)])
  return out
 ca,cb=map(first,[S['physical_coordinate_maps2']['CA'],S['physical_coordinate_maps2']['CB']]);cbi=mat(S['physical_CB_linear_inverse']);assert mul(cb,cbi)==mul(cbi,cb)==eye(9)
 Y=mul(cbi,ca);assert Y==mat(S['YB_in_YA_linear'])==mat(F9['YB_to_YA_linear'])
 assert F9['full9_CA_coordinate_order']==S['source_F9_ids']['CA']
 a=[r[:7]for r in Y[:7]];ai=mat(A['actual_inverse']);assert a==mat(A['actual_full9_CB_zero7_linear_matrix'])==mat(S['A0']);assert mul(a,ai)==mul(ai,a)==eye(7)
 assert all(r[7:]==[F(0),F(0)]for r in Y[:7])
 for family in ['CA','CB']:assert [v[1]for v in S['source_conductor_generators_rows_and_ids'][family][:7]]==S['source_F9_ids'][family][:7]
 # A0^-1 B has linear part z_k on FULL nine coordinates, before ell=0.
 liftfirst=mul(ai,Y[:7]);assert liftfirst==[r+[F(0),F(0)]for r in eye(7)]
 p2=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase2-known-quadratic-closure.json');cols=[pairs(v)for v in p2['all396_known_quadratic_columns135_sparse']];cp=p2['fixed_pivot_columns'];rr=p2['fixed_pivot_rows'];iv=mat(p2['fixed_minor_inverse']);assert len(cp)==len(rr)==len(iv)==131
 free=[i for i in range(135)if i not in rr];assert free==P['unit_NF_rows135']==[44,132,133,134]
 q=[]
 for f in free:
  w=[sum((cols[cp[k]].get(f,F(0))*iv[k][j]for k in range(131)),F(0))for j in range(131)]
  row={f:F(1)}
  for j,c in enumerate(w):
   if c:row[rr[j]]=-c
  q.append(row)
 assert q==[pairs(v)for v in P['literal_fixed_pi4_rows_by135_sparse']]
 assert all(sum((c*col.get(i,F(0))for i,c in row.items()),F(0))==0 for row in q for col in cols)
 assert [[row.get(i,F(0))for i in free]for row in q]==eye(4)
 mon2=list(itertools.combinations_with_replacement(range(9),2));ix={m:i for i,m in enumerate(mon2)};ann=[{45*t+ix[(i,j)]:F(1)}for t in range(3)for i,j in itertools.combinations_with_replacement(range(7),2)]
 assert ann==[pairs(v)for v in A['actual_initial_columns84_by135_sparse']]
 qs=set().union(*map(set,q));ans=set().union(*map(set,ann));assert qs.isdisjoint(ans)
 assert qs=={42,43,44,132,133,134};assert {mon2[i%45]for i in qs}=={(7,7),(7,8),(8,8)}
 projected=[[sum((c*col.get(i,F(0))for i,c in row.items()),F(0))for col in ann]for row in q];assert projected==[[F(0)]*84 for _ in range(4)]
 # Membership is formal: IA IB kills coker(D0->E4), hence Cprime.
 # A0^-1 B in IB has first z_k; zi in IA. No tau-origin substitution.
 save('independent-fullnine-source-lifts-and-annihilator-proof.json',dict(actual_CA_F9=S['source_F9_ids']['CA'],actual_CB_F9=S['source_F9_ids']['CB'],actual_conductor_labels=S['source_conductor_generators_rows_and_ids'],YB_in_YA_full9=Y,A0=a,A0_inverse=ai,all49_two_sided_inverse_entries=True,separate_r_s_zero=True,all63_fullnine_lift_first_entries=liftfirst,ann_initial_columns84=ann,proof='Each conductor scalar gives a factorization Bi -> Ai -> Bi of multiplication by that scalar. After the common P0 tensor base changes each comparison cokernel is killed by its conductor; the composite cokernel C is killed by IA*IB. Cprime is a quotient of C, so each zi*(A0_inverse B)_k e_t is an actual full-R0 relation. Its degree2 initial is zi*zk e_t. The trace tau coordinates are not specialized to their origins.',not_Told_or_nil_annihilator_claim=True))
 save('independent-fixed-pi4-and-disjoint-support-zero.json',dict(fixed_pi4_rows=q,free_NF_indices=free,quadratic_monomial_order=mon2,all396_known_columns_killed=True,all16_embedding_identities=True,ann_target_support=sorted(ans),pi4_support=sorted(qs),all336_projected_entries=projected,added_rank_exactly0=True,known_quad_rank_still131=True,known_quad_quotient_dimension_still4=True,conditional_cubic_h2_le3_not_triggered=True,not_actual_HF2_or_exponent_negative=True))
 assert A['whole_budget_with_supplement_upper']==P['future_whole_QQ_bound']==88555
 assert A['extra_cubic_source_count']==84*9==756 and P['future_cubic_source_count']==44*45+756+3==2739 and P['future_cubic_shape_max']==[27,2766]
 # This replay has 17161 inverse entries, <=2616 sparse first/known data,
 # 4x135 projection workspace, and fixed small formal-source tables.
 slots=17161+sum(map(len,cols))+sum(map(len,ann))+sum(map(len,q))+4*135+4*84+4000;assert slots<88555<150000
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 raw=sum(p.stat().st_size for p in OUT.rglob('*')if p.is_file());assert raw<32<<20 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,scope='SOURCE polynomial lifts and fixed quotient identities only',all84_actual_ann_initial_relations_accepted=True,all336_fixed_projection_entries_zero=True,added_quadratic_rank_exactly0=True,known_quadratic_quotient_still4=True,conditional_cubic_not_triggered=True,no_actual_V2_exponent_nilcube_or_newsource_claim=True,repository_source_history_preserved=True,QQ_slots_upper=slots,primary_whole_bound=88555,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=raw,seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
 print(json.dumps(dict(all_passed=True,added_rank=0,known_quotient_dimension=4,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
