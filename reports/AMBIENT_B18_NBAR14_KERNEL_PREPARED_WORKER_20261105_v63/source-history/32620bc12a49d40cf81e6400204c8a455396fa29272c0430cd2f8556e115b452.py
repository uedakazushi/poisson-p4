#!/usr/bin/env python3
"""Prepared fivecachededgefunctionals and genuine14kernel; beta remainsformal."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];GATE='AMBIENT_B18_NBAR14_KERNEL_SOURCE_PREFLIGHT_20261105_v63';STEM='AMBIENT_B18_NBAR14_KERNEL_CACHED_VALUES_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def rec(v):return {i:Q(c)for i,c in v if Q(c)}
def count(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def guard(v):
 if isinstance(v,Q):assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=2048
 elif isinstance(v,dict):
  for x in v.values():guard(x)
 elif isinstance(v,(list,tuple)):
  for x in v:guard(x)
 return v
def main(peer_path):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={}
 def read(path):p=ROOT/path;pins[path]=sha(p);return json.loads(p.read_text())
 gp=f'reports/{GATE}/report.json';r=read(gp);peer=read(peer_path);assert r['all_passed']and peer['all_passed']and peer['input_sha256'][gp]==sha(ROOT/gp)
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h; pins[p]=h
 ep=f'reports/{GATE}/five-cache-functionals-and-genuine-fourteen-kernel-source.json';e=read(ep);assert sha(ROOT/ep)==r['evidence_sha256'][Path(ep).name]
 arow=rec(e['lambda0_source_row44_sparse']);c0=rec(e['head42_origin_cycle44_sparse']);rho=rec(e['remaining_pure_r_first_cycle_adjoint44_sparse']);s1=[rec(v)for v in e['all4_other_head_pure_r_first_cycles44_sparse']]
 terms=[[[i,c,c0.get(i,Q(0)),c*c0.get(i,Q(0))]for i,c in sorted(arow.items())]]+[[[i,c,v.get(i,Q(0)),c*v.get(i,Q(0))]for i,c in sorted(rho.items())]for v in s1]
 vals=guard([sum((row[-1]for row in term),Q(0))for term in terms]);a=vals[0];b=vals[1:];assert a;xi=Q(-11,16)*a;assert xi
 free=[1,2,3,4,5,7,8,9,10,11,12,13,14,15];K=[[Q(0)]*14 for _ in range(16)]
 for j in range(4):K[j+1][j]=Q(1);K[6][j]=-b[j]/xi
 for j,i in enumerate(free[4:],4):K[i][j]=Q(1)
 assert all(K[i][j]==Q(k==j)for k,i in enumerate(free)for j in range(14))and K[0]==[Q(0)]*14
 knownlong=[Q(0)]+b+[Q(0),xi]+[Q(0)]*9;short=[a]+[Q(0)]*15;assert len(knownlong)==16
 assert all(sum((short[i]*K[i][j]for i in range(16)),Q(0))==0 and sum((knownlong[i]*K[i][j]for i in range(16)),Q(0))==0 for j in range(14))
 det=a*xi;assert det
 connecting=next(p for p in r['input_sha256']if p.endswith('actual-formal-connecting16-and-ten-source-obstructions.json'))
 operator={'nodes':[['actual_external_matrix','CONNECTING16','matrix_roots.actual_ell0_cycles44_by16'],['QQ_matrix_rows',K],['matrix_product',[0,1]]],'shapes':[[44,16],[16,14],[44,14]],'root':2,
  'namespace_bindings':{'CONNECTING16':{'path':connecting,'sha256':pins[connecting]}},'proof':'Accepted16sourcecycles timesactualQQkernelbasis; exactconnecting identifiesthese14barclasseswithNil(D0)/ellNil(D0). Noexpansionofscalarfunctionsor780relations.'}
 data={'five_actual_cached_functional_values_a_b1_b2_b3_b4':vals,'all_five_exact_termwise_cache_products':terms,'xi_equals_minus11a_by16':xi,
  'canonical_long_coordinate_beta_of_head0':'FORMAL_UNEVALUATED','actual_edge_pivot2symbolic_rows':[[a,Q(0)],['beta',xi]],'pivot_determinant_independent_beta':det,
  'actual_edge_kernel_basis16_by14_rows':K,'retained_unit_free_rows14':free,'free_tail_I14_verified':True,'both_edge_rows_kill_all14_independently_of_beta':True,'genuine_Nbar_dimension':14,
  'genuine_Nbar14_barcycle_operator':operator,'canonical_head0_long_edge_not_needed_or_computed':True,'no_new_physical_source_gauge_jet_or_formal_nil_lifts':True,'no_original_nilcube_or_product_values':True}
 guard(data);owned=count([r,peer,e,arow,c0,rho,s1,terms,vals,K,operator,data]);assert owned<=60000,owned
 blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert sum((ROOT/p).stat().st_size for p in pins)+len(blob.encode())<32<<20;ev=out/'actual-five-edge-functionals-and-genuine-Nbar14-cycles.json';ev.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_cached_actual_Nbar14_kernel_and_barcycles','actual_five_functionals':list(map(str,vals)),'actual_xi':str(xi),'genuine_Nbar_dimension':14,'canonical_head0_long_edge_beta_computed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'owned_numeric_slots':owned,'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source-peer-report',required=True);main(p.parse_args().source_peer_report)
