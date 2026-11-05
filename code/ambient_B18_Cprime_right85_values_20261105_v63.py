#!/usr/bin/env python3
"""One cached value job on the complete85 first-linear-syzygy block."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];GATE='AMBIENT_B18_CPRIME_RIGHT85_SOURCE_PREFLIGHT_20261105_v63';STEM='AMBIENT_B18_CPRIME_RIGHT85_VALUES_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main(peer_path):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={}
 def read(path):p=ROOT/path;pins[path]=sha(p);return json.loads(p.read_text())
 gp=f'reports/{GATE}/report.json';r=read(gp);peer=read(peer_path);assert r['all_passed']and peer['all_passed']and peer['input_sha256'][gp]==sha(ROOT/gp)
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h; pins[p]=h
 path=f'reports/{GATE}/complete-right85-and-cached-plane-source-adjoints.json';e=read(path);assert sha(ROOT/path)==r['evidence_sha256'][Path(path).name]
 qp=next(p for p in r['input_sha256']if p.endswith('phase4-all-twenty-actual-quadratic-values.json'));Q2=read(qp)['all44_actual_V2_pureplane_nine_values'];leaves=[Q(c)for col in Q2 for c in col]
 dp=next(p for p in r['input_sha256']if p.endswith('actual-curve-source-and-one-residual-cubic-functional.json'));dual={i:Q(c)for i,c in read(dp)['residual_actual_cubic_dual495_sparse']}
 mons=list(itertools.combinations_with_replacement(range(9),3));mi={m:i for i,m in enumerate(mons)};piv=e['actual_V1_fixed24_source_columns'];values=[];native=[]
 for syz,adj in zip(e['all85_syzygy_columns216_sparse'],e['all85_cached_plane_Q2_adjoints396_sparse']):
  value=sum((Q(c)*leaves[i]for i,c in adj),Q(0));poly={}
  for i,c in syz:
   j,var=divmod(i,9)
   if var<7:continue
   for target in range(3):
    for k,mon in enumerate([(7,7),(7,8),(8,8)]):
     row=165*target+mi[tuple(sorted(mon+(var,)))];poly[row]=poly.get(row,Q(0))+Q(c)*leaves[9*piv[j]+3*target+k]
  poly={i:c for i,c in poly.items()if c};assert value==sum((c*dual.get(i,Q(0))for i,c in poly.items()),Q(0));values.append(value);native.append(poly)
 assert len(values)==85
 for v in values:assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=2048
 nonzero=[i for i,v in enumerate(values)if v];positive=bool(nonzero);data={'all85_actual_residual_cubic_values':values,'all85_actual_plane_cubic_initial_representatives_sparse':[list(p.items())for p in native],
  'all85_adjoint_and_actual_cubic_pairings_equal':True,'all_actual_initial_membership_by_V1_syzygies':True,'source_kernel_sha256':e['kernel_sha256'],'nonzero_indices':nonzero,'first_nonzero_source':nonzero[0]if positive else None,
  'known26_plus_actual_right85_rank':27 if positive else 26,'positive_Cprime_ell3_length9_and_long_exponent2':positive,'zero_is_only_partial_negative':True,'new_physical_source_gauge_or_jet':False}
 blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert sum((ROOT/p).stat().st_size for p in pins)+len(blob.encode())<32<<20;ep=out/'all85-actual-cached-cubic-values-and-source-membership.json';ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_cached_complete_right85_cubic_values','actual_right85_rank':int(positive),'positive_Cprime_ell3_length9_and_long_exponent2':positive,'no_exponent_conclusion_if_all_zero':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source-peer-report',required=True);main(p.parse_args().source_peer_report)
