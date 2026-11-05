#!/usr/bin/env python3
"""One fixed cached curve-column pairing; no source/gauge/jet computation."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
GATE='AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_SOURCE_PREFLIGHT_20261105_v63'
STEM='AMBIENT_B18_CPRIME_CURVE_GA_CUBIC_VALUES_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(v):
 if isinstance(v,Q):return str(v)
 if isinstance(v,dict):return {k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def number_count(v):
 if isinstance(v,(Q,int)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(number_count(x)for x in v.values())
 if isinstance(v,(list,tuple)):return sum(map(number_count,v))
 return 0
def main(source_peer,closure_peer):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));resource.setrlimit(resource.RLIMIT_FSIZE,(32<<20,32<<20));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={}
 def read(path):p=ROOT/path;pins[path]=sha(p);return json.loads(p.read_text())
 gp=f'reports/{GATE}/report.json';r=read(gp);assert r['all_passed'];peer=read(source_peer);cp=read(closure_peer);assert peer['all_passed']and cp['all_passed']
 assert peer['input_sha256'][gp]==sha(ROOT/gp)
 for path,h in r['input_sha256'].items():assert sha(ROOT/path)==h; pins[path]=h
 path=f'reports/{GATE}/actual-curve-source-and-one-residual-cubic-functional.json';e=read(path);assert sha(ROOT/path)==r['evidence_sha256'][Path(path).name]
 closure_path=next(p for p in r['input_sha256']if p.endswith('phase2-complete-streamed-cubic-upper-certificate.json'));assert cp['input_sha256'][closure_path]==sha(ROOT/closure_path);X=json.loads((ROOT/closure_path).read_text());assert X['known_cubic_rank']==26
 dual={i:Q(c)for i,c in e['residual_actual_cubic_dual495_sparse']};curve={i:Q(c)for i,c in e['actual_curve_relation495_sparse']};q27={i:Q(c)for i,c in e['residual_compressed_dual27_sparse']};proj=[{i:Q(c)for i,c in row}for row in X['all495_source_projection_columns27_sparse']]
 assert all(sum((c*col.get(i,Q(0))for i,c in q27.items()),Q(0))==0 for col in [{i:Q(c)for i,c in v}for v in X['selected_original_compressed_columns']])
 terms=[[i,c,dual.get(i,Q(0)),c*dual.get(i,Q(0))]for i,c in sorted(curve.items())];value=sum((v[-1]for v in terms),Q(0));compressed={}
 for i,c in curve.items():
  for j,a in proj[i].items():compressed[j]=compressed.get(j,Q(0))+c*a
 compressed={j:c for j,c in compressed.items()if c};assert value==sum((c*compressed.get(i,Q(0))for i,c in q27.items()),Q(0))
 positive=bool(value);newdet=Q(X['fixed_minor_determinant'])*value
 data={'all_eight_source_pairing_terms':terms,'actual_residual_scalar':value,'same_actual_curve_compressed_column27_sparse':list(compressed.items()),'compressed_and_full495_pairings_equal':True,
  'old_fixed26_minor_determinant':Q(X['fixed_minor_determinant']),'new_27minor_determinant_rows_old26_then25_columns_old26_then_curve':newdet,
  'old2019_source_factorization_sha256':X['all2019_fixed_factorization_sha256'],'old2019_source_projection_sha256':X['all2019_source_and_projection_stream_sha256'],
  'actual_initial_membership_source_report_sha256':sha(ROOT/gp),'known_plus_actual_curve_cubic_rank':27 if positive else 26,'positive_length9_and_long_exponent2':positive,
  'positive_proof':'Schur complement ofold26minor is literalresidualscalar. Nonzero makes27sourcesindependent;allknown2019spanwasexactly26 andnewcurveisa genuineinitial byacceptedcolumnidentity. Thusm³(Cprime/ell³)=0,lengthupper9=knownlower9.',
  'zero_is_only_partial_negative':True,'no_new_physical_source_gauge_or_jet':True}
 peak=number_count([r,peer,cp,e,X,dual,curve,q27,proj,terms,compressed,data]);assert peak<=8000,peak
 for v in [value,newdet]+[p[-1]for p in terms]:assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=2048
 blob=json.dumps(enc(data),separators=(',',':'))+'\n';assert sum((ROOT/p).stat().st_size for p in pins)+len(blob.encode())<32<<20;ep=out/'actual-curve-scalar-and-complete-positive-rank-certificate.json';ep.write_text(blob)
 report={'id':STEM,'all_passed':True,'status':'PROVED_single_actual_curve_initial_pairing','actual_residual_scalar':str(value),'positive_Cprime_ell3_length9_and_long_exponent2':positive,'no_negative_exponent_inference_if_zero':True,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'numeric_slots':peak,'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 rp=out/'report.json';rp.write_text(json.dumps(report,indent=2)+'\n');assert sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in out.iterdir())<32<<20;signal.alarm(0);print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source-peer-report',required=True);p.add_argument('--closure-peer-report',required=True);a=p.parse_args();main(a.source_peer_report,a.closure_peer_report)
