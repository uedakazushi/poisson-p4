#!/usr/bin/env python3
"""One bounded exact cancellation of the already verified seven source maps."""
import argparse, gzip, hashlib, json, resource, signal, time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_LOW_INVARIANT_POLE_CANCELLATION_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s cap')));signal.alarm(60)
 out.mkdir(parents=True,exist_ok=False)
 a,Y,W,B=s.symbols('a Y W B');H=s.Poly(147456*Y-(3072*W-1)**2,a,Y,W,B,domain=s.QQ)
 base=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63'
 review=ROOT/'reports/Q_NEW_FAMILY_LOW_INVARIANT_EXPANSION_20261104_v63_REVIEW/report.json'
 assert json.loads(review.read_text())['all_passed']
 results=[];inputs=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),review]
 for name in ['F4','p0','p1','p2','p3','p4','N12']:
  path=base/(name+'.json.gz');inputs.append(path);d=json.loads(gzip.decompress(path.read_bytes()))
  p=s.Poly.from_dict({tuple(row[1:]):s.Integer(row[0])for row in d['terms']},(a,Y,W,B),domain=s.QQ)
  original=p;k=0;cap=d['denominator_D_power']
  while k<cap:
   q,r=p.div(H)
   if not r.is_zero:break
   p=q;k+=1
  assert p*H**k==original
  raw=dict(name=name,original_terms=len(d['terms']),cancelled_H_power=k,
   denominator_H_power=cap-k,denominator_B_power=d['denominator_B_power'],
   denominator_integer=d['denominator_integer'],source_c_power=d['source_c_power'],
   variables=['a=u/X','Y=X^2','W','B'],H='147456*Y-(3072*W-1)^2',
   terms=[[str(c.p),str(c.q),*ex]for ex,c in p.terms()])
  dest=out/(name+'.json.gz');dest.write_bytes(gzip.compress(json.dumps(raw,separators=(',',':')).encode(),mtime=0))
  results.append({k:v for k,v in raw.items()if k!='terms'}|{'reduced_terms':len(raw['terms']),'exact_remultiplication_verified':True})
  (out/'frontier.json').write_text(json.dumps(results,indent=2)+'\n')
 report=dict(id=ID,all_passed=True,results=results,no_new_points_relation_search_or_degree_increase=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({'all_passed':True,'seconds':report['seconds'],'results':results},indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'reports'/ID);ns=ap.parse_args()
 try:main(ns.output.resolve())
 except (TimeoutError,MemoryError)as e:
  ns.output.mkdir(parents=True,exist_ok=True);(ns.output/'BLOCKED.json').write_text(json.dumps({'complete':False,'constraint':str(e),'no_auto_retry':True})+'\n');raise
