#!/usr/bin/env python3
"""Parse bounded M12 binary quotient/guard witnesses; never resume arithmetic."""
import argparse,json,struct
from pathlib import Path
def parse(path):
 if path.suffix=='.json':
  from fractions import Fraction
  x=json.loads(path.read_text());assert x['version']=='M12SP001';count=0
  for p in list(x['atoms'].values())+list(x['completed_source_values'].values())+x['extras']+[x['current_uncommitted_polynomial']]+x['integer_matrix_prefix']:
   last=None
   for d,c in p:
    q=Fraction(c);assert -32<=d<=216 and(last is None or last<d)and abs(q.numerator).bit_length()<=2048 and q.denominator.bit_length()<=2048;last=d;count+=1
  return dict(version=x['version'],source_phase=x['phase'],cursor=x['cursor'],reason=x['reason'],coefficient_records=count,parse_only_no_resume=True)
 with path.open('rb')as f:
  def get(fmt):
   n=struct.calcsize('<'+fmt);b=f.read(n);assert len(b)==n;return struct.unpack('<'+fmt,b)[0]
  def string():
   n=get('H');assert n<8192;b=f.read(n);assert len(b)==n;return b.decode()
  def z():
   s=get('B');n=get('H');assert s<2 and n<=256;b=f.read(n);assert len(b)==n and(not n or b[0]);v=int.from_bytes(b,'big');assert v.bit_length()<=2048;return -v if s else v
  def p(cap):
   n=get('H');assert n<=cap+1;last=-1
   for _ in range(n):
    d=get('H');v=z();assert last<d<=cap and v;last=d
   return n
  magic=f.read(8);result={'version':magic.decode(),'parse_only_no_resume':True}
  if magic==b'M12BP001':
   k=get('H');n=get('H');assert 0<=k<12 and n==12;h=string();assert len(h)==64;count=p(108);nu=get('I');assert nu==(11-k)**2
   for _ in range(nu):
    i=get('H');j=get('H');assert k<i<12 and k<j<12;count+=p(108)
   result.update(stage=k,input_sha256=h,coefficient_records=count,quotient_updates=nu)
  elif magic==b'M12AS001':
   h=string();reason=string();phase=string();v=[get('I')for _ in range(9)];v=[x-(1<<32)if x>=1<<31 else x for x in v];c=[get('Q')for _ in range(3)]
   for _ in range(4):z()
   n=get('H');assert n==144;count=sum(p(108)for _ in range(n));b=get('H');assert b==5;count+=sum(p(216)for _ in range(b));result.update(input_sha256=h,reason=reason,phase=phase,cursors=v,counters=c,coefficient_records=count)
  elif magic==b'M12FX001':
   a,b,c=z(),z(),z();assert a==-1 and b==0 and c==(1<<2047)-1;assert p(2)==2;result.update(fixture_negative_zero_2048bit_roundtrip=True)
  else:raise ValueError('wrong M12 binary version')
  assert not f.read(1);return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('path',type=Path);a=p.parse_args();print(json.dumps(parse(a.path),indent=2))
