"""Independent stream decoder/readback of saved faces only; no coefficient replay."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_H_FACES_VALUES_READBACK_20261105_v64'
PRIMARY='S61_COL53_H_FACES_VALUES_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def decode(path,H,record,allowed):
 qterms=0;primitive_gcd=0;valuehash=hashlib.sha256();valuehash.update((str(record['terms'])+'\n').encode());witness=None
 with path.open('rb')as f:
  assert f.read(8)==b'Q5CELL1\n';flag=f.read(1)
  def u():
   value=0;shift=0
   while True:
    b=f.read(1);assert b;byte=b[0];value|=(byte&127)<<shift
    if byte<128:return value
    shift+=7;assert shift<=2048
  def i():
   tag=u();b=f.read(tag//2);assert len(b)==tag//2 and(not b or b[0])
   a=int.from_bytes(b,'big');assert a or not tag%2
   return -a if tag%2 else a
  if not record['terms']:
   assert flag==b'\0'and not f.read(1);content=Q(0);shift=[0]*5
  else:
   assert flag==b'\1';num=i();den=i();assert den>0
   content=Q(num,den);assert content and(content.numerator,content.denominator)==(num,den)
   shift=[]
   for j in range(5):a=u();shift.append(a//2 if not a%2 else-(a//2)-1)
   assert shift==record['shift'];n=u();assert n==record['terms'];previous=None
   for k in range(n):
    a=i();assert a
    e=tuple(shift[j]+u()for j in range(5));assert e[0]==H and e in allowed
    order=(sum(e),e)
    if previous is not None:assert previous>order
    previous=order;primitive_gcd=gcd(primitive_gcd,abs(a))
    c=content*a;assert c
    valuehash.update((str(c)+' '+' '.join(map(str,e))+'\n').encode());qterms+=1
    if witness is None or e<tuple(witness['exponent5']):witness={'exponent5':list(e),'coefficient_QQ':str(c),'primitive_integer_coefficient':str(a)}
   assert primitive_gcd==1 and not f.read(1)
  assert valuehash.hexdigest()==record['value_sha256']
  assert str(content)==record['content']and sha(path)==record['sha256']
 return {'H':H,'actual_terms':qterms,'primitive_gcd':primitive_gcd,'content_QQ':str(content),'witness':witness,'value_sha256':valuehash.hexdigest()}
if __name__=='__main__':
 b=ROOT/'reports'/PRIMARY;p=b/'report.json';assert sha(p)=='b8d4309eaddecdd9b362934f7d2f8bde1973487d071ca507b730f426db8fd983'
 r=json.loads(p.read_text());assert r['all_passed']
 f=ROOT/'reports/S61_COL53_H_FACES_SOURCE_20261105_v64/two-complete-positive-H-faces.json';F=json.loads(f.read_text());pins={str(p.relative_to(ROOT)):sha(p),str(f.relative_to(ROOT)):sha(f),'code/s61_col53_H_faces_readback_20261105_v64.py':sha(__file__)}
 ans=[]
 for j,H in enumerate([-3,21]):
  allowed={tuple(e)for e in F['faces'][j]['full5_exponents']};rec=r['faces'][j];path=b/rec['file'];pins[str(path.relative_to(ROOT))]=sha(path)
  ans.append(decode(path,H,rec,allowed))
 assert ans[0]['actual_terms']==411 and ans[1]['actual_terms']==0
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 d={'id':ID,'all_passed':True,'scope':'saved411+0 terms stream readback only; not an independent source coefficient replay','independent_binary_decoder_no_producer_codec_import':True,'all5coords_SOURCE_membership_and_DEGLEX_unique':True,'faces':ans,'input_sha256':pins,'entry_nonzero_primary_value_proof_pending_separate_source_replay_peer':True}
 (out/'report.json').write_text(json.dumps(d,separators=(',',':'))+'\n')
 print(json.dumps(ans,ensure_ascii=False));print('report_sha256',sha(out/'report.json'))
