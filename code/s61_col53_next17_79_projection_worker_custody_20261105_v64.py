"""SOURCE-only immutable worker custody and synthetic packed-projection fixtures."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
from s61_sparse_generic5_source_packed5_profile_v2_20261105_v64 import pack,unpack,Engine
from s61_lex_HL_S2_source_projection_20261105_v64 import inside
ID='S61_COL53_NEXT17_79_PROJECTION_WORKER_CUSTODY_20261105_v64'
G='reports/S61_COL53_NEXT17_79_SOURCE_PROJECTION_GATE_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
if __name__=='__main__':
 p=ROOT/G/'report.json';assert sha(p)=='c51081147c86f7aaafbe51c50b51620abe5259ac71cf54cede2fc767197b4189';r=json.loads(p.read_text());assert r['all_passed']
 pins={str(p.relative_to(ROOT)):sha(p)}
 for f,h in r['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
 for f,h in r['evidence_sha256'].items():assert sha(p.parent/f)==h;pins[str((p.parent/f).relative_to(ROOT))]=h
 tests=0
 for a in[(-3,0,4),(0,-2,-4)]:
  for b in[(-2,0,3),(2,-3,1)]:
   A={(a[0],a[1],a[2]+i,-4+i,10-i)for i in range(5)}
   B={(b[0],b[1],b[2]+i,2-i,-8+i)for i in range(5)}
   AP={pack(e)for e in A if inside(e,a)};BP={pack(e)for e in B if inside(e,b)}
   assert all(pack(unpack(k))==k for k in AP|BP)
   lead=tuple(a[j]+b[j]for j in range(3))
   for cap in[0,1,2]:
    got=set()
    for x in AP:
     e=unpack(x)
     for y in BP:
      f=unpack(y)
      if(e[2]-a[2])+(f[2]-b[2])<=cap:got.add(pack(tuple(e[j]+f[j]for j in range(5))))
    expected={tuple(e[j]+f[j]for j in range(5))for e in A for f in B if inside(e,a)and inside(f,b)and(e[2]-a[2])+(f[2]-b[2])<=cap}
    assert{unpack(k)for k in got}==expected
    assert[unpack(k)for k in Engine._lex_alias(got)]==sorted(expected)
    tests+=1
 # Out-of-range coordinate must fail before packing; never unchecked key sums.
 for e in[(-65,0,0,0,0),(0,192,0,0,0),(0,0,0,-128,0)]:
  try:pack(e)
  except AssertionError:tests+=1
  else:raise AssertionError('packed source overflow accepted')
 for f in[Path(__file__),ROOT/'code/s61_lex_HL_S2_source_projection_20261105_v64.py']:pins[str(f.relative_to(ROOT))]=sha(f)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 d={'id':ID,'all_passed':True,'SOURCE_implementation_and_byte_custody_ONLY':True,'profile_and_signed_values_NOT_run':True,'synthetic_formal_exponent_fixtures':tests,'all5coords_decode_add_check_pack':True,'canonical_tuplelex_hash_alias_verified':True,'VALUE_full5_plus_coeff7_NOT_changed':True,'checkpoint_stream_has_no_list_or_sort_alias':True,'input_sha256':pins}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'));print('fixtures',tests)
