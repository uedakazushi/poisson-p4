#!/usr/bin/env python3
"""SOURCE I/O gate only: same primitive bytes, heapsort aliases and stream readback."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,json,hashlib,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import encode,decode
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell,descending_alias_keys
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STREAM_CODEC_SOURCE_GATE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def main(a):
 assert a.check
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);cp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64/report.json';assert sha(cp)=='f7f9bc898a54b9a457a59ba4146f38b0e5a8ce9f2cbda906379a71db20e5f9cb';r=read(cp);assert r['all_passed'];ap=cp.parent/'actual-cleared-col82-and-commonQ3RE-representation.json';assert sha(ap)==r['evidence_sha256'][ap.name];api=read(ap);inputs={cp:sha(cp),ap:sha(ap),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')}
 for name in ['s61_sparse_generic5_primitive_binary_codec_20261105_v64.py','s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py']:p=ROOT/'code'/name;inputs[p]=sha(p)
 records=[];peak=25000
 for v in api['actual16_col82_numerator_values']:
  p=cp.parent/v['file'];assert sha(p)==v['sha256'];inputs[p]=sha(p);blob=p.read_bytes();assert blob[8];old=decode(blob[9:]);q=Q(old['content']);D={tuple(old['shift'][j]+t[j+1]for j in range(5)):int(t[0])for t in old['terms']};peak=max(peak,25000+32*len(D)+1024);assert peak<150000
  keys=descending_alias_keys(D);assert keys==sorted(D,key=lambda e:(sum(e),e),reverse=True)and all(any(k is x for x in D)for k in keys);del keys
  f=out/f"existing-row-{v['row']}.bin";result=write_cell(f,q,D,old['shift']);assert validate_cell(f,q,D,result)and f.read_bytes()==blob and result['value_sha256']==v['actual_QQ_numerator_value_sha256'];records.append({'row':v['row'],'terms':len(D),'byteidentical_to_original':True,'heap_order_and_key_identity':True,'incremental_all_coefficient_readback':True,'sha256':sha(f)})
 fixtures=[('zero',Q(0),{},[0]*5),('negative_large',Q(-(2**293+1),65537),{(-13,2,0,0,4):1,(0,0,0,0,0):-(2**1651+1),(4,0,2,0,-3):17,(4,0,1,1,-3):-19},[-13,0,0,0,-3])]
 for name,q,D,shift in fixtures:
  f=out/f'fixture-{name}.bin';rec=write_cell(f,q,D,shift);assert validate_cell(f,q,D,rec)
  if D:
   G=0
   for c in D.values():G=gcd(G,abs(c))
   expected={'content':str(q*G),'shift':shift,'terms':[[str(D[e]//G),*[e[j]-shift[j]for j in range(5)]]for e in sorted(D,key=lambda e:(sum(e),e),reverse=True)]};wire=b'Q5CELL1\n\x01'+encode(expected)
  else:wire=b'Q5CELL1\n\x00'
  assert f.read_bytes()==wire;records.append({'fixture':name,'terms':len(D),'legacy_wire_equal':True,'stream_readback':True,'sha256':sha(f)})
 (out/'existing4160-and-fixture-roundtrips.json').write_bytes(enc({'all_passed':True,'records':records,'existing_coefficients':4160,'source_test_only_legacy_arrays_allowed':True,'future_VALUE_no_legacy_encode_decode_arrays':True,'alias_keys_heap_storage':'n references only; custom in-place heapsort has O(1) auxiliary integer scratch; all exponent tuples alias map keys','future_map_plus_heap_plus_scratch_bound':'7n+n+O(1) <=10n+2048'}));base=sum(p.stat().st_size for p in inputs);assert base+sum(p.stat().st_size for p in out.iterdir())+65536<32*1024*1024;report={'id':ID,'all_passed':True,'source_only':True,'no_new_matrix_coefficients':True,'all4160_existing_coefficients_wire_identical':True,'all16_existing_and2_fixtures_stream_readback':True,'SOURCE_owned_upper':peak,'actual_READ_RAW_bytes':base,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');main(p.parse_args())
