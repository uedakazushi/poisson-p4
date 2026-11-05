"""Synthetic Q5 codec/sign/zero and interrupted-segment fixtures; no source query."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
from s61_col53_next17_79_streamed_value_20261105_v64 import recover_committed
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_STREAM_CODEC_FIXTURES_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
def main():
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 cases=[(Q(-7,30),{(-3,0,6,-1,10):12,(-3,0,7,-1,9):-18}),
        (Q(11,37),{(2,-2,0,4,1):-25}),
        (Q(-13,21),{}),(Q(0),{})]
 receipts=[]
 for k,(qi,vals)in enumerate(cases):
  p=out/('fixture-'+str(k)+'.bin');r=write_cell(p,qi,vals);assert validate_cell(p,qi,vals,r)
  body=p.read_bytes()
  if vals:
   d=decode(body[9:]);stored=Q(d['content']);ratio=stored/qi;assert ratio.denominator==1 and ratio.numerator>0
   recovered={tuple(d['shift'][j]+int(t[j+1])for j in range(5)):int(t[0])*ratio.numerator for t in d['terms']}
   assert recovered==vals and(stored>0)==(qi>0)
   if k==0:assert ratio==6 and d['shift']==[-3,0,6,-1,9]
   mathG=ratio.numerator
  else:
   assert body==b'Q5CELL1\n\0'and r['content']=='0'and r['G']==1;mathG=0
   # This branch deliberately performs no qi division, including qi=0.
  receipts.append({'fixture':k,'fixed_q':str(qi),'math_G':mathG,'codec':r})
 # Two literal generations, then an interrupted trailing manifest and payload.
 streamed=out/'stream-cursor';streamed.mkdir()
 payloads=[(out/'fixture-0.bin').read_bytes(),(out/'fixture-2.bin').read_bytes()]
 lines=[];at=0
 for k,b in enumerate(payloads):
  lines.append(enc({'node':k,'stream_offset':at,'bytes':len(b),'payload_sha256':hashlib.sha256(b).hexdigest()}));at+=len(b)
 (streamed/'node-values.q5stream').write_bytes(b''.join(payloads)+b'preserved-incomplete-tail')
 (streamed/'node-values-index.jsonl').write_bytes(b''.join(lines)+b'{"node":2')
 cursor=recover_committed(streamed);assert cursor=={'nodes_committed':2,'last_node':1,'stream_prefix_bytes':at,'manifest_prefix_bytes':sum(map(len,lines))}
 assert (streamed/'node-values.q5stream').read_bytes()[:at]==b''.join(payloads)
 # A complete but hash-invalid third segment is also excluded from the cursor.
 bad=out/'invalid-tail-cursor';bad.mkdir()
 (bad/'node-values.q5stream').write_bytes(b''.join(payloads)+b'x')
 (bad/'node-values-index.jsonl').write_bytes(b''.join(lines)+enc({'node':2,'stream_offset':at,'bytes':1,'payload_sha256':'0'*64}))
 assert recover_committed(bad)==cursor
 pins={f:sha(ROOT/f)for f in['code/s61_col53_next17_79_stream_codec_fixtures_20261105_v64.py','code/s61_col53_next17_79_streamed_value_20261105_v64.py','code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py','code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py']}
 (out/'synthetic-fixtures.json').write_bytes(enc({'cases':receipts,'committed_cursor_after_both_failure_shapes':cursor,'two_generations_byteidentical':True,'actual_shift_may_differ_from_fixed_SOURCE_lead':True}))
 d={'id':ID,'all_passed':True,'synthetic_only_no_source_value_queries':True,'negative_q_content_restoration':True,'actual_zero_q0_and_nonzero_q_mathG0_codecG1':True,'two_generations_offsets_hashes_and_interrupted_cursor':True,'input_sha256':pins,'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()}}
 (out/'report.json').write_bytes(enc(d));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':main()
