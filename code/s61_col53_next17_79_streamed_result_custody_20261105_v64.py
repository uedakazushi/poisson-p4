"""Read-only serialized result custody; not an independent signed source replay."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_STREAMED_RESULT_CUSTODY_20261105_v64'
PRIMARY='S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);src=ROOT/'reports'/PRIMARY;rp=src/'report.json'
 assert sha(rp)=='4b1ed0275771c6f13a8f250b14a71a125277d3b8e6f17fd5995142b8d4dc0c54'
 r=json.loads(rp.read_text());assert r['all_passed']and r['nodes_committed']==1884
 for f,h in r['evidence_sha256'].items():assert sha(src/f)==h
 count=0;at=0;mathzero=0;root=None;tmp=ROOT/'work'/(PRIMARY+'-canonical-cells')
 with(src/'node-values-index.jsonl').open('rb')as m,(src/'node-values.q5stream').open('rb')as s:
  for line in m:
   d=json.loads(line);assert d['stream_offset']==at;body=s.read(d['bytes'])
   assert len(body)==d['bytes']and hashlib.sha256(body).hexdigest()==d['payload_sha256']
   literal=tmp/('node-'+str(d['node'])+'.bin');assert literal.read_bytes()==body
   if not d['actual_nonzero_slots']:
    assert body==b'Q5CELL1\n\0'and d['math_G']==0 and d['codec']['G']==1;mathzero+=1
   if d['node']==2997:root=d
   count+=1;at+=len(body)
  assert not s.read(1)
 assert count==1884 and at==859628 and root['actual_nonzero_slots']==0
 files=[str(rp.relative_to(ROOT)),str((src/'node-values-index.jsonl').relative_to(ROOT)),str((src/'node-values.q5stream').relative_to(ROOT)),str(Path(__file__).relative_to(ROOT)),'research/notes/'+PRIMARY+'.md']
 d={'id':ID,'all_passed':True,'serialization_only_not_independent_signed_source_peer':True,'all1884_two_generation_payloads_byteidentical':True,'complete_cells':count,'empty_cells_mathG0_codecG1':mathzero,'stream_bytes':at,'fixed_needed_projection_zero_only':True,'full_minor_zero_NOT_claimed':True,'input_sha256':{f:sha(ROOT/f)for f in files}}
 (out/'report.json').write_text(json.dumps(d,separators=(',',':'))+'\n');print('report_sha256',sha(out/'report.json'));print('empty_cells',mathzero)
if __name__=='__main__':main()
