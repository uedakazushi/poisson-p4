"""Prepare SOURCE-only continuation from LIVE08, never replay completed45 terms."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from s61_sparse_generic5_source_view_profile_20261105_v64 import Engine,count_numeric,unpack

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 o=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_SOURCE_GATE_20261105_v64';o.mkdir(exist_ok=False)
 old=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_SOURCE_GATE_20261105_v64';g=json.loads((old/'report.json').read_text());assert g['all_passed'];pins=g['input_sha256'].copy()
 for n,h in pins.items():assert sha(ROOT/n)==h
 prior=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_PROFILE_20261105_v64';p8=json.loads((prior/'phase-08.json').read_text());snap=prior/'live-08.q5supports';assert sha(snap)==p8['all_live_snapshot']['sha256']
 with snap.open('rb')as f:
  assert f.read(8)==b'Q5SUPP1\n';head=json.loads(f.read(int.from_bytes(f.read(4),'big')));wanted={x['id']for x in head['maps']};nodes={}
  for rec in head['maps']:
   S=set();lo=[191]*5;hi=[-64]*5;hs={}
   for _ in range(rec['terms']):
    raw=f.read(5);assert len(raw)==5;v=int.from_bytes(raw,'big');assert v not in S;S.add(v);ex=unpack(v)
    for j,x in enumerate(ex):lo[j]=min(lo[j],x);hi[j]=max(hi[j],x)
    hs[ex[0]]=hs.get(ex[0],0)+1
   assert 50000+len(S)+8192<150000
   nodes[rec['id']]=dict(terms_upper=len(S),known_QQ_content=rec['q'],primitive_integer_abs_upper=rec['M'],exponent_box5=[[lo[j],hi[j]]for j in range(5)]if S else[[0,0]]*5,H_slice_sizes=[[h,n]for h,n in sorted(hs.items())],support_sha256=None)
  assert not f.read(1)
 assert len(nodes)==36 and sum(n['terms_upper']for n in nodes.values())==45842
 phase7=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64/phase-07.json';proven=json.loads(phase7.read_text())['nodes']
 for p in sorted(prior.glob('term-8-*.json')):proven.update(json.loads(p.read_text())['nodes'])
 assert len(list(prior.glob('term-8-*.json')))==45
 for k,v in nodes.items():
  for field in ['terms_upper','known_QQ_content','primitive_integer_abs_upper','exponent_box5','H_slice_sizes']:assert v[field]==proven[k][field],(k,field)
 assert p8['accumulator']=='new_593_sum';nodes['new_593_sum']['support_sha256']=p8['boundary_hash']
 state=dict(accumulator='new_593_sum',canonical_accumulator_hash=p8['boundary_hash'],nodes=nodes,snapshot_path=str(snap.relative_to(ROOT)),snapshot_sha256=sha(snap))
 (o/'resume-live08-metadata.json').write_text(json.dumps(state,separators=(',',':'))+'\n')
 # Two-component ID lifecycle fixture, using frozen helper and same retirement.
 e=Engine(metadata_base=10000);e.leaf('F',{'content':'1','shift':[0]*5,'terms':[['1',0,0,0,0,0]]});acc=None
 for component in [9,10]:
  rawname='raw53_component'+str(component);e.leaf(rawname,{'content':'1','shift':[0]*5,'terms':[['1',1,0,0,0,0]]});out='sum'+str(component);e.fused_union_product(out,acc,'F',rawname,'dot'+str(component));acc=out;e.drop(rawname);e.nodes={k:v for k,v in e.nodes.items()if k in e.live};e.ops=[];e.meta=sum(count_numeric(v)+12 for v in e.nodes.values());e.own();assert rawname not in e.nodes and rawname not in e.live
 e.boundary_hash(acc);e.snapshot(o/'two-component-fixture.q5supports')
 worker=ROOT/'code/s61_sparse_generic5_constrained83_col53_last2_view_profile_20261105_v64.py';helper=ROOT/'code/s61_sparse_generic5_source_view_profile_20261105_v64.py';code=worker.read_text();assert 'for component in[9,10]'in code and "rawname='raw53_component'+str(component)"in code and 'eng.drop(rawname)'in code and 'counter=593'in code
 for p in [old/'report.json',old/'resume-views-fused-budget.json',phase7,worker,helper,Path(__file__)]:pins[str(p.relative_to(ROOT))]=sha(p)
 for p in prior.iterdir():
  if p.is_file():pins[str(p.relative_to(ROOT))]=sha(p)
 b=json.loads((old/'resume-views-fused-budget.json').read_text());b.update(id=o.name,worker=str(worker.relative_to(ROOT)),worker_sha256=sha(worker),helper_sha256=sha(helper),resume_phase07=False,resume_phase08=True,resume_snapshot=str(snap.relative_to(ROOT)),resume_snapshot_sha256=sha(snap),remaining_components=[9,10],remaining_terms=48,old45_not_replayed=True)
 # Includes48 receipts, two snapshots/phases, failure snapshot+receipt/report.
 generated=48*49152+2*(665204+65536)+694644+262144+262144+65536
 readraw=sum((ROOT/n).stat().st_size for n in pins);b.update(generated_RAW_upper=generated,actual_READ_files_RAW=readraw,gate_receipt_RAW_reserve=200000,ROOT_receipt_RAW_reserve=65536,fresh_RAW_upper=readraw+generated+200000+65536)
 assert b['fresh_RAW_upper']<33554432
 (o/'resume-views-fused-budget.json').write_text(json.dumps(b,separators=(',',':'))+'\n')
 with(o/'help.log').open('w')as f:subprocess.run([sys.executable,str(worker),'--help'],check=True,stdout=f,stderr=subprocess.STDOUT)
 report=dict(all_passed=True,source_only=True,profile_NOT_run=True,all36_live_profiles_recovered_from_source_receipts=True,all45842_snapshot_keys_and_boxes_verified=True,two_component_ID_lifecycle_fixture_PASS=True,old8_and45_never_replayed=True,remaining_terms=48,fresh_RAW_upper=b['fresh_RAW_upper'],same150k_55math60hard_2GB=True,input_sha256=pins,code_sha256=sha(Path(__file__)),evidence_sha256={p.name:sha(p)for p in o.iterdir()});(o/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(sha(o/'report.json'))
if __name__=='__main__':main()
