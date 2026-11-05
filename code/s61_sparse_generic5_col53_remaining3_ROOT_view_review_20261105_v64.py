"""Independent SOURCE-only review of immutable views and fused support union."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,importlib.util,json,sys,tempfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from s61_sparse_generic5_source_view_profile_20261105_v64 import Engine,unpack
from s61_sparse_generic5_source_packed5_profile_v2_20261105_v64 import pack,count_numeric

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_ROOT_REVIEW_20261105_v64';out.mkdir(exist_ok=False)
 gate=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_SOURCE_GATE_20261105_v64/report.json';g=json.loads(gate.read_text());assert g['all_passed'] and g['profile_not_executed']
 for name,h in g['input_sha256'].items():assert sha(ROOT/name)==h,name
 for name,h in g['evidence_sha256'].items():assert sha(gate.parent/name)==h,name
 bp=gate.parent/'resume-views-fused-budget.json';b=json.loads(bp.read_text());assert b['remaining_components']==[8,9,10] and b['remaining_terms']==93 and b['same_SOURCE_owned_limit']==150000
 assert b['generated_RAW_upper']==8047824 and b['fresh_RAW_upper']==26158253<33554432
 assert b['actual_READ_files_RAW']+b['gate_receipt_RAW_reserve']+b['ROOT_receipt_RAW_reserve']+b['generated_RAW_upper']==b['fresh_RAW_upper']
 assert b['stored_M_bit_guard']==512 and b['positive_M_proven_upper_bits']==404 and b['stored_QQ_num_den_bit_guard']==128 and b['QQ_proven_upper_bits']==80
 assert b['failure_numeric_reserve_ALWAYS_charged']==4096
 # Tiny fixtures only: no original remaining support or coefficient values.
 e=Engine(metadata_base=10000);P={(-1,2,0,1,0),(0,0,1,-1,2),(2,-2,0,0,-1)};T={(0,1,0,0,0),(1,0,-1,0,1)}
 e.save('P',{pack(x)for x in P},Q(2,3),7,{'op':'ROOT_fixture'})
 e.save('T',{pack(x)for x in T},Q(-3,5),11,{'op':'ROOT_fixture'})
 e.copy('V','P',Q(-5,7),(-2,1,2,0,-1));assert e.live['V']['S']is e.live['P']['S']
 V={tuple(x[j]+(-2,1,2,0,-1)[j] for j in range(5))for x in P};assert {unpack(v+e.live['V']['off'])for v in e.live['V']['S']}==V
 e.mul('W','V','T',Q(4,9),(1,-1,0,2,0));W={tuple(x[j]+y[j]+(1,-1,0,2,0)[j]for j in range(5))for x in V for y in T}
 assert {unpack(v)for v in e.live['W']['S']}==W
 C={(0,0,0,0,0),(1,1,1,1,1)};e.save('C',{pack(x)for x in C},Q(5,11),13,{'op':'ROOT_fixture'})
 e.fused_union_product('SUM','C','V','T','DOT');expected=C|{tuple(x[j]+y[j]for j in range(5))for x in V for y in T}
 assert {unpack(v)for v in e.live['SUM']['S']}==expected
 assert 'DOT'not in e.live and e.nodes['DOT']['terms_upper']==len(P)*len(T)
 assert 'SAFE UPPERS' in e.nodes['DOT']['support_profile']
 assert all(e.nodes['SUM']['H_slice_sizes']==[[h,sum(v[0]==h for v in expected)]for h in sorted({v[0]for v in expected})]for _ in [0])
 # Independently derive rational rescaling weights and positive M bound.
 dq=Q(-10,21)*Q(-3,5);aa=Q(5,11);den=aa.denominator*dq.denominator//__import__('math').gcd(aa.denominator,dq.denominator)
 a=int(aa*den);z=int(dq*den);gg=__import__('math').gcd(abs(a),abs(z));M=abs(a//gg)*13+abs(z//gg)*min(len(P),len(T))*7*11
 assert e.live['SUM']['q']==Q(gg,den) and e.live['SUM']['M']==M
 # Failure snapshot stores each backing once plus offsets, exact byte reader.
 shot=e.failure_view_snapshot(out/'fixture-views.q5view')
 with (out/'fixture-views.q5view').open('rb')as f:
  assert f.read(8)==b'Q5VIEW1\n';hdr=json.loads(f.read(int.from_bytes(f.read(4),'big')));blocks=[]
  for n in hdr['blocks']:blocks.append({int.from_bytes(f.read(5),'big')for _ in range(n)})
  assert not f.read(1)
 for r in hdr['maps']:assert {unpack(v+r['o'])for v in blocks[r['b']]}=={unpack(v+e.live[r['id']]['off'])for v in e.live[r['id']]['S']}
 assert next(x for x in hdr['maps']if x['id']=='P')['b']==next(x for x in hdr['maps']if x['id']=='V')['b']
 # Aliased accumulator mutation is rejected before union arithmetic.
 e.copy('ALIAS','SUM')
 try:e.fused_union_product('BAD','SUM','P','T','BAD_DOT')
 except AssertionError:pass
 else:raise AssertionError('shared accumulator accepted')
 e.drop('ALIAS','P','V','T','W');h=e.boundary_hash('SUM')
 canonical=(json.dumps([list(x)for x in sorted(expected)],separators=(',',':'))+'\n').encode();assert h==hashlib.sha256(canonical).hexdigest()
 e.snapshot(out/'fixture-boundary.q5supports')
 worker=ROOT/b['worker'];helper=ROOT/b['helper'];c=worker.read_text();t=helper.read_text()
 assert "for component in[8,9,10]" in c and "counter=188" in c and "eng.fused_union_product" in c and "term-{component}-{ti:02d}" in c
 assert "len(termsaved)==93" in c and "rootdoc['all_passed']and count_numeric(rootdoc)<=2000" in c
 assert "self.owners(c['S'])==1" in t and "self.live.pop(C)" in t and "self.own(1+2048);S.add(v)" in t
 assert "no_dot_support_allocated" in t and "no_set_copy" in t
 report=dict(all_passed=True,source_only=True,no_actual_remaining_support_profile_or_coefficient_VALUE=True,
 all_gate_input_evidence_hashes_current=True,independent_signed5coordinate_views_multiply_fused_union_qM_exact=True,
 unique_backing_storage_and_alias_mutation_rejection=True,lossless_shared_view_snapshot_independently_decoded=True,
 boundary_decoded_tuple_lex_hash_independently_checked=True,old_phase07_only_remaining93_source_bound=True,
 whole_ownership_same150000=True,failure_reserve4096=True,fresh_RAW_upper=b['fresh_RAW_upper'],
 future_scope='ONE SOURCE profile from live07 only; no VALUE; unmaterialized dot cardinality/Hcounts remain upper; endpoint union and snapshot exact. All167x14/87+66/oneP290/F9F10/complements preserved.',
 input_sha256={str(gate.relative_to(ROOT)):sha(gate),str(bp.relative_to(ROOT)):sha(bp),str(worker.relative_to(ROOT)):sha(worker),str(helper.relative_to(ROOT)):sha(helper)},code_sha256=sha(Path(__file__)))
 assert count_numeric(report)<=2000
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert (out/'report.json').stat().st_size<=65536;print(sha(out/'report.json'))
if __name__=='__main__':main()
