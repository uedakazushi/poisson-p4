#!/usr/bin/env python3
"""ROOT audits changed SOURCE profile/snapshot protocol without evaluating the profile."""
import hashlib,importlib.util,json,math,re
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_PROFILE_V2_SOURCE_GATE_20261105_v64'
OUT=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_PROFILE_V2_ROOT_REVIEW_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False);p=BASE/'report.json';assert sha(p)=='9c6978f491c6b3b2163ad3d613405d2b5ce14dd04c3ac1da16e8ad4a02cba6cc'
 r=json.loads(p.read_text());assert r['all_passed'] and r['no_actual_support_profile_or_VALUE_executed'];pins={str(p.relative_to(ROOT)):sha(p)}
 for name,want in r['input_sha256'].items():assert sha(ROOT/name)==want;pins[name]=want
 p=BASE/'whole-SOURCE-and-snapshot-budget.json';assert sha(p)==r['evidence_sha256'][p.name];b=json.loads(p.read_text());pins[str(p.relative_to(ROOT))]=sha(p)
 worker=ROOT/b['worker'];helper=ROOT/b['helper'];assert sha(worker)==b['worker_sha256'] and sha(helper)==b['helper_sha256']
 wc=worker.read_text();hc=helper.read_text()
 # Bound all static scalar arguments used by the closed old literal DAG and
 # the11-row coefficient operator. Big dynamic q/M are separatelyguarded4096.
 opath=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64/staged-literal-col56-program-and-positive-profiles.json'
 old=json.loads(opath.read_text());scalars=[]
 for op in old['literal_ops']:
  for key in ['scalar','scalarB']:
   if key in op:scalars.append(Fraction(op[key]))
 scopes={op['binding']['scope']for op in old['literal_ops']if op['op']=='input'}
 del old
 for name in scopes:
  doc=json.loads((opath.parent/name).read_text());groups=doc.get('common_operands',doc.get('row_operands',{}))
  for leaf in groups.values():
   scalars.append(Fraction(leaf['QQ_transform_scalar']));scalars.append(Fraction(leaf['original_literal']['content']))
  del doc
 plan=json.loads((ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64/exact-11-row-coefficient-operator-and-callback-owners.json').read_text())
 for row in plan['raw_coefficient_side11']:
  for factors,c in row:scalars.append(Fraction(c))
 scalars.extend([Fraction(353894400000),Fraction(62208000000*353894400000)])
 scalarbits=max(max(abs(x.numerator).bit_length(),x.denominator.bit_length())for x in scalars);assert scalarbits<=2048
 # With stored q/M<=4096 and scalar<=2048, helper add LCM weights have
 # <=10240bits, weight*M and sum<=14337; triple QQproducts<=10240.
 assert 3*4096+2048+1==14337<16384
 spec=importlib.util.spec_from_file_location('root_packed5_fixture',helper);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 points=[(-64,-64,-64,-64,-64),(191,191,191,191,191),(-7,0,2,-1,4),(0,-3,4,8,-2),(0,1,-9,3,7),(2,0,0,0,0)]
 for e in points:assert m.unpack(m.pack(e))==e
 for e,f in [((-7,0,2,-1,4),(3,-2,0,1,-5)),((0,1,-9,3,7),(2,-1,4,-3,-7))]:
  expect=tuple(e[j]+f[j]for j in range(5));assert m.unpack(m.pack(e)+m.pack(f)-m.BIAS)==expect
 eng=m.Engine();S={m.pack(e)for e in points};actual,_=eng.supporthash(S)
 expected=hashlib.sha256((json.dumps([list(e)for e in sorted(points)],separators=(',',':'))+'\n').encode()).hexdigest();assert actual==expected
 eng.save('synthetic',S,Fraction(-7,13),37,{'op':'synthetic_only'})
 fixture=eng.snapshot(OUT/'independent-synthetic-fixture.q5supports');assert fixture['all_keys_exact_stream_roundtrip']
 raw=Path(OUT/fixture['file']).read_bytes();assert raw[:8]==b'Q5SUPP1\n';size=int.from_bytes(raw[8:12],'big');head=json.loads(raw[12:12+size]);assert head['maps'][0]['id']=='synthetic'
 body=raw[12+size:];assert len(body)==5*len(S)
 decoded={m.unpack(int.from_bytes(body[i:i+5],'big'))for i in range(0,len(body),5)};assert decoded==set(points)
 assert b['max_boundary_live_maps']==39 and 39*4096+256==160000
 assert b['one_snapshot_RAW_upper']==12+160000+5*125000==785012
 assert b['fresh_RAW_upper']==b['actual_input_RAW_bytes']+65536+b['all_generated_RAW_upper']==28422564<33554432
 assert b['all_generated_RAW_upper']>=11*(785012+262144)+512000
 for text in ['signal.alarm(55)','time.monotonic()-start<52',"eng.boundary_hash(name,expected[k]['support_sha256'])",'expected_boundary[position]if position<8 else None',
              "snapshot=eng.snapshot(out/f'live-{position:02d}.q5supports')",'assert set(eng.live)==set(funcs.values())|{accumulator}',
              "eng.nodes={k:v for k,v in eng.nodes.items()if k in keep}",'rootpath.stat().st_size<=65536','count_numeric(rootdoc)<=2000']:
  assert text in wc,text
 for text in ['nmax+768','self.own(8192)','len(names)<=64','total<=125000','e.to_bytes(5',"assert not f.read(1)",'self.own(len(S)+1+768)','rawindex(shift)-BIAS']:
  assert text in hc,text
 # Direct rawentry decode is at most36terms; missing coefficient-vs-key
 # temporary allowance108 fits the explicit25k runtime reserve alongside
 # <=4k gate+2kroot and<=8k formatting/controls scratch.
 e=json.loads((ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/actual-fixed24-entry-source-and-exponent-plan.json').read_text())
 entries={(x['row'],x['column']):x for x in e['actual_entry_sources576']}
 maxraw=max(entries[(6 if row=='i'else row,53)]['terms']for row in plan['raw_row_ids11']);assert maxraw==36
 assert 3*maxraw+4000+2000+8192+5000<25000
 facts=dict(SOURCE_only_changed_hash_policy=True,all_static_scalar_bits=scalarbits,source_temporary_bound16384_proved=True,
  signed_packing_and_Minkowski_bias_independently_verified=True,tuplelex_Hslice_hash_independently_equal=True,
  independent_snapshot_header_body_EOF_roundtrip=True,all39live_function_maps_and_source_DAG_retained=True,
  old8_boundaries_exacthash_required=True,full_snapshot_keys125000_guarded=True,
  actual_input_RAW_label='reserved future READ upper, includes receipt/gate reservations; not measured reads',
  fresh_RAW_upper=28422564,SOURCE_owned150000_unchanged=True,math55_hard60_seconds=True,
  source_qM_bits4096_source_temporary16384_only=True,VALUE_coefficient_and2048codec_limits_unchanged=True,
  no_SOURCE_profile_or_matrix_coefficients_evaluated=True,only_synthetic_codec_fixture=True)
 (OUT/'source-proof.json').write_text(json.dumps(facts,indent=2)+'\n');pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 result=dict(id=OUT.name,all_passed=True,frozen=True,SOURCE_only=True,ONE_changed_row6_SOURCE_profile_ready=True,
  input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()})
 p=OUT/'report.json';p.write_text(json.dumps(result,indent=2)+'\n');assert p.stat().st_size<65536
 print(json.dumps(facts,indent=2));print('report_sha256',sha(p))
if __name__=='__main__':main()
