#!/usr/bin/env python3
"""Independent fixed homogeneous SOURCE and unrun streaming-worker audit."""
import ast,hashlib,itertools,json,math,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_CUBIC_ELL3_SOURCE_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def plus(v,w,c=F(1)):
 z=dict(v)
 for i,a in w.items():
  z[i]=z.get(i,F(0))+c*a
  if not z[i]:del z[i]
 return z
def sp(v):return [[i,c]for i,c in sorted(v.items())]
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def checkpins(d):
  for p,h in d.get('input_sha256',{}).items():assert sha(ROOT/p)==h;pins[p]=h
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_CUBIC_ELL3_SOURCE_PREFLIGHT_20261105_v63';r=load(primary,h='3e6087d41837068f6fe9f7513f89f154d127476ebeaed0abbafb7e947ed752b9');checkpins(r)
 p=load(primary,'actual-rank1-quadratic-source-and-streamed-cubic-gate.json',h='5c42a527ff3d5f0600c46b5edf33117eadbc22b234699f9ee9758a394e0c156b')
 qp=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW',h='5faf10001ed74c2b9e107d7fa1f4fe873600423affa75305a75a8f1d3fad5bb3');assert qp['all_passed']and qp['actual_Cprime_HF_through2']==[3,3,3]
 V=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63','phase4-all-twenty-actual-quadratic-values.json');L=[pairs(v)for v in load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase1-actual-fullnine-first-V-columns.json')['all44_actual_linear_columns27_sparse']]
 Q4=load('AMBIENT_B18_CPRIME_PRODUCT_ANN_QUOTIENT_SOURCE_PREFLIGHT_20261105_v63');q4=[pairs(v)for v in Q4['literal_fixed_pi4_rows_by135_sparse']];v=list(map(F,V['all20_actual_quadratic_values4'][18]));assert v==[0,F(275,96),F(35,48),F(13,300)]
 free=[0,2,3];q3=[plus(q4[i],q4[1],-v[i]/v[1])for i in free];E=[Q4['unit_NF_rows135'][i]for i in free];assert E==[44,133,134]
 assert q3==[pairs(v)for v in p['exact_NF3_rows_by135_sparse']]and [[q.get(i,F(0))for i in E]for q in q3]==[[F(i==j)for j in range(3)]for i in range(3)]
 m2=list(itertools.combinations_with_replacement(range(9),2));i2={m:i for i,m in enumerate(m2)};m3=list(itertools.combinations_with_replacement(range(9),3));i3={m:i for i,m in enumerate(m3)};K=pairs(V['source_constant_kernel20'][18]);P={}
 for j,c in K.items():
  for target in range(3):
   for k,mon in enumerate([(7,7),(7,8),(8,8)]):P=plus(P,{45*target+i2[mon]:c*F(V['all44_actual_V2_pureplane_nine_values'][j][3*target+k])})
 assert P==pairs(p['actual_quadratic_initial_representative135_sparse'])and [sum((c*P.get(i,F(0))for i,c in row.items()),F(0))for row in q4]==v
 digest=hashlib.sha256()
 for col in L:
  for var in range(9):
   z={}
   for i,c in col.items():t,x=divmod(i,9);z=plus(z,{45*t+i2[tuple(sorted((x,var)))]:c})
   assert all(sum((c*z.get(i,F(0))for i,c in row.items()),F(0))==0 for row in q3);digest.update(json.dumps(enc(sp(z)),separators=(',',':')).encode()+b'\n')
 assert digest.hexdigest()==p['known396_products_sha256']and all(sum((c*P.get(i,F(0))for i,c in row.items()),F(0))==0 for row in q3)
 proj=[]
 for target in range(3):
  for mon in m3:
   d=45*target+i2[mon[:2]];x=mon[2];proj.append({9*j+x:q.get(d,F(0))for j,q in enumerate(q3)if q.get(d,F(0))})
 assert proj==[pairs(v)for v in p['canonical_cubic_projection_columns495_sparse']]
 def embed(j,x):
  target,mi=divmod(E[j],45);return 165*target+i3[tuple(sorted(m2[mi]+(x,)))]
 # Independent ordinary monomial multiplyback for all495 source generators.
 dh=hashlib.sha256()
 for raw in range(495):
  target,j=divmod(raw,165);mon=m3[j];d=45*target+i2[mon[:2]];x=mon[2];lhs={raw:F(1)}
  for k,c in proj[raw].items():j,var=divmod(k,9);lhs=plus(lhs,{embed(j,var):c},-1)
  quad={d:F(1)}
  for j,q in enumerate(q3):quad=plus(quad,{E[j]:q.get(d,F(0))},-1)
  assert all(sum((c*quad.get(i,F(0))for i,c in q.items()),F(0))==0 for q in q3)
  rhs={}
  for i,c in quad.items():t,j=divmod(i,45);rhs=plus(rhs,{165*t+i3[tuple(sorted(m2[j]+(x,)))]:c})
  assert lhs==rhs;dh.update(json.dumps(enc([raw,sp(quad),sp(lhs)]),separators=(',',':')).encode()+b'\n')
 consistency=[plus({9*j+x:F(1)},proj[embed(j,x)],-1)for j in range(3)for x in range(9)]
 assert len(consistency)==27 and p['source_count']==dict(first44_times_S2=1980,actual_quadratic_times_S1=9,ell_cube_targets=3,consistency=27,total=2019)
 old=load('AMBIENT_B18_CPRIME_CUBIC_ELL3_PREPARED_WORKER_20261105_v63',h='888486a7202646c25aa83df8cc8734cb4cfec000f301eeb0888b755aac113d9c');checkpins(old)
 prep=load('AMBIENT_B18_CPRIME_CUBIC_ELL3_PREPARED_WORKER_V2_20261105_v63',h='03426aa5f37e4332eb624620fd05ba99929ae95c4f76d6dda6b6946842ec09a4');checkpins(prep)
 worker=ROOT/prep['prepared_code'];assert sha(worker)==prep['prepared_code_sha256']=='81d350d0ab1ab37ae0783b1bac2ff9306889163d8f94bd628cc4866b9cdcef2d'
 for n,h in old['evidence_sha256'].items():assert sha(ROOT/'reports/AMBIENT_B18_CPRIME_CUBIC_ELL3_PREPARED_WORKER_20261105_v63'/n)==h
 helper=ROOT/'code/ambient_B18_Cprime_cubic_ell3_source_preflight_20261105_v63.py';assert sha(helper)==r['code_sha256']=='13be0d7bb4d1f515649f5a9d72cd668651a5fcd7230b0e2d4db6bdcda593f7ea'
 # Coefficient bounds are proved without generating/ranking2019 columns.
 # Raw multiplication is injective on monomials; each projected entry is
 # bounded by a source column l1 norm times max absolute projection entry.
 qc=[c for row in q3 for c in row.values()];ac=[c for col in L for c in col.values()]+list(P.values());den=math.lcm(*(c.denominator for c in [a*b for a in ac for b in qc]+qc+[F(1)]));assert den==13365000000
 maximum=max(max(sum(map(abs,col.values()))for col in L),sum(map(abs,P.values())),F(2))*max(map(abs,qc));M=math.ceil(den*maximum);assert M==184057500000<=1453443750000
 bits=(math.factorial(27)*M**27).bit_length();derived=bits+2*den.bit_length()+50;assert bits==1104 and derived==1222<1270<2048
 profile=prep['source_size_profile'];assert [profile[k]for k in ['all_input_and_dependency_metadata_upper_including_1742','rank_scratch_upper_including_fixed_layout739','saved_evidence_and_controller_upper']]==[6000,6000,10000]and profile['sum_of_three_disjoint_envelopes']==22000<profile['hard_worker_cap']==24000
 text=worker.read_text();ast.parse(text)
 for token in ['del piv','inverse=[row[n:]for row in aug];del aug','count==2019','count2==2019','peak<=24000','signal.alarm(60)','<=2048',"sourcehash.hexdigest()==backsource.hexdigest()"]:assert token in text
 valueout=ROOT/'reports/AMBIENT_B18_CPRIME_CUBIC_ELL3_VALUES_20261105_v63';assert not valueout.exists();run=subprocess.run(['python',str(worker),'--help'],capture_output=True,text=True,timeout=15,cwd=ROOT);assert run.returncode==0 and not run.stderr and not valueout.exists();(OUT/'independent-help-stdout.txt').write_text(run.stdout)
 save('independent-initial-source-projection-custody-and-positive-criterion.json',dict(actual_K18_pureplane_representative=P,rank1_residue=v,NF3=q3,free2_rows=E,all396_known_quad_products=True,all495_monomial_divisor_multiplybacks=True,divisor_proof_sha256=dh.hexdigest(),all27_consistency_columns=consistency,ann84_redundant_proof='TheirNF4imageis0 andkerNF4isexactknownrank131linearproducts; their756cubicmultiplesalready liein44*S2.',initial_ancestry='The actualV2*K18 andP have identicalNF4. Their omitted nonplane part lies in knownK2, soP itself is an actualquadraticinitial; kerNF3=knownK2+QQP.',compressed_cokernel_proof='Id-EQ maps everyF495generator intoS1kerNF3 by the saved495monomialmultiplybacks. IncludingALL27I-QE andall1980+9+3 projectedknownsources gives the exactknowncubic quotient ofCprime/ell3.',source_count=2019,independent_integer_denominator=den,independent_integer_coefficient_bound=M,independent_minor_bits=bits,independent_derived_bits=derived,disjoint_lifetime_envelopes=[6000,6000,10000],whole_envelope=22000,hard_cap=24000,preserved_old_custody_not_used=True,independent_help_import_only=True,positive_criterion='Only rank27 implies m3(Cprime/ell3)=0 byNakayama and length<=9. ActualHF333 plusDVRrank2/torsion{1,e>=2} give length>=9, hence length9/e2. A nonzero knowncubic quotient is inconclusive.',no_cubic_rank_or_source_values_computed=True))
 for q in [Path(__file__),worker,helper,ROOT/f'research/notes/{ID}.md']:pins[str(q.relative_to(ROOT))]=sha(q)
 raw=sum((ROOT/q).stat().st_size for q in pins)+sum(q.stat().st_size for q in OUT.iterdir()if q.is_file());assert raw<32<<20 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,actual_initial_source_and_compressed_presentation_accepted=True,streaming_worker_custody_and_bounds_pass=True,known_cubic_columns=2019,compressed_rows=27,whole_numeric_envelope=22000,hard_numeric_cap=24000,no_cubic_value_or_rank_run=True,parent_value_GO_still_required=True,input_sha256=pins,evidence_sha256={q.name:sha(q)for q in OUT.iterdir()if q.is_file()and q.name!='report.json'},actual_input_output_RAW=raw,seconds=time.monotonic()-st))
 print(json.dumps(dict(all_passed=True,columns=2019,rows=27,no_rank_run=True,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
