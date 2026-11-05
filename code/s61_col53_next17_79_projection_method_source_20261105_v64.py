"""Pure SOURCE projection proof/fixture and exact request custody, no actual coefficients."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
from s61_lex_HL_S2_source_projection_20261105_v64 import input_project,copy_project,product_project,sum_project,inside
ID='S61_COL53_NEXT17_79_PROJECTION_METHOD_SOURCE_20261105_v64'
PREFLIGHT='reports/S61_COL53_NEXT17_79_FUNCTIONAL_SOURCE_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
if __name__=='__main__':
 p=ROOT/PREFLIGHT/'report.json';assert sha(p)=='359acee48739dc5e99b14bd24d9a65791fd9704510b6cc8ad308b7309b9b870b'
 r=json.loads(p.read_text());assert r['all_passed']and r['SOURCE_PREFLIGHT_only_no_signed_coefficients']
 pins={str(p.relative_to(ROOT)):sha(p)}
 for f,h in r['input_sha256'].items():assert sha(ROOT/f)==h;pins[f]=h
 for f,h in r['evidence_sha256'].items():assert sha(p.parent/f)==h;pins[str((p.parent/f).relative_to(ROOT))]=h
 D=json.loads((p.parent/'chosen-functional-and-thin-source-inputs.json').read_text());target=D['new_functional_request']['exponent5'];corner=D['new_functional_request']['source_positive_minor_corner5']
 assert target==[-5,0,9,-3,19]and corner==[-5,0,7,-3,19]
 assert D['new_functional_request']['target_minus_positive_corner5']==[0,0,2,0,0]
 profiles=D['positive_source_profiles'];a=(-3,0,4);b=tuple(profiles['17,53']['positive_lex_corner5'][:3]);c=tuple(profiles['6,79']['positive_lex_corner5'][:3]);d=tuple(profiles['17,79']['positive_lex_corner5'][:3])
 assert tuple(a[j]+d[j]for j in range(3))==tuple(b[j]+c[j]for j in range(3))==tuple(corner[:3])
 # Exhaustive finite formal-exponent fixture families. These are not source/beta
 # points. Include high-H/negative-L and high-L/negative-s terms.
 checks=0;guardmax=0
 def guard(n):
  global guardmax
  guardmax=max(guardmax,n);assert n<=4096
 def source_set(l,tag):
  h,L,s=l
  S={(h,L,s+i,tag-i,2*i-tag)for i in range(5)}
  S|={(h,L+1,s-9,tag,-7),(h+1,L-8,s-20,-9,tag)}
  return S
 for la in[(0,0,0),(-3,0,4),(2,-4,-7)]:
  for lb in[(0,0,0),(-2,0,3),(1,-6,5)]:
   A=source_set(la,1);B=source_set(lb,3);A2=input_project(A,la,guard);B2=input_project(B,lb,guard)
   lp=tuple(la[j]+lb[j]for j in range(3))
   full={tuple(e[j]+f[j]for j in range(5))for e in A for f in B}
   result,pairs=product_project(A2,B2,la,lb,lp,guard)
   assert result=={e for e in full if inside(e,lp)};checks+=1
   for shift in[(0,0,0,0,0),(2,-3,1,-6,7)]:
    lo=tuple(la[j]+shift[j]for j in range(3));cop=copy_project(A2,la,lo,shift,guard)
    assert cop=={tuple(e[j]+shift[j]for j in range(5))for e in A if inside(e,la)};checks+=1
 for delta in[(0,0,0),(0,0,1),(0,0,2),(0,0,3),(0,1,-20),(1,-20,-30)]:
  la=(-3,0,4);lb=tuple(la[j]+delta[j]for j in range(3));A=source_set(la,1);B=source_set(lb,2)
  lo=min(la,lb);result=sum_project([input_project(A,la,guard),input_project(B,lb,guard)],[la,lb],lo,guard)
  assert result=={e for e in A|B if inside(e,lo)};checks+=1
 helper=ROOT/'code/s61_lex_HL_S2_source_projection_20261105_v64.py'
 note=ROOT/'research/notes'/f'{ID}.md'
 for q in[Path(__file__),helper,note]:pins[str(q.relative_to(ROOT))]=sha(q)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 proof={'exact_projection_statement':'For a priori SOURCE lex-leading triples lA,lB, leading-H,L and relative-s degree<=2 of a product uses only the same projections of both factors; sums shift relative-s offsets and omit only SOURCE-higherH/L or degree>2 children. Actual cancelled faces never shift the supplied SOURCE leads.',
 'proof':'A SOURCE exponent minus its supplied lead is lex-nonnegative in Z^3. In a product that reaches total H-deficit0 both H-deficits are0; then total L-deficit0 forces both L-deficits0; remaining s-deficits are nonnegative and sum<=2. Addition/copy follow the same fixed coordinates. q,t remain untouched.',
 'chosen_minor_leading_triple':corner[:3],'requested_full5_exponent':target,'relative_s_degree':2,'product_degree_pairs':[[0,0],[0,1],[1,0],[0,2],[1,1],[2,0]],'FULL_F11_NOT_REPLACED':True,'actual_zero_faces_or_vertices_never_used_for_pruning':True,'actual_coefficients_evaluated':False,'formal_exponent_fixtures':checks,'synthetic_max_support':guardmax,'future_graph_compilation_lastuse_owned_RAW_and_VALUE_gate_STILL_REQUIRED':True}
 (out/'projection-proof-and-fixed-functional.json').write_bytes(enc(proof))
 report={'id':ID,'all_passed':True,'SOURCE_math_and_implementation_fixture_ONLY':True,'no_actual_minor_coefficient_or_SOURCE_support_run':True,'candidate':[17,79],'functional_exponent5':target,'all5coords_kept':True,'full_F11_and_all_complements_retained':True,'projection_helper_immutable':True,'future_actual_projection_SOURCE_worker_gate_pending':True,'input_sha256':pins,'evidence_sha256':{'projection-proof-and-fixed-functional.json':sha(out/'projection-proof-and-fixed-functional.json')}}
 (out/'report.json').write_bytes(enc(report));print('report_sha256',sha(out/'report.json'));print('formal fixtures',checks)
