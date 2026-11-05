#!/usr/bin/env python3
"""Independent all-source Schur row and genuine barcube-zero proof, no products."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63_REVIEW';PINS={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):PINS[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
def field(p,key):
 PINS[p]=sha(ROOT/p);s=(ROOT/p).read_text();tok=json.dumps(key);assert s.count(tok)==1;i=s.index(':',s.index(tok))+1
 while s[i].isspace():i+=1
 return json.JSONDecoder().raw_decode(s,i)[0]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30));out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 rp='reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63/report.json';r=read(rp);assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p in r['input_sha256']:
  if p.endswith('/report.json'):assert read(p)['all_passed']
 src=read('reports/AMBIENT_B18_NBAR_CUBE_ZERO_PROJECTION_20261105_v63/all-three-Schur-row42-and-structural-bar-cube-zero.json')
 mp='reports/AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63/actual-minimal-Tnil5-849-and-d3-linear-operator.json'
 nodes=field(mp,'operator_nodes')
 for n in nodes:
  if n[0]=='origin_unit_matrix_inverse':assert Q(n[2]['determinant']);n[2]={'determinant':n[2]['determinant']}
 steps=field(mp,'three_exact_unit_Schur_steps');steps=[{k:s[k]for k in['input','unit_rows','unit_columns','retained_rows','target_projection','unit_origin_determinant']}for s in steps]
 point=read('reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase1-actual-m0-and-section.json');E=point['selected_source_columns24']
 assert 42 not in E and nodes[5]==['QQ_matrix_rows',[[str(int(i==j))for j in E]for i in range(44)]]
 assert nodes[8]==['matrix_product',[5,7]]and nodes[6]==['matrix_product',[0,5]]and nodes[7][:2]==['origin_unit_matrix_inverse',6]
 assert all(Q(c)==0 for c in nodes[5][1][42])
 assert nodes[10]==['horizontal_blocks',[1,8,9]]and nodes[9]==['zero']and nodes[13]==['vertical_blocks',[10,12]]
 xp='reports/AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63/phase2-actual-radical-X-firsts.json';PINS[xp]=sha(ROOT/xp);s=(ROOT/xp).read_text();tok=json.dumps('all780_source_columns44_sparse_dual9');assert s.count(tok)==1;i=s.index(':',s.index(tok))+1
 while s[i].isspace():i+=1
 assert s[i]=='[';i+=1;dec=json.JSONDecoder();stream=hashlib.sha256()
 for j in range(780):
  while s[i].isspace()or s[i]==',':i+=1
  col,i=dec.raw_decode(s,i);found=[jets for row,jets in col if row==42];assert len(found)<=1
  assert not found or Q(dict(found[0]).get(0,0))==0
  stream.update(json.dumps([j,found],separators=(',',':')).encode()+b'\n')
 del s,col
 assert stream.hexdigest()==src['all780_row42_jet_records_sha256']
 labels=list(range(128));route=[]
 for st in steps:
  pos=labels.index(42);assert pos not in st['unit_rows']and pos in st['retained_rows']and Q(st['unit_origin_determinant'])
  pn=nodes[st['target_projection']];assert pn[0]=='matrix_sum';sel,neg=pn[1:]
  assert nodes[sel][0]=='row_column_selection'and nodes[sel][2]==st['retained_rows']and nodes[nodes[sel][1]]==['identity']
  assert nodes[neg][0]=='matrix_neg';prod=nodes[nodes[neg][1]];assert prod[0]=='matrix_product'
  C,U,J=prod[1];assert nodes[C]==['row_column_selection',st['input'],st['retained_rows'],st['unit_columns']]
  assert nodes[U][0]=='origin_unit_matrix_inverse'and nodes[nodes[U][1]]==['row_column_selection',st['input'],st['unit_rows'],st['unit_columns']]
  assert nodes[J][0]=='row_column_selection'and nodes[J][2]==st['unit_rows']and nodes[nodes[J][1]]==['identity']
  nxt=st['retained_rows'].index(pos);route.append([pos,nxt]);labels=[labels[k]for k in st['retained_rows']]
 assert route==[[42,18],[18,0],[0,0]]and labels==[42,43,50,68,70]and nodes[59]==['matrix_product',[58,43,28]]
 archive=field('reports/AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63/actual512-origin-values-and-Pi0.json','constant_minimal_projection5_by128_sparse');assert archive[0]==[[42,'1']]
 phi=read('reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63_REVIEW/independent-eight-fixed-values-and-genuine-quartic-dual.json')
 assert phi['actual_complete_dual1650_sparse']==[[329,'1']]
 deep=read('reports/AMBIENT_B18_NBAR_CUBE_DEEP_SOURCE_20261105_v63_REVIEW/independent-actual-product-order-and-whole-function-H-proof.json')
 assert deep['actual_44_cube_orders_lower_bounds'][42]>=5 and deep['m_cube_order5_by_literal_m0_e43_zero']
 data={'all780_origin_X_row42_stream_hash':stream.hexdigest(),'literal_direct_section_row42_zero':True,'raw42_retention_route':route,'three_unit_Schur_projection_rows_at_origin_are_selectors':True,
 'independent_full_Pi0_firstrow128_sparse':[[42,'1']],'archived_row_matches_new_source_derivation':True,'genuine_complete_phi4_1650':[[329,'1']],
 'actual_all_bar_cubes_zero_proof':'Every literal triple of the genuine Nbar14 basis has an actual128 representative (c,H(mc)) of order>=4. The independently proven raw42 component has order>=5. Minimal projection is regular and its constant first row is e42. Thus the first target r4 coefficient is0; this is the COMPLETE actual quartic dual on all1650 monomials, so the m4 class is0. Actual m5Tnil=0 kills higher terms. Multilinearity quantifies all barcubes; no Tnil ring structure is used.',
 'D_finite_intersection_proof':'Let H=tors_D(N). If ell*n∈H, then ell^q(ell*n)=0 for someq, so n∈H. Therefore ellN∩H=ellH. Known N³⊂H and new Nbar³=0 give N³⊂ellH.',
 'conditional_index3':'Only under the separate accepted condition ellH=0 (equiv honest T2 length17) does originalN³=0 follow; knownN²≠0 then gives index3. The condition is not proved here.',
 'original_nil_index_unchanged':[3,7],'cube_scalars_or_new_jets_evaluated':False,'inactive9_or35_tests_unnecessary':True}
 ep=out/'independent-Pi0-row42-and-all-genuine-bar-cubes-zero-proof.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 result={'id':STEM,'all_passed':True,'status':'PROVED_independent_all_genuine_bar_nil_cubes_zero','genuine_Nbar_cube_zero':True,'original_cube_subset_ell_H0':True,'original_cube_vanishing_proved':False,'input_sha256':PINS,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));signal.alarm(0)
if __name__=='__main__':main()
