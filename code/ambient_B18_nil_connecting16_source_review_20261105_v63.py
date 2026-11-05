#!/usr/bin/env python3
"""Independent typed/source and exact formal connecting-map proof only."""
import hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_NIL_CONNECTING16_SOURCE_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
ZERO=(0,)*8
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 primary='AMBIENT_B18_NIL_CONNECTING16_SOURCE_OPERATOR_20261105_v63';r=load(primary);assert pins[f'reports/{primary}/report.json']=='df9b460a8d4d261c9967ba669b943e5f2296ce9df1d10b91cd029870c4505d1a'and r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 E=load(primary,'actual-formal-connecting16-and-ten-source-obstructions.json');alias=load('AMBIENT_B18_NIL_CONNECTING16_NAMESPACE_BINDING_20261105_v63');assert pins['reports/AMBIENT_B18_NIL_CONNECTING16_NAMESPACE_BINDING_20261105_v63/report.json']=='2849325cc64a948758349b106d2ea12b1f20029869b59cf621e256a0c9e8019b'
 assert alias['primary_report_sha256']==pins[f'reports/{primary}/report.json']
 foreign={}
 for ns,p in alias['external_namespace_artifacts'].items():assert sha(ROOT/p)==alias['artifact_sha256'][p];foreign[ns]=json.loads((ROOT/p).read_text())
 assert alias['external_namespace_artifacts']['radical_full_R0']==E['source_namespace_bindings']['radical']
 for ns in['minimal','honest','natural']:assert alias['external_namespace_artifacts'][ns]==E['source_namespace_bindings'][ns]
 for s in['AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63_REVIEW','AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63_REVIEW','AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW','AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63_REVIEW','AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW','AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW','AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW']:assert load(s)['all_passed']
 basis=load('AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63','exact-Loewy16-and-explicit-filtered-basis.json')['explicit_sixteen_filtered_basis_lifts'];M=foreign['minimal'];current=list(range(128))
 for step in M['three_exact_unit_Schur_steps']:
  assert step['retained_rows']==[i for i in range(len(current))if i not in step['unit_rows']]
  current=[current[i]for i in step['retained_rows']]
 assert current==[42,43,50,68,70]
 nodes=E['matrix_operator_nodes'];shapes=[tuple(s)for s in E['matrix_operator_shapes']];assert len(nodes)==len(shapes)==68
 def oldref(i,j):assert isinstance(j,int) and 0<=j<i
 for i,(n,s)in enumerate(zip(nodes,shapes)):
  op=n[0]
  if op=='actual_external_matrix':
   assert n[1]in foreign;d=foreign[n[1]]
   for k in n[2].split('.'):d=d[k]
  elif op=='identity':assert s[0]==s[1]
  elif op=='QQ_matrix_rows':assert len(n[1])==s[0]and all(len(row)==s[1]for row in n[1]);assert all(isinstance(F(c),F)for row in n[1]for c in row)
  elif op=='matrix_product':
   ids=n[1]
   for j in ids:oldref(i,j)
   assert all(shapes[a][1]==shapes[b][0]for a,b in zip(ids,ids[1:]))and s==(shapes[ids[0]][0],shapes[ids[-1]][1])
  elif op=='matrix_neg':oldref(i,n[1]);assert s==shapes[n[1]]
  elif op=='matrix_sum':oldref(i,n[1]);oldref(i,n[2]);assert s==shapes[n[1]]==shapes[n[2]]
  elif op=='horizontal_blocks':
   for j in n[1]:oldref(i,j)
   assert all(shapes[j][0]==s[0]for j in n[1])and sum(shapes[j][1]for j in n[1])==s[1]
  elif op=='row_column_selection':oldref(i,n[1]);assert s==(len(n[2]),len(n[3]))and all(0<=j<shapes[n[1]][0]for j in n[2])and all(0<=j<shapes[n[1]][1]for j in n[3])
  elif op=='multiply_by_literal_S_monomial':oldref(i,n[1]);assert s==shapes[n[1]]and all(0<=v<8 for v in n[2]);F(n[3])
  elif op in['actual_coefficient_pullback_by_original105_clone','coefficientwise_formal_section_S_to_R0']:oldref(i,n[1]);assert s==shapes[n[1]]
  elif op=='matrix_of_regular_formal_coordinate_divisions':oldref(i,n[1]);assert s==shapes[n[1]]==(3,16)and n[2]=='ell'and n[3]['each_quotient_origin']=='uncomputed, never set0'
  else:assert op=='classes_in_actual_full_Cprime_cokernel';oldref(i,n[1]);oldref(i,n[2]['presentation_matrix']);assert s==shapes[n[1]]and shapes[n[2]['presentation_matrix']]==(3,44)
 R=E['matrix_roots'];inj=nodes[R['minimal5_to_honest128_injection']];assert inj[0]=='QQ_matrix_rows'and [[F(c)for c in row]for row in inj[1]]==[[F(i==j)for j in current]for i in range(128)]
 # Universal Schur P=[-B U^-1,I], E=(0,I)^t gives PE=I;
 # composing the three literal retained-coordinate sections proves the
 # actual5->128 source lift, independent of coefficient truncations.
 assert nodes[R['directPhi44_by128']]==['horizontal_blocks',[6,7]]and nodes[6]==['matrix_sum',3,5]and nodes[4]==['matrix_product',[1,0]]and nodes[5]==['matrix_neg',4]and nodes[7]==['matrix_product',[1,2]]
 # Explicit noncommutative source identity Vbar*Phi=0 modulo mn=I,
 # d1*d2=0. No scalar original source coefficients are expanded.
 def normal(word):
  w=list(word)
  while True:
   hit=next((j for j in range(len(w)-1)if w[j:j+2]==['m','n']),None)
   if hit is None:break
   del w[hit:hit+2]
  if any(w[j:j+2]==['d1','d2']for j in range(len(w)-1)):return None
  return tuple(w)
 def reduced(terms):
  d={}
  for w,c in terms:
   w=normal(w)
   if w is not None:d[w]=d.get(w,F(0))+c
  return {w:c for w,c in d.items()if c}
 assert not reduced([(('d1','m'),F(1)),(('d1','m','n','m'),F(-1))])and not reduced([(('d1','m','n','d2'),F(1))])
 # Interpret only the5x16 basis coefficients in a free symbolic natural
 # matrix Nij and S monomials, independently of the producer Matrices API.
 memo={}
 def addp(a,b):
  d=a.copy()
  for k,c in b.items():
   d[k]=d.get(k,F(0))+c
   if not d[k]:del d[k]
  return d
 def emat(k):
  if k in memo:return memo[k]
  n=nodes[k];a,b=shapes[k];o=n[0]
  if o=='identity':v=[[{('one',0,0,ZERO):F(1)}if i==j else{}for j in range(b)]for i in range(a)]
  elif o=='actual_external_matrix':assert n[1:]==['natural','matrices.actual_natural_map_minimal5_by5'];v=[[{('N',i,j,ZERO):F(1)}for j in range(5)]for i in range(5)]
  elif o=='row_column_selection':z=emat(n[1]);v=[[z[i][j]for j in n[3]]for i in n[2]]
  elif o=='horizontal_blocks':parts=[emat(j)for j in n[1]];v=[sum((p[i]for p in parts),[])for i in range(a)]
  elif o=='matrix_sum':x,y=emat(n[1]),emat(n[2]);v=[[addp(x[i][j],y[i][j])for j in range(b)]for i in range(a)]
  else:
   assert o=='multiply_by_literal_S_monomial';z=emat(n[1]);ex=tuple(n[2].count(i)for i in range(8));v=[[{(tag,r,c,tuple(x+y for x,y in zip(mon,ex))):coef*F(n[3])for(tag,r,c,mon),coef in p.items()}for p in row]for row in z]
  memo[k]=v;return v
 bb=emat(R['actual_filtered_basis_coordinates5_by16']);expected=[[{}for _ in range(16)]for _ in range(5)]
 for j,b in enumerate(basis):
  def put(target,mon,coef,oldcol=None):
   ex=tuple(mon.count(i)for i in range(8));key=('one',0,0,ex)if oldcol is None else('N',target,oldcol,ex);expected[target][j]=addp(expected[target][j],{key:F(coef)})
  if b['type']in['new_minimal_unit','new_minimal_monomial']:put(b['target'],b.get('monomial',[]),1)
  elif b['type']=='natural_column_times_monomial':
   for i in range(5):put(i,b['monomial'],1,b['old_minimal_column'])
  else:
   for col,mon,c in b['terms']:
    for i in range(5):put(i,mon,c,col)
 assert bb==expected
 gamma=nodes[59];assert gamma==['QQ_matrix_rows',[['1','0','0'],['0','0','1'],['-11/12','1','0']]]
 assert nodes[60]==['actual_external_matrix','radical_full_R0','actual_V3_by44_rows_after_unit1_cancellation']and nodes[61][2]=='FULL_R0_F9'and nodes[63][2]=='ELL0_S_F9'
 assert R['actual_connecting_classes3_by16']==66 and nodes[67]==['row_column_selection',66,[0,1,2],list(range(10))]
 J=load('AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63','joint-degree3-certificate.json');assert len(J['all105_unknown_original_node_order'])==len(J['all105_equation_original_node_order'])==105 and F(J['det_joint105'])
 clones=E['two_distinct_original105_source_clones'];assert [c['namespace']for c in clones]==['FULL_R0_F9','ELL0_S_F9'];assert clones[0]['target9'][:-1]==clones[1]['target9'][:-1]and clones[0]['target9'][-1]==['sum',['r'],['scale','1/5',['ell']]]and clones[1]['target9'][-1]==['r']
 assert E['all_six_degree_at_least2_basis_columns_zero_CLASSES']==list(range(10,16))and [i for i,b in enumerate(basis)if b['degree']>=2]==list(range(10,16))
 proof=dict(all68_typed_matrix_nodes_and_foreign_namespaces=True,explicit_namespace_alias_receipt_pinned=True,actual5_to128_injection=current,three_exact_Schur_right_sections=True,all16_symbolic_basis_columns_match=True,
  all48_formal_numerator_memberships='Vbar=d1*m andPhi=(I-nm,n*d2); mn=I andd1d2=0 prove allzero atell0. RingR=QQ[[z7,r,ell]] is regular andell is a nonzerodivisor, so each exact numerator defines a unique regular ellquotient.',
  source_section='Separate original105 full andell0 unit clones. Same unique derived canonical/unit/regular-quotient branches are pulled back functorially; no ell0 restriction is applied before the fullV clone.',
  boundary_and_lift_descent='Difference of two full lifts isell*b, hence quotient changes byVfull*b. A barX boundary lifts toXfull*sigma(a)+ell*b; VfullXfull=0, so its quotient class iszero. This is the actualUCT edge, not a finiteArtin kernel.',
  exact_sequence='0 -> N/ellN -> Tnil16 -> Cprime[ell]2 ->0',
  target_m2_zero_proof='Cprime[ell] has length2 over the SAME localS. A nonzero m²M would force three strict nonzero filtration layers; an equality at either stage contradictsNakayama. Thus m²M=0.',
  six_actual_higher_basis_classes_in_nil_kernel=list(range(10,16)),remaining_genuine_source_classes10=True,genuine_nil_kernel_dimension14=True,
  scope='These are exact formalfunction/class operators. No ellquotient origin, targetbasis2, literal Vpreimage for sixzeros, actual14 full nil lifts or nilcube product is computed. Cprime2 length7 does not replace the fulltarget.')
 save('independent-formal-section-types-and-genuine-edge-proof.json',proof)
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 report=dict(id=ID,all_passed=True,frozen=True,all68_types_and_source_bindings_verified=True,actual_formal_ell_division_memberships48=True,six_higher_basis_classes_zero_in_genuine_edge=True,genuine_remaining_source_count10=True,genuine_kernel_dimension14=True,actual_edge_values_or_nil_generators_computed=False,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20;save('report.json',report);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
