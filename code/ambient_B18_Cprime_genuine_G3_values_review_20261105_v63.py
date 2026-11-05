#!/usr/bin/env python3
"""Independent cached normalform/action and genuine target verification."""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def pairs(v):return {int(i):F(c)for i,c in v if F(c)}
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(s,n='report.json',h=None):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p)
  if h:assert sha(p)==h
  return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_GENUINE_G3_VALUES_20261105_v63';r=load(primary,h='0936b3983bcd06853fbd9c2e21871eb856cc0da1a8248344e7b7da471d33d2b8');assert r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 e=load(primary,'actual-Cprime9-actions-and-genuine-G3-target.json',h='f73869c4f307de530a1dd32f181dec5e6b1613bd120f9abb226b757663605dcb')
 assert load('AMBIENT_B18_CPRIME_GENUINE_G3_SOURCE_20261105_v63_REVIEW',h='b8109751d709bde8cf9abb9e0427fc6dfe2894b589c27ad57f1ea1c7bf926b7c')['cached9_action_source_and_genuine_G3_constructor_accepted']
 assert load('AMBIENT_B18_CPRIME_EXACT_TORSION_12_20261105_v63_REVIEW',h='d7ec8b85b6622793fd4d5c06dfe96b1fe4b02d7ddea9ef9be3ebcfd837587f0b')['genuine_target_kerell_mod_imell2_canonical']
 source=load('AMBIENT_B18_CPRIME_GENUINE_G3_SOURCE_PREFLIGHT_20261105_v63','actual-nine-action-and-genuine-G3-source-interface.json');Q=load('AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63','phase4-all-twenty-actual-quadratic-values.json');L=[pairs(v)for v in source['source_L44_columns27_sparse']];cp=source['source_V1_pivot24_columns'];rr=source['source_V1_pivot24_rows'];iv=[[F(c)for c in row]for row in source['source_V1_fixed24_inverse']];q1=[pairs(v)for v in source['first_NF3_rows27_sparse']];q2=[pairs(v)for v in source['quadratic_NF3_rows135_sparse']];E1=source['first_basis_embedding_rows27'];E2=source['quadratic_basis_embedding_rows135'];m2=list(itertools.combinations_with_replacement(range(9),2));i2={m:i for i,m in enumerate(m2)}
 B=[]
 for j,col in enumerate(Q['all44_actual_V2_pureplane_nine_values']):
  p={45*t+i2[mon]:F(col[3*t+k])for t in range(3)for k,mon in enumerate([(7,7),(7,8),(8,8)])};v=[sum((c*p.get(i,F(0))for i,c in row.items()),F(0))for row in q2];nf=list(map(F,Q['all44_projected_NF4_values'][j]));assert v==[nf[0],nf[2]-F(14,55)*nf[1],nf[3]-F(104,6875)*nf[1]];B.append(v)
 def nf(u0,u1,u2):
  coefficients=[sum((c*u1.get(rr[j],F(0))for j,c in enumerate(row)),F(0))for row in iv]
  linear=[sum((c*u1.get(i,F(0))for i,c in row.items()),F(0))for row in q1]
  quadratic=[sum((c*u2.get(i,F(0))for i,c in row.items()),F(0))-sum((a*B[cp[j]][k]for j,a in enumerate(coefficients)),F(0))for k,row in enumerate(q2)]
  return [u0.get(i,F(0))for i in range(3)]+linear+quadratic
 basis=[(t,())for t in range(3)]+[(i//9,(i%9,))for i in E1]+[(i//45,m2[i%45])for i in E2];actions=[]
 for x in range(9):
  cols=[]
  for t,mon in basis:
   m=tuple(sorted(mon+(x,)));cols.append(nf({}, {9*t+m[0]:F(1)}if len(m)==1 else{}, {45*t+i2[m]:F(1)}if len(m)==2 else{}))
  actions.append(cols)
 assert actions==[[list(map(F,col))for col in A]for A in e['actual_R0_actions9_by9_columns']]
 zero=[F(0)]*9;unit=lambda i:[F(j==i)for j in range(9)]
 def apply(A,v):return [sum((c*A[j][i]for j,c in enumerate(v)if c),F(0))for i in range(9)]
 rh=hashlib.sha256()
 for j,v in enumerate(L):
  p={E2[k]:c for k,c in enumerate(B[j])if c};a=nf({},v,p);assert a==zero;rh.update(json.dumps(enc(a)).encode())
 assert rh.hexdigest()==e['all44_actual_relation_reduction_sha256']
 ch=hashlib.sha256()
 for A in actions:
  for C in actions:
   for j in range(9):p=apply(A,C[j]);assert p==apply(C,A[j]);ch.update(json.dumps(enc(p),separators=(',',':')).encode())
 assert ch.hexdigest()==e['all81_commutation_checks_sha256']
 cyclic=[unit(t)for t in range(3)]
 for t,mon in basis[3:]:
  p=unit(t)
  for x in mon:p=apply(actions[x],p)
  cyclic.append(p)
 assert cyclic==[unit(j)for j in range(9)]and all(A[j]==zero for A in actions for j in range(6,9))and all(not any(A[j][:6])for A in actions for j in range(3,6))and all(not any(col[:3])for A in actions for col in A)
 E=actions[8];E_sq=[apply(E,col)for col in E];assert E_sq==[unit(6),zero,unit(8)]+[zero]*6==[list(map(F,col))for col in e['ell_square_columns9']]
 rows=[3,5,6,7,8];cols=[0,2,3,4,5];minor=[[E[j][i]for j in cols]for i in rows];assert all(minor[i][i]==1 and all(minor[i][j]==0 for j in range(i+1,5))for i in range(5))
 S=[[E[j][6+i]for j in range(3)]for i in range(3)];K=unit(1)
 for j in range(3):K[3+j]=-S[j][1]
 kernel=[K,unit(6),unit(7),unit(8)];assert kernel==[list(map(F,v))for v in e['all4_actual_kernel_ell_basis_vectors9']]and all(apply(E,k)==zero for k in kernel)
 assert [[k[i]for i in [1,6,7,8]]for k in kernel]==[unit(i)[:4]for i in range(4)]
 assert [unit(6),unit(8)]==[list(map(F,v))for v in e['image_ell_square_basis_vectors9']]
 G=[K,unit(7)];project=lambda v:[v[1],v[7]];assert G==[list(map(F,v))for v in e['genuine_Cprime_ell_torsion_basis_vectors9']]and [project(k)for k in G]==[[1,0],[0,1]]
 target=[];alpha=[]
 for A in actions:
  ak=apply(A,K);assert apply(E,ak)==zero and not any(ak[:6])and apply(A,G[1])==zero
  a=project(ak);assert a[0]==0;target.append([a,[F(0),F(0)]]);alpha.append(a[1])
 assert target==[[list(map(F,col))for col in A]for A in e['actual_genuine_G3_R0_actions2_by2_columns']]and alpha==list(map(F,e['actual_genuine_G3_action_scalars9']))==[F(0)]*7+[F(-11,16),F(0)]
 save('independent-actual-nine-module-and-genuine-two-target.json',dict(all9_actions=actions,all44_actual_relation_hash=rh.hexdigest(),all81_commutator_hash=ch.hexdigest(),all9_literal_cyclic_vectors_I9=True,strict_three_layer_support=True,ell_rank5_unit_lowertriangular_minor=minor,ell_square_rank2_I2=True,full_kernel4=kernel,image_ell_square2=[unit(6),unit(8)],genuine_target_basis=G,genuine_projection_on_kernel=[1,7],all9_genuine_target_actions=target,actual_genuine_action_scalars=alpha,genuine_R0_module='QQ[[rA]]/(rA²) withallz_i andellacting0; K isacyclicgenerator,rA*K=−11Q1/16. ThisrAistheS cutcoordinate,notDVRell.',genuine_annihilator='(zA0,...,zA6,ell,rA²)',actual_fullCprime_ell_torsion_identified=True,no_finite_kernel_free_block_confusion=True,no_edge_RegDivell_or_original_nil_generators_computed=True,no_new_source_gauge_jet_or_rank_search=True))
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 total=sum((ROOT/p).stat().st_size for p in pins)+sum(p.stat().st_size for p in OUT.iterdir()if p.is_file());assert total<500000 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,actual_Cprime9_actions_and_all44_relations_verified=True,genuine_Cprime_ell_torsion_dimension=2,genuine_module_rA_squared_cyclic=True,genuine_target_action_scalars=list(map(str,alpha)),all81_commutators_and_fixed_ell_blocks_verified=True,no_edge_or_new_source_gauge_jet_values=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},actual_RAW=total,seconds=time.monotonic()-st))
 print(json.dumps(dict(all_passed=True,genuine_dimension=2,rA_action='-11/16',seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
