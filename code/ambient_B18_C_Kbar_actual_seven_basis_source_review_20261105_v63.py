#!/usr/bin/env python3
"""Independent exact seven barclasses, four formal edge entries retained."""
import itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
from ambient_B18_C_connecting_kernel3_split_review_20261105_v63 import c,v,add,sc,mm
ROOT=pathlib.Path(__file__).resolve().parents[1];PRI='AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def numeric_mm(A,B):
 out=[[Q(0)]*len(B[0])for _ in A]
 for j in reversed(range(len(B[0]))):
  for i in reversed(range(len(A))):
   for k in reversed(range(len(B))):out[i][j]+=A[i][k]*B[k][j]
 return out
def main():
 t=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));s=sha(p)
  if h:assert h==s,(k,h,s)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=s
  return json.loads(p.read_text())if parse else None
 pr=pin(ROOT/'reports'/PRI/'report.json','2bdf9a211b43dd0e27ad2f014da5f702af1b102f2ba6e42b6c54188dad08bb85')
 for p,h in pr['input_sha256'].items():
  pin(ROOT/p,h,False)
  if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
 pin(ROOT/'code/ambient_B18_C_Kbar_actual_seven_basis_source_20261105_v63.py',pr['code_sha256'],False)
 en=next(iter(pr['evidence_sha256']));D=pin(ROOT/'reports'/PRI/en,pr['evidence_sha256'][en])
 def src(s):return ROOT/next(p for p in pins if p.endswith(s))
 val=pin(src('AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63/report.json'),'29133a93c9c3f24d1bf228cacd4f49f8fd5233c255a3b85a159598da21cc297d')
 vpeer=pin(src('AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63_REVIEW/report.json'),'bffdb3cb50ea194531b9e237da68e2612872aa0fbdc4edb32f7f9b2db8a5b058');assert vpeer['all_passed']
 known=field(src('actual-M7-stable-surjection-splitting-and-kernel-three-source-basis.json'),'known_values')
 p1,b1,d,a0,a4=map(Q,[known[k]for k in['head1_rg_p1','head1_k2_b1','cyclic_top_d','h0_k1_a0','h4_k1_a4']])
 p2,b2=map(Q,val['head29_A_B_rg_k2'][2:]);p3,b3=map(Q,val['head30_A_B_rg_k2'][2:])
 # Independent coordinates solve each2x2 triangular edge cancellation.
 aU=b3/b1;cU=(p3-aU*p1)/d;aS=(b2-b3/4)/b1;cS=(p2-p3/4-aS*p1)/d
 K=[[Q(0)]*7 for _ in range(11)]
 for j,i in enumerate([6,8,9,10,5]):K[i][j]=1
 K[3][5]=1;K[1][5]=-aU;K[7][5]=-cU
 K[2][6]=1;K[3][6]=-Q(1,4);K[1][6]=-aS;K[7][6]=-cS
 assert K==[[Q(q)for q in row]for row in D['actual_Kbar11_by7']]
 correction=dict(zip(['a3=b3/b1','c3=(p3-a3*p1)/d','aS=(b2-b3/4)/b1','cS=(p2-p3/4-aS*p1)/d'],[aU,cU,aS,cS]))
 assert correction=={k:Q(q)for k,q in D['actual_correction_coefficients'].items()}
 # Honest all-order cycle source, not only the constant56x5 origins.
 five=src('actual-five-head-formal-RegDivell-operator-and-source-binding.json')
 froot=field(five,'matrix_roots');fnodes=field(five,'matrix_operator_nodes');fshape=field(five,'matrix_operator_shapes')
 assert froot['five_formal_cycles_ELL0']==6 and fshape[6]==[56,5]
 assert fnodes[0]==['actual_external_ELL0_matrix','old_unit','typed_m24_by56','ELL0_S_F9']
 assert fnodes[1]==['actual_external_ELL0_matrix','old_unit','typed_n56_by24','ELL0_S_F9']
 assert fnodes[3:7]==[['identity_matrix',56],['matrix_product',1,0],['matrix_difference',3,4],['matrix_product',5,2]]
 E=[[Q(q)for q in row]for row in fnodes[2][1]]
 assert [[i for i in range(56)if E[i][j]]for j in range(5)]==[[22],[23],[29],[30],[55]]
 actual=ROOT/'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/exact-Loewy11-and-source-basis-corollary.json'
 pin(actual,'b38faf02bcee171bffd26dddd4a996aae44f21a8ea7c85d9b0351c8eb150fa27',False)
 basis=field(actual,'actual_filtered_monomial_basis11');assert field(actual,'exact_Loewy_HF')==[5,3,2,1]
 monmap=[[[]for _ in range(11)]for _ in range(5)]
 for j,x in enumerate(basis):monmap[x['target']][j]=[[x['monomial_original8'],'1']]
 assert monmap==D['actual_old11_by_fivehead_polynomial_map5_by11']
 owned=count([pr,D,val,vpeer,known,fnodes,fshape,E,basis,monmap,K,pins])+20000
 assert owned<150000 and raw+1000000<32<<20
 save(out/'positive-before-formal-edge-and-fixed-coordinate-inverse.json',{'pure_SOURCE_no_responses':True,'four_unknown_head_entries_formal':True,'whole_numeric_upper':owned,'no56_scalar_expansion':True})
 edge=[[{}for _ in range(11)]for _ in range(4)]
 edge[0][4]=c(d);edge[2][0]=c(a0);edge[2][4]=c(a4);edge[1][7]=c(d)
 for j,p,b in[(1,p1,b1),(2,p2,b2),(3,p3,b3)]:edge[1][j]=c(p);edge[3][j]=c(b)
 for row,col,name in[(1,0,'p0'),(3,0,'b0'),(1,4,'p4'),(3,4,'b4')]:edge[row][col]=v(name)
 assert mm(edge,[[c(q)for q in row]for row in K])==[[{}for _ in range(7)]for _ in range(4)]
 assert all(K[i][j]==0 for i in[0,4]for j in range(7))
 complement=[4,7,0,1];B=[row+[Q(i==j)for j in complement]for i,row in enumerate(K)]
 # Solve the11 coordinate equations, not a Gaussian inverse/rank search.
 inv=[[Q(0)]*11 for _ in range(11)]
 for j,i in enumerate([6,8,9,10,5]):inv[j][i]=1
 inv[5][3]=1;inv[5][2]=Q(1,4);inv[6][2]=1
 inv[7][4]=1;inv[8][7]=1;inv[9][0]=1;inv[10][1]=1
 inv[8][3]=cU;inv[8][2]=cU/4+cS;inv[10][3]=aU;inv[10][2]=aU/4+aS
 assert B==[[Q(q)for q in row]for row in D['whole11_basis_K_then_S']]
 assert inv==[[Q(q)for q in row]for row in D['explicit_whole11_inverse']]
 identity=[[Q(i==j)for j in range(11)]for i in range(11)]
 assert numeric_mm(inv,B)==identity==numeric_mm(B,inv)
 Sedge=[[row[j]for j in complement]for row in edge]
 # The forced nonzero positions already give det d^2*a0*b1, independent
 # of p0,b0,p4,b4. Symbolic elimination is unnecessary.
 assert Sedge[0]==[c(d),{},{},{}]and Sedge[1][1]==c(d)
 assert Sedge[2][2:]==[c(a0),{}]and Sedge[3][3]==c(b1)
 det=d*d*a0*b1;assert det==Q(11979,2260608800)
 nodes=D['source_operator_nodes'];shapes=D['source_operator_shapes']
 assert nodes[0]==['external_actual_sigma56_by5','five',6,'ELL0_S_F9']
 assert nodes[1]==['actual_Rbar_polynomial_matrix',monmap,'original8_z0_to_z6_r']and nodes[2]==['constant_QQ_matrix',D['actual_Kbar11_by7']]
 assert nodes[3:]==[['matrix_product',0,1],['matrix_product',3,2]]
 assert shapes==[[56,5],[5,11],[11,7],[56,11],[56,7]]and D['actual_barcycles56_by7_root']==4
 for i in[3,4]:
  left,right=nodes[i][1:];assert shapes[left][1]==shapes[right][0]and shapes[i]==[shapes[left][0],shapes[right][1]]
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'actual_four_head_values_and_known_old_p1_b1_d_a0_a4_bound_to_independent_peers':True,
  'all28_symbolic_edge_K_identities_with_p0_b0_p4_b4_formal':True,'both121_coordinate_inverse_products_checked':True,
  'seven_actual_QQ_basis_labels':D['labels'],'same_four_actual_correction_coefficients':{k:str(q)for k,q in correction.items()},
  'complement_edge_determinant_nonzero_independent_uncomputed_heads':str(det),
  'old11_actual_filtered_QQ_basis_and_exact5x11_Rbar_monomial_map_verified':True,
  'actual_full_formal_sigma56_by5_root6_I_minus_nm_E_verified':True,
  'typed56x7_composition_source_shapes_and_root4_verified_without_expansion':True,
  'kernel_and_descent_proof':'The accepted genuine connecting map is Rbar-linear on Told11 with kernel K/ellK7. All seven displayed columns are killed by the full edge, including its four still-formal entries. The explicit11 QQ change of basis proves independence, so seven columns span the full genuine kernel. Compose the actual formal cycle sigma=(I-nbar*mbar)E with the old11 monomial map and K. Rbar-linearity and complete768 relation descent give actual mod-ell cycles, without a naive homology tensor shortcut.',
  'critical_splitting_scope':'B=[K,S] is a QQ vector-space basis of Told11. It does NOT assert an R0-module splitting of originalUCT; S is stable only in the quotientM7, and originalUCT remains nonsplit by the independent e4abd proof.',
  'not_seven_free_R0_generators_or_full_K_lifts_or_K_action_matrix':True,
  'no_new_response_source_gauge_or_jet_and_no_originalH0_nil_inference':True,'owned_upper':owned,'RAW':raw}
 ep=out/'independent-actual-seven-barbasis-formal-edge-and-typed-cycle-source.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'actual_Kbar_QQ_basis7_and_typed56_cycles_PASS':True,'R0_UCT_split_claimed':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'seconds':report['seconds']}))
if __name__=='__main__':main()
