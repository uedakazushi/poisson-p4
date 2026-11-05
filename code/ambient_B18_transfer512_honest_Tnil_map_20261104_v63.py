#!/usr/bin/env python3
"""Actual512 original-column homotopies into honest128 and minimal5.
Source-gated, no producer import/source solve/jet or large scalar expansion.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_honest_nil_T_first_values_review_20261104_v63 as Q
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TRANSFER512_HONEST_TNIL_MAP_20261104_v63';OUT=ROOT/'reports'/ID
class Mat:
 def __init__(self):self.nodes=[];self.shapes=[]
 def add(self,op,shape,*args,**kw):
  assert all(i<len(self.nodes)for i in args);i=len(self.nodes);self.nodes.append(dict(id=i,op=op,args=list(args),shape=shape,**kw));self.shapes.append(shape);return i
 def ext(self,ns,node,shape):return self.add('actual_frozen_matrix_operator',shape,namespace=ns,node=node)
 def mul(self,*args):
  assert all(self.shapes[a][1]==self.shapes[b][0]for a,b in zip(args,args[1:]));return self.add('matrix_product',[self.shapes[args[0]][0],self.shapes[args[-1]][1]],*args)
 def block(self,rr):
  ss=[[self.shapes[j]for j in row]for row in rr];assert all(len({s[0]for s in row})==1 for row in ss)and all(len({row[j][1]for row in ss})==1 for j in range(len(ss[0])));return self.add('block_matrix',[sum(row[0][0]for row in ss),sum(s[1]for s in ss[0])],*(j for row in rr for j in row),blocks=rr)
 def eye(self,n):return self.add('identity',[n,n])
 def zero(self,n,m):return self.add('zero',[n,m])
 def neg(self,a):return self.add('negative',self.shapes[a],a)
 def sub(self,a,b):assert self.shapes[a]==self.shapes[b];return self.add('difference',self.shapes[a],a,b)
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 return x
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 stW='AMBIENT_B18_CONDUCTOR_TRANSFER512_NIL_CYCLES_20261104_v63';stH='AMBIENT_B18_HONEST_NIL_T_OPERATOR_20261104_v63';stM='AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63'
 for st in[stW,stH,stM,'AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63']:
  assert load(st+'_REVIEW')['all_passed']
 W=load(stW,'actual-transfer512-H1-cycles-typed-DAG.json');H=load(stH,'actual-honest-Tnil128-972-formal-operator.json');Min=load(stM,'actual-minimal-Tnil5-849-and-d3-linear-operator.json')
 R=load('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json');cycles=W['actual_512_transfer_cycles'];assert len(cycles)==512
 for c in cycles:assert len(c['literal_H1_cycle44'])==44 and all(W['new_node_types'][i]=='R0'for i in c['literal_H1_cycle44'])
 assert H['actual_S_frame']==Min['actual_S_frame']==['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA']
 assert H['shared_homotopy_template_v1']['not_S_linear']and Min['matrices']['direct_unit_section_n']==8 and Min['matrices']['full_target_projection5_by128']==59
 # Persist a source gate BEFORE constructing any map or new value.
 typed_upper=len(W['new_typed_scalar_nodes'])+len(H['scalar_function_nodes'])+len(H['matrix_operator_nodes'])+len(H['vector_formal_function_operator_nodes'])+len(Min['operator_nodes'])+512*10+512*44+100
 origin_cache_upper=sum(1 if s=='R0'else 4 for s in W['new_node_types'])
 budget=typed_upper+origin_cache_upper+6000+512*13
 assert budget<150000
 save('preflight.json',dict(phase='BEFORE map/value construction',actual_source_cycles512=True,all22528_source_slots_are_original_R0_function_roots=True,shared_frame=H['actual_S_frame'],new_vector_node_bound=5120,new_matrix_node_bound=100,including_frozen_types_source_slots_upper=typed_upper,origin_R0_C4_value_cache_upper=origin_cache_upper,other_cache_and_output_allowance=6000+512*13,total_bound=budget,output_byte_bound=12_000_000,no_new_source_or_jet=True,origin_gate='all W0 zero plus W20 firsts from β0=0 and original C4 p origins; no unknown source callbacks'))
 save('phase.json',dict(phase='sourcegate PASS; building original-column operator'))
 M=Mat();wf=M.add('actual_transfer512_R0_columns_restrict_ell0',[44,512],namespace='transfer',json_pointer='/actual_512_transfer_cycles',field='literal_H1_cycle44',restriction=H['ell0_pullback'])
 m=M.ext('honest',H['matrices']['new_m24_by44'],[24,44]);n=M.ext('honest',H['matrices']['normalized_n44_by24'],[44,24]);d2=M.ext('honest',H['matrices']['d2_new'],[24,84]);phi=M.ext('honest',H['matrices']['Phi44_by128'],[44,128]);nd=M.ext('minimal',Min['matrices']['direct_unit_section_n'],[44,24]);pi=M.ext('minimal',Min['matrices']['full_target_projection5_by128'],[5,128]);mw=M.mul(m,wf)
 vec=[]
 def v(op,dim,*args,**kw):
  assert all(i<len(vec)for i in args);i=len(vec);vec.append(dict(id=i,op=op,dimension=dim,args=list(args),**kw));return i
 columns=[]
 for j in range(512):
  z=v('ORIGINAL_actual_matrix_column_functions',44,matrix=wf,column=j);mv=v('ORIGINAL_actual_matrix_column_functions',24,matrix=mw,column=j);ss=[]
  for sector in range(3):
   part=v('contiguous_vector_slice',8,mv,start=8*sector,length=8)
   ss.append(v('ORIGINAL_complete_Koszul_Hadamard_Euler_lift',28,part,power=[1,1,2][sector],template_namespace='honest',template='shared_homotopy_template_v1',proof='d1*m*W=Vgamma*W=0; same512 actual transfer cycle identity',before_any_normal_form_coefficient_factoring=True))
  b=v('ordered_vector_concatenation',84,*ss);columns.append(v('ordered_vector_concatenation',128,z,b))
 assert len(vec)==5120
 w128=M.add('matrix_from_ORIGINAL_column_vector_operators',[128,512],columns=columns,extend_S_linearly_only_after_columnwise_contraction=True)
 bottom=M.block([[M.zero(84,44),M.eye(84)]])
 J=M.block([[phi],[bottom]])
 phid=M.block([[M.sub(M.eye(44),M.mul(nd,m)),M.mul(nd,d2)]])
 Ji=M.block([[phid],[bottom]])
 wd=M.mul(J,w128);w5=M.mul(pi,wd);w5simple=M.mul(pi,w128)
 sources={name:dict(path=f'reports/{st}/{fn}',sha256=pins[f'reports/{st}/{fn}'])for name,st,fn in[('transfer',stW,'actual-transfer512-H1-cycles-typed-DAG.json'),('honest',stH,'actual-honest-Tnil128-972-formal-operator.json'),('minimal',stM,'actual-minimal-Tnil5-849-and-d3-linear-operator.json')]}
 op=dict(source_artifacts=sources,matrix_nodes=M.nodes,vector_nodes=vec,matrices={'actual_W44':wf,'mW24':mw,'honest_W128':w128,'old128_to_direct128':J,'direct128_to_old128':Ji,'direct_Phi':phid,'direct_W128':wd,'minimal_W5':w5,'same_minimal_W5_without_J_on_these_columns':w5simple},
  original512_source_labels=[dict(index=c['index'],i=c['CA_generator_index'],j=c['CB_raw_index'],nil=c['normalized_nil_generator_index'],basis=c['normalizer_basis_index'])for c in cycles],
  exact_identity='d2*H512=m*W44; Phi_n*[W44;H512]=Phi_ndirect*[W44;H512]=W44',
  section_change_proof='J_n(s,b)=(Phi_n(s,b),b). Both ways map actual(X,h),section,d3 relations; J_ndirect*J_n differs from identity by(n,0), and reverse by(ndirect,0). For these512 lifts J_n fixes each128 vector.',
  actual_map='W→N/ellN→Tnil=coker128; then direct128 canonical section change and minimal5 projection',
  W_mod_ellW_injection_NOT_assumed=True,full_image='(W+ellN)/ellN',no_contraction_Slinearity_assumed=True,no_new_scalar_source_or_jet=True)
 save('actual-transfer512-honest128-minimal5-operator.json',op)
 save('phase.json',dict(phase='typed actual map SAVED before optional cache-origin arithmetic'))
 # Cache-only ORIGIN gate: independent C4 algebra, not a new source jet.
 old=load('AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63','actual-UP-QB-640-typed-operator.json');atom=load('AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63','phase1-actual-atom-firsts.json');P0=load('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase4-native-source-and-R12.json');G0=load('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase5-native-adapter-firsts.json');GL=load('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63','actual-E4-nil-two-generators-and-relations.json')
 X=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase2-actual-radical-X-firsts.json');m0dat=load('AMBIENT_B18_HONEST_NIL_T_FIRST_VALUES_20261104_v63','phase1-actual-m0-and-section.json');lin=load('AMBIENT_B18_HONEST_NIL_T_FIRST_LINEAR_20261104_v63','phase3-all849-minimal-first-relations-and-rank.json')
 atoms={int(i):Q.pairs(v)for i,v in atom['all79_atom_values']};ops=old['operator_scalar_nodes'];oc={}
 def oval(i):
  if i in oc:return oc[i]
  a=ops[i]
  if a[0]=='c':z=F(a[1],a[2])
  elif a[0]=='actual_C0_scalar_function':z=atoms[a[1]].get(0,F(0))
  elif a[0]=='sum':z=sum((F(c,d)*oval(j)for j,c,d in a[1]),F(0))
  elif a[0]=='product':
   z=F(a[1],a[2])
   for j in a[3]:z*=oval(j)
  else:raise AssertionError(a)
  oc[i]=z;return z
 C4=Q.FiniteC([[0]*8],[[0],[0]]);core=[{j:F(c)for j,c in enumerate(P0['kappa9'][i][:4])if F(c)}for i in[0,1,2,6]]
 assert all(not c.get(0,F(0))and not c.get(3,F(0))for c in core)
 var=[{j:F(c)for j,c in enumerate(v[:4])if F(c)}for v in G0['G8'][:5]]+[C4.polynomial(GL['forward_GL4'][s],core)for s in ['h','c','r','l']]
 d=C4.polynomial(R['d_polynomial'],var);qs=[C4.mul(d,C4.polynomial(p,var))for p in R['R9_generators_polynomials9']]
 qcomp={}
 for j,q in enumerate(qs):
  for gamma in range(4):
   value=C4.mul(q,{gamma:F(1)})
   for k,root in enumerate(R['actual_Qrad4_by36_columns'][4*j+gamma]):
    val=value.get(k,F(0));assert root not in qcomp or qcomp[root]==val;qcomp[root]=val
 A0={}
 for i,row in enumerate(R['native_new_A12_by9_nodes']):
  for j,root in enumerate(row):
   v={k:F(c)for k,c in enumerate(X['new_A12_by9_C36'][i][j][:4])if F(c)};assert root not in A0 or A0[root]==v;A0[root]=v
 # Constant Pi0 cleanup bound to the exact18 unit rows and literalthreewedge projection.
 assert Q.sparse(m0dat['m0_24_by44'][7])=={20:F(5,2)}
 assert Min['d3_only_linear']['constant_projection3_by84']==[[str(int(i==j))for i in range(84)]for j in[6,24,26]]
 xp=X['X0_pivot_columns'];E=m0dat['selected_source_columns24'];C=[i for i in range(44)if i not in E];I=[C[i]for i in Min['fixed_minor_data']['compatible_X0_top18_rows_in_C20']];Ji0=[Q.sparse(row)for row in Min['fixed_minor_data']['inverse']]
 Xv=[{int(i):Q.pairs(v)for i,v in col}for col in X['all780_source_columns44_sparse_dual9']]
 hi=[[F(5,2)*Xv[j].get(20,{}).get(i+1,F(0))for j in xp]for i in[0,4,5]];assert hi==[[F(c)for c in row]for row in lin['all54_actual_pivot_h0_readback']]
 pi0=[]
 for target in[42,43]:
  row={target:F(1)}
  for k,j in enumerate(xp):
   for u,a in Ji0[k].items():Q.put(row,I[u],-Xv[j].get(target,{}).get(0,F(0))*a)
  pi0.append(row)
 for wedge,h in zip([6,24,26],hi):
  row={44+wedge:F(1)}
  for k,b in enumerate(h):
   for u,a in Ji0[k].items():Q.put(row,I[u],-b*a)
  pi0.append(row)
 # Free input/large first arrays are not retained during the25088-node origin pass.
 del Xv,X,lin,m0dat,Min,H,vec,M,op
 nodes=W['new_typed_scalar_nodes'];types=W['new_node_types'];values=[]
 for i,a in enumerate(nodes):
  typ=a[0]
  if typ=='R0_QQ':v=F(a[1],a[2])
  elif typ=='R0_actual_external':
   ns,j=a[1:];v=atoms[j].get(0,F(0))if ns=='C0'else oval(j)if ns=='chain'else qcomp[j]
  elif typ=='C_actual_external':assert a[1]=='radical';v=A0[a[2]]
  elif typ=='R0_neg':v=-values[a[1]]
  elif typ=='R0_sum':v=sum((values[j]for j in a[1]),F(0))
  elif typ=='R0_product':
   v=F(1)
   for j in a[1]:v*=values[j]
  elif typ=='C_embed_R0':v=C4.c(values[a[1]])
  elif typ=='C_basis':v={a[1]:F(1)}
  elif typ=='C_sum':v=Q.plus(*(values[j]for j in a[1]))
  elif typ=='C_product':
   v=C4.c(1)
   for j in a[1]:v=C4.mul(v,values[j])
  elif typ=='R0_component_of_C':v=values[a[1]].get(a[2],F(0))
  else:raise AssertionError(typ)
  assert (types[i]=='R0')==isinstance(v,F);values.append(v)
 # Each derivative is computed by the actual β*A5*t term. All omitted
 # derivatives multiply β0 or β0*b0, which were just checked exactly0.
 metadata=[];zero5=[]
 for c in cycles:
  assert all(values[i]==0 for i in c['literal_H1_cycle44'])
  ii,jj,alpha,gamma=[c[k]for k in ['CA_generator_index','CB_raw_index','normalized_nil_generator_index','normalizer_basis_index']];group=[5,7][alpha]
  proots=c['p_four_actual_R0_components'];assert [nodes[i]for i in proots]==[['R0_actual_external','radical',j]for j in R['actual_Qrad4_by36_columns'][4*group+gamma]];pp=[values[i]for i in proots];assert pp[2]==0
  scalar=C4.mul(A0[R['native_new_A12_by9_nodes'][jj][5]],Q.plus(C4.c(pp[2]),Q.scaled({1:F(1)},pp[3]))).get(0,F(0));assert scalar==0
  beta=atoms[W['source_artifacts']['C0']['path']and [7148,7199,7454,7505,7688,7943,8045,505938][ii]];assert not beta.get(0,F(0))
  row20first=[beta.get(k+1,F(0))*scalar for k in range(8)];assert not any(row20first)
  # Only final wedge6/24/26 origins survive Pi0 when z0=0;
  # original homogeneous contraction gives h_i7(0)=∂zi(mz)_7.
  z5=[F(0),F(0)]+[F(5,2)*row20first[i]for i in[0,4,5]];assert not any(z5);zero5.append(z5)
  metadata.append(dict(index=c['index'],all44_actual_source_origins_zero=True,pB_origin=pp[2],pAB_origin=pp[3],actual_row20_first8=row20first,minimal_origin5=z5))
 assert len(zero5)==512
 save('actual512-origin-values-and-Pi0.json',dict(all512=metadata,constant_minimal_projection5_by128_sparse=[[[i,c]for i,c in sorted(row.items())]for row in pi0],minimal_origin_rank=0,all512_columns_zero=True,all25088_actual_source_origins_evaluated=True,C4_nil_generator_origin=[Q.sparse([qs[g].get(k,F(0))for k in range(4)])for g in[5,7]],m0_row7=[20,'5/2'],no_m1_or_new_source_callback_needed_since_all_W0_zero=True))
 save('proof.json',dict(actual_map128_then5=True,all512_ORIGINAL_formal_function_columns_lifted_before_anyNF=True,regularity='d1*mW=0, regularr^n implies last entry vanishes at z0; orderedHadamard and positiveEuler d+2 give d2H=mW',
  old_to_direct_section_comparison='J_n(s,b)=(Phi_n(s,b),b), inverse on cokernels J_ndirect. JndirectJn−I=(n,0)*(−m,d2); reverse difference=(ndirect,0)*(−m,d2). The512 chosen vectors are fixed because Phi_n(z,Hmz)=z.',
  W_image='(W+ellN)/ellN inside honest Tnil',stronger_sufficient_test='length(Tnil/imW)=2 implies N/ellN=imW, N=W by D-Nakayama, and N³=0',
  actual_origin_result='ALL512 minimal origin columns are0. Therefore imW⊂m_S*Tnil and mu(Tnil/imW)=5, so its length is at least5. The sufficient length2 test FAILS.',
  no_false_negative='This does not prove N³ nonzero or refute P=psi(W). K/ellK and the length9 quotient remain unconstructed;2..4 invisibleTor directions retained.',
  no_all_order_Wzero_claim=True,no_new_source_jet_gauge_ranksearch=True))
 for p in[Path(__file__),Path(Q.__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 report=dict(id=ID,all_passed=True,frozen=True,actual_complete_map_shapes=[[44,512],[128,512],[5,512]],typed_vector_nodes=5120,typed_matrix_nodes=len(json.loads((OUT/'actual-transfer512-honest128-minimal5-operator.json').read_text())['matrix_nodes']),all512_columnwise_original_homotopies=True,old_direct128_section_change_explicit=True,all512_minimal_origin_columns_zero=True,origin_rank=0,quotient_minimal_generator_count=5,strong_length2_sufficient_test_FAILED=True,nil_cube_still_OPEN=True,no_new_source_jet_scalar_expansion=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-t)
 save('phase.json',dict(phase='FINAL actual typed map and cache-origin obstruction'));report['evidence_sha256']['phase.json']=sha(OUT/'phase.json');save('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
