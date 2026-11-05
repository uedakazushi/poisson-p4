#!/usr/bin/env python3
"""Independent fixed-source factored R review; no dense inverse/actual alpha."""
import hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field,nc
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='8741d10b4cc1253c086e7f3c89163cc3f87caac2374a62f17402dcbcdf2fb494'
maximum=actual=preactual=0
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def sp(a):return {int(i):Q(q) for i,q in a if Q(q)}
def mark(a):
 global maximum
 maximum=max(maximum,a[0].bit_length(),a[1].bit_length());assert maximum<=134
 return a
def ex(q):q=Q(q);return mark((q.denominator,abs(q.numerator))) if q else (1,0)
def ba(a,b):
 if not a[1]:return b
 if not b[1]:return a
 D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
def bm(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
def note(q):
 global actual
 actual=max(actual,abs(q.numerator).bit_length(),q.denominator.bit_length());assert actual<=134
 return q
def qa(a,b):
 global preactual
 if not a:return b
 if not b:return a
 n=max(max(abs(a.numerator).bit_length()+b.denominator.bit_length(),abs(b.numerator).bit_length()+a.denominator.bit_length())+1,a.denominator.bit_length()+b.denominator.bit_length())
 preactual=max(preactual,n);assert n<=403
 return note(a+b)
def qm(a,b):
 global preactual
 if not a or not b:return Q(0)
 n=max(abs(a.numerator).bit_length()+abs(b.numerator).bit_length(),a.denominator.bit_length()+b.denominator.bit_length())
 preactual=max(preactual,n);assert n<=403
 return note(a*b)
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=True);assert not (out/'report.json').exists();pins={};RAW=0
 def pin(p,h=None,select=None):
  nonlocal RAW
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  k=str(p.relative_to(ROOT))
  if k not in pins:RAW+=p.stat().st_size
  pins[k]=v
  if select:return field(p,select)
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 for p,h in pr['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 E=pin(ROOT/'reports'/PRI/'actual-factored918-rightsolve-source-circuit-and-all-generator-identities.json',
  '878d1b39b00ffab99314d595d579e4bd839a30098f93036b8e54ae53f68db3c5')
 assert sha(ROOT/'code/ambient_B18_ell2_T2_full918_rightsolve_source_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
 D=pin(ROOT/'reports/AMBIENT_B18_ELL2_T2_G2_STRAND_SOURCE_PROFILE_20261105_v63/literal-active-G2-linear-strand-matrix-singleton-units-and-residual-source.json',
  '2230420db5f89f66290a7494106353c2d30023f3b715d2f2e01074aa0cb1da8d')
 I=pin(ROOT/'reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63/pruned46-by41-integer-rank-source-and-live-bound.json',
  'efc72ae140207d2d9eedfa080328932d27a6c34c3a021cf8ef84ad7c9b3e15a2')
 cert=pin(ROOT/'reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63/fixed46-by41-integer-Schur-rank-and-full-strand-kernel-recipes.json',
  '956d3b41691b1088fcf91779b3be32082d4561a75e12a67c210e72df08d74493')
 factor=pin(ROOT/'reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63_REVIEW/independent-fixed-minor-all41-factorizations-and-full-strand-proof.json',
  'ff377febdf22dddf025adfdd83912fe4f60fc940c61936114a1c3f1190f79111')
 assert factor['full918_strand_kernel_equals_known_G3linear242'] and factor['active_rank']==190
 for p,h in [('reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63_REVIEW/report.json','ed6954d48503e2220e1d7ab29526a64ead78387ab00e1012f70b454172ff20e3'),('reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63_REVIEW/report.json','fbc71ad533749bf42d4636c627300d5bd67b9265e6aec2bb9f98b07104c0b4a9')]:assert pin(ROOT/p,h)['all_passed']
 Fraw=pin(ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
  '4a298af981efd32b617db909ec2cfc9dea8ef062d7b545444a514896c946f850',select='full_G2_linear234_by102_columns')
 G1=pin(ROOT/'reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_INPUT_GATE_20261105_v63/actual-full-G1-and186-mixed-adjoint-source-input.json',
  'd10780b841d2a063a317bdaa885d5d652a2ed0fe7aa1b99ad97d5c34b6b033ac',select='exact_full_G1_polynomial_columns')
 F=list(map(sp,Fraw));A=list(map(sp,D['literal_matrix_columns']));assert Fraw==E['actual_full_G2linear_columns']
 assert len(F)==102 and len(A)==252
 nodes=E['scalar_linear_nodes'];outputs=E['full_source_outputs918'];active=E['active_source_outputs252'];block=E['active_unit_block28'];steps=E['active_forward_and_reverse162_source_steps']
 fixtures=list(map(sp,E['all918_source_generator_R_F_columns']));support=sum(map(len,fixtures));assert support<=1711
 va=[0,1,4,5,7,8];vz=[2,3,6];pairs=list(itertools.combinations_with_replacement(range(9),2));pairs.remove((8,8));mi={m:i for i,m in enumerate(pairs)}
 assert E['quadratic_pair_order44']==[list(p) for p in pairs]
 def slot(t,a,b):return 44*t+mi[tuple(sorted((a,b)))]
 assert len(nodes)==2349 and len(outputs)==918 and len(steps)==162
 def key(n):return (n[0],tuple(n[1])) if n[0]=='sum' else tuple(n)
 lookup={key(n):i for i,n in enumerate(nodes)};assert len(lookup)==len(nodes)
 for i,n in enumerate(nodes):
  if n[0]=='scale':assert n[2]<i
  elif n[0]=='sum':assert all(j<i for j in n[1])
  elif n[0]=='input':assert 0<=n[1]<1144
  else:assert n[0]in['constant','unit28_component']
 def sn(q,x):
  if not q or not x:return 0
  if q==1:return x
  return lookup[('scale',str(q),x)]
 def total(seq):
  seq=[j for j in seq if j]
  return 0 if not seq else seq[0] if len(seq)==1 else lookup[('sum',tuple(seq))]
 inp=[lookup[('input',j)] for j in range(1144)]
 m6=[p for p in itertools.combinations_with_replacement(range(6),2) if p!=(5,5)]
 rhs=[inp[slot(t,va[a],va[b])] for t in range(17) for a,b in m6]
 liveR=set(range(340));liveC=set(range(252));xnodes={}
 for old,st in zip(D['fixed_singleton_unit_pivots'],steps):
  row,col,q=old['row340'],old['column252'],Q(old['literal_nonzero_unit'])
  assert st['row340']==row and st['column252']==col and st['pivot']==str(q) and st['kind']==old['singleton_kind']
  assert A[col][row]==q and q and st['saved_rhs_node']==rhs[row]
  rem=[[j,str(A[j][row])] for j in sorted(liveC-{col}) if row in A[j]]
  assert st['remaining_row_coefficients']==rem
  if st['kind']=='row_singleton':
   assert not rem;xnodes[col]=sn(1/q,rhs[row])
   for r in sorted(liveR-{row}):
    if r in A[col]:rhs[r]=total([rhs[r],sn(-A[col][r],xnodes[col])])
  else:assert not(set(A[col])&(liveR-{row}))
  liveR.remove(row);liveC.remove(col)
 sr=cert['selected_residual_rows46'];sc=cert['selected_residual_columns41'];R=I['original340_row_indices46'];C=I['original252_column_indices41'];selected=[C[j] for j in sc]
 U=block['matrix_integer28_by28'];assert U==[[I['residual_integer_matrix46_by41'][i][j] for j in sc] for i in sr]
 assert block['determinant_literal']==factor['selected_integer_determinant'] and int(block['determinant_literal'])
 assert block['input_nodes28']==[sn(Q(I['row_integer_scaling_D46'][i]),rhs[R[i]]) for i in sr]
 assert max(block['input_nodes28'])<min(i for i,n in enumerate(nodes) if n[0]=='unit28_component')
 for k,j in enumerate(selected):xnodes[j]=lookup[('unit28_component',0,k)]
 for j in liveC-set(selected):xnodes[j]=0
 for st in reversed(steps):
  if st['kind']=='column_singleton':xnodes[st['column252']]=sn(1/Q(st['pivot']),total([st['saved_rhs_node']]+[sn(-Q(q),xnodes[j]) for j,q in st['remaining_row_coefficients']]))
 assert active==[xnodes[j] for j in range(252)]
 assert all(outputs[9*i+va[a]]==active[6*i+a] for i in range(42) for a in range(6))
 for i,n in enumerate(nodes):
  if n[0]=='unit28_component':assert n[1]==0 and 0<=n[2]<28
 gauge=[j for j,n in enumerate(outputs) if n==0]
 assert gauge==E['explicit_zero_gauge_source_slots242'] and len(gauge)==242 and sum(bool(j) for j in outputs)==676
 def fixture(j):
  i,v=divmod(j,9);return {slot(t,u,v):q for row,q in F[i].items() for t,u in [divmod(row,9)] if (u,v)!=(8,8)}
 selected_full=[9*(j//6)+va[j%6] for j in selected]
 own=nc([E,D,I,cert,factor,Fraw,G1,pr])+2*(sum(map(len,A))+sum(map(len,F))+support)+4*len(nodes)+4*1144+2*918+20000
 assert own<150000 and RAW<32<<20
 def evaluate(j,positive):
  q=fixture(j);y=fixtures[j];zero=(1,0) if positive else Q(0);vals=[]
  for n in nodes:
   op=n[0]
   if op=='constant':assert n[1]=='0';v=zero
   elif op=='input':v=ex(q.get(n[1],0)) if positive else q.get(n[1],Q(0))
   elif op=='unit28_component':v=ex(y.get(selected_full[n[2]],0)) if positive else y.get(selected_full[n[2]],Q(0))
   elif op=='scale':v=bm(ex(n[1]),vals[n[2]]) if positive else qm(Q(n[1]),vals[n[2]])
   else:
    v=zero
    for k in reversed(n[1]):v=ba(v,vals[k]) if positive else qa(v,vals[k])
   vals.append(v)
  for r in range(28):
   v=zero
   for col,c in enumerate(U[r]):
    if positive:v=ba(v,bm(ex(c),ex(y.get(selected_full[col],0))))
    else:v=qa(v,qm(Q(c),y.get(selected_full[col],Q(0))))
   if positive:ba(v,vals[block['input_nodes28'][r]])
   else:assert v==vals[block['input_nodes28'][r]]
  result={k:vals[node] for k,node in enumerate(outputs) if (vals[node][1] if positive else vals[node])}
  if not positive:assert result==y
  target={}
  for k,c in result.items():
   src,var=divmod(k,9)
   for row,a in F[src].items():
    t,u=divmod(row,9)
    if(u,var)==(8,8):continue
    z=slot(t,u,var)
    if positive:target[z]=ba(target.get(z,zero),bm(ex(a),c))
    else:
     target[z]=qa(target.get(z,Q(0)),qm(a,c))
     if not target[z]:del target[z]
  if positive:
   for z,a in q.items():ba(target.get(z,zero),ex(a))
  else:assert target==q
  return target
 # The supplied source witnesses are known literal inputs. Bound all DAG,
 # block and product verifications before any signed QQ fixture replay.
 for j in range(918):evaluate(j,True)
 assert maximum<=pr['positive_known_fixture_prior_bits']==134
 save(out/'phase-positive-before-signed-known-fixtures.json',{'known_fixed_fixture_positive_bound':maximum,
  'source_prior134_and403_valid':True,'whole_own_numeric_upper':own,'producer_whole66378_bound_scope':'accepted producer allocation; independent reader has separate owned envelope',
  'no_dense_inverse_or_actual_alpha_K_rho_b1_Z2':True,'RAW_read':RAW})
 hashes=[]
 for j in range(918):
  target=evaluate(j,False);h=hashlib.sha256(json.dumps(sorted(target.items()),default=str,separators=(',',':')).encode()).hexdigest()
  assert h==E['all918_F_R_F_identity_fingerprints'][j];hashes.append(h)
 for j in range(252):
  src,var=divmod(j,6);slotj=9*src+va[var];q=fixture(slotj)
  active_input={20*t+k:q[slot(t,va[a],va[b])] for t in range(17) for k,(a,b) in enumerate(m6) if slot(t,va[a],va[b]) in q}
  assert active_input==A[j]
  y={k:fixtures[slotj].get(9*(k//6)+va[k%6],Q(0)) for k in range(252)};image={}
  for k,c in y.items():
   if not c:continue
   for row,a in A[k].items():
    image[row]=qa(image.get(row,Q(0)),qm(a,c))
    if not image[row]:del image[row]
  assert image==A[j]
 assert actual==44 and preactual<=403
 proof={'all162_literal_singleton_units_and_forward_reverse_DAG_bindings':True,
  'fixed28_unit_matrix_selected_original_rows_columns_and_literal_det_peer_bound':True,
  'all918_times28_unit_equations_independently_checked':True,
  'all2349_scalar_nodes_linear_topological_block_not_selfdependent':True,
  'all918_DAG_outputs_equal_saved_source_witnesses':True,
  'all918_F_R_F_full1144_coordinate_identities':True,'all252_A_R_A_identities':True,
  'all918_identity_hashes_sha256':hashlib.sha256(''.join(hashes).encode()).hexdigest(),
  'ordinary_coefficient_pairs44_ell2_omission_and_flattenings_verified':True,
  'source_rank676_kernel242_accepted_from_same_strand_source_peers':True,
  'formal_covector_proof':'For source alpha annihilating kerF, v-RFv lies in kerF by FRF=F, hence alpha RF=alpha. Only true-image RHS are lifted.',
  'no_arbitrary_alpha_numeric_RHS_or_K_rho_b1_Z2_response_bound':True,
  'unit28_adjugate_recipe_remains_factored_no_dense_inverse':True,
  'positive_known_fixture_bound':maximum,'actual_known_fixture_max_bits':actual,'actual_preop_max_bits':preactual,
  'whole_own_numeric_upper':own,'hard_numeric':150000,'RAW_read':RAW,
  'hard_RAW':32<<20,'hard_seconds':60,'hard_memory_bytes':2_000_000_000}
 ep=out/'independent-factored-full918-source-DAG-all-identities-and-known-fixture-proof.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'all918_source_identities_and_factored_rightsolve_accepted':True,
  'alpha_K_rho_b1_Z2_remain_formal':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds'],'positive':maximum,'actual':actual,'owned':own}))
if __name__=='__main__':main()
