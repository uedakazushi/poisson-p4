#!/usr/bin/env python3
"""Factored source rightsolve and918 fixed generator identities, alpha formal."""
from fractions import Fraction as Q
from pathlib import Path
import itertools,json,math,resource,signal,time
from ambient_B18_ell2_T2_G2_strand_rank_input_review_20261105_v63 import sha,count,fields,sp
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_FULL918_RIGHTSOLVE_SOURCE_20261105_v63'
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);pins={}
 def rd(p):pins[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def fld(p,w):pins[p]=sha(ROOT/p);return fields(ROOT/p,w)
 def put(n,d):
  (o/n).write_text(json.dumps(d,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in o.iterdir())<33554432
 accepted=[('reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63_REVIEW/report.json','ed6954d48503e2220e1d7ab29526a64ead78387ab00e1012f70b454172ff20e3'),('reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63_REVIEW/report.json','fbc71ad533749bf42d4636c627300d5bd67b9265e6aec2bb9f98b07104c0b4a9')]
 accepted_data=[]
 for p,h in accepted:
  accepted_data.append((p,rd(p)));assert accepted_data[-1][1]['all_passed'] and pins[p]==h
 pb='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_SOURCE_PROFILE_20261105_v63/';D=rd(pb+'literal-active-G2-linear-strand-matrix-singleton-units-and-residual-source.json');A=list(map(sp,D['literal_matrix_columns']));pivs=D['fixed_singleton_unit_pivots'];assert len(A)==252 and len(pivs)==162
 rb='reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_INPUT_GATE_20261105_v63/';I=rd(rb+'pruned46-by41-integer-rank-source-and-live-bound.json')
 cert=rd('reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63/fixed46-by41-integer-Schur-rank-and-full-strand-kernel-recipes.json');factor=rd('reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63_REVIEW/independent-fixed-minor-all41-factorizations-and-full-strand-proof.json')
 assert cert['active_G2_linear_strand_rank']==190 and factor['active_rank']==190
 det=int(factor['selected_integer_determinant']);fac=[[Q(int(x),det) for x in c] for c in factor['all41_original_column_factor_numerators_common_det']]
 sr=cert['selected_residual_rows46'];sc=cert['selected_residual_columns41'];R=I['original340_row_indices46'];C=I['original252_column_indices41'];selected=[C[j] for j in sc]
 F=list(map(sp,fld('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['full_G2_linear234_by102_columns'])['full_G2_linear234_by102_columns']));assert len(F)==102
 G1=fld('reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_INPUT_GATE_20261105_v63/actual-full-G1-and186-mixed-adjoint-source-input.json',['exact_full_G1_polynomial_columns'])['exact_full_G1_polynomial_columns'];L1=[]
 for col in G1[:17]:
  v={}
  for head,ex,x in col:
   if sum(ex)==1:v[9*head+ex.index(1)]=Q(x)
  L1.append(v)
 inherited=set(p for p,h in accepted)
 for path,rep in accepted_data:
  for k,h in rep['input_sha256'].items():
   if k in pins:assert pins[k]==h;inherited.add(k)
  for k,h in rep['evidence_sha256'].items():
   p=str(Path(path).parent/k)
   if p in pins:assert pins[p]==h;inherited.add(p)
 assert set(pins)<=inherited,('unbound_fixed_source',set(pins)-inherited)
 va=[0,1,4,5,7,8];vz=[2,3,6];pairs=list(itertools.combinations(range(3),2));m9=[p for p in itertools.combinations_with_replacement(range(9),2) if p!=(8,8)];ix9={m:i for i,m in enumerate(m9)};m6=[p for p in itertools.combinations_with_replacement(range(6),2) if p!=(5,5)]
 def outslot(target,a,b):return 44*target+ix9[tuple(sorted((a,b)))]
 # Literal scalar linear DAG, plus ONE fixed28 scalar-unit solve BLOCK.
 nodes=[['constant','0']];index={('constant','0'):0}
 def node(v):
  k=json.dumps(v,separators=(',',':'))
  if k not in index:index[k]=len(nodes);nodes.append(v)
  return index[k]
 inp=[node(['input',j]) for j in range(1144)]
 def scale(q,x):
  if not q or not x:return 0
  if q==1:return x
  return node(['scale',str(q),x])
 def total(v):
  v=[x for x in v if x]
  if not v:return 0
  if len(v)==1:return v[0]
  return node(['sum',v])
 active_inputs=[inp[outslot(t,va[a],va[b])] for t in range(17) for a,b in m6];assert len(active_inputs)==340
 rhs=active_inputs.copy();xnodes={};liveR=set(range(340));liveC=set(range(252));steps=[];source_rows=[{j:A[j][i] for j in range(252) if i in A[j]} for i in range(340)]
 for p in pivs:
  i,j=p['row340'],p['column252'];q=Q(p['literal_nonzero_unit']);kind=p['singleton_kind'];after=sorted(liveC-{j});row=[(k,source_rows[i][k]) for k in after if k in source_rows[i]]
  if kind=='row_singleton':
   assert not row;xnodes[j]=scale(1/q,rhs[i])
   for k in sorted(liveR-{i}):
    if k in A[j]:rhs[k]=total([rhs[k],scale(-A[j][k],xnodes[j])])
  else:assert kind=='column_singleton' and not (set(A[j])&(liveR-{i}))
  steps.append({'row340':i,'column252':j,'kind':kind,'pivot':str(q),'saved_rhs_node':rhs[i],'remaining_row_coefficients':[[k,str(x)] for k,x in row]});liveR.remove(i);liveC.remove(j)
 U=[[I['residual_integer_matrix46_by41'][i][j] for j in sc] for i in sr]
 input28=[scale(Q(I['row_integer_scaling_D46'][i]),rhs[R[i]]) for i in sr]
 block={'kind':'fixed_QQ_unit28_solve','matrix_integer28_by28':U,'determinant_literal':str(det),'input_nodes28':input28,'output_component_count':28,'inverse_recipe':'adj(U)/det(U),neverexpanded','source_nested_minor_certificate_path':'reports/AMBIENT_B18_ELL2_T2_G2_STRAND_RANK_20261105_v63/fixed46-by41-integer-Schur-rank-and-full-strand-kernel-recipes.json','source_committed_minor_bits':306,'source_product_difference_bits':613,'arbitrary_numeric_RHS_application_NOT_authorized':True}
 for k,j in enumerate(selected):xnodes[j]=node(['unit28_component',0,k])
 for j in liveC-set(selected):xnodes[j]=0
 for s in reversed(steps):
  if s['kind']=='column_singleton':xnodes[s['column252']]=scale(1/Q(s['pivot']),total([s['saved_rhs_node']]+[scale(-Q(q),xnodes[k]) for k,q in s['remaining_row_coefficients']]))
 active_out=[xnodes[j] for j in range(252)];outputs=[0]*918
 for j,n in enumerate(active_out):i,a=divmod(j,6);outputs[9*i+va[a]]=n
 def ps(j,a,b):return inp[outslot(j,a,b)]
 def qs(i,head,a,b):return inp[outslot(17+3*i+head,a,b)]
 for i,z in enumerate(vz):
  for j in range(17):
   for a in va:outputs[9*(42+17*i+j)+a]=scale(-1,ps(j,z,a))
   outputs[9*(42+17*i+j)+z]=scale(-1,ps(j,z,z))
 for i,j in pairs:
  for a in range(17):outputs[9*(42+17*i+a)+vz[j]]=scale(-1,ps(a,vz[i],vz[j]))
 for ip,(i,j) in enumerate(pairs):
  for head in range(3):
   for a in va:
    v=[scale(L1[k].get(9*head+a,Q(0)),outputs[9*(42+17*i+k)+vz[j]]) for k in range(17)];v.append(scale(-1,qs(i,head,vz[j],a)));outputs[9*(93+3*ip+head)+a]=total(v)
 zero_recipes=[(0,0,0,0,1,-1),(0,1,0,1,1,-1),(0,2,1,0,2,1),(1,0,0,0,2,-1),(1,1,2,0,1,1),(1,2,0,2,2,-1),(2,1,1,1,2,-1),(2,2,1,2,2,-1)]
 for head in range(3):
  for pair,coefvar,qi,a,b,sgn in zero_recipes:outputs[9*(93+3*pair+head)+vz[coefvar]]=scale(sgn,qs(qi,head,vz[a],vz[b]))
 assert len(outputs)==918 and sum(x!=0 for x in outputs)==676
 # Sparse projection R_A A: use the accepted ALL41 factors, not a newinverse.
 # Forward row-singletons force originalsource coordinates. Column-singletons
 # are recovered in reverse by the literal unchanged row equation.
 def active_projection(basis,mode):
  exact=mode=='exact';one=Q(1) if exact else (1,1);zero=Q(0) if exact else (1,0);v={}
  for s in steps:
   if s['kind']=='row_singleton' and s['column252']==basis:v[basis]=one
  if basis in C:
   for j,q in zip(selected,fac[C.index(basis)]):
    if q:v[j]=q if exact else ex(q)
  for s in reversed(steps):
   if s['kind']!='column_singleton':continue
   j=s['column252'];p=Q(s['pivot']);q=one if j==basis else zero
   for k,a in s['remaining_row_coefficients']:
    a=Q(a)/p
    if exact:q=ck(q+ck(a*((Q(1) if k==basis else Q(0))-v.get(k,Q(0)))))
    else:q=ba(q,bm(ex(a),ba(one if k==basis else zero,v.get(k,zero))))
   if exact:
    if q:v[j]=q
    else:v.pop(j,None)
   elif q[1]:v[j]=q
  return v
 maximum=0
 def mark(a):
  nonlocal maximum
  maximum=max(maximum,a[0].bit_length(),a[1].bit_length());assert maximum<=384;return a
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else (1,0)
 def ba(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def bm(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 posA=[active_projection(j,'positive') for j in range(252)]
 # A source-basis fixture supplies an exact column F_i times one variable.
 def fixture(slot):
  i,v=divmod(slot,9);q={}
  for k,x in F[i].items():
   row,u=divmod(k,9)
   if (u,v)==(8,8):continue
   z=outslot(row,u,v);assert z not in q;q[z]=x
  return q
 def circuit_fixture(slot,mode,pa):
  exact=mode=='exact';q=fixture(slot);q=q if exact else {i:ex(x) for i,x in q.items()};zero=Q(0) if exact else (1,0)
  i,v=divmod(slot,9);active=6*i+va.index(v) if i<42 and v in va else None
  aout=pa[active] if active is not None else {};assert (active is not None) or not any(outslot(j,va[a],va[b]) in q for j in range(17) for a,b in m6)
  vals=[]
  for n in nodes:
   op=n[0]
   if op=='constant':val=zero
   elif op=='input':val=q.get(n[1],zero)
   elif op=='unit28_component':
    # The accepted activeprojection evaluates the same source unit solve
    # on all252originalAcolumns, including source eliminations/restorations.
    # Its selected28 value is exactly the residual selected coordinate.
    val=aout.get(selected[n[2]],zero)
   elif op=='scale':
    if exact:val=ck(Q(n[1])*vals[n[2]])
    else:val=bm(ex(n[1]),vals[n[2]])
   else:
    assert op=='sum';val=zero
    for k in n[1]:val=ck(val+vals[k]) if exact else ba(val,vals[k])
   vals.append(val)
  for i in range(28):
   check=zero
   for j,c in enumerate(U[i]):
    if exact:check=ck(check+ck(Q(c)*aout.get(selected[j],zero)))
    else:check=ba(check,bm(ex(c),aout.get(selected[j],zero)))
   if exact:assert check==vals[input28[i]]
   else:ba(check,vals[input28[i]])
  y={j:vals[k] for j,k in enumerate(outputs) if vals[k] if exact} if exact else {j:vals[k] for j,k in enumerate(outputs) if vals[k][1]}
  return y
 # The selected coordinate ofactiveprojection afterreverse neverchanges the
 # residualselected entries, so the unitblock oracle is source-faithful.
 pos=[];support=0
 for j in range(918):
  y=circuit_fixture(j,'positive',posA);pos.append(y);support+=len(y);q={}
  for k,c in y.items():
   src,var=divmod(k,9)
   for row,x in F[src].items():
    target,u=divmod(row,9)
    if (u,var)==(8,8):continue
    z=outslot(target,u,var);q[z]=ba(q.get(z,(1,0)),bm(ex(x),c))
  for z,x in fixture(j).items():ba(q.get(z,(1,0)),ex(x))
 constant=count([D,A,pivs,I,cert,factor,fac,R,C,F,G1,L1,steps,U,block,nodes,outputs,pins]);owned_prior=constant+count([posA,pos])+8*support+4*len(nodes)+4*1144+2*918+12000;assert owned_prior<150000
 raw_inputs=sum((ROOT/p).stat().st_size for p in pins);output_numeric_upper=count([steps,U,block,nodes,outputs])+4*support+12000;raw_prior=raw_inputs+256*output_numeric_upper+2000000;assert raw_prior<33554432
 put('phase-positive-before-signed-source-fixtures.json',{'phase':'complete_source_and_positive_bounds_before_knownQQ_fixture_arithmetic','nodes':len(nodes),'outputs':918,'nonzero_gauge_outputs676':True,'positive_fixture_support':support,'positive_fixture_committed_bits':maximum,'temporary_prior_bits':3*maximum+1,'owned_numeric_prior':owned_prior,'literal_inputs_RAW':raw_inputs,'whole_RAW_prior':raw_prior,'active_fixedunit_minor_guard306_613':True,'arbitrary_expanded_R_or_numeric_RHSheight_NOT_claimed':True,'no_alpha_K_rho_b1_Z2_or_response_values':True})
 actual=0
 def ck(q):
  nonlocal actual
  actual=max(actual,abs(q.numerator).bit_length(),q.denominator.bit_length());assert actual<=maximum;return q
 exactA=[active_projection(j,'exact') for j in range(252)]
 for j,y in enumerate(exactA):
  assert not (set(y)&(set(C)-set(selected)))
  q={}
  for k,c in y.items():
   for i,x in A[k].items():
    q[i]=ck(q.get(i,Q(0))+ck(c*x))
    if not q[i]:del q[i]
  assert q==A[j]
 fixtures=[];fingerprint=[]
 for j in range(918):
  y=circuit_fixture(j,'exact',exactA);q={}
  for k,c in y.items():
   src,var=divmod(k,9)
   for row,x in F[src].items():
    target,u=divmod(row,9)
    if (u,var)==(8,8):continue
    z=outslot(target,u,var);q[z]=ck(q.get(z,Q(0))+ck(x*c))
    if not q[z]:del q[z]
  assert q==fixture(j)
  fixtures.append([[k,str(x)] for k,x in sorted(y.items())]);fingerprint.append(sha_bytes(json.dumps(sorted(q.items()),default=str,separators=(',',':')).encode()))
 assert sum(x!=0 for x in outputs)==190+306+102+54+24==676
 data={'input_shape1144':[26,44],'output_shape918':[102,9],'variable_order9':['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell'],'quadratic_pair_order44':m9,'actual_full_G2linear_columns':[[[i,str(x)] for i,x in sorted(v.items())] for v in F],'scalar_linear_nodes':nodes,'active_unit_block28':block,'active_forward_and_reverse162_source_steps':steps,'active_source_outputs252':active_out,'full_source_outputs918':outputs,'explicit_zero_gauge_source_slots242':[j for j,n in enumerate(outputs) if not n],'all918_source_generator_R_F_columns':fixtures,'all918_F_R_F_identity_fingerprints':fingerprint,'all252_active_A_R_A_identity_checked':True,'all918_times28_actual_unit_block_equations_checked':True,'all918_F_R_F_identity_checked':True,'rank_decomposition':[190,306,102,54,24],'kernel_decomposition':[62,126,51,3],'all1144_output_compatibilities_on_image':'F*R*F=F onALL918originalsourcebasisgenerators impliesF*R=P_image andall1144coordinatesonANYtrueimageinput. No arbitraryoutput is asserted liftable.','alpha_Kfull_FORMAL':'Foranysourcecovectoralphaannihilatingtheknown242G3kernel,alpha*R*F=alpha. Thusalpha*R is anoutputcovector;actualalphavaluesandresponse/support/cancellation/heightremainseparategates. UseKfullwithPsi,notKU.','no_expanded676_by1144_inverse':'onefixed28unitadjugateblockplusliteral162singletonand306/102/54/24tensorrecipes','fixed_block_numeric_bound':'unitmatrixentriesknown;itsorderedminors<=306bits/productdifference<=613. NoexpandedinverseorapplicationtoarbitrarynumericRHSisevaluated.','no_K_alpha_rho_b1_Z2_upper_or_higher_source_responses':True}
 ep='actual-factored918-rightsolve-source-circuit-and-all-generator-identities.json';put(ep,data)
 report={'id':STEM,'all_passed':True,'SOURCE_factored_linearoperator_only':True,'all918_generator_identities_checked':True,'all1144_image_coordinate_compatibilities_implied':True,'source_rank676_kernel242':True,'alpha_remains_FORMAL':True,'positive_known_fixture_prior_bits':maximum,'known_fixture_actual_bits':actual,'fixed_unit_minor_bound306_temp613':True,'owned_numeric_prior':owned_prior,'whole_RAW_prior':raw_prior,'no_actual_K_alpha_rho_b1_Z2_upper_or_higher_response':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep:sha(o/ep)},'seconds':time.monotonic()-start};put('report.json',report);print(json.dumps({'all_passed':True,'sha256':sha(o/'report.json'),'nodes':len(nodes),'source_owned':owned_prior,'positive':maximum,'actual':actual}))
def sha_bytes(b):
 import hashlib
 return hashlib.sha256(b).hexdigest()
if __name__=='__main__':main()
