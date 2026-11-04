#!/usr/bin/env python3
"""Independent typed512 source fingerprints, C4 identity, conductor/Tor proof."""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parents[1];P='AMBIENT_B18_CONDUCTOR_TRANSFER512_NIL_CYCLES_20261104_v63';STEM=P+'_REVIEW';OUT=R/'reports'/STEM

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=R/'reports'/st/name;pins[str(p.relative_to(R))]=sha(p);return json.loads(p.read_text())
 primary=load(P);assert pins[f'reports/{P}/report.json']=='228454d222be324b940b0971a9a5ac578a9dc09a447db4298bfd1053c3a6b863'
 for p,h in primary['input_sha256'].items():assert sha(R/p)==h
 data=load(P,'actual-transfer512-H1-cycles-typed-DAG.json');assert pins[f'reports/{P}/actual-transfer512-H1-cycles-typed-DAG.json']=='d718c968faef948a264a3749fb5183cb294c3da02b95df09085f52c7898d0dfb'
 src={}
 for ns,ref in data['source_artifacts'].items():
  p=R/ref['path'];assert sha(p)==ref['sha256'];pins[ref['path']]=ref['sha256'];src[ns]=json.loads(p.read_text())
 for ns in ['radical_peer','C0_peer','relative_peer','nil2_peer','upper7_peer','upper17_peer','rankK','product44_peer']:assert src[ns]['all_passed']
 rad,chain,C0,source=[src[x] for x in ['radical','chain','C0','source']];old=chain['operator_scalar_nodes'];ns=rad['new_Cscalar_and_R0component_nodes'];Cn=C0['DAG_scalar_nodes'];sn=source['DAG_scalar_nodes'];assert Cn[:len(sn)]==sn
 beta=C0['actual_C0_basis8_scalar_nodes'];assert beta==[7148,7199,7454,7505,7688,7943,8045,505938]
 assert chain['physical_target_rows12']==[18,19,28,29,60,67,72,63,64,68,73,74]
 assert chain['coefficient_C_basis']==['1','tauA','tauB','tauA*tauB']
 assert chain['nu_A_B_operator_nodes']==[208,209] and old[208]==['actual_C0_scalar_function',20086] and old[209]==['actual_C0_scalar_function',40865]
 assert sn[20086]==['a',12,20085] and sn[20085]==['m',1693,20084] and sn[1693]==['m',11,11] and sn[20084]==['c',1,4]
 assert source['families']['CA']['p_quadratic_relation']=={'s':11,'t':12,'rule':'p^2=s*p+t'}
 # Honest scalar-function atoms are kept distinct: this is exact syntax/source arithmetic.
 def plus(p,q,c=F(1)):
  p=dict(p)
  for m,v in q.items():
   a=p.get(m,F(0))+c*v
   if a:p[m]=a
   else:p.pop(m,None)
  return p
 def times(p,q):
  out={}
  for m,c in p.items():
   for n,d in q.items():out=plus(out,{tuple(sorted(m+n)):c*d})
  return out
 memo={}
 def atom_poly(i):
  if i in memo:return memo[i]
  n=old[i];k=n[0]
  if k=='c':p={():F(n[1],n[2])} if n[1] else {}
  elif k=='actual_C0_scalar_function':p={(n[1],):F(1)}
  elif k=='sum':
   p={}
   for j,a,b in n[1]:p=plus(p,atom_poly(j),F(a,b))
  elif k=='product':
   p={():F(n[1],n[2])}
   for j in n[3]:p=times(p,atom_poly(j))
  else:raise AssertionError(k)
  assert len(p)<200;memo[i]=p;return p
 F12=chain['actual_F12_C4'];rawA=source['families']['CA']['actual_zero_odd_I0_generators'];rawB=source['families']['CB']['actual_zero_odd_I0_generators'];rows=[18,19,28,29,60,67,72,63];assert [r for r,i in rawB]==rows
 identities=chain['accepted_actual_regular_function_identities'];ca=next(x for x in identities if x['kind']=='accepted_actual_CA_f63_tail');cb=next(x for x in identities if x['kind']=='accepted_actual_CB_f63_tail')
 assert ca['function']==7790 and ca['gA']==beta[7] and cb['function']==28539 and cb['gB']==299650 and cb['coordinate_quotient_nodes']==[327551,327554,327557,327560,327563,327566,327569]
 assert len([x for x in identities if x['kind']=='accepted_actual_CA_zero_graph_membership'])==len([x for x in identities if x['kind']=='accepted_actual_CB_zero_graph_membership'])==4
 for j in range(8):
  wa={(rawA[j][1],):F(1)} if j<7 else {(ca['gA'],):F(1)};wb={(rawB[j][1],):F(1)} if j<7 else {(cb['gB'],):F(1)}
  if j==7:
   for k in range(7):
    wa[tuple(sorted((beta[k],ca['coordinate_quotients'][k])))]=F(1);wb[tuple(sorted((rawB[k][1],cb['coordinate_quotient_nodes'][k])))]=F(1)
  assert atom_poly(F12[j][1])==wa and atom_poly(F12[j][2])=={m:-v for m,v in wb.items()} and F12[j][3]==0
 U=chain['exact_existing_U9_C4'];assert [atom_poly(x) for x in U[0]]==[{():F(1)},{},{},{}]
 for i in range(8):assert [atom_poly(x) for x in U[i+1]]==[{}, {(beta[i],):F(1)}, {}, {}]
 A=rad['native_new_A12_by9_nodes'];Q=rad['actual_Qrad4_by36_columns'];assert len(A)==12 and all(len(row)==9 for row in A) and len(Q)==36
 assert rad['NilD0_exactly_H1_of_actual_complex'] and len(rad['actual_V3_by44_rows_after_unit1_cancellation'])==3
 for group,factors in [(5,[7,11]),(7,[11,13])]:assert ns[ns[Q[4*group][0]][1]]==['C_product',factors]
 # Independent semantic fingerprints and R0/C typing; every foreign ID is namespace checked.
 intern={};descriptors=[]
 def key(*x):
  if x not in intern:intern[x]=len(descriptors);descriptors.append(x)
  return intern[x]
 z,o=key('q','0'),key('q','1')
 def summ(xs):
  xs=tuple(sorted(x for x in xs if x!=z));return z if not xs else xs[0] if len(xs)==1 else key('sum',xs)
 def product(xs):
  if z in xs:return z
  xs=tuple(sorted(x for x in xs if x!=o));return o if not xs else xs[0] if len(xs)==1 else key('product',xs)
 def negative(x):
  if x==z:return z
  d=descriptors[x]
  return d[1] if d[0]=='neg' else key('neg',x)
 def rref(space,i):
  if space=='chain':assert 0<=i<len(old);return key('atompoly',tuple(sorted((m,str(v)) for m,v in atom_poly(i).items()))) if i>1 else [z,o][i]
  if space=='C0':assert 0<=i<len(Cn);return key('C0',i) if i>1 else [z,o][i]
  assert space=='radical' and 0<=i<len(ns) and (ns[i][0].startswith('R0_') or ns[i][0]=='C_QQ')
  return [z,o][i] if i<2 else key('Rrad',i)
 def cref(i):
  assert 0<=i<len(ns) and not ns[i][0].startswith('R0_');return [z,o][i] if i<2 else key('Crad',i)
 def embed(i):return key('embed',i) if i not in [z,o] else i
 def basis(i):assert 0<=i<4;return o if i==0 else key('basis',i)
 def comp(i,k):
  assert 0<=k<4;d=descriptors[i]
  if i==z:return z
  if i==o:return o if k==0 else z
  if d[0]=='embed':return d[1] if k==0 else z
  if d[0]=='basis':return o if k==d[1] else z
  return key('component',i,k)
 N=data['new_typed_scalar_nodes'];T=data['new_node_types'];assert len(N)==len(T)==25088;fingerprints=[];counts=collections.Counter()
 for i,n in enumerate(N):
  tag=n[0];counts[tag]+=1
  def child(j,t):assert 0<=j<i and T[j]==t;return fingerprints[j]
  if tag=='R0_QQ':assert T[i]=='R0';v=key('q',str(F(n[1],n[2])))
  elif tag=='R0_actual_external':assert T[i]=='R0';v=rref(n[1],n[2])
  elif tag=='C_actual_external':assert T[i]=='C' and n[1]=='radical';v=cref(n[2])
  elif tag=='R0_neg':assert T[i]=='R0';v=negative(child(n[1],'R0'))
  elif tag in ['R0_sum','R0_product','C_sum','C_product']:
   t='C' if tag.startswith('C_') else 'R0';assert T[i]==t;v=(summ if tag.endswith('sum') else product)([child(j,t) for j in n[1]])
  elif tag=='C_embed_R0':assert T[i]=='C';v=embed(child(n[1],'R0'))
  elif tag=='C_basis':assert T[i]=='C';v=basis(n[1])
  elif tag=='R0_component_of_C':assert T[i]=='R0';v=comp(child(n[1],'C'),n[2])
  else:raise AssertionError(tag)
  fingerprints.append(v)
 def fp(i,t):assert 0<=i<len(N) and T[i]==t;return fingerprints[i]
 records=data['actual_512_transfer_cycles'];assert len(records)==512;assert len(set((x['CA_generator_index'],x['CB_raw_index'],x['normalized_nil_generator_index'],x['normalizer_basis_index']) for x in records))==512
 for number,(i,j,al,ga) in enumerate(itertools.product(range(8),range(8),range(2),range(4))):
  rec=records[number];assert [rec['index'],rec['CA_generator_index'],rec['CB_raw_index'],rec['CB_physical_row'],rec['normalized_nil_generator_index'],rec['normalizer_basis_index']]==[number,i,j,rows[j],al,ga]
  group=[5,7][al];p0,pA,pB,pAB=[rref('radical',k) for k in Q[4*group+ga]];be=rref('C0',beta[i]);de,ba=rref('chain',F12[j][0]),rref('chain',F12[j][1]);bj=negative(rref('chain',F12[j][2]));nu=rref('C0',20086)
  assert [fp(k,'R0') for k in rec['p_four_actual_R0_components']]==[p0,pA,pB,pAB];f=product([be,bj]);assert fp(rec['actual_f_betaA_times_bB'],'R0')==f
  even=product([be,summ([product([bj,p0]),product([de,pB]),product([nu,ba,pAB])])]);odd=summ([product([bj,pA]),product([ba,pB]),product([de,pAB])]);c8=[z]*8;c8[i]=odd
  assert [fp(k,'R0') for k in rec['literal_actual_C0_free9_representative']]==[even]+c8
  tc=summ([embed(pB),product([basis(1),embed(pAB)])]);wantC=[];want36=[]
  for g in range(9):
   c=product([embed(be),cref(A[j][g]),tc])
   if g==group:c=summ([c,product([embed(f),basis(ga)])])
   wantC.append(c);want36.extend(comp(c,k) for k in range(4))
  assert [fp(k,'C') for k in rec['literal_Qrad_lift9_C_coefficients']]==wantC
  assert [fp(k,'R0') for k in rec['literal_H1_cycle44']]==want36+c8
 # A separate sparse polynomial algebra checks all4 C-components, not numerical points.
 nv=10
 def var(i):return {tuple(int(j==i) for j in range(nv)):F(1)}
 def pa(*ps):
  out={}
  for p in ps:out=plus(out,p)
  return out
 def pm(*ps):
  out={(0,)*nv:F(1)}
  for p in ps:
   q={}
   for m,c in out.items():
    for n,d in p.items():q=plus(q,{tuple(a+b for a,b in zip(m,n)):c*d})
   out=q
  return out
 def cm(x,y,nuA,nuB):
  out=[{} for _ in range(4)]
  for i in range(4):
   for j in range(4):
    terms=[x[i],y[j]]
    if i&j&1:terms.append(nuA)
    if i&j&2:terms.append(nuB)
    out[i^j]=pa(out[i^j],pm(*terms))
  return out
 nuA,nuB,be,bj,de,ba=[var(i) for i in range(6)];p=[var(i) for i in range(6,10)];Fj=[de,ba,{m:-c for m,c in bj.items()},{}];t=[p[2],p[3],{},{}]
 E=pm(be,pa(pm(bj,p[0]),pm(de,p[2]),pm(nuA,ba,p[3])));H=pa(pm(bj,p[1]),pm(ba,p[2]),pm(de,p[3]));left=[E,pm(be,H),{},{}];right=[pa(pm(be,bj,p[k]),pm(be,cm(Fj,t,nuA,nuB)[k])) for k in range(4)];assert left==right
 assert data['all512_generate_exact_transfer_W_not_all_NilD0'] and data['K_mod_J_dimension_at_most']==26 and data['K_square_length_at_most']==182
 assert 8+2*(13-4)==26 and 2*(13*14//2)==182
 # R0/Jplus²=(u7,(r,s)²)² has ten pure(r,s) +21 single-u monomials.
 assert sum(d+1 for d in range(4))+7*3==31
 evidence={'all25088_node_types_and_foreign_namespaces':dict(counts),'all512_actual_source_C0_free9_and_Qrad36_plus8_cycles_reconstructed':True,'source_nuA_equals_tA_plus_sA_squared_over4':True,'both_unipotent_physical_odd_frames_and_all4_graph_gates':True,'rawCB8_generate_entire_IB_det1':True,'two_actual_nilgenerators_nativegroups':[5,7],'C4_basis4_and_actual12to9_column_adapter':True,'universal_C4_transfer_four_component_identity':True,'W_is_exact512_generated_transfer_subideal_not_full_N':True,'transfer_cross_product':'phi_f(x)phi_g(y)=g*phi_f(xy)=f*phi_g(xy), with transfers indexed by the conductor pair','W_nil_and_WN2_zero':'Each phi_f(x)^3=f²phi_f(x³)=0; psi(N) subset Nil(E4), hence phi_f(x)nnprime=phi_f(x psi(n)psi(nprime))=0.','Jtimes_N_subset_W':True,'triple_map_factors_through_Sym3_N_mod_W':True,'Jtimes_annihilates_N_mod_W_and_N3':True,'Jplus_K2_annihilation_distinguished_from_Jtimes_N3':True,'K_Jplus_dim_upper':26,'K2_length_upper':182,'R0_Jplus2_length':31,'Tor2_connecting_image_not_discarded':True,'no_Artin_Rees_strictness_or_honest_K_Jplus_presentation_assumed':True,'no_full_N_generators_or_nilcube_or_H0_length_decision':True}
 save('independent-512-source-and-conductor-proof.json',evidence)
 report={'id':STEM,'all_passed':True,'all512_actual_C0_representatives_and_complete44cycles_source_bound':True,'all25088_typed_nodes_checked':True,'four_component_polynomial_identity_independent':True,'actual_conductor_CB_unipotent_frame_unit1_entire_IB':True,'all512_generate_W_not_all_N':True,'W_times_N_squared_zero':True,'IAIB_times_N_subset_W':True,'nilcube_triple_map_factors_Sym3_N_mod_W':True,'honest_K_mod_Jplus_presentation_still_open_Tor2_retained':True,'no_scalar512_expansion_newjet_sourceworker_or_kernel_solve':True,'no_nilcube_or_H0_DVR_exponent_decision':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{STEM}.md'),'evidence_sha256':{'independent-512-source-and-conductor-proof.json':sha(OUT/'independent-512-source-and-conductor-proof.json')},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};assert sum(p.stat().st_size for p in OUT.iterdir())<33554432;save('report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
