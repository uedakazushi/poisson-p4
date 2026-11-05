#!/usr/bin/env python3
"""Independent ONE549 frame/ring and empty-positive-path native proof."""
import itertools,json,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
FRAME='AMBIENT_B18_NATIVE_ONE549_FRAME_SOURCE_20261105_v63';POS='AMBIENT_B18_NATIVE_ONE549_POSITIVE_INPUT_20261105_v63';ID='AMBIENT_B18_NATIVE_ONE549_SOURCE_POSITIVE_20261105_v63_REVIEW'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
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
 def packet(st,h):
  r=pin(ROOT/'reports'/st/'report.json',h)
  for p,q in r['input_sha256'].items():
   pin(ROOT/p,q,False)
   if p.endswith('_REVIEW/report.json'):assert json.loads((ROOT/p).read_text())['all_passed']
  en=next(iter(r['evidence_sha256']));return r,pin(ROOT/'reports'/st/en,r['evidence_sha256'][en])
 fr,S=packet(FRAME,'6fe41836acc116ca03338d59ca91515351ce61f61fd9d66b27bb68538892f66b')
 pr,P=packet(POS,'62d66bc4406244a863719b5e0ec9cdb2b780c5b2d61475c39ecc6796512f1b4f')
 pin(ROOT/'code/ambient_B18_native_one549_frame_source_20261105_v63.py',fr['code_sha256'],False)
 pin(ROOT/'code/ambient_B18_native_one549_positive_input_20261105_v63.py',pr['code_sha256'],False)
 def src(s):return ROOT/next(p for p in pins if p.endswith(s))
 original=src('combined-OLD2-D-and-four-pivot-source-factor-inventory.json')
 c0=field(original,'original619_c0_source_combinations');rec=next(x for x in field(original,'kernel617_source_recipes')if x['free_case619']==549)
 a=field(src('actual-known-bottom-first-rank-and-corrected-source-kernel.json'),'all_zero_bottom_and_Y0_source_combinations')[549]
 assert c0[549]==a['original768_combination']==S['literal_full_original768_source']==[[298,'1/4'],[649,'1']]
 assert rec==S['same617_lower_basis_recipe']and rec['subtract_actual_basis563609']==['0','0']
 aa=a['actual_a1_27_times9'];assert aa==S['full_a1_27times9_32_entries_UNEVALUATED']and len(aa)==32
 assert [[i,q]for i,q in aa if i%9 in[7,8]]==[[61,'-79/24'],[62,'-11/30']]
 chart=src('pruned-known-Cfull224-input-and-positive-source-envelopes.json')
 assert field(chart,'pivot27_source_indices')[6]==648
 for k,sourcekey in [('Gamma_two_current56','Gamma_two'),('P27_columns56_current','P27_columns56'),('Pi3_columns56_current','Pi3_columns56')]:assert S[k]==field(chart,sourcekey)
 assert S['D_two_by10_by56']==field(original,'D_two_by10_by56')
 assert S['D_ordinary_pair_order10']==field(original,'pair_order10')==[[0,7],[0,8],[1,7],[1,8],[4,7],[4,8],[5,7],[5,8],[7,7],[7,8]]
 op=field(src('actual-complete-Res4-syzygy-operator.json'),'native_S12_by32')
 assert [x for x in op['koszul28_in_first8_rows']if x[1]==2]==[[0,2,3,1],[3,2,0,-1]]
 assert S['full_original_source_S12_by32']==op
 mu=field(src('actual-T2-mu-input-source-and-cache-gate.json'),'native_S0_columns')
 assert mu[9]==[]==S['full_intrinsic_Res4_S0_column9_exactempty']and mu[8]==[[15,'-1/2']]==S['comparison_pivot648_mask0_S0_column8_nonzero']
 B=list(map(tuple,S['ordinary_basis_J18']));bi={m:i for i,m in enumerate(B)};assert len(B)==len(bi)==18
 canon=src('canonical33-to-B17-C68-G8-source-input-and-response-budget.json')
 assert B==list(map(tuple,field(canon,'basis_B17')))+[(0,0,0,0,0,0,0,2,1)]
 assert S['canonical_branch']==field(canon,'canonical30_branch')and S['trace_basis_order']==field(canon,'free4_basis_mask_order')
 table=src('exact-existing-raw33-multiplication-table.json');rb=[]
 for z,r,e in field(table,'basis33'):
  q=[0]*9
  if z:q[z-1]=1
  q[7]=r;q[8]=e;rb.append(tuple(q))
 images=[bi.get(e)for e in rb]
 assert [images.index(i)for i in range(18)]==S['raw33_selector18']==[0,7,6,3,2,8,9,22,23,20,21,14,15,12,13,24,25,30]
 for i,e in enumerate(B):
  for j,q in enumerate(e):
   if q:f=list(e);f[j]-=1;assert tuple(f)in bi
 rawtable=field(table,'table_columns_for_each_left_basis')
 for i,j in itertools.product(range(33),repeat=2):
  v={images[k]:Q(q)for k,q in rawtable[i][j]if images[k]is not None and Q(q)}
  e=tuple(x+y for x,y in zip(rb[i],rb[j]));k=bi.get(e)if images[i]is not None and images[j]is not None else None
  assert v==({k:Q(1)}if k is not None else{})
 nf=src('actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json')
 mons=field(nf,'monomials660_degree_le4');NFs=field(nf,'all660_actual_raw33_normal_forms');assert len(mons)==len(NFs)==660
 for e,v in zip(mons,NFs):
  p={images[k]:Q(q)for k,q in v if images[k]is not None and Q(q)};k=bi.get(tuple(e))
  assert p==({k:Q(1)}if k is not None else{})
 g68=src('actual-G8-C68-and-monic-second-values.json');g28=src('actual32-G8-r2ell-and-two-monic-cached-coefficients.json')
 G=[v+q for v,q in zip(field(g68,'full_G8_C68'),field(g28,'actual_G8_r2ell_eight_by4'))]
 assert G==S['already_accepted_G8_J18_C72_images_COPIED_NOT_REEVALUATED']
 assert field(g68,'G8_order')==field(g28,'source_G8_order')==S['G8_order']
 nu=[[q[k]for k in S['raw33_selector18']]for q in field(canon,'actual_nu_A_B33_cached_UNEVALUATED')]
 assert nu==S['already_accepted_nu_two_J18_images']and[v[17]for v in nu]==field(g28,'actual_nu_A_B_r2ell')==['0','0']
 assert G[0][:4]==['0','0','0','1/2']and G[3][:4]==['0','0','0','0']and all(Q(v[0])==0 for v in nu)
 ns=[{i for i,q in enumerate(v)if Q(q)}for v in nu]
 assert all(tuple(x+y for x,y in zip(B[i],B[j]))not in bi for a,b in itertools.product(ns,repeat=2)for i in a for j in b)
 nr=field(src('phase2-actual-N-B11.json'),'N_C_B44');ne=field(src('phase2-actual-T-N-G8-C16.json'),'N');first=field(src('phase5-native-adapter-firsts.json'),'N')
 N=[]
 for i in range(12):
  row=[]
  for j in range(12):
   b=nr[i][j];e=ne[i][j];f=first[i][j]
   for ri,fi in[(0,0),(3,1),(4,2),(5,5),(6,6),(1,8)]:assert b[4*ri:4*ri+4]==f[4*fi:4*fi+4]
   assert b[:8]==e[:8]
   row.append([b[:4],b[12:16],b[16:20],b[20:24],b[24:28],b[4:8],e[8:12],b[8:12],e[12:16]])
  N.append(row)
 assert N==S['N_known12_by12_by9_by4_cached_operands_ONLY']and S['N_known_J18_positions']==[0,1,2,3,4,5,6,15,16]
 owned=count([fr,S,pr,P,c0,a,aa,op,mu,B,rawtable,mons,NFs,G,nu,nr,ne,first,N,pins])+30000
 assert owned<150000 and raw+2000000<32<<20
 save(out/'positive-before-independent-empty-path-check.json',{'all_source_arrays_and_metadata_charged':True,'whole_numeric_upper':owned,'no_signed_response_job':True,'no_new_N2_N3_or_physical_coefficients':True})
 # A Boolean semiring source calculation is independent of the producer's
 # positiveD/U arithmetic: a missing path certifies exact(1,0) envelope.
 SS=[set()for _ in B]
 for row,g in[(0,3),(3,0)]:
  for base,mask in itertools.product(range(18),range(4)):
   if not Q(G[g][4*base+mask]):continue
   if not mask&1:SS[base].add(4*row+(mask^1))
   else:
    for nb in ns[0]:
     e=tuple(x+y for x,y in zip(B[base],B[nb]))
     if e in bi:SS[bi[e]].add(4*row+(mask^1))
 assert not SS[0]
 for v,support in zip(P['positive_intrinsic_source9_columns18'],SS):assert {i for i,d,u in v if int(u)}>set() if False else {i for i,d,u in v if int(u)}==support
 primitive=P['native_primitive_input_from_frozen_source'];assert primitive['basis18']==S['ordinary_basis_J18']and primitive['N_known_positions']==S['N_known_J18_positions']and primitive['N_known12_by12_by9_by4']==N
 assert primitive['G8_C72']==G and primitive['nu_two_J18']==nu and primitive['Gamma_two56']==S['Gamma_two_current56']and primitive['D_two_by10_by56']==S['D_two_by10_by56']
 Gamma=[{i for i,q in row if i<48 and Q(q)}for row in S['Gamma_two_current56']]
 D=[[{i for i,q in enumerate(row)if i<48 and Q(q)}for row in ch]for ch in S['D_two_by10_by56']]
 known=S['N_known_J18_positions'];cache={}
 def support(kind,ch,pair,target):
  key=(kind,ch,pair,target)
  if key in cache:return cache[key]
  cov=Gamma[ch]if kind==0 else D[ch][pair];z=set()
  for source,im,physical,kpos,mask in itertools.product(range(12),range(4),range(12),range(9),range(4)):
   if not Q(N[source][physical][kpos][mask])or 4*physical+(mask^im)not in cov:continue
   collision=mask&im;nb=known[kpos]
   if not collision:
    if nb==target:z.add(4*source+im)
   elif collision in[1,2]:
    if any(tuple(x+y for x,y in zip(B[nb],B[vb]))==B[target]for vb in ns[collision-1]):z.add(4*source+im)
   else:assert collision==3 # every monic pair vanishes structurally above
  cache[key]=z;return z
 cubic=[(0,17),(5,16),(6,15),(15,6),(16,5)];assert P['native_Gamma_cubic5_term_order']==[list(x)for x in cubic]
 firstbase={0:1,1:2,4:3,5:4,7:5,8:6};all_empty=0
 for ch in range(2):
  for j,(l,s)in enumerate(cubic):
   assert not(support(0,ch,-1,l)&SS[s])and P['native_Gamma_cubic_two_by5_D_U'][ch][j]==['1','0'];all_empty+=1
  for p,(a,b)in enumerate(S['D_ordinary_pair_order10']):
   terms=[(0,7+p)]+[(firstbase[x],firstbase[y])for x,y in([(a,b)]if a==b else[(a,b),(b,a)])]
   assert len(terms)==len(P['native_D_second10_unique_first_decomposition_D_U'][ch][p])
   for j,(l,s)in enumerate(terms):
    assert not(support(1,ch,p,l)&SS[s])and P['native_D_second10_unique_first_decomposition_D_U'][ch][p][j]==['1','0'];all_empty+=1
 assert P['positive_native_Gamma_plus_D_two_D_U']==[['1','0'],['1','0']]
 for key,v in P['positive_covectorL_row_cache']:
  assert {i for i,d,u in v if int(u)}==cache[tuple(key)]
  assert all(int(d)>0 and int(u)>0 for i,d,u in v)
 assert len(cache)==len(P['positive_covectorL_row_cache'])==68
 peer=pin(src('AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json'),'10844d8720465054d6b9f7a7aa408dd95dd613a2e857eb52e50cd5a71a847e92');assert peer['all_passed']
 vv=pin(src('actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json'));position=vv['original768_native23_sources'].index(648)
 assert position==2 and vv['Gamma_Xp2_native23_twochannels_rr_rell'][position]==[['0','0'],['0','0']]
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'literal_ONE549_original768_combo_298quarter_plus649_and_all32_a1_checked':True,'native_source649_Koszul2_mask1_full_S0_zero_pivot648_mask0_nonzero_distinguished':True,
  'all660_nonhomogeneous_NFs_and1089_products_independently_verified_for_J18':True,'same_canonical_C68_C28_and_actual_N_first_rr_re_frame_checked':True,
  'all_monic_pair_products_zero_by_exact_support':True,'independent_Boolean_path_semiring_all68_positive_row_supports_and18_S_supports_match':True,
  'all_positive_native_piece_intersections_empty':all_empty,'two_native_total_envelopes':[['1','0'],['1','0']],
  'unknown_N_factor_proof':'For every D second coefficient unknownL2 multiplies the ENTIRE Res4 S0 column9=0, not just a weighted kernel. For Gamma pure r2ell unknownL3 likewise multiplies S0; unknown nonplaneL2 has a z factor and cannot project to pureplane in the verified monomial J18 quotient. Known Lrr/Lre contain their full varying-monic carries. Nothing sets N2/N3 themselves to0.',
  '648_outer_coupling':'Only plane a1 positions61,62 (pivot27position6 row648) enter Gamma pureplane. Both fullL2S0 target values rr/re on both channels are0 by independentROOT92cache. Its individual S0 isNONZERO and was not zero-filled.',
  'old_quarter_and_all55_FIRST_lane_SCOPE':'OLD298quarter and all other full32 a1/55 FIRST carries remain in the separate A0/V62 source peer lane. This report by itself does NOT claim wholeONE549 or all549/617 upper functional zero.',
  'no_signed_native_response_or_new_N2N3_source_physical_gauge_jet':True,'no_worker_needed_for_this_zero_native_subcomponent':True,
  'owned_upper':owned,'RAW':raw}
 ep=out/'independent-one549-frame-J18-and-empty-positive-native-path-certificate.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_POSITIVE_ONLY':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'native_subcomponent_exact_zero_proved_without_signed_job':True,'whole_ONE549_claimed':False,'seconds':time.monotonic()-t}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owned,'empty_piece_paths':all_empty,'seconds':report['seconds']}))
if __name__=='__main__':main()
