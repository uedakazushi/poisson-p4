#!/usr/bin/env python3
"""Independent cache-frame and monomial-keyed positive whole605 SOURCE audit."""
import hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_ELL2_T2_UPPER605_POSITIVE_WHOLE_20261105_v63';ID=BASE+'_REVIEW'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
 t0=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pins={};raw=0;bits=0;pre=0;cursor={};parsepeak=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=ROOT/p;k=str(p.relative_to(ROOT));z=sha(p)
  if h:assert z==h,(k,h,z)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=z;return json.loads(p.read_text())if parse else None
 def packet(st,h=None):
  r=pin(f'reports/{st}/report.json',h);assert r['all_passed']
  for p,z in r.get('input_sha256',{}).items():pin(p,z,False)
  name=next(iter(r['evidence_sha256']));return r,pin(f'reports/{st}/{name}',r['evidence_sha256'][name])
 wr,W=packet(BASE,'9431e686c987488602817964dcaa3bee47177a2684c776314478a16dfea45a8c')
 cr,S=packet('AMBIENT_B18_ELL2_T2_UPPER605_OLD_FIRST_CACHE_SOURCE_20261105_v63','1edbeec47f225e128a1a55f50847a8a12de7e64b039b63030cfc0d77fd190548')
 for r,fn in[(wr,'ambient_B18_ell2_T2_upper605_positive_whole_20261105_v63.py'),(cr,'ambient_B18_ell2_T2_upper605_oldfirst_cache_source_20261105_v63.py')]:pin('code/'+fn,r['code_sha256'],False)
 def path(s):return next(p for p in pins if p.endswith(s))
 # Stage the original FIRST sources before loading the positive/native books.
 cf=pin(path('pruned-known-Cfull224-input-and-positive-source-envelopes.json'))
 old=field(ROOT/path('actual-first640-X-mod-m2.json'),'source_columns640_rows56_sparse_dual9')
 nat=field(ROOT/path('actual-native128-X-mod-m2.json'),'source_columns128_rows48_sparse_dual9')
 parsepeak=count([wr,W,cr,S,cf,old,nat,pins])+18000;assert parsepeak<150000
 piv=S['point27_source_indices'];assert piv==cf['pivot27_source_indices']and len(piv)==27 and S['Gamma_two']==cf['Gamma_two']
 L0=dict(S['native_L0_selected_columns']);Le=dict(S['native_Le_selected_columns']);need0=set();neede=set()
 recs={(r['pivot27'],r['variable9']):r for r in S['individual_FIRST_column_recipes162']}
 assert set(recs)==set(itertools.product(range(27),[0,1,4,5,7,8]))
 for(p,var),rec in sorted(recs.items(),reverse=True):
  orig=piv[p];d=rec['first_column_recipe'];assert rec['original768']==orig
  if var==8:
   kind,arg=cf['pivot27_ell_recipes'][p]
   if kind=='literal':assert d=={'kind':'literal56','column':arg}
   else:
    assert kind=='native'and arg==orig-640 and d=={'kind':'native_ell_SUM','native128_source':arg,'Se':cf['native_Se'][arg],'S0':cf['native_S0'][arg]}
    need0.update(k for k,x in d['Se']);neede.update(k for k,x in d['S0'])
  else:
   dual=8 if var==7 else var+1;src=old[orig]if orig<640 else nat[orig-640]
   want=[[row,dict(v)[dual]]for row,v in src if Q(dict(v).get(dual,'0'))]
   assert d['kind']=='literal56'and d['column']==want and d['original_old8_dual_index']==dual
   if orig>=640:assert d['native_bottom8_exact0']and all(row<48 for row,x in want)
   if var==7:assert want==cf['pivot27_r_columns'][p]
 assert set(L0)==need0 and set(Le)==neede
 assert all(v==cf['native_L0'][k]for k,v in L0.items())and all(v==cf['native_Le'][k]for k,v in Le.items())
 del old,nat,cf,L0,Le
 j=pin(path('617-direct-old-known-first-formal-Jacobian-positive-input.json'))
 oldsource=pin(path('combined-OLD2-D-and-four-pivot-source-factor-inventory.json'))
 assert S['shared_literal_C0_program532']==j['literal_C0_program']and S['shared_full_typed_program757']==oldsource['typed_program']
 assert S['Gamma_two']==j['Gamma_two']and S['D_two_by10_by56']==oldsource['D_two_by10_by56']
 assert S['OLD_Gamma_source_columns618']==j['pruned_original_old_columns_roots']and S['OLD_D_source_columns618']==oldsource['D_old_root_columns']
 assert S['raw619_original768_constant_sources']==oldsource['original619_c0_source_combinations']==j['original619_source768_combinations']
 rec617={r['free_case619']:r['subtract_actual_basis563609']for r in oldsource['kernel617_source_recipes']}
 assert all(rec617[i]==ab for i,ab in S['corrected605_recipes_basis563609'])
 known=field(ROOT/path('actual-known-bottom-first-rank-and-corrected-source-kernel.json'),'all_zero_bottom_and_Y0_source_combinations')
 parsepeak=max(parsepeak,count([wr,W,cr,S,j,oldsource,known,pins])+18000);assert parsepeak<150000
 for i,r in enumerate(known):assert r['original768_combination']==S['raw619_original768_constant_sources'][i]and r['full56_Xgamma1_equals_Xp0a1']and r['actual_a1_27_times9']==S['raw619_actual_a1_27times9'][i]
 assert len(known)==619;del known
 boundary={i:(kind,v)for i,kind,v in j['literal_boundary_firsts_point_r_ell_and_membership']}
 f79raw=field(ROOT/path('phase1-independent-full9-first640-source-proof.json'),'first640');f79={i:dict(v)for i,v in f79raw['all79_atom_values_full10']};del f79raw
 physicalraw=field(ROOT/path('actual-source105-ell2-raw33-values.json'),'all_requested164_root_values');phys={i:dict(v)for i,v in physicalraw};del physicalraw
 sel=S['raw33_B21_selector'];assert sel==[0,7,6,5,4,3,2,1,8,9,22,23,20,21,14,15,12,13,24,25,30]
 actual=extension=0
 for fn,mode,v in S['actual_or_formal_B21_boundary96']:
  kind,first3=boundary[fn];assert len(v)==21 and not Q(v[0])and list(map(Q,[v[0],v[8],v[9]]))==list(map(Q,first3))
  if fn in phys:
   assert mode=='actual_physical164_B21'and list(map(Q,v))==[Q(phys[fn].get(k,'0'))for k in sel];actual+=1
   if fn in f79:assert list(map(Q,v[:10]))==[Q(f79[fn].get(k,'0'))for k in range(10)]
  else:
   assert mode=='FULL9firsts_with_arbitrary_higher_extension_ONLY_after_whole605_certificate'and fn in f79
   assert list(map(Q,v))==[Q(f79[fn].get(k,'0'))for k in range(10)]+[Q(0)]*11;extension+=1
 assert actual==extension==48 and len(boundary)==96
 ring=pin(path('ONE508-source-closure-B21-and-exact-carry-zeros.json'));assert S['ordinary_B21_basis']==ring['ordinary_basis_B21']and ring['raw33_to_B21_selector21']==sel
 rbasis=pin('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json')
 nfp='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json';pin(nfp,'4ba0113ce973c01e156188d22c3d593e83c53ddd301ad731ee13feb52aceb20f',False)
 mons=field(ROOT/nfp,'monomials660_degree_le4');nfs=field(ROOT/nfp,'all660_actual_raw33_normal_forms')
 b=list(map(tuple,S['ordinary_B21_basis']));bi={m:i for i,m in enumerate(b)};assert len(b)==len(bi)==21
 rawmon=[]
 for z,r,e in rbasis['basis33']:
  m=[0]*9
  if z:m[z-1]=1
  m[7]=r;m[8]=e;rawmon.append(tuple(m))
 images=[bi.get(m)for m in rawmon];assert [images.index(i)for i in range(21)]==sel
 for i,k in itertools.product(range(33),repeat=2):
  got={images[j]:Q(x)for j,x in rbasis['table_columns_for_each_left_basis'][i][k]if images[j]is not None and Q(x)}
  target=bi.get(tuple(a+c for a,c in zip(rawmon[i],rawmon[k])))if images[i]is not None and images[k]is not None else None
  assert got==({target:Q(1)}if target is not None else{})
 for m,v in zip(mons,nfs):
  got={images[j]:Q(x)for j,x in v if images[j]is not None and Q(x)};target=bi.get(tuple(m));assert got==({target:Q(1)}if target is not None else{})
 assert len(mons)==len(nfs)==660
 parsepeak=max(parsepeak,count([wr,W,cr,S,j,oldsource,boundary,f79,phys,ring,rbasis,mons,nfs,b,bi,rawmon,images,pins])+18000);assert parsepeak<150000
 del j,oldsource,boundary,f79,phys,ring,rbasis,mons,nfs,rawmon,images
 n92=pin(path('actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json'))
 direct=pin(path('AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json'),'10844d8720465054d6b9f7a7aa408dd95dd613a2e857eb52e50cd5a71a847e92');assert direct['all_passed']
 assert S['accepted_native92_original_sources23']==n92['original768_native23_sources']==piv[4:]
 assert S['accepted_native92_Gamma_plane_targets23_twochannels_rr_rell']==n92['Gamma_Xp2_native23_twochannels_rr_rell'];del n92,direct
 br,BO=packet('AMBIENT_B18_ELL2_T2_UPPER605_BOOLEAN_PROFILE_20261105_v63','02ed9a5d789ae48a3d28c68076437de6250979f0ed8f671a186f3c3731b7b9da')
 nr,NF=packet('AMBIENT_B18_NATIVE_SHARED72_FINE_POSITIVE_20261105_v63','072c20b07c37929a97f49da4b0f8235f846d30a4d05f9ebada0d5c858edffa76')
 peer=pin('reports/AMBIENT_B18_NATIVE_SHARED72_SOURCE_FINE_20261105_v63_REVIEW/report.json','0384d2623fab486eb0dbd1a4df239bd04e6a7ff778f7464c454ec324c541683e');assert peer['all_passed']
 oldpeer=pin('reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63_REVIEW/report.json','36a453a118ec96340d5fd72c57e65ddebbaeeea108c53040f34bc3b5471cee8a');assert oldpeer['all_passed']
 inv=pin('reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63/complete605-shared-source-linearity-and-global24-opaque-frontier.json')
 assert S['cases605']==inv['cases605']==W['all605_case_coverage_checked']and S['opaque_cases12_kept_outside_value_scope']==inv['opaque_cases12']
 qpath=path('known-Psi-point-correction-and-full-b1-source-annihilation.json');Qrows=field(ROOT/qpath,'Q_two6')
 qreport=pin('reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/report.json');assert qreport['all_passed']and qreport['evidence_sha256'][pathlib.Path(qpath).name]==pins[qpath]
 # Everything retained now is an immutable source or positive envelope. The
 # signed original source/frame comparisons above computed no new response.
 parsed=count([wr,W,cr,S,br,BO,nr,NF,peer,oldpeer,inv,Qrows,qreport,b,bi,piv,recs,pins])
 owner=max(parsepeak,2*parsed+50000);assert owner<150000 and raw+4_000_000<32*1024*1024
 save(out/'positive-before-independent-whole605-semiring.json',{'staged_literal_frame_peak':parsepeak,'whole_numeric_upper':owner,'RAW_input_plus_output_upper':raw+4000000,'physical48_actual_and48_FULL9first_formal_higher_only_after_corrected605':True,'source_committed_hard384_temporary1153':True,'native_inherited354_1063_retained':True,'no_signed_response_or_new_sourcejet':True})
 def ck(x):
  nonlocal bits
  bits=max(bits,abs(x).bit_length());assert bits<=384,(cursor,bits);return x
 def mulint(x,y):
  nonlocal pre
  pre=max(pre,abs(x).bit_length()+abs(y).bit_length());assert pre<=1153;return ck(x*y)
 def addint(x,y):
  nonlocal pre
  pre=max(pre,max(abs(x).bit_length(),abs(y).bit_length())+1);assert pre<=1153;return ck(x+y)
 def lcm(x,y):return mulint(x//math.gcd(x,y),y)
 def ex(x):v=Q(x);return(ck(v.denominator),ck(abs(v.numerator)))if v else(1,0)
 def prod(a,c):return(mulint(a[0],c[0]),mulint(a[1],c[1]))if a[1]and c[1]else(1,0)
 def sumdu(vs,maximum=False):
  vs=[v for v in vs if v[1]]
  if not vs:return(1,0)
  d=1
  for a,u in reversed(vs):d=lcm(d,a)
  us=[mulint(u,d//a)for a,u in reversed(vs)]
  if maximum:return d,max(us)
  u=0
  for x in us:u=addint(u,x)
  return d,u
 zero=(0,)*9
 def combine(terms):
  groups={}
  for m,v in terms:
   if v[1]:groups.setdefault(m,[]).append(v)
  return {m:sumdu(vs)for m,vs in groups.items()}
 def psum(terms):return combine((m,prod(a,v))for a,p in reversed(terms)for m,v in reversed(list(p.items())))
 def pmul(a,c):
  return combine((tuple(x+y for x,y in zip(m,n)),prod(v,w))for m,v in reversed(list(a.items()))for n,w in reversed(list(c.items()))if tuple(x+y for x,y in zip(m,n))in bi)
 leaves={fn:{b[i]:ex(x)for i,x in enumerate(v)if Q(x)}for fn,mode,v in S['actual_or_formal_B21_boundary96']}
 assert sum(map(len,leaves.values()))==152
 C=dict(leaves)
 for i,d in S['shared_literal_C0_program532']:
  cursor={'phase':'positive_literal532','node':i}
  if d[0]=='c':v={zero:ex(Q(d[1],d[2]))}if d[1]else{}
  elif d[0]=='n':v=dict(C[d[1]])
  elif d[0]=='a':v=psum([((1,1),C[d[1]]),((1,1),C[d[2]])])
  else:assert d[0]=='m';v=pmul(C[d[1]],C[d[2]])
  C[i]=v
 ccells=sum(len(C[i])for i,d in S['shared_literal_C0_program532']);assert ccells==119
 boolC=[[i,sorted(bi[m]for m in C[i])]for i,d in S['shared_literal_C0_program532']]
 assert hashlib.sha256(json.dumps(sorted(boolC)).encode()).hexdigest()==BO['Boolean_program_fingerprints']['C0']
 G=[{i:ex(x)for i,x in row if Q(x)}for row in S['Gamma_two']];DD=[[{i:ex(x)for i,x in enumerate(row)if Q(x)}for row in ch]for ch in S['D_two_by10_by56']]
 gc=dict(S['OLD_Gamma_source_columns618']);dc=dict(S['OLD_D_source_columns618']);old3={}
 for old,roots in gc.items():old3[old]=[sumdu([prod(G[ch][row],C[node].get(b[20],(1,0)))for row,(kind,node)in roots if kind=='C0'and row in G[ch]])for ch in range(2)]
 del C,boolC
 T={}
 for i,d in S['shared_full_typed_program757']:
  cursor={'phase':'positive_typed757','node':i}
  if d[0]=='actual_C0_scalar_function':v=dict(leaves[d[1]])
  elif d[0]=='sum':v=psum([(ex(Q(a,den)),T[k])for k,a,den in d[1]])
  else:
   assert d[0]=='product';v={zero:ex(Q(d[1],d[2]))}if d[1]else{}
   for k in d[3]:v=pmul(v,T[k])
  T[i]=v
 tcells=sum(map(len,T.values()));assert tcells==324
 boolT=[[i,sorted(bi[m]for m in v)]for i,v in sorted(T.items())];assert hashlib.sha256(json.dumps(boolT).encode()).hexdigest()==BO['Boolean_program_fingerprints']['typed']
 old2={}
 for old,roots in gc.items():
  for ch in range(2):old3[old][ch]=sumdu([old3[old][ch]]+[prod(G[ch][row],T[node].get(b[20],(1,0)))for row,(kind,node)in roots if kind=='typed'and row in G[ch]])
 for old,roots in dc.items():old2[old]=[sumdu([prod(DD[ch][pair][row],T[node].get(b[10+pair],(1,0)))for row,node in roots for pair in reversed(range(10))if row in DD[ch][pair]])for ch in range(2)]
 del T,boolT
 L0={k:{i:ex(x)for i,x in row if Q(x)}for k,row in S['native_L0_selected_columns']};Le={k:{i:ex(x)for i,x in row if Q(x)}for k,row in S['native_Le_selected_columns']}
 first={}
 for key,rec in sorted(recs.items(),reverse=True):
  d=rec['first_column_recipe']
  if d['kind']=='literal56':v={i:ex(x)for i,x in d['column']if Q(x)}
  else:
   groups={}
   for name,mat in [('S0',Le),('Se',L0)]:
    for k,x in reversed(d[name]):
     for i,y in mat[k].items():groups.setdefault(i,[]).append(prod(ex(x),y))
   v={i:sumdu(vs)for i,vs in groups.items()}
  first[key]=v
 assert sum(map(len,first.values()))==1295
 dpairs=S['fixed_pair_order10'];assert dpairs==[[0,7],[0,8],[1,7],[1,8],[4,7],[4,8],[5,7],[5,8],[7,7],[7,8]]
 ftable={};bterms=[]
 for ch,p,var in itertools.product(range(2),range(27),[0,1,4,5,7,8]):
  terms=[]
  for pair,(a,c)in reversed(list(enumerate(dpairs))):
   targets=([c]if var==a else[])+([a]if var==c and a!=c else[])
   for other in targets:
    terms.extend(prod(x,first[p,other][i])for i,x in DD[ch][pair].items()if i in first[p,other])
  ftable[ch,p,var]=sumdu(terms);bterms.append([ch,p,var,len(terms)])
 assert bterms==BO['all324_FIRST_table_operational_term_counts']
 assert BO['actual_support_counts']=={'boundary96_scalar_cells':152,'C0_532_scalar_cells':119,'typed757_scalar_cells':324,'FIRST162_scalar_cells':1295,'FIRST324_scalar_dot_terms':sum(r[3]for r in bterms)}
 # Each case is a FORMAL raw extension until the certified563/609 correction.
 xp=S['accepted_native92_Gamma_plane_targets23_twochannels_rr_rell'];i92={orig:i for i,orig in enumerate(S['accepted_native92_original_sources23'])}
 raw619={};unions={}
 def uni(name,ch,z):unions[name,ch]=sumdu([unions.get((name,ch),(1,0)),z],True)
 for case in reversed(range(619)):
  c0=S['raw619_original768_constant_sources'][case];a1=S['raw619_actual_a1_27times9'][case];channels=[]
  for ch in range(2):
   cursor={'phase':'positive_raw619_formal_extension','case':case,'channel':ch}
   v=[sumdu([prod(ex(x),old3[col][ch])for col,x in c0 if col<640]),sumdu([prod(ex(x),old2[col][ch])for col,x in c0 if col<640])]
   v.append(sumdu([prod(ex(x),ftable[ch,p,var])for slot,x in a1 for p,var in[divmod(slot,9)]if var in[0,1,4,5,7,8]]))
   v.append(sumdu([prod(ex(x),ex(xp[i92[piv[p]]][ch][0 if var==8 else 1]))for slot,x in a1 for p,var in[divmod(slot,9)]if p>=4 and var in[7,8]]))
   channels.append(v)
   for name,z in zip(['raw_OLD3','raw_OLD2','raw_FIRST','raw_nativeOuter'],v):uni(name,ch,z)
  raw619[case]=channels
 ncase={i:[tuple(map(int,z))for z in channels]for i,channels in NF['all605_native_sum_fine_twochannels']};assert set(ncase)==set(S['cases605'])
 complete={};nonzero=[]
 for case,(a,beta)in reversed(S['corrected605_recipes_basis563609']):
  channel=[]
  for ch in range(2):
   v=[sumdu([raw619[case][ch][k],prod(ex(a),raw619[563][ch][k]),prod(ex(beta),raw619[609][ch][k])])for k in range(4)]
   for name,z in zip(['OLD3','OLD2','FIRST','nativeOuter'],v):uni(name,ch,z)
   uni('native',ch,ncase[case][ch]);z=sumdu(v+[ncase[case][ch]]);uni('upper',ch,z);channel.append(z)
  complete[case]=channel
  if any(u for d,u in channel):nonzero.append(case)
 Qdense=[dict(row)for row in Qrows]
 for target in range(6):uni('six',target,sumdu([prod(ex(Qdense[ch].get(target,'0')),unions['upper',ch])for ch in range(2)]))
 expected={(name,ch):tuple(map(int,z))for name,ch,z in W['positive_common_envelopes']};assert unions==expected
 nonzero=sorted(nonzero);wanted=[525,529,530,550,552,554,558,560,562,566,567,569,573,576,577,589,595,599,603,612,617]
 assert nonzero==W['structural_positive_nonzero_case_labels']==wanted
 zero584=sorted(set(S['cases605'])-set(nonzero));assert len(zero584)==584
 # Positive Q for each605 separately, solely to verify all six classes are
 # controlled by the complete common envelope. No actual class is evaluated.
 for case,zs in complete.items():
  for target in range(6):
   z=sumdu([prod(ex(Qdense[ch].get(target,'0')),zs[ch])for ch in range(2)]);d,u=unions['six',target]
   assert d%z[0]==0 and mulint(z[1],d//z[0])<=u
   if case in zero584:assert z==(1,0)
 assert bits==W['positive_whole_operation_bits']==wr['positive_whole_operation_bits']==257
 assert W['future_value_committed_prior']==wr['future_value_committed_prior']==354 and W['future_value_temporary_prior']==1063
 prior554=pin('reports/AMBIENT_B18_ELL2_T2_UPPER554_VALUES_V2_20261105_v63_REVIEW/report.json');assert prior554['all_passed']and prior554['actual_upper_two']==['0','0']and prior554['actual_full_six']==['0']*6
 active20=sorted(set(nonzero)-{554});assert len(active20)==20
 observed=parsed+count([leaves,G,DD,gc,dc,old3,old2,L0,Le,first,ftable,bterms,raw619,unions,xp,i92,ncase,complete,Qdense,expected,nonzero,zero584,prior554,pins])+18000
 assert observed<owner and all(sha(ROOT/p)==h for p,h in pins.items())
 e={'all_passed':True,'literal_cache_physical48_B21_and48_only79_FULL9firsts_checked':True,'arbitrary_higher_extensions_safe_ONLY_after_complete_corrected605_not_raw563609':True,
  'all162_FIRST_recipes_rebuilt_from_original640_native128_and_actual_L0Se_LeS0':True,'all619_full_a1_and_original_source_corrections_checked':True,
  'B21_all1089_products_and660_NFs_independently_checked':True,'full532_and757_SOURCE_graphs_keep_all_FIRSTcross':True,'positive_cells_boundary152_C0119_typed324_FIRST1295':True,
  'all324_FIRST_table_terms_verified_ordinary_rr_once':True,'accepted_native92_all23_source_targets_bound_to10844':True,'shared72_fine_actual_restriction0384_pinned':True,
  'method':'Monomial9-tuples (rather than coefficient indices), reversed term collection and two-pass integer LCM/positive scaled sums; percase corrections before common max union. No producer imports or signed response.',
  'all_raw619_FORMAL_component_envelopes_and_corrected605_COMPLETE_bounds_recomputed':True,'all_common_D_U_envelopes_match_exactly':True,'all605_six_Q_positive_bounds_divide_and_dominate':True,
  'source_zero_cases584':zero584,'source_zero_two_and_six_classes_for_these584':True,'positive_nonzero_labels21_NOT_actual_nonzero':nonzero,
  'case554_separate_accepted_actual_zero_peer':pins['reports/AMBIENT_B18_ELL2_T2_UPPER554_VALUES_V2_20261105_v63_REVIEW/report.json'],'remaining20_need_separate_value_GO':active20,
  'OLD4_outer_scope':'Complete exact12-function linearity/formal cancellation36a453 used onlyafter605 source correction. No raw609carry absent claim; its formal response and basis563 remain in every bound.',
  'scope':'SOURCE positive and structural zero proof only. No actual605/20 responses, rank/HF/length/H0 conclusion. FullPsi/comparison correction Q kept; native intermediate354 goes through fine outputs before upper/Q. WholeVALUE code/custody/ROOTGO separate.',
  'staged_source_parse_peak':parsepeak,'whole_owned_upper':owner,'observed_held_upper':observed,'RAW_input_plus_output_upper':raw+4000000,'positive_committed_bits':bits,'temporary_precheck_bits':pre,'inherited_native_union354_1063':True}
 ep=out/'independent-cache-frame-full605-positive-bound-and584-sourcezero-proof.json';save(ep,e)
 rep={'id':ID,'all_passed':True,'frozen':True,'SOURCE_POSITIVE_ONLY_no605_or20_signed_responses':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t0};save(out/'report.json',rep);signal.alarm(0)
 print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':owner,'observed':observed,'bits':bits,'pre':pre,'source_zeros':584,'positive_nonzeros':21,'seconds':rep['seconds']}))
if __name__=='__main__':main()
