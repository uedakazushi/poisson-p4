#!/usr/bin/env python3
"""Independent native617 kernel restriction and shared72 positive SOURCE audit."""
import hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_NATIVE_SHARED72_SOURCE_FINE_20261105_v63_REVIEW'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pins={};raw=0;bits=0;pre=0;cursor={}
 def read(p,h=None,parse=True):
  nonlocal raw
  p=ROOT/p;k=str(p.relative_to(ROOT));z=sha(p)
  if h:assert h==z,(k,h,z)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=z;return json.loads(p.read_text())if parse else None
 def packet(st,h=None):
  r=read(f'reports/{st}/report.json',h);assert r['all_passed']
  for p,z in r.get('input_sha256',{}).items():read(p,z,False)
  name=next(iter(r['evidence_sha256']));e=read(f'reports/{st}/{name}',r['evidence_sha256'][name]);return r,e
 sr,S=packet('AMBIENT_B18_NATIVE_SHARED72_INPUT_SOURCE_20261105_v63','4f01895013eb59d8d1d5721eff04d8131e9a73f52f2449f9613731e311c704e6')
 fr,F=packet('AMBIENT_B18_NATIVE_SHARED72_FINE_POSITIVE_20261105_v63','072c20b07c37929a97f49da4b0f8235f846d30a4d05f9ebada0d5c858edffa76')
 kr,K=packet('AMBIENT_B18_NATIVE_POINT_KERNEL105_SOURCE_20261105_v63','15b6ff370efd02939a1fac290ca6c24c1237ab3046f6a82984d22b1071e8d534')
 mr,M=packet('AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63','faec4ff020e7f1af70cbf58d825ee09c85aeb31000d540b0f261c48ade511754')
 br,B=packet('AMBIENT_B18_NATIVE_KERNEL105_FUNCTIONAL_SOURCE_20261105_v63','be68ef59b96e555a8f9e8b0f6a5592648072f2904cb9b1fc943e8dd630ba5d95')
 nr=read('reports/AMBIENT_B18_NATIVE_ONE554_NUMERIC_SOURCE_20261105_v63_REVIEW/report.json','f9581e3b3c18c9f98856317d9cf72da89f7330d3149b45ddffb82d4cd9bcb8a9');assert nr['all_passed']
 pp='reports/AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63/ONE554-positive-four-native-and-nonzero644-outer-input.json';P=read(pp,nr['input_sha256'][pp])
 shared=read('reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63_REVIEW/report.json','36a453a118ec96340d5fd72c57e65ddebbaeeea108c53040f34bc3b5471cee8a');assert shared['all_passed']
 sd=read('reports/AMBIENT_B18_ELL2_T2_UPPER605_SHARED_SOURCE_20261105_v63/complete605-shared-source-linearity-and-global24-opaque-frontier.json')
 meta=read('reports/AMBIENT_B18_NATIVE_KERNEL617_TO105_SOURCE_20261105_v63/metadata-distinct-basis-support-and-source-uniqueness.json','c07cdf2f87eb389532ee1dc3156f443cb86177d889c00234d705d024943d0faa')
 for r,fn in [(sr,'ambient_B18_native_shared72_input_source_20261105_v63.py'),(fr,'ambient_B18_native_shared72_fine_positive_20261105_v63.py'),(kr,'ambient_B18_native_point_kernel105_source_20261105_v63.py'),(mr,'ambient_B18_native_kernel617_to105_source_20261105_v63.py'),(br,'ambient_B18_native_kernel105_functional_source_20261105_v63.py')]:read('code/'+fn,r['code_sha256'],False)
 p=P['full_primitive_input_J18_N_G_nu_GammaD'];assert B['same_J18_with_r2ell_primitive_input']==p
 pruned=S['native_primitive_known_cache'];assert pruned['N_known_first8_by12_by9_by4']==p['N_known12_by12_by9_by4'][:8]
 for k in ['basis18','N_known_positions','G8_C72','nu_two_J18','Gamma_two56','D_two_by10_by56']:assert pruned[k]==p[k]
 assert pruned['N_known_positions']==[0,1,2,3,4,5,6,15,16]
 def src(suffix):return next(k for k in pins if k.endswith(suffix))
 op_path=src('actual-complete-Res4-syzygy-operator.json');op=field(ROOT/op_path,'native_S12_by32')
 mu_path=src('actual-T2-mu-input-source-and-cache-gate.json');MU=field(ROOT/mu_path,'native_S0_columns')
 oldpath=src('combined-OLD2-D-and-four-pivot-source-factor-inventory.json');c619=field(ROOT/oldpath,'original619_c0_source_combinations');recipes=field(ROOT/oldpath,'kernel617_source_recipes')
 assert len(c619)==619 and all(len(v)==len({i for i,x in v})for v in c619)
 expectedmaps={r['free_case619']:r for r in M['all617_source_maps']};assert set(expectedmaps)=={r['free_case619']for r in recipes}==set(range(619))-{563,609}
 pairs=list(itertools.combinations(range(8),2));want=[[i,k,j,1]for k,(i,j)in enumerate(pairs)]+[[j,k,i,-1]for k,(i,j)in enumerate(pairs)]
 assert sorted(op['koszul28_in_first8_rows'])==sorted(want)
 assert op['free4_zero_source_columns']==[[8,28,1],[9,29,1],[10,30,1],[11,31,1]]
 base=list(map(tuple,p['basis18']));bi={m:i for i,m in enumerate(base)};table=[[bi.get(tuple(a+b for a,b in zip(x,y)))for y in base]for x in base]
 parsed=count([sr,S,fr,F,kr,K,mr,M,br,B,nr,P,shared,sd,meta,p,pruned,op,MU,c619,recipes,expectedmaps,pairs,want,base,bi,table,pins])
 prior=parsed+50000;assert prior<150000 and raw+4_000_000<32*1024*1024
 save(out/'positive-before-independent-kernel-and-fine-paths.json',{'all_input_numeric_labels_and_duplicate_maps':parsed,'whole_owned_upper':prior,'RAW_input_plus_output_upper':raw+4000000,'integer_commit_hard':384,'integer_temporary_hard':1153,'fine_SOURCE_only_no_signed_responses':True,'point_override_permitted_only_after_full48_kernel_proof':True})
 def ck(x):
  nonlocal bits
  bits=max(bits,abs(x).bit_length());assert bits<=384,(cursor,bits);return x
 def prod(x,y):
  nonlocal pre
  b=abs(x).bit_length()+abs(y).bit_length();pre=max(pre,b);assert b<=1153;return ck(x*y)
 def plus(x,y):
  nonlocal pre
  b=max(abs(x).bit_length(),abs(y).bit_length())+1;pre=max(pre,b);assert b<=1153;return ck(x+y)
 def lcm(x,y):return prod(x//math.gcd(x,y),y)
 def q(x):v=Q(x);ck(v.numerator);ck(v.denominator);return v
 def qm(x,y):prod(x.numerator,y.numerator);prod(x.denominator,y.denominator);v=x*y;ck(v.numerator);ck(v.denominator);return v
 def qa(x,y):
  a=prod(x.numerator,y.denominator);b=prod(y.numerator,x.denominator);plus(a,b);prod(x.denominator,y.denominator);v=x+y;ck(v.numerator);ck(v.denominator);return v
 def madd(v,k,x):
  z=qa(v.get(k,Q(0)),x)
  if z:v[k]=z
  else:v.pop(k,None)
 def axpy(v,w,a):
  for k,x in w.items():madd(v,k,qm(a,x))
 G=[[q(x)for x in row]for row in p['G8_C72']];nu=[[q(x)for x in row]for row in p['nu_two_J18']]
 assert all(not row[i]for row in nu for i in [0,5,6])
 assert all(table[i][j]is None for x in nu for y in nu for i,a in enumerate(x)if a for j,b in enumerate(y)if b)
 c=[row[3]for row in G];assert all(row[:3]==[0,0,0]for row in G)
 assert c==list(map(Q,['1/2','25/8','0','0','0','5','-125/2','0']))
 # Full48x128 point source matrix, from original signs and actual point trace.
 S0=[]
 for native in range(128):
  col,mask=divmod(native,4);v={}
  if col<28:
   i,j=pairs[col]
   for row,g,sgn in [(j,i,-1),(i,j,1)]:
    if not mask and c[g]:madd(v,4*row+3,sgn*c[g])
  else:v={4*(col-20)+mask:Q(1)}
  assert v=={i:q(x)for i,x in MU[native]}
  S0.append(v)
 complement=K['fixed23_source_complement'];outrows=K['fixed23_output_rows'];assert len(complement)==len(outrows)==23
 diag=[]
 for j,src in enumerate(complement):
  col=[S0[src].get(i,Q(0))for i in outrows];assert col==[(-c[0]if j<7 else Q(1))if k==j else Q(0)for k in range(23)];diag.append(col[j])
 det=Q(1)
 for x in diag:det=qm(det,x)
 assert det==Q(-1,128)==Q(K['fixed23_minor_determinant'])
 ownK=[]
 for rec in K['explicit105_kernel_basis_columns']:
  src=rec['free_native_source'];assert src not in complement;v={src:Q(1)}
  for j,row in enumerate(outrows):
   if S0[src].get(row):madd(v,complement[j],-S0[src][row]/diag[j])
  image={}
  for col,x in v.items():axpy(image,S0[col],x)
  assert not image and v=={i:q(x)for i,x in rec['column']};ownK.append(v)
 free=[r['free_native_source']for r in K['explicit105_kernel_basis_columns']]
 assert len(free)==len(set(free))==105 and set(free)|set(complement)==set(range(128))and not set(free)&set(complement)
 assert M['kernel105_free_native_indices']==free
 # Source correction and reconstruction precede restriction to605.
 rawc=[{i-640:q(x)for i,x in row if i>=640 and Q(x)}for row in c619];maps={};allcoords=0
 for rec in reversed(recipes):
  case=rec['free_case619'];v=dict(rawc[case])
  for b,x in zip([563,609],rec['subtract_actual_basis563609']):axpy(v,rawc[b],-q(x))
  image={}
  for src,x in v.items():axpy(image,S0[src],x)
  assert not image and all(not v.get(i)for i in range(112,128))
  coords={i:v[src]for i,src in enumerate(free)if v.get(src)};re={}
  for i,x in reversed(list(coords.items())):axpy(re,ownK[i],x)
  assert re==v
  e=expectedmaps[case];assert e['subtract563609']==rec['subtract_actual_basis563609']
  assert v=={i:q(x)for i,x in e['literal_native128_constant_combination']}and coords=={i:q(x)for i,x in e['kernel105_coordinates']}
  maps[case]=(v,coords);allcoords+=len(coords)
 good=sd['cases605'];assert len(good)==605 and good==sorted(set(maps)-set(meta['opaque_second_case_indices12']))
 ids=sorted({i for case in good for i in maps[case][1]});assert len(ids)==72 and sum(len(maps[case][1])for case in good)==120
 assert ids==S['kernel72_original105_indices']==meta['opaque_free605_basis_indices']
 restrictions={r['kernel105']:r for r in S['kernel72_restrictions']};assert set(restrictions)==set(ids)
 for i,rec in restrictions.items():assert ownK[i]=={src:q(x)for src,x in rec['original_native_source_combination']}
 raw75=sorted({src for i in ids for src in ownK[i]});actualraw={src for case in good for src in maps[case][0]}
 assert len(raw75)==75 and len(actualraw)==72 and set(raw75)-actualraw=={4,8,12}and raw75==S['raw75_source_indices']
 assert all(S0[i]for i in[4,8,12])
 programs={r['native128']:r for r in S['raw75_original_Koszul_program']};assert set(programs)==set(raw75)
 for src,rec in programs.items():
  col,mask=divmod(src,4);i,j=pairs[col]
  assert rec=={'native128':src,'inputmask':mask,'Koszul_entries':[[i,j,1],[j,i,-1]]}
 mapping={r['case619']:{i:q(x)for i,x in r['kernel105_coordinate_pairs']}for r in S['all605_constant_native_coordinate_recipes']}
 assert set(mapping)==set(good)and all(mapping[case]==maps[case][1]for case in good)
 # Signed SOURCE arithmetic stops here. Remaining arithmetic is positive D/U.
 R={tuple(k):{i:(int(d),int(u))for i,d,u in row}for k,row in P['reused68_positive_covectorL_rows']};assert len(R)==68
 def ex(x):v=q(x);return(v.denominator,abs(v.numerator))if v else(1,0)
 def path(xs):
  n=d=1
  for den,u in xs:
   if not u:return None
   n=prod(n,u);d=prod(d,den)
  return n,d
 def accumulate(factory,size):
  ds=[1]*size;us=[0]*size
  for i,n,d in factory():ds[i]=lcm(ds[i],d)
  for i,n,d in factory():assert ds[i]%d==0;us[i]=plus(us[i],prod(n,ds[i]//d))
  return[(d,u)if u else(1,0)for d,u in zip(ds,us)]
 def sumdu(vs):return accumulate(lambda:((0,u,d)for d,u in reversed(vs)if u),1)[0]
 EG=[[ex(x)for x in row]for row in G];EN=[[ex(x)for x in row]for row in nu]
 # Independent common-class priors before enumerating any72 fine envelopes.
 def book(vs):
  d=1
  for den,u in vs:d=lcm(d,den)
  return d,max([d]+[prod(u,d//den)for den,u in vs])
 classes={'R68':[z for row in R.values()for z in row.values()],'G8':[z for row in EG for z in row],
  'nu':[(1,0)]+[z for row in EN for z in row],'kernelweights':[ex(x)for i in ids for x in ownK[i].values()]}
 fineD=1;fineU=1<<20;finebook={}
 for name,vs in classes.items():d,u=book(vs);finebook[name]=list(map(str,(d,u)));fineD=prod(fineD,d);fineU=prod(fineU,u)
 fineprior=max(fineD.bit_length(),fineU.bit_length());assert fineprior==242 and 34*48*3*2*4*19<1<<20
 mappingrows={}
 for case,coords in mapping.items():
  d=1
  for x in coords.values():d=lcm(d,x.denominator)
  u=0
  for x in coords.values():u=plus(u,prod(abs(x.numerator),d//x.denominator))
  mappingrows[case]=(d,u)
 composed=max([fineprior]+[max(prod(fineD,d).bit_length(),prod(fineU,u).bit_length())for d,u in mappingrows.values()if u]);assert composed==254
 phasepath='reports/AMBIENT_B18_NATIVE_SHARED72_FINE_POSITIVE_20261105_v63/phase-before-positive-S-and-fine-contractions.json'
 phase=read(phasepath,fr['evidence_sha256']['phase-before-positive-S-and-fine-contractions.json']);assert phase['positive_classes']==finebook and phase['positive_common_D']==str(fineD)and phase['positive_common_U_paths2pow20']==str(fineU)
 cub=[(0,17),(5,16),(6,15),(15,6),(16,5)];dpairs=[(0,7),(0,8),(1,7),(1,8),(4,7),(4,8),(5,7),(5,8),(7,7),(7,8)];vf={0:1,1:2,4:3,5:4,7:5,8:6}
 expected={r['kernel105']:r for r in F['native_kernel72_Gamma_D_sum_fine']};fp=dict(F['all72_positive_S18by48_fingerprints_AFTER_point_override']);point={i:v for i,v in F['raw_positive_point_bounds_BEFORE_exact_source_override']};values={};nonzeropoint=0;allpieces=0
 enc=lambda z:list(map(str,z))
 for i in reversed(ids):
  cursor={'phase':'positive_shared72','kernel105':i}
  def spaths():
   for src,w in reversed(list(ownK[i].items())):
    op=programs[src];im=op['inputmask']
    for row,g,sgn in reversed(op['Koszul_entries']):
     for gb in reversed(range(18)):
      for mask in reversed(range(4)):
       x=EG[g][4*gb+mask]
       if not x[1]:continue
       collision=mask&im;om=mask^im
       if not collision:
        n,d=path((ex(w),x));yield 48*gb+4*row+om,n,d
       elif collision in[1,2]:
        for vb,y in reversed(list(enumerate(EN[collision-1]))):
         dest=table[gb][vb]
         if y[1]and dest is not None:
          n,d=path((ex(w),x,y));yield 48*dest+4*row+om,n,d
  SS=accumulate(spaths,18*48)
  pospoint=[[k,*enc(z)]for k,z in enumerate(SS[:48])if z[1]];assert pospoint==point[i];nonzeropoint+=bool(pospoint)
  # The independent FULL48 kernel identity above, and ONLY that identity,
  # permits the exact point coefficient0. Keep every positive higher cell.
  SS[:48]=[(1,0)]*48
  encoded=[[[col,*enc(SS[48*b+col])]for col in range(48)if SS[48*b+col][1]]for b in range(18)]
  assert hashlib.sha256(json.dumps(encoded,separators=(',',':')).encode()).hexdigest()==fp[i]
  def dot(row,b):return sumdu([(prod(d,SS[48*b+k][0]),prod(u,SS[48*b+k][1]))for k,(d,u)in row.items()if u and SS[48*b+k][1]])
  piecebook=[];channels=[]
  for ch in range(2):
   cg=[dot(R[0,ch,-1,rb],sb)for rb,sb in cub];ds=[]
   for pi,(a,b)in enumerate(dpairs):
    terms=[(0,7+pi)]+[(vf[x],vf[y])for x,y in([(a,b)]if a==b else[(a,b),(b,a)])]
    ds.append([dot(R[1,ch,pi,rb],sb)for rb,sb in terms])
   A=sumdu(cg);D=sumdu([z for row in ds for z in row]);total=sumdu([A,D]);values[i,ch]=total
   channels.append([enc(A),enc(D),enc(total)]);piecebook.append([list(map(enc,cg)),[[enc(z)for z in row]for row in ds]]);allpieces+=len(cg)+sum(map(len,ds))
  assert channels==expected[i]['two_channels_Gamma_D_sum_DU']
  assert hashlib.sha256(json.dumps(piecebook,separators=(',',':')).encode()).hexdigest()==expected[i]['fine68_piece_DU_sha256']
 mapped={case:v for case,v in F['all605_native_sum_fine_twochannels']}
 for case in reversed(good):
  for ch in range(2):
   # Keep original denominator multiplication, including weights, even when
   # the final positive rational quotient could be reduced further.
   def paths():
    for i,w in reversed(list(mapping[case].items())):
     z=path((ex(w),values[i,ch]))
     if z:n,d=z;yield 0,n,d
   v=accumulate(paths,1)[0];assert enc(v)==mapped[case][ch]
 # Independently reconstruct coarse basis341 then605354 source priors.
 def flat(x):
  if isinstance(x,list):
   for y in x:yield from flat(y)
  else:yield ex(x)
 broad={'known_partial_N9':list(flat(p['N_known12_by12_by9_by4'])),'canonical_G8':list(flat(p['G8_C72'])),'monics':list(flat(p['nu_two_J18'])),
  'Gamma_and_D_scalar_coefficients':[ex(z)for row in p['Gamma_two56']for k,z in row]+[ex(z)for ch in p['D_two_by10_by56']for row in ch for z in row],
  'native_kernel105_weights':[ex(x)for col in ownK for x in col.values()]}
 bd=1;bu=1<<28
 for name,vs in broad.items():
  d,u=book(vs);h=B['positive_operand_classes'][name];assert(str(d),str(u),len(vs))==(h['D'],h['U_with_one'],h['operand_count']);bd=prod(bd,d);bu=prod(bu,u)
 assert max(bd.bit_length(),bu.bit_length())==341 and(str(bd),str(bu))==(S['native_basis_common_D'],S['native_basis_common_U'])
 assert 3*2*12*16*9*18*19*11==B['direct_paths_per_target_upper']<1<<28
 coarse={case:v for case,v in S['all605_positive_native_sum_envelopes_same_two_channels']};coarsebits=341
 for case,(d,u)in mappingrows.items():
  z=(prod(bd,d),prod(bu,u))if u else(1,0);assert[enc(z),enc(z)]==coarse[case];coarsebits=max(coarsebits,max(x.bit_length()for x in z))
 assert coarsebits==354 and 3*354+1==1063
 assert F['new_fine_positive_prior']==fineprior and F['new_coordinate_composition_positive_prior']==composed and F['source_committed_union_prior']==354 and F['source_temporary_union_prior']==1063
 assert nonzeropoint>0 and allpieces==72*68
 observed=parsed+count([S0,diag,ownK,free,rawc,maps,ids,restrictions,raw75,actualraw,programs,mapping,R,G,nu,EG,EN,classes,finebook,mappingrows,phase,values,SS,encoded,expected,fp,point,mapped,broad])+10000
 assert observed<prior and all(sha(ROOT/k)==h for k,h in pins.items())
 e={'all_passed':True,'full48x128_point_matrix_rebuilt_from_original_Koszul_and_actual_trace_G8':True,'fixed23_unit_block_determinant':'-1/128',
  'all105_kernel_columns_derived_by_fixed_diagonal_complement_no_rank_search':True,'all105_full48_kernel_checks_and_free_I105':True,
  'all617_literal_full48S0cN0_and128_reconstructions_verified':True,'all619_source_index_uniqueness_verified':True,'raw563609_corrections_preserved':True,
  'actual72_kernel_indices':ids,'actual605_coordinate_nonzeros120':True,'formal_raw75_indices':raw75,'formal_extras4_8_12_point_NONZERO':True,
  'no_individual_raw_auxiliary_declared_actual':'Only each constant kernel restriction is an actual projected functional. Missing N2mixed/N3 remain formal; kernel S0 removes whole highest terms and nonplane unknownL2 cannot enter pureplane in verified J18 quotient.',
  'point_positive_S0_nonzero_columns':nonzeropoint,'point_override_ONLY_after_exact_full48_kernel_proof':True,'all_higher_S_cells_and_single_nu_FIRSTcross_retained':True,
  'all72_complete_S_fingerprints_and4896_fine_DU_pieces_checked':allpieces,'all72_Gamma_D_total_DU_and605_bothchannel_DU_match':True,'denominator_certificates_match_exactly_not_just_ratios':True,
  'same_literal_frame_for_68_R_rows_as_independent_f958_peer':True,'fine_prior242_mapping254_and_inherited341_354_union1063_reproduced':True,
  'fine_method':'Reverse source/mask/base loops; two-pass raw-denominator LCM and scaled positive integer numerator sums. No producer imports or signed native/605 response values.',
  'scope':'SOURCE/fine positive only. OLD605 proof is separate36a453; FIRST/Psi/fullcomparison/native92 outer remain. Fine outputs must precede Q; whole immutable code/custody/newROOTGO required.',
  'whole_owned_upper':prior,'observed_held_upper':observed,'RAW_input_plus_output_upper':raw+4000000,'committed_bits_observed':bits,'temporary_precheck_bits_observed':pre}
 ep=out/'independent-full617-kernel72-and-all-fine605-DU-source-proof.json';save(ep,e)
 rep={'id':ID,'all_passed':True,'frozen':True,'SOURCE_POSITIVE_ONLY_no_native_or605_responses':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start};save(out/'report.json',rep);signal.alarm(0)
 print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'owned':prior,'observed':observed,'bits':bits,'pre':pre,'positive_nonzero_points':nonzeropoint,'seconds':rep['seconds']}))
if __name__=='__main__':main()
