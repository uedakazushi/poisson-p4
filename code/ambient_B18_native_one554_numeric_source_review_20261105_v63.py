#!/usr/bin/env python3
"""Independent absolute-path/LCM audit of the native554 SOURCE envelopes."""
import hashlib,itertools,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import field
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_NATIVE_ONE554_NUMERIC_SOURCE_20261105_v63_REVIEW'
FST='AMBIENT_B18_NATIVE_ONE554_FRAME_SOURCE_20261105_v63'
PST='AMBIENT_B18_NATIVE_ONE554_POSITIVE_INPUT_20261105_v63'
BST='AMBIENT_B18_NATIVE_ONE549_POSITIVE_INPUT_20261105_v63'
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(55);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};raw=0;cursor={};maxcommit=0;maxpre=0
 def pin(p,h=None,parse=True):
  nonlocal raw
  p=pathlib.Path(p);k=str(p.relative_to(ROOT));z=sha(p)
  if h:assert z==h,(k,h,z)
  if k not in pins:raw+=p.stat().st_size
  pins[k]=z
  return json.loads(p.read_text())if parse else None
 def packet(st,h=None):
  r=pin(ROOT/'reports'/st/'report.json',h);assert r['all_passed']
  for p,z in r.get('input_sha256',{}).items():pin(ROOT/p,z,False)
  if 'code_sha256'in r:
   cs=[p for p in (ROOT/'code').glob('*.py')if p.stem==st.lower()]
   # Producer capitalization differs; the exact code bytes are pinned below.
  en=next(iter(r.get('evidence_sha256',{})),None)
  return r,pin(ROOT/'reports'/st/en,r['evidence_sha256'][en])if en else None
 fr,S=packet(FST,'cb03c3d5cfbc9cba86f9c15280b303a032dafa08fcddaad9ae782168171b13f5')
 pr,P=packet(PST,'e9dbeecb12be941d83a97b701ce7b723c6db06200bdc698d8233300d5b6e11bd')
 br,BP=packet(BST,'62d66bc4406244a863719b5e0ec9cdb2b780c5b2d61475c39ecc6796512f1b4f')
 basefr,BS=packet('AMBIENT_B18_NATIVE_ONE549_FRAME_SOURCE_20261105_v63','6fe41836acc116ca03338d59ca91515351ce61f61fd9d66b27bb68538892f66b')
 inherited=pin(ROOT/'reports/AMBIENT_B18_NATIVE_ONE549_SOURCE_POSITIVE_20261105_v63_REVIEW/report.json','b02f0d2118aab1499477d6c05603e8e0a7f8738b83cc702795fd96224a234e18');assert inherited['all_passed']
 hr,H=packet('AMBIENT_B18_NATIVE_ONE554_COMMON_DEN_HEIGHT_20261105_v63')
 hs,_=packet('AMBIENT_B18_NATIVE_ONE554_HEIGHT_SCOPE_20261105_v63')
 for st,r,fn in [(FST,fr,'ambient_B18_native_one554_frame_source_20261105_v63.py'),(PST,pr,'ambient_B18_native_one554_positive_input_20261105_v63.py'),('height',hr,'ambient_B18_native_one554_common_den_height_20261105_v63.py')]:pin(ROOT/'code'/fn,r['code_sha256'],False)
 src=P['full_primitive_input_J18_N_G_nu_GammaD'];assert src==BP['native_primitive_input_from_frozen_source']
 for a,b in [('basis18','ordinary_basis_J18'),('N_known_positions','N_known_J18_positions'),('N_known12_by12_by9_by4','N_known12_by12_by9_by4_cached_operands_ONLY'),('G8_C72','already_accepted_G8_J18_C72_images_COPIED_NOT_REEVALUATED'),('nu_two_J18','already_accepted_nu_two_J18_images'),('Gamma_two56','Gamma_two_current56'),('D_two_by10_by56','D_two_by10_by56')]:assert src[a]==BS[b]
 assert src['Gamma_two56']==S['current_Gamma_two56']and src['D_two_by10_by56']==S['current_D_two_by10_by56']
 assert P['parameterized_four_native_source_records']==S['parameterized_four_native_source_records']
 assert P['reused68_positive_covectorL_rows']==BP['positive_covectorL_row_cache']
 def source(suffix):return ROOT/next(p for p in pins if p.endswith(suffix))
 nativeop=field(source('actual-complete-Res4-syzygy-operator.json'),'native_S12_by32')
 mu=field(source('actual-T2-mu-input-source-and-cache-gate.json'),'native_S0_columns')
 a=field(source('actual-known-bottom-first-rank-and-corrected-source-kernel.json'),'all_zero_bottom_and_Y0_source_combinations')[554]
 combo=[[289,'25/4'],[379,'-1'],[643,'29/4'],[645,'25'],[651,'-25/2'],[670,'1']]
 assert a['original768_combination']==S['literal_full_original768_source']==combo
 assert a['actual_a1_27_times9']==S['all29_a1_kept_UNEVALUATED']==P['all29_a1_kept_UNEVALUATED']
 assert len(S['all29_a1_kept_UNEVALUATED'])==29
 assert [[i,q]for i,q in a['actual_a1_27_times9']if i%9 in[7,8]]==S['only_plane_a1_pivot5_original644']==[[52,'-175/3'],[53,'-377/60']]
 assert S['D_ordinary_pair_order10']==BS['D_ordinary_pair_order10']==[[0,7],[0,8],[1,7],[1,8],[4,7],[4,8],[5,7],[5,8],[7,7],[7,8]]
 for rec,(j,w)in zip(S['parameterized_four_native_source_records'],combo[2:]):
  c,im=divmod(j-640,4)
  assert(rec['original768'],rec['native128'],rec['S32_source'],rec['input_trace_mask'],rec['constant_source_coefficient'])==(j,j-640,c,im,w)
  terms=[['G8',row,g,sgn]for row,col,g,sgn in nativeop['koszul28_in_first8_rows']if col==c]
  assert terms==rec['intrinsic_source_entries']and len(terms)==2
  assert mu[j-640]==rec['full_point_S0_column48']==[]
 assert mu[4]==S['outer_pivot644_mask0_FULL_S0_column_nonzero']==[[11,'-1/2']]
 vr=pin(source('AMBIENT_B18_NATIVE92_ROOT_DIRECT_REVIEW_20261105_v63/report.json'),'10844d8720465054d6b9f7a7aa408dd95dd613a2e857eb52e50cd5a71a847e92');assert vr['all_passed']
 v=pin(source('actual92-Gamma-nativepivot23-plane2-and322-full-product-pieces.json'))
 assert v['original768_native23_sources'].index(644)==1
 outer=v['Gamma_Xp2_native23_twochannels_rr_rell'][1]
 assert outer==S['outer_full_GammaXp644_rr_rell_actual_cached_twochannels']==P['actual644_GammaXp_rr_rell_twochannels_UNPAIRED']and any(Q(q)for row in outer for q in row)
 # Before any independent envelope arithmetic: charge all decoded labels,
 # duplicated input arrays, future rational copies and dense D/U cells.
 parsed=count([fr,S,pr,P,br,BP,basefr,BS,inherited,hr,H,hs,src,nativeop,mu,a,combo,v,vr,pins])
 prior=parsed+2*count(src)+68*48*3+18*48*3+18000
 assert prior<150000 and raw+2_000_000<32*1024*1024
 save(out/'positive-before-independent-numeric-paths.json',{'all_parsed_numeric_labels_count':parsed,'whole_owned_upper':prior,'RAW_input_plus_output_upper':raw+2_000_000,'committed_integer_hard':384,'temporary_integer_hard':1153,'method':'two-pass path LCM then reverse positive scaled-numerator sum; no signed responses','no_new_source_physical_IFT_gauge_jet':True})
 def before(*v):
  nonlocal maxpre
  b=max(v,default=0);maxpre=max(maxpre,b)
  if b>1153:save(out/'guard-prefix.json',{'cursor':cursor,'prebits':b,'no_capraise':True});raise RuntimeError('1153 temporary STOP')
 def ck(x):
  nonlocal maxcommit
  b=abs(x).bit_length();maxcommit=max(maxcommit,b)
  if b>384:save(out/'guard-prefix.json',{'cursor':cursor,'committed_bits':b,'no_capraise':True});raise RuntimeError('384 committed STOP')
  return x
 def mul(x,y):before(abs(x).bit_length()+abs(y).bit_length());return ck(x*y)
 def add(x,y):before(max(abs(x).bit_length(),abs(y).bit_length())+1);return ck(x+y)
 def lcm(x,y):return mul(x//math.gcd(x,y),y)
 def q(x):z=Q(x);ck(z.numerator);ck(z.denominator);return z
 G=[[q(x)for x in row]for row in src['G8_C72']];nu=[[q(x)for x in row]for row in src['nu_two_J18']]
 N=[[[[q(x)for x in masks]for masks in bases]for bases in row]for row in src['N_known12_by12_by9_by4']]
 basis=list(map(tuple,src['basis18']));bi={m:i for i,m in enumerate(basis)};known=src['N_known_positions'];assert known==[0,1,2,3,4,5,6,15,16]
 # The inherited independent monomial quotient proof covers all660 NFs and
 #1089 products; additionally recheck the exact monic-pair support here.
 for x,y in itertools.product(nu,repeat=2):
  assert all(tuple(a+b for a,b in zip(basis[i],basis[j]))not in bi for i,p in enumerate(x)if p for j,z in enumerate(y)if z)
 Gamma=[{i:q(z)for i,z in row if i<48 and Q(z)}for row in src['Gamma_two56']]
 DD=[[{i:q(z)for i,z in enumerate(row)if i<48 and Q(z)}for row in ch]for ch in src['D_two_by10_by56']]
 # Independently evaluate the full point Res4(S) columns, using actual
 # trace collisions and signs, without any functional contraction.
 for rec in S['parameterized_four_native_source_records']:
  point=[Q(0)]*48;im=rec['input_trace_mask']
  for _,r,g,sgn in rec['intrinsic_source_entries']:
   for mask,z in enumerate(G[g][:4]):
    collision=mask&im;z*=sgn
    if collision==3:z*=nu[0][0]*nu[1][0]
    elif collision:z*=nu[collision-1][0]
    point[4*r+(mask^im)]+=z
  assert not any(point)
 # Raw products retain their original denominators, hence match D divisibility
 # certificates as well as rational absolute bounds. No cancellation/reduction.
 def product_path(xs):
  n=d=1
  for z in xs:
   if not z:return None
   n=mul(n,abs(z.numerator));d=mul(d,z.denominator)
  return n,d
 def accumulate(factory,size):
  ds=[1]*size;us=[0]*size
  for i,n,d in factory():ds[i]=lcm(ds[i],d)
  for i,n,d in factory():assert ds[i]%d==0;us[i]=add(us[i],mul(n,ds[i]//d))
  return [(d,u)if u else(1,0)for d,u in zip(ds,us)]
 def rpaths(key):
  kind,ch,pair,target=key;cov=Gamma[ch]if kind==0 else DD[ch][pair]
  for outrow,z in sorted(cov.items(),reverse=True):
   physical,om=divmod(outrow,4)
   for im in reversed(range(4)):
    mask=om^im;collision=mask&im
    if collision==3:continue
    for r in reversed(range(12)):
     for kp in reversed(range(9)):
      x=N[r][physical][kp][mask];nb=known[kp]
      if not x:continue
      if collision:
       for vb,y in reversed(list(enumerate(nu[collision-1]))):
        if y and tuple(a+b for a,b in zip(basis[nb],basis[vb]))==basis[target]:
         n,d=product_path((z,x,y));yield 4*r+im,n,d
      elif nb==target:
       n,d=product_path((z,x));yield 4*r+im,n,d
 expectedrows={tuple(k):{i:(int(d),int(u))for i,d,u in row}for k,row in P['reused68_positive_covectorL_rows']}
 rows={};rowmax=0
 for key,e in sorted(expectedrows.items(),reverse=True):
  cursor={'phase':'all68_numeric_L_rows','key':key};got=accumulate(lambda:rpaths(key),48)
  assert got==[e.get(i,(1,0))for i in range(48)],('L',key)
  rows[key]=got;rowmax=max(rowmax,max(max(d.bit_length(),u.bit_length())for d,u in got))
 assert len(rows)==68
 def spaths():
  for rec in reversed(S['parameterized_four_native_source_records']):
   w=q(rec['constant_source_coefficient']);im=rec['input_trace_mask']
   for _,r,g,sgn in reversed(rec['intrinsic_source_entries']):
    for base in reversed(range(18)):
     for mask in reversed(range(4)):
      x=G[g][4*base+mask]
      if not x:continue
      collision=mask&im;om=mask^im
      if not collision:
       n,d=product_path((w,x));yield base*48+4*r+om,n,d
      elif collision in[1,2]:
       for vb,y in reversed(list(enumerate(nu[collision-1]))):
        dest=bi.get(tuple(a+b for a,b in zip(basis[base],basis[vb])))
        if y and dest is not None:n,d=product_path((w,x,y));yield dest*48+4*r+om,n,d
 SS=accumulate(spaths,18*48)
 for k,row in enumerate(P['positive_weighted_intrinsic_S18']):
  e={i:(int(d),int(u))for i,d,u in row};assert SS[k*48:(k+1)*48]==[e.get(i,(1,0))for i in range(48)]
 assert SS[:48]==[(1,0)]*48
 newmax=max(max(d.bit_length(),u.bit_length())for d,u in SS)
 def sum_du(terms):
  def gen():
   for n,d in terms:yield 0,n,d
  return accumulate(gen,1)[0]
 def dot(row,base):
  return sum_du([(mul(u,SS[48*base+i][1]),mul(d,SS[48*base+i][0]))for i,(d,u)in reversed(list(enumerate(row)))if u and SS[48*base+i][1]])
 vf={0:1,1:2,4:3,5:4,7:5,8:6};cub=[(0,17),(5,16),(6,15),(15,6),(16,5)];fine=[];totals=[];absouter=[]
 for ch in range(2):
  channel=[]
  for j,(rb,sb)in enumerate(cub):
   z=dot(rows[(0,ch,-1,rb)],sb);assert list(map(str,z))==P['positive_Gamma_weighted_cubic_two_by5'][ch][j];channel.append(z)
  for p,(a,b)in enumerate(S['D_ordinary_pair_order10']):
   ts=[(0,7+p)]+[(vf[x],vf[y])for x,y in([(a,b)]if a==b else[(a,b),(b,a)])]
   for j,(rb,sb)in enumerate(ts):
    z=dot(rows[(1,ch,p,rb)],sb);assert list(map(str,z))==P['positive_D_weighted_second_two_by10_unique_first_decomp'][ch][p][j];channel.append(z)
  assert len(channel)==34;fine+=channel;oo=[]
  for weight,x in zip(['-377/60','-175/3'],outer[ch]):
   n,d=product_path((q(weight),q(x)));oo.append((d,n))
  assert[list(map(str,z))for z in oo]==P['positive_outer644_two_by2'][ch];absouter.append(oo)
  total=sum_du([(u,d)for d,u in channel+oo if u]);assert list(map(str,total))==P['positive_native_plus_absouter_two'][ch];totals.append(total)
  newmax=max(newmax,max(max(d.bit_length(),u.bit_length())for d,u in channel+oo+[total]))
 assert len(fine)==68 and rowmax==145 and newmax==111
 # Rebuild an independent common denominator and positive direct-path prior.
 classes={'known_partial_N9':src['N_known12_by12_by9_by4'],'canonical_G8':src['G8_C72'],'monics':src['nu_two_J18'],
  'Gamma_and_D_scalar_coefficients':[z for row in src['Gamma_two56']for i,z in row]+[z for ch in src['D_two_by10_by56']for row in ch for z in row],
  'native_c0_constants':[r['constant_source_coefficient']for r in S['parameterized_four_native_source_records']]}
 def leaves(x):
  if isinstance(x,list):
   for y in x:yield from leaves(y)
  else:yield q(x)
 info={}
 for name,x in classes.items():
  d=1;n=0
  for z in leaves(x):d=lcm(d,z.denominator);n+=1
  u=d
  for z in leaves(x):u=max(u,mul(abs(z.numerator),d//z.denominator))
  info[name]={'D':str(d),'U_including_scalar_one':str(u),'D_bits':d.bit_length(),'U_bits':u.bit_length(),'operands':n}
 assert info==H['scalar_operand_classes']
 d=u=1
 for x in info.values():d=mul(d,int(x['D']));u=mul(u,int(x['U_including_scalar_one']))
 u=mul(u,1<<28);paths=4*2*12*4*4*9*18*(1+18)*11
 assert paths==52005888==H['explicit_native_direct_paths_upper']and paths<1<<28
 assert str(d)==H['common_D']and str(u)==H['common_absolute_scaled_U_with_path2pow28']
 assert max(d.bit_length(),u.bit_length())==339 and 3*339+1==1018
 primitivebits=max(max(abs(z.numerator).bit_length(),z.denominator.bit_length())for x in classes.values()for z in leaves(x))
 assert primitivebits==hs['existing_primitive_N_G_nu_GammaD_operand_bits152']==152
 assert hs['reused68_operator_row_D_U_bits145']==rowmax and hs['new_WEIGHTED_S_and_contraction_positive_operations_bits111']==newmax
 assert maxcommit<=384 and maxpre<=1153
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 evidence={'all_passed':True,'source_frame_full6_c0_and29_a1_current_Gamma_D_verified':True,'same_J18_N_source_by_physical_transpose12_and_canonical_G8_nu_inputs_bound':True,
  'all_four_full_masked_S0_columns_recomputed_zero':True,'outer644_mask0_full_S0_nonzero':[ [11,'-1/2'] ],'individual_L2S0_retained':True,
  'independent_method':'Reverse physical-output and source loop order, two-pass integer path denominators LCM then positive scaled-numerator sum. No producer import, positive Fraction summation, or signed response.',
  'all_L_rows':68,'all_L_D_U_cells':3264,'all_S_D_U_cells':864,'all_fine_native_D_U_components':68,'all_outer_ABS_products':4,
  'all_exact_D_U_pairs_match_not_only_rational_quotients':True,'all_monic_pair_supports_zero_independently_rechecked':True,
  'numeric_row_envelope_bits':rowmax,'new_S_fine_outer_total_bits':newmax,'primitive_operand_bits':primitivebits,
  'direct_path_classes':info,'direct_paths_upper':paths,'common_D':str(d),'common_positive_scaled_U':str(u),'native_C72_committed_prior339_temporary1018':True,
  'proper_partial_paths_covered_by_scalar_one_in_each_U':True,'source_unknown_N2N3':'Remain formal. Full weighted S0=0 removes their whole kernel contractions; nonplane unknown terms have z factors in the independently proved J18 monomial quotient. No coefficient of N2/N3 was assigned zero.',
  'scope':'Native SOURCE/envelope only. OLD/FIRST/Q assembly and VALUE code/custody are separate. No signed native/outer/whole upper response, HF, rank, or case554zero/nonzero claim.',
  'whole_owned_upper':prior,'RAW_input_plus_output_upper':raw+2000000,'actual_max_committed_bits':maxcommit,'actual_max_temporary_precheck_bits':maxpre}
 ep=out/'independent-all3264-L864-S68-fine-numeric-source-envelope-proof.json';save(ep,evidence)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_POSITIVE_ONLY':True,'input_sha256':pins,'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},'whole_ONE554_VALUE_GO_claimed':False,'seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':prior,'RAW':raw+2000000,'rowbits':rowmax,'newbits':newmax,'seconds':report['seconds']}))
if __name__=='__main__':main()
