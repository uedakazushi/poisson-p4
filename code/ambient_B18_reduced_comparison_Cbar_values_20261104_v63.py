#!/usr/bin/env python3
"""One bounded existing-cache value/rank job; no source or gauge worker."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_REDUCED_COMPARISON_CBAR_VALUES_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,list):return [enc(y)for y in x]
 if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
 return x
def main(out):
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
 files={
 'B':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
 'coords':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json',
 'cache_peer':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63_REVIEW/report.json',
 'Cbar':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63/actual-Cbar-minimal-resolution-and-three-targets.json',
 'Cbar_peer':'reports/AMBIENT_B18_CBAR_LITERAL_RESOLUTION_20261104_v63_REVIEW/report.json',
 'nil2_peer':'reports/AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW/report.json',
 'origin':'reports/AMBIENT_B18_REDUCED_COMPARISON_ORIGIN_20261104_v63/report.json',
 }
 expected={
 'B':'d63f10af5050e2437b5df135a462a664c31dc7dccd771d376136ba2d841a2427',
 'coords':'c8de9c27b4d881bb236081a3558ef9408668357fd38f4a28fe3434d079207ac3',
 'cache_peer':'bb6e7e7c5327528f2be31b1d29db5891461d12adb319daa184b6c6addd63d693',
 'Cbar':'f88c4d8e177e7d00010205e21b327ad2af36c4539a2c12f81c59341608800159',
 'Cbar_peer':'7e17f85c19492d3e81dc502b38f33819e6045045bb63b131d160019adec3e378',
 'nil2_peer':'02888d3c4b987327f6df2d8faf1f2adffa3f9c01af47e27aff8c5c0982e75472',
 'origin':'0524a73bc87ff7c12ffa04ccd5137aca5d0f15173a725e76f50dbde3f6229674',
 }
 for k,p in files.items():assert sha(ROOT/p)==expected[k],k
 for k in ['cache_peer','Cbar_peer','nil2_peer','origin']:assert read(ROOT/files[k])['all_passed']
 B,C,L=[read(ROOT/files[k])for k in ['B','coords','Cbar']]
 assert B['B11_variable_order8']==L['actual_variables']==['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA']
 assert B['B11_exact_ideal']=='(rA^3, all z_i*z_j, all z_i*rA^2, zA2,zA3,zA6)'
 assert B['C_B44_basis']=='B11_basis in given order, each multipliedby(1,tauA,tauB,tauA*tauB)'
 assert L['actual_three_targets']==['tauA','tauB','tauA*tauB'] and L['r_lengths']==[3,3,2]
 mons=[tuple(x)for x in B['B11_basis_monomials8']];mi={m:i for i,m in enumerate(mons)}
 shifts=list(map(Q,L['CB_shift_coefficients7']));assert [i for i,c in enumerate(shifts)if c]==[0,4,5]
 assert len(C['kappa9_C_B44'])==9 and all(len(v)==44 for v in C['kappa9_C_B44'])
 # Source/annihilator metadata is complete before this finite value worker starts.
 pre={'all9_actual_proper_coordinates_saved_exact_C44':True,'B11_basis_count':11,'C44_dimension':44,
 'missing_z2_z3_z6_annihilate_ALL_three_sectors':True,
 'all_B11_defining_relations_in_three_sector_annihilators':True,
 'actual_nil_generator_count':2,'basis_multiples':4,'new_image_vectors':8,
 'max_image_closure_columns':8*(1+8+36),'no_unknown_coefficients_set_zero':True,
 'new_source_gauge_jet_or_point':False,'visible_term_product_upper_bound':50000,'QQ_slot_upper_bound':30000,
 'max_output_bytes':32000000}
 (out/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
 def bmul(x,y):
  z=[Q(0)]*11
  for i,a in enumerate(x):
   if not a:continue
   for j,b in enumerate(y):
    if not b:continue
    m=tuple(u+v for u,v in zip(mons[i],mons[j]));k=mi.get(m)
    if k is not None:z[k]+=a*b
  return z
 nu=[]
 for sp in B['nu_A_B_actual_B11']:
  v=[Q(0)]*11
  for i,c in sp:v[i]=Q(c)
  nu.append(v)
 def cmul(x,y):
  z=[Q(0)]*44
  for a,b in itertools.product(range(4),repeat=2):
   v=bmul(x[a::4],y[b::4])
   if a&b&1:v=bmul(v,nu[0])
   if a&b&2:v=bmul(v,nu[1])
   for j,c in enumerate(v):z[4*j+(a^b)]+=c
  return z
 basis44=[]
 for i in range(44):v=[Q(0)]*44;v[i]=Q(1);basis44.append(v)
 for i,j in itertools.product(range(44),repeat=2):
  want=[Q(0)]*44
  for k,c in B['C_B44_multiplication_basis_products'][i][j]:want[k]=Q(c)
  assert cmul(basis44[i],basis44[j])==want
 # Actual B11 -> three cyclic R0 quotient ring maps, not evaluation at new points.
 sector_maps=[]
 for sec,N in enumerate([3,3,2]):
  cols=[]
  for m in mons:
   c=Q(1);degree=m[7]
   for i in range(7):
    if m[i]:
     c*=(-shifts[i])**m[i]if sec==1 else Q(0)
     degree+=2*m[i]
   v=[Q(0)]*N
   if degree<N:v[degree]=c
   cols.append(v)
  sector_maps.append(cols)
 def proj(v,sec):return [sum(v[j]*sector_maps[sec][j][k]for j in range(11))for k in range([3,3,2][sec])]
 def rmul(x,y):
  z=[Q(0)]*len(x)
  for i,a in enumerate(x):
   for j,b in enumerate(y):
    if i+j<len(z):z[i+j]+=a*b
  return z
 for i,j,sec in itertools.product(range(11),range(11),range(3)):
  x=[Q(0)]*11;y=x.copy();x[i]=Q(1);y[j]=Q(1)
  assert proj(bmul(x,y),sec)==rmul(proj(x,sec),proj(y,sec))
 def add(*vs):return [sum(x)for x in zip(*vs)]
 def scale(v,c):return [Q(c)*x for x in v]
 kappa=[list(map(Q,v))for v in C['kappa9_C_B44']]
 a,d,x,t=[kappa[i]for i in [0,1,2,6]]
 h=add(a,scale(d,-1),x,scale(t,-2));cc=scale(add(d,x),Q(1,2));ll=scale(add(d,scale(x,-1),a),Q(1,3))
 nilgens=[cmul(d,h),cmul(cmul(d,cc),ll)]
 def to_cbar(v):return proj(v[1::4],0)+proj(v[2::4],1)+proj(v[3::4],2)
 images=[];Cvalues=[]
 for alpha,gamma in itertools.product(range(2),range(4)):
  v=cmul(nilgens[alpha],basis44[gamma]);Cvalues.append(v);images.append(to_cbar(v))
 assert all(v[i]==0 for v in images for i in [0,3,6])
 # Literal Cbar coordinate actions from its exact three source sectors.
 lens=[3,3,2];offsets=[0,3,6]
 actions=[]
 for var in range(8):
  cols=[]
  for j in range(8):
   v=[Q(0)]*8
   for sec,(N,off)in enumerate(zip(lens,offsets)):
    if off<=j<off+N:
     deg=j-off+(1 if var==7 else 2)
     coef=Q(1)if var==7 else (-shifts[var]if sec==1 else Q(0))
     if deg<N:v[off+deg]=coef
   cols.append(v)
  actions.append(cols)
 def act(cols,v):return [sum(cols[j][i]*v[j]for j in range(len(v)))for i in range(len(cols[0]))]
 closure=[];labels=[]
 for j,v in enumerate(images):closure.append(v);labels.append([j,[]])
 for j,v in enumerate(images):
  for i in range(8):closure.append(act(actions[i],v));labels.append([j,[i]])
  for i,k in itertools.combinations_with_replacement(range(8),2):closure.append(act(actions[i],act(actions[k],v)));labels.append([j,[i,k]])
 assert len(closure)==360
 normalized={};selected=[]
 def reduce(v):
  v=v.copy()
  for pivot,b in sorted(normalized.items()):
   if v[pivot]:q=v[pivot];v=[x-q*y for x,y in zip(v,b)]
  return v
 for j,v in enumerate(closure):
  w=reduce(v)
  if any(w):
   p=next(i for i,c in enumerate(w)if c);normalized[p]=scale(w,1/w[p]);selected.append(j)
 assert all(not any(reduce(v))for v in closure)
 assert all(not any(reduce(act(A,b)))for A in actions for b in normalized.values())
 free=[i for i in range(8)if i not in normalized];q=len(free)
 def qmap(v):return [reduce(v)[i]for i in free]
 qactions=[]
 for A in actions:
  cols=[]
  for i in free:v=[Q(0)]*8;v[i]=Q(1);cols.append(qmap(act(A,v)))
  qactions.append(cols)
 assert all(qmap(act(A,b))==[Q(0)]*q for A in actions for b in normalized.values())
 for i,j,k in itertools.product(range(8),range(8),range(q)):
  e=[Q(0)]*q;e[k]=Q(1)
  assert act(qactions[i],act(qactions[j],e))==act(qactions[j],act(qactions[i],e))
 def rank(cols,n):
  piv={}
  for v in cols:
   w=v.copy()
   for p,b in sorted(piv.items()):
    if w[p]:c=w[p];w=[x-c*y for x,y in zip(w,b)]
   if any(w):p=next(i for i,c in enumerate(w)if c);piv[p]=scale(w,1/w[p])
  return len(piv)
 m1=[c for A in qactions for c in A]
 m2=[act(A,c)for A in qactions for c in m1]
 m3=[act(A,c)for A in qactions for c in m2]
 assert not any(any(v)for v in m3)
 hfun=[q-rank(m1,q),rank(m1,q)-rank(m2,q),rank(m2,q)];assert hfun[0]==3
 assert q==4 and free==[0,3,4,6] and all(not any(c)for A in qactions[:7]for c in A)
 lam=-qactions[7][0][2]/qactions[7][1][2];assert lam==Q(11,12)
 split=[[Q(1),lam,Q(0),Q(0)],[Q(0),Q(0),Q(0),Q(1)],[Q(0),Q(1),Q(0),Q(0)],[Q(0),Q(0),Q(1),Q(0)]]
 assert rank(split,4)==4
 assert act(qactions[7],split[0])==act(qactions[7],split[1])==[Q(0)]*4
 assert act(qactions[7],split[2])==split[3] and act(qactions[7],split[3])==[Q(0)]*4
 result={
 'all9_proper_C44_values_from_immutable_phase3':True,
 'no_new_source_or_gauge_solution':True,
 'B11_to_three_sector_ring_maps':sector_maps,
 'all363_sector_basis_product_identities':True,'all1936_C44_basis_products_match_source':True,
 'actual_nil_generators_C44':nilgens,'all8_nil_basis_multiple_values_C44':Cvalues,
 'all8_actual_nil_image_vectors_in_Cbar8':images,
 'Cbar8_literal_coordinate_actions_columns':actions,
 'all360_R0_submodule_closure_columns':closure,'all360_labels':labels,
 'image_pivots':sorted(normalized),'selected_original_closure_columns':selected,
 'normalized_triangular_image_basis':[[p,b]for p,b in sorted(normalized.items())],
 'image_rank':len(normalized),'image_invariant_under_all8_actions':True,
 'quotient_original_basis_indices':free,'quotient_length':q,'quotient_Hilbert_function':hfun,
 'all8_quotient_coordinate_actions_columns':qactions,'quotient_actions_commute':True,
 'exact_decomposition_basis_columns_in_quotient':split,
 'exact_decomposition':'k plus k plus k[[rA]]/(rA²), all zA0..6 act zero',
 'right_exact_cokernel_basechange_no_kernel_Tor_drop':True,
 'actual_Cprime_bar_is_not_generic_DVR_Cprime':True,
 'nil_D0_full_generators_KJK_or_cube_decided':False,
 }
 ep=out/'actual-eight-nil-images-and-reduced-Cbar-cokernel.json';ep.write_text(json.dumps(enc(result),indent=2)+'\n')
 assert ep.stat().st_size<32000000
 rels=list(files.values())+[f'research/notes/{STEM}.md',str(Path(__file__).relative_to(ROOT))]
 report={'id':STEM,'all_passed':True,'status':'PROVED_actual_reduced_comparison_ell0_cokernel_from_existing_B11',
 'actual_nil_image_rank':len(normalized),'actual_reduced_Cbar_length':q,'actual_reduced_Cbar_HF':hfun,
 'exact_reduced_Cbar_decomposition':'k plus k plus k[[rA]]/(rA²), all zA0..6 zero',
 'quotient_actions_and_NF_invariance_verified':True,
 'new_source_gauge_jet_point_GB_or_780_expansion':False,
 'nil_D0_generators_KJK_or_cube_decided':False,
 'input_sha256':{p:sha(ROOT/p)for p in rels},'evidence_sha256':{ep.name:sha(ep),'preflight.json':sha(out/'preflight.json')},
 'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k]for k in ['id','all_passed','actual_nil_image_rank','actual_reduced_Cbar_length','actual_reduced_Cbar_HF','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output.resolve())
