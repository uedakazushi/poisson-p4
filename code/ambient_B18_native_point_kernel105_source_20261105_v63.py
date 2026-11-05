#!/usr/bin/env python3
"""Explicit constant native point kernel and finite-jet determinacy source."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_NATIVE_POINT_KERNEL105_SOURCE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};maximum=0
 def read(p):q=ROOT/p;pins[p]=sha(q);return json.loads(q.read_text())
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 p='reports/AMBIENT_B18_NATIVE_CANONICAL68_VALUES_20261105_v63/actual-G8-C68-and-monic-second-values.json';g=read(p)
 op=read('reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json')
 assert read('reports/AMBIENT_B18_NATIVE_CANONICAL68_DIRECT_REVIEW_20261105_v63/report.json')['all_passed']
 assert read('reports/AMBIENT_B18_NATIVE_ONE549_SOURCE_POSITIVE_20261105_v63_REVIEW/report.json')['all_passed']
 def ck(q):
  nonlocal maximum
  maximum=max(maximum,abs(q.numerator).bit_length(),q.denominator.bit_length());assert maximum<=32;return q
 c=[]
 for v in g['full_G8_C68']:assert v[:3]==['0','0','0'];c.append(F(v[3]))
 assert c==list(map(F,['1/2','25/8','0','0','0','5','-125/2','0'])) and c[0]==F(1,2)
 pairs=[(i,j) for i in range(8) for j in range(i+1,8)];pairindex={p:i for i,p in enumerate(pairs)}
 want=[]
 for k,(i,j) in enumerate(pairs):want.extend([[i,k,j,1],[j,k,i,-1]])
 assert op['native_S12_by32']['koszul28_in_first8_rows']==want
 assert op['native_S12_by32']['free4_zero_source_columns']==[[8,28,1],[9,29,1],[10,30,1],[11,31,1]]
 # Prior32 coversconstantfractions(ci/c0)andatmosttwoinputconstantproducts;
 # no function responses or point-rank search is performed.
 basis=[]
 for k in range(28):
  for mask in [1,2,3]:basis.append({'kind':'nil_input_mask','free_native_source':4*k+mask,'column':[[4*k+mask,'1']]})
 for i in range(1,8):
  for j in range(i+1,8):
   v={4*pairindex[i,j]:F(1)}
   a=ck(-c[i]/c[0]);b=ck(c[j]/c[0])
   if a:v[4*pairindex[0,j]]=a
   if b:v[4*pairindex[0,i]]=b
   # FixedQQ Koszulidentitycheck, no source jet or Gaussian elimination.
   image=[F(0)]*8
   for src,q in v.items():
    k,mask=divmod(src,4);assert mask==0;a0,b0=pairs[k]
    image[a0]=ck(image[a0]+ck(c[b0]*q));image[b0]=ck(image[b0]-ck(c[a0]*q))
   assert image==[F(0)]*8
   basis.append({'kind':'scalar_mask_Koszul_kernel','free_native_source':4*pairindex[i,j],
    'column':[[k,str(v)] for k,v in sorted(v.items())]})
 assert len(basis)==105
 free=[v['free_native_source'] for v in basis];units=[4*pairindex[0,j] for j in range(1,8)]+list(range(112,128))
 assert len(set(free))==105 and len(set(units))==23 and set(free).isdisjoint(units) and set(free)|set(units)==set(range(128))
 rows=[4*j+3 for j in range(1,8)]+list(range(32,48));det=ck((-c[0])**7);assert det==F(-1,128)
 ep='constant-kernel105-basis-unit23-minor-and-native-functional-determinacy.json'
 e={'actual_C4_G8_point_equals_c_times_tauAtauB':list(map(str,c)),'point_monic_relations':'tauA²=tauB²=0; fromsamecanonicalbranch',
  'explicit105_kernel_basis_columns':basis,'fixed23_source_complement':units,'fixed23_output_rows':rows,'fixed23_minor_block':'diag(-c0 I7, I16)','fixed23_minor_determinant':str(det),
  'all21_scalar_mask_kernel_identities_checked':True,'all84_nonzero_mask_images_zero_by_tauAtauB_nilproduct':True,
  'completeness_proof':'EverykerS0 coefficienton105free sources determinesa linearcombinationofthese105columns; subtractit, remainderissupportedin23unit sources, andfixedminor forcesremainder0.',
  'native_X0_same_kernel':'Res4(N0^t)unit, acceptedfullnormalizerSOURCE; nooldB/nativecancellationinference.',
  'finite_native_functional_determinacy_theorem':{
   'hypotheses':'constant cN inker(fullRes4 S0), regularR0coefficients, fixedΓ/DwithDonlythe10nonzzpairs andΓcubicr²ell, samevaryingmonictracebasis.',
   'Dsecond':'D[L0*S2+alluniqueordinaryLfirst*Sfirst]; unknownL2*S0cN=0.',
   'Gamma_cubic':'Γ[L0*S3r²ell+Lr*S2rell+Lell*S2rr+Lrr*Sell+Lrell*Sr]; unknownL3*S0cN=0 andnonplaneL2 cannotenterpureplane.',
   'data_needed':'CanonicalG8/nu imageJ18_with_r2ell plusN0,activefirst,ell,rr/rell only. Allalreadyfinitecache-bound; Nmixed2/N3remainformalbutprojectedfunctionalindependent.',
   'no_annihilation_claim':'Thisprovesdetermination, NOTthatanyofthe105functionalvaluesarezero; fulla1/OLD/FIRST/outercorrectionsremainseparate.'},
  'no_native_functional_values_or_605_response_job':True,'source_prior32_temporary97':True,'actual_constant_QQ_bits':maximum}
 dump(ep,e);r={'id':STEM,'all_passed':True,'pure_constant_source_and_determinacy_only':True,'no_point_rank_search_or_functional_value':True,
  'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',r);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'constant_bits':maximum,'seconds':r['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
