#!/usr/bin/env python3
"""Exact native Rell4 cache/source input, no adapter reconstruction yet."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_first_old_cache_20261105_v63 import mul
from ambient_B18_ell2_raw33_quotient_20261105_v63 import mv,srec
from ambient_B18_native128_first_20261104_v63 import PATHS as OLD_PATHS
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 raw=read('reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json');basis=raw['actual33_basis'];actions=[list(map(H.parse,cols)) for cols in raw['coordinate_actions9_sparse_columns']]
 images=[]
 for s,r,e in basis:
  v=[F(0)]*4
  if s==0 and r<=1 and e<=1:v[r+2*e]=F(1)
  images.append(v)
 def project(v):return [sum((c*images[i][j] for i,c in v.items()),F(0)) for j in range(4)]
 for i,(s,r,e) in enumerate(basis):
  for j in range(33):
   v={j:F(1)}
   for k in ([s-1] if s else [])+[7]*r+[8]*e:v=mv(actions[k],v)
   assert project(v)==mul(images[i],images[j])
 for k,action in enumerate(actions):
  var=[F(0)]*4
  if k==7:var[1]=F(1)
  if k==8:var[2]=F(1)
  for j,v in enumerate(action):assert project(v)==mul(var,images[j])
 mbase='reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/'
 ms=read(mbase+'report.json');peer=read('reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json');assert ms['all_passed'] and peer['all_passed']
 gam=read(mbase+'phase2-canonical30-exact-values.json');proper=read(mbase+'phase3-actual-kappa9-and-normalized-nilframe.json');nu=read(mbase+'phase1-source150-readback-fingerprints.json')['nu_A_B33']
 assert gam['all30_constraints_exact'] and gam['all30_gauge_B11_projections_exact'] and proper['actual_same_CA_frame_and_fiber_constraints_exact']
 def scalar(p):return project({i:F(c) for i,c in enumerate(p) if F(c)})
 def C132(p):
  z=[F(0)]*16
  for i,c in enumerate(map(F,p)):
   if c:
    b,mask=divmod(i,4)
    for j,a in enumerate(images[b]):z[4*j+mask]+=a*c
  return z
 ev=list(map(scalar,gam['actual_even15_B33']));od=list(map(scalar,gam['actual_odd15_B33']));actual=list(map(C132,proper['actual14_C132']));kap=list(map(C132,proper['kappa9_C132']));nu=list(map(scalar,nu))
 src=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json');phi={i:project(H.parse(p)) for i,p in zip(src['joint105_unknown_order'],src['all105_unknown_final_values'])}
 b=read('reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json');bm=[tuple(x) for x in b['basis56']];indices=[]
 for r,e in [(0,0),(1,0),(0,1),(1,1)]:v=[0]*9;v[7]=r;v[8]=e;indices.append(bm.index(tuple(v)))
 assert all(phi[i]==[H.parse(p).get(k,F(0)) for k in indices] for i,p in zip(b['unknown105_order'],b['all105_cached_values_B56']))
 firstell=read('reports/AMBIENT_B18_ELL2_NATIVE_GAUGE_FIRST_ELL_20261105_v63/actual30-firstell-and-native-S1ell.json');assert [v[2] for v in ev+od]==list(map(F,firstell['actual_gauge30_first_ell']));assert [v[2] for v in nu]==list(map(F,firstell['actual_nuA_nuB_first_ell']))
 oldpath='reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/phase3-actual-proper-G8-B11.json';old=read(oldpath)
 assert [v[:8] for v in actual]==[list(map(F,v[:8])) for v in old['actual14_C_B44']]
 assert [v[:8] for v in kap]==[list(map(F,v[:8])) for v in old['kappa9_C_B44']]
 physical=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json');F12=[]
 for rec in physical['physical12_coefficients']:
  d=project(H.parse(rec['delta']));a=project(H.parse(rec['bA']));b=project(H.parse(rec['bB']));p=[F(0)]*16
  for j in range(4):p[4*j]=d[j];p[4*j+1]=a[j];p[4*j+2]=-b[j]
  F12.append(p)
 assert [r['physical_row'] for r in physical['physical12_coefficients']]==[18,19,28,29,60,67,72,63,64,68,73,74]
 assert all(v[0]==v[1]==v[2]==0 for v in nu),'monicsstartsecondonRell4'
 refs={k:OLD_PATHS[k] for k in ['source','interface','proper','equations','first5','second5','membership','preparation','normalizer','Q_operator']}
 for p in refs.values():hs[p]=H.sha(ROOT/p)
 slots={'cached_selected_input_coefficients':105*4+30*4+14*16+9*16+12*16+2*4,
  'source33raw_and_MS_full_input_arrays_upper':15000,
  'full17_telescope_core5_R12_T_N_matrix_and_scratch':16*(17*12+17*5+12*5+5*5+7*5+7*5+6*12*12),
  'all619_native_source_folds_and_remainingvalues':10000,'output_adapter_and_intrinsic_S2_upper':16*12*12+48*128,
  'source_equations_refs_and_metadata':10000}
 assert sum(slots.values())<150000
 evidence={'Rell4_basis':['1','rA','ell','rA*ell'],'raw33_to_Rell4_columns':[[str(c) for c in v] for v in images],
  'all1089_raw33_products_and297_coordinate_intertwiners_verified':True,
  'actual105_Rell4':[[i,list(map(str,v))] for i,v in sorted(phi.items())],
  'actual_even15_Rell4':[list(map(str,v)) for v in ev],'actual_odd15_Rell4':[list(map(str,v)) for v in od],
  'actual14_C16':[list(map(str,v)) for v in actual],'actual_kappa9_C16':[list(map(str,v)) for v in kap],
  'actual_physicalF12_C16':[list(map(str,v)) for v in F12],'actual_nu_A_B_Rell4':[list(map(str,v)) for v in nu],
  'all105_B56_Rell4_overlap_exact':True,'all30_firstell_overlap_exact':True,'all14_kappa9_constant_rfirst_overlap_exact':True,
  'normalizer_source_references':refs,
  'worker_required_checks':['full17offintersection','fullfirst5monomialtelescopes','dH12=R12F12','TH12=(G8,0)','NF12=d(G8,0)','N0/rfirstcachedmatch','G8firstellpeeredmatch','all619nativeS0zero','all619L1rOmegaEll_and_L1ellOmegaR_retained'],
  'native_formula':'GammaL0*S2rell*cN +GammaL1r*OmegaEll+GammaL1ell*OmegaR, plusalreadyfrozenoldknownparts(includingnegativeSp1a1)',
  'actual_N1ell_S2rell_or_nativevalues_not_computed':True,
  'whole_slot_bound':sum(slots.values()),'slot_breakdown':slots,
  'limits':{'seconds':60,'memory_bytes':2_000_000_000,'raw_artifact_bytes':32*1024*1024,'QQ_coefficients':150000,'bits':384},
  'no_new_physical_implicit_gauge_or_higher_jet':True}
 ep=out/'actual-native-Rell4-cache-source-input-and-adapter-gate.json';H.dump(ep,evidence)
 report={'id':STEM,'all_passed':True,'actual_native_adapter_or_secondvalues_started':False,'whole_slot_bound':sum(slots.values()),
  'all_source33_and_B56_and_firstell_frames_overlap':True,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),
  'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','whole_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
