#!/usr/bin/env python3
"""Cached C72 frame/monic/Nfirst source bundle, no G8second or70 response."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_C72_FRAME_SOURCE_20261105_v63'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 base='reports/AMBIENT_B18_ELL2_T2_NATIVE70_A18_SOURCE_20261105_v63';r=read(base+'/report.json');assert r['all_passed'] and r['all70_fullC4_point_kernel_verified']
 e=read(base+'/all70-C4-kernels-and-full-raw33-A18-ring-map-source.json');Pi=list(map(sp,e['raw33_to_A18_columns']));mon=list(map(tuple,e['A18_basis_monomials']));mi={m:i for i,m in enumerate(mon)}
 assert all(not v or list(v.values())==[Q(1)] for v in Pi)
 copied={next(iter(v)):i for i,v in enumerate(Pi) if v};assert len(copied)==18 and {i for v in Pi for i in v}==set(range(18))
 ms='reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63'
 assert read(ms+'/report.json')['all_passed'] and read(ms+'_REVIEW/report.json')['all_passed']
 proper=read(ms+'/phase3-actual-kappa9-and-normalized-nilframe.json',['actual14_C132','kappa9_C132'])
 nuold=read(ms+'/phase1-source150-readback-fingerprints.json',['nu_A_B33'])['nu_A_B33']
 gauge=read(ms+'/phase2-canonical30-exact-values.json',['all30_constraints_exact','all30_gauge_B11_projections_exact','fixed_J30_determinant']);assert gauge['all30_constraints_exact'] and gauge['all30_gauge_B11_projections_exact'] and Q(gauge['fixed_J30_determinant'])==Q(4,9)
 def scalar(p):assert len(p)==33;return [str(Q(p[copied[i]])) for i in range(18)]
 def trace(p):assert len(p)==132;return [str(Q(p[4*copied[i]+mask])) for i in range(18) for mask in range(4)]
 actual=list(map(trace,proper['actual14_C132']));kap=list(map(trace,proper['kappa9_C132']));nu=list(map(scalar,nuold));assert all(Q(p[0])==0 for p in nu)
 oldsource=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase4-native-source-and-R12.json',['actual14','kappa9'])
 oldadapter=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase5-native-adapter-firsts.json',['N','G8'])
 oldnu=read('reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase1-source-firsts.json',['nu_first9'])['nu_first9']
 ellgate=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_INPUT_GATE_20261105_v63/actual-native-Rell4-cache-source-input-and-adapter-gate.json',['actual14_C16','actual_kappa9_C16','actual_nu_A_B_Rell4'])
 elladapter=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase2-actual-T-N-G8-C16.json',['N','actualG8'])
 # OldC36 base order=(1,z0..6,r); C16=(1,r,ell,r*ell).
 oldindices=list(range(9));rellindices=[0,mi[(7,)],mi[(8,)],mi[(7,8)]]
 def restrict(p,indices):return [p[4*i+mask] for i in indices for mask in range(4)]
 for current,old,ell in [(actual,oldsource['actual14'],ellgate['actual14_C16']),(kap,oldsource['kappa9'],ellgate['actual_kappa9_C16'])]:
  assert len(current)==len(old)==len(ell)
  for a,b,c in zip(current,old,ell,strict=True):assert list(map(Q,restrict(a,oldindices)))==list(map(Q,b)) and list(map(Q,restrict(a,rellindices)))==list(map(Q,c))
 for a,b,c in zip(nu,oldnu,ellgate['actual_nu_A_B_Rell4'],strict=True):assert list(map(Q,a[:9]))==list(map(Q,b)) and [Q(a[i]) for i in rellindices]==list(map(Q,c))
 # Actual N0/N1 are combined as coefficient functions in the SAME trace
 # basis; no second coefficient or alternative gauge is invented here.
 N36=oldadapter['N'];N16=elladapter['N'];Nfirst=[]
 for i in range(12):
  row=[]
  for j in range(12):
   a,b=list(map(Q,N36[i][j])),list(map(Q,N16[i][j]));assert len(a)==36 and len(b)==16 and a[:4]==b[:4] and a[32:36]==b[4:8]
   row.append(list(map(str,a+b[8:12])))
  Nfirst.append(row)
 Gfirst=[]
 for a,b in zip(oldadapter['G8'],elladapter['actualG8'],strict=True):
  a,b=list(map(Q,a)),list(map(Q,b));assert a[:4]==b[:4] and a[32:36]==b[4:8];Gfirst.append(list(map(str,a+b[8:12])))
 # Res4(N0 transpose) in k[tA,tB]/(tA^2,tB^2), coefficient-wise.
 L0=[{} for _ in range(48)]
 for intrinsic in range(12):
  for source in range(12):
   n=list(map(Q,Nfirst[intrinsic][source][:4]))
   for inp in range(4):
    for mask,x in enumerate(n):
     if x and not(mask&inp):L0[4*intrinsic+inp][4*source+(mask^inp)]=x
 mug=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['L0_Res4_Ntranspose_columns'])
 assert L0==list(map(sp,mug['L0_Res4_Ntranspose_columns']))
 assert all(Q(p[0])==0 for p in kap)
 K=read('reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',['source9_variable_order','selected5_pivot_order','five_actual_rational_roots'])
 prep=read('reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json',['generators3'])
 inds=[K['source9_variable_order'].index(name) for name in K['selected5_pivot_order']];assert len(inds)==5 and len(prep['generators3'])==3
 program={'core_kappa_indices4':[0,1,2,6],'D_unit':'1-kappa[1]^2','P_unit':'5+kappa[1]^2','five_pivot_kappa_indices':inds,'five_pivot_order':K['selected5_pivot_order'],'five_actual_rational_roots':K['five_actual_rational_roots'],'regular_CI_generators3':prep['generators3'],
 'generator_order8':'kappa[pivot_i]-root_i (i0..4),thenregularCI3','not_radical9':True,'native_source_type':'Koszul28_in_first8_rows_plus_free4_zero_rows;12by32 beforeRes4'}
 owned=count([e,proper,nuold,actual,kap,nu,oldsource,oldadapter,ellgate,elladapter,Nfirst,Gfirst,L0,K,prep])+8000;assert owned<150000
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-C72-cache-frame-and-nine-direction-Nfirst-source-bundle.json'
 ev={'A18_basis_monomials':mon,'trace_basis_mask_order4':['1','tauA','tauB','tauA_tauB'],'base_mask_index':'4*base_index+mask',
 'raw33_to_A18_coeff_copy_source_indices18':[copied[i] for i in range(18)],'actual14_C72':actual,'actual_kappa9_C72':kap,'actual_nu_A_B_A18':nu,
 'actual_N12_by12_point_and_first40':Nfirst,'actual_G8_point_and_first40':Gfirst,
 'all_actual14_kappa9_old8_and_Rell4_overlaps_exact':True,'all_nu_old8_and_Rell4_overlaps_exact':True,'all_N0_and_Nrfirst_C36_C16_exact':True,'Res4_N0transpose_equals_same_L0':True,
 'regular_G8_source_program':program,'G8_second_values_not_evaluated':True,
 'conditional_extension_scope':'FutureformalNextensionmustuseallactualN0/N1coefficientsaboveandSAMEactualnu.The70fullC4kernelproofthenkillsANYdifferenceinN2,m_base^2*m_base=0. Individualpointpivotsnoteligible; retainactualLeS0firstcarry.',
 'source_unit_proof':'nu0=0andkappa_scalar_point0,soC72maxidealJ=(base,tauA,tauB)hasJ^7=0; D/P scalarunits1,5. Theirinversepowerscanusebinomialseriesin(d^2)throughpower3,since(d^2)^4=0. No implicit or gauge solve needed.',
 'remaining_before_values':['sourceprogramallpolynomialtermsandfiniteC72tablewithpositiveheight/ownedbound','actualG8C72evaluationmatchesG8point/all9firstandavailableC16S2rell','weightedintrinsicS0/full70kernelsandNextensionusedonlythere','old70primitiveF2/sourcefirstcrossandall27pointfirstcarries','independentwholegateandROOTGO'],
 'no_native70_response_or_N2_or_new_source_coefficients':True,'source_owned_bound':owned}
 ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'actual_C72_source_frame_and_all9_first_bindings_verified':True,'regular_G8_program_not_evaluated':True,'native70_values_not_evaluated':True,'source_owned_bound':owned,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'owned':owned,'Nfirst_scalar_fields':12*12*40,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
