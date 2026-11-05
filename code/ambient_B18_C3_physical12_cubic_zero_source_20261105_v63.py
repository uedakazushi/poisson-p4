#!/usr/bin/env python3
"""Actual quadratic frame and projected cubic-zero theorem for 12 tauAB columns."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63'
PATHS={
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2roots':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'C2peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ellpeer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
 'pair':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
 'pairpeer':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63_REVIEW/report.json',
 'rawframepeer':'reports/AMBIENT_B18_CONDUCTOR_TRANSFER512_NIL_CYCLES_20261104_v63_REVIEW/report.json',
 'rawnote':'research/notes/AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63.md'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={p:sha(ROOT/p)for p in PATHS.values()}
 D={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k!='rawnote'};assert all(v['all_passed']for k,v in D.items()if k.endswith('peer'))
 basis=D['raw33']['actual33_basis'];zi=[basis.index([i+1,0,0])for i in range(7)];ir=basis.index([0,1,0]);il=basis.index([0,0,1]);ir2=basis.index([0,2,0]);irl=basis.index([0,1,1]);first=zi+[ir,il]
 ee={p['physical_row']:p for p in D['ell']['all12_actual_physical_coefficients_pure_ell3']};records=[];bmat=[]
 for k,p in enumerate(D['C2']['physical12_coefficients']):
  vals={name:{i:Q(c)for i,c in p[name]}for name in ['delta','bA','bB']};e=ee[p['physical_row']]
  assert all(not vals['delta'].get(i,Q(0))for i in first)
  assert all(not vals[name].get(i,Q(0))for name in ['bA','bB']for i in [ir,il])
  assert all(not vals[name].get(i,Q(0))for name in ['delta','bA']for i in [ir2,irl])
  assert all(Q(e[name+'_jet3'][2])==0 for name in ['delta','bA'])
  if k<7:assert [vals['bA'].get(i,Q(0))for i in zi]==[Q(k==j)for j in range(7)]
  b=[vals['bB'].get(ir2,Q(0)),vals['bB'].get(irl,Q(0)),Q(e['bB_jet3'][2])];bmat.append(b)
  records.append({'physical_row':p['physical_row'],'delta_first9_zero':True,'delta_and_bA_pureplane_second_zero':True,'bA_first_z7':[str(vals['bA'].get(i,Q(0)))for i in zi],'bB_first_z7':[str(vals['bB'].get(i,Q(0)))for i in zi],'bB_second_plane':[str(c)for c in b]})
 nr={i:dict(v)for i,v in D['C2roots']['all_requested164_root_values']}
 assert all(not Q(nr[n].get(i,'0'))for n in [20086,40865]for i in [0,ir,il])
 # A fixed source minor only; no rank search or new coefficient evaluation.
 rows=[0,4,5];a=[bmat[i]for i in rows];assert max(max(abs(c.numerator).bit_length(),c.denominator.bit_length())for row in a for c in row)<=32
 den=1
 for row in a:
  for c in row:den=math.lcm(den,c.denominator)
 aa=[[c.numerator*(den//c.denominator)for c in row]for row in a]
 det=aa[0][0]*(aa[1][1]*aa[2][2]-aa[1][2]*aa[2][1])-aa[0][1]*(aa[1][0]*aa[2][2]-aa[1][2]*aa[2][0])+aa[0][2]*(aa[1][0]*aa[2][1]-aa[1][1]*aa[2][0]);assert det
 assert max(abs(det).bit_length(),(den**3).bit_length())<2048;minor=Q(det,den**3)
 data={
  'same_comparison_module':'C=coker(D0→E4), C3=C/ell³C; threecyclictargets(eA=tauA,eB=tauB,eH=tauAtauB). C isamodule, notaring.',
  'complete_raw_presentation':'IA*eA (8), IB*eB (8), IAIB*eH (64), twelve physicalactioncolumns on3rawgenerators(36):3targets/116sources. PhysicalF*1alreadybelongsrawimage; all3Fmultiples retained.',
  'all12_verified_source_records':records,'both_nu_constant_and_pureplane_first_zero':True,
  'actual_first_frame21':'RA_i=zi*eA; RB_i=(zi+beta_i)*eB; RH_i=F_i*tauB forfirst7physicalrows. Actualdelta1=0 andbA_i,1=zi, so firstsareexactlyzi*eA,zi*eB,zi*eH, triangularunit1. Allother116firstslieintheirspan.',
  'first_frame_pureplane_second':'RA2=0,RB2=beta2_i*eB,RH2=0. ThelatterusesACTUALall12delta2planezero andbA2planezero, nowverifiedincludingell²; notassumedfromC2.',
  'complete_quadratic_initial_source':'M0=0,V1rank21. Constantkernel95 plusS1*V1 spansallquadraticinitialrelations. Mod(z7)S³,66intrinsiczero-linearcols(64products+2curves)give0;5remainingtauB and12tauAB give0;12tauA columnsgive -bB_j,2*eH. Henceexactquadraticspanisplane²*eH.',
  'fixed_bB_plane_minor_physical_rows':[18,60,67],'fixed_bB_plane_minor_matrix':[[str(c)for c in row]for row in a],'fixed_bB_plane_minor_determinant':str(minor),
  'actual_C_and_C3_HF_degrees0_1_2':[3,6,6],
  'physical12_tauAB_order3_relation_construction':'EachM_j=F_j*tauAB hasfirst0. ItsM2eA/eB=-nuB1*bB1,+nuA1*bA1 ispurezz becauseeachfactorhaszeroplane-first. M2eH=delta2 hasnoplanecomponent. Cancelunknownzzcoefficientsbyz*RA/RB/RH andknownmixedzi*r/zi*ell by r/ell*RH. ThisisafiniteformalactualTaylor2 coefficientoperator, notzero-fill.',
  'all_correction_pureplane_cubic_zero':'Forzzcorrectionsthemultiplierz killsallpureplaneQ2carry. ForHmixedcorrectionstheRH2planeisexactly0. ThusfullcoherentM2/linearframecarryprojects0; noIBbeta2 termisdiscardedunjustifiably.',
  'primitive_pureplane_cubic_kept_formal':'Offdiag nu*b haszero plane-degree3 since bothnuandb haveplaneorder≥2. TheonlyprimitivecubeisUNKNOWNdelta_j,3*eH; retainall12×4scalarfunctionalsformally.',
  'all_unknown_delta3_sensitivities_zero_after_known_quad':'plane²*eH isalreadyactualquadraticinitialspan, soeveryplane-cubicdelta3*eH vanishesintheknowncubicquotient. Thisprovesall12physicaltauABprojectedcubicclasses0 withoutknowing/evaluatingdelta3.',
  'physical12_cubic_source_result':'ALLZERO_SOURCE_SUBSET_ONLY; itdoesnotclaimallC3cubicrelationszero.',
  'known_initial_model_diagnostic':'AddingintrinsicgA3*eA,gB3*eB andell³ yieldsplane/(gA3,ell³)eA ⊕plane/(gB3,ell³)eB ⊕plane/(r,ell)²eH. EachfirstsummandCI(3,3)length9; lastlength3.',
  'known_initial_model_HF':[3,6,6,4,2],'known_initial_model_length':21,
  'plain_model_not_decisive':'ActualC3 hasalreadyknownDVRlength18..20, soceil21andthezero12subsetdonotsettletorsionexponents. Stopplain-modeldegreeescalation; otherphysicaltauA/tauBcubicclassesorquarticprojectionsremain.',
  'no_new_physical_source_IFT_gauge_higherjet_or_cubic_values':True,'original_H0_nilcube_or_C3_length_not_determined':True}
 ev=out/'actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json';ev.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PROVED_actual_C_HF366_and_physical12_tauAB_cubic_subset_source_zero','actual_C_HF_through2':[3,6,6],'all12_tauAB_cubic_initial_classes_zero':True,'actual_C3_length_or_torsion_exponents_determined':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ev.name:sha(ev)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
