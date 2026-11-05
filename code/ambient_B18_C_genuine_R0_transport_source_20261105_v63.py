#!/usr/bin/env python3
"""Genuine C[ell] R0 transport: free corrections, cachedHB identities, no new dots."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63'
PATHS={
 'T7':'reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63/actual-filtered23-Dfree16-quotient7-and-genuine4-source.json',
 'T7_peer':'reports/AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63_REVIEW/report.json',
 'exact':'reports/AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63/report.json',
 'exact_peer':'reports/AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63_REVIEW/report.json',
 'formal':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/full-HB6-by2-affine-source-operator-and-old27-leaves.json',
 'formal_peer':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63_REVIEW/report.json',
 'HB':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63/phase2-full-six-coordinate-HB-matrix-and-fixed-rank-certificate.json',
 'HB_peer':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_VALUES_20261105_v63_REVIEW/report.json',
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ell_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decode(x):return{tuple(k):Q(v)for k,v in x}
def encode(x):return[[list(k),str(v)]for k,v in sorted(x.items())]
def add(a,b):
 v=a.copy()
 for k,q in b.items():
  v[k]=v.get(k,Q(0))+q
  if not v[k]:del v[k]
 return v
def scale(a,c):return{k:q*c for k,q in a.items()if q*c}
def mm(A,B):return[[sum((A[i][k]*B[k][j]for k in range(len(B))),Q(0))for j in range(len(B[0]))]for i in range(len(A))]
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['exact']['actual_comparison_C_torsion_exponents']==[1,1,2,3]
 f=d['formal'];P=[f['P3_A'],f['P3_B']];nf=f['full_NF4'];ratios=[]
 for side in range(2):
  a1=-Q(nf[3*side][5*side+1]);a2=a1*a1-Q(nf[3*side][5*side]);a3=-Q(nf[3*side+2][5*side+1]);ratios.append([a1,a2,a3])
 # Derived ACTION operators, never paired with the stored scalar leaves.
 action_ops=[];Y2=[]
 for side in range(2):
  a1,a2,a3=ratios[side]
  alpha=[add(scale(decode(P[side][k][0]),a1),scale(decode(P[side][k][1]),Q(-1)))for k in range(3)]
  for j in[1,2]:
   c=add(scale(decode(P[side][j][0]),a2),scale(decode(P[side][j][2]),Q(-1)))
   dd=add(scale(decode(P[side][j][0]),a3),scale(decode(P[side][j][3]),Q(-1)))
   response=add(add(alpha[j-1],scale(alpha[j],a1)),scale(c,Q(-1)))
   old=decode(f['full_HB_6_by2_affine_operators'][3*side][j-1]);assert response==scale(old,Q(-1))
   assert all(k[0]!='mu03'for k in response)
   action_ops.append({'side':side,'low_generator':j-1,'operator':encode(response),'identically_minus_cached_HB_operator':True})
   Y2.append({'side':side,'low_generator':j-1,'r_ell_free_coefficient_c':encode(c),'ell2_free_coefficient_d':encode(dd)})
 basis=d['raw33']['actual33_basis'];p2=[basis.index([0,2,0]),basis.index([0,1,1])]
 phys={r['physical_row']:r for r in d['C2']['physical12_coefficients']};ells={r['physical_row']:r for r in d['ell']['all12_actual_physical_coefficients_pure_ell3']}
 flags=[]
 for i,record in enumerate(d['frame']['all12_verified_source_records'][:7]):
  row=record['physical_row'];assert record['delta_first9_zero']and record['delta_and_bA_pureplane_second_zero']
  assert list(map(Q,record['bA_first_z7']))==[Q(i==j)for j in range(7)]
  for key,ellkey in[('delta','delta_jet3'),('bA','bA_jet3')]:
   vals=dict(phys[row][key]);assert all(Q(vals.get(j,'0'))==0 for j in p2)and Q(ells[row][ellkey][2])==0
  flags.append([i,row,'delta2_plane0','bAi2_plane0','bAi1=zi'])
 oldHB=[[Q(v)for v in row]for row in d['HB']['actual_full_HB6_by2_matrix']]
 M=[[-q for q in oldHB[i]]for i in[0,3]];assert M==[[Q(3,5),Q(-6)],[Q(121,80),Q(-121,8)]]
 R=[[Q(0)]*4 for _ in range(4)]
 for i in range(2):R[2+i][:2]=M[i]
 assert mm(R,R)==[[0]*4 for _ in range(4)]and all(M[i][1]==-10*M[i][0]for i in range(2))and M[0][0]
 C=[[Q(1),0,Q(10),0],[0,0,Q(1),0],[0,Q(3,5),0,0],[0,Q(121,80),0,Q(1)]]
 I=[[Q(1),Q(-10),0,0],[0,0,Q(5,3),0],[0,Q(1),0,0],[0,0,Q(-121,48),Q(1)]]
 ident=[[Q(i==j)for j in range(4)]for i in range(4)];assert mm(C,I)==ident==mm(I,C)
 J=[[Q(0)]*4 for _ in range(4)];J[1][0]=1;assert mm(mm(I,R),C)==J
 data={
  'scope':'ACTUAL genuineG=C[ell]=kerE_C4/imE_C4³, transportedusingtheprimitiveD-freeF. NotnaiveR0actionontheDquotientC4/F.',
  'source_genuine_basis':['g_rH=q^-1(t1-alpha_A*t3-alpha_B*t4)','g_ellH=q^-1(t2-beta_A*t3-beta_B*t4)','g_Ahigh=q^-1(t5)','g_Bhigh=q^-1(t6)'],
  'full_quartic_free_top_source':'m4C4 isspanAr ell³,Br ell³ andiscontainedinimE³; allm5C4=0.',
  'free_degree2_correction_operators_UNPAIRED':Y2,
  'actual_representative_construction':'Forlowg_jchoosek_j as SOURCEkerE7coordinate. Modm4,ell*k_j hasonlyfreecubic terms c_A*r ell²A+d_A*ell³A+ Banalog; subtract Y2_j=c_A*r ellA+d_A*ell²A+B. Remainderell*(k_j-Y2_j) liesinm4C4, sochooseY3_j∈span(r ell²A,r ell²B)⊂Fdegree3 withellY3_j equaltheremainder. Thenk_j-Y2_j-Y3_j∈kerE23. TheunknownY3_j isRETAINED,neverzero-filled. Highg representatives t5/t6 similarlycorrectbyFdegree3.',
  'unknown_Y3_action_invariance':'ForANYY∈Fdegree3,mY⊂m4C4⊂imE³, soallR0maximalidealactionsonGareindependentoftheseunknowncorrections. Thisiswhyfullquartic/freeoffdiagonalvaluesarenotneeded.',
  'literal_zH_source_checks':flags,
  'z_action_proof':'RHi=(0,delta_i,bAi),bAi1=zianddelta_i1=0. Theverifiedplane2partsdelta_i2andbAi2are0,soalltheirquadratictermscontainz. Sinceeachz_j*rawtarget∈m²C4 bythefirst21frame, theequationziH=-delta_iB-(bAi-zi)H impliesziH∈m³C4. Hencezi*(rHorellH)∈m4; zi*anydegree2correctionalso∈m4. AllvanishinGmodimE³. No mixedsourceleafisassignedzero.',
  'high_socle_proof':'Bothhighgenuineclasseshaverepresentativesofdegree3; multiplyingbym landsinm4C4⊂imE³. Thusallm killthem.',
  'r_transport_proof':'Usingthecorrectedlowrepresentatives, q(rg_rH)=alpha0+a1*alpha1-c1 onAandb1onB, q(rg_ellH)=alpha1+a1*alpha2-c2. Free d*ell² termsandunknownY3 contributeonlyF/imE³anddrop. Substitutingalpha_k=a1*Pk30-Pk21,c_k=a2*Pk30-Pk12 givesEXACTLYminusNF4(ellP0-rP1) andminusNF4(ellP1-rP2), checkedcoefficientwiseagainstthefrozenaffineoperators.',
  'four_r_action_operator_identities':action_ops,
  'actual_r_action_low2_to_high2':[[str(x)for x in row]for row in M],
  'actual_r_matrix4':[[str(x)for x in row]for row in R],
  'actual_z0_through_z6_and_ell_actions':'all4x4zero; inoriginalr,scoordinates,sactionequalsractionbecauseell=5(s-r)',
  'fixed_r_rank1_witness':'3/5, secondlowcolumn=-10first','r_squared_zero':True,
  'actual_R0_module_type':'R0/(z0,...,z6,ell,r²) ⊕ (R0/m)^2',
  'actual_mu_R0_genuine_C_ell':3,'genuine_C_ell_m2_annihilated':True,'genuine_C_ell_socle_dimension':3,
  'canonical_J2_plus_two_k_basis':['g_rH','r*g_rH','10*g_rH+g_ellH','g_Bhigh'],
  'basis_change4':[[str(x)for x in row]for row in C],'basis_change_inverse4':[[str(x)for x in row]for row in I],
  'both_basis_change_products_identity':True,'canonical_r_Jordan_action':[[str(x)for x in row]for row in J],
  'no_new_r_action_scalar_dots':'ThefouractionsareexactnegativesofALREADYpeeredHBnumbers; noaffineleafvaluepairing/newcoefficientqueryisrunhere.',
  'retained_limits':['Numericalalpha/betaE7coordinatesstillbelongtotheirfourdotgate; fullC4representativesY3arenotcomputed.',
   'ThisdoesnotmakeF R0stableanddoesnotgiveall23C4actiontables.',
   'Originalcomparisonconnectingvalues/Kbar7/originalH0andnilcuberemainuncomputed.'],
  'source_profile':{'degree':4,'numeric_envelope':50000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'genuine-four-R0-action-transport-and-cached-HB-operator-identities.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_ACTUAL_GENUINE_C_ELL_R0_TYPE_J2_PLUS_K2_NO_NEW_ACTION_DOTS',
  'actual_genuine_C_ell_R0_mu':3,'actual_genuine_C_ell_m2_zero':True,'new_scalar_action_dots_evaluated':False,
  'original_connecting_or_H0_values_computed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
