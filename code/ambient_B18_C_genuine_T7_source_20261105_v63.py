#!/usr/bin/env python3
"""Actual filtered23 basis, primitiveD-free16, and genuine4 via aDquotient7."""
import hashlib,json,math,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_GENUINE_T7_SOURCE_20261105_v63'
PATHS={
 'exact':'reports/AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63/report.json',
 'exact_peer':'reports/AMBIENT_B18_C_EXACT_TORSION1123_20261105_v63_REVIEW/report.json',
 'value':'reports/AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63/all-full-tauAB12-values-and-fixed-known-rank-C4-certificates.json',
 'value_peer':'reports/AMBIENT_B18_C_FULL_TAUAB12_VALUES_20261105_v63_REVIEW/report.json',
 'formal':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/full-HB6-by2-affine-source-operator-and-old27-leaves.json',
 'formal_peer':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63_REVIEW/report.json',
 'complete':'reports/AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63/complete-full-C-cubic155-operator-with-four-formal-ell3-callbacks.json',
 'complete_peer':'reports/AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json',
}
ZERO=(0,)*4
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def const(q):return{ZERO:Q(q)}if q else{}
def var(i):return{tuple(int(i==j)for j in range(4)):Q(1)}
def p_add(a,b):
 v=a.copy()
 for k,q in b.items():
  v[k]=v.get(k,Q(0))+q
  if not v[k]:del v[k]
 return v
def p_scale(a,c):return{k:q*c for k,q in a.items()if q*c}
def p_mul(a,b):
 v={}
 for x,c in a.items():
  for y,d in b.items():
   k=tuple(a+b for a,b in zip(x,y));v[k]=v.get(k,Q(0))+c*d
 return{k:q for k,q in v.items()if q}
def matmul(A,B):
 return[[sum_poly(p_mul(A[i][k],B[k][j])for k in range(len(B)))for j in range(len(B[0]))]for i in range(len(A))]
def sum_poly(xs):
 v={}
 for x in xs:v=p_add(v,x)
 return v
def decode(x):return{tuple(k):Q(v)for k,v in x}
def encode(x):return[[list(k),str(v)]for k,v in sorted(x.items())]
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['exact']['actual_comparison_C_torsion_exponents']==[1,1,2,3]
 assert d['value']['C4_rank4_certificate_if_positive']['rows_ordered']==[0,3,2,5]
 f=[]
 for label,target,rpow in[('A',0,0),('B',1,0),('rA',0,1),('rB',1,1)]:
  for e in range(4):f.append({'label':label+'*ell^'+str(e),'raw_target':target,'r_power':rpow,'ell_power':e,'m_degree':rpow+e})
 t=[{'label':label,'raw_target':target,'r_power':r,'ell_power':e,'m_degree':r+e}
    for label,target,r,e in[('H',2,0,0),('rH',2,1,0),('ellH',2,0,1),('r²A',0,2,0),('r²B',1,2,0),('r²ellA',0,2,1),('r²ellB',1,2,1)]]
 assert dict(sorted(Counter(x['m_degree']for x in f+t).items()))=={0:3,1:6,2:6,3:6,4:2}
 assert len({(x['raw_target'],x['r_power'],x['ell_power'])for x in f+t})==23
 eF=[[int(j%4<3 and i==j+1)for j in range(16)]for i in range(16)]
 projection=[[int(j==16+i)for j in range(23)]for i in range(7)]
 formal=d['formal'];PA=formal['P3_A'];PB=formal['P3_B'];nf3=d['complete']['full_no_cut_NF6']
 a1=-Q(nf3[0][0]);b1=-Q(nf3[3][4]);assert(a1,b1)==(Q(8,25),Q(96,275))
 ops=[]
 for label,P,k,c in[('alpha_A',PA,1,a1),('alpha_B',PB,1,b1),('beta_A',PA,2,a1),('beta_B',PB,2,b1)]:
  op=p_add(p_scale(decode(P[k][0]),c),p_scale(decode(P[k][1]),Q(-1)))
  assert all(key[0]=='cache'for key in op)
  ops.append({'label':label,'operator_old_C2_only':encode(op)})
 # Here polys use four ABSTRACT variablesalphaA,alphaB,betaA,betaB.
 # No signed affine response is evaluated.
 E=[[{}for _ in range(7)]for _ in range(7)]
 E[2][0]=const(1);E[5][1]=var(0);E[6][1]=var(1);E[5][2]=var(2);E[6][2]=var(3);E[5][3]=const(1);E[6][4]=const(1)
 K=[[{}for _ in range(4)]for _ in range(7)]
 K[1][0]=const(1);K[3][0]=p_scale(var(0),Q(-1));K[4][0]=p_scale(var(1),Q(-1))
 K[2][1]=const(1);K[3][1]=p_scale(var(2),Q(-1));K[4][1]=p_scale(var(3),Q(-1));K[5][2]=const(1);K[6][3]=const(1)
 assert all(not x for row in matmul(E,K)for x in row)
 assert[[E[i][j]for j in[0,3,4]]for i in[2,5,6]]==[[const(int(i==j))for j in range(3)]for i in range(3)]
 assert[[K[i][j]for j in range(4)]for i in[1,2,5,6]]==[[const(int(i==j))for j in range(4)]for i in range(4)]
 assert all(not x for row in matmul(matmul(E,E),E)for x in row)
 known={tuple(k):Q(v)for k,v in formal['known27_old_cached_leaf_values']};profiles=[];maxH=0
 for op in ops:
  U=0;D=1
  for key,w in op['operator_old_C2_only']:
   q=Q(w);v=known[tuple(key)];U+=(abs(q.numerator)+q.denominator-1)//q.denominator*((abs(v.numerator)+v.denominator-1)//v.denominator)
   D=math.lcm(D,q.denominator*v.denominator)
  h=max((U*D).bit_length(),D.bit_length());maxH=max(maxH,h);profiles.append([op['label'],str(U),str(D),h])
 data={
  'scope':'SameactualC4/sourcebasis; D-linear primitive-free complement andgenuineell-torsion targetonly. NOTanR0/algebrasplitting.',
  'raw_targets':['A=tauA','B=tauB','H=tauA*tauB'],
  'literal_filtered23_basis_F16_then_t7':f+t,
  'basis_proof':'Fullactualgrade0..3basesare3,6,6,6(no-cutNF6). Ingrade4theactual4sourcecutminorrows0,3,2,5 isunit,soleavingliteralrows1,4 givesAr ell³,Br ell³ asactualbasis. Withm5C4zero,the23monomialsformactualfilteredbasis,notjustanuppercover.',
  'F16_E_shift_matrix':eF,
  'free_generator_labels':['A','B','rA','rB'],
  'E3_free_generator_images':['ell³A','ell³B','r ell³A','r ell³B'],
  'primitive_D_free_proof':'These4E³imagesareindependent:2literalgrade3coordinatesell³A/B,and2literalgrade4coordinatesr ell³A/B. SinceCtoriskilledbyell³,E³onC4identifiesC/Tor modellwithimE³. Hencethechosen4sourceclassesprojecttoabasisofC/Tor modell;D-NakayamagivestheirfullDspanF≅D4mappingisomorphicallytoC/Tor. ThusC=F⊕Tor_D C asDmoduleswiththisparticularF, FNOTassertedR0stable.',
  'T7':'C4/F16≅Tor_D(C) asDmodules, actualbasis(H,rH,ellH,r²A,r²B,r²ellA,r²ellB)',
  'actual_D_linear_quotient_q7_by23':projection,
  'four_ell_coefficients_UNPAIRED':ops,
  'known_old27_leaf_values_copied':formal['known27_old_cached_leaf_values'],
  'coefficient_formula':'q(ell*rH)=q(-P1³),q(ell²H)=q(-P2³). ModF,Ar³=-a1*r²ellA,Br³=-b1*r²ellB,while r ell²/ell³ belongF. Hencea1*PA30−PA21andb1*PB30−PB21. AllP4/sourceδ4 lieinm4C4⊂F. Everyz-containingcubicisquadraticmonomialtimesze_sector∈m2C4,henceinm4⊂F; Hplanecubicsalsoinm4byHquadratics. No unknownfourthcoeff isassignedzero.',
  'formal_variables':['alpha_A','alpha_B','beta_A','beta_B'],
  'symbolic_E7_rows7_columns7':[[encode(x)for x in row]for row in E],
  'symbolic_kernel7_by4':[[encode(x)for x in row]for row in K],
  'all28_EK_products_identically0':True,'fixed_E_rank3_rows':[2,5,6],'fixed_E_rank3_columns':[0,3,4],
  'fixed_kernel_freeI4_rows':[1,2,5,6],
  'E2_rank1_gate':'TheactualknownDtype1123forces(beta_A,beta_B)!=0; numericalfour-coefficientjobmustverifythisliteralnonzeroifrun. Noresponseiscomputedhere.',
  'genuine_target_isomorphism':'q inducesker(E23)/im(E23³)→ker(E7). ItskernelisF16[ell]=ell³F16=imE23³; itissurjectivebecauseC=F⊕Ctor. ThuscanonicalgenuineC[ell]isrepresentedbythefourkerE7classesrelative to thefixedsourceF, withoutnaivefiniteker orArtinell division.',
  'transport_gate':'F isonlyD-stable:r(rA)=r²A isoutsideF. DoNOTapplynaivequotientr/zactionsorclaimqR0-linear. FullR0actionsaretransportedviauniqueactualtorsionlift, stillseparateSOURCEgate. ThefourkerE7classesarequotientcoordinates, notexpandedC4/formalnil-lifts.',
  'connecting_gate':'ActualoldcomparisonRegDivell sourcecanbeprojectedbythisD-linearq; coefficient/source-orderandR0transportmustbeprovedbeforeedgevalues. No connectingvalue iscomputed.',
  'four_future_cached_value_profiles':profiles,'future_coefficient_height_bound':maxH,'future_before_operation_bound':2*maxH+1,
  'fixed_source_profile':{'max_degree':4,'numeric_envelope':40000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW':33554432},
  'no_new_physical_gauge_IFT_jet_or_scalar_responses':True,
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-filtered23-Dfree16-quotient7-and-genuine4-source.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_PRIMITIVE_DFREE16_T7_E7_GENUINE4_FOUR_CACHE_RESPONSES_UNPAIRED',
  'actual_filtered23_basis_proved':True,'actual_D_quotient7_genuine_target_proved':True,'four_scalar_response_values_computed':False,
  'R0_action_or_connecting_values_computed':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
