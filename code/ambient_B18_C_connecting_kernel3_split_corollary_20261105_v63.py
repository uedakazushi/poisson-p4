#!/usr/bin/env python3
"""Actual quotient-level edge split and kernel J2+k, retaining 8 unknowns."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_KERNEL3_SPLIT_20261105_v63'
PATHS={
 'value':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/phase2-full-coherent-head23-and-low-two-values.json',
 'value_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'seven':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63/honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
 'seven_peer':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63_REVIEW/report.json',
 'G':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63/genuine-four-R0-action-transport-and-cached-HB-operator-identities.json',
 'G_peer':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63_REVIEW/report.json',
 'dichotomy':'reports/AMBIENT_B18_C_CONNECTING_KERNEL3_DICHOTOMY_SOURCE_20261105_v63/actual-thirteen-parameter-edge-and-exact-kernel3-module-dichotomy.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def c(q):return{'QQ':Q(q)}if q else{}
def v(k):return{k:Q(1)}
def add(a,b):
 out=a.copy()
 for k,q in b.items():
  out[k]=out.get(k,Q(0))+q
  if not out[k]:del out[k]
 return out
def scale(a,q):return{k:x*q for k,x in a.items()if x*q}
def mul(a,b):
 if not a or not b:return{}
 assert set(a)<={'QQ'}or set(b)<={'QQ'},'Only constant times affine needed for actualsource identities'
 return scale(b,a.get('QQ',Q(0)))if set(a)<={'QQ'}else scale(a,b.get('QQ',Q(0)))
def mm(A,B):return[[total(mul(A[i][k],B[k][j])for k in range(len(B)))for j in range(len(B[0]))]for i in range(len(A))]
def total(xs):
 out={}
 for x in xs:out=add(out,x)
 return out
def encode(A):return[[[[k,str(q)]for k,q in sorted(x.items())]for x in row]for row in A]
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};d={}
 for k,p in PATHS.items():
  d[k]=json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 actual=d['value'];assert actual['simple_target_rank2_witness_nonzero']
 rg,b1=map(Q,actual['head23_canonical_J2_k2_socle_coordinates']);low=[[Q(x)for x in r]for r in actual['fivehead_low_rH_ellH_rows']]
 dd=Q(actual['actual_cyclic_top_chi'][4]);a0,a4=low[0][1],low[4][1]
 assert(dd,a0,b1,rg,a4)==(Q(-36,1681),Q(2,25),Q(1331,9216),Q(-55495,968256),Q(-4248,210125))
 assert a0*b1==Q(1331,115200)
 M=[[{}for _ in range(7)]for _ in range(4)]
 M[0][4]=c(dd);M[1][6]=c(dd);M[2][0]=c(a0);M[2][4]=c(a4)
 for j in[0,2,3,4]:M[1][j]=v('p'+str(j));M[3][j]=v('b'+str(j))
 M[1][1]=c(rg);M[3][1]=c(b1)
 Rm=[[c(Q(q))for q in r]for r in d['seven']['actual_r_action7']];Rg=[[c(Q(q))for q in r]for r in d['G']['canonical_r_Jordan_action']]
 assert mm(M,Rm)==mm(Rg,M)
 sel=[4,6,0,1];Smat=[[c(i==j)for j in sel]for i in range(7)]
 assert mm(Rm,Smat)==mm(Smat,Rg)
 target=[[M[i][j]for j in sel]for i in range(4)]
 det=dd*dd*a0*b1;assert det==Q(11979,2260608800)
 b3_over_b1=scale(v('b3'),1/b1)
 gamma_u=scale(add(v('p3'),scale(v('b3'),-rg/b1)),1/dd)
 Bs=add(v('b2'),scale(v('b3'),Q(-1,4)));Ps=add(v('p2'),scale(v('p3'),Q(-1,4)))
 Bs_over_b1=scale(Bs,1/b1);gamma_s=scale(add(Ps,scale(Bs,-rg/b1)),1/dd)
 K=[[{}for _ in range(3)]for _ in range(7)]
 K[3][0]=c(1);K[1][0]=scale(b3_over_b1,Q(-1));K[6][0]=scale(gamma_u,Q(-1));K[5][1]=c(1)
 K[2][2]=c(1);K[3][2]=c(Q(-1,4));K[1][2]=scale(Bs_over_b1,Q(-1));K[6][2]=scale(gamma_s,Q(-1))
 assert all(not x for row in mm(M,K)for x in row)
 Rk=[[c(0)for _ in range(3)]for _ in range(3)];Rk[1][0]=c(1)
 assert mm(Rm,K)==mm(K,Rk)
 # After deleting literal selector columnsS and their rows, K on source
 # rows2,3,5 is a fixed determinant1 matrix. No unknown inverse needed.
 rem=[[K[i][j]for j in range(3)]for i in[2,3,5]]
 assert rem==[[{}, {},c(1)],[c(1),{},c(Q(-1,4))],[{},c(1),{}]]
 outdata={
  'scope':'ACTUAL R0-linear splitting ONLY of honest quotientM7→genuineG4, notTold11→G, nottheKbar7extension.',
  'known_values':{'cyclic_top_d':str(dd),'h0_k1_a0':str(a0),'h4_k1_a4':str(a4),'head1_rg_p1':str(rg),'head1_k2_b1':str(b1)},
  'remaining_eight_uncomputed_actual_high_coords':['p0','b0','p2','b2','p3','b3','p4','b4'],
  'actual_edge4_by7_with_eight_formal_coords':encode(M),
  'literal_source_S_basis':['h4','r*h4','h0','h1'],'S_selector7_by4':encode(Smat),
  'all28_source_S_action_identities_verified':True,
  'S_stability_proof':'InM7, r*h4=rh4,r²h4=0,rh0=rh1=0,andallz/ell0. ThusS≅J2+k² actualR0submodule. ThesearetheliteralacceptedM7relations,notassumedactionsonTold.',
  'actual_target_frame':['delta(h4)','r*delta(h4)','delta(h0)','delta(h1)'],
  'delta_S_matrix4':encode(target),'delta_S_determinant':str(det),
  'determinant_proof':'Expansionfirstrow/secondcolumn givesd²*a0*b1 independentofall8unknownhighcoords; allthreeknownunitsaredirectactualcachedvalues.',
  'actual_quotient_split':'M7=S⊕ker(delta), S→G4isR0-lineariso. Sectionis(δ|S)^−1 withimageS. No uncomputedhighcoordinateisassignedzero.',
  'kernel_M_basis7_by3_with_four_formal_coords':encode(K),
  'kernel_basis_labels':['u=h3−(b3/b1)h1−((p3−b3*p1/b1)/d)*rh4','r*u=rh3','s=(h2−h3/4)−((b2−b3/4)/b1)h1−(((p2−p3/4)−(b2−b3/4)*p1/b1)/d)*rh4'],
  'all12_actual_delta_K_products_formally_zero':True,'all21_kernel_action_identities':True,
  'actual_kernel_module_type':'J2(r) plus k','kernel_free_minor_on_rows2_3_5':'1',
  'five_known_actual_Kbar_classes':['z5*h4','z5*r*h4','r²*h4','r³*h4','r*h3'],
  'known_five_Kbar_proof':'ITold4wasalreadyanactualfourdimensionalkernel. r*h3 mapsr*delta(h3)=0sincehead3lowtwo0andmkillsthehighgenuinecoords. ItsimageinM7isrh3!=0,soitisindependentoftheIToldfour.',
  'honest_Kbar_basis_gate':'Thefiveknownclassesplusu,s give anactualKbar7basis withfouruncomputedtrueparametersp2,b2,p3,b3. TheyarebarclassesinthehonestUCTkernel; fullformalK-liftsareNOTexpanded.',
  'Told_non_adoption':'TheliteralS_liftspan(h4,rh4,h0,h1)inToldisnotshownstable; r²h4isnonzerothere. ThisquotientsplitdoesNOTsplit0→Kbar7→Told11→G4,doesNOTsplit0→ITold4→Kbar7→(J2+k),anddoesNOTdeterminetheR0typeofKbar7.',
  'no_original_H0_torsion_exponents_length_or_nilcube_decision':True,
 }
 owned=S.count(d)+S.count(outdata)+10000;assert owned<50000;outdata['whole_source_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-M7-stable-surjection-splitting-and-kernel-three-source-basis.json';ep.write_text(json.dumps(outdata,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'ACTUAL_QUOTIENT_M7_SPLIT_G4_PLUS_J2K_WITH_EIGHT_HIGHS_FORMAL',
  'actual_kernel_M3_module':'J2 plus k','actual_M7_sequence_splits_R0':True,'actual_Kbar7_extension_split_claimed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
