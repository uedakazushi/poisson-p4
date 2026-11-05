#!/usr/bin/env python3
"""Pure/source thirteen-parameter edge and exact length3-kernel dichotomy."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_KERNEL3_DICHOTOMY_SOURCE_20261105_v63'
PATHS={
 'seven':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63/honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
 'seven_peer':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63_REVIEW/report.json',
 'G':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63/genuine-four-R0-action-transport-and-cached-HB-operator-identities.json',
 'G_peer':'reports/AMBIENT_B18_C_GENUINE_R0_TRANSPORT_SOURCE_20261105_v63_REVIEW/report.json',
 'adjoint':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63/actual-three-row-sigma2-source-adjoint-and-coherent-fourth-frontier.json',
 'adjoint_peer':'reports/AMBIENT_B18_C_CONNECTING_SIGMA2_ADJOINT_SOURCE_20261105_v63_REVIEW/report.json',
 'plane':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63/all56-plane-orders-surviving-sigma2-adjoint-and-U9-namespace-correction.json',
 'plane_peer':'reports/AMBIENT_B18_C_CONNECTING_PLANE_ORDER_SOURCE_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
 v=a.copy()
 for k,q in b.items():
  v[k]=v.get(k,Q(0))+q
  if not v[k]:del v[k]
 return v
def sc(a,c):return{k:q*c for k,q in a.items()if q*c}
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'seven':['actual_r_action7','source7_basis','fixed_I_Told_source_basis'],
  'G':['canonical_r_Jordan_action'],
  'adjoint':['actual_sigma0_56_by5_sparse'],
  'plane':['actual_plane2_H_three_monomials_rows3_columns56']}
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 for j,row in[(1,23),(2,29),(3,30)]:assert d['adjoint']['actual_sigma0_56_by5_sparse'][j]==[[row,'1']]
 H2=[[Q(q)for q in row]for row in d['plane']['actual_plane2_H_three_monomials_rows3_columns56']]
 assert all(H2[k][col]==0 for k in[1,2]for col in[23,29,30])
 Rm=[[Q(q)for q in row]for row in d['seven']['actual_r_action7']];Rg=[[Q(q)for q in row]for row in d['G']['canonical_r_Jordan_action']]
 labels=['d']+[f'p{i}'for i in range(5)]+['a0','a4']+[f'b{i}'for i in range(5)];assert len(labels)==13
 M=[[{}for _ in range(7)]for _ in range(4)]
 def var(k):assert k in labels;return{k:Q(1)}
 M[0][4]=var('d');M[1][6]=var('d')
 for j in range(5):M[1][j]=var(f'p{j}');M[3][j]=var(f'b{j}')
 M[2][0]=var('a0');M[2][4]=var('a4')
 for i in range(4):
  for j in range(7):
   lhs={};rhs={}
   for k in range(7):lhs=add(lhs,sc(M[i][k],Rm[k][j]))
   for k in range(4):rhs=add(rhs,sc(M[k][j],Rg[i][k]))
   assert lhs==rhs
 simple=[[var('a0'),{},{}],[var('b0'),var('b1'),add(var('b2'),sc(var('b3'),Q(-1,4)))]]
 data={
  'scope':'Actual SOURCE constraints on remainingM7→genuineG4 and universal module-kernel classification; no edge coefficients or originalH0 computed.',
  'actual_source_module':'J2²+k³, source7=(h0,h1,h2,h3,h4,rh3,rh4)',
  'actual_target_module':'J2+k², canonical4=(g,rg,k1,k2), allz/ell0,r²0',
  'low_source_proof':'V1_plane=0, so lowgenuinecoordinates areH2_(rell/ell²)*sigma0. Forsourceheads1,2,3 sigma0=e23,e29,e30 andBOTH H2rows inthosecolumnsareliteral0. Thusbothlowcoordinatesonallthreeheadsare0,withoutan affineleafpairing.',
  'actual_cyclic_top':'chi=(0,0,0,0,d), d!=0bysurjectivity; oldchi parameterc=0SOURCEproved',
  'actual_simple_k1_top':'a=(a0,0,0,0,a4), a0!=0bysurjectivityofthetargettoprows; no numericalvaluechosen',
  'thirteen_unknown_parameters':labels,
  'actual_edge_matrix4_by7_with_formal_parameters':[[[[k,str(q)]for k,q in sorted(p.items())]for p in row]for row in M],
  'all28_r_equivariance_entries_identically_zero':True,
  'source_three_simple_basis':['h0','h1','h2-h3/4'],
  'simple_target_projection':'ann(r:G)/(rG)=span(k1,k2), nottheentireG/rGcyclictop',
  'simple_source_matrix2_by3':[[[[k,str(q)]for k,q in sorted(p.items())]for p in row]for row in simple],
  'unconditional_simple_source_rank_minimum':1,
  'exact_kernel_dichotomy':'rank(simple2×3)=2 iff (b1,b2-b3/4)!=(0,0) iffker(M7→G4)≅J2+k. rank1 iffboth0 iffker≅k³. Surjectivityrulesoutrank0.',
  'proof':'AcyclicJ2targetisprojectiveoverA=k[r]/r²,soits surjection splits; equivalentlytop2fromthetwoJ2sourcechains hasonecyclicconstraint. Thetargettwo simpletopsreceiveA=imagesofthethree simplesourcegenerators. If rankA2, thetwo simpletargets split fromsourcek³; remainingJ2²+k mapsontoJ2 andkernelJ2+k. IfrankA1, surjectivityforces theone survivingkernel direction of thecyclictopconstraint to map nontrivially to simpletarget/imA. HenceeveryactualkernelvectormusthavezeroJ2top,so r kills thelength3kernel, givingk³. No diagonalization/full17valuecomputationneeded.',
  'cheap_selected_witness':'simpleminor(h0,h1)=a0*b1, independentoftheunknownb0. Sincea0!=0,SAMEhead23canonicalk2!=0alonecertifiestheJ2+kcase. Ifb1=0,keepb2-b3/4unknown; doNOTinferk³.',
  'kernel_scope':'Thislength3moduleisquotientofactualK/ellKbyITold4. Thehonestexact0→ITold4→Kbar7→kernel3→0isretained; notKbar/I Kbar. TheextensionandoriginalH0torsionexponentsarenotdeterminedbythisdichotomy.',
  'values_not_evaluated':True,'selective_source_receipts':receipts,
 }
 owned=S.count(d)+S.count(data)+12000;assert owned<50000;data['whole_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-thirteen-parameter-edge-and-exact-kernel3-module-dichotomy.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PURE_SOURCE_THIRTEEN_PARAMETERS_KERNEL_J2K_OR_K3_NO_EDGE_VALUES',
  'remaining_actual_edge_parameter_count':13,'low_head1_2_3_structural_zero':True,'actual_kernel_case_selected':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
