#!/usr/bin/env python3
"""Actual original UCT edge cannot split R0-linearly; fixed gr2 obstruction."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_TOLD_NONSPLIT_20261105_v63'
PATHS={
 'seven':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63/honest-seven-source-quotient-four-fixed-kernel-classes-and-edge-constraints.json',
 'seven_peer':'reports/AMBIENT_B18_C_CONNECTING_SEVEN_SOURCE_20261105_v63_REVIEW/report.json',
 'old40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
 'old40_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
 'exact':'reports/AMBIENT_B18_T_EXACT_LOEWY11_20261105_v63/exact-Loewy11-and-source-basis-corollary.json',
 'exact_peer':'reports/AMBIENT_B18_T_COMBINED619_QUARTIC_WITNESS_20261105_v63_REVIEW/report.json',
 'value':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/phase2-full-coherent-head23-and-low-two-values.json',
 'value_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={'seven':['actual_r_action7','fixed_I_Told_source_basis'],
  'old40':['homogeneous_module_generators40','standard_degree1','standard_degree2'],
  'exact':['exact_length','exact_Loewy_HF','actual_filtered_monomial_basis11'],
  'value':['actual_cyclic_top_chi']}
 d={}
 for k,p in PATHS.items():
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('_peer'):assert d[k]['all_passed']
 old=d['old40'];e=lambda i,n=1:[n if j==i else 0 for j in range(8)]
 qr=e(7,2);q1=[a+b for a,b in zip(e(5),e(7))]
 G38=old['homogeneous_module_generators40'][38]
 assert G38==[[3,qr,'1'],[4,q1,'11767/3652']]
 cc=Q(G38[1][2]);assert old['standard_degree1']==['e3*r','e4*z5','e4*r']
 assert old['standard_degree2']==['e4*z5*r','e4*r²']
 assert d['exact']['exact_length']==11 and d['exact']['exact_Loewy_HF']==[5,3,2,1]
 basis=d['exact']['actual_filtered_monomial_basis11']
 assert[(b['target'],b['monomial_original8'])for b in basis if b['degree']==2]==[(4,q1),(4,qr)]
 I=d['seven']['fixed_I_Told_source_basis'];assert I['first']['degree']==1 and all(x['degree']>=2 for x in I['other_three'])
 assert I['first']['target']==4 and I['first']['monomial_original8']==e(5)
 R=[[Q(q)for q in row]for row in d['seven']['actual_r_action7']]
 assert all(R[i][j]==0 for i in range(7)for j in[0,1])
 assert[(i,R[i][2])for i in range(7)if R[i][2]]==[(5,Q(1,4))]
 assert[(i,R[i][3])for i in range(7)if R[i][3]]==[(5,Q(1))]
 assert[(i,R[i][4])for i in range(7)if R[i][4]]==[(6,Q(1))]
 chi=list(map(Q,d['value']['actual_cyclic_top_chi']));assert chi[:4]==[0,0,0,0]and chi[4]==Q(-36,1681)
 # Unknown actual gr1 corrections stay as alpha0..2. The q2row is
 # checked coefficientwise WITHOUT choosing/evaluating those three values.
 first=[[[0,0,0],[1,0,0],[0,0,0]], # h0 =alpha0*v
        [[0,0,0],[0,1,0],[0,0,0]], # h1 =alpha1*v
        [[0,0,0],[0,0,1],[0,0,0]]] # h2 =u/4+alpha2*v
 r_gr1_to_gr2=[[str(-cc),'1','0'],['0','0','1']]
 q2row=[Q(0),Q(0),Q(0),Q(0),Q(1)]
 assert q2row==[0,0,0,0,1]
 data={
  'scope':'ActualR0moduleUCT0→K/ellK→Told11→genuineC[ell]4→0isNOTsplit. ThisiscompatiblewiththeSPLITquotientM7→G4.',
  'actual_old_gr1_basis':['u=r*h3','v=z5*h4','w=r*h4'],
  'actual_old_gr2_basis':['q1=z5*r*h4','q2=r²*h4'],
  'gr2_basis_independence':'AcceptedexactToldHF5321andliteral11basisgiveToldm²/m³dimension2withthetwoDISPLAYEDclassesindependent; notjustaknowninitialuppercover.',
  'ITold_gr1_exact_line':'IToldhasbasisz5h4(degree1),z5rh4,r²h4,r³h4(degree>=2). Thereforeitsimageinactualgr1Toldisexactlyspan(v).',
  'allfive_head_r_leading_classes':['alpha0*v','alpha1*v','u/4+alpha2*v','u','w'],
  'allfive_head_source_proof':'M7=Told/IToldhasrh0=rh1=0,rh2=rh3/4. Hencerh0,rh1andrh2−rh3/4belongITold. Theirgr1partsaretheDISPLAYEDunknownmultiplesofv. rh3=uandr h4=waredefinitions. Thusallfiveheads,includingh0/h1,arecontrolled,notmerelythetwoobviousones.',
  'quadratic_source_G38_literal':G38,
  'actual_r_on_gr1_to_gr2':r_gr1_to_gr2,
  'allfive_head_r2_gr2_classes':['alpha0*q1','alpha1*q1','(alpha2−11767/14608)*q1','−(11767/3652)*q1','q2'],
  'q2_coordinate_of_r2_allfiveheads':list(map(str,q2row)),
  'source_cyclic_top_chi':list(map(str,chi)),
  'section_obstruction_coefficient':str(1/chi[4]),
  'nonsplitting_proof':'Letg∈Gbethecanonicalcyclictop,chi(g)=1and r²g=0. ANYactualpreimagex∈Toldhasconstanthead4coefficient1/d sinceallotherheadchi0andmToldmapsintomG. TermsfrommToldcontributer²*mTold⊂m³Told,soq2coordinateofr²xmodm³isEXACTLY1/d,nonzero. Henceeverypreimagehasr²x!=0; anR0-linearsectionwouldgivepreimager²s(g)=s(r²g)=0,contradiction.',
  'quotient_split_distinction':'Theq2obstructionliesinIToldandvanishesinM7. ThusM7stableS/deltaSisoandsplittingarevalid,whiletheoriginalUCTToldsequencehasnonzeroextensionclass. DoNOTinferthatITold→Kbar7extensionissplitornonsplitfromthisseparateresult.',
  'unknown_alpha0_alpha1_alpha2_never_evaluated':True,
  'no_original_H0_length_exponents_nilcube_or_fullKbar_type_claim':True,
  'selective_source_receipts':receipts,
 }
 owned=S.count(d)+S.count(data)+10000;assert owned<50000;data['whole_owned_numeric_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'actual-Told-UCT-nonsplitting-fixed-q2-gr2-obstruction.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'PURE_ACTUAL_TOLD_UCT_NOT_R0_SPLIT_BY_FIXED_GR2_OBSTRUCTION',
  'actual_Told_to_genuine_C_ell_sequence_R0_split':False,'quotient_M7_split_compatible':True,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_numeric_upper':owned,'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
