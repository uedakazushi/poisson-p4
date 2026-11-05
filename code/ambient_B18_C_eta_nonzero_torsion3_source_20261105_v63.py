#!/usr/bin/env python3
"""Conditional pure source: one eta!=0 gives full HF366544... and type1123."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_ETA_NONZERO_TORSION3_SOURCE_20261105_v63'
PATHS={
 'line':'reports/AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63/complete-full-C-cubic155-operator-with-four-formal-ell3-callbacks.json',
 'line_peer':'reports/AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json',
 'reduction':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json',
 'reduction_peer':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/report.json',
 'obstruction':'reports/AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63/full-C-cubic-obstruction-and-five-graph03-source-zeros.json',
 'obstruction_peer':'reports/AMBIENT_B18_C_FULL_GRAPH03_FORCED_ZERO_20261105_v63_REVIEW/report.json',
 'cut':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63/one-CUT1-residual-and-formal-coordinate-rank-certificate.json',
 'cut_peer':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Q(0))
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 line=d['line'];assert line['complete_formal_operator_shape']==[6,155]
 nf3=[[Q(x)for x in row]for row in line['full_no_cut_NF6']]
 nf4=[[Q(0)]*15 for _ in range(6)];curves=[]
 for s in range(2):
  # normalized g=r³+a1r²ell+a2r ell²+a3ell³ from the reviewed fullNF3.
  a=[Q(1)]+[-nf3[3*s+k][4*s]for k in range(3)];curves.append(a)
  b=[[a[1]**2-a[2],-a[1],1,0,0],
     [a[1]*a[2]-a[3],-a[2],0,1,0],
     [a[1]*a[3],-a[3],0,0,1]]
  for k in range(3):
   assert dot(b[k],a+[Q(0)])==dot(b[k],[Q(0)]+a)==0
   nf4[3*s+k][5*s:5*s+5]=list(map(Q,b[k]))
 # The actual initial relation ell³(A-B/4) follows conditionally from
 # eta!=0 and exact kernel(NF3). Its r/ell multiples are actual initials.
 rpoly=[Q(0)]*15;epoly=rpoly.copy();rpoly[3]=1;rpoly[8]=Q(-1,4);epoly[4]=1;epoly[9]=Q(-1,4)
 cols=[[dot(row,v)for row in nf4]for v in[rpoly,epoly]]
 assert cols==[[0,1,0,0,Q(-1,4),0],[0,0,1,0,0,Q(-1,4)]]
 assert cols[0][1]*cols[1][2]-cols[0][2]*cols[1][1]==1
 assert d['cut']['actual_C3_HF_if_positive']==[3,6,6,4,0]
 assert d['cut']['cut3_torsion_pattern_if_positive']==[1,1,2,'e>=3']
 assert d['obstruction']['actual_nonzero_r_times_ell2_classes']['A']==['0','1','0','0']
 data={
  'scope':'CONDITIONAL pure theorem for FULL actual comparison C. Assumption: at least one of the FOUR free eta_j=delta_ell³ is nonzero. No eta value is evaluated here.',
  'same_source_and_initial_membership':'The reviewed complete full155 source gives actual cubic initial class eta*ell³(A-B/4). Under eta!=0, invert that fixed nonzero QQ constant and use exact known NF3 kernel to obtain the actual homogeneous initial relation ell³(A-B/4). Its r/ell multiples are actual quartic initials. This is NOT a claim that arbitrary upper-cover relations are actual.',
  'normalized_full_curve_cubics':[[str(x)for x in a]for a in curves],
  'full_known_quartic_NF6_rows15_columns':[[str(x)for x in row]for row in nf4],
  'two_actual_conditional_quartic_columns':[[str(x)for x in col]for col in cols],
  'fixed_rank2_rows':[1,2],'fixed_rank2_minor':'1',
  'conditional_full_HF_proof':'The no-cut quartic upper cover hasdim6, and these two actual initials haveindependentimages, soHF4<=4. Reductionm^(n+1)=ell*m^n plusrank_D C4 givesmu_D(m4C)=HF4>=4. HenceHF4=4 andm4C isD-free4. Multiplicationbyell givesallhigherHF4; HF3=5 by the complete eta theorem.',
  'conditional_full_C_HF':[3,6,6,5,'4 for every degree>=4'],
  'conditional_Hilbert_series':'(3+3t-t^3-t^4)/(1-t)',
  'conditional_exponent_atmost4':'ell4C⊂m4C isD-torsionfree, soell4Tor_D(C)=0. The accepted type112e,e>=3 now leavesonlye3ore4.',
  'F_torsion_structure':'F=m3C hasrank4/mu5, andellF=m4C isD-free. ThusTor_D(F)=k. Moreoverm*Tor_D(F)⊂mF∩Tor_D(C)=m4C∩Tor_D(C)=0.',
  'exclude_e4_steps':[
   'Assumee=4. The accepted type(1,1,2,4) givesell³Tor_D(C)≅k. ItliesinF andtherefore equalsTor_D(F).',
   'C/m4C haslength3+6+6+5=20, whileC/ell³C haslength19andm4C⊂ell³C. Henceell³C/m4C haslength1.',
   'Tor_D(F)=k liesinell³C andintersectsm4C=0. Thereforeell³C=m4C⊕Tor_D(F), asD-modules only.',
   'Foranyx writeell³x=ell*y+q withy∈F,q∈Tor_D(F). Put t=ell²x−y. Thenell*t=q.',
   'r*t∈F andell*(r*t)=r*q=0, since m*Tor_D(F)=0. Thusr*t∈Tor_D(F)⊂ell³C.',
   'Consequentlyr ell²x=r*y+r*t∈m4C+ell³C=ell³C. Thiscontradictsthe actualnonzero r ell²A inC3.'
  ],
  'conditional_exact_comparison_C_DVR_type':{'free_rank':4,'torsion_exponents':[1,1,2,3],'torsion_length':7},
  'conditional_comparison_C_torsion_m_annihilator':'m4*Tor_D(C)=0 sinceitliesinD-free m4C. No m3-annihilation isclaimed.',
  'conditional_next_target':'C/ell4C haslength23, andactualC[ell] iscanonicallykerE/imE³ onthatquotient. Thisdoesnotconstruct23actions/basis/edge.',
  'no_original_H0_nilcube_or_total_exhaustion_transfer':True,
  'all_eta_values_UNREAD':True,'allzero_eta_scope':'The conditional theorem doesnotapply ifall4eta=0; doNOTinferexponent3oranyfinitehigherquotientlengthfromthisalone.',
  'source_profile':{'degree':4,'QQ_matrix_slots':90,'numeric_envelope':25000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'conditional-eta-quartic-rank2-and-exclusion-of-exponent4.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_CONDITIONAL_NONZERO_ETA_IMPLIES_EXACT_C_TORSION1123_NO_VALUES',
  'actual_eta_values_evaluated':False,'unconditional_remaining_exponent_determined':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
