#!/usr/bin/env python3
"""Pure C4 length/HF and a full-HB residual criterion for exponent3."""
import hashlib,itertools,json,signal,time
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C4_HB_RESIDUAL_SOURCE_20261105_v63'
PATHS={
 'reduction':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63/report.json',
 'reduction_peer':'reports/AMBIENT_B18_C_FULL_LINEAR_REDUCTION_SOURCE_20261105_v63_REVIEW/report.json',
 'eta':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63/report.json',
 'eta_peer':'reports/AMBIENT_B18_C_FULL_ETA4_CACHE_VALUES_20261105_v63_REVIEW/report.json',
 'cut':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63/one-CUT1-residual-and-formal-coordinate-rank-certificate.json',
 'cut_peer':'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63_REVIEW/report.json',
 'HB':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json',
 'HB_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json',
 'HB_values':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63/phase2-actual-RIGHT2-quartics-two-factors-and-rank.json',
 'HB_values_peer':'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'full_source':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/report.json',
 'full_operator':'reports/AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63/full-HB6-by2-affine-source-operator-and-old27-leaves.json',
}
ZERO=(0,)*6
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def const(q):return{ZERO:Q(q)}if q else{}
def var(i):return{tuple(int(j==i)for j in range(6)):Q(1)}
def add(a,b):
 v=a.copy()
 for k,q in b.items():
  v[k]=v.get(k,Q(0))+q
  if not v[k]:del v[k]
 return v
def mul(a,b):
 v={}
 for x,c in a.items():
  for y,d in b.items():
   k=tuple(p+q for p,q in zip(x,y));v[k]=v.get(k,Q(0))+c*d
 return{k:q for k,q in v.items()if q}
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()};assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 assert d['reduction']['unconditional_m3C_equals_ell_m2C_proved']
 assert d['eta']['actual_full_C_HF_degree3']==6
 assert d['cut']['actual_C3_HF_if_positive']==[3,6,6,4,0]and d['cut']['cut3_torsion_pattern_if_positive']==[1,1,2,'e>=3']
 assert d['HB_values']['actual_quartic_matrix_rows_A4_B4_columns_HB0_HB1'][0]==['-3/5','6']
 assert d['full_operator']['rows']==['A*r²ell²','A*r ell³','A*ell4','B*r²ell²','B*r ell³','B*ell4']
 # Ordered rows [0,j,2,5] for j=1 or4, columns HB0,HB1,ell4A,ell4B.
 # The last four HB entries p,q,s,t stay arbitrary, not zero-filled.
 x,y,p,q,s,t=[var(i)for i in range(6)]
 mat=[[const(Q(-3,5)),const(6),{},{}],[x,y,{},{}],[p,q,const(1),{}],[s,t,{},const(1)]]
 det={}
 for perm in itertools.permutations(range(4)):
  sign=(-1)**sum(perm[i]>perm[j]for i in range(4)for j in range(i+1,4));v=const(sign)
  for i in range(4):v=mul(v,mat[i][perm[i]])
  det=add(det,v)
 expected=add(mul(const(Q(-3,5)),y),mul(const(-6),x));assert det==expected
 data={
  'scope':'C4=C/ell4C, actual original comparison module over sameR0,D=QQ[[ell]]. Puretheorem, no new pivot03/HB values.',
  'unconditional_C4_m5_zero_proof':'C3m4=0 givesm4C⊂ell³C. Fullreductionm5C=ell*m4C yieldsm5C⊂ell4C. Thereforem5(C4)=0 WITHOUTnewcoefficientorclosurejob.',
  'unconditional_C4_prefix_HF':[3,6,6,6],
  'unconditional_C4_length_values':[23,24],
  'unconditional_C4_HF_values':[[3,6,6,6,2,0],[3,6,6,6,3,0]],
  'length_proof':'The acceptedDtypeisD4+tors(1,1,2,e>=3). ThuslengthC4=16+1+1+2+min(4,e)=23ife3and24ife>=4. Sinceell4C⊂m4C,theHFprefixthrough3isfull3666,length21. Combinedwithm5C4zero,thisgivesexactHF4values2or3.',
  'full_quartic_NF6_order':d['full_operator']['rows'],
  'actual_C4_added_source_initials':{'ell4A_NF':[0,0,1,0,0,0],'ell4B_NF':[0,0,0,0,0,1],
   'ell4H_NF':'0inthisknownquarticNF becauseHplane2alreadyzero; no other fullC relation isadded.'},
  'conditional_HB_residual_test':'Letw=HB1+10HB0. Ifw_1!=0 orw_4!=0 (A*r ell³orB*r ell³), choosej1or4. Theliteral4sourcecolumnsHB0,HB1,ell4A,ell4B onorderedrows[0,j,2,5]have determinant(-3/5)*w_j.',
  'universal_matrix_rows0_j_2_5_columnsHB0_HB1_cutA_cutB':[['-3/5','6','0','0'],['x','y','0','0'],['p','q','1','0'],['s','t','0','1']],
  'universal_determinant_monomials_x_y_p_q_s_t':[[list(k),str(v)]for k,v in sorted(det.items())],
  'universal_determinant':'-3*y/5-6*x=(-3/5)*(y+10*x), independentofp,q,s,t',
  'positive_completion_proof':'A nonzero displayedminor removesatleast4dimensionsfromtheknownfullquarticNF6inC4, soHF4(C4)<=2. Theactualprefix21andm5C4zero giveC4length<=23, whiletheacceptedDVRtypegiveslower23. HenceC4length23/HF36662andexacte3.',
  'conditional_exact_comparison_C_DVR_type':{'free_rank':4,'torsion_exponents':[1,1,2,3],'torsion_length':7},
  'negative_scope':'Ifbothw1andw4arezero,thiscriterionisnegativeonly. FullHBmayhaverank2throughitsell4coordinates; thene3or4remains. IfHBknownrank1,thewholequarticspaceisstillunknown. DoNOTinferexacte>=4orHF4C4=3fromnegativeknownsources.',
  'no_full_C_ell3H_relation_assumed':True,
  'no_new_source_jet_or_coefficient_values':True,
  'no_original_H0_nilcube_or_exact_nilindex_transfer':True,
  'fixed_source_profile':{'matrix_size':[4,4],'universal_symbol_count':6,'max_polynomial_degree':4,'numeric_envelope':25000,'hard_numeric':150000,'hard_seconds':60,'hard_memory_bytes':2147483648,'hard_RAW':33554432},
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'C4-finite-HF-and-two-HB-residual-coordinate-criterion.json';ep.write_text(json.dumps(data,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'PURE_C4_HF36662_OR36663_AND_CONDITIONAL_HB_RESIDUAL_E3_NO_VALUES',
  'unconditional_C4_m5_zero_proved':True,'unconditional_C4_length_23_or24_proved':True,'actual_HB_residual_values_evaluated':False,'remaining_C_exponent_determined':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
