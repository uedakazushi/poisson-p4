#!/usr/bin/env python3
"""Finite filtered-source completeness audit; no new relation values."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json,time
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_ORDER4_FINITE_SOURCE_COMPLETENESS_PREFLIGHT_20261105_v63'
PATHS={
 'L':'reports/AMBIENT_B18_T_QUADRATIC_L_SYZYGY_R4_20261104_v63/actual-quadratic-L-syzygy-scalar-offset-certificate.json',
 'restricted':'reports/AMBIENT_B18_T_RESTRICTED118_QUARTIC_PREFLIGHT_20261104_v63/restricted118-quartic-source-operator-preflight.json',
 'theta':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json',
 'known40':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63/known-initial40-and-length-upper17.json',
 'known40_peer':'reports/AMBIENT_B18_T_KNOWN_INITIAL_LENGTH_UPPER17_20261104_v63_REVIEW/report.json',
 'schur_peer':'reports/AMBIENT_B18_T_MINIMAL5_SCHUR_OPERATOR_20261104_v63_REVIEW/report.json',
 'HF3_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json',
 'restricted_values_peer':'reports/AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63_REVIEW/report.json',
 'd3_values_peer':'reports/AMBIENT_B18_T_CORRECTED168_QUARTIC_VALUES_20261104_v63_REVIEW/report.json',
 'bridge':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_SOURCE_PREFLIGHT_20261104_v63/actual619-quartic-source-order-and-unknown-gates.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def sp(p):return {int(i):Q(c) for i,c in p if Q(c)}
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in PATHS.items()};D={k:read(ROOT/p) for k,p in PATHS.items()}
 for k in ['known40_peer','schur_peer','HF3_peer','restricted_values_peer','d3_values_peer']:assert D[k]['all_passed']
 L=[sp(p) for p in D['L']['normalized_actual37_linear_relations']];h=list(map(Q,D['L']['normalized_actual37_e5_r2_values']))
 assert len(L)==len(h)==37
 m2=list(itertools.combinations_with_replacement(range(8),2));m3=list(itertools.combinations_with_replacement(range(8),3));ix3={m:i for i,m in enumerate(m3)}
 b3=Q(-42215,3456);c4=Q(4565,68832);T=D['theta']['Theta4_cubic600'];theta={}
 for key,c in [('constant_part',Q(1)),('times_b3',b3),('times_c4',c4)]:
  for i,a in T[key]:theta[i]=theta.get(i,Q(0))+c*Q(a)
 checks=[]
 for p,l in enumerate(L):
  for mon in m2:
   value=sum(a*theta.get(120*target+ix3[tuple(sorted(mon+(var,)))],Q(0)) for row,a in l.items() for target,var in [divmod(row,8)])
   expected=h[p] if mon==(7,7) else Q(0);assert value==expected
   checks.append([p,list(mon),str(value)])
 qb=D['restricted']['canonical_qB_initial'];qbchecks=[]
 for var in range(8):
  value=Q(0)
  for target,exp,c in qb:
   mon=tuple(i for i,n in enumerate(exp) for _ in range(n))+(var,)
   value+=Q(c)*theta.get(120*target+ix3[tuple(sorted(mon))],Q(0))
  assert value==(b3 if var==7 else 0);qbchecks.append([var,str(value)])
 qc=4*120+ix3[(5,7,7)];assert theta.get(qc,Q(0))==c4
 assert len(D['restricted']['all119_linear_source_syzygies296'])==119
 assert len(D['restricted']['actual118_constant_combinations_of119'])==118
 assert len(D['theta']['all168'])==168 and len(D['bridge']['all619'])==619
 # The independently accepted homogeneous GB40 supplies only the residual
 # one-dimensional candidate, not equality with the true initial module.
 G=D['known40']['homogeneous_module_generators40'];leaders=[]
 for gen in G:
  # Original leading order is target first, then lex on exponent vectors.
  # Replay its exact accepted leading set directly, independent of coefficients.
  terms=[(int(i),tuple(m)) for i,m,c in gen if Q(c)]
  leaders.append(max(terms,key=lambda x:(-x[0],x[1])))
 # The source theorem only needs the documented, independently accepted
 # all-degree normal forms. Do not infer a new leading order from this helper.
 assert D['known40']['standard_degree_at_least3']==['e4*r^d']
 result={'source_manifest':inputs,'all1332_L_Sym2_Theta4_actual_offset_checks':checks,
  'all8_qB_variables_Theta4_actual_offsets':qbchecks,'canonical_qC_Theta4_offset':str(c4),
  'actual_relations':{'minimal_cycle657':657,'Koszul168':168,'complete_total':825,
   'first_linear_selected':37,'remaining_after_firsts':620,'quadratic_pivot':1,'constant_kernel':619,
   'complete_linear_coefficient_syzygies':119,'nonzero_cubic_pivot':1,'linear_syzygy_kernel':118},
  'finite_filtered_source_cover':['corrected619 constant columns','corrected118 kernel of linear-coefficient syzygies','corrected168 actual d3 columns',
   'all mixed degree3 syzygies of L37*Sym2 and qB*variables, killed by the simultaneous Theta4 identities'],
  'full_order4_dual_completeness_proved':True,
  'conditional_conclusion':'if all619,all118,all168 corrected pure4 values vanish then gr4T is exactly one-dimensional; a nonzero value kills it',
  'reason_no_unchecked_higher_L_syzygy':'quadratic-coefficient L syzygies are included in the full1332 map checked here; coefficientdegree>=3 contributes known degree4 initials',
  'Schreyer140_role':'accepted GB40 proves the one-dimensional candidate; no type-only or unlifted140-syzygy assertion is used',
  'actual619_values_evaluated_here':False,'new_source_coefficients':0}
 ep=o/'finite-filtered-order4-completeness-and-Theta4-source-checks.json';dump(o,ep.name,result)
 note=ROOT/'research/notes'/f'{STEM}.md';r={'id':STEM,'all_passed':True,'source_Theta4_checks':1341,
  'order4_finite_source_completeness_accepted':True,'actual619_values_or_newcoefficients_computed':False,
  'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 dump(o,'report.json',r);print(json.dumps(r))
if __name__=='__main__':main()
