#!/usr/bin/env python3
"""Independent all3960 factor identities and allorder descent proof; no values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','cert':f'reports/{A}/actual3960-rho-factor-identities-or-obstructions.json',
'input':'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63/actual660-monomial-NFs-and-combined-rho-raw33-factor-input.json',
'input_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_INPUT_GATE_20261105_v63_REVIEW/report.json',
'raw':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
'raw_peer':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63_REVIEW/report.json',
'adjoint_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def dot(a,b):return sum((c*b.get(i,Q(0)) for i,c in a.items()),Q(0))
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in ['input_peer','raw_peer','adjoint_peer'])
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 I=D['input'];C=D['cert'];nf=list(map(sp,I['all660_actual_raw33_normal_forms']));weights=[[sp(p) for p in q] for q in I['two_three_component_original_Taylor_functionals']];bars=[[sp(p) for p in q] for q in I['candidate_two_three_component_raw33_covectors']]
 assert [[enc(p) for p in q] for q in bars]==C['actual_two_three_component_bar_rho33'];checks=0
 for i in range(2):
  for target in range(3):
   for j,p in enumerate(nf):assert dot(bars[i][target],p)==weights[i][target].get(j,Q(0));checks+=1
 assert checks==3960 and C['failures']==[] and C['combined_two_rho_factor_raw33']
 omitted=[tuple(w.count(i) for i in range(9)) for d in range(5) for w in itertools.combinations_with_replacement(range(9),d) if w.count(8)>=2];assert len(omitted)==55 and [list(m) for m in omitted]==C['omitted55_ell_degree_ge2_raw33_zero_and_rho_identically_zero']
 actions=[list(map(sp,A)) for A in D['raw']['coordinate_actions9_sparse_columns']]
 for m in omitted:
  p={0:Q(1)}
  for i,n in enumerate(m):
   for _ in range(n):p=mv(actions[i],p)
  assert not p
 ep=o/'independent-all3960-component-identities-and715-monomial-descent.json';dump(o,ep.name,{'inputs':inputs,'all3960_factor_identities_verified':True,'all55_omitted_ell2_monomials_killed':True,'all715_degree_le4_monomials_accounted':True,'actual_two_three_component_bar_rho33':C['actual_two_three_component_bar_rho33'],
'formal_allorder_proof':'TheacceptedactualcoordinatequotientR0/(IAIB,ell²) israw33andannihilatesm5. Originalrho is a finiteTaylorfunctional supportedindegree<=4, soalsoannihilatesm5. All715 ordinarymonomials throughdegree4,including55ell²terms,haveequalvaluesunderoriginalrhoandcandidatebar_rho composedwiththeactualquotient. Theircontinuouslinear mapsagreeonALLformalfunctions. Thisiscomponentwiseforallthreephysicaltargetrows, notmerelyalongonecurve.',
'actual_invariant_source_formula':'ForactualPhi=V56*Wv, thepeered286-coordinateadjointgives(alpha,beta)=bar_rho*(Phi inraw33). ActualprimitivePhi muststillbeboundtosavedphysical12/monicsbeforeanydot.',
'A23_graph_callbacks_remain_separate':True,'actual_Phi_dot_invariants_M1_or_graph_values_computed':False,'new_source_or_implicit_or_gauge_coefficients':0})
 report={'id':STEM,'all_passed':True,'both_combined_functionals_descend_actual_raw33_allorders':True,'all3960_identities_and715_monomials_verified':True,'actual_invariant_or_M1_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
