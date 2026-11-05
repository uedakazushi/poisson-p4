#!/usr/bin/env python3
"""Independent exact polynomial-plane G1 and104-slot nuisance-aware adjoint gate."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual-plane-G1-degree23-and-twelve-adjoint-source-input.json',
'graph':'reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json',
'graph_peer':'reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63_REVIEW/report.json',
'negative_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_INVARIANT_ADJOINT6_20261105_v63_REVIEW/report.json',
'T':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
'T_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
'active':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
'res':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
'res_peer':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json',
'actions':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
'actions_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def poly(p):return {(i,tuple(m)):Q(c) for i,q in p for m,c in q if Q(c)}
def shift(p,m):return {(i,(a+m[0],b+m[1])):c for(i,(a,b)),c in p.items() if b+m[1]<2}
def encode(p):return [[i,list(m),str(c)] for(i,m),c in sorted(p.items())]
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in D if k.endswith('peer'))
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['gate'];actions=[list(map(sp,M)) for M in D['actions']['C2_coordinate_actions9_sparse_columns']];lifts=list(map(poly,D['res']['physical_lift_m0_columns']));eta=D['res']['physical_target3_basis_positions'];basis=[{i:Q(1)} for i in range(14)]
 # These are actual monomial lifts, not merely a formal module type.
 for j,p in enumerate(lifts):
  value={}
  for(i,(r,e)),c in p.items():
   q={eta[i]:c}
   for _ in range(r):q=mv(actions[7],q)
   for _ in range(e):q=mv(actions[8],q)
   value=add(value,q)
  assert value==basis[j]
 active=D['active']['active_variables'];plane=[]
 for k in D['T']['active_F1_free17']:
  block,j=divmod(k,14);variable=8 if block==0 else active[block-1]
  p=shift(lifts[j],(1,0)) if variable==7 else shift(lifts[j],(0,1)) if variable==8 else {}
  for i,c in actions[variable][j].items():p=add(p,scale(lifts[i],-c))
  assert all(sum(m) in [2,3] for i,m in p);plane.append(p)
 plane.extend({} for _ in range(9));assert list(map(encode,plane))==E['exact_plane_G1_polynomial_columns']
 assert [sorted({sum(m) for i,m in p}) for p in plane]==E['plane_degrees_per_source_column']
 columns=[{} for _ in range(12)]
 for j,p in enumerate(plane):
  for(i,(r,e)),c in p.items():
   degree=r+e
   for var,(dr,de) in enumerate([(1,0),(0,1)]):
    if e+de<2:columns[(degree-2)*6+2*i+e+de]=add(columns[(degree-2)*6+2*i+e+de],{2*j+var:c})
   if degree==2:
    for var,(dr,de) in enumerate([(2,0),(1,1)]):
     if e+de<2:columns[6+2*i+e+de]=add(columns[6+2*i+e+de],{52+2*j+var:c})
 assert list(map(enc,columns))==E['twelve_source_adjoint_columns104']
 inv=[]
 for p in D['graph']['two_invariant_functionals_on_all234_F1_coefficients']:
  q={}
  for k,c in sp(p).items():j,var=divmod(k,9);assert var in [7,8];q[2*j+var-7]=c
  inv.append(q)
 assert list(map(enc,inv))==E['two_invariant_functionals_on104_with_F2_zero'] and all(k<52 for p in inv for k in p)
 assert E['whole_twelve_adjoint_slot_bound']==10858<150000 and E['twelve_adjoint_solve_not_started']
 ep=o/'independent-exact-plane-G1-and-full104-adjoint-input.json';dump(o,ep.name,{'inputs':inputs,'all14_actual_physical_monomial_lifts_evaluate_to_C2_basis':True,'all26_exact_plane_G1_polynomial_columns_reconstructed':True,'actual_plane_degree_support':[2,3],'all12_adjoint_columns104_and_two_F1_only_targets_verified':True,
'allorder_current_frame_proof':'Rawtargetprojectionp0 uses theactual14monomiallifts, p0eta=I3. Theaccepted11polynomialsourcecorrections satisfyd1*a11=K=I-eta*p0, hencep0*d1*a11=0. Subsequent selectedd2 corrections are exactsourcecycles, sod1*d2=0 kills them. Free17imagesarethereforeEXACTLYp0*d1raw*Efree17 in thecurrent26frame. Theother9zero-variablecolumns vanishonplane. NoinverseTaylorcoefficientisintroduced.',
'exact_source_coefficient_equations':'Setz7=0 andF=M*Wv. F0=0,G1plane1=0. Degree3:Phi3=G1plane2*F1. Degree4:Phi4=G1plane3*F1+G1plane2*F2. Thefull52F2slotsarenuisanceunknowns, neverzero-assigned.',
'required_future_adjoint_check':'Each12-entryrho must represent its target onall104slots, includingzero coefficientsonALL52F2slots. Only thenmayPhi3/Phi4functionalvalues be requested.',
'source_availability':'Phi3is genuinecommon105plane3 alreadycached. Phi4r4/r3ellremainsanuncomputedsourcecallback; B56truncationdoesnotfillitwith0.',
'no_new_Phi4_M1_M2_source_or_adjoint_values_computed':True,'whole_twelve_adjoint_slot_bound':10858})
 report={'id':STEM,'all_passed':True,'exact_current_frame_plane_G1_and104_slot_adjoint_gate_accepted':True,'actual_Phi4_M1_M2_or_adjoint_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
