#!/usr/bin/env python3
"""Independent full polynomial G1 and selected mixed coefficient closure audit."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,itertools,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_MIXED_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','gate':f'reports/{A}/actual-full-G1-and186-mixed-adjoint-source-input.json',
'graph':'reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json',
'plane':'reports/AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63/actual-plane-G1-degree23-and-twelve-adjoint-source-input.json',
'plane_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63_REVIEW/report.json',
'T':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
'active':'reports/AMBIENT_B18_ELL2_C2_SECOND_BETTI_ACTIVE_PREFLIGHT_20261105_v63/active-five-d3-source-split-and-budget.json',
'res':'reports/AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63/actual-C2-resolution-and11-unit-source-cancellation.json',
'action':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
'action_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def mon(word):return tuple(word.count(i) for i in range(9))
def encode(p):return [[i,list(m),str(c)] for(i,m),c in sorted(p.items())]
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['plane_peer']['all_passed'] and D['action_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['gate'];actions=[list(map(sp,M)) for M in D['action']['C2_coordinate_actions9_sparse_columns']];lifts=[]
 for col in D['res']['physical_lift_m0_columns']:
  p={}
  for target,v in col:
   for(r,e),c in v:
    if Q(c):p[target,(0,)*7+(r,e)]=Q(c)
  lifts.append(p)
 full=[]
 for k in D['T']['active_F1_free17']:
  block,j=divmod(k,14);var=8 if block==0 else D['active']['active_variables'][block-1];p={}
  for(target,m),c in lifts[j].items():
   n=list(m);n[var]+=1
   if n[8]<2:p[target,tuple(n)]=c
  for i,c in actions[var][j].items():p=add(p,scale(lifts[i],-c))
  full.append(p)
 for z in D['active']['zero_actions']:
  for target in range(3):full.append({(target,mon((z,))):Q(1)})
 assert len(full)==26 and list(map(encode,full))==E['exact_full_G1_polynomial_columns']
 eta=D['res']['physical_target3_basis_positions']
 for j,p in enumerate(full):
  assert all(1<=sum(m)<=3 and sum(m[:7])<=1 for target,m in p)
  linear={target*9+m.index(1):c for(target,m),c in p.items() if sum(m)==1};assert linear==sp(D['graph']['actual_G1linear3_by26_columns'][j])
  plane={(target,(m[7],m[8])):c for(target,m),c in p.items() if not any(m[:7])};assert plane=={(target,tuple(m)):Q(c) for target,m,c in D['plane']['exact_plane_G1_polynomial_columns'][j]}
  value={}
  for(target,m),c in p.items():
   q={eta[target]:c}
   for k,e in enumerate(m):
    for _ in range(e):q=mv(actions[k],q)
   value=add(value,q)
  assert not value
 mons2=[mon(w) for w in itertools.combinations_with_replacement(range(9),2) if w!=(8,8)];m3=[mon((z,7,7)) for z in range(7)]+[mon((z,7,8)) for z in range(7)]+[mon((7,7,7)),mon((7,7,8))];m4=[mon((7,7,7,7)),mon((7,7,7,8))]
 outputs=[(d,i,m) for d,ms in [(2,mons2),(3,m3),(4,m4)] for i in range(3) for m in ms];assert [[d,i,list(m)] for d,i,m in outputs]==E['output186_order'];oi={k:i for i,k in enumerate(outputs)};cols=[{} for _ in outputs]
 for j,p in enumerate(full):
  for(target,m),c in p.items():
   for var in range(9):
    n=list(m);n[var]+=1;k=(sum(n),target,tuple(n))
    if k in oi:cols[oi[k]]=add(cols[oi[k]],{9*j+var:c})
   for i,w in enumerate([(7,7),(7,8)]):
    n=list(m)
    for var in w:n[var]+=1
    k=(sum(n),target,tuple(n))
    if k in oi:cols[oi[k]]=add(cols[oi[k]],{234+2*j+i:c})
 assert list(map(enc,cols))==E['actual186_source_adjoint_columns286']
 assert D['graph']['two_invariant_functionals_on_all234_F1_coefficients']==E['two_invariant_functionals286']
 # Exhaustively check that every omitted F2/F3 monomial cannot contribute
 # to the chosen186 outputs; no hidden nuisance coefficient was setzero.
 plane2={mon((7,7)),mon((7,8))};omitted2=[m for m in mons2 if m not in plane2];all3=[mon(w) for w in itertools.combinations_with_replacement(range(9),3) if w.count(8)<2];checks=0
 for p in full:
  for(target,m),c in p.items():
   for n in omitted2:
    q=tuple(a+b for a,b in zip(m,n));assert (sum(q),target,q) not in oi;checks+=1
   if sum(m)==1:
    for n in all3:
     q=tuple(a+b for a,b in zip(m,n));assert (sum(q),target,q) not in oi;checks+=1
 assert sum(E['slot_breakdown'].values())==118371==E['whole_known_array_slot_bound']<150000 and E['no_adjoint_solve_or_new_source_values']
 ep=o/'independent-full-G1-and-closed286-by186-source-system.json';dump(o,ep.name,{'inputs':inputs,'all26_exact_current_frame_full_G1_polynomials_and_module_relations_verified':True,'all_linear_and_plane_overlaps_verified':True,'all186_adjoint_columns_on286_domain_reconstructed':True,'omitted_F2_F3_monomial_checks':checks,
'complete_nuisance_scope_proof':'G1haspositiveorder,linearpartonlyz. Phi2usesonlyF1. InselectedPhi3(zi*r²,zi*r*ell,plane),G1lin*F2withanyzhasatleast2zanddoesnotcontribute. InplanePhi4,G1lin*F3andanyF2withzhavezandvanish. All234F1and52planeF2coefficientsremainindependent. F0=0 is thepreviouslyverifiedM0Wv condition.',
'source_values_not_assumed':'Phi2isincachedB56. MixedPhi3mustbecomparedfunctionallywithselected135/source33beforevalues. Phi4r4/r3ellremainsformal. Truncationdoesnotassignunavailablecallbacks0.',
'whole_known_array_slot_bound':118371,'no_new_physical_source_implicit_gauge_M1_M2_or_adjoint_values':True})
 report={'id':STEM,'all_passed':True,'exact_fullG1_and_closed286_domain186_output_gate_accepted':True,'actual_adjoint_or_source_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds']}))
if __name__=='__main__':main()
