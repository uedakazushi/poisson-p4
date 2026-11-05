#!/usr/bin/env python3
"""Independent explicit286-coordinate adjoints and180-output obstructions."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,mv
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','cert':f'reports/{A}/actual-two-full286-mixed-source-adjoint-certificates.json','input':'reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_INPUT_GATE_20261105_v63/actual-full-G1-and186-mixed-adjoint-source-input.json','input_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_INPUT_GATE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def dot(a,b):return sum((c*b.get(i,Q(0)) for i,c in a.items()),Q(0))
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['input_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 S=D['input'];C=D['cert'];cols=list(map(sp,S['actual186_source_adjoint_columns286']));psi=list(map(sp,S['two_invariant_functionals286']));assert C['output186_order']==S['output186_order'] and len(cols)==186
 for i,r in enumerate(C['degree23_only_results']):
  q=sp(r['obstruction286']);assert not r['representation_exists'] and all(not dot(q,p) for p in cols[:180]) and dot(q,psi[i])==Q(r['nonzero_value'])!=0
 used=set();supports=[]
 for i,r in enumerate(C['complete_results']):
  rho=sp(r['rho']);assert r['representation_exists'] and mv(cols,rho)==psi[i]
  assert all(k<234 for k in mv(cols,rho)) # all52 nuisanceF2 coefficients cancel exactly
  used.update(rho);supports.append(len(rho))
 assert supports==[11,16]
 inventory=[[i,S['output186_order'][i]] for i in sorted(used)];bydegree={d:[r for r in inventory if r[1][0]==d] for d in [2,3,4]};assert list(map(len,bydegree.values()))==[6,6,4]
 assert all(r[1][1]==1 for r in bydegree[2]+bydegree[3]) and {r[1][1] for r in bydegree[4]}=={0,1}
 ep=o/'independent-full286-adjoints-and-minimal-Phi-callback-inventory.json';dump(o,ep.name,{'inputs':inputs,'both180_output_nonexistence_obstructions_verified':True,'both_explicit186_covectors_verified_on_all286_domain_slots':True,'all52_F2_nuisance_coefficients_exactly_zero':True,'actual_rho_support_counts':supports,'exact_selected_source_callback_inventory':inventory,
'source_value_identity':'alpha=rho0*Phi,beta=rho1*Phi, withPhi2all/selectedPhi3/Phi4plane outputs inacceptedorder. TheidentityisvalidforanyactualformalF=M*Wv withF0=0 becauseallnuisancecoordinateswereincluded.',
'what_is_not_decided':'Neitheralpha norbeta norcanonicalM1/graphvaluesareevaluated. Failurewith180outputsonlyisamethodlimitation, notproofthatnewphysicalsourcecoefficientjobsarenecessary. RequiredmixedPhi3/Phi4callbacksmustfirstbecoherentlyboundtocachesornewauthorizedfunctionalinputs.',
'new_source_or_adjoint_solve_or_Phi_values_computed':False})
 report={'id':STEM,'all_passed':True,'both_full286_adjoint_certificates_and180_obstructions_accepted':True,'actual_invariant_M1_graph_or_Phi_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'rho_support':supports}))
if __name__=='__main__':main()
