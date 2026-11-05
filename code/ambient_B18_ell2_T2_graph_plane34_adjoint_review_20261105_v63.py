#!/usr/bin/env python3
"""Independent verification of saved104-coordinate adjoint obstructions."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_ADJOINT_20261105_v63';STEM=A+'_REVIEW'
P={'primary':f'reports/{A}/report.json','cert':f'reports/{A}/actual-two-plane-source-adjoint-certificates.json','input':'reports/AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63/actual-plane-G1-degree23-and-twelve-adjoint-source-input.json','input_peer':'reports/AMBIENT_B18_ELL2_T2_GRAPH_PLANE34_INPUT_GATE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def sp(p):return {i:Q(c) for i,c in p if Q(c)}
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def dot(a,b):return sum((c*b.get(i,Q(0)) for i,c in a.items()),Q(0))
def main():
 t=time.monotonic();o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False);inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['input_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 cols=list(map(sp,D['input']['twelve_source_adjoint_columns104']));psi=list(map(sp,D['input']['two_invariant_functionals_on104_with_F2_zero']));cert=D['cert'];verified=[]
 for key,n in [('degree3_only_results',6),('complete_results',12)]:
  for i,r in enumerate(cert[key]):
   d=sp(r['obstruction104']);assert not r['representation_exists'] and all(dot(c,d)==0 for c in cols[:n]) and dot(psi[i],d)==Q(r['nonzero_value'])!=0
   verified.append({'outputs':n,'invariant':i,'obstruction':r['obstruction104'],'nonzero_pairing':r['nonzero_value']})
 assert all(any(i>=52 for i,c in r['obstruction104']) for r in cert['complete_results'])
 ep=o/'independent104-nuisance-aware-plane-adjoint-obstructions.json';dump(o,ep.name,{'inputs':inputs,'all_four_saved_obstruction_identities':verified,'all12_full104_rows_checked':True,'complete_obstructions_have_nonzero_F2_nuisance_component':True,'scope':'NeitherPhi3onlynorPhi3+Phi4planeoutputsdeterminestheseinvariantsintheunconstrained104-coefficientmodel. Thisdoesnotprovenewhigherphysicalvaluesnecessary: additionalmixedzsourceconstraintsmaycutdownambiguity.','actual_M1_M2_Phi3_Phi4_or_graph_values_not_evaluated':True})
 report={'id':STEM,'all_passed':True,'all_four_plane_adjoint_nonexistence_certificates_verified':True,'actual_invariant_or_M1_M2_values_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};dump(o,'report.json',report);print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json')}))
if __name__=='__main__':main()
