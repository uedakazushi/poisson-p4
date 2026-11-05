#!/usr/bin/env python3
"""One-variable module structure only; no ungraded S-splitting claim."""
import hashlib,json,resource,signal,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_TOLD_TNIL_COKERNEL_R_CHAIN_20261105_v63';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def mm(A,B):return [[sum(A[i][k]*B[k][j]for k in range(5))for j in range(5)]for i in range(5)]
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 D=load('AMBIENT_B18_TOLD_TNIL_COKERNEL_GRADED_R2_20261105_v63');assert D['all_passed']and D['actual_cokernel_HF']==[3,1,1]
 F=load('AMBIENT_B18_TOLD_TNIL_ZERO3_FIRST_MAP_20261105_v63','actual-three-first-columns-and-fixed-quotient-proof.json');assert F['quotient4_by40_after_origin_products'][1]=={'7':'1'}
 mini=load('AMBIENT_B18_HONEST_NIL_T_MINIMAL_OPERATOR_20261104_v63','actual-minimal-Tnil5-849-and-d3-linear-operator.json')
 assert mini['minimal_target_labels'][0]=={'kind':'original_free44_source','index':42}
 J=[[0]*5 for _ in range(5)];J[1][0]=1;J[2][1]=1;J2=mm(J,J);assert J2[2][0]==1 and not any(x for row in mm(J2,J)for x in row)
 # Each zi action can only send the three top generators q,a',b' to r²q.
 # Every elementary slot therefore has product0 with r and all other zi.
 E=[]
 for c in[0,3,4]:
  M=[[0]*5 for _ in range(5)];M[2][c]=1;E.append(M)
 zero=[[0]*5 for _ in range(5)]
 assert all(mm(J,M)==mm(M,J)==zero for M in E)
 assert all(mm(M,N)==zero for M in E for N in E)
 evidence=dict(actual_module='Q=coker(Told -> Tnil)',actual_Q_HF=[3,1,1],actual_Q_length=5,
  literal_q='image of actual minimal free44 row42',r_variable='rA, not DVR ell',
  proof='q,rq,r²q have successive nonzero graded classes. Choose two other top lifts a,b with zero e42 component. Then ra,rb lie in m²Q=k*r²q. Write ra=c_a*r²q and rb=c_b*r²q and replace a,b by a-c_a*rq,b-c_b*rq. They are killed by r and retain independent top classes. The five vectors form a k-basis.',
  one_variable_module='Q is noncanonically k[[rA]]/(rA³) direct sum k direct sum k',
  basis_order=['q','rA*q','rA²*q','a-c_a*rA*q','b-c_b*rA*q'],actual_r_action_in_this_existing_basis=J,
  zi_image_containment='Each zAi*Q is contained in m²Q=k*rA²q',
  consequences=['rA*zAi*Q=0','zAi*zAj*Q=0','rA³*Q=0'],
  full_S_action_template='zi is an unknown linear combination of E_(2,0),E_(2,3),E_(2,4)',
  all_three_template_slot_products_and_r_commutators_checked=True,
  number_of_unevaluated_S_action_scalars=21,
  not_S_module_direct_sum=True,not_ell_DVR_exponent_statement=True,
  no_actual_c_a_c_b_or21_scalars_computed=True,no_source_or_jet=True)
 p=OUT/'actual-Q-one-variable-chain-and-unknown-S-action-template.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 for p2 in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p2.relative_to(ROOT))]=sha(p2)
 rep=dict(id=ID,all_passed=True,frozen=True,actual_Q_one_variable_rA_Jordan_lengths=[3,1,1],actual_Q_HF=[3,1,1],full_S_action_unknown_scalars=21,no_S_module_split_or_ell_DVR_exponents_claim=True,no_new_source_scalar_or_jet=True,input_sha256=pins,evidence_sha256={p.name:sha(p)},seconds=time.monotonic()-st)
 (OUT/'report.json').write_text(json.dumps(rep,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
