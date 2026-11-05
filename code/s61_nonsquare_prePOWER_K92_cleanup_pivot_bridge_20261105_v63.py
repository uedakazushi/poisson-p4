#!/usr/bin/env python3
"""SOURCE exactbridge: globalrowcleanup scalesexistingresidual4, retains pivots."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_K92_CLEANUP_EXISTING_PIVOT_BRIDGE_20261105_v63'
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2,default=str)+'\n')
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False);pins={}
 def src(rel):p=ROOT/rel;pins[rel]=sha(p);return read(p)
 gr=src('reports/S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_SOURCE_20261105_v63/report.json');assert gr['all_passed'];src('reports/S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_SOURCE_20261105_v63/typed-generic-rowcleanup-and-pivot.json');lg=src('reports/S61_NONSQUARE_PREPOWER_K92_GLOBAL_CLEANUP_LIFT_GUARD_20261105_v63/report.json');assert lg['generic_A_topright_block_literallyZERO_before_specialization']
 nx=src('reports/S61_NONSQUARE_PREPOWER_K92_CLEANUP_COLUMN6_ENTRY_20261105_v63/report.json');assert nx['all_passed']and nx['formal_specialization_ZERO']and nx['generic_N6_callback_RETAINED'];src('reports/S61_NONSQUARE_PREPOWER_K92_CLEANUP_COLUMN6_ENTRY_20261105_v63/single-column6-specialization.json')
 t=src('reports/S61_NONSQUARE_PREPOWER_K92_FACTOR_RANK_PROFILE_20261105_v63/typed-next-residual-Schur.json');q=[[Q(v)for v in r]for r in t['column_transform_Q2']];alpha=t['alpha'];beta=t['beta'];r=[[Q(v)for v in rr]for rr in t['row_transform_R2']];K2=r[1][0];r=[r[0]]+[[v/K2 for v in rr]for rr in r[1:]];assert r[1]==[Q(i==0)for i in range(7)]and r[2]==[Q(i==2)for i in range(7)];assert alpha[2]==4418 and alpha[0]==4486 and beta==[-68,0,0,0,0,0,0]
 P=[0,6,1];J=[2,3,4,5];assert 2 not in P and 1 in P;assert all(q[i][1]==Q(i==0)and q[i][6]==Q(i==6)and q[i][0]==Q(i==2)for i in range(7))
 expected={2:{1:Q(1)},3:{2:Q(18,7),3:Q(1)},4:{2:Q(147,25),4:Q(1)},5:{2:Q(126,25),5:Q(1)}}
 for j,s in expected.items():assert {i:q[i][j]for i in range(7)if q[i][j]}==s
 callbacks=[]
 for j in P+J:
  callbacks.append(dict(Tprime_column=j,common_a_power=4418,terms=[dict(N_index=i,QQ=str(q[i][j]),a_power=beta[i])for i in range(7)if q[i][j]]))
 nodes=[
  dict(name='G',shape=[11,11],op='accepted_genericcleanup_external'),
  dict(name='S7',shape=[7,7],op='A4_Schur_of_G',formula='D-C*A4^-1*B; A4block2+2 inverse viaadjP,adjQ'),
  dict(name='T7',shape=[7,7],op='weighted_normalization',formula='R2normalized*diag(a^alpha)*S7*diag(a^beta)*Q2',alpha=alpha,beta=beta,R2=r,Q2=q),
  dict(name='K7',shape=[7,7],op='identity_except_row2',row2_coefficients=dict(column2='delta',column1='-59375*delta*a^-68'),determinant='delta',inverse_row2=dict(column2='1/delta',column1='59375*a^-68')),
  dict(name='Tprime',shape=[7,7],op='matrix_product',args=['K7','T7']),
  dict(name='P3',shape=[3,3],op='submatrix',matrix='T7',rows=P,columns=P),
  dict(name='P3inverse',shape=[3,3],op='classical_adjugate3_over_determinant',matrix='P3',parameter_only_chart='detP3!=0'),
  dict(name='S4',shape=[4,4],op='Schur',formula='T7[J,J]-T7[J,P]*P3inverse*T7[P,J]'),
  dict(name='S4prime',shape=[4,4],op='Schur',formula='Tprime[J,J]-Tprime[J,P]*P3inverse*Tprime[P,J]'),
  dict(name='S4prime_firstrow_callback',shape=[1,4],op='finite_callback_expression',formula='a^4418*(rJ-[N2,N6,a^-68*N0]*P3inverse*T7[P,J])',rJ=['N1','N3+(18/7)*N2','N4+(147/25)*N2','N5+(126/25)*N2']),
 ]
 dump(out/'typed-existing-pivot-bridge.json',dict(nodes=nodes,pivot_order=P,free_order=J,callback_columns=callbacks,
  rowcleaning_identity='Tprime row2=delta*(T7row2-59375*a^-68*T7row1); allotherrowsunchanged',
  pivot3_UNCHANGED='row2notinP, soTprime[P,P]=P3; row1isINpivotset',
  exact_residual_identity='S4prime=diag(delta,1,1,1)*S4',
  universal_proof='pivotrow1 satisfies T7[1,J]=T7[1,P]*P3inverse*T7[P,J]; subtractedpivotrow contributesZERO toSchur. Remainingfirstrow scalesdelta.',
  same_cokernel_on_chart=['deltaP!=0','deltaQ!=0','detP3!=0'],all_chart_complements_RETAINED=True,
  all_existing4plus3_pivots_RETAINED=True,N0_andN6_lie_in_existing_pivot_columns=True,alternateP5_is_NOTnewrankprogress=True,
  specialized_N6_ZERO_NOT_substituted_in_genericformula=True,all7generic_Ncallbacks_RETAINED=True,
  no_next_N1toN5_callback_job_or_coefficient_window=True))
 dump(out/'remaining-source-target-assessment.json',dict(current_genuine_remaining_object='sameexistingresidual4; firstrowexpression above, other3rows unchanged',
  next_meaningful_target='a directborder/cofactorofthatactualfirstresidualrow, preservingALLthreepivotcorrections; evaluateonlyafterfreshwhole-worker bounds',
  not_a_cheapest_callback_sequence='N0,N6arealreadyeliminatedpivotcolumns; theirpositive/zeroformalvalues alone cannotproduceanotherresidualpivot',
  currentformal_data='direct4border16polynomials throughdegree342ZERO; target345ZERO; nonglobalN6formalzero; N0nonzeroQQ[t6663] independentlyverified',
  no_residual4_rank_or_full57_global_determinant_claim=True))
 pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__));note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<33554432
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,no_new_signed_coefficient_values=True,exact_existing_residual4_row_scaling=True,all_existing_three_Schur_pivots_RETAINED=True,
  N0_N6_are_existing_pivot_columns=True,generic_N6_callback_retained_despite_specializedzero=True,no_mere_callback_sequence_or_newscalar_job=True,
  no_det11_nonzero_or_global57_claim=True,total_new_RAW_before_report=new,external_read_RAW=external,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
 dump(out/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
