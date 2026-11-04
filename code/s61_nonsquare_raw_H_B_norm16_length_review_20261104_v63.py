#!/usr/bin/env python3
"""One independent cofactor00 + PID exact-length corollary audit."""
import argparse,gzip,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_LENGTH_20261104_v63_REVIEW'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 ib=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63';cb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_CANONICAL_20261104_v63';peer=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_20261104_v63_REVIEW';lb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_LENGTH_20261104_v63';rb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63';eb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_EXCEPTION_SQUAREFREE_20261104_v63';qb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_Q_LINEAR_BATCH_20261104_v63'
 assert sha(peer/'report.json')=='14ec6552b7f64be98c26936f7324e282d9652b59f8a475ef7dae1568d5b3cba3'and read(peer/'report.json')['all_passed']
 c=read(ib/'certificate.json');n=read(cb/'native-certificate.json.gz');D=int(c['determinant_cost_bounds_only']['matrix_uniform_QQ_denominator']);matrix=[]
 for row in c['matrix_rows_source_order']:
  rr=[]
  for p in row:
   z={}
   for m,s in p:
    assert m[1:]==[0,0]and m[0]>=0;q=Q(s)*D;assert q.denominator==1;z[m[0]]=q.numerator
   rr.append(z)
  matrix.append(rr)
 v0={}
 for m,s in n['only16_solution_cofactors'][0]:
  assert m[1:]==[0,0]and m[0]>=0;q=Q(s);assert q.denominator==1;v0[m[0]]=q.numerator
 assert v0
 inp=out/'cofactor-source-integer.txt'
 with inp.open('w')as f:
  for p in [p for row in matrix for p in row]+[v0]:
   f.write(str(len(p))+'\n')
   for k,v in sorted(p.items(),reverse=True):f.write(str(k)+' '+str(v)+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_norm16_cofactor_review_20261104_v63.cpp';exe=work/'cofactor-check'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','cofactor-review',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 native=read(out/'native-certificate.json');assert native['cofactor00_equals_saved_V0']and {int(k):int(v)for k,v in native['cofactor00']}==v0
 # Exact det identity follows from universal adjugate identity plus prior fullNV.
 assert native['all_Bareiss_divisions_reverse_multiplied']and native['exact_divisions']==sum(k*k for k in range(1,15))==1015
 delta=n['fflu_solver_delta'];degree=max(m[0]for m,v in delta);valuation=min(m[0]for m,v in delta);assert degree==148 and valuation==31
 old=read(rb/'certificate.json');exc=read(eb/'certificate.json');H=old['curve_cover_monic_H_d4'];assert H==exc['source_finite_cover_inputs']['H'];assert max(m[2]for m,v in H)==4
 assert [(m,Q(v))for m,v in H if m[2]==4]==[([0,0,4],Q(1))]and all(m[1]>=0 and m[2]>=0 for m,v in H)
 # The only Laurent powers in H are a, which is already a unit in the base.
 Qcover=exc['source_finite_cover_inputs']['Q'];assert [(m,Q(v))for m,v in Qcover if m[0]==2]==[([2,0,0],Q(1))]and max(m[0]for m,v in Qcover)==2 and all(min(m)>=0 for m,v in Qcover)
 qreport=read(qb/'report.json');assert qreport['all_five_linear_Q_components_closed_including_previous']and sum(qreport['untested_factor_degrees'])==32
 lreport=read(lb/'report.json');assert lreport['exact_old_projection_length']==117 and lreport['exact_old_Hcover_length']==468 and lreport['updated_geometric_Bbad_bound']==724
 proof=dict(one_cofactor00_independently_proved=True,no_full16_determinant_or_solver_reexecution=True,det_identity='adj(N)*(N V)=det(N)V=delta adj(N)e0. Firstcomponent gives det(N)V0=delta cofactor00. Sincecofactor00=V0!=0 inQQ[a], det(N)=delta.',source_M_really_multiplication_by_k1='Bound to prior independent16source division checks; h1monic16, no negative a powers in M.',N_equals_D_M_with_D_QQunit=True,PID_Smith_length=148,a_primary_length=31,exact_a_local_projection_QQlength=117,monic_H_over_aunit_d_degree=4,exact_old_H_cover_QQlength=468,remaining_Q_e_degree_sum=32,remaining_Q_cover_atmost_32_times_2_times_4=256,total_GEOMETRIC_Bbad_points_atmost724=True,imported_other_containment_and_unit_gates='ExistingEXCEPTION_SQUAREFREEfinitecover and acceptedQ_LINEAR_BATCH/Lclosed results. No new proof ofBbad existence, emptiness, reducedness, orWHOLEBbadlength.',not_global_scheme_length724=True,source_no_negative_a_exponents_M=True,only_H_Laurent_a_allowed=True)
 (out/'proof-certificate.json').write_text(json.dumps(proof,indent=2)+'\n')
 failure=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM16_COFACTOR_COMPILE_FAILURE_20261104_v63'
 assert read(failure/'failure.json')['cofactor_computation_started']is False and sha(inp)==sha(failure/'cofactor-source-integer.txt')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',peer/'report.json',ib/'report.json',ib/'certificate.json',cb/'report.json',cb/'native-certificate.json.gz',lb/'report.json',lb/'certificate.json',rb/'report.json',rb/'certificate.json',eb/'report.json',eb/'certificate.json',qb/'report.json',failure/'failure.json',failure/'compile.log',failure/'failed-source.cpp']
 report=dict(id=ID,all_passed=True,one_own_Bareiss15_cofactor00=True,all1015_exact_divisions_remultiplied=True,actual_det_equals_delta_independently_proved=True,exact_QQ_lengths148_117_468=True,geometric_Bbad_bound724_only=True,no_new_point_prime_GB_full16det_or_authorworker=True,native_seconds=native['native_seconds'],input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},executable_sha_not_published=sha(exe),seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();main(a.output,a.work)
