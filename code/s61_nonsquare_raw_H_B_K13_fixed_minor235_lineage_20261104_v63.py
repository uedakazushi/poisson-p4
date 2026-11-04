#!/usr/bin/env python3
"""Source-only guard: fixedminor235 is the already-used g, so its norm is zero."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_FIXED_MINOR235_LINEAGE_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):return {tuple(m):Q(v)for m,v in p}
def mul(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():
   k=tuple(i+j for i,j in zip(m,n));z[k]=z.get(k,Q(0))+c*d
 return{k:v for k,v in z.items()if v}
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 b=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';p=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_INPUT_20261104_v63';n=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PAIR_NORM_20261104_v63';r=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63';cost=ROOT/'reports/S61_NONSQUARE_RAW_H_B_MINORS_COST_20261104_v63'
 c=read(b/'certificate.json');native=read(b/'native-certificate.json');index=list(itertools.combinations(range(7),3)).index((2,3,5));assert index==26 and c['minor35_to_distinct_generator_index'][index]==15
 rec=native['records'][index];assert rec['rows7']==[2,3,5]and native['all35_source_scalings_verified'];minor=poly(rec['prepared_integer_minor']);gen=poly(c['localized_ideal_generators_monic_QQ_after_namedunit_and_QQ_only_scaling'][15]);assert poly(read(p/'certificate.json')['g'])==gen
 k=max(gen);ratio=minor[k]/gen[k];assert minor=={m:ratio*v for m,v in gen.items()}
 norm=read(n/'native-certificate.json');res=read(r/'certificate.json');assert mul(poly(res['common_G']),poly(res['oldh_quotient_h1']))==poly(norm['remaining_h_coprime_U'])
 pre=read(cost/'certificate.json');source_scaling=dict(row_units=[pre['known_row_unit_divisors'][i]for i in[2,3,5]],column_units=pre['known_column_unit_divisors'],column_QQ_multipliers=pre['integer_column_multipliers_QQ_units'])
 proof=dict(fixed_rows=[2,3,5],actual_B35_index=26,localized_generator_index=15,prepared_minor_equals_QQunit_times_oldg=str(ratio),all59_source_coefficients_proportional=True,old_pair_g_source_exactly_same_generator=True,old_norm_source='AcceptedPAIR_NORM: clearedNorm(g)=QQ/a/FULLUknownunitfactor*oldh; oldh=G*h1.',whole_product_oldh_G_h1_checked=True,norm_zero_on_closed_Ka13_H4_cover=True,norm_zero_reason='h1(a,t)=0 onKa13 byacceptedprojectionisomorphism, sooldh=G*h1=0. Numeric/a clears areunits. Norm(QQunit*g)=QQunit^4*Norm(g)=0. UunitonwholeH4coveralso makesactualsourcepreconditioningunits.',source_unit_scaling=source_scaling,no_Norm4_recomputation=True,no_alternate_minor_or_Gaussian_or_new_point_prime=True,not_actual_bad_existence_or_empty=True,authorization_not_consumed_by_known_zero_recomputation=True)
 (out/'source-lineage-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',b/'report.json',b/'certificate.json',b/'native-certificate.json',p/'report.json',p/'certificate.json',n/'report.json',n/'native-certificate.json',r/'report.json',r/'certificate.json',cost/'report.json',cost/'certificate.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_NAMED_SATURATION_20261104_v63/report.json']
 report=dict(id=ID,all_passed=True,fixed_minor_is_old_g_exact=True,known_Norm4_zero_no_recomputation=True,no_alternate_job_run=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'source-lineage-proof.json':sha(out/'source-lineage-proof.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
