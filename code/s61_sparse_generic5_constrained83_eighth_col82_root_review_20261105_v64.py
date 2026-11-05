#!/usr/bin/env python3
"""Independent source-only review of ONE col82 convolution bounds and custody."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(65536),b''):h.update(x)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def num(x):
 if isinstance(x,bool):return 0
 if isinstance(x,(float,int)):return 1
 if isinstance(x,dict):return sum(num(v)for v in x.values())
 if isinstance(x,list):return sum(map(num,x))
 return 0
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 target=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64';assert not target.exists(),'source review must precede emission'
 gp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_IMPLEMENTATION_GATE_20261105_v64/report.json';assert sha(gp)=='00170ebc42ae16908f3fe3178386d23b15b5adc5697d7f5704a9d37a89d5260b'
 gate=read(gp);ap=gp.parent/'implementation-custody-storage-bit-addendum.json';a=read(ap)
 for name,h in gate['evidence_sha256'].items():assert sha(gp.parent/name)==h
 for name,h in a['actual_READ_input_sha256'].items():assert sha(ROOT/name)==h
 worker=ROOT/'code/s61_sparse_generic5_constrained83_eighth_col82_coefficients_20261105_v64.py';assert sha(worker)==a['worker_sha256']
 assert a['source_custody_counts']=={'raw':64,'explicit':4,'delta':1,'reuse':80,'M':16}
 assert len(a['all165_leaf_custody_checks'])==165 and all(z['literal_source_identity']for z in a['all165_leaf_custody_checks'])
 index=read(ROOT/a['row_index_file']);commonraw=read(ROOT/index['common_operand_file'])['common_operands'];assert len(index['row_records16'])==16
 rows=[];owned=0;tuples=0;total_maxbits=0
 # Coefficient-free Minkowski supports and integer l1 coefficient bounds.
 for ref in index['row_records16']:
  rp=ROOT/ref['row_operand_file'];pp=ROOT/ref['plan_file'];assert sha(rp)==ref['row_operand_sha256'] and sha(pp)==ref['plan_sha256']
  raw=read(rp);plan=read(pp);sources={**commonraw,**raw['row_operands']};profiles={}
  for k,s in sources.items():
   p=s['original_literal'];shift=[p['shift'][j]+s['transform_monomial_exponent5'][j]for j in range(5)]
   exps={tuple(shift[j]+z[j+1]for j in range(5))for z in p['terms']};assert len(exps)==len(p['terms'])
   profiles[k]=(Q(p['content'])*Q(s['QQ_transform_scalar']),exps,max([abs(int(z[0]))for z in p['terms']]+[0]))
  weights=[];D=1
  for p in plan['products']:
   q=Q(p['full_parameter_monomial_prefactor']['QQ'])
   for name in p['factor_labels']:q*=profiles[name][0]
   assert q==Q(p['source_product_content']);D=D*q.denominator//gcd(D,q.denominator);weights.append(q)
  assert D==int(plan['integer_common_D']);M=0;union=set();peak=0;rowtuples=0
  for p,q in zip(plan['products'],weights):
   pp0=[profiles[k]for k in p['factor_labels']];sizes=[len(z[1])for z in pp0];N=prod(sizes);assert N==p['literal_tuple_count']
   # Fix a coefficient exponent. Select all factors except the largest;
   # its unique exponents admit <=one term for each choice of the others.
   primitive_bound=prod(z[2]for z in pp0)*(N//max(sizes));weight=q.numerator*(D//q.denominator)
   assert primitive_bound==int(p['primitive_product_integer_abs_upper']) and weight==int(p['exact_scalar_integer_weight_after_D'])
   M+=abs(weight)*primitive_bound
   support={(0,0,0,0,0)}
   for _,S,_ in pp0:
    nxt={tuple(x+y for x,y in zip(e,f))for e in support for f in S};peak=max(peak,len(support)+len(nxt)+len(union));support=nxt
   shift=p['full_parameter_monomial_prefactor']['exponent5'];support={tuple(e[j]+shift[j]for j in range(5))for e in support};assert len(support)==p['exact_product_support_cardinality'];union|=support;rowtuples+=N
  assert union==set(map(tuple,plan['exact_source_support_sorted'])) and M==int(plan['integer_numerator_abs_upper'])
  n=len(union);assert 40000+28*n+128<=plan['future_owned_extra_including40kmeta']<=106064
  # Source/support checking owns no convolution value arrays. Count every
  # integer/content exponent in the literal profiles and all live metadata.
  live=16000+num(a)+num(gate)+num(index)+num(commonraw)+num(raw)+num(plan)+7*sum(len(z[1])for z in profiles.values())+6*peak+num(rows)
  owned=max(owned,live);assert live<150000
  assert M.bit_length()+2<=a['family_partial_integer_bit_upper']<2048
  total_maxbits=max(total_maxbits,M.bit_length()+2);tuples+=rowtuples
  rows.append(dict(row=ref['row'],support=n,partial_integer_bits=M.bit_length()+2,source_tuples=rowtuples,owned=live))
  del profiles,sources,raw,plan,pp0,support,union
 assert tuples==2728788 and a['coefficient_owned_limit']==106064 and a['scalar_unreduced_temporary_bit_upper']==124
 assert a['actual_READ_input_RAW_bytes']==sum((ROOT/p).stat().st_size for p in a['actual_READ_input_sha256'])
 assert a['fresh_RAW_upper']==a['actual_READ_input_RAW_bytes']+a['binary_output_upper']+a['other_artifact_upper']+a['prepared_GATE_HELP_input_upper']<32*1024*1024
 # Manual formula audit: stage8 entry S8_ij=S7_ij - S7_i62*A52,j/p8,
 # p8=32768 H^3 E/(273375 t^2), S7_i62=M_i/(Q^3(q+1)).
 # Clear Q^3(q+1)E: every old product receives E and the new summand is
 # -(273375/32768)*t^2/H^3*M_i*A52,j, exactly the frozen final product.
 formula=dict(pivot='32768*H^3*E/(273375*t^2)',E='32*s0+15*q0+35',denominator='Q^3*(q0+1)*E',last_prefactor='-273375*t^2/(32768*H^3)',all_old_numerator_terms_prefixed_by_E=True,no_extra_P290_scale=True,all169_rows_87extra_top_RHS_66_constraints_and_other15_columns_retained=True)
 out=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_ROOT_SOURCE_REVIEW_20261105_v64';out.mkdir(exist_ok=False)
 result=dict(id=out.name,all_passed=True,source_only=True,no_coefficient_emitted=True,source_tuple_count=tuples,source_owned_peak=owned,source_partial_integer_bit_upper=total_maxbits,value_owned_bound=106064,fresh_RAW_upper=a['fresh_RAW_upper'],rows=rows,formula_review=formula,input_sha256={str(gp.relative_to(ROOT)):sha(gp),str(ap.relative_to(ROOT)):sha(ap),str(worker.relative_to(ROOT)):sha(worker),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},limitations='One16 col82 source only; not a full inverse, chart-complement exhaustion, or selected-square cokernel equivalence.',seconds=time.monotonic()-t)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:result[k]for k in ['all_passed','source_owned_peak','source_tuple_count','source_partial_integer_bit_upper','seconds']}));print(sha(out/'report.json'))
if __name__=='__main__':main()
