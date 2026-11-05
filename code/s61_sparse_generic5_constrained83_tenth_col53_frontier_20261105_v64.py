#!/usr/bin/env python3
"""Read-only failed SOURCE boundary and a universal factored alternative."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,re,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_TENTH_COL53_SOURCE_FRONTIER_20261105_v64'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_TENTH_COL53_STAGED_SOURCE_20261105_v64'
CODE='code/s61_sparse_generic5_constrained83_tenth_col53_staged_source_20261105_v64.py'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def main(a):
 assert a.freeze
 start=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);old=ROOT/'reports'/OLD;inputs={}
 def load(p):
  inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 z=load(old/'partial-source-program.json');failure=load(old/'failure.json');pre=load(old/'packed5-and-whole-SOURCE-preflight.json')
 assert failure['forecast']==164259 and failure['limit']==150000 and failure['context']['output']=='row6_tenth_A'and not failure['completed_rows']
 assert len(z['ops'])==len(z['profiles'])==186
 assert sha(ROOT/CODE)=='df29eb616d6a5b52349ea8af0fac0ed902bc12eb3bf0c09cf177a05976453760'
 inputs[CODE]=sha(ROOT/CODE)
 for p,h in z['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]=h
 def read(stem,name):
  p=ROOT/'reports'/stem/name;inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 E={(v['row'],v['column']):v for v in read('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','actual-fixed24-entry-source-and-exponent-plan.json')['actual_entry_sources576']}
 scope=read('S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64','actual-row35-ninth-unit-and-full168x15-operator.json')
 scope10=read('S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_UNIT_SOURCE_V2_20261105_v64','actual-F10-pivot-and-full167x14-typed-operator.json')
 cost=read('S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64','remaining15-literal-source-cost-comparison.json')
 api=read('S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64','actual-cleared-col82-and-commonQ3RE-representation.json')
 vapi=read('S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64','actual15-col56-and-once-shared-row35.json')
 idx=read('S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_UNIT_COL82_SOURCE_20261105_v64','row-stream-input-and-plan-index16.json')
 common=load(old/'common-source-operands.json')['common_operands']
 base=25000+sum(map(numeric,[E,scope,scope10,cost,api,vapi,idx,common]))+4096
 meta=numeric(z['profiles'])+numeric(z['ops']);n=z['profiles']
 names=['Q','P7','E','S','Delta','A52_53','shared_M2_31_sum_a_sum_b_sum_c','shared_M2_83_sum_a_sum_b_sum_c','shared_M2_104_sum_a_sum_b_sum_c','shared_G','shared_B51','shared_P7F9','shared_UF10','row35_N8_sum_m','shared58_ninth_output','row6_ninth_output','Y6']
 resident=sum(n[k]['terms_upper']for k in names);twice_new=failure['forecast']-base-meta-resident-256
 assert(base,meta,resident)==(44728,5291,13314)and twice_new==2*50335
 product_count=twice_new//2;actual_source_before_alias=base+meta+resident+product_count+256
 assert actual_source_before_alias==113924<150000
 # Formal polynomial identity only, with all source factor symbols distinct.
 labels=['cp','c9','c10','bi','b58','aij','a58j','a35j','ai56','a5856','a3556']
 N=len(labels);zero=(0,)*N
 def sym(i):return{tuple(1 if j==i else 0 for j in range(N)):1}
 def add(A,B,c=1):
  C=dict(A)
  for e,v in B.items():
   C[e]=C.get(e,0)+c*v
   if not C[e]:del C[e]
  return C
 def mul(A,B):
  C={}
  for e,c in A.items():
   for f,d in B.items():
    k=tuple(e[j]+f[j]for j in range(N));C[k]=C.get(k,0)+c*d
  return{k:v for k,v in C.items()if v}
 cp,c9,c10,bi,b58,aij,a58j,a35j,ai56,a5856,a3556=map(sym,range(N))
 nij=add(mul(cp,aij),mul(mul(c9,bi),a35j),-1)
 n58j=add(mul(cp,a58j),mul(mul(c9,b58),a35j),-1)
 ni56=add(mul(cp,ai56),mul(mul(c9,bi),a3556),-1)
 n5856=add(mul(cp,a5856),mul(mul(c9,b58),a3556),-1)
 T=mul(c10,n5856)
 left=add(mul(T,nij),mul(mul(c10,ni56),n58j),-1)
 bracket=add(add(mul(T,aij),mul(mul(c10,ni56),a58j),-1),mul(mul(mul(c9,c10),add(mul(b58,ai56),mul(bi,a5856),-1)),a35j))
 right=mul(cp,bracket);assert left==right
 evidence={
  'all_passed':True,'read_only_diagnosis_no_Minkowski_replay':True,
  'failed_SOURCE_program_completed':False,'completed_profiles_and_ops':186,
  'boundary':{'base_numeric':base,'registered_metadata':meta,'prior_live_supports':resident,'new_positive_support_count_reconstructed_from_guard':product_count,'SOURCE_before_full_hash_alias':actual_source_before_alias,'SOURCE_full_alias_forecast':164259,'guard':150000,'full_alias_not_allocated':True,'no_new_matrix_coefficients':True},
  'future_VALUE_warning':{'same_declared_positive_profile_only_product_7numeric_forecast':7*product_count,'actual_coefficient_support_NOT_computed':True,'hash_stream_change_alone_does_NOT_authorize_VALUE':True},
  'formal_source_identity':{'symbols':labels,'all_universal_coefficients_checked':True,'left_terms':len(left),'equation':'N10=P7*F9*[U*F10*N8_ij-C10*Ni56*N8_58j+C9*C10*(B58*N8_i56-Bi*N8_5856)*N8_35j]','substitutions':['cp=P7*F9','c9=62208000000*s0*q0*t0^2','c10=353894400000*H*t0','U*F10=C10*N9_58,56','N9_k,l=cp*N8_k,l-c9*Bk*N8_35,l'],'source_binding':'The accepted old col56 staged SOURCE + all-coefficient VALUE review pins each Ni56 to this SAME N8 callback; internal N8_i56 was not serialized as a new numeric cache.','not_a_GCD_or_numeric_division':True},
  'possible_changed_representation_NOT_adopted_for_VALUE':{'reduced_denominator':'U^2*F10','numerator':'P7 times the displayed bracket','historical_F9_chart_and_complement_must_still_remain':True,'no_F9_irreducibility_or_coprimality_asserted':True,'fresh_support_lifetime_height_RAW_and_ROOT_go_REQUIRED':True},
  'all167x14_87plus66_RHS_oneP290_other13_and_complements_kept':True
 }
 ep=out/'exact-source-boundary-and-universal-factored-alternative.json';ep.write_bytes(enc(evidence))
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:inputs[str(p.relative_to(ROOT))]=sha(p)
 raw=sum((ROOT/p).stat().st_size for p in inputs)+ep.stat().st_size+65536;assert raw<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'old_SOURCE_guard_preserved':True,'not_a_completed_col53_profile':True,'no_Minkowski_or_coefficient_retry':True,'full_alias_forecast':164259,'SOURCE_without_alias':113924,'declared_product_positive_support':50335,'future_VALUE_positive_profile_product_alone_forecast':352345,'universal_formal_factor_identity_checked':True,'future_representation_source_only_no_VALUE_GO':True,'READ_and_output_RAW_upper':raw,'input_sha256':inputs,'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start}
 (out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
