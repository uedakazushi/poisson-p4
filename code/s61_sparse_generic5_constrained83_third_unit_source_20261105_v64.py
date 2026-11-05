#!/usr/bin/env python3
"""Third actual typedunit and SOURCE bounds for ONLY21 smallcol48differences."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_THIRD_UNIT_SOURCE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
def gate(stem,h):
 p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
p,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');ep=p.parent/'actual-fixed24-entry-source-and-exponent-plan.json';assert sha(ep)==r['evidence_sha256'][ep.name];inputs[ep]=sha(ep);v=read(ep);E={(z['row'],z['column']):z for z in v['actual_entry_sources576']}
p,tr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_UNIT2_FULL175_TRANSPORT_20261105_v64','a2eaa1fb659bda0892a20bb9c4483db5d2f48d3781a6ff8d6769dbc7de530647');tp=p.parent/'full175x22-and-all-RHS-consistency-transport.json';assert sha(tp)==tr['evidence_sha256'][tp.name];inputs[tp]=sha(tp);t=read(tp)
p,ur=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR_COL81_SOURCE_FRONTIER_20261105_v64','47d939f1227ebdb6c088ffade5443f17d2399b32387c550bf147659c91bbc558');up=p.parent/'actual-unchanged-columns-and-six-source-units.json';assert sha(up)==ur['evidence_sha256'][up.name];inputs[up]=sha(up);u=read(up);assert u['post_unit2_columns_exactly_unchanged_all175']==[48,50,76]
def get(z):
 if not z['terms']:return{'content':'0','shift':[0]*5,'terms':[]}
 ref=z['source'];assert ref['kind']=='unchanged20_column_literal';path=ROOT/ref['file'];h=r['input_sha256'][str(path.relative_to(ROOT))];assert sha(path)==h;inputs[path]=h
 with path.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
 assert hashlib.sha256(blob).hexdigest()==ref['payload_sha256'];p=decode(blob);assert len(p['terms'])==z['terms'];return p
def mon(p):
 assert len(p['terms'])==1;z=p['terms'][0];return Q(p['content'])*int(z[0]),[p['shift'][j]+z[j+1]for j in range(5)]
p50,e50=mon(get(E[83,50]));p48,e48=mon(get(E[83,48]));ratio=p48/p50;shift=[x-y for x,y in zip(e48,e50)];assert(p50,e50)==(Q(8,405),[0,0,0,0,2]);assert(ratio,shift)==(Q(256,3),[1,0,0,0,-1])
rows=[i for i in t['selected22_actual_row_ids']if i!=83];cols=[j for j in t['residual22_column_ids']if j!=50];full=[i for i in t['full_residual175_row_ids']if i!=83];top=[i for i in t['remaining109_top_row_ids']if i!=83];assert len(rows)==len(cols)==21 and len(full)==174 and len(top)==108 and len(t['all66_original_formal_zero_constraint_row_ids'])==66
plans=[]
for row in rows:
 operands=[]
 for col,sgn in [(48,1),(50,-1)]:
  z=E[row,col];p=get(z)
  if not z['terms']:continue
  q=Q(p['content'])*sgn*(ratio if col==50 else 1);me=[min(p['shift'][j]+z[j+1]for z in p['terms'])+(shift[j]if col==50 else 0)for j in range(5)];ma=[max(p['shift'][j]+z[j+1]for z in p['terms'])+(shift[j]if col==50 else 0)for j in range(5)];mb=max(abs(int(z[0])).bit_length()for z in p['terms']);operands.append({'column':col,'terms':z['terms'],'QQ_content_with_sign':str(q),'primitive_integer_bit_upper':mb,'lower5_after_transform':me,'upper5_after_transform':ma,'source':z['source']})
 D=1
 from math import gcd
 for op in operands:q=Q(op['QQ_content_with_sign']);D=D*q.denominator//gcd(D,q.denominator)
 M=sum(abs(Q(op['QQ_content_with_sign']).numerator)*(D//Q(op['QQ_content_with_sign']).denominator)*((1<<op['primitive_integer_bit_upper'])-1)for op in operands);U=sum(z['terms']for z in operands)
 lower=[min(z['lower5_after_transform'][j]for z in operands)for j in range(5)]if operands else[0]*5;upper=[max(z['upper5_after_transform'][j]for z in operands)for j in range(5)]if operands else[0]*5;eb=sum(max(1,((upper[j]-lower[j]).bit_length()+6)//7)for j in range(5));cap=2048+(D.bit_length()+7)//8+U*((M.bit_length()+7)//8+eb+16)
 plans.append({'row':row,'column':48,'formula':'oldAi48-(256H/(3t0))*oldAi50','operand_refs':operands,'support_union_cardinality_upper':U,'integer_common_denominator':str(D),'integer_numerator_abs_upper':str(M),'actual_exponent_lower5':lower,'actual_exponent_upper5':upper,'binary_RAW_upper':cap,'whole_owned_upper':30000})
operator={'all_passed':True,'source_only':True,'old_full175_operator':str(tp.relative_to(ROOT)),'old_full175_sha256':sha(tp),'third_pivot_row':83,'third_pivot_column':50,'actual_pivot_QQ_coefficient':str(p50),'actual_pivot_exponent5':e50,'actual_inverse_QQ_coefficient':str(1/p50),'actual_inverse_exponent5':[-x for x in e50],'old_t0_unit_no_extra_chart':True,'full174_row_ids':full,'top108_row_ids':top,'selected21_row_ids':rows,'residual21_column_ids':cols,'old87_other_top_RHS':t['other87_top_RHS_consistency_rows'],'old66_zero_RHS':t['all66_original_formal_zero_constraint_row_ids'],'all174x21_typed_roots':'N_i_j=M2_i_j-(oldA_i50/p)*M2_83_j, for every listed174row and21column. M2 is the exact previousfull175callback, not a newlyexpandedmatrix.','RHS':'h3_i=h2_i-(oldA_i50/p)*h2_83; all87other_top+66constraints kept.','substitution':'x50=(h2_83-sum_j M2_83_j*xj)/p, then restoreoldx51,x49 and frozenearlierunknowns.','all_actual_refs_and_single_P290_scale_kept':True,'cokernel':'coker(full175x22) isomorphiccoker(full174x21) by ONEactualunit. Full174coker=top108coker direct_sumK5^66; selected21square isNOTsamecoker.','column48_simplification_ratio_QQ':str(ratio),'column48_simplification_ratio_exponent5':shift,'column48_exact_formula':'N_i48=oldA_i48-(256H/(3t0))*oldA_i50, since bothcolumns48/50 wereunchangedbyfirst2units.','row17_oldAi50_terms':E[17,50]['terms'],'row17_oldAi50_is_zero':E[17,50]['terms']==0,'row17_unit_persistence_unproved_until_difference':True,'old_postunit2_col81_values_are_only_input_not_currentpostunit3_values':True,'other20_current_columns_stay_typed_no_coefficients_updated':True,'no_fixed24_rank24_or_inverse_claim':True,'chart':'Old baseunits*P290*detA*detB*H*s0*q0*t0. Thirdpivotonlyt0alreadyunit; noextra divisor.','no_new_polynomial_values':True}
(out/'actual-third-unit-and-full174x21-typed-operator.json').write_bytes(enc(operator));(out/'col48-21-source-support-bit-storage-plans.json').write_bytes(enc({'all_passed':True,'source_only':True,'all21_exact_difference_plans':plans,'ratio_QQ':str(ratio),'ratio_exponent5':shift,'no_coefficient_updates':True}))
base=sum(p.stat().st_size for p in inputs);files=sum(p.stat().st_size for p in out.iterdir());future=base+files+sum(p['binary_RAW_upper']for p in plans)+2*1024*1024;owned=25000+numeric(v)+numeric(t)+numeric(u)+numeric(plans)+numeric(operator)+1024;assert owned<150000 and future<32*1024*1024
report={'id':ID,'all_passed':True,'source_only':True,'full_residual_shape':[174,21],'selected_subblock_shape':[21,21],'third_actual_t0_unit_verified':True,'row17_Ai50_not_zero':E[17,50]['terms']>0,'future_col48_owned_upper':30000,'future_col48_RAW_upper':future,'source_owned_upper':owned,'no_polynomial_update':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
