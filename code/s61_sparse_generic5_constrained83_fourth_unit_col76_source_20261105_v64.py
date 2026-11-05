#!/usr/bin/env python3
"""SOURCE only: fourth fixedtypedpivot and one smallest cleared-column gate."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
import hashlib,json,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_UNIT_COL76_SOURCE_20261105_v64'
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
def evidence(p,r,name):
 e=p.parent/name;assert sha(e)==r['evidence_sha256'][name];inputs[e]=sha(e);return e,read(e)
p,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');_,v=evidence(p,r,'actual-fixed24-entry-source-and-exponent-plan.json');E={(z['row'],z['column']):z for z in v['actual_entry_sources576']};del v
p,third=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_THIRD_UNIT_SOURCE_20261105_v64','46803a74176f4e299cbbeb561d16c4bfc8d42a23f9677b4adea383e0c6650002');thirdp,t=evidence(p,third,'actual-third-unit-and-full174x21-typed-operator.json')
p,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR3_COL48_COEFFICIENTS_20261105_v64','999fb1eb5b4b7253ef4a123de1cf40238531662f2a104f7719db48ba75a2d1c0');_,api=evidence(p,cr,'actual-residual-col48-and-structural-frontier.json');c48={z['row']:{**z,'path':str((p.parent/z['file']).relative_to(ROOT))}for z in api['actual21_residual_col48_entries']};nu={i:z['terms']for i,z in c48.items()};del api
p,nr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR_COL81_COEFFICIENTS_20261105_v64','0f05a5001b12f4b109801d60ca680416e13fa3243fb06b492af7c88be5b00f1c');_,na=evidence(p,nr,'actual-residual-col81-and-structural-frontier.json');n81={z['row']:z['terms']for z in na['actual22_residual_col81_entries']};del na
p,dr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64','0ca442b4eb828e535bdb527357826c6f4cd7edb48f23735f08742b798ada9371');_,dv=evidence(p,dr,'whole-cached-rows11-44-exact-comparison.json');nd={z['column']:z['difference_terms_count']for z in dv['all24_complete_actual_row_comparisons']};del dv;gc.collect()
def oldpoly(row,col):
 z=E[row,col]
 if not z['terms']:return{'content':'0','shift':[0]*5,'terms':[]},{'kind':'sourcezero','row':row,'column':col}
 ref=z['source'];assert ref['kind']=='unchanged20_column_literal';path=ROOT/ref['file'];h=r['input_sha256'][str(path.relative_to(ROOT))];assert sha(path)==h;inputs[path]=h
 with path.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
 assert hashlib.sha256(blob).hexdigest()==ref['payload_sha256'];a=decode(blob);assert len(a['terms'])==z['terms'];return a,ref
def cpoly(row):
 z=c48[row];path=ROOT/z['path'];assert sha(path)==z['sha256'];inputs[path]=sha(path);blob=path.read_bytes();assert blob[:8]==b'Q5CELL1\n'
 if not z['terms']:assert blob[8]==0;return{'content':'0','shift':[0]*5,'terms':[]},{'path':z['path'],'sha256':z['sha256']}
 assert blob[8]==1;a=decode(blob[9:]);assert len(a['terms'])==z['terms'];return a,{'path':z['path'],'sha256':z['sha256']}
pivot,pivotref=cpoly(31);actual={tuple(pivot['shift'][j]+z[j+1]for j in range(5)):Q(pivot['content'])*int(z[0])for z in pivot['terms']};assert actual=={(2,0,0,1,-1):Q(256,1215),(2,0,0,0,-1):Q(256,6075)}
rows=[i for i in t['selected21_row_ids']if i!=31];cols=[j for j in t['residual21_column_ids']if j!=48];full=[i for i in t['full174_row_ids']if i!=31];top=[i for i in t['top108_row_ids']if i!=31];assert len(rows)==len(cols)==20 and len(full)==173 and len(top)==107
# No expanded M2 values are invented. Countmajorants use explicit original
# source operands, exceptcol81 has accepted actualstage2counts and76unchanged.
def bound2(i,j):
 if j==81:return n81[i]
 if j==76:return E[i,j]['terms']
 return E[i,j]['terms']+E[i,49]['terms']*(E[44,j]['terms']+nd[j])+E[i,51]['terms']*nd[j]
compare=[]
for j in cols:
 cost=sum(2*bound2(i,j)+2*E[i,50]['terms']*bound2(83,j)+nu[i]*bound2(31,j)+nu[i]*E[31,50]['terms']*bound2(83,j)for i in rows);compare.append({'column':j,'literal_tuple_product_majorant':cost})
compare.sort(key=lambda z:(z['literal_tuple_product_majorant'],z['column']));assert compare[0]=={'column':76,'literal_tuple_product_majorant':1068};assert E[31,76]['terms']==0 and E[83,76]['terms']==9 and E[31,50]['terms']==2
bundle={};profiles={}
def leaf(label,a,ref,scale=Q(1),shift=(0,0,0,0,0)):
 n=len(a['terms']);content=Q(a['content'])*scale;S={tuple(a['shift'][j]+z[j+1]+shift[j]for j in range(5))for z in a['terms']};assert len(S)==n;M=max([abs(int(z[0]))for z in a['terms']]+[0]);bundle[label]={'original_literal':a,'source':ref,'QQ_transform_scalar':str(scale),'transform_monomial_exponent5':list(shift)};profiles[label]={'terms':n,'content':content,'support':S,'max_primitive_integer_abs':M};return label
leaf('Q',{'content':'1','shift':[0]*5,'terms':[['5',0,0,0,1,0],['1',0,0,0,0,0]]},{'kind':'explicit5q0plus1'})
a,ref=oldpoly(31,50);leaf('P',a,ref);a,ref=oldpoly(83,76);leaf('R',a,ref)
for i in rows:
 a,ref=oldpoly(i,76);leaf('W'+str(i),a,ref);a,ref=oldpoly(i,50);leaf('V'+str(i),a,ref);a,ref=cpoly(i);leaf('U'+str(i),a,ref,Q(6075,256),(-2,0,0,0,1))
plans=[];alltuples=0;peak=25000+numeric(E)+numeric(t)+numeric(bundle)+sum(z['terms']for z in profiles.values())*6;assert peak<150000
for row in rows:
 recipes=[(['Q','W'+str(row)],Q(1),(0,0,0,0,0)),(['U'+str(row),'P','R'],Q(405,8),(0,0,0,0,-2)),(['Q','V'+str(row),'R'],Q(-405,8),(0,0,0,0,-2))];active=[]
 for labels,preq,preshift in recipes:
  if any(profiles[k]['terms']==0 for k in labels):continue
  content=preq
  for k in labels:content*=profiles[k]['content']
  n=prod(profiles[k]['terms']for k in labels);alltuples+=n;active.append((labels,content,preshift,n))
 D=1
 for _,q,_,_ in active:D=D*q.denominator//gcd(D,q.denominator)
 union=set();M=0;products=[]
 for labels,q,preshift,n in active:
  S={(0,0,0,0,0)};stages=[]
  for k in labels:
   T={tuple(x[j]+y[j]for j in range(5))for x in S for y in profiles[k]['support']};stages.append({'input_support':len(S),'next_factor':k,'next_terms':profiles[k]['terms'],'exact_next_support':len(T)});S=T
  S={tuple(e[j]+preshift[j]for j in range(5))for e in S};union|=S
  assignments=n//max(profiles[k]['terms']for k in labels);integer_bound=assignments*prod(profiles[k]['max_primitive_integer_abs']for k in labels);weight=q.numerator*(D//q.denominator);M+=abs(weight)*integer_bound;products.append({'factor_labels':labels,'source_product_content':str(q),'prefactor_monomial_exponent5':list(preshift),'exact_scalar_integer_weight_after_D':str(weight),'literal_tuple_count':n,'per_output_monomial_assignment_upper':assignments,'coefficient_product_integer_abs_upper':str(integer_bound),'support_stages':stages,'exact_product_support_cardinality':len(S)})
 U=len(union);lower=[min(e[j]for e in union)for j in range(5)]if union else[0]*5;upper=[max(e[j]for e in union)for j in range(5)]if union else[0]*5;eb=sum(max(1,((upper[j]-lower[j]).bit_length()+6)//7)for j in range(5));raw=2048+(D.bit_length()+7)//8+U*((M.bit_length()+7)//8+eb+16);owned=25000+16*(sum(profiles[k]['terms']for k in ['Q','P','R','U'+str(row),'V'+str(row),'W'+str(row)])+U+max([p['exact_product_support_cardinality']for p in products]+[0]))+2048;peak=max(peak,25000+numeric(E)+numeric(t)+numeric(bundle)+numeric(plans)+6*len(union)+512);assert owned<60000 and peak<150000
 plans.append({'row':row,'column':76,'value_representation':'cleared_F_over_common_Q','cleared_formula':'Q*W+(405/(8t0^2))*U*P*R-(405/(8t0^2))*Q*V*R','integer_common_D':str(D),'integer_numerator_abs_upper':str(M),'products':products,'exact_source_support_sorted':[list(e)for e in sorted(union,key=lambda e:(sum(e),e),reverse=True)],'exact_union_support_cardinality':U,'expected_shift_lower5':lower,'upper5':upper,'uncompressed_binary_RAW_upper':raw,'future_owned_extra_including25kmeta':owned,'structural_source_zero':not active,'structurally_unchanged_actual_col76':profiles['U'+str(row)]['terms']==profiles['V'+str(row)]['terms']==0})
assert alltuples==1068
operator={'all_passed':True,'source_only':True,'fourth_pivot_row':31,'fourth_pivot_column':48,'actual_pivot_source':pivotref,'actual_pivot_factorization':'(256/6075)*H^2*(5q0+1)/t0','typed_inverse':'(6075/256)*t0/(H^2*(5q0+1))','added_chart':'Q=5q0+1 nonzero; its complement remainsOPEN','old_chart':'oldbaseunits*P290*detA*detB*H*s0*q0*t0','full173_row_ids':full,'top107_row_ids':top,'selected20_row_ids':rows,'residual20_column_ids':cols,'all87_other_top_RHS':t['old87_other_top_RHS'],'all66_zero_constraint_RHS':t['old66_zero_RHS'],'all173x20_RF_roots':'S4_i_j=S3_i_j-(S3_i48/p4)*S3_31_j; S3 is the frozen exactfull174callback. Everyactualreference is retained.','RHS':'h4_i=h3_i-(S3_i48/p4)*h3_31. Retain all87other_top and66zero constraints.','substitution':'x48=(h3_31-sum_j S3_31_j*xj)/p4; then restorex50,x51,x49, oneoldP290scale, tenconstraints and49units tooriginalnormal42.','full_cokernel':'coker(full174x21) isomorphic coker(full173x20) onaddedQchart; coker(full173)=coker(top107x20) direct_sumK5^66. Selected20square isNOTsamecoker.','no_actual_fixed24_basis_rank24_or_fullinverse_claim':True,'selected_next_source_column':76,'column76_cleared_identity':'Let Ui=(6075/256)*t0/H^2*S3_i48, P=oldA31,50, R=oldA83,76, Vi=oldAi50, Wi=oldAi76. Then Q*S4_i76=Q*Wi+(Ui*P-Q*Vi)*R/p3, p3=8t0^2/405. This usesoldA31,76=0 andcol76unchangedafterfirst2.','new_col76_representation_not_column_operation':'Futurevalues are numeratorF/commonQ. No actualnewQcolumnscale has been applied; never multiplyP290again.','all_other19_new_columns_stay_typed':True,'no_fourth_pivot_polynomial_updates':True}
(out/'actual-fourth-unit-and-full173x20-typed-operator.json').write_bytes(enc(operator));(out/'source-column-comparison-and-cost-model.json').write_bytes(enc({'all_passed':True,'source_only':True,'all20_fixed_candidate_costs':compare,'minimum_in_this_literal_tuple_majorant_model':76,'not_a_field_rank_or_value_search':True,'no_claim_minimal_after_uncomputed_cancellations':True}));(out/'col76-fixed-source-operands.json').write_bytes(enc({'all_passed':True,'parameter_order':['H','L','s0','q0','t0'],'source_operands':bundle}));(out/'all20-cleared-col76-support-height-plans.json').write_bytes(enc({'all_passed':True,'source_only':True,'common_denominator_polynomial':'5q0+1','all20_source_plans':plans,'total_expanded_primitive_tuple_majorant':alltuples,'no_polynomial_values_computed':True}))
base=sum(p.stat().st_size for p in inputs);sourcebytes=sum(p.stat().st_size for p in out.iterdir());future=base+sourcebytes+sum(p['uncompressed_binary_RAW_upper']for p in plans)+2*1024*1024;ownfuture=max(p['future_owned_extra_including25kmeta']for p in plans);assert future<32*1024*1024 and ownfuture<60000
report={'id':ID,'all_passed':True,'source_only':True,'fourth_actual_typed_pivot_verified':True,'added_chart_Q5q0plus1':True,'full_residual_shape':[173,20],'selected_subblock_shape':[20,20],'next_source_column':76,'source_tuple_majorant':alltuples,'future_coeff_owned_upper':ownfuture,'future_coeff_RAW_upper':future,'source_owned_peak':peak,'no_new_polynomial_coefficients_or_value_worker':True,'future_immutable_implementation_gate_REQUIRED':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
