#!/usr/bin/env python3
"""SOURCE only: fifthfixedunit, knownQfactors, fixedcheapestclearedcolumn."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
import argparse,json,hashlib,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_UNIT_COL81_SOURCE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def plus(A,B,c=Q(1)):
 R=dict(A)
 for e,q in B.items():
  R[e]=R.get(e,Q(0))+c*q
  if not R[e]:del R[e]
 return R
def times_Q(A):
 R=dict(A)
 for e,q in A.items():
  f=list(e);f[3]+=1;f=tuple(f);R[f]=R.get(f,Q(0))+5*q
 return{e:q for e,q in R.items()if q}
def divide_known_Q(A):
 if not A:return{},{}
 m=min(e[3]for e in A);groups={}
 for e,q in A.items():groups.setdefault(e[:3]+e[4:],{})[e[3]-m]=q
 T={};R={}
 for key,p in groups.items():
  for j in range(max(p),0,-1):
   q=p.get(j,Q(0))/5
   if q:T[key[:3]+(j-1+m,)+key[3:]]=q;p[j]=p.get(j,Q(0))-5*q;p[j-1]=p.get(j-1,Q(0))-q
  if p.get(0,Q(0)):R[key[:3]+(m,)+key[3:]]=p[0]
 assert plus(times_Q(T),R)==A;return T,R

def primitive(A):
 if not A:return{'content':'0','shift':[0]*5,'terms':[]}
 D=1
 for q in A.values():D=D*q.denominator//gcd(D,q.denominator)
 ints={e:(q*D).numerator for e,q in A.items()};G=0
 for z in ints.values():G=gcd(G,abs(z))
 lo=[min(e[j]for e in A)for j in range(5)];return{'content':str(Q(G,D)),'shift':lo,'terms':[[str(ints[e]//G),*[e[j]-lo[j]for j in range(5)]]for e in sorted(A,key=lambda e:(sum(e),e),reverse=True)]}

def main(a):
 assert a.build
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')};parsed=0;decoded_terms=0
 def load(p):
  nonlocal parsed
  z=read(p);parsed+=numeric(z);inputs[p]=sha(p);assert 25000+parsed<150000;return z
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;z=load(p);assert z['all_passed'];return p,z
 def ev(p,r,name):
  f=p.parent/name;assert sha(f)==r['evidence_sha256'][name];return f,load(f)
 p,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');_,v=ev(p,r,'actual-fixed24-entry-source-and-exponent-plan.json');E={(z['row'],z['column']):z for z in v['actual_entry_sources576']};del v
 p,fr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_UNIT_COL76_SOURCE_20261105_v64','349d1606feb479950eb690241bfdb95c49586731465b9f456e6656d20ea467cc');op4,t=ev(p,fr,'actual-fourth-unit-and-full173x20-typed-operator.json');_,oldbundle=ev(p,fr,'col76-fixed-source-operands.json');oldbundle=oldbundle['source_operands']
 p,nr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR_COL81_COEFFICIENTS_20261105_v64','0f05a5001b12f4b109801d60ca680416e13fa3243fb06b492af7c88be5b00f1c');n81p,n81api=ev(p,nr,'actual-residual-col81-and-structural-frontier.json');n81={z['row']:{**z,'path':str((p.parent/z['file']).relative_to(ROOT))}for z in n81api['actual22_residual_col81_entries']};del n81api
 p,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_COEFFICIENTS_20261105_v64','1c6e8c908858d66248a04499b3671fc399760054088ef65b8a265d4b32f7eeba');c76p,c76api=ev(p,cr,'actual-cleared-col76-and-commonQ-representation.json');c76={z['row']:{**z,'path':str((p.parent/z['file']).relative_to(ROOT))}for z in c76api['actual20_col76_numerator_values']};del c76api
 p,dr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64','0ca442b4eb828e535bdb527357826c6f4cd7edb48f23735f08742b798ada9371');_,dv=ev(p,dr,'whole-cached-rows11-44-exact-comparison.json');nd={z['column']:z['difference_terms_count']for z in dv['all24_complete_actual_row_comparisons']};del dv
 p,pr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_COL76_FIFTH_UNIT_SOURCE_20261105_v64','97b0aaa448bcbea32baa675f374d0337aacb6c11918896d44949b473686df5e6');_,proof5=ev(p,pr,'two-exact-identities-and-candidate-existing-chart-unit.json');assert proof5['no_added_divisor_over_current_Qchart']
 def getbin(z):
  nonlocal decoded_terms
  path=ROOT/z['path'];assert sha(path)==z['sha256'];inputs[path]=sha(path);b=path.read_bytes();assert b[:8]==b'Q5CELL1\n';v=decode(b[9:])if b[8]else{'content':'0','shift':[0]*5,'terms':[]};assert len(v['terms'])==z['terms'];decoded_terms+=len(v['terms']);return v
 # Factor ONLY the knownlinearQ in existingcurrentcol76representations.
 C={};factorproofs=[]
 for i,z in c76.items():
  F=poly(getbin(z));T,R=divide_known_Q(F);assert not R;literal=primitive(T);C[i]=literal;factorproofs.append({'row':i,'original_current_numerator':{'path':z['path'],'sha256':z['sha256']},'exact_current_value_quotient':literal,'Q_remultiply_equals_original_F':True,'new_future_column_value_NOT_computed':True})
 current14=poly(C[14]);current35=poly(C[35]);p5dict=plus(current14,current35,Q(6));assert p5dict=={(1,0,0,1,0):Q(10,6075),(1,0,0,0,0):Q(2,6075)}
 # Fixedsource P factor makes Qdivisibility globalallrows, not an observation.
 P=poly(oldbundle['P']['original_literal']);P0,rem=divide_known_Q(P);assert not rem and P0=={(1,0,0,0,0):Q(-1,2025)}
 T31=poly(getbin(n81[31]));T0,remT=divide_known_Q(T31);assert not remT
 G=plus(poly(getbin(n81[14])),poly(getbin(n81[35])),Q(6));Gq,Gr=divide_known_Q(G);assert Gr and len(G)==9
 rows=[i for i in t['selected20_row_ids']if i!=14];cols=[j for j in t['residual20_column_ids']if j!=76];full=[i for i in t['full173_row_ids']if i!=14];top=[i for i in t['top107_row_ids']if i!=14];assert len(rows)==len(cols)==19 and len(full)==172 and len(top)==106
 def n2(i,j):
  if j==81:return n81[i]['terms']
  return E[i,j]['terms']+E[i,49]['terms']*(E[44,j]['terms']+nd[j])+E[i,51]['terms']*nd[j]
 nu={i:len(oldbundle['U'+str(i)]['original_literal']['terms'])for i in rows};nc={i:len(C[i]['terms'])for i in rows};comparison=[]
 for j in cols:
  cost=sum(2*n2(i,j)+nu[i]*n2(31,j)+2*nu[i]*n2(83,j)+2*E[i,50]['terms']*n2(83,j)+nc[i]*(n2(14,j)+n2(35,j))for i in rows);comparison.append({'column':j,'uniform_source_tuple_majorant':cost})
 comparison.sort(key=lambda z:(z['uniform_source_tuple_majorant'],z['column']));assert comparison[0]['column']==81
 bundle={};profiles={}
 def leaf(label,literal,source,scale=Q(1),shift=(0,0,0,0,0)):
  n=len(literal['terms']);content=Q(literal['content'])*scale;S={tuple(literal['shift'][j]+z[j+1]+shift[j]for j in range(5))for z in literal['terms']};assert len(S)==n;M=max([abs(int(z[0]))for z in literal['terms']]+[0]);bundle[label]={'original_literal':literal,'source':source,'QQ_transform_scalar':str(scale),'transform_monomial_exponent5':list(shift)};profiles[label]={'n':n,'q':content,'S':S,'M':M}
 leaf('Q',{'content':'1','shift':[0]*5,'terms':[['5',0,0,0,1,0],['1',0,0,0,0,0]]},{'kind':'explicit5q0plus1'})
 leaf('T',getbin(n81[31]),{'kind':'actual_stage2_col81','row':31,'path':n81[31]['path'],'sha256':n81[31]['sha256']});leaf('R',getbin(n81[83]),{'kind':'actual_stage2_col81','row':83,'path':n81[83]['path'],'sha256':n81[83]['sha256']})
 leaf('G',primitive(G),{'kind':'fixed_source_row_sum','rows':[14,35],'QQ_coefficients':['1','6'],'source_terms':[{'path':n81[i]['path'],'sha256':n81[i]['sha256']}for i in(14,35)],'Q_division_has_nonzero_remainder':True})
 for i in rows:
  leaf('W'+str(i),getbin(n81[i]),{'kind':'actual_stage2_col81','row':i,'path':n81[i]['path'],'sha256':n81[i]['sha256']});leaf('C'+str(i),C[i],{'kind':'known_Q_quotient_existing_col76','row':i,'original_numerator_path':c76[i]['path'],'original_numerator_sha256':c76[i]['sha256'],'Q_remultiply_witness':True})
  for prefix in('U','V'):
   z=oldbundle[prefix+str(i)];leaf(prefix+str(i),z['original_literal'],z['source'],Q(z['QQ_transform_scalar']),tuple(z['transform_monomial_exponent5']))
 plans=[];tuplesall=0;maxbuffers=0
 for i in rows:
  recipes=[(['Q','W'+str(i)],Q(1),(0,0,0,0,0)),(['U'+str(i),'T'],Q(-1),(0,0,0,0,0)),(['Q','U'+str(i),'R'],Q(-1,40),(1,0,0,0,-2)),(['Q','V'+str(i),'R'],Q(-405,8),(0,0,0,0,-2)),(['C'+str(i),'G'],Q(-6075,2),(-1,0,0,0,0))];active=[]
  for labels,pre,shift in recipes:
   if any(profiles[k]['n']==0 for k in labels):continue
   q=pre
   for k in labels:q*=profiles[k]['q']
   n=prod(profiles[k]['n']for k in labels);tuplesall+=n;active.append((labels,q,shift,n))
  D=1
  for _,q,_,_ in active:D=D*q.denominator//gcd(D,q.denominator)
  union=set();M=0;products=[]
  for labels,q,shift,n in active:
   S={(0,0,0,0,0)};stages=[]
   for k in labels:
    B={tuple(a[j]+b[j]for j in range(5))for a in S for b in profiles[k]['S']};stages.append({'previous_support':len(S),'factor':k,'factor_terms':profiles[k]['n'],'exact_next_support':len(B)});maxbuffers=max(maxbuffers,len(S)+len(B)+len(union));S=B
   S={tuple(e[j]+shift[j]for j in range(5))for e in S};union|=S;assign=n//max(profiles[k]['n']for k in labels);h=assign*prod(profiles[k]['M']for k in labels);w=q.numerator*(D//q.denominator);M+=abs(w)*h;products.append({'factor_labels':labels,'source_product_content':str(q),'prefactor_monomial_exponent5':list(shift),'exact_scalar_integer_weight_after_D':str(w),'literal_tuple_count':n,'per_output_assignment_upper':assign,'product_primitive_integer_abs_upper':str(h),'support_stages':stages,'exact_product_support_cardinality':len(S)})
  U=len(union);lower=[min(e[j]for e in union)for j in range(5)]if union else[0]*5;upper=[max(e[j]for e in union)for j in range(5)]if union else[0]*5;eb=sum(max(1,((upper[j]-lower[j]).bit_length()+6)//7)for j in range(5));raw=2048+(D.bit_length()+7)//8+U*((M.bit_length()+7)//8+eb+16);ninput=sum(profiles[k]['n']for k in ['Q','T','R','G','W'+str(i),'C'+str(i),'U'+str(i),'V'+str(i)]);owned=25000+16*(ninput+U+max([p['exact_product_support_cardinality']for p in products]+[0]))+2048;assert owned<60000
  plans.append({'row':i,'column':81,'value_representation':'new_S5_col81=Laurent_numerator_over_commonQ','integer_common_D':str(D),'integer_numerator_abs_upper':str(M),'products':products,'exact_source_support_sorted':[list(e)for e in sorted(union,key=lambda e:(sum(e),e),reverse=True)],'exact_union_support_cardinality':U,'expected_shift_lower5':lower,'upper5':upper,'uncompressed_binary_RAW_upper':raw,'future_owned_extra_including25kmeta':owned,'structural_source_zero':not active})
 refined=sum(2*n81[i]['terms']+nu[i]*n81[31]['terms']+2*nu[i]*n81[83]['terms']+2*E[i,50]['terms']*n81[83]['terms']+nc[i]*len(G)for i in rows);assert tuplesall==refined
 operator={'all_passed':True,'source_only':True,'row_operation':'row14 <- row14+6row35, determinant1 on ALLcolumns/RHS','fifth_pivot_row_label':14,'fifth_pivot_column':76,'actual_pivot_value':'(2H/6075)*(5q0+1)','actual_typed_inverse':'6075/(2H*(5q0+1))','no_new_chart_divisor':True,'chart':'Previousoldbaseunits*P290*detA*detB*H*s0*q0*t0 and Q=5q0+1 nonzero; QcomplementOPEN.','full172_row_ids':full,'top106_row_ids':top,'selected19_row_ids':rows,'residual19_column_ids':cols,'all87_other_top_RHS':t['all87_other_top_RHS'],'all66_zero_constraints':t['all66_zero_constraint_RHS'],'all172x19_source_roots':'S5_i_j=S4_i_j-(S4_i76/p5)*(S4_14_j+6*S4_35_j). For all172rows and19cols. Everyold source callback retained.','RHS':'h5_i=h4_i-(S4_i76/p5)*(h4_14+6*h4_35). No RHSQdivisibility is assumed.','substitution':'x76=((h4_14+6h4_35)-sum_j(S4_14_j+6S4_35_j)*xj)/p5, thenrestorex48,x50,x51,x49, oneP290scale,tenconstraints,and49units tooriginalnormal42.','full_cokernel':'coker(full173x20) isomorphiccoker(full172x19) onexistingchart; full172coker=top106coker direct_sumK5^66. Selected19square isNOTsamecoker.','commonQ_factor_global_proof':'P=oldA31,50=Q*P0 withP0=-H/2025. Thus everycurrentF_i76=Q*Ci, Ci=Wi+(Ui*P0-Vi)*R76/p3 Laurent for ALLrows, not onlyselected. Current20Ci are exactknownQquotient witnesses.','all_current_col76_Qquotient_representations_preserved':True,'global_commonQ_matrix_representation':'For selectedorcallbackrows: Q*S5_i_j=Q*M2_i_j-Ui*M2_31_j+(Q/p3)*(Ui*P0-Vi)*M2_83_j-J5*Ci*(M2_14_j+6*M2_35_j), J5=6075/(2H). RHS transforms remain separate.','next_fixed_column':81,'next_cleared_formula':'Q*S5_i81=Q*W_i-U_i*T-(H/(40t0^2))*Q*U_i*R-(405/(8t0^2))*Q*V_i*R-(6075/(2H))*C_i*G. W=M2_i81,T=M2_3181,R=M2_8381,G=M2_1481+6M2_3581.','commonQ_is_representation_NOT_new_columnoperation':True,'P290onlyonce_no_newscale':True,'other18_newcolumns_stay_typed':True,'no_new_actual_fifth_col81_or_othercolumn_values':True,'no_fixed24_basis_rank19_or_completeinverse_claim':True}
 factor={'all_passed':True,'source_only':True,'current_col76_numerator_Qquotient_witnesses20':factorproofs,'P_equals_Q_times_P0':{'P_literal':oldbundle['P']['original_literal'],'P0_exact':primitive(P0),'Q_remultiply_checked':True},'T31_equals_Q_times_T0':{'source':n81[31],'T0_exact':primitive(T0),'Q_remultiply_checked':True},'G81_Q_not_divisible':{'G_exact':primitive(G),'known_Q_quotient_part':primitive(Gq),'known_Q_remainder_NONZERO':primitive(Gr),'G_equals_Qquotient_plus_remainder_checked':True},'not_RF_GCD':'Only division by explicitlinear5q0+1 plus literalremultiply; no unknownfieldgcd.','no_futurecolumn_values_computed':True}
 (out/'fifth-unit-and-full172x19-typed-operator.json').write_bytes(enc(operator));(out/'existing-commonQ-factor-and-remainder-witnesses.json').write_bytes(enc(factor));(out/'fixed19-column-cost-comparison.json').write_bytes(enc({'all_passed':True,'source_only':True,'all19_uniform_costs':comparison,'minimum_in_this_model':81,'col81_refined_using_literal9termG':refined,'no_fieldrank_or_otherrow_weight_search':True}));(out/'next-col81-fixed-source-operands.json').write_bytes(enc({'all_passed':True,'parameter_order':['H','L','s0','q0','t0'],'source_operands':bundle}));(out/'all19-cleared-col81-support-height-plans.json').write_bytes(enc({'all_passed':True,'source_only':True,'common_denominator_polynomial':'5q0+1','all19_source_plans':plans,'literal_tuple_majorant':tuplesall,'no_new_values':True}))
 # WholeSOURCE bound includes everyparsedJSON, alloutputobjects and redundant
 # decodedliteralbuffers, even though most were released before supportwork.
 sourceupper=25000+parsed+numeric(operator)+numeric(factor)+numeric(bundle)+numeric(plans)+12*decoded_terms+20*maxbuffers+4096;assert sourceupper<150000
 base=sum(p.stat().st_size for p in inputs);sourcebytes=sum(p.stat().st_size for p in out.iterdir());future=base+sourcebytes+sum(p['uncompressed_binary_RAW_upper']for p in plans)+2*1024*1024;futureown=max(p['future_owned_extra_including25kmeta']for p in plans);assert future<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'fifth_typed_unit_existing_chart_verified':True,'full_residual_shape':[172,19],'selected_subblock_shape':[19,19],'next_fixed_column':81,'fixed_model_col81_minimum':True,'knownQfactor_all_existing20_verified':True,'G81_knownQ_remainder_nonzero':True,'source_tuple_majorant':tuplesall,'future_coeff_owned_upper':futureown,'future_coeff_RAW_upper':future,'whole_source_numeric_upper_including_parse':sourceupper,'no_new_fifth_column_values_or_valueworker':True,'future_immutable_implementation_gate_REQUIRED':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
