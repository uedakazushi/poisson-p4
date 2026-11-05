#!/usr/bin/env python3
"""SOURCE-onlysixth typedunit and row-streamed18-column74 bounds, no values."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
import argparse,json,hashlib,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_UNIT_COL74_SOURCE_20261105_v64'
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
def primitive(A):
 if not A:return{'content':'0','shift':[0]*5,'terms':[]}
 D=1
 for q in A.values():D=D*q.denominator//gcd(D,q.denominator)
 Z={e:(q*D).numerator for e,q in A.items()};G=0
 for z in Z.values():G=gcd(G,abs(z))
 lo=[min(e[j]for e in A)for j in range(5)]
 return{'content':str(Q(G,D)),'shift':lo,'terms':[[str(Z[e]//G),*[e[j]-lo[j]for j in range(5)]]for e in sorted(A,key=lambda e:(sum(e),e),reverse=True)]}
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def mono(p):
 assert len(p['terms'])==1;z=p['terms'][0]
 return Q(p['content'])*int(z[0]),tuple(p['shift'][j]+z[j+1]for j in range(5))
def mm(A,B):return A[0]*B[0],tuple(x+y for x,y in zip(A[1],B[1]))
def mi(A):return 1/A[0],tuple(-x for x in A[1])
def neg(A):return-A[0],A[1]

def main(a):
 assert a.build
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'rows').mkdir();inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')};parsed=0;decoded=0
 def load(p):
  nonlocal parsed
  z=read(p);parsed+=numeric(z);inputs[p]=sha(p);assert 25000+parsed<150000;return z
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=load(p);assert r['all_passed'];return p,r
 def ev(p,r,name):
  e=p.parent/name;assert sha(e)==r['evidence_sha256'][name];return e,load(e)
 p,sr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');_,v=ev(p,sr,'actual-fixed24-entry-source-and-exponent-plan.json');E={(z['row'],z['column']):z for z in v['actual_entry_sources576']};del v
 p,fr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_UNIT_COL81_SOURCE_20261105_v64','d668d98214a0f755afe96990eb143f355c90ecaf8df7a9c9a988eeb4509b11fe');op5,t=ev(p,fr,'fifth-unit-and-full172x19-typed-operator.json');_,oldbundle=ev(p,fr,'next-col81-fixed-source-operands.json');oldbundle=oldbundle['source_operands']
 p,vr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIFTH_COL81_COEFFICIENTS_20261105_v64','6dd3081ad8db9e4ef33801391590c5a78b4359c5d8d0957c48e39f6ad495b10a');_,api=ev(p,vr,'actual-cleared-col81-and-commonQ-representation.json');Fi={z['row']:{**z,'path':str((p.parent/z['file']).relative_to(ROOT))}for z in api['actual19_col81_numerator_values']};del api
 p,ur=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_COL81_SIXTH_UNIT_SOURCE_20261105_v64','35f38b9a19991a566524e986657f4ce3c55a2ba1ba3e4082a932b9ad56d7a30e');_,u6=ev(p,ur,'actual-row104-unit-and-fullrow-RHS-alias.json');assert u6['zero_elimination_coefficients']=={'stage3_V104':True,'stage4_U104':True,'stage5_C104':True}
 p,dr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64','0ca442b4eb828e535bdb527357826c6f4cd7edb48f23735f08742b798ada9371');_,dv=ev(p,dr,'whole-cached-rows11-44-exact-comparison.json');deltas={z['column']:z['difference_terms']for z in dv['all24_complete_actual_row_comparisons']};nd={j:len(z)for j,z in deltas.items()};delta74=deltas[74];delta51=deltas[51];del dv,deltas;gc.collect()
 rows=[i for i in t['selected19_row_ids']if i!=104];cols=[j for j in t['residual19_column_ids']if j!=81];full=[i for i in t['full172_row_ids']if i!=104];top=[i for i in t['top106_row_ids']if i!=104];assert len(rows)==len(cols)==18 and len(full)==171 and len(top)==105
 nu={i:len(oldbundle['U'+str(i)]['original_literal']['terms'])for i in rows};nc={i:len(oldbundle['C'+str(i)]['original_literal']['terms'])for i in rows};nv={i:len(oldbundle['V'+str(i)]['original_literal']['terms'])for i in rows}
 def n2(i,j):return E[i,j]['terms']+E[i,49]['terms']*(E[44,j]['terms']+nd[j])+E[i,51]['terms']*nd[j]
 costs=[]
 for j in cols:
  c=sum(3*n2(i,j)+2*nu[i]*n2(31,j)+3*nu[i]*n2(83,j)+3*nv[i]*n2(83,j)+2*nc[i]*(n2(14,j)+n2(35,j))+Fi[i]['terms']*n2(104,j)for i in rows);costs.append({'column':j,'literal_tuple_majorant_clearQ2':c})
 costs.sort(key=lambda z:(z['literal_tuple_majorant_clearQ2'],z['column']));assert costs[0]=={'column':74,'literal_tuple_majorant_clearQ2':101405};j=74
 # Read only authenticatedsmall cells, but charge the ENTIREoldmatrix READ
 # set andhashit bystream once; no former2.9MBbudget isreused.
 checked_files=set()
 def oldpoly(k,c):
  nonlocal decoded
  z=E[k,c]
  if not z['terms']:return{'content':'0','shift':[0]*5,'terms':[]},{'kind':'sourcezero','row':k,'column':c}
  ref=z['source'];assert ref['kind']=='unchanged20_column_literal';path=ROOT/ref['file'];expected=sr['input_sha256'][str(path.relative_to(ROOT))]
  if path not in checked_files:assert sha(path)==expected;inputs[path]=expected;checked_files.add(path)
  with path.open('rb')as f:f.seek(ref['offset']);blob=f.read(ref['bytes'])
  assert hashlib.sha256(blob).hexdigest()==ref['payload_sha256'];z=decode(blob);decoded+=len(z['terms']);return z,ref
 a0,ra=oldpoly(44,49);b0,rb=oldpoly(44,51);A=mono(a0);B=mono(b0);D=mono(primitive({tuple(z[1:]):Q(z[0])for z in delta51}));assert A==(Q(-1,324000),(0,0,1,0,1))and D==(Q(8,675),(1,0,0,1,0))
 inner_pref=[(Q(1),(0,0,0,0,0)),neg(mi(A)),mm(B,mm(mi(A),mi(D))),neg(mi(D))]
 common={};rowbundles={};profiles={}
 def leaf(label,p,src,scale=Q(1),shift=(0,0,0,0,0),iscommon=False,row=None):
  n=len(p['terms']);S={tuple(p['shift'][h]+z[h+1]+shift[h]for h in range(5))for z in p['terms']};assert len(S)==n;data={'original_literal':p,'source':src,'QQ_transform_scalar':str(scale),'transform_monomial_exponent5':list(shift)};profiles[label]={'n':n,'q':Q(p['content'])*scale,'S':S,'M':max([abs(int(z[0]))for z in p['terms']]+[0])}
  if iscommon:common[label]=data
  else:rowbundles.setdefault(row,{})[label]=data
 leaf('Q',{'content':'1','shift':[0]*5,'terms':[['5',0,0,0,1,0],['1',0,0,0,0,0]]},{'kind':'explicitQ5qplus1'},iscommon=True)
 leaf('Q2',{'content':'1','shift':[0]*5,'terms':[['25',0,0,0,2,0],['10',0,0,0,1,0],['1',0,0,0,0,0]]},{'kind':'explicitknownQsquare'},iscommon=True)
 z,ref=oldpoly(44,j);leaf('S',z,ref,iscommon=True);leaf('Delta',primitive({tuple(z[1:]):Q(z[0])for z in delta74}),{'kind':'acceptedliteralrow11minus8row44','column':j,'report':'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64/report.json'},iscommon=True)
 shared=[31,83,14,35,104]
 for k in shared:
  for c in(j,49,51):
   z,ref=oldpoly(k,c);leaf(f'A{k}_{c}',z,ref,iscommon=True)
 for i in rows:
  for c in(j,49,51):
   label=f'A{i}_{c}'
   if label in common:continue
   z,ref=oldpoly(i,c);leaf(label,z,ref,row=i)
  for prefix in('U','V','C'):
   z=oldbundle[prefix+str(i)];leaf(prefix+str(i),z['original_literal'],z['source'],Q(z['QQ_transform_scalar']),tuple(z['transform_monomial_exponent5']),row=i)
  z=Fi[i];path=ROOT/z['path'];assert sha(path)==z['sha256'];inputs[path]=sha(path);blob=path.read_bytes();value=decode(blob[9:])if blob[8]else{'content':'0','shift':[0]*5,'terms':[]};assert len(value['terms'])==z['terms'];decoded+=len(value['terms']);leaf('F'+str(i),value,{'kind':'actualfifthcol81numerator','row':i,'path':z['path'],'sha256':z['sha256'],'value_representation':'S5_i81=F_i/Q'},row=i)
 # Sourceprovenance/sharedutility sizes are partofMETA; VALUEphase MUSTstream
 # one rowbundle/plan andcompile/droprawliterals, notparseallrowbundles.
 profile_n=sum(z['n']for z in profiles.values());post_base=25000+numeric(E)+numeric(t)+numeric(oldbundle)+numeric(common)+sum(numeric(z)for z in rowbundles.values())+6*profile_n+4096;assert post_base<150000
 (out/'common-source-operands.json').write_bytes(enc({'all_passed':True,'parameter_order':['H','L','s0','q0','t0'],'common_operands':common}))
 unitproof={'a_source':ra,'a_monomial':{'QQ':str(A[0]),'exponent5':list(A[1])},'b_source':rb,'b_monomial':{'QQ':str(B[0]),'exponent5':list(B[1])},'d_known_delta51':{'QQ':str(D[0]),'exponent5':list(D[1])},'all_four_inner_prefactors':[{'QQ':str(q),'exponent5':list(e)}for q,e in inner_pref]}
 (out/'old-unit2-monomial-source-prefactors.json').write_bytes(enc(unitproof))
 peak=max(25000+parsed,post_base);alltuples=0;index=[];caps=0;maxowned=0;maxbuffers=0
 for i in rows:
  # Q*H5: sevenouterterms, thenEXACTfour-source expansionofeachM2_kj.
  outer=[(i,['Q2'],(Q(1),(0,0,0,0,0))),(31,['Q','U'+str(i)],(Q(-1),(0,0,0,0,0))),(83,['Q2','U'+str(i)],(Q(-1,40),(1,0,0,0,-2))),(83,['Q2','V'+str(i)],(Q(-405,8),(0,0,0,0,-2))),(14,['Q','C'+str(i)],(Q(-6075,2),(-1,0,0,0,0))),(35,['Q','C'+str(i)],(Q(-18225),(-1,0,0,0,0))),(104,['F'+str(i)],(Q(-765450000),(-1,0,0,-1,0)))]
  active=[]
  for k,olabels,opre in outer:
   inners=[([f'A{k}_{j}'],inner_pref[0]),([f'A{k}_49','S'],inner_pref[1]),([f'A{k}_49','Delta'],inner_pref[2]),([f'A{k}_51','Delta'],inner_pref[3])]
   for ilabels,ipre in inners:
    labels=olabels+ilabels
    if any(profiles[l]['n']==0 for l in labels):continue
    pref=mm(opre,ipre);q=pref[0]
    for l in labels:q*=profiles[l]['q']
    nt=prod(profiles[l]['n']for l in labels);alltuples+=nt;active.append((k,labels,pref,q,nt))
  dc=1
  for _,_,_,q,_ in active:dc=dc*q.denominator//gcd(dc,q.denominator)
  union=set();M=0;products=[]
  for k,labels,pref,q,nt in active:
   # Forecast all source supportsets before each update, including oldS/newT,
   # accumulatorunion, allrowplans alreadyondisk, neverresidentall18plans.
   literal_upper=prod(profiles[l]['n']for l in labels);forecast=post_base+numeric(products)+numeric(index)+20*(len(union)+literal_upper)+4096;assert forecast<150000
   S={(0,0,0,0,0)};stages=[]
   for l in labels:
    T={tuple(x[h]+y[h]for h in range(5))for x in S for y in profiles[l]['S']};maxbuffers=max(maxbuffers,len(S)+len(T)+len(union));stages.append({'factor':l,'previous_support':len(S),'factor_terms':profiles[l]['n'],'exact_next_support':len(T)});S=T
   S={tuple(e[h]+pref[1][h]for h in range(5))for e in S};union|=S;assign=nt//max(profiles[l]['n']for l in labels);h=assign*prod(profiles[l]['M']for l in labels);w=q.numerator*(dc//q.denominator);M+=abs(w)*h;products.append({'source_M2_row':k,'factor_labels':labels,'full_parameter_monomial_prefactor':{'QQ':str(pref[0]),'exponent5':list(pref[1])},'source_product_content':str(q),'exact_scalar_integer_weight_after_D':str(w),'literal_tuple_count':nt,'per_output_assignment_upper':assign,'primitive_product_integer_abs_upper':str(h),'support_stages':stages,'exact_product_support_cardinality':len(S)})
  U=len(union);lower=[min(e[h]for e in union)for h in range(5)]if union else[0]*5;upper=[max(e[h]for e in union)for h in range(5)]if union else[0]*5;eb=sum(max(1,((upper[h]-lower[h]).bit_length()+6)//7)for h in range(5));raw=2048+(dc.bit_length()+7)//8+U*((M.bit_length()+7)//8+eb+16);nshared=sum(profiles[l]['n']for l in common);nrow=sum(profiles[l]['n']for l in rowbundles[i]);owned=25000+16*(nshared+nrow+U+max([p['exact_product_support_cardinality']for p in products]+[0]))+2048;assert owned<150000;maxowned=max(maxowned,owned)
  plan={'row':i,'column':j,'value_representation':'new_S6_i74=Laurent_numerator_over_Q_squared','integer_common_D':str(dc),'integer_numerator_abs_upper':str(M),'products':products,'exact_source_support_sorted':[list(e)for e in sorted(union,key=lambda e:(sum(e),e),reverse=True)],'exact_union_support_cardinality':U,'expected_shift_lower5':lower,'upper5':upper,'uncompressed_binary_RAW_upper':raw,'future_owned_extra_including25kmeta':owned,'future_common_and_row_leaf_terms':nshared+nrow,'no_B104j_Qdivisibility_assumed':True}
  rp=out/'rows'/f'row-{i}-source-operands.json';pp=out/'rows'/f'row-{i}-support-height-plan.json';rp.write_bytes(enc({'all_passed':True,'row':i,'row_operands':rowbundles[i]}));pp.write_bytes(enc(plan));index.append({'row':i,'row_operand_file':str(rp.relative_to(ROOT)),'row_operand_sha256':sha(rp),'plan_file':str(pp.relative_to(ROOT)),'plan_sha256':sha(pp),'source_tuple_majorant':sum(p['literal_tuple_count']for p in products),'binary_RAW_upper':raw,'future_owned_upper':owned});caps+=raw
  peak=max(peak,post_base+numeric(products)+numeric(index)+numeric(plan)+20*maxbuffers+4096);assert peak<150000;del union,products,plan
 assert alltuples==101405
 operator={'all_passed':True,'source_only':True,'sixth_pivot_row':104,'sixth_pivot_column':81,'actual_unit':'H*q0*(5q0+1)/765450000','typed_inverse':'765450000/(H*q0*(5q0+1))','no_added_chart_factor':True,'current_row104_all19_and_arbitrary_RHS_alias':'M2_row104/h2_row104 fromacceptedzeroV/U/C. DoNOTreplacebyrawA104; M2includesfirsttwoactualsourceSchurs.','full171_row_ids':full,'top105_row_ids':top,'selected18_row_ids':rows,'residual18_column_ids':cols,'all87_other_top_RHS':t['all87_other_top_RHS'],'all66_zero_constraints':t['all66_zero_constraints'],'all171x18_RF_roots':'S6_i_j=S5_i_j-(S5_i81/p6)*M2_104_j; ALLoriginalsourcecallbacks andnormalcoordinateRHS retained.','RHS':'h6_i=h5_i-(S5_i81/p6)*h2_104. No QdivisibilityofanyRHSassumed.','substitution':'x81=(h2_104-sum_j M2_104_j*xj)/p6, thenrestorex76(row14+6row35),x48,x50,x51,x49,oneP290scale,tenconstraintsand49units.','full_cokernel':'full172x19 unit-equivalenttofull171x18; full171coker=top105coker direct_sumK5^66. Selected18square isNOTsamecoker.','next_fixed_column':j,'global_clear_Q2_identity':'Q^2*S6_i_j=Q*H5_i_j-(765450000/(H*q0))*F5_i81*M2_104_j; sourceH5commonQ andnewF5acceptednumerators.','all_28_source_products_per_row':'ExpandsevenouterM2coefficientsandEXACTM2_kj=A_kj-A_k49*S/a+A_k49*b*Delta/(a*d)-A_k51*Delta/d. NoM2_104j auxiliaryvalueorQdivisibilityinvented.','knownQfactor_policy':'UsecommonQ^2conservatively. No hiddenRFgcdorunprovedQcancel. TheQcomplement remainsOPEN.','representation_NOT_new_columnoperation':True,'P290onlyonce_all_other17_columns_stay_typed':True,'no_new_actual_sixth_column_values':True,'no_fixed24_basis_rank18_globalinverse_or_exhaustion_claim':True,'chart':'Allpreviousbaseunits*P290*detA*detB*H*s0*q0*t0 andQnonzero; no newchartfactor.','VALUEstreaming_obligation':'ParseONLYcommonleafs+one rowleafbundle+one rowplan; doNOTloadall18supportplans orfulloldmatrix coefficients.'}
 (out/'sixth-unit-and-full171x18-typed-operator.json').write_bytes(enc(operator));(out/'fixed18-source-cost-comparison.json').write_bytes(enc({'all_passed':True,'source_only':True,'all18_uniform_clearQ2_source_costs':costs,'selected_minimum_in_this_model':j,'no_field_rank_or_otherrow_weightsearch':True}));(out/'row-stream-input-and-plan-index18.json').write_bytes(enc({'all_passed':True,'source_only':True,'common_operand_file':str((out/'common-source-operands.json').relative_to(ROOT)),'common_operand_sha256':sha(out/'common-source-operands.json'),'row_records18':index,'all18_source_tuple_majorant':alltuples,'VALUE_must_stream_one_row':True}))
 base=sum(p.stat().st_size for p in inputs);sourcebytes=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());future=base+sourcebytes+caps+2*1024*1024;assert future<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'sixth_typed_unit_existing_chart':True,'full_residual_shape':[171,18],'selected_subblock_shape':[18,18],'next_fixed_column':74,'source_tuple_majorant':alltuples,'common_Q_denominator_power':2,'future_coeff_owned_upper':maxowned,'future_coeff_RAW_upper':future,'whole_source_numeric_upper':peak,'old14MBmatrix_actualREAD_correctly_charged':True,'no_new_sixth_column_values_or_valueworker':True,'future_immutable_implementation_gate_REQUIRED':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
