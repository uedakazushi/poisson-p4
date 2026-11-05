"""SOURCE-only existence of one more generic pivot; no updated matrix values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,re,resource,time
ROOT=Path(__file__).resolve().parents[1]
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
ID='S61_COL53_NONZERO_PIVOT_CHART_SOURCE_20261105_v64'
INPUTS={
 'reports/S61_COL53_H_FACES_VALUES_20261105_v64_REVIEW/report.json':'d1708919e26ff8d484855c987d432cc4e5245cab8ae58164c5862bb44a348c26',
 'reports/S61_COL53_H_FACES_VALUES_20261105_v64/report.json':'b8d4309eaddecdd9b362934f7d2f8bde1973487d071ca507b730f426db8fd983',
 'reports/S61_COL53_H_FACES_VALUES_READBACK_20261105_v64/report.json':'c38b8696b7278e5cbf175b15cd0a068f453db6feb031241f677860d5eb2befbc',
 'reports/S61_COL53_H_FACES_SOURCE_CODE_20261105_v64_REVIEW/report.json':'7c6e258796e3a5cd5f4b6a5d612c56c2def02b4ee03fd587fbefa11a7ff13783',
 'reports/S61_COL53_H_FACES_SOURCE_20261105_v64/report.json':'5e66838fb845def468a9a1c883fd08a5836f21fdfd3a88de8f4d62109013d718',
 'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64/report.json':'e2968d3bfdb21d674d4aad2669fbb4eb87f78c6b815ff3a796c5cb8030c123a7',
 'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_UNIT_SOURCE_V2_20261105_v64/report.json':'ad9dc2bc80fa50e1bf0946f5bf515450a04acad56818cb23851120111be8824b',
 'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_ROOT_SOURCE_REVIEW_20261105_v64/report.json':'560a97ddf59fb740ca67807bb17114034f6879d337844bc6ed1dbc73a1e9633d'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(x):return(json.dumps(x,separators=(',',':'))+'\n').encode()
if __name__=='__main__':
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};loaded=[]
 def read(f,h=None):
  p=ROOT/f
  if h:assert sha(p)==h
  pins[f]=sha(p);v=json.loads(p.read_text());loaded.append(v);return v
 docs={}
 for f,h in INPUTS.items():docs[f]=read(f,h);assert docs[f]['all_passed']
 peer=docs[next(iter(INPUTS))]
 assert peer['independent_actual_signed_SOURCE_VALUE']and peer['fixed_entry_nonzero']and peer['all411_minH_coefficients_and_zero_maxH_face_PASS']
 source=docs['reports/S61_COL53_H_FACES_SOURCE_20261105_v64/report.json']
 graphpath='reports/S61_COL53_H_FACES_SOURCE_20261105_v64/compact-face-graph.json'
 graph=read(graphpath,source['evidence_sha256']['compact-face-graph.json'])
 assert graph['root']==1202 and len(graph['nodes'])==1203
 assert graph['nodes'][1202]['name']=='new_1013_sum'
 # This is the ENTIRE unpruned finite original expression, not the411-term face.
 full_bracket={'file':graphpath,'sha256':pins[graphpath],'root_index':1202,'root_name':'new_1013_sum','all_nodes':1203,'use_ALL_original_args_not_Hmode_tied_subgraphs':True,'not_the_Hminus3_face_polynomial':True}
 readback=docs['reports/S61_COL53_H_FACES_VALUES_READBACK_20261105_v64/report.json']
 witness=readback['faces'][0]['witness']
 assert witness['exponent5']==[-3,0,6,-1,10]and Q(witness['coefficient_QQ'])==Q(241815,134217728)
 assert readback['faces'][1]['actual_terms']==0
 sidepath='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
 side=docs[sidepath+'/report.json'];owners=read(sidepath+'/exact-11-row-coefficient-operator-and-callback-owners.json',side['evidence_sha256']['exact-11-row-coefficient-operator-and-callback-owners.json'])
 assert owners['prefactor']=='P7/(U^2*F10)'and owners['stage8_vector_identity_exact']and owners['all167x14_and87plus66_RHS_P290once_historicalF9_F10_complements']
 assert next(v for v in owners['all_function_owners15']if v['row']==6)['N9_i56_actual_binary']['sha256']=='add8d8119182cec56feb05b38587ec42ac0c442f867e6b60a44a158613dd1019'
 # Universal antisymmetric identity, independently expanded in 12 free symbols.
 names=['P7','F9','C9','C10','Aij','Bi','Wj','A58j','Ai56','A5856','W56','B58']
 n=len(names)
 def V(s):return{tuple(1 if j==names.index(s)else 0 for j in range(n)):1}
 def add(a,b,sign=1):
  d=dict(a)
  for e,c in b.items():
   d[e]=d.get(e,0)+sign*c
   if not d[e]:del d[e]
  return d
 def mul(a,b):
  d={}
  for e,c in a.items():
   for f,v in b.items():
    g=tuple(e[j]+f[j]for j in range(n));d[g]=d.get(g,0)+c*v
  return{e:c for e,c in d.items()if c}
 P7,F9,C9,C10,Aij,Bi,Wj,A58j,Ai56,A5856,W56,B58=(V(s)for s in names)
 cp=mul(P7,F9);Ni56=add(mul(cp,Ai56),mul(mul(C9,Bi),W56),-1)
 N5856=add(mul(cp,A5856),mul(mul(C9,B58),W56),-1)
 UF10=mul(C10,N5856)
 N9ij=add(mul(cp,Aij),mul(mul(C9,Bi),Wj),-1)
 N958j=add(mul(cp,A58j),mul(mul(C9,B58),Wj),-1)
 lhs=add(mul(UF10,N9ij),mul(mul(C10,Ni56),N958j),-1)
 bracket=add(add(mul(UF10,Aij),mul(mul(C10,Ni56),A58j),-1),mul(mul(mul(C9,C10),add(mul(B58,Ai56),mul(Bi,A5856),-1)),Wj))
 rhs=mul(cp,bracket);assert lhs==rhs
 # P7 is already the old Q*R unit, not a new chart factor.
 common=read('reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64/common-source-operands.json')
 p7=common['common_operands']['P7']
 assert p7['QQ_transform_scalar']=='1'and p7['transform_monomial_exponent5']==[0]*5
 assert p7['original_literal']=={'content':'1','shift':[0]*5,'terms':[['5',0,0,0,2,0],['6',0,0,0,1,0],['1',0,0,0,0,0]]}
 typedpath='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_TENTH_UNIT_SOURCE_V2_20261105_v64'
 typed=docs[typedpath+'/report.json']
 T=read(typedpath+'/actual-F10-pivot-and-full167x14-typed-operator.json',typed['evidence_sha256']['actual-F10-pivot-and-full167x14-typed-operator.json'])
 assert T['actual_p10']=='F10/(353894400000*H*t0*F9)'and T['P290_applied_once_no_double_scaling']
 assert T['common_representation_ONLY']['U']=='Q^4*R^2*E'
 assert T['common_representation_ONLY']['D10']=='U^2*F9*F10'
 full=T['full167_row_ids'];top=T['top101_row_ids'];sel=T['selected14_row_ids'];cols=T['residual14_column_ids'];extra=T['all87_other_top_RHS'];zero=T['all66_zero_constraints']
 assert list(map(len,[full,top,sel,cols,extra,zero]))==[167,101,14,14,87,66]
 assert 6 in sel and 53 in cols and set(full)==set(top)|set(zero)and not set(top)&set(zero)
 assert set(top)==set(sel)|set(extra)and not set(sel)&set(extra)
 match=re.search(r'previous(\d+) typed units.*?(\d+) constraint units and(\d+) old units',T['inverse_substitution']);assert match
 previous,constraint,old=map(int,match.groups());typed_count=previous+1;total=typed_count+constraint+old
 assert(typed_count,constraint,old,total)==(10,10,49,69)
 assert len(cols)+total==83 and len(full)+total==236
 after={'full_row_ids':[x for x in full if x!=6],'top_row_ids':[x for x in top if x!=6],'selected_row_ids':[x for x in sel if x!=6],'column_ids':[x for x in cols if x!=53],'extra87_RHS':extra,'zero66_consistency_RHS':zero}
 assert list(map(len,[after['full_row_ids'],after['top_row_ids'],after['selected_row_ids'],after['column_ids']]))==[166,100,13,13]
 # No matrix/RHS multiplication or scalar inverse is evaluated here.
 existence={'all_passed':True,'SOURCE_ONLY_no_new_matrix_coefficients_or_RHS_evaluated':True,'current_counts_from_manifest':{'full':[167,14],'top':[101,14],'selected':[14,14],'extras':87,'zero_constraints':66,'already_completed_typed_units':typed_count,'old_constraint_units':constraint,'old_prefix_units':old,'total_original83_units':total},'one_prospective_pivot':{'row':6,'column':53,'typed_ordinal':typed_count+1,'original83_ordinal':total+1,'full_bracket':full_bracket,'nonzero_witness':witness,'actual_entry':'p = P7*F11_full_DAG/(U^2*F10)','inverse_exists_on_new_chart':'U^2*F10/(P7*F11_full_DAG)','P7':'(5*q0+1)*(q0+1)=Q*R','U':'Q^4*R^2*E','F11_full_DAG_NOT_face_polynomial':True},'chart':{'coefficient_field':'QQ(H,L,s0,q0,t0)','parameter_only':True,'old_units_and_complements':T['old_parameter_units'],'retain_F10':T['actual_F10'],'retain_historical_F9':True,'new_open':'D(F11_full_DAG) inside the SAME old chart','new_and_ALL_old_complements_OPEN':True,'nonempty_open_proof':'product of old nonzero base factors and new nonzero Laurent polynomial is nonzero in the parameter-domain','no_absolute_irreducibility_factor_or_coprimality_claim':True},'prospective_shapes_and_ids_ONLY':after,'formal_Schur_existence_not_coefficient_update':{'m_i_j':'M10[i,j] - M10[i,53]*M10[6,j]/p for every remaining full row and column','RHS':'h10[i] - M10[i,53]*h10[6]/p; ALL87 extra +66 consistency RHS retained','inverse_substitution':'x53=(h10[6]-sum_j M10[6,j]*xj)/p, followed by exact10typed+10constraint+49old restoration with oneP290','cokernel':'on this chart, coker(full167x14)=coker(full166x13); full166coker=top100coker direct_sum K5^66; selected13 square NOT samecokernel','no_actual_typed_update_executed':True},'universal_Plucker_identity_terms':len(lhs),'full_matrix_rank_inverse_exhaustion_NOT_proved':True}
 owned=50000+sum(count_numeric(v)for v in loaded)+count_numeric(after)+10000;assert owned<150000
 raw=sum((ROOT/p).stat().st_size for p in pins)+262144;assert raw<32*1024*1024
 note=ROOT/'research/notes'/f'{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
 (out/'actual-full-bracket-unit-chart-and-manifest-counts.json').write_bytes(enc(existence))
 report={'id':ID,'all_passed':True,'SOURCE_ONLY_no_new_matrix_coefficient_or_RHS_update':True,'independent_actual_value_peer_PINNED':True,'one_more_generic_K5_pivot_exists':True,'prospective_original83_ordinal':total+1,'prospective_typed_ordinal':typed_count+1,'prospective_shape':[166,13],'actual_updated_matrix_NOT_computed':True,'all_complements_and_87plus66_RHS_retained':True,'source_numeric_upper':owned,'RAW_upper':raw,'input_sha256':pins,'evidence_sha256':{'actual-full-bracket-unit-chart-and-manifest-counts.json':sha(out/'actual-full-bracket-unit-chart-and-manifest-counts.json')},'seconds':time.monotonic()-start}
 (out/'report.json').write_bytes(enc(report));print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
