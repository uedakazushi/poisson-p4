#!/usr/bin/env python3
"""Metadata-only decoded-label and causal producer VALUE budget clarification."""
import hashlib,json,pathlib,re,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
ID='S61_COL53_VERTEX_PRODUCER_BUDGET_CLARIFICATION_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def legacy(x):
 if isinstance(x,bool):return 0
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',x) else 0
 if isinstance(x,list):return sum(map(legacy,x))
 if isinstance(x,dict):return sum(map(legacy,x.values()))
 return 0
def main():
 started=time.monotonic();out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def read(p,h=None):
  s=sha(p)
  if h:assert s==h
  pins[str(p.relative_to(ROOT))]=s;return json.loads(p.read_text())
 gatepath=ROOT/'reports/S61_COL53_VERTEX_SOURCE_GATE_V2_20261105_v64/report.json'
 gate=read(gatepath,'6ad0a98af6dd6539bf2d9cdf16ec8e431ee70df66dbafef540c6aac32c73a1b8')
 gp=gatepath.parent/'graph.json';G=read(gp,gate['evidence_sha256']['graph.json'])['nodes']
 valuepath=ROOT/'reports/S61_COL53_VERTEX_VALUES_V2_20261105_v64/report.json';read(valuepath,'ef9cb1c8175ce59df1cd95274bcc92df719036a61b4d7d41d080237e4b25ea20')
 rp=ROOT/'reports/S61_COL53_VERTEX_SOURCE_VALUES_V2_20261105_v64_REVIEW/report.json';read(rp,'d7e8134bbed2dc1483bef0933d9c892223222583619f72583d564ea2dbb66831')
 old=sum(legacy({k:v for k,v in n.items() if k!='leads'}) for n in G)
 decoded=legacy(G);missing=sum(5 for n in G for v in n['leads'] if v is not None)
 assert old==20856 and decoded==60956 and missing==40100 and decoded-old==missing
 assert gate['owned_upper']==old+len(G)*20+40000==84916
 docs={};base=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
 inputs=[n for n in G if n['op']=='input'];N=max(n['literal_profile']['terms_upper'] for n in inputs);assert N==1762
 for scope in {n['source']['binding']['scope'] for n in inputs if 'binding' in n['source']}:
  path=base/scope;d=read(path,gate['input_sha256'][str(path.relative_to(ROOT))]);docs[scope]=legacy(d)
 assert max(docs.values())==6517
 # Numeric strings counted as THREE slots; all coefficient-map exponent
 # keys counted, although the old genericcounter ignored tuplekeys.
 old_p_and_coeff=14*N+8
 max_decode_overlap=old_p_and_coeff+(8*N+8)
 max_new_coeff_overlap=(8*N+8)+6*N+6*N
 max_binding_doc_overlap=old_p_and_coeff+max(docs.values())
 temporary=max(max_decode_overlap,max_new_coeff_overlap,max_binding_doc_overlap)
 assert temporary==38780
 reserve=temporary+1220;assert reserve==40000
 revised=decoded+20*len(G)+reserve;assert revised==125016<150000
 report={'id':ID,'all_passed':True,'metadata_only':True,'target_arithmetic_rerun':False,
  'old_primary_gate_and_value_bytes_unchanged':True,'independent_reader118477_is_separate_layout':True,
  'old84916_scope':'Pre-JSONsourcegraph counter omitted tuplevertexlabels. It is not the literal decoded VALUE layout. This clarification replaces only that explanatory wholecount, notthehardcap,sourceprooforvalues.',
  'preJSON_counter20856':old,'decoded_counter60956':decoded,'decoded_vertex_labels40100':missing,
  'value_coefficient_array_reserve':'1203*10 QQ coefficients, twointegerlimbs each, 24060slots; aliases into graph are notduplicates.',
  'coefficient_array_upper':20*len(G),'max_input_polynomial_terms':N,'binding_doc_numeric_counts':docs,
  'literal_poly_plus_coefficient_dictionary_per_term_slots':[8,6],
  'causal_overlap_bounds':{'old_poly_and_old_coeff_plus_new_decode':max_decode_overlap,'new_poly_old_coeff_new_coeff_comprehension':max_new_coeff_overlap,'old_poly_old_coeff_plus_newbindingdoc':max_binding_doc_overlap},
  'causal_temporary_max':temporary,'controller_result_metadata_extra':1220,'unchanged_reserved_temporary_slots':40000,
  'producer_decoded_VALUE_whole_positive_upper':revised,'unchanged_hard_limit':150000,
  'scope':'Metadata/ownershipclarification only. Novertexcoefficient,entry,face,rank,parameterpointorprime recomputed. Thetenzero values and theirindependentpeer are unchanged.',
  'input_sha256':pins,'code_sha256':sha(pathlib.Path(__file__)),'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert sha(gatepath)==pins[str(gatepath.relative_to(ROOT))] and sha(valuepath)==pins[str(valuepath.relative_to(ROOT))]
 print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'revised_whole':revised}))
if __name__=='__main__':main()
