#!/usr/bin/env python3
"""Compile only immutable scalar metadata and supplied pure-ell source2 inputs."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63'
PATHS={
 'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
 'gate164':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_SOURCE_PREFLIGHT_20261104_v63/actual-Artin20-source-and-cost-preflight.json',
 'Q20':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json',
 'Q20peer':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW/report.json',
 'Q20code':'code/ambient_B18_Cprime_quadratic20_values_v3_20261105_v63.py',
 'old33':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'C2peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def field(text,key):
 token=json.dumps(key);assert text.count(token)==1;i=text.index(':',text.index(token))+1
 while text[i].isspace():i+=1
 v,j=json.JSONDecoder().raw_decode(text,i);return v,[i,j]
def count(v):
 if isinstance(v,(int,Q)):return 1
 if isinstance(v,str):
  try:Q(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 if isinstance(v,dict):return sum(map(count,v.values()))
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 return 0
def main():
 t=time.monotonic();signal.alarm(60);out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);pins={p:sha(ROOT/p)for p in PATHS.values()}
 gate=json.loads((ROOT/PATHS['gate164']).read_text());qtext=(ROOT/PATHS['Q20']).read_text();mon,_=field(qtext,'base27_monomials');values,vint=field(qtext,'all105_cached_composite_values');nu,_=field(qtext,'nu_A_B27')
 assert mon.index([])==0;ind=[mon.index(m)for m in [[],[8],[8,8]]];phi=[[i,[v[j]for j in ind]]for i,v in values];del values,qtext
 assert {i for i,v in phi}==set(gate['unknown105_order'])and len(phi)==105
 order=gate['needed_closure_order'];roots=gate['requested_root_nodes'];assert len(order)==18584 and len(roots)==164 and order==sorted(order)
 for k in ['Q20peer','C2peer']:assert json.loads((ROOT/PATHS[k]).read_text())['all_passed']
 oldtext=(ROOT/PATHS['old33']).read_text();oldroots,rint=field(oldtext,'all_requested164_root_values');del oldtext
 b33,bint=field((ROOT/PATHS['raw33']).read_text(),'actual33_basis');i0=b33.index([0,0,0]);iell=b33.index([0,0,1]);expected=[]
 for i,v in oldroots:
  z=dict(v);expected.append([i,[z.get(i0,'0'),z.get(iell,'0')]])
 assert {i for i,v in expected}==set(roots);del oldroots,b33
 text=(ROOT/PATHS['source']).read_text();token='"DAG_scalar_nodes"';assert text.count(token)==1;pos=text.index(':',text.index(token))+1
 while text[pos].isspace():pos+=1
 assert text[pos]=='[';start=pos;pos+=1;decoder=json.JSONDecoder();records=[];want=0;seen=set();unknown={i for i,v in phi};oldid=0
 while want<len(order):
  while text[pos].isspace()or text[pos]==',':pos+=1
  node,pos=decoder.raw_decode(text,pos)
  if oldid==order[want]:
   if oldid in unknown:rec=[oldid,'input']
   else:
    op=node[0];assert op in ['c','a','m','n','i']
    if op!='c':assert all(j in seen for j in node[1:])
    rec=[oldid]+node
   records.append(rec);seen.add(oldid);want+=1
  oldid+=1
 assert set(roots).issubset(seen);del text,seen,unknown,node
 metadata=count(records);assert metadata<=3*18584
 inverses=[[r[0],r[2],'1']for r in records if r[1]=='i'];assert inverses==gate['unit_inverse_center_checks']
 # Keep only supplied inputs/metadata. NO scalar DAG node is evaluated here.
 compiled={'operator_records_originalID_op_args':records,'provided105_originalID_pure_ell_coefficients':phi,'requested_root_nodes164':roots,'equation105_order':gate['equation105_order'],'CA_target_positions':gate['target_subtraction_positions'],'actual_target_pure_ell':'z7=0,rA=0,sA=ell/5','expected164_constant_and_ell_coefficients_from_C2':expected,'nu_roots_CA_CB':gate['nu_roots_CA_CB'],'expected_nu_pure_ell3_from_Q20':[[v[j]for j in ind]for v in nu],'physical12_coefficient_roots':gate['physical12_coefficient_roots'],'unit_inverse_centers6':inverses,'source_field_custody':{'DAG_stream_interval':[start,pos],'only_needed_nodes_retained':True,'Q20_105_value_field_bytes':vint,'C2_164_root_field_bytes':rint,'raw33_basis_field_bytes':bint},'all105_unknown_inputs_are_cache_values':True,'scalar_readback_values_not_evaluated':True}
 cp=out/'compiled-original164-DAG-and-provided-pure-ell-inputs.json';cp.write_text(json.dumps(compiled,separators=(',',':'))+'\n')
 profile={'source_operator_records':len(records),'operator_numeric_metadata':metadata,'provided_phi_coefficients':315,'provided_phi_labels':105,'all_retained_DAG_values_dense_upper':3*18584,'dense_ID_value_array':'None control slots and implicit indices, no secondnumeric ID map; nodes index byoriginalID','release_before_evaluation':['needed_closure_order18584','full original source JSON text','source gate otherfields','all27coefficient arrays beyond the supplied315'],
 'future_owned_budget':{'compiled_all_metadata_and_supplied_inputs':count(compiled),'all18584_value_triples':3*18584,'output_all164_triples_plus_labels':656,'output_physical12_triples_and_nu':150,'controller_peer_and_serialization_reserve':10000},
 'source_bit_bound':2048,'bit_proof':'Every intermediate is restriction of SAME Q20 B27 guardedscalar DAG toR0/(z7,rA,ell³). Thisisanalgebrahomomorphism; units stayunits, and Q20 sourceworker guardchecked≤2048 onEACHvalue. No largerbitbound/new solve isassumed.',
 'QQ_metadata_hardcap':150000,'future_static_cap':130000,'future_seconds':60,'future_address_bytes':2<<30,'future_input_plus_output_bytes':2<<20,
 'required_readbacks':['all105equations onprovided105inputs exact, CA targetpos56ell/5','all164constant+ell firstcoeffs matchpriorC2','bothnu pureell3 matchQ20','all6 unitinverse centers1 andproducts1','all36 physical12pureellcoefficients preservefirst overlaps'],
 'scope':'Onlyemitexistingpureell secondrootvalues. NoIFT/gauge/higherjet/source63guess; noC3quadraticrank/cubicrank/exponents yet.',
 'pre_arithmetic_delta_scope_correction':'C2proveddelta2_r²=delta2_rell=0, notdelta2_ell². The latterremainsprovided-cache-readbackunknown untilthisONEjob. Nofullplanezero/upper21/cubicprimitivecancellation waspromoted.'}
 assert sum(profile['future_owned_budget'].values())<130000
 gp=out/'pure-ell-readback-custody-height-and-whole-budget.json';gp.write_text(json.dumps(profile,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_GATE_only_existing164_pure_ell2_scalar_readbacks','new_IFT_gauge_higherjet':False,'readback_values_not_evaluated':True,'future_owned_numeric_upper':sum(profile['future_owned_budget'].values()),'future_seconds':60,'source_bit_bound':2048,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{p.name:sha(p)for p in [cp,gp]},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
