#!/usr/bin/env python3
"""Single-root38648 metadata/positive bound for coherent HEAD23; no values."""
import hashlib,json,math,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEAD23_CACHE_SOURCE_20261105_v63'
PATHS={
 'head':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63/full-coherent-head23-two-socle-cache-polynomials-and-low-two-row-source.json',
 'head_report':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_SOURCE_20261105_v63/report.json',
 'original':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json',
 'original_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63_REVIEW/report.json',
 'joint':'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63/compiled-six-even-roots-and-provided-joint105-plane3-cache.json',
 'joint_peer':'reports/AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63_REVIEW/report.json',
 'oldroots':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'oldroots_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ell_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
}
MON=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)];MI={m:i for i,m in enumerate(MON)}
MAXH=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilabs(q):return(abs(q.numerator)+q.denominator-1)//q.denominator
def ht(x):return max(abs(x[0].numerator).bit_length(),x[0].denominator.bit_length())if isinstance(x[0],Q)else (x[0]*x[1]).bit_length()+1
def track(x):
 global MAXH
 MAXH=max(MAXH,ht(x));assert MAXH<=2048
 return x
def bounds(q):return track((ceilabs(q),q.denominator))
def plus(a,b):
 if isinstance(a[0],Q)and isinstance(b[0],Q):return track((a[0]+b[0],1))
 if isinstance(a[0],Q):a=bounds(a[0])
 if isinstance(b[0],Q):b=bounds(b[0])
 return track((a[0]+b[0],math.lcm(a[1],b[1])))
def times(a,b):
 if isinstance(a[0],Q)and isinstance(b[0],Q):return track((a[0]*b[0],1))
 if isinstance(a[0],Q):a=bounds(a[0])
 if isinstance(b[0],Q):b=bounds(b[0])
 if not a[0]or not b[0]:return(0,1)
 return track((a[0]*b[0],a[1]*b[1]))
def zero():return[(Q(0),1)]+[(0,1)]*5
def conv(a,b):
 v=zero()
 for i,x in enumerate(MON):
  for j,y in enumerate(MON):
   k=MI.get((x[0]+y[0],x[1]+y[1]))
   if k is not None:v[k]=plus(v[k],times(a[i],b[j]))
 return v
def main():
 started=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 original=S.selective(ROOT/PATHS['original'],['operator_records_originalID_op_args','requested_root_nodes164'],receipts)
 omap={r[0]:r for r in original['operator_records_originalID_op_args']};wanted=set()
 def visit(i):
  if i in wanted:return
  wanted.add(i);r=omap[i]
  if r[1]not in['c','input']:
   for j in r[2:]:visit(j)
 visit(38648);records=[omap[i]for i in sorted(wanted)]
 assert 38648 in original['requested_root_nodes164']and len(records)==2194
 source_phase1=S.count(original)+len(omap)+len(wanted)+S.count(records)+12000
 assert source_phase1<=150000
 del original,omap,wanted
 joint=S.selective(ROOT/PATHS['joint'],['provided105_originalID_r_s_coefficients10','provided_monomials10','all630_known105_lower_overlap_checked','B6_monomials','literal_source_substitution'],receipts)
 assert joint['provided_monomials10']==[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]]
 assert joint['B6_monomials']==[list(m)for m in MON]and joint['all630_known105_lower_overlap_checked']
 peers={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k.endswith('_peer')};assert all(p['all_passed']for p in peers.values())
 assert json.loads((ROOT/PATHS['head_report']).read_text())['all_passed']
 given=joint['provided105_originalID_r_s_coefficients10'];assert len(given)==105
 provided_bounds={};conversion_receipts=[]
 for i,arr in given:
  q=list(map(Q,arr));v=zero();v[0]=(q[0],1)
  # Giveninputsubstitutiononly: no requested DAG root is evaluated.
  v[1]=plus(bounds(q[1]),bounds(q[2]));v[2]=times((Q(1,5),1),bounds(q[2]));v[3]=times((Q(1,5),1),plus(bounds(q[4]),times((Q(2),1),bounds(q[5]))));v[4]=times((Q(1,25),1),bounds(q[5]));v[5]=times((Q(1,25),1),plus(bounds(q[8]),times((Q(3),1),bounds(q[9]))))
  provided_bounds[i]=v;conversion_receipts.append([i,[[str(a),str(b)]for a,b in v]])
 old=S.selective(ROOT/PATHS['oldroots'],['all_requested164_root_values'],receipts)
 oldroot=dict(dict(old['all_requested164_root_values'])[38648]);del old
 ell=S.selective(ROOT/PATHS['ell'],['all164_actual_pure_ell_root_values3'],receipts)
 ellroot=dict(ell['all164_actual_pure_ell_root_values3'])[38648];del ell
 expected=[oldroot.get(0,'0'),oldroot.get(8,'0'),ellroot[1],oldroot.get(25,'0'),ellroot[2]]
 assert Q(expected[0])==Q(ellroot[0])
 uses=Counter()
 for r in records:
  if r[1]not in['c','input']:uses.update(r[2:])
 uses[38648]+=1;cache={};peak=0;inverses=[];constants=[];majorhash=hashlib.sha256()
 for r in records:
  i,op=r[:2]
  if op=='input':v=provided_bounds[i]
  elif op=='c':v=zero();v[0]=(Q(r[2],r[3]),1);constants.append(r)
  elif op=='a':v=[plus(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
  elif op=='n':v=cache[r[2]].copy();v[0]=(-v[0][0],1)
  elif op=='m':v=conv(cache[r[2]],cache[r[3]])
  elif op=='i':
   x=cache[r[2]];assert x[0][0]==1
   v=zero();v[0]=(Q(1),1)
   for j,mon in enumerate(MON[1:],1):
    for k,other in enumerate(MON[1:],1):
     tail=(mon[0]-other[0],mon[1]-other[1]);t=MI.get(tail)
     if t is not None and t<j:v[j]=plus(v[j],times(x[k],v[t]))
   conv(x,v);conv(v,x);inverses.append([i,r[2],'1'])
  else:raise AssertionError(op)
  cache[i]=v;peak=max(peak,len(cache));majorhash.update((json.dumps([i,[[str(a),str(b)]for a,b in v]],separators=(',',':'))+'\n').encode())
  if op not in['input','c']:
   for j in r[2:]:
    uses[j]-=1
    if not uses[j]:del cache[j]
  if not uses[i]:del cache[i]
 assert peak==120 and inverses==[[21179,21178,'1']]
 root_bound=cache[38648][5];root_H=MAXH;root_pre=2*root_H+1
 head=json.loads((ROOT/PATHS['head']).read_text());known={tuple(k):Q(v)for k,v,loc in head['all_input_values_already_cached']}
 leaf_bounds={k:bounds(q)for k,q in known.items()};leaf_bounds[('NEW_bB63_lambda12',)]=root_bound
 output_polys=head['head1_genuine_high_A_B_UNPAIRED']+head['head1_canonical_socle_coordinates_rg_k2_UNPAIRED']+[p for pair in head['five_head_low_rows_rH_ellH_UNPAIRED']for p in pair]
 post_bounds=[]
 for poly in output_polys:
  v=(0,1)
  for mon,w in poly:
   term=(Q(w),1)
   for label in mon:term=times(term,leaf_bounds[tuple(label)])
   v=plus(v,term)
  post_bounds.append(v)
 # Simple-target witness is low_ellH(h0)*canonical_k2(head1), with
 # the otheroffdiagonal retainedformal (it doesnotenterdeterminant).
 witness_bound=times(post_bounds[5],post_bounds[3]);post_H=MAXH;pre=max(root_pre,2*post_H+1);assert pre<=2048
 numeric_components={'compiled_records_provided105_lowerchecks':S.count(records)+S.count(given)+S.count(expected)+S.count(inverses),
   'two_dependency_counts_and_originalID_index':3*len(records),'lastuse6QQvalues_with_copy':12*(peak+2),
   'head_polynomials_input_values_locators':S.count(head),'outputs_factors_peer_controller_serialization':10000}
 future=sum(numeric_components.values());assert future<=75000
 source_phase2=S.count(records)+S.count(joint)+S.count(provided_bounds)+S.count(conversion_receipts)+S.count(cache)+2*len(records)+S.count(head)+S.count(peers)+12000
 assert max(source_phase1,source_phase2)<150000
 compiled={'root_original_ID':38648,'root_semantics':'physical CB odd coefficient bB63, notdelta63, notnewgB','operator_records_originalID_op_args':records,
  'provided105_originalID_r_s_coefficients10':given,'provided_monomials10':joint['provided_monomials10'],'B6_monomials':[list(m)for m in MON],
  'B6_equations':['r²=0','ell³=0'],'source_substitution':'s=r+ell/5, SAME provided common105 branch; ordinary coefficients',
  'expected_root_lower_B6_order_1_r_ell_rell_ell2':expected,'unit_inverse_centers':inverses,
  'all630_old_input_lower_provenance':'TheENTIRE supplied105×10 table is exactly the independently peered RIGHT2 cache source table; itsfirst6degree≤2entries/630previousoverlapreceiptareinheritedwithouta newsolve.',
  'requested_root_value_not_evaluated':True}
 profile={'scope':'ONE existing-cache root38648 thirdcoefficient, then coherent two HEAD23 highs andfiveheadlowtwo rows. Nootherheadhigherresponse/newphysical/IFT/gauge/jet.',
  'original_records':len(records),'lastuse_peak':peak,'literal_n_d_constants':constants,
  'all_node_positive_majorant_sha256':majorhash.hexdigest(),'root_stored_and_diagnostic_height_majorant':root_H,'root_before_operation_bound':root_pre,
  'root_lambda_majorant_U_D':[str(q)for q in root_bound],'post_all_partial_height_majorant':post_H,'post_before_operation_bound':2*post_H+1,
  'head_output_polynomial_majorants_U_D':[[str(a),str(b)]for a,b in post_bounds],'simpletarget_minor_majorant_U_D':[str(a)for a in witness_bound],
  'whole_before_operation_bound':pre,'future_owned_components':numeric_components,'future_static_owned_upper':future,
  'source_owned_phase_upper':[source_phase1,source_phase2],'static_cap':75000,'hard_owned_numeric':100000,'global_owned_numeric':150000,
  'hard_QQ_bits':2048,'seconds':60,'memory_bytes':2147483648,'fresh_RAW_cap':2000000,'global_RAW_cap':33554432,
  'mandatory_checks':['literalroot38648sameoriginalDAG','all105providedinputbranchesexact','all5rootlowercoeffsoverlap','unitcenter1andbothinverseproducts',
    'strictparseandprecheckBEFOREeachQQoperation','HEAD23lowtwo0,chi0/1/2/3constraints','allpostpolynomialtermsretained','canonicalrg/k2transformationexact'],
  'positive_simpletarget_minor_scope':'Onlyiflow_ellH(h0)*k2(head1)!=0,actualker(M7→G4)RtypeJ2+k. Other17unknownedgecoefficientsnotcomputed. Ifzero,thiswitnessonlynegative,notkerk³.',
  'no_requestedroot_or_headvalue_evaluated':True,'selective_source_custody':receipts}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);cp=out/'compiled-single-bB63-root-and-provided-joint105-cache.json';cp.write_text(json.dumps(compiled,separators=(',',':'))+'\n')
 pp=out/'single-root-head23-coherent-positive-height-and-whole-budget.json';pp.write_text(json.dumps(profile,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_SINGLE_CACHED_ROOT38648_PLUS_HEAD23_POLYNOMIAL_GATE_NO_VALUES',
  'future_owned_numeric_upper':future,'whole_before_operation_bound':pre,'actual_root_or_edge_values_computed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)for p in[cp,pp]},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
