#!/usr/bin/env python3
"""Compile five graph-cubic cache readbacks and a rigorous majorant, no readbacks."""
import hashlib
import json
import math
import re
import signal
import time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C3_GRAPH5_CACHED_SOURCE_V2_20261105_v63'
PATHS={
 'compiled164':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
 'Q20':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2roots':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'raw33':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'frame':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
 'frame_peer':'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json',
}
MONS=[(0,0),(1,0),(0,1),(1,1),(0,2),(1,2)]
MIDX={m:i for i,m in enumerate(MONS)}
PAIRS=[[(a,b)for a,u in enumerate(MONS)for b,v in enumerate(MONS)
         if(u[0]+v[0],u[1]+v[1])==m]for m in MONS]


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):return json.dumps(x,separators=(',',':'),sort_keys=True)
def ncount(v):
 if isinstance(v,(int,float,Q)):return 1
 if isinstance(v,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',v)))
 if isinstance(v,(list,tuple)):return sum(map(ncount,v))
 if isinstance(v,dict):return sum(map(ncount,v.values()))
 return 0
def ceil_abs(q):return(abs(q.numerator)+q.denominator-1)//q.denominator
def zero():return[(0,1)]*6
def input_bound(v,oldmon):
 result=zero()
 for(nr,ns),s in zip(oldmon,v):
  q=Q(s)
  for ld in range(ns+1):
   m=(nr+ns-ld,ld)
   if m in MIDX:
    ix=MIDX[m];t=q*math.comb(ns,ld)/5**ld;u,d=result[ix]
    result[ix]=(u+ceil_abs(t),math.lcm(d,t.denominator))
 result[0]=(Q(v[0]),1)
 return result
def bound_term(a,b):
 if isinstance(a[0],Q):
  if not a[0]:return(0,1)
  if isinstance(b[0],Q):return(a[0]*b[0],1)
  return(ceil_abs(a[0])*b[0],a[0].denominator*b[1])
 if isinstance(b[0],Q):return bound_term(b,a)
 if not a[0]or not b[0]:return(0,1)
 return(a[0]*b[0],a[1]*b[1])
def bound_sum(a,b):
 if isinstance(a[0],Q):return(a[0]+b[0],1)
 return(a[0]+b[0],math.lcm(a[1],b[1]))


def main():
 started=time.monotonic();signal.alarm(60)
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 pins={p:sha(ROOT/p)for p in PATHS.values()}
 d={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()}
 assert all(v['all_passed']for k,v in d.items()if k.endswith('_peer'))
 input_numeric=ncount(d);assert input_numeric<100000
 c,j=d['compiled164'],d['joint'];recs={r[0]:r for r in c['operator_records_originalID_op_args']}
 rows=[63,64,68,73,74]
 meta=[r for r in c['physical12_coefficient_roots']if r['physical_row']in rows]
 assert[r['physical_row']for r in meta]==rows
 roots=[r[n]for r in meta for n in['aA_node','aB_node']]
 wanted=set()
 def visit(i):
  if i in wanted:return
  wanted.add(i);r=recs[i]
  if r[1]not in['input','c']:
   for a in r[2:]:visit(a)
 for i in roots:visit(i)
 order=sorted(wanted);records=[recs[i]for i in order]
 uses=Counter()
 for r in records:
  if r[1]not in['input','c']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 qvals=dict(d['Q20']['all105_cached_composite_values']);qm=d['Q20']['base27_monomials']
 qi=[qm.index(m)for m in[[],[7],[8],[7,7],[7,8],[8,8]]]
 phi=dict(j['all105_unknown_jets']);assert set(phi)==set(qvals)==set(i for i,v in c['provided105_originalID_pure_ell_coefficients'])
 known_overlap=[]
 for i,v in j['all105_unknown_jets']:
  a=list(map(Q,v));lower=[a[0],a[1]+a[2],a[2]/5,a[3]+a[4]+a[5],(a[4]+2*a[5])/5,a[5]/25]
  assert lower==[Q(qvals[i][k])for k in qi]
  known_overlap.append(i)
 assert len(known_overlap)==105 and j['all105_equations_zero_through_degree3']
 # Future readback checks: all five lower B6 coefficients of each root.
 # This does not compute the remaining r*ell^2 coefficient.
 oldroots={i:dict(v)for i,v in d['C2roots']['all_requested164_root_values']}
 ellroots=dict(d['ell']['all164_actual_pure_ell_root_values3'])
 b33=d['raw33']['actual33_basis'];inds=[b33.index(m)for m in[[0,0,0],[0,1,0],[0,0,1],[0,1,1]]]
 expected=[[i,[oldroots[i].get(k,'0')for k in inds]+[ellroots[i][2]]]for i in roots]
 firstframe={r['physical_row']:r for r in d['frame']['all12_verified_source_records']}
 c2_byrow={r['physical_row']:r for r in d['C2']['physical12_coefficients']}
 graph_premises=[]
 nonzz2=[b33.index([i+1,1,0])for i in range(7)]+[b33.index([i+1,0,1])for i in range(7)]
 nonzz2 += [b33.index([0,2,0]),b33.index([0,1,1])]
 for r in rows:
  vals=dict(c2_byrow[r]['delta'])
  assert all(Q(x)==0 for x in firstframe[r]['bA_first_z7'])
  assert all(Q(vals.get(i,'0'))==0 for i in nonzz2)
  assert all(Q(vals.get(b33.index(m),'0'))==0 for m in[[0,3,0],[0,2,1]])
  graph_premises.append([r,True,True,True])
 # A majorant, never a positive-order scalar-DAG readback. Each positive
 # component stores (integer ceil absolute bound, denominator divisor).
 # The point component is an already fixed source-origin constant.
 provided={i:input_bound(v,j['monomials'])for i,v in phi.items()}
 cache={};peak=0;max_u=max_d=max_h=max_center=0;hasher=hashlib.sha256();inverses=[]
 for r in records:
  i,op=r[:2]
  if op=='input':v=provided[i]
  elif op=='c':v=zero();v[0]=(Q(r[2],r[3]),1)
  elif op=='a':v=[bound_sum(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
  elif op=='n':v=cache[r[2]].copy();v[0]=(-v[0][0],1)
  elif op=='m':
   v=zero();v[0]=(cache[r[2]][0][0]*cache[r[3]][0][0],1)
   for m in range(1,6):
    for a,b in PAIRS[m]:v[m]=bound_sum(v[m],bound_term(cache[r[2]][a],cache[r[3]][b]))
  elif op=='i':
   x=cache[r[2]];assert x[0][0]==1;inverses.append([i,r[2],'1'])
   v=zero();v[0]=(Q(1),1)
   for m in range(1,6):
    for a,b in PAIRS[m]:
     if a:v[m]=bound_sum(v[m],bound_term(x[a],v[b]))
  else:raise AssertionError(op)
  max_center=max(max_center,abs(v[0][0].numerator).bit_length(),v[0][0].denominator.bit_length())
  for u,den in v[1:]:
   ub,db=u.bit_length(),den.bit_length();max_u=max(max_u,ub);max_d=max(max_d,db);max_h=max(max_h,ub+db)
   assert ub+db<=2048
  hasher.update((enc([i,[[str(u),str(den)]for u,den in v]])+'\n').encode())
  cache[i]=v;peak=max(peak,len(cache))
  if op not in['input','c']:
   for a in r[2:]:
    uses[a]-=1
    if not uses[a]:del cache[a]
  if not uses[i]:del cache[i]
 height=max(max_h,max_center);operation_height=2*height+1
 assert operation_height<2048 and peak==198 and len(records)==9823
 # Count the SOURCE's live coefficient and metadata copies conservatively.
 source_owned=input_numeric+len(order)*3+peak*13+12000
 assert source_owned<150000
 compiled={
  'B6_monomials':[list(m)for m in MONS],
  'B6_equations':['r^2','ell^3'],
  'operator_records_originalID_op_args':records,
  'provided105_originalID_r_s_coefficients10':j['all105_unknown_jets'],
  'provided_monomials10':j['monomials'],
  'literal_source_substitution':'s=r+ell/5; then quotient by r^2 and ell^3',
  'root_pairs_physical5':meta,
  'expected10_lower_B6_coefficients_order_1_r_ell_rell_ell2':expected,
  'all630_known105_lower_overlap_checked':True,
  'both_original_inverse_centers':inverses,
  'readback_coefficients_not_evaluated':True,
 }
 cp=out/'compiled-ten-even-roots-and-provided-joint105-plane3-cache.json'
 cp.write_text(enc(compiled)+'\n')
 profile={
  'source_scope':'Five actual physical F_graph*tauB cubic relations, rows63,64,68,73,74; no whole C3 presentation/rank or exponent claim.',
  'source_graph_premises_row_bAfirst0_delta2nonzz0_delta3_r3_r2ell0':graph_premises,
  'full_cubic_carry_proof':'These five columns have first0; the delta2 eB component has no nonzz terms, while bA2 eH MAY contain mixed z*r and z*ell terms. Clear delta2/bA2 zz terms by z times the ACTUAL first frame, retaining z in every beta2 carry. Clear bA2 mixed H terms by r/ell times RH_i; its eB plane2 is actual delta_i2_plane=0, so this full mixed correction has zero eB plane3. Primitive bA3 eH and its H corrections die modulo actual plane2 eH; the remaining eB cubic is exactly delta3 eB.',
  'least_unknown_functional':'Existing C2 delta3 r^3 and r^2*ell vanish on these five rows. Mod ell^3, the only remaining cubic function is delta3_(r*ell^2). It equals (old delta_(r*s^2)+3 old delta_(s^3))/25 after the SAME source substitution.',
  'unknowns_are_readbacks_not_new_physics':'All105 source coefficients through3 are already supplied and independently reviewed. New job would only compose ten existing finite scalar-DAG roots; no IFT, gauge, point, source solve or new jet.',
  'same105_common_frame_ids_verified':True,
  'all630_known_input_overlap_originalIDs':known_overlap,
  'source_equations_through3_inherited_from_joint_peer':True,
  'B6_projection_is_an_algebra_homomorphism_and_m4_zero':True,
  'majorant_semantics':'For each positive coefficient, |q|<=integer U and denominator(q) divides D. Addition U sums/D lcm; convolution sums product U/D product-lcm. Exact origin constants retain all zero/one factors. Unit inverse coefficients recur by increasing total degree, origin1. No positive-order scalar readback is computed.',
  'every_node_majorant_sha256':hasher.hexdigest(),
  'max_positive_abs_ceil_bits':max_u,'max_denominator_divisor_bits':max_d,
  'max_coefficient_height_bound':height,'max_future_QQ_temporary_precheck_bound':operation_height,
  'node_count':len(records),'provided105_inputs_used':sum(r[1]=='input'for r in records),'last_use_peak_nodes':peak,
  'root10_final_majorants':[[i,[[str(u),str(den)]for u,den in cache[i]]]for i in roots],
  'source_owned_numeric_upper_including_metadata':source_owned,
  'future_owned_budget':{
    'operator_numeric_metadata':ncount(records),
    'provided105_full10_jets_and_labels':ncount(j['all105_unknown_jets']),
    'future_converted105_six_values_and_labels':735,
    'dependency_count_dense_integer_slots':len(records),
    'lastuse_value_buffers':6*peak+12,
    'retained_expected_outputs_and_source_frame_metadata':2000,
    'controller_peer_copy_output_serialization_reserve':20000},
  'future_static_cap':75000,'future_hard_numeric':100000,
  'future_hard_seconds':60,'future_hard_memory_bytes':2147483648,
  'future_RAW_input_plus_output_bytes':1000000,'hard_RAW_bytes':33554432,'hard_QQ_bits':2048,
  'future_mandatory_checks':['all630 provided105 source2 overlaps rechecked','all10 roots const/r/ell/rell andell2 overlap prior C2/pureell cache','both source inverse products exact in B6','all5 graph delta lower five coefficients zero; only r*ell2 can remain','all10 roots source IDs and row ordering preserved','strict QQ precheck before every operation and root-only STOP on failure'],
  'nonzero_consequence':'One nonzero delta_(r*ell2) gives an actual nonzero eB cubic class; multiplication by r gives nonzero eB r^2*ell^2 in the CI(3,3) summand, hence known C3 upper19 only.',
  'zero_consequence':'All five zero excludes only these five graph-cubic sources. It does not settle C3 or unknown quartic/coherent corrections.',
 }
 assert sum(profile['future_owned_budget'].values())<75000
 ep=out/'graph5-source-carry-cache-binding-and-operational-height-profile.json';ep.write_text(json.dumps(profile,indent=2)+'\n')
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_five_graph_cubic_cache_readbacks_NOT_EVALUATED','new_source_gauge_or_jet':False,'readback_values_evaluated':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)for p in[cp,ep]},'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)


if __name__=='__main__':main()
