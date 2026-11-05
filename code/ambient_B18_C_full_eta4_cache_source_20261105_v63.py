#!/usr/bin/env python3
"""Source-only full-C four free-row ell^3 readbacks; no requested root is evaluated."""
import hashlib,json,math,re,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_FULL_ETA4_CACHE_SOURCE_20261105_v63'
PATHS={
 'theorem':'reports/AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63/report.json',
 'theorem_peer':'reports/AMBIENT_B18_C_FULL_CUBIC_LINE_SOURCE_20261105_v63_REVIEW/report.json',
 'original':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json',
 'old_compiled':'reports/AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63/compiled-eight-even-roots-and-provided-joint105-plane3-cache.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ell_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
}
ROWS=[19,28,29,72]
MAX_BOUND_H=0

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):return json.dumps(x,separators=(',',':'),sort_keys=True)
def count(x):
 if isinstance(x,(int,float,Q)):return 1
 if isinstance(x,str):return int(bool(re.fullmatch(r'-?\d+(?:/\d+)?',x)))
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,dict):return sum(map(count,x.values()))
 return 0
def ceilabs(q):return(abs(q.numerator)+q.denominator-1)//q.denominator
def bh(a):
 if isinstance(a[0],Q):return max(abs(a[0].numerator).bit_length(),a[0].denominator.bit_length())
 return 0 if not a[0] else a[0].bit_length()+a[1].bit_length()
def track(a):
 global MAX_BOUND_H
 MAX_BOUND_H=max(MAX_BOUND_H,bh(a));assert MAX_BOUND_H<=2048
 return a
def plus(a,b):
 if isinstance(a[0],Q)and isinstance(b[0],Q):return track((a[0]+b[0],1))
 assert not isinstance(a[0],Q)and not isinstance(b[0],Q)
 return track((a[0]+b[0],math.lcm(a[1],b[1])))
def times(a,b):
 if isinstance(a[0],Q):
  if isinstance(b[0],Q):return track((a[0]*b[0],1))
  if not a[0]or not b[0]:return track((0,1))
  return track((ceilabs(a[0])*b[0],a[0].denominator*b[1]))
 if isinstance(b[0],Q):return times(b,a)
 if not a[0]or not b[0]:return track((0,1))
 return track((a[0]*b[0],a[1]*b[1]))
def zero():return[(Q(0),1)]+[(0,1)]*3
def conv(a,b):
 v=zero();v[0]=times(a[0],b[0])
 for n in range(1,4):
  for k in range(n+1):v[n]=plus(v[n],times(a[k],b[n-k]))
 return v

def main():
 started=time.monotonic();signal.alarm(60)
 pins={p:sha(ROOT/p)for p in PATHS.values()}
 th=json.loads((ROOT/PATHS['theorem']).read_text());tp=json.loads((ROOT/PATHS['theorem_peer']).read_text())
 assert th['all_passed']and tp['all_passed']and pins[PATHS['theorem']]in tp['input_sha256'].values()
 # Phase1 holds two immutable DAG encodings only; no root coefficients.
 original=json.loads((ROOT/PATHS['original']).read_text());g=json.loads((ROOT/PATHS['old_compiled']).read_text())
 omap={r[0]:r for r in original['operator_records_originalID_op_args']}
 records=g['operator_records_originalID_op_args'];meta=g['root_pairs_physical4']
 assert[r['physical_row']for r in meta]==ROWS
 roots=[r[n]for r in meta for n in['aA_node','aB_node']]
 assert roots==[20117,40895,20132,40910,20135,40913,20159,40937]
 assert all(omap[r[0]]==r for r in records)
 wanted=set()
 def visit(i):
  if i in wanted:return
  wanted.add(i);r=omap[i]
  if r[1]not in['c','input']:
   for a in r[2:]:visit(a)
 for i in roots:visit(i)
 assert sorted(wanted)==[r[0]for r in records]
 original_inputs=original['provided105_originalID_pure_ell_coefficients']
 original_equations=original['equation105_order']
 # A map/set controller counts its integer keys, as well as all parsed
 # source arrays. Lists below borrow immutable records; they do not clone QQ.
 phase1=count(original)+count(g)+len(omap)+len(wanted)+len(roots)+2000
 assert phase1<=140000
 del original,omap,wanted
 j=json.loads((ROOT/PATHS['joint']).read_text());jp=json.loads((ROOT/PATHS['joint_peer']).read_text())
 ell=json.loads((ROOT/PATHS['ell']).read_text());ep=json.loads((ROOT/PATHS['ell_peer']).read_text())
 assert jp['all_passed']and ep['all_passed']
 assert j['all105_equations_zero_through_degree3']and j['all105_equation_original_node_order']==original_equations
 assert j['monomials']==[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2],[3,0],[2,1],[1,2],[0,3]]
 assert j['all105_unknown_jets']==g['provided105_originalID_r_s_coefficients10']
 old=dict(original_inputs);provided=[];bounds={};overlaps=0
 for i,v in j['all105_unknown_jets']:
  a=[Q(v[k])/5**n for n,k in enumerate([0,2,5,9])]
  assert a[:3]==list(map(Q,old[i]));overlaps+=3
  provided.append([i,list(map(str,a))])
  bounds[i]=[(a[0],1)]+[track((ceilabs(q),q.denominator))for q in a[1:]]
 assert overlaps==315
 expected_all=dict(ell['all164_actual_pure_ell_root_values3'])
 expected=[[i,expected_all[i]]for i in roots]
 phase2=count(g)+count(j)+count(ell)+count(jp)+count(ep)+count(th)+count(tp)+count(provided)+count(bounds)+15000
 assert phase2<=140000
 del j,ell,old,bounds
 provided_bounds={i:[(Q(v[0]),1)]+[track((ceilabs(Q(q)),Q(q).denominator))for q in v[1:]]for i,v in provided}
 uses=Counter()
 for r in records:
  if r[1]not in['input','c']:uses.update(r[2:])
 for i in roots:uses[i]+=1
 cache={};peak=0;inverses=[];mh=hashlib.sha256();constants=[]
 for r in records:
  i,op=r[:2]
  if op=='input':v=provided_bounds[i]
  elif op=='c':v=zero();v[0]=track((Q(r[2],r[3]),1));constants.append([i,r[2],r[3]])
  elif op=='a':v=[plus(a,b)for a,b in zip(cache[r[2]],cache[r[3]])]
  elif op=='n':v=cache[r[2]].copy();v[0]=track((-v[0][0],1))
  elif op=='m':v=conv(cache[r[2]],cache[r[3]])
  elif op=='i':
   x=cache[r[2]];assert x[0][0]==1
   v=zero();v[0]=track((Q(1),1))
   for n in range(1,4):
    for k in range(1,n+1):v[n]=plus(v[n],times(x[k],v[n-k]))
   # Positive absolute-majorant convolution also bounds ALL summands and
   # partial sums in BOTH future exact inverse-product diagnostics.
   conv(x,v);conv(v,x);inverses.append([i,r[2],'1'])
  else:raise AssertionError(op)
  mh.update((enc([i,[[str(a),str(b)]for a,b in v]])+'\n').encode())
  cache[i]=v;peak=max(peak,len(cache))
  if op not in['input','c']:
   for a in r[2:]:
    uses[a]-=1
    if not uses[a]:del cache[a]
  if not uses[i]:del cache[i]
 # Diagnostic delta subtraction is bounded as addition of absolute values.
 deltas=[]
 for r in meta:
  a,b=cache[r['aA_node']],cache[r['aB_node']]
  delta=[track((a[0][0]-b[0][0],1))]+[plus(u,v)for u,v in zip(a[1:],b[1:])]
  deltas.append([r['physical_row'],[[str(u),str(d)]for u,d in delta]])
 assert len(records)==7945 and peak==184
 op_bound=2*MAX_BOUND_H+1;assert op_bound<=2048
 future_budget={
  'compiled_original_records_and_source_labels':count(records)+count(meta)+count(provided)+count(expected)+count(original_inputs)+count(inverses),
  'dependency_counter_integer_values_and_keys':2*len(records),
  'converted105_four_coeffs_and_keys':525,
  'lastuse_four_QQ_buffers':4*peak+16,
  'numeric_controls_profile_peer_and_output_copy_reserve':18000,
 }
 future=sum(future_budget.values());assert future<=75000
 source_phase3=count(g)+count(provided)+count(provided_bounds)+count(expected)+count(cache)+2*len(records)+12000
 source_peak=max(phase1,phase2,source_phase3);assert source_peak<150000
 compiled={
  'parameter':'QQ[ell]/(ell^4), r=0 and s=ell/5 in SAME actual common105 source frame',
  'coefficient_order':['1','ell','ell^2','ell^3'],
  'operator_records_originalID_op_args':records,
  'provided105_originalID_pure_ell4':provided,
  'provided105_originalID_pure_ell3_prior':original_inputs,
  'root_pairs_physical4':meta,
  'expected8_roots_prior_pure_ell3':expected,
  'unit_inverse_centers':inverses,
  'requested_root_coefficients_not_evaluated':True,
 }
 profile={
  'scope':'NEW full-C free eta_j=delta_(ell^3), NOT the old alpha_j=delta_(r*ell^2). Four rows19/28/29/72 in one future job.',
  'source_substitution':'r=0,s=ell/5; input degree k is old(0,k)/5^k. This is a ring homomorphism through degree3.',
  'source_uniqueness':'Same provided105 original IDs/equations/targets, reviewed joint105 fixed unit Jacobian; all315 lower pure-ell input coefficients agree. This binds source branches/gauge, not a new solve.',
  'all315_prior_input_overlap':True,'all24_prior_root_overlaps_required':True,
  'literal_all_n_d_constants':constants,
  'node_count':len(records),'last_use_peak_nodes':peak,'minimal_original_DAG_closure_verified':True,
  'majorant_semantics':'Positive coefficient |q|<=U and den(q)|D. All sum/product and inverse recurrence partials tracked. Both inverse-product diagnostic convolutions and delta subtraction tracked as well. Origin constants exact; no positive root coefficient computed.',
  'every_node_majorant_sha256':mh.hexdigest(),
  'max_stored_or_diagnostic_partial_height_majorant':MAX_BOUND_H,
  'max_operational_before_add_mul_bound':op_bound,
  'root_delta_final_majorants':deltas,
  'source_owned_phase_bounds':[phase1,phase2,source_phase3],
  'source_owned_numeric_peak_bound':source_peak,
  'future_owned_numeric_components':future_budget,'future_owned_numeric_upper':future,
  'future_static_numeric_cap':75000,'future_hard_numeric':100000,
  'global_hard_numeric':150000,'hard_QQ_bits':2048,'hard_seconds':60,'hard_memory_bytes':2147483648,
  'future_RAW_input_output_cap_bytes':1000000,'global_RAW_cap_bytes':33554432,
  'future_mandatory_checks':['all315 known105 input overlaps','all24 prior root coefficients','both unit inverse products in QQ[ell]/ell^4','all four delta coefficients of degree0..2 zero','every numerator and denominator parsed losslessly','precheck before every QQ add/mul/inverse diagnostic/delta','same original IDs and ordered physical rows'],
  'positive_scope':'One nonzero eta_j invokes the complete full-C cubic source theorem and gives actual full HF3=5. This does NOT determine the remaining DVR exponent or original H0/nilcube.',
  'allzero_scope':'All four eta_j=0 exhaust the complete extra cubic operator and give actual full HF3=6. This does NOT determine the remaining DVR exponent or original H0/nilcube.',
  'no_IFT_gauge_physical_solve_or_new_jet':True,
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 cp=out/'compiled-eight-free-even-roots-and-pure-ell4-inputs.json';cp.write_text(enc(compiled)+'\n')
 pp=out/'full-eta4-source-frame-and-whole-operational-profile.json';pp.write_text(json.dumps(profile,indent=2)+'\n')
 assert count(compiled)*2+count(profile)+12000<150000
 report={'id':STEM,'all_passed':True,'status':'SOURCE_ONLY_NEW_FULL_C_ETA4_ELL3_NO_READBACK_VALUES',
  'requested_values_evaluated':False,'input_sha256':pins,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'evidence_sha256':{p.name:sha(p)for p in[cp,pp]},
  'seconds':time.monotonic()-started}
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));signal.alarm(0)

if __name__=='__main__':main()
