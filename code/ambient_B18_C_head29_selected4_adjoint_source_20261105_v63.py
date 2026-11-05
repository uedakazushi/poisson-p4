#!/usr/bin/env python3
"""SOURCE only: one origin adjoint for delta63_r²ell², no fourth response."""
import hashlib,json,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_HEAD29_SELECTED4_ADJOINT_SOURCE_20261105_v63'
PATHS={
 'original':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json',
 'original_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63_REVIEW/report.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'availability':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_CACHE_AVAILABILITY_20261105_v63/literal-raw33-plane9-map-and-one-selected-fourth-head29-gate.json',
}
MAXPRE=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def h(q):return max(abs(q.numerator).bit_length(),q.denominator.bit_length())
def before(n):
 global MAXPRE
 MAXPRE=max(MAXPRE,n);assert n<=2048,'SOURCE_QQ_PREOP_STOP'
def add(a,b):before(2*max(h(a),h(b))+1);return a+b
def mul(a,b):before(h(a)+h(b));return a*b
def divide(a,b):assert b;before(h(a)+h(b));return a/b
def put(d,k,q):
 if not q:return
 d[k]=add(d.get(k,Q(0)),q)
 if not d[k]:del d[k]
def sparse_mm(A,B):
 out=[]
 for row in A:
  v={}
  for k,a in row.items():
   for j,b in B[k].items():put(v,j,mul(a,b))
  out.append(v)
 return out
def closure(roots,omap):
 wanted=set();todo=list(roots)
 while todo:
  i=todo.pop()
  if i in wanted:continue
  wanted.add(i);r=omap[i]
  if r[1]not in['c','input']:todo.extend(r[2:])
 return wanted
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 original=S.selective(ROOT/PATHS['original'],['operator_records_originalID_op_args','provided105_originalID_pure_ell_coefficients','equation105_order'],receipts)
 omap={r[0]:r for r in original['operator_records_originalID_op_args']};roots=[20144,40922];need=closure(roots,omap)
 g_records=[omap[i]for i in sorted(need)]
 source_phase0=S.count(original)+len(omap)+len(need)+S.count(g_records)+10000;assert source_phase0<150000
 phase={'phase':'frozen_metadata_before_originAD','input_sha256':pins,'g_roots':roots,'g_closure_records':len(g_records),
  'numeric_source_preflight_upper':source_phase0,'hard_numeric':150000,'hard_QQ_bits':2048,'hard_seconds':60,'hard_memory':2147483648,
  'fourth_response_RHS_or_lower_coefficients_NOT_evaluated':True}
 (out/'phase.json').write_text(json.dumps(phase,indent=2)+'\n')
 # Preserve the old full records by reference only while deriving the
 # eventual rho-supported closure; origin arithmetic uses justg_records.
 joint=S.selective(ROOT/PATHS['joint'],['all105_unknown_original_node_order','all105_equation_original_node_order',
   'all105_unknown_jets','fixed105_Jacobian','fixed105_Jacobian_inverse','CA_F9_targets'],receipts)
 for k in['original_peer','joint_peer']:assert json.loads((ROOT/PATHS[k]).read_text())['all_passed']
 U=joint['all105_unknown_original_node_order'];F=joint['all105_equation_original_node_order'];UI={i:j for j,i in enumerate(U)}
 assert F==original['equation105_order']and len(U)==len(F)==105 and joint['CA_F9_targets']==['0']*7+['r','s']
 given={i:Q(v[0])for i,v in joint['all105_unknown_jets']};assert set(given)==set(U)
 assert given=={i:Q(v[0])for i,v in original['provided105_originalID_pure_ell_coefficients']}
 # Cull the large immutable record table from the arithmeticphase count;
 # only omap retains shared original records, countedexplicitly below.
 points={};unit_origins=[]
 for r in g_records:
  i,op=r[:2]
  if op=='input':v=given[i]
  elif op=='c':before(max(abs(r[2]).bit_length(),abs(r[3]).bit_length()));v=Q(r[2],r[3])
  elif op=='a':v=add(points[r[2]],points[r[3]])
  elif op=='n':v=-points[r[2]]
  elif op=='m':v=mul(points[r[2]],points[r[3]])
  elif op=='i':assert points[r[2]]==1;v=divide(Q(1),points[r[2]]);unit_origins.append([i,r[2],'1'])
  else:raise AssertionError(op)
  points[i]=v
 assert points[20144]-points[40922]==0
 weights={20144:Q(1),40922:Q(-1)};grad={}
 for r in reversed(g_records):
  i,op=r[:2];w=weights.get(i,Q(0))
  if not w:continue
  if op=='input':put(grad,UI[i],w)
  elif op=='c':pass
  elif op=='a':put(weights,r[2],w);put(weights,r[3],w)
  elif op=='n':put(weights,r[2],-w)
  elif op=='m':put(weights,r[2],mul(w,points[r[3]]));put(weights,r[3],mul(w,points[r[2]]))
  elif op=='i':put(weights,r[2],-divide(w,mul(points[r[2]],points[r[2]])))
 J=[{i:Q(q)for i,q in enumerate(row)if Q(q)}for row in joint['fixed105_Jacobian']]
 Inv=[{i:Q(q)for i,q in enumerate(row)if Q(q)}for row in joint['fixed105_Jacobian_inverse']]
 identity=[{i:Q(1)}for i in range(105)];assert sparse_mm(J,Inv)==identity==sparse_mm(Inv,J)
 rho=sparse_mm([grad],Inv)[0];assert sparse_mm([rho],J)[0]==grad
 needed_roots=roots+[F[i]for i in sorted(rho)];future_need=closure(needed_roots,omap)
 future_records=[omap[i]for i in sorted(future_need)]
 uses=Counter(needed_roots)
 for r in future_records:
  if r[1]not in['c','input']:uses.update(r[2:])
 live=set();peak=0
 for r in future_records:
  i,op=r[:2];live.add(i);peak=max(peak,len(live))
  # Root RHScoefficient canbestreamedatanode'screation; no extra nine-vector
  # needsremainaliveonlybecauseitcontributestothefinalrhosum.
  if i in needed_roots:
   uses[i]-=needed_roots.count(i)
  if op not in['c','input']:
   for j in r[2:]:
    uses[j]-=1
    if not uses[j]:live.remove(j)
  if not uses[i]:live.remove(i)
 arithmetic_phase=S.count(original)+len(omap)+S.count(joint)+S.count(points)+S.count(weights)+S.count(grad)+S.count(rho)+S.count(J)+S.count(Inv)+S.count(g_records)+10000
 # These are source-custody aliases, not QQ copies. Full immutable records
 # were already counted once inoriginal; do notcharge their shared records
 # asecondtimeing_records (only listindex metadatabelow).
 arithmetic_phase-=S.count(g_records);arithmetic_phase+=len(g_records)
 assert arithmetic_phase<150000,('SOURCE_ARITHMETIC_OWNERSHIP_STOP',arithmetic_phase)
 data={
  'scope':'ONE selectedIFT4 sourceadjoint atfrozenorigin ONLY. No R4 RHS,lowerdegree4functionalevaluation,truePhi4 ordelta63_r²ell²responseiscomputed.',
  'actual_g_original_root_pair':roots,'g_semantics':'delta63=aA63−aB63 inSAMEactualfull105/F9 branch, notgA/f63odd',
  'unknown105_original_order':U,'equation105_original_order':F,
  'origin_gradient_g_wrt105_sparse':[[i,str(q)]for i,q in sorted(grad.items())],
  'fixed_rho_equation_order_sparse':[[i,str(q)]for i,q in sorted(rho.items())],
  'fixed_rho_original_equation_nodes':[[F[i],str(q)]for i,q in sorted(rho.items())],
  'all105_rho_J_equals_grad_exact':True,'both_fixed105_inverse_products_verified':True,
  'no_new_inverse_or_source_solve':'TheexactstoredJ105inverseisapplied,notrecomputedbynumericIFT. OriginreverseADusesonlyexisting105originvaluesandoriginalrationalnodes.',
  'g_origin_AD_closure_original_records':g_records,'origin_unit_inverse_checks':unit_origins,
  'selected4_formula':'delta63_(r²ell²)=g_lower4_(r²ell²)−Σrho_i*F_i(lowerPhi≤3)_(r²ell²). g_lower4andRHS4areUNCOMPUTED. LinearCA_F9targetshavezeroorder4.',
  'future_B9':'QQ[r,ell]/(r³,ell³); onlyorder4monomialr²ell². LowerPhiinputsgiventhrough3, slotr²ell²set0ONLYforlower-term/RHSassembly, neverassumedthetrueimplicitPhi4zero.',
  'future_needed_original_records':future_records,'future_functional_roots':needed_roots,
  'future_given105_plane_jets10':joint['all105_unknown_jets'],
  'future_only_ONE_implicit_scalar_response':'No105fourthvector/fulljet/gauge4/newphysicalsolve. Exactfixedrhoandalllowercarrymustbecombined beforevalue.',
  'future_stream_lastuse_nodes':peak,
  'source_profiles':{'phase0_owned_upper':source_phase0,'originAD_owned_upper':arithmetic_phase,'max_origin_AD_precheck_bits':MAXPRE},
  'future_positive_source_code_gate_OPEN':'NeedB9every-nodemajorantandallpostcoherentoperator/minorheight,ownershipincludingmetadata/copies,RAW,preparedworker/custody/freshSOURCEpeer/ROOTONEGO BEFOREanyRHS4/truefourthvalue.',
  'all_new_fourth_values_still_formal':True,'no_originalH0_or_nilclaim':True,
  'selective_source_receipts':receipts,
 }
 ep=out/'one-selected-delta63-fourth-origin-adjoint-and-unrun-B9-source-closure.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_ONE_ORIGIN_RHO_FOR_SELECTED_DELTA63_FOURTH_NO_RESPONSE',
  'rho_support':len(rho),'rho_zero':not bool(rho),'future_original_records':len(future_records),'future_stream_lastuse_nodes':peak,
  'true_fourth_or_RHS4_response_evaluated':False,'higher_value_GO_closed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_source_owned_upper':max(source_phase0,arithmetic_phase),'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
