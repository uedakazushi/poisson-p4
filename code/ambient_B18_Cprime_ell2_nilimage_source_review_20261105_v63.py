#!/usr/bin/env python3
"""Independent finite33 SOURCE/frame/budget and exact module-method peer.
No new gauge, physical source, jet, nilimage or connecting value is evaluated.
"""
import collections,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_ELL2_NILIMAGE_SOURCE_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def encode(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):encode(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [encode(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(encode(v),indent=2)+'\n')
 stem='AMBIENT_B18_CPRIME_ELL2_NILIMAGE_SOURCE_PREFLIGHT_20261105_v63';r=load(stem);assert pins[f'reports/{stem}/report.json']=='81da59ebf78a9f30b906784ea7b47ac57ca8a4fa91eca5ccb010d8e62522c2a7'and r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/n)==h
 P=load(stem,'source33-nilimage-Cprime2-and-formal-edge-gate.json')
 S=load('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json');T=load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-source105-ell2-raw33-values.json');A=load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63','exact-existing-raw33-multiplication-table.json');C2=load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-physical12-and-C2-right-exact-quotient.json')
 assert load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW')['all_passed']
 assert T['all105_equations_exactly_zero']and T['finite_actual_function_uniqueness']and C2['actual_C2_dimension']==14
 ids=set(T['joint105_unknown_order']);assert len(ids)==105;nodes=S['DAG_scalar_nodes'];roots=[i for f in S['families'].values()for pair in f['complete_normalized_physical_w_coefficients'][0]for i in pair]+[f['actual_tracefree_nu']for f in S['families'].values()];need=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in need:continue
  need.add(i)
  if i in ids:continue
  n=nodes[i];assert n[0]in['c','a','m','n','i']
  if n[0]!='c':todo.extend(n[1:])
 assert len(need)==P['source_ancestor_count']
 uses=collections.Counter(roots)
 for i in need:
  if i not in ids and nodes[i][0]!='c':uses.update(nodes[i][1:])
 live=set();peak=0
 for i in sorted(need):
  live.add(i);peak=max(peak,33*len(live))
  if i not in ids and nodes[i][0]!='c':
   for j in nodes[i][1:]:
    uses[j]-=1
    if not uses[j]:live.remove(j)
 assert peak==13299==P['dense_source_stream_peak_QQ']
 interface=load('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63','actual-source-physical-firstjet-interface.json');proper=load('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63','literal-proper-source14-preflight.json');J=load('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase2-fixed30-Jacobian.json')
 assert len(interface['actual_group15_positions'])==15 and len(proper['CA_fiber10_fixed_coefficient_indices'])==10 and len(proper['CA_free_fiber9_coefficient_indices'])==9
 assert interface['retained_source9_order']==['a','d','v2','v3','v4','v5','v6','v7','v8']
 a0=[list(map(F,v))for v in J['A0']];a1=[list(map(F,v))for v in J['A1_eta_delta']];ji=[list(map(F,v))for v in J['J30_inverse']];jj=[a0[i]+[F(0)]*15 for i in range(15)]+[a1[i]+a0[i]for i in range(15)];assert F(J['det_J30'])==F(4,9)
 assert all(sum((jj[i][k]*ji[k][j]for k in range(30)),F(0))==F(i==j)and sum((ji[i][k]*jj[k][j]for k in range(30)),F(0))==F(i==j)for i in range(30)for j in range(30))
 assert load('AMBIENT_B18_HONEST_NIL_T_GAUGE120_20261105_v63_REVIEW')['all_passed']
 nil=load('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63','actual-E4-nil-two-generators-and-relations.json');assert load('AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63_REVIEW')['all_passed']
 assert nil['old_variable_order']==['a','d','x','t']
 expected={'h':{(1,0,0,0):F(1),(0,1,0,0):F(-1),(0,0,1,0):F(1),(0,0,0,1):F(-2)},'c':{(0,1,0,0):F(1,2),(0,0,1,0):F(1,2)},'l':{(1,0,0,0):F(1,3),(0,1,0,0):F(1,3),(0,0,1,0):F(-1,3)}}
 assert all({tuple(e):F(c)for e,c in nil['forward_GL4'][k]}==v for k,v in expected.items())
 # Independent monomial filtration and all1089 plane homomorphism checks.
 basis=A['basis33'];degrees=[int(z!=0)+r+e for z,r,e in basis];assert degrees==A['basis_degrees']and max(degrees)==4 and degrees[0]==0 and min(degrees[1:])==1
 table=A['table_columns_for_each_left_basis'];assert all(degrees[k]>=degrees[i]+degrees[j]for i,row in enumerate(table)for j,col in enumerate(row)for k,c in col if F(c))
 plane=[(0,0),(1,0),(0,1),(1,1)];index={m:i for i,m in enumerate(plane)};mp=[{index[r,e]:F(1)}if not z and(r,e)in index else{}for z,r,e in basis]
 def put(v,k,c):
  v[k]=v.get(k,F(0))+c
  if not v[k]:del v[k]
 def mul(x,y):
  out={}
  for i,c in x.items():
   for j,d in y.items():
    mon=tuple(a+b for a,b in zip(plane[i],plane[j]));k=index.get(mon)
    if k is not None:put(out,k,c*d)
  return out
 for i in reversed(range(33)):
  for j in reversed(range(33)):
   p={}
   for k,c in table[i][j]:
    for n,d in mp[k].items():put(p,n,F(c)*d)
   assert p==mul(mp[i],mp[j])
 assert [[[i,str(c)]for i,c in p.items()]for p in mp]==P['formal_plane4_quotient_maps33']
 # Trace nilpotence proof: with base m^5=0, a surviving tauA^u*tauB^v
 # times a base monomial of degree t needs t+floor(u/2)+floor(v/2)<=4.
 maxdegree=max(t+u+v for t in range(5)for u in range(11)for v in range(11)if t+u//2+v//2<=4);assert maxdegree==10
 rank=load('AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63');assert rank['actual_generic_DVR_rank_Cprime']==2 and rank['actual_DVR_torsion_minimal_generators']==2
 assert load('AMBIENT_B18_NIL_DVR_FINITENESS_20261104_v63_REVIEW')['all_passed']
 assert P['expected_Cprime2_interval_from_rank2_and_tors_mu2']==[6,8]
 for e1,e2 in itertools.product(range(1,4),repeat=2):
  length=4+min(e1,2)+min(e2,2);assert 6<=length<=8 and((length==6)==(e1==e2==1))
 budget=P['hard_budget'];assert sum(budget.values())==118256==P['whole_saved_live_QQ_upper'];assert budget['gauge30_base33']==30*33 and budget['four_gauge_residual_correction_stages']==4*2*30*33 and budget['physical_source_DAG_dense_peak']==peak
 assert P['new_derived_gauge_coefficient_upper']==990 and P['hard_seconds']==60 and P['hard_total_output_bytes']==32<<20 and P['future_derived_coefficient_bit_guard']==2048
 proof=dict(source105_closed_original_DAG_ancestors=len(need),source_dense_peak_QQ=peak,all30_fixed_unit_Jacobian_left_right_inverse_identities=1800,
  same_frame='same actual15 group positions, CA flag5+fixed fiber10; retained kappa9[a,d,v2,...,v8], oldcore(a,d,x,t)=(kappa0,kappa1,kappa2,kappa6)',
  actual_nilgens='d*h,d*c*l, h=a-d+x-2t,c=(d+x)/2,l=(a+d-x)/3; retain all5 d*u relations and fullactual proper values',
  all1089_plane_homomorphism_products=True,base_m5_zero_from_degree_filtration=True,trace_maximal_ideal_power11_zero=True,
  right_exact_quotient_proof='NE is generated overE4 by2 nilgens; E4 is generated overR0 by4 tracebasis classes, hence eightimages generate its image. C2 already factors throughR0/(IA*IB,ell²), so33basis closure spans the actual image. Cprime/ell²=(C/ell²)/im(NE), by right-exactness only; no H1 basechange.',
  length6_target_proof='Cprime=D^2 plus torsion with2 elementary divisors. lengthCprime2=4+min(e1,2)+min(e2,2). Only length6 forces both e=1. Then genuineCprime[ell] injects Cprimebar and equals ker(ell:Cprime2)/imell there. Length7/8 leaves genuine target open.',
  formal_edge_proof='For a FORMAL ell-independent lift sigma(cbar), Vfull*sigma(cbar) is divisible byell by the exactcycle identity. Its quotient class is independent of the lift modulo imVfull and lies inCprime[ell]. Its reduction is (partialell Vfull)|ell0*cbar. No division occurs inArtin33.',
  remaining_gates=['Before values pin accepted B11 gauge/proper cache and C_special all26 ell0 projection; verify actual150 physical, all30equations, gaugeorigin+ell0overlap, nil2image overlap and Cprimebar length4/actions.','Keep source m1/regular-division and cycle-first sensitivities open until separately canceled or supplied. No14nilcycle or actualconnecting matrix is yet built.','Enforce peak118256,2048bit and60s/2GB/32MiB guards; stop rather than enlarge if crossed.'],
  source_gate_passes_with_mandatory_runtime_checks=True,whole_saved_live_QQ_upper=118256,new_physical_source_or_jet_required=False,canonical990_values_NOT_computed=True,Cprime2_length_and_actual14kernel_NOT_computed=True)
 save('independent-frame-unit-plane-and-right-exact-method-proof.json',proof)
 for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,source_frame_budget_gate_passes=True,whole_QQ_upper=118256,actual_gauge_nilimage_or_connecting_values_computed=False,finite_H1_or_ell_division_shortcut_rejected=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20;save('report.json',result);print(json.dumps({k:v for k,v in result.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
