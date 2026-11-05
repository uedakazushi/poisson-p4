#!/usr/bin/env python3
"""SOURCE: a socle-insensitive, lower-cache-only delta63 finite algebra."""
import hashlib,json,math,signal,time
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_B63_SAFE_OUTPUT_SOURCE_20261105_v63'
PATHS={
 'gradient':'reports/AMBIENT_B18_C_HEAD29_SELECTED4_ADJOINT_SOURCE_20261105_v63/one-selected-delta63-fourth-origin-adjoint-and-unrun-B9-source-closure.json',
 'gradient_peer':'reports/AMBIENT_B18_C_HEAD29_SELECTED4_ADJOINT_SOURCE_20261105_v63_REVIEW/report.json',
 'gradient_counter':'reports/AMBIENT_B18_C_HEAD29_SELECTED4_ADJOINT_SOURCE_20261105_v63/whole-source-counter-completion-no-math-rerun.json',
 'head29':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63/coherent-head29-high-operator-and-three-selected-plane4-source-gate.json',
 'head29_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63_REVIEW/report.json',
 'Q20':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63/phase1-common105-and-114-source-equations.json',
 'Q20_peer':'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63_REVIEW/report.json',
 'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
 'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def basis():
 plane=[m for d in range(5)for m in combinations_with_replacement([7,8],d)]
 mixed=[(i,)+m for i in range(7)for d in range(3)for m in combinations_with_replacement([7,8],d)]
 zz=list(combinations_with_replacement([0,4,5],2))
 return sorted(plane+mixed+zz,key=lambda m:(len(m),m))
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 d={}
 fields={'gradient':['actual_g_original_root_pair','origin_gradient_g_wrt105_sparse','fixed_rho_equation_order_sparse'],
  'Q20':['base27_monomials','all105_cached_composite_values'],
  'joint':['all105_unknown_original_node_order','all105_equation_original_node_order','all105_unknown_jets','CA_F9_targets']}
 for k,p in PATHS.items():
  if k in fields:d[k]=S.selective(ROOT/p,fields[k],receipts)
  elif k.endswith('_peer'):
   d[k]=json.loads((ROOT/p).read_text());assert d[k]['all_passed']
 assert d['gradient']['actual_g_original_root_pair']==[20144,40922]
 assert d['gradient']['origin_gradient_g_wrt105_sparse']==d['gradient']['fixed_rho_equation_order_sparse']==[]
 assert d['joint']['CA_F9_targets']==['0']*7+['r','s']
 mon=basis();mi={m:i for i,m in enumerate(mon)};assert len(mon)==63 and mon[0]==()
 for m in mon:
  for i in set(m):assert tuple(list(m[:m.index(i)])+list(m[m.index(i)+1:]))in mi
 table=[[mi.get(tuple(sorted(a+b)),-1)for b in mon]for a in mon]
 # Downward closure means the complement is a monomial ideal.  These
 # index-only checks also verify the finite product implementation.
 assert all(table[0][i]==table[i][0]==i for i in range(63))
 variables=[mi[(i,)]for i in range(9)];comm=0
 for i in variables:
  for j in variables:
   for k in range(63):
    a=table[i][k];b=table[j][k]
    assert (-1 if a<0 else table[j][a])==(-1 if b<0 else table[i][b]);comm+=1
 socle=[i for i in range(63)if all(table[j][i]<0 for j in variables)]
 unknown=[i for i,m in enumerate(mon)if (all(j>=7 for j in m)and len(m)==4)
  or(len(m)==3 and sum(j<7 for j in m)==1)or(len(m)==2 and all(j<7 for j in m))]
 assert socle==unknown and len(unknown)==32
 known=[i for i in range(63)if i not in set(unknown)];assert len(known)==31
 plane10=[m for m in mon if all(i>=7 for i in m)and len(m)<=3]
 b27=[tuple(m)for m in d['Q20']['base27_monomials']]
 assert len(plane10)==10 and len(b27)==27 and set(b27)<=set(mon)
 maps=[]
 for target in[plane10,b27]:
  tm={m:i for i,m in enumerate(target)};proj=[tm.get(m,-1)for m in mon];checks=0
  for i,a in enumerate(mon):
   for j,b in enumerate(mon):
    p=table[i][j];lhs=-1 if p<0 else proj[p]
    rhs=tm.get(tuple(sorted(a+b)),-1)if proj[i]>=0 and proj[j]>=0 else -1
    assert lhs==rhs;checks+=1
  maps.append({'target_monomials':[list(m)for m in target],'projection63':proj,'all_product_checks':checks})
 # Only supplied INPUT coefficients are converted/compared here.  No
 # rational source output is composed and no implicit solution is found.
 joint=dict(d['joint']['all105_unknown_jets']);old=dict(d['Q20']['all105_cached_composite_values'])
 U=d['joint']['all105_unknown_original_node_order'];assert set(U)==set(joint)==set(old)and len(U)==105
 r_s=[(a,b)for n in range(4)for b in range(n+1)for a in[n-b]]
 lower=[];overlap=0;conversion_terms=0
 for node in U:
  arr=list(map(Q,joint[node]));v={}
  assert len(arr)==10
  for q,(a,b)in zip(arr,r_s):
   for ell in range(b+1):
    m=(7,)*(a+b-ell)+(8,)*ell
    v[m]=v.get(m,Q(0))+q*Q(math.comb(b,ell),5**ell);conversion_terms+=1
  w=list(map(Q,old[node]));assert len(w)==27
  for m,q in zip(b27,w):
   if all(i>=7 for i in m):assert v[m]==q;overlap+=1
   else:v[m]=q
  assert set(v)=={mon[i]for i in known}
  lower.append([node,[str(v.get(m,Q(0)))for m in mon]])
 assert overlap==630
 # This table is a lower extension, not an assertion that unprovided true
 # implicit coefficients are zero.  An arbitrary difference lives in the
 # 32-dimensional socle.  The formal rational chain rule then applies.
 owned=S.count(d)+S.count(table)+S.count(maps)+S.count(lower)+S.count(mon)+S.count(known)+S.count(unknown)+12000
 assert owned<150000
 data={
  'scope':'SOURCE lower-cache sufficiency for entire delta63 in B63. No delta63 new coefficient, odd root, head29/head30 response or source solve evaluated.',
  'same_regular_R0_variables':['z0','z1','z2','z3','z4','z5','z6','r','ell'],
  'basis63_monomials':[list(m)for m in mon],
  'basis_counts':{'plane_le4':15,'zi_times_plane_le2':42,'zz_active_0_4_5':6},
  'product63_indices_minus1_zero':table,'all_variable_commutations_on_basis':comm,
  'socle32_indices':socle,'known31_indices':known,
  'socle_checks':{'maximalideal_times_unknown32_zero':32*9,'socle_squared_zero':32*32},
  'trusted_quotient_homomorphisms':maps,
  'provided105_lower_extension63':lower,
  'all630_provided_plane_le2_overlaps_exact':True,'input_conversion_term_count':conversion_terms,
  'lower_extension_warning':'Zero entries on socle32 encode a chosen lower extension ONLY. True Phi4, zi-plane2 and activezz input coefficients remain arbitrary formal unknowns. No equation105 is claimed in whole B63 for this lower extension.',
  'safe_output_theorem':'Let Phi and L have the same residue Phi0 and agree on known31. Then Delta=Phi-L lies in Soc(B63)^105. For every original rational DAG with origin-unit denominators, f(Phi)-f(L)=sum_i(df/dPhi_i)(Phi0)*Delta_i, because m*Soc=0 and Soc²=0. The actual f=20144-40922 has ALL105 derivatives0 by independent forward-AD peer, hence its ENTIRE B63 image is identical for Phi and L.',
  'formal_chain_rule_induction':{'add':'linear identity','multiply':'(ab)(Phi)-(ab)(L)=a0*Delta_b+b0*Delta_a, Delta_a*Delta_b=0','inverse':'(a^-1)(Phi)-(a^-1)(L)=-a0^-2*Delta_a; a0 unit','rational_nodes':'Finite induction applies to original c/a/m/n/i records without evaluating any requested coefficient.'},
  'odd_output_scope':'For bA/bB/nu or other nonzero-gradient functions the full B63 output is NOT asserted true. Only the trusted plane10 or B27 quotient is valid. Any later worker must tag/extract precisely those outputs and reject untrusted slots.',
  'needed_delta_outputs':'Head29 plane4; head30 plane4, zi-plane2 and activezz2. All belong to the safe delta63 image, including adapted-curve and ORIGINAL-z Hadamard contractions.',
  'head30_source_gate_still_OPEN':'Coherent full head30 q7/H/m1/m2 source polynomial and primitive root availability/majorants/code/custody/freshpeer/ROOT ONEGO remain to be frozen before values.',
  'all_new_output_values_still_uncomputed':True,'no_new_IFT_source_gauge_or_physical_jet':True,
  'whole_SOURCE_numeric_upper':owned,'selective_source_receipts':receipts,
 }
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 ep=out/'finite63-socle-safe-output-theorem-and-actual-lower-input-overlaps.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_B63_SOCLE_SAFE_ENTIRE_DELTA63_NO_OUTPUT_VALUES',
  'safe_basis_dimension':63,'unknown_socle_dimension':32,'all630_input_overlaps':True,
  'actual_delta_or_head_responses_evaluated':False,'higher_value_GO_closed':False,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{ep.name:sha(ep)},'whole_source_numeric_upper':owned,'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
