#!/usr/bin/env python3
"""Source/bound gate for TWO cached pairings only; no pairing evaluated."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
from ambient_B18_ell2_raw33_quotient_20261105_v63 import srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_ALPHA29_BETA29_INPUT_GATE_20261105_v63'
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=H.sha(ROOT/p);return H.selected(p,k)
 graph=read('reports/AMBIENT_B18_ELL2_T2_LINEAR_GRAPH4_SOURCE_20261105_v63/actual23-linear-graph-four-common-contractions-and-six-adjoint-input.json');H0=list(map(H.parse,graph['H0_free3_columns']));assert H0==[{23:F(1)},{29:F(1)},{30:F(1)}];assert graph['two_gauge_invariants']=='alpha=mu[0,r], beta=mu[0,ell]-mu[1,r]';W=H.parse(graph['actual_H0v_source56']);assert W=={23:F(-3375,11),29:F(-50),30:F(25,2)}
 M=list(map(H.parse,field('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json','Mnew0_columns')));assert all(not M[j] for j in [23,29,30])
 source='reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63/actual3960-rho-factor-identities-or-obstructions.json';Rho=list(map(lambda row:list(map(H.parse,row)),field(source,'actual_two_three_component_bar_rho33')))
 assert read('reports/AMBIENT_B18_ELL2_T2_GRAPH_RAW33_FACTOR_20261105_v63_REVIEW/report.json')['all_passed'];assert read('reports/AMBIENT_B18_ELL2_T2_GRAPH_MIXED_ADJOINT_20261105_v63_REVIEW/report.json')['all_passed']
 phy=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json');p=phy['physical12_coefficients'][7];assert p['physical_row']==63 and p['roots']['physical_row']==63 and p['roots']['aA_node']==20144 and p['roots']['aB_node']==40922 and p['roots']['bB_node']==38648 and p['roots']['delta_operator']==['subtract',20144,40922]
 assert read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed'];V29=[H.parse(p['delta']),{}, {i:-x for i,x in H.parse(p['bB']).items()}];V30=[{},H.parse(p['delta']),H.parse(p['bA'])]
 # These are previous cached vector VALUES, not new rho pairings.
 pW='reports/AMBIENT_B18_ELL2_T2_GRAPH_CACHE_VALUES_20261105_v63/phase1-actual-two-invariants-from-raw33-primitive.json';vw=read(pW);VW=list(map(H.parse,vw['primitive_V_Wv_full33']));aW=list(map(F,vw['actual_two_invariants_alpha_beta']))
 p23='reports/AMBIENT_B18_ELL2_T2_ACTIVE4_CACHE_BUNDLE_20261105_v63/phase1-actual-alpha23-beta23-and-primitive23.json';v23=read(p23);V23=list(map(H.parse,v23['actual_primitive_V23_raw33']));a23=list(map(F,v23['actual_alpha23_beta23']))
 for stem in ['GRAPH_CACHE_VALUES','ACTIVE4_CACHE_BUNDLE']:
  assert read('reports/AMBIENT_B18_ELL2_T2_'+stem+'_20261105_v63_REVIEW/report.json')['all_passed']
 maximum=0;checks=[]
 def mark(b):
  nonlocal maximum
  maximum=max(maximum,b[0].bit_length(),b[1].bit_length());return b
 def exact(x):x=F(x);return mark((x.denominator,abs(x.numerator)))
 def prod(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 dotbounds=[];derivebounds=[]
 for ch,row in enumerate(Rho):
  b=(1,0)
  for t,rr in enumerate(row):
   for i,c in rr.items():b=plus(b,prod(exact(c),exact(V29[t].get(i,F(0)))))
  dotbounds.append(b);d=plus(exact(aW[ch]),prod(exact(F(3375,11)),exact(a23[ch])));d=plus(d,prod(exact(F(50)),b));derivebounds.append(prod(exact(F(2,25)),d))
 # Prior for all99 exact source-vector linear identities; no rho dot here.
 for t in range(3):
  for i in set(VW[t])|set(V23[t])|set(V29[t])|set(V30[t]):
   b=exact(VW[t].get(i,F(0)))
   for c,v in [(F(-3375,11),V23),(F(-50),V29),(F(25,2),V30)]:b=plus(b,prod(exact(c),exact(v[t].get(i,F(0)))))
 assert maximum<=384
 H.dump(out/'phase.json',{'phase':'source_vector_binding_only_no_pairings','positive_prior_bits':maximum,'no_alpha29_beta29_value_evaluated':True})
 bits=0
 def ck(x):
  nonlocal bits
  bits=max(bits,abs(x.numerator).bit_length(),x.denominator.bit_length());assert bits<=maximum;return x
 for t in range(3):
  for i in range(33):
   q=F(0)
   for c,v in [(F(-3375,11),V23),(F(-50),V29),(F(25,2),V30)]:q=ck(q+ck(c*v[t].get(i,F(0))))
   assert q==VW[t].get(i,F(0));checks.append([t,i])
 data={'primitive_V29_raw63_mask1':[srec(v) for v in V29],'primitive_V30_raw63_mask2_source_check_only':[srec(v) for v in V30],'two_universal_bar_rho33':[[srec(v) for v in row] for row in Rho],
  'known_alpha_beta_on_Wv':list(map(str,aW)),'known_alpha23_beta23':list(map(str,a23)),'source_Wv_coefficients':[[j,str(x)] for j,x in sorted(W.items())],
  'recovered30_formula':'u30=(2/25)*(uW+(3375/11)*u23+50*u29), separatelyforalphaandbeta. Thisfollowsall99source-vectoridentitiesanduniversalrho-linearity;noextraV30pairingneeded.',
  'two_new_dot_positive_bounds':[list(map(str,b)) for b in dotbounds],'two_recovered30_positive_bounds':[list(map(str,b)) for b in derivebounds],
  'universal_source_condition':'M0e29=M0e30=M0e230; F=M e_j hasF0=0, V e_j=G1F. UniversalALL286functionalidentityincludingF2nuisanceandall3960raw33descentappliesbeforechoosingj. No assumptionofcanonicalmu1values.',
  'all99_Wv_source_vector_identities_checked':True,'no_pairing_or_new_function_values':True,'source_raw63_roots':p['roots'],'source_order':'R0(zA0..6,rA,ell), ell=5(sA-rA); physicaltrace variablesaretheacceptedtracefreecache,notF9sA'}
 future=2*(count(data)+4)+5000;raw=2*len(json.dumps(data).encode())+1000*(2*math.ceil(maximum*math.log10(2))+40)+1_000_000;owned=count([graph,M,phy,vw,v23,Rho,data])+8000;ready=future<60000 and owned<60000 and raw<32*1024*1024 and 3*maximum+1<=1153
 data.update({'whole_numeric_bound':future,'whole_RAW_bound':raw,'source_committed_bit_bound':maximum,'source_temporary_bit_bound':3*maximum+1});ep=out/'two-pairing-source-binding-and-positive-height-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':ready,'no_pairing_or_new_function_value_evaluated':True,'new_pairing_count_only':2,'all99_Wv_source_supports_verified':True,'source_committed_bit_bound':maximum,'source_temporary_bit_bound':3*maximum+1,'known_vector_check_maximum_bits':bits,'whole_numeric_bound':future,'whole_RAW_bound':raw,'source_owned_bound':owned,'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:H.sha(ep)},'seconds':time.monotonic()-start};H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','new_pairing_count_only','source_committed_bit_bound','source_temporary_bit_bound','whole_numeric_bound','whole_RAW_bound','source_owned_bound']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
