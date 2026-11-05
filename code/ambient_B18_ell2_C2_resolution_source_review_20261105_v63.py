#!/usr/bin/env python3
"""Independent C2 Loewy receipts, polynomial units and resolution proof."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,math,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv,mp,eye
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_C2_RESOLUTION_SOURCE_PREFLIGHT_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','filtration':f'reports/{A}/actual-C2-Loewy-filtration.json',
 'operator':f'reports/{A}/actual-C2-resolution-and11-unit-source-cancellation.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'C2_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'raw':'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def poly(p):return {tuple(m):Q(c) for m,c in p if Q(c)}
def vector(v):return {int(i):poly(p) for i,p in v if p}
def pe(p):return [[list(m),str(c)] for m,c in sorted(p.items())]
def ve(v):return [[i,pe(p)] for i,p in sorted(v.items()) if p]
def times(a,b):
 z={}
 for m,c in a.items():
  for n,d in b.items():
   k=tuple(x+y for x,y in zip(m,n))
   if k[-1]<2:z=add(z,{k:c*d})
 return z
def vadd(a,b,c=Q(1)):
 out={i:dict(p) for i,p in a.items()}
 for i,p in b.items():
  q=add(out.get(i,{}),scale(p,c))
  if q:out[i]=q
  else:out.pop(i,None)
 return out
def pmatmul(A,B):
 out=[[{} for j in B[0]] for i in A]
 for i,row in enumerate(A):
  for k,a in enumerate(row):
   for j,b in enumerate(B[k]):
    if a and b:out[i][j]=add(out[i][j],times(a,b))
 return out
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'ACTUAL14_ACTIONS_AND_TYPED_POLYNOMIAL_GATES','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items()};assert D['primary']['all_passed'] and D['C2_peer']['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 actions=[[sp(p) for p in M] for M in D['C2']['C2_coordinate_actions9_sparse_columns']];E=actions[8];R=actions[7]
 assert all(mp(M,N)==mp(N,M) for M in actions for N in actions) and all(not p for p in mp(E,E))
 current=eye(14);dims=[14]
 for stage in D['filtration']['stages']:
  columns=[mv(M,p) for M in actions for p in current]
  assert list(map(enc,columns))==stage['input_columns']
  fixed={}
  for receipt in stage['fixed_rank_receipts']:
   v=columns[receipt['source_column']]
   for k,c in receipt['subtract']:v=add(v,scale(fixed[k],-Q(c)))
   k=receipt['pivot'];assert k not in fixed and v.get(k)==Q(receipt['divide_by'])
   v=scale(v,1/Q(receipt['divide_by']));assert enc(v)==receipt['normalized'];fixed[k]=v
  for p in columns:
   v=dict(p)
   for k,b in sorted(fixed.items()):v=add(v,scale(b,-v.get(k,Q(0))))
   assert not v
  current=[fixed[k] for k in sorted(fixed)];assert list(map(enc,current))==stage['basis'] and len(current)==stage['rank'];dims.append(len(current))
 assert dims==[14,11,5,1,0] and D['filtration']['actual_Loewy_HF']==[3,6,4,1]
 op=D['operator'];rawindices=D['C2']['C2_basis_raw45_indices'];raw=D['raw'];branch=raw['canonical_branch_quotients6'][0]['basis_plane_monomials'];physical=[]
 for q in rawindices:
  if q<6:sector=0;a,e=branch[q]
  elif q<12:sector=1;a,e=branch[q-6]
  else:
   b,a,e=raw['actual33_basis'][q-12];assert b==0;sector=2
  physical.append((sector,a,e))
 eta=[rawindices.index(i) for i in [0,6,12]];assert eta==op['physical_target3_basis_positions']==[0,5,10]
 nonphysical=[i for i in range(14) if i not in eta];one={(0,0):Q(1)};ell={(0,1):Q(1)};r={(1,0):Q(1)}
 m0=[{s:{(a,e):Q(1)}} for s,a,e in physical]
 assert list(map(ve,m0))==op['physical_lift_m0_columns']
 for j,(s,a,e) in enumerate(physical):
  v={eta[s]:Q(1)}
  for _ in range(a):v=mv(R,v)
  if e:v=mv(E,v)
  assert v=={j:Q(1)}
 Avec=list(map(vector,op['source_a11_columns']));K=list(map(vector,op['target_K11_columns']));anchors=op['source_anchors11']
 assert len(Avec)==len(K)==len(anchors)==11 and len(set(anchors))==11
 for j,(a,k) in enumerate(zip(Avec,K)):
  expected=vadd({nonphysical[j]:one},{eta[physical[nonphysical[j]][0]]:{physical[nonphysical[j]][1:]:Q(1)}},-1)
  assert k==expected
  actual={}
  for source,p in a.items():
   if source<14:idx=source;v=ell;matrix=E
   else:assert 112<=source<126;idx=source-112;v=r;matrix=R
   actual=vadd(actual,{idx:times(v,p)})
   actual=vadd(actual,{i:scale(p,-c) for i,c in matrix[idx].items()})
  assert actual==k
  constants={i:p.get((0,0),Q(0)) for i,p in a.items() if p.get((0,0),Q(0))}
  assert constants=={anchors[j]:Q(-1)}
 Qm=[[Avec[j].get(i,{}) for j in range(11)] for i in anchors];Qi=[[poly(p) for p in row] for row in op['selected_Q_inverse11']]
 assert [[pe(p) for p in row] for row in Qm]==op['selected_Q11']
 I=[[one if i==j else {} for j in range(11)] for i in range(11)]
 N=[[add(Qm[i][j],I[i][j]) for j in range(11)] for i in range(11)]
 assert not any(p for row in pmatmul(N,N) for p in row) and op['selected_N_nilpotence_index']==2
 assert pmatmul(Qm,Qi)==pmatmul(Qi,Qm)==I
 # Verify the physical target splitting as two inverse polynomial maps.
 target=K+[{i:one} for i in eta]
 def target_inverse(v):
  out={i:v[j] for i,j in enumerate(nonphysical) if j in v}
  for j,p in v.items():
   for s,q in m0[j].items():out=vadd(out,{11+s:times(p,q)})
  return out
 assert all(target_inverse(v)=={j:one} for j,v in enumerate(target))
 complement=[i for i in range(126) if i not in anchors];assert complement==op['source_complement115_indices'] and len(complement)==115
 def embed(p):return {(0,)*7+(a,e):c for (a,e),c in p.items()}
 G1=[]
 for source in complement:
  if source<14:var=8;j=source;matrix=E
  else:var=(source-14)//14;j=(source-14)%14;matrix=actions[var]
  exponent=tuple(int(i==var) for i in range(9));v={}
  for s,p in m0[j].items():v[s]=times(embed(p),{exponent:Q(1)})
  for i,c in matrix[j].items():
   for s,p in m0[i].items():
    q=add(v.get(s,{}),scale(embed(p),-c))
    if q:v[s]=q
    else:v.pop(s,None)
  assert all((0,)*9 not in p for p in v.values());G1.append(v)
 assert [{i:p for i,p in v.items() if p} for v in G1]==list(map(vector,op['reduced_G1_3_by115_poly_columns']))
 ranks=[14*sum(math.comb(8,j) for j in range(min(8,i)+1)) for i in range(4)]
 assert ranks==[14,126,518,1302] and op['after11_unit_cancellation']==[3,115,518,1302]
 ep=o/'independent-Loewy3641-polynomial11units-and-resolution-proof.json'
 dump(o,ep.name,{'inputs':inputs,'fixed_filtration_receipts_and_all_image_spans_verified':True,'actual_Loewy_dimensions':dims,'actual_HF':[3,6,4,1],
  'all14_physical_r_ell_monomial_lifts_reproduce_actual_actions':True,'eta3_positions':eta,
  'all11_global_polynomial_d1_a_equals_K_verified':True,'source_anchors11':anchors,'selected_minor_has_origin_minus_identity':True,
  'N_squared_zero_and_global_two_sided_inverse_verified':True,'target_basis_two_sided_polynomial_unit_verified':True,
  'all3_by115_G1_polynomial_entries_verified':True,'full_resolution_first_ranks':ranks,'after11units':[3,115,518,1302],
  'G2_complete_typed_formula':'P_complement*(I-a11*Q11_inverse*P_anchors)*d2; source clearing does not change complement coordinates',
  'G3_complete_typed_formula':'original total periodic/Koszul d3 into unchanged F2 rank518',
  'exact_periodic_proof':'u+ell*w maps u+Lw; kernel augmentation=(ellI-L)w; ker(ellI-L)=(ellI+L)w and conversely, sinceL²0',
  'exact_Koszul_total_proof':'all8 Ai commute withL andeachother; completed monic translation xI-A onC2[[x]] is regular for any coefficient module, noA-flatness ofC2 required',
  'source_unit_cancellation_proof':'source basis[a,E_C-a*pi_nonphysical*d1*E_C] unit byQ11; target[K,eta] unit; d1d2=0 killsfirst11 transformed d2 coordinates',
  'Vext_optional_exact_kernel_proof':'in target[K,eta] coordinates [etaV,K] equalsblock[[0,I11],[V,0]], hence samekerV andsamecokerC2',
  'coarse_T2_if_source_maps_are_constructed':{'generators':574,'relations':2185,'required':['d1m=V','Vn=d1','d2h=I-mn','768cyclelifts']},
  'actual_comparison_m_n_h_or_T2_H0_values_not_computed':True,'new_source_coefficients':0})
 note=ROOT/'research/notes'/f'{STEM}.md';report={'id':STEM,'all_passed':True,
  'actual_C2_Loewy_HF':[3,6,4,1],'actual_C2_Loewy_dimensions':dims,
  'full_periodic_Koszul_resolution_and11_source_unit_cancellations_accepted':True,'reduced_first_ranks':[3,115,518,1302],
  'actual_comparison_maps_T2_H0_not_constructed':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),
  'note_sha256':sha(note),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
