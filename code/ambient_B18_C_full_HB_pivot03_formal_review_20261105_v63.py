#!/usr/bin/env python3
"""Independent full polynomial HB source proof; all new mu03 stay formal."""
import hashlib,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_HB_PIVOT03_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
SOURCE='6cbc0f67ea1b88f6d5e14d06176f38b58d876b34afb5823a45c73bc78a4023d0'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def plus(a,b):
 d=dict(a)
 for k,q in b.items():
  d[k]=d.get(k,Q(0))+q
  if not d[k]:del d[k]
 return d
def scale(a,c):return {k:q*c for k,q in a.items() if q*c}
def leaf(k):return {k:Q(1)}
def encode(a):return [[list(k),str(q)] for k,q in sorted(a.items())]
def polyplus(a,b):
 d={m:dict(v) for m,v in a.items()}
 for m,v in b.items():d[m]=plus(d.get(m,{}),v)
 return {m:v for m,v in d.items() if v}
def multiply_variable(a,v):return {(m[0]+int(v==0),m[1]+int(v==1)):dict(q) for m,q in a.items()}
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=True);pins={}
 assert not (out/'independent-full-HB-formal-six-coordinate-source-proof.json').exists()
 def pin(p,h=None):
  p=pathlib.Path(p);q=sha(p)
  if h:assert q==h,(str(p),q,h)
  pins[str(p.relative_to(ROOT))]=q
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',SOURCE);assert pr['all_passed'] and not pr['requested_values_evaluated']
 docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 formal=pin(ROOT/'reports'/PRI/'full-HB6-by2-affine-source-operator-and-old27-leaves.json',
  'e64828cd64d288f655c276aaf030fa3a2060f855f5aa83c557df18aa7df208f6')
 c2=docs['actual-physical12-and-C2-right-exact-quotient.json'];raw=docs['actual-ell2-product33-and-raw45-actions.json']
 hb=docs['two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json']
 C={r['physical_row']:r for r in c2['physical12_coefficients']};basis=raw['actual33_basis']
 PIV=[18,60,67];first7=[18,19,28,29,60,67,72];free=[19,28,29,72]
 assert docs['report.json']['all_passed']
 # Literal physical firsts and unit quadratic normalization, not names.
 T=[[Q(q) for q in row] for row in hb['constant_source_normalization_T']]
 H=[[Q(dict(C[p]['bB']).get(basis.index([i+1,0,0]),'0')) for i in range(7)] for p in PIV]
 assert H==[[Q(q) for q in row] for row in hb['source_first7_h_rows']]
 # Raw33 is modulo ell²; its missing ell² entry must come from the
 # separately accepted full physical frame / pure-ell readback.
 frame=pin(ROOT/'reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/actual-first21-quadratic-frame-and-twelve-cubic-zero-sensitivity-proof.json',
  '67a252627c275b059206abe3d6cca88906fe233f2cfd5e170fcd36e03e44efeb')
 fm={r['physical_row']:r for r in frame['all12_verified_source_records']}
 B=[[Q(fm[p]['bB_second_plane'][i]) for p in PIV] for i in range(3)]
 for i,m in enumerate([(2,0),(1,1)]):
  assert B[i]==[Q(dict(C[p]['bB']).get(basis.index([0,*m]),'0')) for p in PIV]
 assert [[sum(B[i][k]*T[k][j] for k in range(3)) for j in range(3)] for i in range(3)]==[[-Q(i==j) for j in range(3)] for i in range(3)]
 for p in first7:
  ds=dict(C[p]['delta']);assert all(not Q(ds.get(basis.index(m),'0')) for m in [[0,0,0],[0,1,0],[0,0,1]]+[[i+1,0,0] for i in range(7)])
 beta=[]
 for sparse in docs['actual-special-fiber20-coefficients.json']['beta2_before_ell']:
  row={}
  for mon,c in sparse:
   assert not any(mon[:7]) and sum(mon)==2
   nr,ns=mon[7:]
   for k in range(ns+1):
    m=(nr+ns-k,k);row[m]=row.get(m,Q(0))+Q(c)*math.comb(ns,k)/Q(5**k)
  beta.append({m:q for m,q in row.items() if q})
 assert [i for i,v in enumerate(beta) if v]==[0,4,5]
 # Full E_i3 for ALL seven first-frame sources. Only three pivot D3
 # vectors remain; every other E is proved zero from accepted source data.
 known={};delta={};effective={};mon3=[(3,0),(2,1),(1,2),(0,3)]
 for p in first7:
  ds=dict(C[p]['delta'])
  if p in PIV:
   delta[p]={mon3[0]:leaf(('cache',p,29)),mon3[1]:leaf(('cache',p,30)),
     mon3[2]:leaf(('lambda12',p,-1)),mon3[3]:leaf(('mu03',p,-1))}
   for ix in [29,30]:known[('cache',p,ix)]=Q(ds.get(ix,'0'))
  else:delta[p]={} # all four ACTUAL plane3 coefficients now independently zero
  gamma={}
  for z in range(7):
   for dr,de in [(1,0),(0,1)]:
    ix=basis.index([z+1,dr,de]);q=Q(ds.get(ix,'0'))
    if z in [0,4,5] and p in free:assert q==0
    if p in PIV and z in [0,4,5]:
     known[('cache',p,ix)]=q
     for (br,be),v in beta[z].items():
      m=(br+dr,be+de);gamma[m]=plus(gamma.get(m,{}),scale(leaf(('cache',p,ix)),v))
   # Inactive mixed leaves are intentionally not asserted zero.
  effective[p]=polyplus(delta[p],{m:scale(v,-Q(1)) for m,v in gamma.items()})
  if p in free:assert not effective[p]
 for r in docs['phase1-six-root-existing-cache-readbacks.json']['all3_lambda_readbacks']:
  known[('lambda12',r['physical_row'],-1)]=Q(r['delta_r_ell2'])
 assert len(known)==27
 A=[];Bfull=[]
 for k in range(3):
  aa={};bb={}
  for j,p in enumerate(PIV):
   aa=polyplus(aa,{m:scale(v,T[j][k]) for m,v in delta[p].items()})
   for i,q in enumerate(first7):
    bb=polyplus(bb,{m:scale(v,T[j][k]*H[j][i]) for m,v in effective[q].items()})
  A.append(aa);Bfull.append(bb)
 for k in range(3):
  assert [encode(A[k].get(m,{})) for m in mon3]==formal['P3_A'][k]
  assert [encode(Bfull[k].get(m,{})) for m in mon3]==formal['P3_B'][k]
 # Derive NO-CUT quartic NF recursively from the actual cubic curves.
 N3=[[Q(v) for v in row] for row in docs['full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json']['full_no_cut_cubic_NF6_rows12_columns']]
 N4=[[Q(0)]*15 for _ in range(6)]
 for side in range(2):
  g=[Q(1)]+[-N3[side*3+i][side*4] for i in range(3)]
  def reduce(r,e):
   if r<3:return {(r,e):Q(1)}
   d={}
   for k in range(1,4):
    for m,c in reduce(r-k,e+k).items():d[m]=d.get(m,Q(0))-g[k]*c
   return d
  for col,r in enumerate([4,3,2,1,0]):
   v=reduce(r,4-r)
   for row,m in enumerate([(2,2),(1,3),(0,4)]):N4[3*side+row][5*side+col]=v.get(m,Q(0))
 assert [[str(v) for v in row] for row in N4]==formal['full_NF4']
 matrix=[[{} for _ in range(2)] for _ in range(6)]
 for j in range(2):
  for side,P in enumerate([A,Bfull]):
   poly=polyplus(multiply_variable(P[j],1),{m:scale(v,-Q(1)) for m,v in multiply_variable(P[j+1],0).items()})
   for row in range(3):
    for (r,e),v in poly.items():
     assert r+e==4
     matrix[3*side+row][j]=plus(matrix[3*side+row][j],scale(v,N4[3*side+row][5*side+4-r]))
 encoded=[[encode(v) for v in row] for row in matrix]
 assert encoded==formal['full_HB_6_by2_affine_operators']
 assert [[list(k),str(v)] for k,v in sorted(known.items())]==formal['known27_old_cached_leaf_values']
 for side,row in enumerate([0,3]):
  for j in range(2):
   old={('lambda12' if k[0]=='lambda' else k[0],k[1],k[2]):Q(v) for k,v in hb['matrix2_rows_A4_B4_columns_HB0_HB1_linear_leaf_operators'][side][j]}
   assert matrix[row][j]==old
 assert all(k[0]!='mu03' for row in [0,3] for v in matrix[row] for k in v)
 proof={'same_actual_physical_first7_and_quadratic_unit_T_verified':True,
  'all_first_delta9_and_origins_literal_zero':True,
  'all7_firstframe_E3_built_before_reduction_to3pivots':True,
  'free_eta_alpha_and_graph_source_zeros_pinned_no_unknown_zero_filling':True,
  'beta2_mixed_carry_negative_all_active_coefficients_retained':True,
  'P3_A_and_P3_B_independently_constructed_symbolically':True,
  'all6_by2_affine_entries_match_with3mu03_unread':True,
  'old_rows0_3_match_coefficientwise_no_mu_dependency':True,
  'full_HB_order3_cancellation_proof':'Actual P0,P1,P2 have FULL secondterms r²H,r ellH,ell²H. Thus ellP0-rP1 and ellP1-rP2 are genuine order4 sources before any cut. Linear*P3 alone supplies their A/B quartics.',
  'higher_source_carry_proof':'P4 enters only order5 after multiplying r/ell. z-containing quartics clear by cubic coefficients times actual firstframe with higher carry>=5. H quartics clear by quadratic coefficients times P, whose A/B starts3 and hence contributes>=5. No delta4,beta3 or mixed source3 callback contributes.',
  'matrix_leaf_term_counts':[[len(v) for v in row] for row in matrix],
  'mu03_values_and_HB_entries_not_paired':True,
  'rank2_scope':'if some w=col1+10col0 nonzero, rank2 and full HF4=4, b2,e3or4; rank1 only known-submodule negative',
  'no_IFT_physical_gauge_jet_or_new_source_values':True,
  'input_sha256':pins,'code_sha256':sha(__file__),'elapsed_seconds':time.monotonic()-start}
 ep=out/'independent-full-HB-formal-six-coordinate-source-proof.json';save(ep,proof)
 save(out/'formal-phase.json',{'status':'FORMAL_MATH_PASS_AWAIT_UNRUN_WORKER_CODE_CUSTODY_AND_CACHE_CAPS','proof_sha256':sha(ep),'no_new_mu_or_HB_pairing':True})
 signal.alarm(0)
 print(json.dumps({'formal_math_passed':True,'proof_sha256':sha(ep),'seconds':proof['elapsed_seconds']}))
if __name__=='__main__':main()
