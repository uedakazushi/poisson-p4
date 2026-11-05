#!/usr/bin/env python3
"""Independent full quartic witness and DVR embedding/tail proof."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_FULL_QUARTIC_BOUND_SOURCE_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='deed64d50a1fd9271cfddb7256257b8d4c769292c0c138f21f6e0c3c09c0ff1f'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={}
 def pin(p,h=None):
  p=pathlib.Path(p);v=sha(p)
  if h:assert v==h,(str(p),v,h)
  pins[str(p.relative_to(ROOT))]=v
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 docs={}
 for p,h in pr['input_sha256'].items():
  d=pin(ROOT/p,h);docs[pathlib.Path(p).name]=d
  if p.endswith('_REVIEW/report.json'):assert d['all_passed']
 D=pin(ROOT/'reports'/PRI/'full-quartic-known-HB-witness-and-m2-torsion-shape.json',
   '3185353ccd54d99a99246f90153ac02a7defca391b50be16d23ececa943830ad')
 pin(ROOT/'code/ambient_B18_C_full_quartic_bound_source_20261105_v63.py',pr['code_sha256'])
 pin(ROOT/'research/notes'/f'{PRI}.md',pr['note_sha256'])
 R=docs['full-C-ell-reduction-and-conditional-graph03-exponent3-proof.json']
 N3=[[Q(q) for q in row] for row in R['full_no_cut_cubic_NF6_rows12_columns']]
 N4=[[Q(0)]*15 for _ in range(6)]
 for t in range(2):
  g=[Q(1)]+[-N3[3*t+k][4*t] for k in range(3)]
  def reduce(i,j):
   if i<3:return {(i,j):Q(1)}
   d={}
   for k in range(1,4):
    for m,c in reduce(i-k,j+k).items():d[m]=d.get(m,Q(0))-g[k]*c
   return d
  for c,i in enumerate([4,3,2,1,0]):
   v=reduce(i,4-i)
   for k,mon in enumerate([(2,2),(1,3),(0,4)]):N4[3*t+k][5*t+c]=v.get(mon,Q(0))
 assert [[str(q) for q in row] for row in N4]==D['known_full_quartic_NF6_rows15_columns']
 oldA=list(map(Q,['44/625','-8/25','1','0','0']))
 oldB=list(map(Q,['6252/75625','-96/275','1','0','0']))
 assert N4[0][:5]==oldA and N4[3][5:10]==oldB
 HB=docs['two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json']
 assert HB['HB_sources']==[{'label':'HB0','formula':'ell*P0-r*P1'},{'label':'HB1','formula':'ell*P1-r*P2'}]
 Hval=docs['phase2-actual-RIGHT2-quartics-two-factors-and-rank.json']
 M=[[Q(q) for q in row] for row in Hval['actual_quartic_matrix_rows_A4_B4_columns_HB0_HB1']]
 assert M==[[-Q(3,5),Q(6)],[-Q(121,80),Q(121,8)]]
 # Actual full P0,1 secondterms r²H,r ellH give order3 cancellation
 # before any ell³ cut. Unknown ell³ coefficients in P3 survive in the
 # full quartic rows1/2/4/5; the old rows0/3 annihilate them literally.
 assert oldA[3:]==[0,0] and oldB[3:]==[0,0] and M[0][0]
 embeddingproof=[
  'Set M=m²C and T=Tor_D C. Finite-D reduction gives rankD M4, muD M=HF2=6, and muD ellM=HF3=6.',
  'Thus S=Tor_D M has exactly two cyclic summands D/ell^a,D/ell^b; both a,b>=2 because none disappears after multiplication by ell.',
  'S=M intersect T embeds in T of actual type(1,1,2,e). Hence ell²S embeds in ell²T, which is CYCLIC (or zero). A submodule of a cyclic module over the DVR D is cyclic.',
  'Therefore at most one of a,b exceeds2, so after sorting a=2 and b>=2. This gives HF4=muD ell²M<=5 independently of the HB calculation.',
  'Together with rankD4 this proves HF4 in{4,5}. No quartic completeness or vanishing of unknown ell³/ell⁴ coordinates is assumed.',
  'All decompositions here are D-module structure theorem decompositions, neither R0-module nor algebra splittings.'
 ]
 exponentproof=[
  'ell²T subset M because ell²C subset m²C; it is torsion, hence ell²T subset S. The largest exponent e-2 of ell²T is at most b. Thus e<=b+2.',
  'Actual C3m4=0 gives ell²M=m4C subset ell³C. For a finite D-module, Tor_D(ell^n C)=ell^n Tor_D C; this follows from any D-free plus torsion decomposition, without R0 compatibility.',
  'Taking torsion gives ell²S subset ell³T. These cyclic modules have lengths b-2 and e-3 (e>=3), so b-2<=e-3. Thus e>=b+1.',
  'Consequently e is b+1 or b+2. Both parameters are still unknown; in particular HF4=4 means b=2 and only constrains e to3 or4.'
 ]
 # Finite sanity receipts for the general invariant-factor proof, not its
 # replacement: direct ell-power dimensions and generating series formula.
 receipts=[]
 for b in range(2,9):
  hs=[3,6]+[4+int(j<4)+int(j<b+2) for j in range(2,b+5)]
  assert hs[2:4]==[6,6] and hs[4] in[4,5] and hs[-1]==4
  diffs=[hs[0]]+[hs[k]-hs[k-1] for k in range(1,len(hs))]
  expected=[0]*len(hs);expected[0]=3;expected[1]=3;expected[4]-=1;expected[b+2]-=1
  assert diffs==expected
  receipts.append({'b':b,'possible_e':[b+1,b+2],'finite_HF_prefix':hs})
 proof={'scope':'same actual full C, not C3 or original H0',
  'fresh_recursive_no_cut_NF4':[[str(q) for q in row] for row in N4],
  'full_HB0_order3_cancellation_without_cut_and_old_row0_3_projection':True,
  'actual_nonzero_full_quartic_coordinate':{'source':'HB0','row':0,'coefficient':'-3/5'},
  'independent_DVR_embedding_proof_of_HF4_4_or5':embeddingproof,
  'noncanonical_D_type_m2C':{'free_rank':4,'torsion_exponents':[2,'b>=2']},
  'all_higher_HF_formula':'h_j=4+1_{j<4}+1_{j<b+2} for every j>=2; h0=3,h1=6',
  'Hilbert_series':'(3+3t-t^4-t^(b+2))/(1-t)',
  'additional_e_b_refinement_proof':exponentproof,'remaining_e_cases':['b+1','b+2'],
  'finite_sanity_receipts_not_general_proof':receipts,
  'no_actual_b_or_e_value_and_no_quartic_completeness_claim':True,
  'unknown_full_quartic_coordinates_not_assigned_zero':True,
  'no_new_source_values_physical_gauge_IFT_or_jet':True}
 ep=out/'independent-full-HB-projection-DVR-embedding-and-parameterized-tail.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
  'code_sha256':sha(__file__),'evidence_sha256':{ep.name:sha(ep)},
  'full_HF4_4_or5_and_m2_type2_b_accepted':True,'additional_e_in_bplus1_bplus2_proved':True,
  'b_and_e_values_open':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
