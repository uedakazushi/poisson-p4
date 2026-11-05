#!/usr/bin/env python3
"""Full-nine Vfirst SOURCE/frame and bounded known-closure method review.
No actual Vfirst, rank, new source/gauge or higher coefficient is evaluated.
"""
import hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(tuple,list)):return [enc(x)for x in v]
 return v
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={}
 def load(s,n='report.json'):
  p=ROOT/'reports'/s/n;pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),indent=2)+'\n')
 stem='AMBIENT_B18_CPRIME_LINEAR_INITIAL_SOURCE_PREFLIGHT_20261105_v63';r=load(stem);assert pins[f'reports/{stem}/report.json']=='00ab656d92a89bfaf434e117de770403beb2a127f0ff73e20a28a0642c5f8925'and r['all_passed']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/stem/n)==h
 P=load(stem,'full-nine-first-source-and-linear-initial-completion-gate.json');A=load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63','exact-existing-raw33-multiplication-table.json');T=load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-source105-ell2-raw33-values.json');K=load('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63','phase3-actual-kappa9-and-normalized-nilframe.json');C5=load('AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63','actual4variable7core-polynomials.json');R=load('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63','actual-radical-normalizer-comparison-operator.json')
 assert load('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW')['actual_Cprime2_length']==7
 assert load('AMBIENT_B18_RADICAL_COMPARISON_OPERATOR_20261104_v63_REVIEW')['all_passed']
 assert load('AMBIENT_B18_REDUCED_COMPARISON_GENERIC_RANK2_20261104_v63_REVIEW')['all_passed']
 basis=A['basis33'];sel=[i for i,d in enumerate(A['basis_degrees'])if d==1];assert len(sel)==9
 mapping={0:0}
 for i in sel:
  z,rr,e=basis[i];var=z-1 if z else 7 if rr else 8;mapping[i]=var+1
 assert set(mapping.values())==set(range(10))
 # The raw33 basis reverses z labels. Use the EXPLICIT canonical-variable
 # permutation before forming any native/negative-U column.
 assert [mapping[i]for i in sel]==[7,6,5,4,3,2,1,8,9]
 mp=[{mapping[i]:F(1)}if i in mapping else{}for i in range(33)]
 def put(v,k,c):
  v[k]=v.get(k,F(0))+c
  if not v[k]:del v[k]
 def times(a,b):
  z={}
  for i,c in a.items():
   for j,d in b.items():
    if not i or not j:put(z,i+j,c*d)
  return z
 for j in reversed(range(33)):
  for i in reversed(range(33)):
   v={}
   for k,c in A['table_columns_for_each_left_basis'][i][j]:
    for n,d in mp[k].items():put(v,n,F(c)*d)
   assert v==times(mp[i],mp[j])
 assert P['actual_first_quotient_basis33_indices']==[0]+sel and P['actual_first_quotient_basis_labels']==[basis[i]for i in[0]+sel]
 assert P['actual_R0_coordinate_order9']==['zA0','zA1','zA2','zA3','zA4','zA5','zA6','rA','ell']
 assert len(K['kappa9_C132'])==9 and all(len(v)==132 and F(v[0])==0 for v in K['kappa9_C132'])
 assert C5['selected5_pivot_order']==['v4','v3','v5','v7','v8']and C5['source9_variable_order']==['a','d','v2','v3','v4','v5','v6','v7','v8']
 for k in C5['selected5_pivot_order']:
  p=C5['five_actual_rational_roots'][k];assert isinstance(p['D_power'],int)and isinstance(p['P_power'],int)and p['D_power']>=0 and p['P_power']>=0
 assert len(R['R9_generators_polynomials9'])==9 and R['actual_C_variable_rules']['GL4_variable_input_order']==['a','d','x','t']
 # Check the genuine gA-first input rather than silently setting it0.
 roots={i:{j:F(c)for j,c in v}for i,v in T['all_requested164_root_values']};assert 7790 in roots
 gfirst={mapping[i]:c for i,c in roots[7790].items()if i in mapping and c};assert not gfirst
 # In base m²=0, trace maximal ideal^5=0. Therefore full proper rational
 # roots inC40 are evaluated exactly by finite unit inverses; τ constants
 # MUST be retained before the base-linear trace coefficient is read.
 maxdeg=max(t+u+v for t in range(2)for u in range(5)for v in range(5)if t+u//2+v//2<=1);assert maxdeg==4
 m2=list(itertools.combinations_with_replacement(range(9),2));m3=list(itertools.combinations_with_replacement(range(9),3));assert len(m2)==45 and len(m3)==165
 assert P['known_quadratic_closure_shape']==[135,396] and 44*9==396 and 44*45==1980
 # A canonical quadratic divisor leaves<=h2 entries for each of495
 # cubic monomials. ALL consistency I-qE columns are retained.
 for h2 in range(4):
  assert 9*h2<=27
  cubic_live=135*h2+495*h2+27*27+2*27*27+27*27+54
  assert cubic_live<=4860<5000
 # Audit the repaired, still UNRUN worker only for its storage contract.
 worker=ROOT/'code/ambient_B18_Cprime_linear_initial_values_20261105_v63.py';wh=sha(worker);assert wh=='1be2cda944e1dfc0a608e0ac03035168438f4240abc279ad34739146c474b028';text=worker.read_text()
 assert 'def cubic_stream():'in text and 'del qcols,iv'in text and 'del piv3'in text and 'reduced=[]'not in text
 assert text.count('in cubic_stream()')==1 and 'in enumerate(cubic_stream())'in text
 assert "for f,k in enumerate(free):"in text and "yield ['consistency',f,var]"in text and 'all_compressed_columns_never_materialized'in text
 budget=P['budget'];assert sum(budget.values())==85555 and budget['future_q_le3_and_cubic_stream']==5000 and P['whole_saved_live_QQ_upper']==85555
 gap=ROOT/'work/marked-square-v63/linear-initial-source-stream-gap/gap.json';oldgap=json.loads(gap.read_text());assert oldgap['possible_retained_QQ_upper']==54270 and oldgap['worker_sha256']!=wh
 proof=dict(all1089_fullnine_first_quotient_products=True,canonical_variable_map_from_raw33=mapping,actual_kappa_full_C40_required=True,actual_core_order_and5residual_order=C5['selected5_pivot_order'],known_unit_D_P_residues=[1,5],trace_maximal_ideal_power5_zero=True,gA_first_exactly_zero_from_existing_source=True,
  complete_first_source='FullVgamma is the actual3x44 presentation and haszero origin because Cprime needs3 minimal generators (accepted Cprimebar). Thus its44 actual fullnine first columns exhaust degree1 initial relations. Compute fullproper polynomials/rational roots andtrace carries BEFORE coefficient extraction.',
  positive_upper_proof='Knownquadratics=L1*S1 are only a submodule of trueinitial relations. If h2<=3, use canonical quadraticdivisor q,E with Id-Eq inknowncubic=L1*S2. QuotientF495/K3 is exactly F_(9h2)/(qK3+(I-qE)). Include all3ell³targets. If the known third quotient is0, actualm3=m4 andNakayama killsm3; lengthCprime/ell³ <=3+h1+h2.',
  lower_and_completion='Actual DVR rank2 withtorsion{1,e>=2} gives lengthCprime/ell³=6+1+min(e,3)>=9. If the proved known upper is9, then length9 ande=2. Failed or larger upper is inconclusive, not an exponent negative.',
  fixed_storage_gap_preserved=oldgap,repaired_worker_sha256=wh,real_streaming_contract='At most27 selectedoriginal andnormalized compressed columns; hash bothpasses, delete quadratic135 scratch andnormalizedbasis before small inverse, no full2010columns. Quadraticdivisor405+canonicalq3_1485+basis729+aug1458+inverse729+scratch54=4860<5000.',
  whole_saved_live_QQ_upper=85555,source_gate_passes_after_stream_repair=True,no_actual_Vfirst_rank_or_exponent_value=True,
  mandatory_runtime_checks=['all44_Vgamma_origins0 andactualninevariablelabels','allunitinverse andd-coordinate identities','every fixedminor inverse andallfactorization hashes','allI-qE consistency columns andall3ell³targets','if third quotient0 assert knownupper>=9; hard60s/2GB/32MiB andcoefficient/liveguard'] )
 save('independent-fullnine-frame-known-initial-proof-and-stream-repair.json',proof)
 for p in[Path(__file__),worker,gap,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 rep=dict(id=ID,all_passed=True,frozen=True,fullnine_C40_source_frame_verified=True,conditional_cubic_stream_storage_repair_verified=True,whole_QQ_upper=85555,no_actual_first_rank_or_new_source_values=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st)
 assert sum(p.stat().st_size for p in OUT.iterdir())<32<<20;save('report.json',rep);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
