#!/usr/bin/env python3
"""Independent reversed old-C2 dots and fixed genuine-target E7/kernel checks."""
import hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
import ambient_B18_C_full_eta4_cache_values_review_20261105_v63 as ar
from ambient_B18_C3_alpha4_cached_source_review_20261105_v63 import nc
ROOT=pathlib.Path(__file__).resolve().parents[1]
PRI='AMBIENT_B18_C_GENUINE_T7_CACHED_VALUES_20261105_v63';ID=PRI+'_REVIEW'
PRIMARY='c93b6250e8681fb40647f668bf01312822b6f504c9fcd10680fce059cd606fec'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None):
  nonlocal RAW
  p=pathlib.Path(p);s=sha(p)
  if h:assert s==h,(str(p),s,h)
  key=str(p.relative_to(ROOT))
  if key not in pins:RAW+=p.stat().st_size
  pins[key]=s
  return json.loads(p.read_text()) if p.suffix=='.json' else None
 pr=pin(ROOT/'reports'/PRI/'report.json',PRIMARY);assert pr['all_passed']
 docs={p:pin(ROOT/p,h) for p,h in pr['input_sha256'].items()}
 peer=next(v for p,v in docs.items() if p.endswith('_REVIEW/report.json'))
 assert peer['D_only_primitive_split_genuine_T7_source_and_code_PASS']
 D=next(v for p,v in docs.items() if p.endswith('actual-filtered23-Dfree16-quotient7-and-genuine4-source.json'))
 C2=next(v for p,v in docs.items() if p.endswith('actual-physical12-and-C2-right-exact-quotient.json'))
 E=pin(ROOT/'reports'/PRI/'four-cached-coefficients-actual-E7-and-genuine-four-target-coordinates.json',
  '65703acebd05d111a30f051dbebb13a8085573ae11afeb5960d7421d0eb697d0')
 assert sha(ROOT/'code/ambient_B18_C_genuine_T7_cached_values_20261105_v63.py')==pr['code_sha256']
 assert sha(ROOT/'research/notes'/f'{PRI}.md')==pr['note_sha256']
 owned=sum(nc(x) for x in docs.values())+nc(E)+nc(pr)+12000
 assert owned<30000 and RAW<600000
 save(out/'pre-arithmetic-independent-four-dots-bound.json',{'owned_upper':owned,'RAW_read':RAW,'source_coefficient_bound43_pre87':True,'hard_bits':2048,'hard_value':30000,'hard_seconds':60,'hard_memory':2<<30})
 delta={r['physical_row']:dict(r['delta']) for r in C2['physical12_coefficients']}
 keys={tuple(k) for op in D['four_ell_coefficients_UNPAIRED'] for k,w in op['operator_old_C2_only']}
 assert len(keys)==24 and all(k[0]=='cache' for k in keys)
 ar.premax=ar.valmax=0;leaves={k:ar.parse(delta[k[1]].get(k[2],'0')) for k in keys};dots={};products={}
 for op in D['four_ell_coefficients_UNPAIRED']:
  value=Q(0)
  for k,w in reversed(op['operator_old_C2_only']):
   p=ar.mul(ar.parse(w),leaves[tuple(k)]);value=ar.add(value,p);products[(op['label'],tuple(k))]=p
  dots[op['label']]=value
 assert {k:str(v) for k,v in dots.items()}==E['all4_actual_cached_coefficients']
 assert products=={(label,tuple(k)):ar.parse(p) for label,k,w,q,p in E['all60_termwise_products']}
 assert len(products)==60 and dots=={'alpha_A':-Q(85,6),'alpha_B':Q(8315,432),'beta_A':Q(275,3),'beta_B':-Q(359975,1728)}
 a,b,c,d=[dots[k] for k in ['alpha_A','alpha_B','beta_A','beta_B']]
 M=[[Q(0)]*7 for _ in range(7)];M[2][0]=1;M[5][1]=a;M[6][1]=b;M[5][2]=c;M[6][2]=d;M[5][3]=1;M[6][4]=1
 K=[[Q(0)]*4 for _ in range(7)];K[1][0]=1;K[3][0]=-a;K[4][0]=-b;K[2][1]=1;K[3][1]=-c;K[4][1]=-d;K[5][2]=1;K[6][3]=1
 def mm(A,B):
  output=[]
  for row in A:
   nr=[]
   for j in range(len(B[0])):
    v=Q(0)
    for k in reversed(range(len(B))):v=ar.add(v,ar.mul(row[k],B[k][j]))
    nr.append(v)
   output.append(nr)
  return output
 assert [[str(q) for q in row] for row in M]==E['actual_E7']
 assert [[str(q) for q in row] for row in K]==E['actual_kernel7_by4']
 assert mm(M,K)==[[0]*4 for _ in range(7)]
 assert [[M[i][j] for j in [0,3,4]] for i in [2,5,6]]==[[int(i==j) for j in range(3)] for i in range(3)]
 assert [[K[i][j] for j in range(4)] for i in [1,2,5,6]]==[[int(i==j) for j in range(4)] for i in range(4)]
 M2=mm(M,M);assert c and M2[5][0]==c and M2[6][0]==d
 assert sum(bool(q) for row in M2 for q in row)==2
 assert [[str(q) for q in row] for row in M2]==E['actual_E2']
 assert mm(M2,M)==[[0]*7 for _ in range(7)] and ar.premax<=87
 proof={'all24_old_C2_leaf_bindings_and60_reversed_products':True,
  'four_actual_coefficients':{k:str(v) for k,v in dots.items()},
  'actual_E7_and_kernel7_by4_match':True,'all28_EK_products_and_fixed_I3_freeI4':True,
  'actual_E_ranks':[7,3,1,0],'beta_vector_literal_nonzero':True,'E3_zero':True,
  'actual_genuine_C_ell_dimension':4,
  'source_genuine_target_isomorphism_bound_to_D_only_primitiveF':True,
  'not_naive_R0_quotient_action_and_no_connecting_values':True,
  'no_actual_full_C4_or_formal_torsion_lifts_computed':True,
  'own_max_before_operation_bits':ar.premax,'own_max_stored_QQ_bits':ar.valmax,'owned_upper':owned,
  'no_new_source_physical_gauge_IFT_or_jet':True,'no_original_H0_nil_transfer':True}
 ep=out/'independent-four-old-C2-dots-and-actual-genuine-E7-kernel.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'all_four_coefficients_E7_and_genuine_kernel_verified':True,
  'R0_transport_and_connecting_values_separate':True,'elapsed_seconds':time.monotonic()-start}
 save(out/'report.json',report);signal.alarm(0)
 print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json'),'seconds':report['elapsed_seconds']}))
if __name__=='__main__':main()
