"""Independent SOURCE review: row permutation/inverses, full worker, no ring values."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
R=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=R/'reports/S61_TWO_TAIL_H_CHANGED15_ROOT_SOURCE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
t=time.monotonic();pins={}
def load(s):
 p=R/s;pins[s]=sha(p);return json.loads(p.read_text())
c=load('reports/S61_TWO_TAIL_H_CHANGED15_CUSTODY_20261105_v63/report.json')
assert pins[next(iter(pins))]=='38216fcd90439f1efcdbf2acbc9e236f63092dce7c01e49c8c767b7f6e833755'
for s,h in c['input_sha256'].items():assert sha(R/s)==h;saved=pins.setdefault(s,h)
p='reports/S61_TWO_TAIL_H_CHANGED15_CUSTODY_20261105_v63/'
s=load(p+'genuine15-source-identity-and-frame.json');b=load(p+'source-and-whole-budget.json')
for n,h in c['evidence_sha256'].items():assert sha(R/p/n)==h
worker=R/'code/s61_two_tail_H_changed15_axis_20261105_v63.py';assert sha(worker)==c['worker_sha256']=='8e610f36d998118eaddad98d18ba33d3d6efa7054ad85e4694559fd3c4218d9a'
assert s['target_original32_columns']==[22,11,13,25,15,19,14,18,28,30,29,17,27,1,2]
assert s['genuine_current15_rows']==[0,1,3,6,2,11,12,9,7,13,10,14,4,5,8]
assert s['column_unit_product']=='B^2/L^4' and s['anchor_NOT_inverted']
assert s['Laplace_index_sum']==sum(range(12,15))*2==78 and s['Laplace_sign']==1
A=s['all15_row_operator'];B=s['all15_row_operator_inverse'];n=15
# Laurent polynomial multiplication in QQ[a,a^-1], with a separate dictionary
# accumulation, rather than the producer's permutation of a prior receipt.
def entry_prod(a,b):
 z={}
 for i,u in a.items():
  for j,v in b.items():
   k=int(i)+int(j);z[k]=z.get(k,Q(0))+Q(u)*Q(v)
 return {k:v for k,v in z.items()if v}
checks=0
for left,right in [(B,A),(A,B)]:
 for i in range(n):
  for j in range(n):
   z={}
   for k,ap in left[i].items():
    for e,v in entry_prod(ap,right[int(k)].get(str(j),{})).items():z[e]=z.get(e,Q(0))+v
   z={e:v for e,v in z.items()if v};assert z==({0:Q(1)}if i==j else {});checks+=1
rows=c['axis_input_polynomials'];assert len(rows)==16
assert [(x['row'],x['column'])for x in rows]==[(i,j)for i in range(4)for j in range(4)]
terms=sum(len(x['coefficients'])for x in rows);assert terms==85
assert min(a for x in rows for a,v in x['coefficients'])==0 and max(a for x in rows for a,v in x['coefficients'])==9
assert max(abs(int(v)).bit_length()for x in rows for a,v in x['coefficients'])==188
assert all(int(x['D'])%32003 for x in rows)
assert b['scalar_product_upper']==2288 and b['owned_prior']==2231
assert b['hard_products']==10000 and b['hard_owned']==4096 and b['field_temporary_bits']==30
r=dict(id=out.name,all_passed=True,SOURCE_only_no_determinant_or_ring_values=True,input_sha256=pins,
 all450_actual_Laurent_row_inverse_products=checks,all16_axis_input_cells_source_bound=True,ROM_bits=188,axis_terms=85,
 manually_reviewed_full_worker_sha256=sha(worker),
 independent_mathematical_argument='After A Ri-ci R_anchor for last three rows, old12 cells vanish; determinant scales A^3. Bottom block AV-cu has det A^2 detPair by universal 4x4 block identity. Both Laplace selected row and column sums are39, giving sign+. Cancel A^2 in original polynomial domain, never localize at A. Thus A D15new=detPair D12. All column and row normalization factors are U6 units; combined quotient nonzero must be tested, not separate factors.',
 code_review=['subset recurrence row=|mask|-1 and sign=row+column_position valid over any commutative ring','monic degree4 reduction convolution16+reduction12; no field/reduced assumption','only QQ denominator inverses modulo known prime; a inverse formula validated in ring','mandatory D12*detPair uses independently Newton-verified existing D12 remainder','input85 coefficients188bit distinguished from field temporary30bit','hard products10000/owned4096/RAW4MiB/45soft60hard2GB unchanged; complete failure prefix retained'],
 one_fixed_combined_value_job_authorized=True,output_stem='S61_TWO_TAIL_H_CHANGED15_AXIS_ENTRY_20261105_v63',
 no_other_prime_axis_minor_fallback_authorized=True,independent_signed_postpeer_required=True,
 no_H0_wholeH_Fitting_or_global_exhaustion_claim=True,code_sha256=sha(Path(__file__)),seconds=time.monotonic()-t)
q=out/'report.json';q.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n');print(sha(q))
