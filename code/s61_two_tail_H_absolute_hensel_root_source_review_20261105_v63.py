"""Independent coefficient-content and finite absolute-factor SOURCE review."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;return json.loads(q.read_text())
s=load('reports/S61_TWO_TAIL_H_ABSOLUTE_HENSEL_SOURCE_20261105_v63/report.json','6b1ebfdb365ee0b24e6ebcd14e07b4e4a4e7e95feeecf28c6e8ed5c3b1fa4c99');assert s['all_passed']
p=R/'reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/factor-sparse.txt';pins[str(p.relative_to(R))]=sha(p)
assert pins[str(p.relative_to(R))]=='241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c'
lines=p.read_text().splitlines();assert int(lines[0])==len(lines)-1==94
F={}
for line in lines[1:]:
 a,e,d,c=map(int,line.split());F[a,e]=F.get((a,e),0)+c
F={k:v for k,v in F.items()if v};assert len(F)==42
assert max(a for a,e in F)==6 and max(e for a,e in F)==9 and max(a+e for a,e in F)==9
assert {e:c for (a,e),c in F.items()if a==6}=={0:17578125*55,1:17578125*66}
# Remainder of 6^9 P under 6e+5, using integer Horner evaluation,
# independently of producer long division. This is coefficient algebra only.
p0={e:c for (a,e),c in F.items()if a==0};rem=sum(c*(-5)**e*6**(9-e)for e,c in p0.items())
assert rem==-1774213652736 and rem%32003==24204
assert 17578125%32003 and 55%32003 and 66%32003 and rem%32003
out=R/'reports/S61_TWO_TAIL_H_ABSOLUTE_HENSEL_ROOT_SOURCE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_theorem_and_content_only_no_roots_Hensel_or_candidate_values=True,input_sha256=pins,
 d1_terms42_checked=True,a_degree6_preserved=True,e_degree9=True,formal_content_remainder=rem,formal_content_remainder_mod32003=24204,
 finite_proof_review=['aLC linear L=55+66e; content must divide L; nonzero remainder gives a-primitivity over algebraic closure','squarefree seed with4Fp linear plus1irreducible quadratic splits into6distinct roots in singleFp2','simple-root uniqueness lifts all6 roots through e9 in Fp2[[e]]; full six-root product must reconstruct Fbar/LC mod e10','any geometric factor has positive a-degree, leading coefficient scalar*1 or scalar*L, and e-degree<=9','after LC normalization every factor is the unique seed-root subset product, so coefficients belong to Fp2','31proper subsets containing seed0 times2LC allocations cover every unordered factorization; complete e<=18,a<=6 products decide62 candidates','absolute modp irreducibility plus a-primitivity and preserved a-degree6 excludes char0 geometric factors via Gauss valuation normalization','globalH geometric a-content excludes d and primitive bilinear Lglobal by accepted explicit witnesses; preserved d1 a-degree6 then excludes factors ofH'],
 implementation_requirements=['fresh immutable complete worker/custody plus full code audit before ONE GO','Fp2 Karatsuba sums reduced moduloP before multiplication to preserve30bit temporary bound','all root enumeration/seed/Legendre/inverse/normalize/series/subsets/full62candidateproducts counted in hard2million baseproducts','only onecandidate G/H/fullproduct resident, all previous candidates represented by mismatch witness/hash','all roots/subsets/currentcandidate/pending scalar/cursors/failure prefixes preserved under24000owned/4MiB/45soft60hard2GB'],
 proposed_worker_preparation_authorized=True,factor_values_NOT_authorized=True,no_irreducibility_result_yet=True,no_new_prime_or_beta_physical_point=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p))
