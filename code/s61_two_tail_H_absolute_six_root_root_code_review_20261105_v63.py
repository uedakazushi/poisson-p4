"""ROOT immutable six-root worker CODE/custody review, no root/series values."""
from pathlib import Path
import json,hashlib
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
cpath='reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_CUSTODY_20261105_v63/report.json'
c=load(cpath,'2df082f7682db14be7b0cd7df73f75952265e551afb41a39b80dd3a8a36c2846')
w='code/s61_two_tail_H_absolute_six_root_worker_20261105_v63.py';wh='76dc04a1b1d24513307f1a0fb4895cfff98223ed50ee08190b7d4c6a8922787e'
assert sha(R/w)==c['worker_sha256']==wh;pins[w]=wh
for p,h in c['runtime_input_sha256'].items():assert sha(R/p)==h;pins[p]=h
load('reports/S61_TWO_TAIL_H_ABSOLUTE_HENSEL_ROOT_SOURCE_REVIEW_20261105_v63/report.json','523ac8adfaec65208092de73dcb724693b2b084f4c07df8a4062a29f51865113')
for name,h in c['evidence_sha256'].items():assert sha(R/Path(cpath).parent/name)==h
out=R/'reports/S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ROOT_CODE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,CODE_SOURCE_custody_only_no_factor_root_Hensel_values=True,input_sha256=pins,
 manually_reviewed_complete200line_worker=True,prior_independent_content_and_absolute_factor_theorem_accepted=True,
 code_review=['fixedP32003 trialdivision2..178 provesprime; sourcecontent andalpha/55 denominators areunits','seedmonicquartic exactdivision/product with monicquadratic, all4roots fromexhaustivealgebraicfieldenumeration; Legendre nonresidue provesquadraticfield','Fp2conjugation/norm/inversion formula and explicitunitcheck correct; Karatsuba sums reducebefore30bit products','sixsimple rootseries update coefficientn usingpreviouscoefficients andfixedderivativeinverse; fullrootjetschecked','linearfactors andall64subsets usefixeda-major/e<10 representation; allocate55+66e thencomplete sixproductverified','odd masks1..61 giveexact31unorderedproper subsets; twoLCallocations =>62FULLa<=6/e<=18products, firstmismatchwitnessretainedonlyafterentireproductcomputed','matchedpairmeansonlyfixedimagefactor; no match meansPRIMARYabsoluteirred requiringindependentpostpeer','everyfield product counts hard2M, prospectivepending+loop/scalar/poly/series state retained; 24kowner and4MiB streamingRAW/45soft60hard2GB unchanged'],
 hard_base_products=2000000,hard_owned_numeric=24000,hard_RAW_bytes=4194304,field_temp_bits=30,soft_seconds=45,hard_seconds=60,address_space_bytes=2000000000,
 one_fixed_value_job_authorized=True,output_stem='S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ENTRY_20261105_v63',
 independent_signed_identity_and_geometric_theorem_postpeer_required=True,no_extra_prime_point_factor_farm_or_retry=True,
 no_irreducibility_result_or_fullrank_Fitting_or_exhaustion_claim_yet=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p))
