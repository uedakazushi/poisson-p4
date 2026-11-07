"""ROOT R7 cached-four-root finite absolute-factor test SOURCE/code review."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
cp='reports/S61_TWO_TAIL_R7_ABSOLUTE_REUSE_CUSTODY_20261105_v63/report.json'
c=load(cp,'fb2936e1129d156686dcf34813753a7e31b8bb45f37259a13e1c082acdde1794')
w='code/s61_two_tail_R7_absolute_reuse_worker_20261105_v63.py';wh='e20b9b0d553231ecff2e26925bebb4c17c7d85a54086136cbf1dc970144e5b5f';assert sha(R/w)==c['worker_sha256']==wh;pins[w]=wh
for p,h in c['runtime_input_sha256'].items():assert sha(R/p)==h;pins[p]=h
f=json.loads((R/'reports/S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63/canonical-factor-candidates.json').read_text())
poly={tuple(k):int(v)for obj in f['factors']if obj['total_degree']==7 for k,v in obj['terms']}
assert {k:v for k,v in poly.items()if k[0]==4}=={(4,0,1):11250000,(4,0,2):-3515625,(4,1,2):9281250}
assert all(d==0 for a,e,d in poly if a==0)
P={e:v for (a,e,d),v in poly.items()if a==0};rem=sum(v*(-5)**e*6**(6-e)for e,v in P.items());assert rem==33802272 and rem%32003==7104
out=R/'reports/S61_TWO_TAIL_R7_ABSOLUTE_ROOT_CODE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_CODE_only_no_R7_root_or_candidate_values=True,input_sha256=pins,
 complete_worker_manually_reviewed=True,global_geometric_a_content_unit_checked=True,Fp_e_content_remainder=rem,Fp_e_content_remainder_mod32003=7104,
 theorem='R7aLC=140625*d*(80-25d+66ed), nonzero a0 independentd excludesbothd andprimitivebilinearL overC. d1 andmodp preservea4. Fourdistinctcachedroots inverifiedFp2 with L*product4=Rbar mod e7 giveuniqueHensel factors; anygeometricfactor hasleading1orL andedegree<=6, so7unorderedproperseedsubsets*2allocations=14FULLa4/e12products suffice forabsoluteirreducibility. No APIcompletenessassumption.',
 code_audit=['exactR7d1/fulloldQ35coefficientbinding+old e7..9zero, fourcachedjets only; no seed orHenselconstruction','fixedprime/Legendre verifyquadraticfield, foursimple roots fromdistinctconstants+unitLC55; all28R7rootresiduals andwholeLproduct checked','Fp2reducedKaratsuba30bit, all16subsetproducts precision7; oddmasks1..13 twoLCallocations cover14FULL65cellproducts','allcompletedroot/subset/currentcandidate/fieldpending/polywork/cursors savedifguardSTOP, no automaticretry','100kbaseproducts/12kowned/2MiB/45soft60hard2GB unchanged'],
 one_fixed_cached_R7_value_job_authorized=True,output_stem='S61_TWO_TAIL_R7_ABSOLUTE_REUSE_ENTRY_20261105_v63',
 independent_finite_identity_and_geometric_theorem_postpeer_required=True,no_new_prime_beta_source_root_construction_or_factor_farm=True,
 no_R7_irreducibility_result_or_badFitting_or_fullrank_claim_yet=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p))
