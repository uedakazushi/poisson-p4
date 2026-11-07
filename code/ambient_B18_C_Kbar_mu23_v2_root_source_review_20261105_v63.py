"""ROOT complete muK2/3 cache-only SOURCE/V2 control review, no scalar pairing."""
from pathlib import Path
import json,hashlib
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
peer='reports/AMBIENT_B18_C_KBAR_MU23_SOURCE_CODE_20261105_v63_REVIEW/report.json'
p=load(peer,'ecc7c58ac4a9ca1e131c2905ef4cb91a0f0f233738361bae15d3a1d20ed4112d')
source='reports/AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63/report.json'
s=load(source,'332217ddd9dbdde76278f4b9d5e3f44240b69f86cab3821c75a91a699009fd7e')
custody='reports/AMBIENT_B18_C_KBAR_MU23_VALUE_CONTROL_V2_CUSTODY_20261105_v63/report.json'
c=load(custody,'c8aadfb38eac43bd8a7391ee35267fb3a8726ce339313fecbb2d67dc8addafcc')
w='code/ambient_B18_C_Kbar_mu23_cached_values_v2_20261105_v63.py';wh='94d50c2461dc1f9d77e62ccbf547acd9c97f0bcd93b3bc829a2829e4ab77b7e7';assert sha(R/w)==wh;pins[w]=wh
compiled='reports/AMBIENT_B18_C_KBAR_MU23_CACHE_SOURCE_20261105_v63/complete-old-Pg-second-cache-source-and-fourteen-action-positive-gate.json'
assert sha(R/compiled)==p['input_sha256'][compiled]=='317b2c5e491e7f8369e73a7578982731415ac45d77444219308dcc5cc34b0127';pins[compiled]=sha(R/compiled)
for path,h in p['input_sha256'].items():assert sha(R/path)==h;pins[path]=h
out=R/'reports/AMBIENT_B18_C_KBAR_MU23_V2_ROOT_SOURCE_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_CODE_only_no_rr_t37_or14action_pairing=True,input_sha256=pins,
 root_full_original_math_worker_and_all_V2_control_diffs_reviewed=True,
 independent_SOURCE_domain_review_accepted=True,
 scope_checks=['originalKbar/rKbar hasbasisU,S,q2; all14ziU/S projectedtoq2 plusknownziq2 commutation exhaust maximalidealimage','37actualnormalizedR37 lifts/original768weights andliteralpiq2 index179 checked independently; noM3','pure-r firstSchurYs1head4all84 andXp1raw55all27zero leavesRAWPg_r2 only','21oldADphase3 scalarcacheinputs allopaqueSECONDgradients empty; no newIFT/gauge/source derivative','independent two-rowz5U/S F2serialization correction restoresnonzero q1inrKbar, q2projection0/compiledweights unchanged','37mandatory origin/rreadback then37rrdots and14same37/direct21rootfolds, lambdaS² witness minor checked','if anyactionnonzero originalmuR0K=2 byNakayama withknownlower2; ifall14zero fullmaximalidealquotientdim3 givesmu3','mu_DK7 andoriginalH0/nilunchanged; nofullformalKlifts constructed'],
 control_checks=['allSCIlists initialized andaliaslive scientificstate; perrow/dot/pendingparse+lookup+add+mul operands/computedcandidate captured','caughtSIGALRM/fullSTOPprefix withpartials/lower/t/actions/terms/cache; fixtureconfirmsfailurecandidate preserved','independentfailureupper10202<=static14076; privatehard60000/RAW800k unchanged','mandatoryexternalAS2GB andwall60 runner, noautomaticretry/capraise'],
 committed_prior=136,temporary_prior=273,static_numeric_prior=14076,private_hard_numeric=60000,RAW_cap=800000,
 one_fixed_cache37plus14_value_job_authorized=True,output_stem='AMBIENT_B18_C_KBAR_MU23_CACHED_VALUES_20261105_v63',
 independent_signed_postpeer_required=True,no_extra_Gamma12_dots_or_new_physical_derivatives_authorized=True,
 no_mu_result_H0_or_nil_conclusion_yet=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');print(sha(p))
