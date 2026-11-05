"""ROOT whole-code/SOURCE review for the one fixed streamed minor projection."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
worker='code/s61_col53_next17_79_streamed_value_20261105_v64.py';wh='a946d904a20a78bc7f9a211ee010620ac0280035bce6565eb7b729ecbaffe808'
gate='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64/report.json';gh='c5823b8df2afda1ff7c4285889da3f45e06ace7d00d89fc27c72133464fa98e3'
assert sha(R/worker)==wh;pins[worker]=wh
g=load(gate,gh)
peer='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_GATE_20261105_v64_REVIEW/report.json'
p=load(peer,'012f7d578a5240351e40451479f4caeb22997b85247a9e2956de00f572cab27a')
assert p['input_sha256'][worker]==wh and p['input_sha256'][gate]==gh
load('reports/S61_COL53_NEXT17_79_STREAMED_VALUE_CUSTODY_20261105_v64/report.json','e75cf2b3576812d20022139e41f7096b701830e04703265f0255094dc9d6975f')
load('reports/S61_COL53_NEXT17_79_STREAMED_PHASE_COMPANION_20261105_v64/report.json','2577d4a15012a9ae99c806b93c55fe7f00c70d5c29ec25584eaa141d1aec94e9')
for f,h in g['input_sha256'].items():assert sha(R/f)==h;pins[f]=h
assert g['whole_numeric_forecast']==138535 and g['fresh_RAW_upper']==26155445
bpath=str(Path(gate).parent/'streamed-whole-lifetime-and-raw-budget.json');b=json.loads((R/bpath).read_text());assert sha(R/bpath)==g['evidence_sha256'][Path(bpath).name];pins[bpath]=sha(R/bpath)
assert b['global_numeric_guard']==150000 and b['hard_VALUE_bits']==2048 and b['temporary_guard_bits']==1921
assert b['ALL_generated_RAW_upper']==23968000 and b['one_generation_RAW_upper']==9268352
out=R/'reports/S61_COL53_NEXT17_79_STREAMED_ROOT_SOURCE_REVIEW_20261105_v64';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_code_and_custody_only_no_chosen_values=True,input_sha256=pins,
 manually_reviewed_complete_worker_and_helpers=True,independent_whole_peer_accepted=True,
 normal_success_owned_upper=138535,legitimate_failure_prefix_owned_upper=138658,corrected_child_load_upper=126982,anchor_stop_with_retained_children_upper=128892,global_guard=150000,
 fresh_metadata_exact=113478,fixed_metadata_base=118000,all1884_dependency_integer_weight_and_positive_partial_bounds_independently_verified=True,
 fresh_RAW_upper=26155445,all_generated_RAW_upper=23968000,two_complete_codec_generations_counted=True,
 code_audit=['fixed qi polynomial recurrence and positive integer qStored/qi restoration, including negative qi','actual empty mathG0 versus existing codec sentinelG1 flag0, no empty ratio division','all five coordinates and7 numeric charge per VALUE map term; SOURCE lead never recomputed from actual zero','complete1884 fixed graph/source129/rootdegree2; fullF11 remains the DAG input, old face only independent anchor check','sparse sum/copy/multiply integer weights, exact projection HL ties and relative S degree','canonical gcd/shift/content, both codec generations byteidentical, per-segment hash/index recovery','input/control/helper hashes pinned and independent parser/sign/empty/cursor fixtures accepted','metadata-only companion closes childdecode plus2048 and legitimate failurecodec/retainedchild charges below unchanged150k','all completed tempcells/stream/hashmanifest and partial accumulator retained on STOP; no automatic retry'],
 one_fixed_signed_projection_job_authorized=True,output_stem='S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64',
 no_new_beta_prime_or_source_or_matrix_RHS_update_authorized=True,all_restoration_P290_once_and_F9_F10_F11_complements_kept=True,
 independent_signed_postpeer_required=True,no_pivot_or_generic_exhaustion_result_yet=True,code_sha256=sha(Path(__file__)))
q=out/'report.json';q.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n');assert q.stat().st_size<=65536;print(sha(q))
