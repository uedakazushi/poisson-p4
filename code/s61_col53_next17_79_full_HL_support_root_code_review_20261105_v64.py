"""ROOT complete leading-HL SOURCE-support code review; no support convolution."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
stem='reports/S61_COL53_NEXT17_79_FULL_HL_SUPPORT_GATE_20261105_v64/';gate=stem+'report.json'
g=load(gate,'2171851bb814c112b93ab35950be9e82f92f05812a5d5dde26563fe4141699da')
load('reports/S61_COL53_NEXT17_79_FULL_HL_SUPPORT_CUSTODY_20261105_v64/report.json','a12114ca0c3b50b092388ccf952ef3022fdd8ed480a2e8102848038235025588')
w='code/s61_col53_next17_79_full_HL_support_20261105_v64.py';wh='beed82ad471b6855d452e4cfb4f3a3d1fc890605af7c8d154e0b9d0af9e82c61';assert sha(R/w)==wh;pins[w]=wh
for p,h in g['input_sha256'].items():assert sha(R/p)==h;pins[p]=h
for p,h in g['evidence_sha256'].items():assert sha(R/stem/p)==h
b=json.loads((R/stem/'whole-SOURCE-support-stream-budget.json').read_text());assert b['worker_sha256']==wh and b['whole_SOURCE_upper']==138331 and b['fresh_RAW_upper']==5468713
assert b['explicit_scalar_inventory']==110983
out=R/'reports/S61_COL53_NEXT17_79_FULL_HL_SUPPORT_ROOT_CODE_REVIEW_20261105_v64';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_code_review_only_no_support_convolutions_or_signed_values=True,input_sha256=pins,
 manually_reviewed_complete127line_worker_and_accepted_packed5_leaf_helpers=True,
 independent_metadata_correction='Conservatively add3768 parsed numeric mode keys to explicit110983:114751 remains below fixed118000, with3249additional slack;4096source control reserve already present. Frozen worker/gate unchanged.',
 conservative_metadata_bound=114751,fixed_metadata_base=118000,
 code_review=['same1884fullSLPnodes; completeleadingHLface noS cutoff, rootS7..38/t19 remainallfull5coordinates','SOURCEonly packed40 injectivekeys; eachproduct unpacksbothkeys/addsall5coordinates/shift thenpackrange+roundtrip beforeinsert; futureVALUE7 unchanged','streamedexactpositive setunions/Minkowskiproducts withonlytwochildsets+output, no coefficient evaluation','per-nodepositivecardUPPER asserted, beforeinsert prospectiveowner150k guard andleafdecode16152 charged','children explicitlyreleasedbeforeoutputsort/atomicstreamsegment;thinSHAindex/prefix maintained','originalsignedlowbandcells are untouched; newSOURCEsets needed because completeHLface supports were never computed','old C<=2signedcells remainreusable forfuturehigherband extension, no permissionfornewVALUE','allcompletedsegments/partialcurrentset/cursors retainedonSTOP, noautomaticretry'],
 SOURCE_numeric_upper=138331,global_guard=150000,fresh_RAW_upper=5468713,pair_count_upper=763161,
 one_SOURCE_support_profile_authorized=True,output_stem='S61_COL53_NEXT17_79_FULL_HL_SUPPORT_20261105_v64',
 signed_VALUE_authorized=False,no_face_zero_nonzero_or_newpivot_claim=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');assert p.stat().st_size<=65536;print(sha(p))
