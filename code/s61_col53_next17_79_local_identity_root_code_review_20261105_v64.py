"""ROOT fresh read-only signed-postpeer CODE gate, no node arithmetic."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parents[1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();pins={}
def load(p,h):
 q=R/p;assert sha(q)==h;pins[p]=h;d=json.loads(q.read_text());assert d['all_passed'];return d
stem='reports/S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64_LOCAL_IDENTITY_POSTPEER_GATE/'
gate=stem+'report.json';gh='9eb3ed13ecab356c8001366e636d769e27def438a299aa8bd7500367e9034e95';g=load(gate,gh)
w='code/s61_col53_next17_79_local_identity_postpeer_20261105_v64.py';wh='492200c950fdcdf6e2d6894f3f5dcc6c41cb421934cce7c17803c3f9b06298c5'
assert sha(R/w)==g['code_sha256']==wh;pins[w]=wh
for p,h in g['input_sha256'].items():assert sha(R/p)==h;pins[p]=h
for p,h in g['evidence_sha256'].items():assert sha(R/stem/p)==h
b=json.loads((R/stem/'positive-local-identity-prior-and-custody.json').read_text());assert b['whole_numeric_prior']==143793 and b['whole_RAW_upper']==27389132
out=R/'reports/S61_COL53_NEXT17_79_LOCAL_IDENTITY_ROOT_CODE_REVIEW_20261105_v64';out.mkdir(exist_ok=False)
r=dict(id=out.name,all_passed=True,SOURCE_CODE_only_no_saved_node_arithmetic=True,input_sha256=pins,
 full233line_worker_manually_reviewed=True,
 code_audit=['own bounded canonical Q5parser reconstructs fixedqi integer coefficients includingnegativeqi/nonprimitivecontent/actualempty','both1884savedgenerations bytecompared and sourcepinnedtreehash verified','independenttrueleaf binding and rational/monomialtransforms match originalphysical inputs, no producerevaluator imported','all1884COPY/ADD/MUL localidentities reversedtermorder fromsavedchildren exactdictmatch; proofbytopologicalinduction certifiesroot','samefull5coords/fixedSOURCElead/degree2/sourceC129 andsourceM824/value7 retained','prospectiveownership includes currentoracle, twochildren, expectedmap, singleliteralinputtransient21000 andall1884reportrows inMETA','scalarproducts/rawqi restore/codec bounds1921 underhard2048, 143793<150000/RAW27389132<32MiB','allcompletedlocalproofrows andpartialexpectedmap preserved ifSTOP, no automaticretry/newsource/matrixupdate'],
 hard_numeric_guard=150000,whole_numeric_prior=143793,whole_RAW_prior=27389132,temporary_bits=1921,hard_VALUE_bits=2048,
 one_read_only_signed_postpeer_authorized=True,output_stem='S61_COL53_NEXT17_79_STREAMED_VALUE_20261105_v64_LOCAL_IDENTITY_REVIEW',
 no_new_coefficient_functional_or_wholeminor_rank_claim_authorized=True,code_sha256=sha(Path(__file__)))
p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');assert p.stat().st_size<=65536;print(sha(p))
