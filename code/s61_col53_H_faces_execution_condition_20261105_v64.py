"""Metadata receipt only: close ROOT's explicit conditional ONE two-face authorization."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_H_FACES_ROOT_EXECUTION_CONDITION_20261105_v64'
PINS={'reports/S61_COL53_H_FACES_SOURCE_20261105_v64/report.json':'5e66838fb845def468a9a1c883fd08a5836f21fdfd3a88de8f4d62109013d718','reports/S61_COL53_H_FACES_VALUE_GATE_20261105_v64/report.json':'fe261a87001601168e3f52eae66fe4bd5e1e80ddda810836de7b656736ac09fe','code/s61_col53_H_faces_values_20261105_v64.py':'d7e9a87cda9384a31f4899e887d2a3816e4530fa815d94786708b6abf0934357'}
PEER='reports/S61_COL53_H_FACES_SOURCE_CODE_20261105_v64_REVIEW/report.json'
PEER_SHA='7c6e258796e3a5cd5f4b6a5d612c56c2def02b4ee03fd587fbefa11a7ff13783'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
if __name__=='__main__':
 p=ROOT/PEER;assert sha(p)==PEER_SHA
 r=json.loads(p.read_text());assert r['all_passed']and r['frozen']and r['actual_signed_face_values_evaluated']is False
 for f,h in PINS.items():assert r['input_sha256'][f]==h and sha(ROOT/f)==h
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 doc={'id':ID,'all_passed':True,'metadata_only_not_independent_value_proof':True,'ROOT_review_and_conditional_authorization_received':True,'ROOT_instruction':'ONE twoHface VALUE only after Mark fresh SOURCE+CODE+MODEplan+OWN/RAW PASS pins exact5e668/fe261/d7e9; no actual coefficients before that condition; preserve atomic faces/failure prefixes; no automatic replay','ROOT_instruction_condition_now_fulfilled':True,'independent_Mark_report':{'file':PEER,'sha256':PEER_SHA,'all_passed':True},'H_faces_only':[-3,21],'fixed_row':6,'fixed_column':53,'whole_owned_upper':126262,'META':98171,'five_coordinate_term_charge':7,'fresh_RAW_upper':17661238,'generated_RAW_upper':2097152,'commit_bits':640,'temp_bits':1921,'stored_VALUE_bits':2048,'hard_seconds':60,'signal_seconds':55,'memory_bytes':2000000000,'no_other_H_or_entryzero_rank_inverse_claim':True,'input_sha256':{**PINS,PEER:PEER_SHA,'code/s61_col53_H_faces_execution_condition_20261105_v64.py':sha(__file__)}}
 b=(json.dumps(doc,separators=(',',':'))+'\n').encode();assert len(b)<=65536
 (out/'report.json').write_bytes(b);print('report_sha256',sha(out/'report.json'))
