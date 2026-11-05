"""SOURCE stop byte custody and streaming numeric census; never replay compiler."""
from pathlib import Path
import hashlib,json,re
ROOT=Path(__file__).resolve().parents[1]
ID='S61_COL53_NEXT17_79_COMPILE_STOP_CUSTODY_20261105_v64'
BASE='reports/S61_COL53_NEXT17_79_COMPILED_SOURCE_20261105_v64'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(d):return(json.dumps(d,separators=(',',':'))+'\n').encode()
if __name__=='__main__':
 b=ROOT/BASE;failure=json.loads((b/'failure.json').read_text())
 assert failure['all_passed']is False and '157050'in failure['error']
 p=b/'partial-compiled-source.json';text=p.read_text()
 token=re.compile(r'"(?:\\.|[^"\\])*"|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?|true|false|null')
 integer_string=re.compile(r'[+-]?[0-9]+(?:/[0-9]+)?')
 count=0;tokens=0
 for m in token.finditer(text):
  t=m.group();tokens+=1
  if t[0]=='"':
   if integer_string.fullmatch(t[1:-1]):count+=3
  elif t not in['true','false','null']:count+=1
 # No json.loads of the entire prefix is used in this SOURCE census.
 files=[p,b/'failure.json',ROOT/'code/s61_col53_next17_79_compile_source_20261105_v64.py',ROOT/'research/notes/S61_COL53_NEXT17_79_COMPILED_SOURCE_20261105_v64.md',Path(__file__)]
 pins={str(f.relative_to(ROOT)):sha(f)for f in files}
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 d={'id':ID,'all_passed':True,'STOP_preserved_not_compiler_completion':True,'SOURCE_only_no_support_signed_coefficients':True,'actual_stop':failure,'full_prefix_bytes':p.stat().st_size,'stream_numeric_token_census_upper':count,'stream_tokens_seen':tokens,'prospective_fresh_prefix_parse_plus25k_runtime':count+25000,'counts_keys_numerical_strings_too_so_conservative':True,'no_full_prefix_parse_or_Source_arithmetic_replay':True,'immutable_executed_compiler_file':str(files[2].relative_to(ROOT)),'input_sha256':pins}
 (out/'report.json').write_bytes(enc(d));print(json.dumps({k:v for k,v in d.items()if k not in['input_sha256','actual_stop']}));print('report_sha256',sha(out/'report.json'))
