#!/usr/bin/env python3
"""Narrow immutable V2/V3 import repair; --help only, never main."""
import ast,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CPRIME_QUADRATIC20_IMPORT_CORRECTION_20261105_v63_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 base=ROOT/'reports/AMBIENT_B18_CPRIME_QUADRATIC20_IMPORT_CORRECTION_20261105_v63';p=base/'report.json';assert sha(p)=='8676ff5f578def96f2289bacb09a4e30037e670e4ba7bf1ced5cdb49e63b14cc';d=json.loads(p.read_text());pins[str(p.relative_to(ROOT))]=sha(p)
 for p,h in d['input_sha256'].items():assert sha(ROOT/p)==h;pins[p]=h
 for p,h in d['evidence_sha256'].items():assert sha(base/p)==h;pins[str((base/p).relative_to(ROOT))]=h
 v2=ROOT/'code/ambient_B18_Cprime_quadratic20_values_v2_20261105_v63.py';v3=ROOT/d['corrected_code'];assert sha(v3)==d['corrected_code_sha256']=='16e83e2bc7d50236eb65440465d76eff09388d17d5fb02bf7e1783dcc6541a14'
 before=v2.read_text();after=v3.read_text();expected=before.replace(' import fixed,sp,rec\n',' import fixed,sp\n').replace('def qq(v):','def rec(v):return {i:Q(c)for i,c in v if Q(c)}\ndef qq(v):').replace(" for p,h in report['input_sha256'].items():", " correction=read('reports/AMBIENT_B18_CPRIME_QUADRATIC20_IMPORT_CORRECTION_20261105_v63/report.json');assert correction['all_passed']and correction['corrected_code_sha256']==sha(Path(__file__))\n for p,h in report['input_sha256'].items():")
 assert expected==after
 fn=next(n for n in ast.parse(after).body if isinstance(n,ast.FunctionDef)and n.name=='rec');module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]));ns={'Q':F};exec(compile(module,'<pure-sparse-decoder-only>','exec'),ns)
 cases=[[],[[0,'0'],[2,'-7/3'],[19,'4/8']],[[1,'9'],[5,'-111/37'],[7,'0']]]
 assert all(ns['rec'](v)=={i:F(c)for i,c in v if F(c)}for v in cases)
 valueout=ROOT/'reports/AMBIENT_B18_CPRIME_QUADRATIC20_VALUES_20261105_v63';assert not valueout.exists()
 run=subprocess.run(['python',str(v3),'--help'],capture_output=True,text=True,timeout=15,cwd=ROOT);assert run.returncode==0 and not run.stderr and not valueout.exists();(OUT/'independent-help-stdout.txt').write_text(run.stdout);(OUT/'independent-help-stderr.txt').write_text(run.stderr)
 for p in [Path(__file__),v2,v3,ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id=ID,all_passed=True,frozen=True,narrow_local_decoder_import_repair_only=True,exact_V2_preserved=True,V3_sha256=sha(v3),all_mathematical_formulas_and_caps_unchanged=True,independent_import_only_help_pass=True,no_value_main_or_source_arithmetic_run=True,source_review_e7b7_remains_frozen=True,authorized_one_V3_value_gate_ready=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name!='report.json'},seconds=time.monotonic()-st)
 (OUT/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
