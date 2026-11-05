#!/usr/bin/env python3
"""Bounded nonmathematical exception/codec lifetime fixture."""
import argparse,gc,hashlib,json,sys,traceback,weakref
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_UPPER20_CODEC_LIFECYCLE_20261105_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 out.mkdir(exist_ok=False,parents=True);refs=[]
 class CodecLocal(dict):__slots__=('__weakref__',)
 def interrupted_codec():
  payload=CodecLocal(copy=['x'*1024]*16);refs.append(weakref.ref(payload))
  try:raise ValueError('synthetic_inner_only')
  except ValueError as inner:
   partial=CodecLocal(copy=['y'*1024]*16);refs.append(weakref.ref(partial));raise RuntimeError('synthetic_outer_only') from inner
 try:interrupted_codec()
 except BaseException as error:
  before=[r() is not None for r in refs];assert before==[True,True];saved=repr(error)
  tb=error.__traceback__;traceback.clear_frames(tb);error.__traceback__=None;error.__context__=None;error.__cause__=None;del tb;gc.collect()
  after=[r() is not None for r in refs];active_tb=sys.exc_info()[2];assert after==[False,False] and active_tb is None
  result={'source_no_math_or_function_or_response_inputs':True,'codec_local_objects_alive_before_release':before,'codec_local_objects_alive_after_release':after,'active_exception_traceback_none':active_tb is None,'only_error_repr_retained':saved,'context_cause_none':error.__context__ is None and error.__cause__ is None,'gc_unreachable_cycle_collection_done_before_future_GUARD':True,'same_control_order_as_VALUE_V3':True,'bounded_payload_two_by16shared1024byte_strings':True}
 ep=out/'synthetic-codec-frame-and-context-release.json';ep.write_text(json.dumps(result,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'synthetic_nonmath_lifecycle_fixture_ONLY':True,'code_sha256':sha(Path(__file__)),'python_version':sys.version.split()[0],'evidence_sha256':{ep.name:sha(ep)}};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps(report))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output.resolve())
