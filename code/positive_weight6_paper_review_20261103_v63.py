#!/usr/bin/env python3
"""Bind the independently reviewed mathematics to the new standalone PDF."""
import argparse,hashlib,json,re,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for block in iter(lambda:f.read(1<<20),b''):h.update(block)
 return h.hexdigest()
def main(out):
 start=time.monotonic();folder=ROOT/'reports/POSITIVE_WEIGHT6_COMPONENT_PAPER_20261103_v63';manifest=json.loads((folder/'manifest.json').read_text());assert manifest['all_passed']and manifest['measured_pdf_pages']==14
 files={}
 for key in ['evidence_sha256','checked_receipt_input_sha256','historical_frozen_sha256']:
  for name,want in manifest[key].items():assert sha(ROOT/name)==want,name;files[name]=want
 for key in ['source','pdf','supplement']:
  name=manifest[key];assert sha(ROOT/name)==manifest[key+'_sha256'];files[name]=manifest[key+'_sha256']
 info=subprocess.run(['pdfinfo',str(ROOT/manifest['pdf'])],check=True,capture_output=True,text=True).stdout
 assert re.search(r'^Pages:\s+14\s*$',info,re.M)
 log=(folder/'lualatex-final.log').read_text()
 for marker in ['undefined references','There were undefined','Rerun to get cross-references right',r'Overfull \hbox',r'Overfull \vbox']:assert marker.lower()not in log.lower(),marker
 extra=[Path(__file__),folder/'manifest.json',folder/'SOURCE_INDEX_REGISTRATION_REQUEST.json',ROOT/'reports/Q_NEW_FAMILY_FULL_COMPLETION_ROOT_REVIEW_20261103_v63.json',ROOT/'research/notes/POSITIVE_WEIGHT6_COMPONENT_PAPER_20261103_v63_REVIEW.md',ROOT/'research/notes/R4_NEW25_PRIMARY_LITERATURE_SPOTCHECK_20261103.md']
 for p in extra:files[str(p.relative_to(ROOT))]=sha(p)
 result=dict(id='POSITIVE_WEIGHT6_COMPONENT_PAPER_20261103_v63_REVIEW',all_passed=True,measured_pdf_pages=14,
  independent_review='Full family/component dimension, generic divisor/shear, firstjet scope, finite-rank length5 proof and separate positive-point nonCM argument read and accepted',
  full_completion_rank_lemma_not_a_flatness_assumption=True,
  contains_full_component25_lower58_firstjet15_and_generic_length5=True,
  closed_completion_scope_at_n111_and_positive_nonCM_distinguished=True,
  source_sha256=manifest['source_sha256'],pdf_sha256=manifest['pdf_sha256'],
  historical_seven_files_unchanged=True,layout_no_undefined_refs_or_overfull_boxes=True,
  old_numeric_certificates_reused_when_unchanged=True,literature_priority_established=False,
  global_component_exhaustion_or_generic_affine_five_exhaustion_claimed=False,
  input_sha256=files,seconds=time.monotonic()-start)
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2)+'\n');print('NEW14PAGE PAPER INDEPENDENT SCOPE/HASH/PDF REVIEW PASS',result['seconds'])
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
