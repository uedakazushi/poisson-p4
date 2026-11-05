#!/usr/bin/env python3
"""Strict inspection ONLY of exactK92AB001; never resumes arithmetic."""
import argparse,hashlib,json
from pathlib import Path
def main(p):
 d=json.loads(p.read_text());assert d['version']=='K92AB001'and d['no_rank_or_determinant_result_claim'];assert len(d['matrix'])==121 and len(d['five_buffers'])==5
 terms=0
 for arr,cap in[(d['matrix'],89),(d['five_buffers'],177)]:
  for pp in arr:
   old=-1
   for g,v in pp:assert old<g<cap and int(v)and abs(int(v)).bit_length()<=1738;old=g;terms+=1
 assert sorted(d['row_permutation'])==sorted(d['column_permutation'])==list(range(11))and d['sign']in[-1,1]
 assert d['next_operation_NOT_committed']and abs(int(d['current_quotient_coefficient'])).bit_length()<=864
 print(json.dumps(dict(state_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),version='K92AB001',all_exactintegers_parse=True,coefficient_terms=terms,stage=d['stage'],phase=d['phase'],cursor=[d['row'],d['column'],d['multiply_next_i'],d['multiply_next_j'],d['division_shift'],d['division_next_j']],NO_resume_or_rank_claim=True),indent=2))
if __name__=='__main__':
 a=argparse.ArgumentParser(description=__doc__);a.add_argument('--state',type=Path,required=True);main(a.parse_args().state)
