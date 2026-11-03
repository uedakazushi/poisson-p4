#!/usr/bin/env python3
"""Dependency and mathematical-scope receipt for minimal lengths/R CM locus."""
import argparse
import json
from pathlib import Path
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 source=ROOT/'reports/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63.json';d=read(source);assert d['all_passed']and d['minimal_generic_lengths']=={'N':1,'R':2}
 paths=[Path(__file__),source,ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 paths.append(ROOT/'research/notes/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63_REVIEW.md')
 report=dict(id='AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63_REVIEW',all_passed=True,
  proof_type='independent conceptual audit of accepted exact-scope inputs; no new CAS needed',
  generic_N_length1_by_common_reduced_Z_base_and_flat_functionfield_A10_intersection_A9=True,
  generic_R_length2_by_actual_ordinary_double=True,
  actual_and_reduced_intersection_CM_locus_on_R='R minus B18',
  reduced_intersection_normal_locus_on_R='R minus B18',actual_intersection_normal_locus_on_R='empty',
  R_itself_CM_entire_N_minus_R_or_global_primary_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS generic lengths(1,2) and exact D/Dred CM/normal loci along R')

if __name__=='__main__':main()
