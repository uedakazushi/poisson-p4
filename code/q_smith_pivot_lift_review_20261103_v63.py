#!/usr/bin/env python3
"""Exact all75 nonlinear Smith row-combination ideal-membership audit."""
import argparse
import json
from pathlib import Path
import sympy as s
from q_smith_global_model_review_20261103_v63 import ROOT,read,source_polynomials,plus,scaled,scalar,sha

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63';frames=read(folder/'frames.json.gz');source=read(folder/'source109.json.gz')
 values=source_polynomials(source['nodes']);equations=[values[i]for i in source['all350_residual_nodes']];count=0
 for block in frames['Jacobi_complement_Smith_blocks']:
  U=block['left'];rows=block['rows']
  for i,text in enumerate(block['diagonal']):
   if not s.sympify(text):continue
   record=next(p for p in frames['source_Smith_pivots']if p['column']==block['columns'][i])
   expected=plus(*(scaled(equations[r],scalar(U[i][j]))for j,r in enumerate(rows)if s.sympify(U[i][j])))
   assert expected==values[record['source_node']];count+=1
 assert count==75
 paths=[Path(__file__),ROOT/'code/q_smith_global_model_review_20261103_v63.py',folder/'frames.json.gz',folder/'source109.json.gz',ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63_REVIEW.json']
 report=dict(id='Q_SMITH_PIVOT_LIFT_20261103_v63_REVIEW',all_passed=True,
  all75_nonlinear_transformed_pivots_equal_exact_U_linear_combinations_of_original350=True,
  coefficient_ring='QQ[a,w108]',no_samples_or_truncation=True,
  complements_frozen_global_model_review=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths})
 with args.output.open('x')as f:json.dump(report,f,indent=2);f.write('\n')
 print('PASS all75 nonlinear Smith pivot lifts in original350 ideal')

if __name__=='__main__':main()
