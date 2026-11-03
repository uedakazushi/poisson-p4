#!/usr/bin/env python3
"""Exact coefficient support/span checks and honest bounded-solver receipts."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--pencil-replay',type=Path,required=True)
 ap.add_argument('--attempts-replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 paths=[Path(__file__)];packets=[]
 for stem,replay in [('S61_AFFINE_QUADRATIC_PENCIL_AUDIT',args.pencil_replay),('S61_AFFINE_NEW_PROJECTION_ATTEMPTS',args.attempts_replay)]:
  directory=ROOT/f'reports/{stem}_20261103_v62';p=directory/'report.json';d=load(p);fresh=load(replay)
  ignore={'seconds','checked_at_utc'}
  assert {k:v for k,v in d.items()if k not in ignore}=={k:v for k,v in fresh.items()if k not in ignore}
  paths.append(p)
  for file,h in d['input_sha256'].items():q=ROOT/file;assert sha(q)==h,file;paths.append(q)
  for file,h in d['evidence_sha256'].items():q=directory/file;assert sha(q)==h,file;paths.append(q)
  packets.append((directory,d))
 directory,pencil=packets[0];p=32003
 with gzip.open(ROOT/'reports/S61_AFFINE_SMALL_PARAMETER_PROJECTION_20261103_v62/restricted-original134.json.gz','rt')as f:encoded=json.load(f)
 mons=list(itertools.combinations_with_replacement(range(45),2));index={m:i for i,m in enumerate(mons)};exact=defaultdict(Fraction)
 for row,data in enumerate(encoded):
  for mon,coefficients in data:
   if len(mon)!=2:continue
   for exp,c in coefficients:
    assert sum(exp)<=1;k=0 if not sum(exp)else exp.index(1)+1;exact[k,row,index[tuple(mon)]]+=Fraction(c)
 rows=[{i for(k,i,j),v in exact.items()if k==r and v}for r in range(6)];assert list(map(len,rows[1:]))==[64]*5
 arrays=np.load(directory/'coefficient-pencils-modp.npz');C=np.zeros((6,134,1035),dtype=np.int64)
 for(k,i,j),v in exact.items():assert v.denominator%p;C[k,i,j]=v.numerator*pow(v.denominator,-1,p)%p
 assert np.array_equal(C,arrays['C'])
 variable={j for(k,i,j),v in exact.items()if k and v};cols=[j for j in range(1035)if j not in variable];assert len(cols)==25
 span=load(directory/'constant-column-exact-span.json');A=s.Matrix(134,25,lambda i,j:s.Rational(exact[0,i,cols[j]]));basis=span['basis_columns'];B=A[:,basis];Q=s.Matrix(span['rational_span_coefficients'])
 assert B*Q==A and B[span['row_basis'],:].det()!=0 and len(basis)==16
 selected=arrays['pivot_columns'];inv=arrays['inverse'];assert np.array_equal(inv@C[0][:,selected]%p,np.eye(134,dtype=np.int64))
 N=np.array([inv@C[k][:,selected]%p for k in range(1,6)]);assert np.array_equal(N,arrays['N'])
 traces=[int(np.trace(M)%p)for M in N];assert traces==pencil['selected_minor_logarithmic_derivatives_mod32003']and all(traces)
 for _,d in packets:
  assert not d.get('affine_dimension_or_length_result_obtained',False)and not d.get('full_affine_dimension_or_length_result_obtained',False)
 attempts=packets[1][1]['attempts'];assert [a['exit_code']for a in attempts]==[124,14,124]
 assert all(not a['exact_unit_or_eliminant_obtained']for a in attempts)
 paths += [ROOT/'research/notes/S61_AFFINE_PREPROCESSING_BLOCKERS_20261103_v62_REVIEW.md',
  ROOT/'reports/S61_AFFINE_MONIC_DIMENSION_CRITERION_20261103_v62_REVIEW.json',
  ROOT/'reports/S61_AFFINE_FINITE_WRONSKI_CRITERION_20261103_v62_REVIEW.json']
 result=dict(id='S61_AFFINE_PREPROCESSING_BLOCKERS_20261103_v62_REVIEW',all_passed=True,
  fresh_subsecond_pencil_and_receipt_producers_match=True,all134_exact_QQ_quadratic_coefficients_reconstructed=True,
  exact_QQ_nonzero_row_support_counts_for_C1_through_C5=[64]*5,exact_parameter_free_columns=25,exact_QQ_constant_column_rank=16,
  selected_modular_inverse_and_each_N_trace_independently_checked=True,
  selected_minor_nonconstant_in_each_parameter_accepted=True,universal_constant_minor_impossibility_claimed=False,
  bounded_projection_failures_exit_codes=[124,14,124],q32_probe_memory_failure_honestly_retained=True,
  failed_solver_jobs_repeated=False,full_unit_generic_finiteness_or_length_conclusion_obtained=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS exact pencil limitations and bounded projection blockers; no unit claim')

if __name__=='__main__':main()
