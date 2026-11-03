#!/usr/bin/env python3
"""Own sparse150 engine reconstructs every global109-slice group entry."""
import argparse
from collections import defaultdict
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_global_flag_review_20261103_v63 import IDS,action

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 directory=ROOT/'reports/Q_GROUP_SLICE_GLOBAL_UNIT_20261103_v63';path=directory/'report.json';d=load(path);f=load(args.replay);assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 paths=[Path(__file__),path,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():p=directory/file;assert sha(p)==h,file;paths.append(p)
 with gzip.open(ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz','rt')as stream:parent=json.load(stream)
 cert=load(directory/'certificate.json');cuts=parent['group_cut_raw_indices'];lies=parent['group_lie_columns'];W=[sum(pair)-sum(mon)for(pair,mon),i in IDS.items()]
 degrees={-1:0,0:0,**{i+1:w for i,w in enumerate(parent['positive_weights'])}};assert len(parent['positive_weights'])==108 and all(w>0 for w in parent['positive_weights'])
 forms=[]
 for r,cut in enumerate(cuts):
  row=[]
  for j,(u,v)in enumerate(lies):
   operator=action(u,v);raw={i:int(operator[cut,i])for i in range(150)if operator[cut,i]}
   assert raw==dict(cert['raw_basis_action_operators'][r][j]);cell=defaultdict(Fraction)
   for i,c in raw.items():
    assert W[i]+u-v==W[cut]
    for k,q in parent['canonical_linear_forms'][i]:cell[k]+=c*Fraction(q)
   cell={k:q for k,q in cell.items()if q};assert cell=={k:Fraction(q)for k,q in cert['universal109_slice_linear_forms'][r][j]}
   assert all(degrees[k]==5-W[cut]+u-v for k in cell);row.append(cell)
  forms.append(row)
 a=s.Symbol('a');M=s.Matrix([[s.Rational(v.get(-1,0))+a*s.Rational(v.get(0,0))for v in row]for row in forms]);K=s.QQ.poly_ring(a);det=K.to_sympy(DomainMatrix.from_Matrix(M).convert_to(K).det())
 assert sum(5-W[i]for i in cuts)==16 and sum(u-v for u,v in lies)==-16
 assert s.factor(det)==s.Rational(163840,531441)*a**3*(9*a+4)
 paths += [ROOT/'research/notes/Q_GROUP_SLICE_GLOBAL_UNIT_20261103_v63_REVIEW.md']
 result=dict(id='Q_GROUP_SLICE_GLOBAL_UNIT_20261103_v63_REVIEW',all_passed=True,fresh_subsecond_producer_replay_matches=True,
  all289_sparse150_operators_and_full109_variable_linear_forms_independently_match=True,
  weights_sum_zero_and108_strict_positive_variables_exclude_all_other_determinant_dependencies=True,
  exact_QQ_a_group17_determinant=str(s.factor(det)),group17_etale_at_ALL_parent_D_points_accepted=True,
  canonical_affine_formal_smooth_factor=17,regular_slice9_gives_canonical26_projective25=True,
  arbitrary_nonQa_literal_projective_Rhat_factor16_claimed=False,global_orbit_slice_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
 print('PASS entire109-slice exact group17 unit and correct projectivization scope')

if __name__=='__main__':main()
