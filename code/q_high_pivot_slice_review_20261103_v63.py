#!/usr/bin/env python3
"""Independent complete closed-slice ideal and rank2 Jacobian audit."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from q_smith_global_model_review_20261103_v63 import ROOT,read,sha
from q_local_q51_review_20261103_v63 import evaluate

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 folder=ROOT/'reports/Q_GLOBAL_HIGH_PIVOT_CLOSED_SLICE_20261103_v63';public=folder/'report.json';d=read(public);f=read(args.replay/'report.json');assert{k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
 packed=folder/'certificate.json.gz';assert packed.read_bytes()==(args.replay/'certificate.json.gz').read_bytes();c=read(packed);paths=[Path(__file__),public,packed,ROOT/'code/q_local_q51_review_20261103_v63.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py']
 for file,h in d['input_sha256'].items():p=ROOT/file;assert sha(p)==h,file;paths.append(p)
 for file,h in d['evidence_sha256'].items():assert sha(folder/file)==h
 parent=read(ROOT/'reports/Q_SMITH_GLOBAL_FINITE_MODEL_20261103_v63/model37.json.gz');R,a,u,v=ring('a,s,t',s.QQ);point=[a]+[R.zero]*36;point[35]=u;point[36]=v
 values=evaluate(parent['nodes'],point,R);raw=[values[i]for i in parent['raw150_reconstruction_nodes']];res=[values[i]for i in parent['all350_residual_nodes']]
 assert raw==[R.from_expr(s.sympify(p))for p in c['all150_raw_coefficients']]and res==[R.from_expr(s.sympify(p))for p in c['all350_residuals']]
 generators=[a*u,a*v];forward=[[R.from_expr(s.sympify(p))for p in row]for row in c['forward350_by2']];reverse=[[R.from_expr(s.sympify(p))for p in row]for row in c['reverse2_by350']]
 assert len(forward)==350 and len(reverse)==2
 assert all(f==row[0]*generators[0]+row[1]*generators[1]for f,row in zip(res,forward))
 assert all(g==sum((p*f for p,f in zip(row,res)),R.zero)for g,row in zip(generators,reverse))
 assert sum(bool(p)for p in res)==15 and all(mon[1]+mon[2]<=1 for p in raw for mon in p)
 X=s.symbols('x0:5');ss,tt=s.symbols('s t');H=X[0]*X[2]**2/3+X[1]**2*X[2]-s.Rational(4,3)*X[0]*X[1]*X[3]-s.Rational(2,9)*X[0]**2*X[4]-ss*X[3]*X[4]**2+tt*X[4]**3/3
 assert s.expand(H-s.sympify(c['a_zero_cubic_H']))==0
 P=s.zeros(5);pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
 for k,(i,j)in enumerate(pairs):
  P[i,j]=s.expand(sum(raw[15*k+n].as_expr().subs(s.Symbol('a'),0)*X[p]*X[q]for n,(p,q)in enumerate(mons)));P[j,i]=-P[i,j]
 assert P[2,3]==s.diff(H,X[4])and P[2,4]==-s.diff(H,X[3])and P[3,4]==s.diff(H,X[2])
 assert all(P[i,j]==0 for i,j in pairs if(i,j)not in[(2,3),(2,4),(3,4)])
 for i,j,k in itertools.combinations(range(5),3):assert s.expand(sum(P[i,l]*s.diff(P[j,k],X[l])+P[j,l]*s.diff(P[k,i],X[l])+P[k,l]*s.diff(P[i,j],X[l])for l in range(5)))==0
 assert all(s.expand(sum(s.diff(P[i,j],X[i])for i in range(5)))==0 for j in range(5))and P[3,4]!=0
 paths += [ROOT/'research/notes/Q_GLOBAL_HIGH_PIVOT_CLOSED_SLICE_20261103_v63_REVIEW.md']
 result=dict(id='Q_GLOBAL_HIGH_PIVOT_CLOSED_SLICE_20261103_v63_REVIEW',all_passed=True,
  fresh_producer_and_entire_certificate_match=True,
  independent_complete_parent_DAG_QQ_a_s_t_pullback_matches_all150_and350=True,
  all350_forward_and2_reverse_polynomial_ideal_certificates=True,
  exact_closed_slice_ideal='(a*s,a*t)=(a) intersection (s,t)',reduced_no_embedded_primes=True,
  independent_a_zero_rank2_actual_Jacobian_tensor_and_divergence_all_Jacobi=True,
  existing_Rat113_membership_by_exact_dx0_dx1_dH_definition=True,
  closed_slice_a_nonflat_only_full_parent_nonflatness_inference=False,
  full_Qzero_branches_or_new_rank2_classification_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
 with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS complete350 closed ideal(as,at), reduced decomposition and actual Rat113 rank2 plane')

if __name__=='__main__':main()
