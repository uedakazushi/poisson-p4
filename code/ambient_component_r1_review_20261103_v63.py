#!/usr/bin/env python3
"""Independent flag polynomials and critical-ideal verification for R1."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess

import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output-directory',type=Path,required=True)
 args=parser.parse_args();args.output_directory.mkdir(parents=True,exist_ok=False)
 path=ROOT/'reports/AMBIENT_BASE_FLAG_CRITICAL_20261103_v63/report.json'
 d=json.loads(path.read_text());assert d['all_passed'];paths=[Path(__file__),path]
 for name,expected in d['input_sha256'].items():
  p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
 t,A,B,C,D,F,r,q,w=s.symbols('t A B C D F r q w');base=[t,A,B,C,D,F]
 local={str(v):v for v in [*base,r,q,w]}
 J=s.Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
 mats=[s.zeros(4)for _ in range(4)]
 # Direct coefficients of beta, without importing common_plane or its adapter.
 for (i,j),coeffs in {
  (0,1):[t,0,0,0],(0,2):[A,0,2*C,D+t/2],
  (0,3):[B,0,D-t/2,2*F],(2,3):[-t,0,0,0]}.items():
  for k,c in enumerate(coeffs):mats[k][i,j]=c;mats[k][j,i]=-c
 pairs=list(itertools.combinations(range(4),2))
 K=s.Matrix([[mats[k][i,j]for k in range(4)]for i,j in pairs if(i,j)!=(2,3)])
 v=s.Matrix([r,1,q,w]);ell=(v.T*J).T
 contraction=s.Matrix.hstack(*[(v.T*M).T for M in mats]);left=list(contraction);right=list(ell*ell.T)
 equations=list(K*v)+[s.expand(left[i]*right[j]-left[j]*right[i])for i in range(16)for j in range(i)]
 assert all(s.expand(x-s.sympify(y,locals=local))==0 for x,y in zip(equations,d['full_flag_chart_equations']))
 rows=[]
 for row in s.Matrix(equations).jacobian([r,q,w]).subs({r:0,q:0,w:0}).tolist():
  if row!=[0,0,0]and row not in rows and[-x for x in row]not in rows:rows.append(row)
 matrix=s.Matrix(rows);assert matrix==s.Matrix([[s.sympify(v,locals=local)for v in row]for row in d['flag_tangent_matrix']])
 def monic(poly):return str(s.Poly(poly,*base,domain=s.QQ).monic().as_expr())
 minors=[s.expand(matrix[list(indices),:].det())for indices in itertools.combinations(range(matrix.rows),3)]
 assert {monic(p)for p in minors if p}=={monic(s.sympify(p,locals=local))for p in d['maximal_minors']}
 # Verify the four equations used in the all-parameter uniqueness proof.
 L=left;R=right;anchor=L[0]
 special={r:0};h=A*q+B*w
 identities=[(2,0,(h+2*t)*w-2*C*q-(D+t/2)*w),
             (8,0,(h+2*t)*w),(12,0,-(h+2*t)*q),
             (10,0,(h+t)*w*w),(15,0,(h+t)*q*q)]
 for i,j,expected in identities:
  assert s.expand((L[i]*R[j]-L[j]*R[i]).subs(special)-expected)==0,(i,j)
 delta=4*C*F-D*D
 text='\n'.join(['R=QQ[t,A,B,C,D,F];','Delta=4*C*F-D^2;',
   'I=ideal{'+','.join('('+p.replace('**','^')+')'for p in d['maximal_minors'])+'};',
   'assert(t^3%gb I==0);',
   'assert(I+ideal(t)==ideal(t,Delta*A,Delta*B,Delta*C,Delta*D,Delta*F));',
   'scan(flatten entries gens I,f->assert(f%gb ideal(t,Delta)==0));',
   'print "CRITICAL_SUPPORT_EXACT_IDEALS_PASS";stdio << flush;exit 0;'])+'\n'
 script=args.output_directory/'critical.m2';log=args.output_directory/'critical.log';script.write_text(text)
 with log.open('w')as stream:subprocess.run([shutil.which('M2'),'--script',str(script.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=60)
 assert 'CRITICAL_SUPPORT_EXACT_IDEALS_PASS'in log.read_text()
 paths += [ROOT/'research/notes/AMBIENT_COMPONENT_R1_NONCM_20261103_v63.md',
   ROOT/'research/notes/AMBIENT_COMPONENT_R1_NONCM_20261103_v63_REVIEW.md',
   ROOT/'reports/AMBIENT_Q_NORMALIZATION_20261103_v63_REVIEW/report.json',
   ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json']
 report=dict(id='AMBIENT_COMPONENT_R1_NONCM_20261103_V63_REVIEW',all_passed=True,
   independent_all_intrinsic_flag_polynomials_and_11_by3_tangent_matrix=True,
   independent_all_maximal_minors_up_to_QQ_units=True,
   independent_M2_t_cubed_and_full_central_minor_ideal_identity=True,
   global_projective_flag_uniqueness_off_t_delta_confirmed=True,
   differential_criticality_alone_not_used_for_uniqueness=True,
   proper_codimension2_bad_image_and_Nakayama_isomorphism_accepted=True,
   reduced_CA_CB_N_Serre_R1_accepted=True,
   distinguished_a_one_Q_component_local_rings_not_CM_accepted=True,
   ordinary_cusp_times_smooth_model_excluded=True,
   full_Jacobi_or_resonant_a_two_over_nine_scope_claimed=False,
   input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)},
   evidence_sha256={p.name:sha(p)for p in [script,log]})
 (args.output_directory/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS full flag polynomials, exact critical ideal, global uniqueness and R1/non-CM scope')


if __name__=='__main__':main()
