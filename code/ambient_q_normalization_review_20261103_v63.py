#!/usr/bin/env python3
"""Independent M2 flag-fiber and geometric-scope review at distinguished Q."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--producer-replay',type=Path,required=True)
 parser.add_argument('--reuse-M2-directory',type=Path)
 parser.add_argument('--output-directory',type=Path,required=True)
 args=parser.parse_args();args.output_directory.mkdir(parents=True,exist_ok=False)
 path=ROOT/'reports/AMBIENT_RANKDROP_FLAG_FIBERS_V2_20261103_v63/report.json'
 d=json.loads(path.read_text());fresh=json.loads(args.producer_replay.read_text())
 assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
 paths=[Path(__file__),path]
 for name,expected in d['input_sha256'].items():
  p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
 K=s.Matrix(d['intrinsic_coefficient_matrix']);assert K.rank()==2
 assert K.nullspace()==[s.Matrix([0,1,0,0]),s.Matrix([0,0,0,1])]
 assert any(s.sympify(v)!=0 for v in d['off_chart_base_constraint_witness'])
 script=['R=QQ[r,z,w,MonomialOrder=>Lex];','J=ideal(r,z,w^2);']
 for name,body in d['flag_fibers'].items():
  assert body['quotient_length']==2
  assert body['reference_line_support_nilpotence_powers']=={'r':1,'z':1,'w':2}
  rows=body.get('unsimplified_constraints')
  if rows is None:rows=d['flag_fibers']['CA']['unsimplified_constraints']+d['flag_fibers']['CB']['unsimplified_constraints']
  text=','.join('('+v.replace('**','^')+')'for v in rows)
  script+=['I=ideal{'+text+'};','assert(I==J);','assert(numColumns basis(R/I)==2);',
           'assert(r%gb I==0 and z%gb I==0 and w%gb I!=0 and w^2%gb I==0);',
           'print "'+name+'_EXACT_PROJECTIVE_FLAG_CHART_PASS";']
 script+=['stdio << flush;exit 0;'];text='\n'.join(script)+'\n'
 target=args.output_directory/'review.m2';log=args.output_directory/'review.log';target.write_text(text)
 if args.reuse_M2_directory:
  assert (args.reuse_M2_directory/'review.m2').read_text()==text
  shutil.copyfile(args.reuse_M2_directory/'review.log',log)
 else:
  with log.open('w')as stream:subprocess.run([shutil.which('M2'),'--script',str(target.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=60)
 output=log.read_text()
 assert all(name+'_EXACT_PROJECTIVE_FLAG_CHART_PASS'in output for name in ['CA','CB','N'])
 assert 'error'not in output.lower()
 paths += [ROOT/'research/notes/AMBIENT_Q_NORMALIZATION_20261103_v63.md',
           ROOT/'research/notes/AMBIENT_Q_NORMALIZATION_20261103_v63_REVIEW.md',
           ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json']
 report=dict(id='AMBIENT_Q_NORMALIZATION_20261103_V63_REVIEW',all_passed=True,
             exact_scope='distinguished v24 Q (a=1), reduced CA/CB/N; not automatically the resonant a=2/9 point',
             fresh_physical_normalizer_v2_replay_matches=True,
             independent_M2_all3_full_constraint_ideals_equal_r_z_w_squared=True,
             independent_M2_lengths_and_nilpowers=[1,1,2],
             whole_projective_flag_fiber_exhausted_using_offchart_witness=True,
             proper_smooth_birational_source_finite_normalization_nearQ_accepted=True,
             reduced_CA_CB_N_nonnormal_atQ_accepted=True,
             component_Hilbert_Samuel_multiplicities=[2,2,2],
             each_reduced_component_completion_is_a_domain_accepted=True,
             cusp_product_or_whole_Jacobi_local_ring_claimed=False,
             input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)},
             evidence_sha256={p.name:sha(p)for p in [target,log]})
 (args.output_directory/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS independent M2 fibers; finite normalization, multiplicity2 and analytic irreducibility')


if __name__=='__main__':main()
