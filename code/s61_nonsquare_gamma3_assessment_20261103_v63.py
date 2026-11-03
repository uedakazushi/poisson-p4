#!/usr/bin/env python3
"""Read-only binding/scope receipt for the proposed all-U free Gamma3 quotient."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import tarfile
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_NONSQUARE_GAMMA3_FREE_QUOTIENT_ASSESSMENT_20261103_v63'
NOTES=[
 'S61_BINARY_TRANSVERSE_ALGEBRA_20261101',
 'S61_BINARY_FORMAL_RING_20261101','S61_BINARY_FORMAL_RING_REVIEW_20261101',
 'S61_BINARY_UNIFORM_CUBIC_20261101',
 'S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102',
 'S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102',
 'S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102',
 'S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_REVIEW_20261102',
 'S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101',
 'S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62',
 'S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62',
 'S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102',
 'S61_NONSQUARE_J0_EXISTING_DEFECT_AUDIT_20261103_v63',
]

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def assess(output):
 start=time.monotonic();base=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62'
 source=json.loads((base/'source-report.json').read_text());assert source['all_passed']
 cols=source['normal_monomials'];weights=(1,1,1,1,1,2,1,1)
 assert len(cols)==91 and len(set(map(tuple,cols)))==91
 expected={ex for ex in itertools.product(range(4),repeat=8) if sum(x*w for x,w in zip(ex,weights))==3}
 assert set(map(tuple,cols))==expected
 assert sum(ex[5]==0 for ex in cols)==84 and sum(ex[5]==1 for ex in cols)==7
 with tarfile.open(base/'matrices.tar.gz') as archive:
  f=archive.extractfile('full');h=hashlib.sha256()
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
  fullsha=h.hexdigest()
 assert fullsha==source['full_matrix_sha256']=='dae59904925260605628d925ae0565de084c9af0ae1bccf1774f1af067830ad6'
 rowpath=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62/report.json'
 row=json.loads(rowpath.read_text());assert row['all_passed'] and row['global_free_row_module_rank']==57
 assert not row['free_cokernel_claimed'] and not row['all_point_cubic_rank57_claimed']
 genericpath=ROOT/'reports/S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.json'
 generic=json.loads(genericpath.read_text());assert generic['all_passed'] and generic['generic_weighted_ranks']==[6,57,239]
 formalpath=ROOT/'reports/S61_BINARY_FORMAL_RING_20261101.json'
 formal=json.loads(formalpath.read_text());assert formal['all_passed'] and formal['associated_graded_free_ranks']==[1,7,23,34]
 assert formal['completed_ring_finite_free_over']=='C[[v0,v1,v3]]'
 horizontalpath=ROOT/'reports/S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102.json'
 horizontal=json.loads(horizontalpath.read_text());assert horizontal['all_passed'] and horizontal['weighted_algebra_free_rank_over_Kd']==65
 assert 'K=QQ(a,e)' in horizontal['scope']
 x,y,z=s.symbols('x y z');m=s.Matrix([-y,x]);p=s.Matrix([[x,y]])
 assert p*m==s.zeros(1,1) and m.subs({x:1,y:0}).rank()==1 and m.subs({x:0,y:0}).rank()==0
 assert p.subs({x:0,y:0})==s.zeros(1,2)
 files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',base/'source-report.json',base/'report.json',base/'matrices.tar.gz',rowpath,genericpath,formalpath,horizontalpath]
 files += [ROOT/'research/notes'/f'{name}.md' for name in NOTES]
 result=dict(id=STEM,all_passed=True,kind='read-only structural assessment; no new all-U rank theorem',
   actual_full_matrix_sha256=fullsha,actual_shape=[63,91],weighted_columns=91,light_cubics=84,heavy_light=7,
   global_row_module_free_rank=57,existing_full_cubic_defect_codimension_at_least=2,
   all_U_free_Gamma3_quotient_premise_found=False,nonexistence_of_such_quotient_proved=False,
   conditional_determinant_argument_valid=True,locally_free_quotient_suffices=True,
   missing_finite_certificate='P(34x91)*M=0 and P*T=I34 over only known-unit denominators, or independent all-U flat thickening',
   counterexample=dict(relation=['-y','x'],quotient_image_ideal='(x,y)',defect='V(x,y) in A3',
      padded_shape=[57,91],padded_generic_cokernel_rank=34,padded_special_cokernel_rank=35,
      torsion_free_cokernel=True,double_dual_free=True,map_to_double_dual_surjective=False),
   no_GB_or_existing_certificate_rerun=True,
   input_sha256={str(p.relative_to(ROOT)):sha(p) for p in dict.fromkeys(files)},seconds=time.monotonic()-start)
 with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print('PASS',STEM,result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assess(a.output)
