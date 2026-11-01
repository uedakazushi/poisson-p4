#!/usr/bin/env python3
"""Finite coverage audit of every denominator of the pure double-root orbit."""
import argparse,json,hashlib
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def verify():
 stems=['S61_NONBINARY_DOUBLE_ROOT_OPEN','S61_NONBINARY_DOUBLE_ROOT_ARTIFACTS','S61_NONBINARY_DOUBLE_ROOT_ENDPOINT','S61_NONBINARY_DOUBLE_ROOT_TORSION_DIVISORS','S61_NONBINARY_DOUBLE_ROOT_MIXED_JUNCTIONS','S61_NONBINARY_RESONANCE_EXCLUSIONS'];reports={};paths=[Path(__file__)]
 for stem in stems:
  path=ROOT/f'reports/{stem}_20261101.json';r=json.loads(path.read_text());assert r['all_passed'];reports[stem]=r;paths.append(path)
  for rel,sha in r.get('input_sha256',{}).items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==sha,(stem,rel)
 A,B,U,V=s.symbols('u14 u49 u26 u35');cover=[(U,3,'torsion_divisors'),(V,3,'torsion_divisors'),(B,5,'mixed_junctions'),(A-2*B,5,'minimal_endpoint'),(A-B,3,'torsion_divisors'),(A+B,3,'artifacts'),(2*A-B,0,'proper_beta_image'),(5*A-7*B,3,'torsion_divisors'),(5*A-4*B,0,'proper_beta_image'),(5*A-3*B,5,'mixed_junctions'),(5*A-B,0,'proper_beta_image'),(20*A*A-59*A*B+65*B*B,3,'artifacts'),(11*A**3-39*A*A*B+51*A*B*B-7*B**3,3,'artifacts')]
 r=reports[stems[0]];minor=s.sympify(r['jacobian_minor']);factors=s.factor_list(minor)[1]
 norm=lambda f:s.Poly(f,A,B,U,V,domain=s.QQ).monic().as_expr()
 assert {norm(f) for f,n in factors}=={norm(f) for f,n,k in cover}
 assert all(norm(s.sympify(f)) in {norm(g) for g,n,k in cover} for f in r['denominator_factors'])
 assert r['symmetric_zw2_orbit_verified'] and r['completed_ring_rank']==3 and r['homogenizer_zero']
 return dict(id='S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101',all_passed=True,coverage=[dict(factor=str(f),generic_local_length=n if n else None,method=k) for f,n,k in cover],number_of_distinct_QQ_factors=13,conclusion='Over sufficiently general beta every pure z^2*w or z*w^2 R0 orbit point is removed by t-saturation.',limitations=['Mixed tails, flag/binary loci and remaining minimal strata are separate.','No full affine-fiber finiteness, five-point exhaustion, or component upper bound.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,factors=r['number_of_distinct_QQ_factors'])))
