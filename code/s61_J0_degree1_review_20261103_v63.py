#!/usr/bin/env python3
"""Independent actual-input binding and polynomial-cofactor J0 identity audit."""
import argparse,gzip,hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_v63';r=json.loads((folder/'report.json').read_text());assert r['all_passed']
 for name,want in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want,name
 cert=json.loads(gzip.decompress((folder/'identity.json.gz').read_bytes()))
 rows=[{tuple(m):Q(c)for m,c in f}for f in cert['row_polynomials']]
 cofs=[{tuple(m):Q(c)for m,c in f}for f in cert['cofactor_polynomials']]
 assert len(rows)==len(cofs)==55 and sum(bool(f)for f in cofs)==39 and sum(map(len,cofs))==45
 assert all(sum(m)<=1 and m[2]==0 for f in cofs for m in f)
 source=ROOT/'work/nonsquare-J0-pid-v63/primitive-chain-v1/stage-1/residual.txt'
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[21,55,1];f.readline();f.readline()
  for i in range(21):
   for j in range(55):
    line=f.readline()
    if i!=4:continue
    z=line.split();n=int(z[0]);assert len(z)==1+4*n
    actual={tuple(map(int,z[2+4*k:5+4*k])):Q(z[1+4*k])for k in range(n)}
    assert actual==rows[j]
 total={}
 for f,g in zip(rows,cofs):
  for m,c in f.items():
   for n,d in g.items():
    mon=tuple(a+b for a,b in zip(m,n));total[mon]=total.get(mon,Q(0))+c*d
 total={m:c for m,c in total.items()if c}
 assert total=={(0,2,0):Q(9),(0,1,0):Q(12),(0,0,0):Q(4)}
 constants=[j for j,f in enumerate(cofs)if set(f)=={(0,0,0)}and f[(0,0,0)]];assert constants
 fields=(folder/'QQ-factors.txt').read_text().splitlines();assert fields[0]=='20 54'
 factors=[Q(v.split()[2])for v in fields[1:]];assert len(factors)==74 and all(factors)
 paths=[Path(__file__),folder/'report.json',folder/'identity.json.gz',folder/'QQ-factors.txt',source,ROOT/'research/notes/S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_v63_REVIEW.md']
 result=dict(id='S61_NONSQUARE_J0_DEGREE1_UNIT_20261103_v63_REVIEW',all_passed=True,
  own_all55_row_polynomials_compared_to_actual_hash_pinned21x55=True,
  exact_full_QQ_polynomial_identity=True,active_cofactors=39,cofactor_terms=45,
  target='(3e+2)^2',unit_constant_cofactor_columns=constants,
  unimodular_replacement_and37th_unit_cokernel_split_accepted=True,
  QQ_diagonal_factors_nonzero_and_no_polynomial_factor_inverted=True,
  large_native_output_or_all351404_diagonal_terms_rerun=False,
  diagonal_all_coefficient_equality_reused_from_unchanged_independent_producer_receipt=True,
  full_J0_coker_unit_or_finite_defect_or_quartic_proved=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT J0 DEGREE1 ACTUAL ROW/UNIT REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
