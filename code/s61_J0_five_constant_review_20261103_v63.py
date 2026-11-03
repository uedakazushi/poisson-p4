#!/usr/bin/env python3
"""Independent Fraction-dictionary verification of the second J0 unit-row certificate."""
import argparse,gzip,hashlib,json,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main(out):
 out.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 folder=ROOT/'reports/S61_NONSQUARE_J0_SECOND_UNIT_20261103_v63'
 report=json.loads((folder/'report.json').read_text());assert report['all_passed']
 for name,want in report['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==want,name
 cert=json.loads(gzip.decompress((folder/'identity.json.gz').read_bytes()))
 selected={int(j):{tuple(m):Q(c)for m,c in poly}for j,poly in cert['source_polynomials'].items()}
 source=ROOT/'work/nonsquare-J0-pid-v63/primitive-span-v1/residual.txt'
 with source.open()as f:
  assert list(map(int,f.readline().split()))==[22,56,1];f.readline();f.readline()
  for i in range(22):
   for j in range(56):
    line=f.readline()
    if i!=13 or j not in selected:continue
    z=line.split();n=int(z[0]);assert len(z)==1+4*n
    actual={}
    for k in range(n):
     assert int(z[4+4*k])==0
     actual[tuple(map(int,z[2+4*k:4+4*k]))]=Q(z[1+4*k])
    assert actual==selected[j]
 delta={}
 for j,multiplier in zip([7,9,26,36,46],[5757,2400,12320,-22945,8225]):
  for mon,c in selected[j].items():delta[mon]=delta.get(mon,Q(0))+multiplier*c
 delta={m:c for m,c in delta.items()if c}
 def mul(f,g):
  h={}
  for (a,b),c in f.items():
   for (d,e),v in g.items():h[a+d,b+e]=h.get((a+d,b+e),Q(0))+c*v
  return {m:c for m,c in h.items()if c}
 scalar=Q(cert['nonzero_QQ_scalar']);assert scalar
 target={(2,0):scalar}
 target=mul(target,{(1,0):Q(100),(0,2):Q(18),(0,1):Q(-9),(0,0):Q(-2)})
 for _ in range(2):target=mul(target,{(0,1):Q(3),(0,0):Q(2)})
 assert target==delta=={tuple(m):Q(c)for m,c in cert['full_delta_terms']}
 with (folder/'row13.witness').open()as f:
  assert int(f.readline())==56;witness={}
  for j in range(56):
   z=f.readline().split()
   if z[0]=='0':continue
   assert z[0]=='1'and z[2:]==['0','0','0'];witness[j]=Q(z[1])
 assert witness=={j:Q(c)/scalar for j,c in zip([7,9,26,36,46],[5757,2400,12320,-22945,8225])}
 paths=[Path(__file__),folder/'report.json',folder/'identity.json.gz',folder/'row13.witness',ROOT/'reports/S61_NONSQUARE_J0_PID_20261103_v63/report.json',source,ROOT/'research/notes/S61_NONSQUARE_J0_SECOND_UNIT_20261103_v63_REVIEW.md']
 result=dict(id='S61_NONSQUARE_J0_SECOND_UNIT_20261103_v63_REVIEW',all_passed=True,
  independent_arithmetic='fractions.Fraction coefficient dictionaries; no SymPy or native polynomial engine',
  all_five_source_polynomials_compared_to_hash_pinned_actual22x56=True,
  exact_QQ_unit_identity='sum_5_primitive_integer_coefficient_columns = c*a^2*C*(3e+2)^2',
  normalized_constant_replacement_coefficient=str(Q(5757)/scalar),
  rank_preserving21x55_cokernel_split_accepted=True,
  native_large_residual_regenerated=False,
  all_point_full_rank_or_J0_defect_finiteness_or_quartic_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('INDEPENDENT J0 FIVE-CONSTANT UNIT REVIEW PASS',result['seconds'])

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
