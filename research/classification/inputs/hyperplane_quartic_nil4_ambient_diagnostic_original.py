from pathlib import Path
import json,hashlib,time
source_path=Path('code/classification_two_hyperplanes_cubic_component_review.py')
source=source_path.read_text();prefix=source[:source.index('report_path=')]
exec(prefix)
P=tensor({(1,2):add(term(1,1,2),term(-2,2,2)),(1,3):add(term(2,1,3),term(1,2,4)),(1,4):add(term(-1,1,4),term(4,2,4)),(2,3):term(2,2,3),(2,4):term(1,2,4),(3,4):term(2,3,4)})
assert all(not p for p in strict_jacobi(P))
D=[[0]*150 for _ in row_basis]
for col,(ij,m) in enumerate(column_basis):
 variation=tensor({ij:{m+(0,)*(N-5):F(1)}})
 mixed=tensor({kl:add(P.get(kl,{}),mul(term(1,19),variation.get(kl,{})))for kl in pairs})
 for inds,p in projective_jacobi(mixed).items():
  for exponent,a in p.items():
   if exponent[19]==1:
    assert a.denominator==1 and not any(exponent[5:19])
    D[row_index[inds,exponent[:5]]][col]=int(a)
cert=modular_minor(D)
print(json.dumps(dict(D=cert['rank_mod_prime'],tangent_upper=125-cert['rank_mod_prime'],seconds=time.monotonic()-started)),flush=True)
r=dict(task='E_reduced_nonnormal',point='twice type20 topQ: xy-2y²,2xz+yw,-xw+4yw,2yz,yw,2zw',
 arithmetic_source=str(source_path),arithmetic_source_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
 Dshape=[350,150],Dmatrix_sha256=matrix_hash(D),rank_lower_certificate=cert,projective_tangent_upper_bound=125-cert['rank_mod_prime'],
 scope='Single exact fullambient tangent diagnostic; no coverage assertion without source uniqueness and knowncomponent lower bound.',seconds=time.monotonic()-started)
Path('reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TYPE20_TOP_AMBIENT_DIAGNOSTIC_20261107.json').write_text(json.dumps(r,indent=2)+'\n')
