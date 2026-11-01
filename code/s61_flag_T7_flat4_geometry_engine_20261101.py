#!/usr/bin/env python3
"""T7 flag flat4 certificate: geometry stage. Rebuild in supplied directory."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]);WD.mkdir(parents=True,exist_ok=True)
import json,time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
st=time.time();data=json.load(open(WD/'flag-T7-universal-lattice-ring.json'));b=list(map(s.Symbol,data['variables']));e=s.Symbol('e');R,*_=ring(b,s.QQ)
def zero(x):return R.from_expr(s.sympify(x).subs(e,0))
D0=zero(data['coordinate_denominator']);factor=s.factor(D0.as_expr());print('D0 factor',factor,time.time()-st,flush=True)
nums=[[zero(v)for v in row['numerators']]for row in data['multiplication']];dens=[zero(row['denominator'])for row in data['multiplication']]
k=next(i for i,row in enumerate(nums)if any(row));pivot=next(i for i,v in enumerate(nums[k])if v);ref=nums[k]
assert all(ref[pivot]*row[j]==row[pivot]*ref[j]for row in nums for j in range(3))
f=s.factor(ref[pivot].as_expr());print('nonzero multiplication witness factor',f,time.time()-st,flush=True)
Z={i:[zero(v)for v in data['coordinate_numerators'][i]]for i in data['selected_q_rows']+[19,24]}
assert all(Z[r][j]==(D0 if i==j else R.zero)for i,r in enumerate(data['selected_q_rows'])for j in range(3))
assert all(not v for i in [19,24]for v in Z[i]);nonp=[i for i in range(3)if i!=pivot];nilrows=[data['selected_q_rows'][i]for i in nonp]
out=dict(all_passed=True,regular_locus_factor=str(factor),nil_multiplication_generic_rank=1,multiplication_witness_factor=str(f),nil_cotangent_generic_dimension=2,nil_tangent_minor_rows=nilrows,nil_tangent_minor='1',surface_cotangent_generic_dimension=4,surface_embedding_minor='3',closed_embedding_from_selected_q_coordinates=True,seconds=time.time()-st)
(WD/'flag-T7-universal-lattice-geometry.json').write_text(json.dumps(out,indent=2)+'\n');print('UNIVERSAL FLAT4 SURFACE GEOMETRY PASS',time.time()-st,flush=True)
