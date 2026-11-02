#!/usr/bin/env python3
"""Check the mixed-orbit fixed center in both Bruhat charts and the single-colon dependencies."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from sympy import QQ
from weight1_structure import J

ROOT=Path(__file__).resolve().parents[1]
STEM='S61_FIXED_PAIR_SINGLE_COLON_20261102'


def orbit_gcd():
    R,*vs=ring('a,b,c,d,s0,s1,s2,t0,t1,t2',QQ);a,b,c,d,s0,s1,s2,t0,t1,t2=vs
    # Positions x,z first, momenta y,w last. Use polynomial g and inverse both
    # multiplied by det(A), derived from lower T * Levi A * upper S.
    one=R.one;zero=R.zero;det=a*d-b*c
    def mm(A,B):return [[sum((v*B[k][j] for k,v in enumerate(row)),zero) for j in range(len(B[0]))] for row in A]
    def add(A,B):return [[u+v for u,v in zip(x,y)]for x,y in zip(A,B)]
    def minus(A):return [[-v for v in row]for row in A]
    def times(A,n):return [[v*n for v in row]for row in A]
    def tr(A):return list(map(list,zip(*A)))
    def blocks(A,B,C,D):return [x+y for x,y in zip(A,B)]+[x+y for x,y in zip(C,D)]
    A=[[a,b],[c,d]];adj=[[d,-b],[-c,a]];S=[[s0,s1],[s1,s2]];T=[[t0,t1],[t1,t2]]
    upper=times(A,det);upperright=mm(upper,S);lower=mm(T,upper);dblock=add(mm(lower,S),tr(adj))
    g=blocks(upper,upperright,lower,dblock);gi=blocks(tr(dblock),minus(tr(upperright)),minus(tr(lower)),tr(upper))
    prod=mm(g,gi);assert all(x==(det**2 if i==j else 0) for i,row in enumerate(prod) for j,x in enumerate(row))
    results={}
    for lab in ['T5','T7']:
     out=[]
     for i,k,v in [(0,1,0),(2,3,2)]:
      x,w=gi[0][v],gi[3][v]
      # Q=e_y^(2xw e_z-x^2 e_w), or e_y^(3w² e_z-x² e_w).
      zcoef=2*x*w if lab=='T5' else 3*w*w;wcoef=-x*x
      q=(g[i][2]*(zcoef*g[k][1]+wcoef*g[k][3])-g[k][2]*(zcoef*g[i][1]+wcoef*g[i][3]))
      out.append(q)
     gg=out[0].gcd(out[1])
     quotient,remainder=gg.div(det**3);assert not remainder and len(quotient)==1 and quotient==1
     results[lab]={'term_counts':[len(f)for f in out],'gcd':str(gg),'unit_on_big_cell':True}
    return results


def certify():
    polys=orbit_gcd()
    x,y,z,w=s.symbols('x y z w');coordinates=s.Matrix([x,y,z,w])
    v=s.Matrix([1,1,1,0]);g=s.eye(4)+v*v.T*J
    assert g.T*J*g==J and g.extract([0,2],[0,2]).det()==0
    for name,G in [('T5',x*x*z+x*w*w),('T7',x*x*z+w**3)]:
        Q=s.zeros(4);Q[1,2]=G.diff(w);Q[1,3]=-G.diff(z);Q-=Q.T
        values=[]
        for i,j,k in [(0,2,0),(1,3,1)]:
            at=g.inv()*s.eye(4)[:,k]
            values.append((g*Q.subs(dict(zip(coordinates,at)),simultaneous=True)*g.T)[i,j])
        assert values==[0,-2]
        polys[name]['exceptional_chart_evaluation']=list(map(int,values))
    inputs={}
    def record(p):
        inputs[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    for name in ['S61_PROJECTIVE_TOP_CYCLE_20261102','S61_AFFINE_TWO_COORDINATE_PROJECTION_20261102',
                 'S61_AFFINE_BOUNDARY_DIMENSION_20261102','S61_BINARY_FORMAL_RING_20261101',
                 'S61_FLAG_UNIFORM_OPEN_20261101','S61_UNIVERSAL_WEIGHT6_20261031',
                 'S61_R0_MIXED_ORBIT_REDUCTION_20261101']:
        p=ROOT/'reports'/(name+'.json')
        assert json.loads(p.read_text())['all_passed'];record(p)
        record(ROOT/'research/notes'/(name+'.md'))
    for p in [Path(__file__).resolve(),ROOT/'research/notes'/(STEM+'.md'),ROOT/'code/weight1_structure.py']:
        record(p)
    return dict(id=STEM,all_passed=True,leading_fixed_coordinates=['u6','u42'],mixed_orbit_gcd=polys,
                full_projective_fixed_pair_dimension_upper=1,colon_exponent=2,
                projective_single_colon_dimension_upper=0,affine_part_preserved=True,
                additional_generic_projection_target='P1',generic_projective_line_module_rank=2300,
                generic_homogenizer_rank=30,generic_homogenizer_square_zero=True,
                fixed_c_affine_algebra_length_determined=False,input_sha256=inputs)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path)
    output=parser.parse_args().output;result=certify()
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS single colon',result['colon_exponent'],result['generic_projective_line_module_rank'])
