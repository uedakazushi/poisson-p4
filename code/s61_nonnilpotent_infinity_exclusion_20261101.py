#!/usr/bin/env python3
"""Finite exact rank stratification excluding all nonnilpotent infinity generically.

The archived v03/v04/v05 normal-form exhaustion is an explicit dependency.
Rank-drop values in the two-parameter families are fully enumerated.
"""
from __future__ import annotations
import argparse,hashlib,json
from itertools import combinations
from pathlib import Path
import sympy as s
from s61_model import S61,X
from s61_nil22_infinity_20261101 import tensor,beta_matrix,polynomial_minor,sm
from weight1_structure import J,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def log_tensor(A,B,C,D,E,P=0,R=0):
    x,y,z,w=X[1:]
    return tensor([A*x*y+B*z*w,E*x*z,(C-E)*x*w+P*y*y,(-B-E)*y*z+R*x*x,(B-C+E)*y*w,C*x*y+D*z*w])
def connected_blocks(M):
    columns=list(range(M.cols));blocks=[]
    while columns:
        cs={columns.pop(0)};rs=set()
        while True:
            old=(len(cs),len(rs));rs|={i for i in range(M.rows) if any(M[i,j]!=0 for j in cs)};cs|={j for j in columns if any(M[i,j]!=0 for i in rs)}
            if old==(len(cs),len(rs)):break
        columns=[j for j in columns if j not in cs];blocks.append((sorted(rs),sorted(cs),M.extract(sorted(rs),sorted(cs))))
    return blocks

def rank_one_zeros(B,variables):
    d,e=variables
    assert B!=s.zeros(*B.shape)
    assert all(s.expand(B[i,j]*B[k,l]-B[i,l]*B[k,j])==0 for i,k in combinations(range(B.rows),2) for j,l in combinations(range(B.cols),2))
    vals=[s.factor(v) for v in B if v!=0];gcd=s.gcd_list(vals);assert gcd.is_Rational and gcd!=0
    factors=set()
    for v in vals:
        for f,power in s.factor_list(v,d,e)[1]:
            assert s.Poly(f,d,e).total_degree()==1;factors.add(f)
    points=set()
    for f,h in combinations(sorted(factors,key=str),2):
        A,b=s.linear_eq_to_matrix([f,h],[d,e])
        if A.det()==0:continue
        point=tuple(A.inv()*b)
        if all(s.expand(v.subs(dict(zip(variables,point,strict=True))))==0 for v in vals):points.add(point)
    return points,{'rank_at_generic_point':1,'all_2_by_2_minors_zero':True,'entry_gcd':str(gcd),'linear_factors':list(map(str,sorted(factors,key=str))),'zero_points':[[str(v) for v in p] for p in sorted(points)]}

def equal_tensor(tt,ee,hh,aa=0,bb=0,cc=0,dd=0,ff=0,gg=0):
    x,y,z,w=X[1:]
    return tensor([(1-ee)*x*y-ff*x*w-gg*y*z+(1-hh)*z*w,tt*x*z,(-2*cc-gg)*x*y+(ee-tt)*x*w+aa*y*z+gg*z*w,ff*x*y+bb*x*w-(1+tt-hh)*y*z+(-2*dd-ff)*z*w,cc*y*y+(1+tt-ee-hh)*y*w+dd*w*w,ee*x*y+ff*x*w+gg*y*z+hh*z*w])

def run():
    m=S61();x,y,z,w=X[1:];d,e=s.symbols('d e');families=[]
    for label,Q,count in [('G1',log_tensor(1,0,0,d,e),27),('G2',log_tensor(1,d/2,0,d/2,e),23)]:
        M=beta_matrix(m,Q);blocks=connected_blocks(M);assert len(blocks)==12;points=set();saved=[]
        for rows,cols,B in blocks:
            zeros,info=rank_one_zeros(B,(d,e));points|=zeros;info.update({'rows':rows,'columns':cols,'matrix':sm(B)});saved.append(info)
        exceptions=[]
        for point in sorted(points):
            sub=dict(zip((d,e),point,strict=True));rank=sum(B.subs(sub).rank() for _,_,B in blocks);assert rank>=10
            exceptions.append({'d':str(point[0]),'e':str(point[1]),'rank':rank})
        assert len(exceptions)==count
        # Diagonal Cartan fixes every member, including all exceptional points.
        torus=s.diag(s.Symbol('a'),1/s.Symbol('a'),s.Symbol('b'),1/s.Symbol('b'));v=s.Matrix([x,y,z,w]);pulled=torus.inv()*Q.subs(dict(zip(v,torus*v,strict=True)),simultaneous=True)*torus.inv().T
        assert (pulled-Q).applyfunc(s.cancel)==s.zeros(4)
        families.append({'label':label,'Q':sm(Q),'blocks':saved,'rank_drop_points':exceptions,'rank_away_from_list':12,'uniform_rank_bound':10,'torus_dimension':2,'generic_beta_incidence_bound':14,'exceptional_beta_incidence_bound':14})
    # The A=0, B!=0 boundary of G2 is not covered by the A=1 chart.
    boundary_Q=log_tensor(0,s.Rational(1,2),0,s.Rational(1,2),e)
    boundary_M=beta_matrix(m,boundary_Q);boundary_blocks=[];boundary_roots=set();boundary_generic_rank=0
    for rows,cols,B in connected_blocks(boundary_M):
        vals=[s.factor(v) for v in B if v!=0]
        if not vals:
            boundary_blocks.append({'rows':rows,'columns':cols,'matrix':sm(B),'rank':0});continue
        assert all(s.expand(B[i,j]*B[h,l]-B[i,l]*B[h,j])==0 for i,h in combinations(range(B.rows),2) for j,l in combinations(range(B.cols),2))
        gcd=s.factor(s.gcd_list(vals));roots=s.solve(gcd,e) if gcd.has(e) else []
        assert all(root.is_Rational for root in roots)
        assert all(s.degree(f,e)==1 for f,power in s.factor_list(gcd,e)[1])
        boundary_generic_rank+=1;boundary_roots.update(roots)
        boundary_blocks.append({'rows':rows,'columns':cols,'matrix':sm(B),'rank':1,'entry_gcd':str(gcd),'rank_drop_roots':list(map(str,roots))})
    assert boundary_generic_rank==12 and boundary_roots=={-s.S.One,-s.Rational(1,2),s.S.Zero,s.Rational(1,2)}
    boundary_exceptional=[]
    for root in sorted(boundary_roots):
        rank=boundary_M.subs(e,root).rank();assert rank==11
        boundary_exceptional.append({'e':str(root),'rank':rank})
    aa,bb=s.symbols('a b');torus=s.diag(aa,1/aa,bb,1/bb);v=s.Matrix([x,y,z,w])
    assert (torus.inv()*boundary_Q.subs(dict(zip(v,torus*v,strict=True)),simultaneous=True)*torus.inv().T-boundary_Q).applyfunc(s.cancel)==s.zeros(4)
    boundary={'label':'G2_A_zero','Q':sm(boundary_Q),'blocks':boundary_blocks,'generic_rank':12,'rank_drop_points':boundary_exceptional,'uniform_rank_bound':11,'torus_dimension':2,'beta_incidence_bound':13}
    constant={}
    for name,params in [('T1',[1,0,0,3,1,1,0]),('T2',[s.Rational(1,2),0,s.Rational(1,2),3,1,1,0]),('T3',[1,s.Rational(3,2),0,s.Rational(3,2),-s.Rational(1,2),1,0]),('T4',[-1,0,2,3,1,1,3])]:constant[name]=log_tensor(*params)
    for dd in [s.S.Zero,s.Rational(1,2)]:
        aa,bb,cc=1,0,-dd/2;constant['one_zero_nil_'+str(dd)]=tensor([(1-dd)*x*y,(cc+dd)*x*z-aa*x*w,-bb*x*z-cc*x*w,(-cc-dd)*y*z+aa*y*w,bb*y*z+cc*y*w,dd*x*y])
    eqdata=[('C1',1,0,1,0,0,1,0,0,0),('C2',1,s.Rational(1,2),1,0,0,1,0,0,0),('D1',-1,0,1,0,0,0,1,0,0),('D2',-1,0,s.Rational(1,2),0,0,0,1,0,0),('A1',s.Rational(1,2),0,1,1,0,0,0,0,0),('A2',s.Rational(1,2),s.Rational(1,2),1,1,0,0,0,0,0),('B1',-s.Rational(1,2),0,1,0,1,0,0,0,0),('B2',-s.Rational(1,2),0,s.Rational(1,2),0,1,0,0,0,0),('Kplus',s.Rational(1,2),0,s.Rational(1,2),-s.Rational(1,2),0,1,0,0,-s.Rational(1,2)),('Kminus',-s.Rational(1,2),s.Rational(1,2),1,0,s.Rational(1,2),0,1,-s.Rational(1,2),0)]
    for name,*args in eqdata:constant[name]=equal_tensor(*args)
    fixed=[]
    for label,Q in constant.items():
        M=beta_matrix(m,Q);witness=minor_witness(M);assert witness['rank']>=11
        fixed.append({'label':label,'Q':sm(Q),'beta_matrix':sm(M),'minor':witness,'beta_image_dimension_upper_bound':26-witness['rank']})
    a=s.Symbol('a');mixed=[]
    for label,bb,cc in [('Z1',0,s.Rational(1,4)),('Z2',0,s.Rational(1,2)),('Z3',s.Rational(1,2),s.Rational(1,2))]:
        Q=tensor([(1-bb)*x*y+(s.Rational(1,2)-cc)*z*z,bb*x*z/2,(1-a-2*cc)*x*z+bb*x*w/2,-bb*y*z/2,a*y*z-bb*y*w/2,bb*x*y+cc*z*z]);M=beta_matrix(m,Q);rr,co,det=polynomial_minor(M,{a:0},[a]);assert len(co)>=11 and det!=0 and not det.has(a)
        A=s.Matrix([sum(s.diff(Q[i,j],X[j+1]) for j in range(4)) for i in range(4)]).jacobian(list(X[1:]));assert A!=s.zeros(4) and A.T*J+J*A==s.zeros(4)
        v=s.Matrix([x,y,z,w]);lie=Q.applyfunc(lambda f:sum((A*v)[i]*s.diff(f,v[i]) for i in range(4)))-A*Q-Q*A.T;assert lie.applyfunc(s.expand)==s.zeros(4)
        mixed.append({'label':label,'Q':sm(Q),'beta_matrix':sm(M),'minor_rows':rr,'minor_columns':co,'constant_minor':str(det),'uniform_rank_bound':len(co),'nonzero_modular_stabilizer':sm(A),'beta_image_dimension_upper_bound':26-len(co)})
    dependencies=[ROOT/'legacy/v03/proofs/regular_semisimple_theorem.txt',ROOT/'legacy/v04/semisimple_boundary.tex',ROOT/'legacy/v05/proofs/mixed_modular_theorem.txt']
    inputs=[Path(__file__),ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_infinity_flags.py',ROOT/'code/s61_model.py',m.source/'symplectic_slice_reduced_equations.json']+dependencies
    return {'id':'S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101','all_passed':True,'field':'QQ','log_families':families,'G2_A_zero_boundary':boundary,'fixed_representatives':fixed,'mixed_lines':mixed,'universal_beta_image_dimension_upper_bound':15,'historical_point_exhaustion_dependencies':[str(p.relative_to(ROOT)) for p in dependencies],'scope':'Every nonnilpotent modular point of the134-leading locus lies over a proper closed subset of W16. Thus excluded over sufficiently general beta. Does not assert exclusion over all Delta!=0, nor any affine fiber bound.','proof_note':'research/notes/S61_NONNILPOTENT_INFINITY_EXCLUSION_20261101.md','input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run();args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'log_rank_drop_counts':[27,23],'fixed_representatives':len(r['fixed_representatives']),'mixed_lines':3,'beta_image_dimension_bound':15,'output':str(args.output)}))
if __name__=='__main__':main()
