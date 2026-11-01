#!/usr/bin/env python3
"""Exact tangent and quadratic obstruction at a nonsquare binary infinity point."""
import sys,json,time
from pathlib import Path
import sympy as s
from itertools import combinations_with_replacement
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize

def verify():
    start=time.time();model=S61();from s61_model import X,J
    x,y,z,w=X[1:];B=X[0]**2*J;B[2,3]+=x*x+w*w;B[3,2]-=x*x+w*w
    pt=model.normalize(B)['parameters']
    for j,i in enumerate([i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]):pt[i]=s.Integer(j+2)
    pt[11]=0;pt[12]=0
    beta={i:pt[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    two,polys=specialize(model,beta,True);q0=s.Matrix([pt[i] for i in two]);k=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=k]
    A=s.zeros(135,45);B=s.zeros(135,1)
    for row,f in enumerate(polys):
     h=0
     for mon,c in f.items():
      if len(mon)==2:
       i,j=mon;h+=c*q0[i]*q0[j];A[row,i]+=c*q0[j];A[row,j]+=c*q0[i]
      elif len(mon)==1:B[row,0]+=c*q0[mon[0]]
     assert h==0
    Jac=A[:,cols].row_join(B);tan=s.Matrix.hstack(*Jac.nullspace());co=s.Matrix.hstack(*Jac.T.nullspace()).T
    assert A.rank()==34 and Jac.rank()==35 and tan.shape==(45,10) and tan[-1,:]==s.zeros(1,10)
    n=tan.cols;rows={i:list(tan[cols.index(i),:]) if i!=k else [0]*n for i in range(45)};tr=list(tan[-1,:]);mons=list(combinations_with_replacement(range(n),2));cf=[]
    for f in polys:
     out={}
     for mon,c in f.items():
      ls=[rows[i] for i in mon]+[tr]*(2-len(mon))
      for i in range(n):
       for j in range(n):
        ij=tuple(sorted((i,j)));out[ij]=out.get(ij,0)+c*ls[0][i]*ls[1][j]
     cf.append([out.get(mon,0) for mon in mons])
    proj=(co.to_DM()*s.Matrix(cf).to_DM()).to_Matrix();_,piv=proj.T.to_DM().rref();basis=proj[list(piv),:];v=s.symbols('v:'+str(n));eq=[s.factor(sum(c*v[i]*v[j] for c,(i,j) in zip(row,mons))) for row in basis.tolist()]
    assert len(piv)==5
    H=s.Matrix.vstack(*[s.hessian(f,v) for f in eq]);assert H.rank()==7
    result=dict(point=list(map(str,pt)),chart_fixed_q=k,chart_value=str(q0[k]),jacobian=[[str(x) for x in row] for row in Jac.tolist()],tangent=[[str(x) for x in row] for row in tan.tolist()],cokernel_rows=[[str(x) for x in row] for row in co[list(piv),:].tolist()],quadratic_monomials=mons,quadratic_obstruction_basis=[[str(x) for x in row] for row in basis.tolist()],equations=list(map(str,eq)),seconds=time.time()-start)

    from primitive_jet_elliptic_curve_20261101 import primitive_beta
    zz,uu,form=primitive_beta(model);form=form.subs({uu[i]:beta[j] for i,j in enumerate([j for j,w in enumerate(model.reduced['parameter_weights']) if w==1])})
    a,b,d=s.symbols('a b d');v1=s.Matrix([a,1,0,b]);v2=s.Matrix([-b,0,1,d])
    sec=s.Matrix([(v1.T*form.subs(dict(zip(zz,v1)))*v2)[0],(v1.T*form.subs(dict(zip(zz,v2)))*v2)[0]])
    assert sec.subs({a:0,b:0,d:0})==s.zeros(2,1)
    jsec=sec.jacobian([a,b,d]).subs({a:0,b:0,d:0});assert jsec.rank()==2
    result.update(id='S61_BINARY_QUADRATIC_OBSTRUCTION_20261101',all_passed=True,field='QQ',quadratic_span_rank=5,common_hessian_rank=7,canonical_curve_derivative=[[str(x) for x in row] for row in jsec.tolist()],canonical_curve_smooth=True)
    return result


def main():
    import argparse,hashlib
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=verify();root=Path(__file__).resolve().parents[1]
    inputs=[Path(__file__),root/'code/s61_model.py',root/'code/s61_flag_transverse_obstruction_20261101.py',root/'code/primitive_jet_elliptic_curve_20261101.py',root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    result['source_sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,tangent_dimension=10,quadratic_obstructions=5)))
if __name__=='__main__':main()
