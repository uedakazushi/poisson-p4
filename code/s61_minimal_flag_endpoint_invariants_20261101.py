#!/usr/bin/env python3
"""Intrinsic equations for the four remaining minimal flag endpoint strata."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis
ROOT=Path(__file__).resolve().parents[1]
enc=lambda M:[[str(v)for v in row]for row in M.tolist()]
def run():
    st=time.time();old=json.loads((ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text());table=json.loads((ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json').read_text());b=s.symbols('b0:15');t=s.Symbol('t');H=s.Matrix(old['flag']['beta_kernel']);D=s.Matrix(old['frame_change']['beta_derivative']);stab=s.Matrix(old['frame_change']['stabilizer']);m=S61();z,u,B=primitive_beta(m);v=s.Matrix([B[i,j].diff(z[k])for(i,j),k in COORDS]);M=v.jacobian(u);rows=list(M.T.to_DM().rref()[1]);iv=M[rows,:].inv();BB=B.subs(dict(zip(u,H*s.Matrix(b))));finite={}
    for j in[3,4]:
        A=sum((stab[i,j]*sp_basis()[i]for i in range(10)),s.zeros(4));g=s.eye(4)+t*A;assert A*A==s.zeros(4)and g*J[1:,1:]*g.T==J[1:,1:];C=g.T*BB.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g;nv=(iv*s.Matrix([C[i,j].diff(z[k])for(i,j),k in COORDS])[rows,:]).applyfunc(s.expand);finite[j]=nv[1:,:]
    double={b[0]:0,b[1]:0};Ddouble=D.subs(double);ratios=[-s.S.One,-s.Integer(4),s.Rational(7,5)];ratio_cert=[]
    for r in ratios:
        f=b[3]-r*b[12];der=s.Matrix([[s.diff(f,x)for x in b]])*Ddouble;assert der==s.Matrix([[0,-f,0,0,0,0]]);moved=finite[4].subs(double);assert s.expand(moved[2]-b[2]-t*(b[12]-2*b[3]))==0;assert all(moved[i]==b[i].subs(double)for i in[0,1,3,12]);parameter=-b[2]/(b[12]-2*b[3]);assert s.cancel(moved[2].subs(t,parameter))==0;ratio_cert.append(dict(ratio=str(r),derivatives=enc(der),normalization_parameter=str(parameter),normalization_open='b12 !=0',normalization_coefficient=str(1-2*r)))
    delta=5*b[3]-4*b[12];sub={b[0]:0,b[3]:s.Rational(4,5)*b[12]};Ds=D.subs(sub);deltader=s.Matrix([[s.diff(delta,x)for x in b]])*D.subs(b[0],0);assert deltader==s.Matrix([[0,-delta,0,0,0,0]])
    P=b[1]*t*t-s.Rational(6,5)*b[12]*t+b[14]-3*b[5];Q=b[6]-(b[14]+b[5])*t+s.Rational(9,5)*b[12]*t*t-b[1]*t**3;flow=finite[3].subs(sub);assert s.expand(flow[14]-3*flow[5]-P)==0 and s.expand(flow[6]-Q)==0
    F=s.factor(125*s.resultant(P,Q,t)/b[1]);weights=[]
    for j in range(6):
        derivative=s.expand(sum(s.diff(F,x)*Ds[i,j]for i,x in enumerate(b)));q=s.cancel(derivative/F);assert not q.free_symbols;weights.append(q)
    assert weights==[0,-2,0,0,0,-2];disc=s.factor(s.discriminant(F,b[6]));assert disc==16*(50*b[1]*b[5]-9*b[12]**2)**2*(-25*b[1]*b[14]+75*b[1]*b[5]+9*b[12]**2);poly=s.Poly(F,b[6]);content=s.gcd_list(poly.all_coeffs());assert content==1;assert s.degree(-25*b[1]*b[14]+75*b[1]*b[5]+9*b[12]**2,b[14])==1
    puresub={b[4]:0,b[6]:0,b[14]:3*b[5]};assert s.expand(F.subs(puresub))==0
    # Both target equations contain the exact source H10, whose orbit closure has dimension16.
    pure=s.Matrix(table['strata']['pure_cubic']['beta']);assert pure[0]==0 and pure[1]==0 and s.expand((H*s.Matrix(b)-pure)[0])==0;bs=dict(zip(b,pure[1:,:]));assert s.expand(delta.subs(bs))==0 and s.expand(F.subs(bs))==0
    deps=['reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json','code/s61_model.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py'];paths=[Path(__file__)]+[ROOT/f for f in deps]
    return dict(id='S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101',all_passed=True,finite_beta_flows={str(j):list(map(str,v))for j,v in finite.items()},double_boundary_ratio_cases=ratio_cert,pure_cubic=dict(delta=str(delta),P=str(P),Q=str(Q),quartic=str(F),weights=list(map(str,weights)),discriminant_in_b6=str(disc),primitive_content='1',absolute_irreducibility_reason='Primitive quadratic in b6 with nonsquare linear factor in b14 in its discriminant.',source_containment=True),conclusion='The three double-boundary strata are b0=b1=0,b3=r*b12; the pure-cubic stratum is b0=delta=F=0.',scope=['Intrinsic exact equations for the four flag endpoint closures; no saturation-removal assertion.','Uses connectedness of the projective flag stabilizer and the already-certified16D source images.'],proof_note='research/notes/S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.time()-st)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
    with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
