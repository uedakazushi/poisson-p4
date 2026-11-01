#!/usr/bin/env python3
"""Five-variable normal form for the free-M0 flag endpoint stratum."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,J,X
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from weight1_structure import COORDS,sp_basis,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def run():
    start=time.time()
    m=S61();old=json.load(open(ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json'));b=s.symbols('b0:15');t=s.Symbol('t');H=s.Matrix(old['flag']['beta_kernel']);stab=s.Matrix(old['frame_change']['stabilizer']);D=s.Matrix(old['frame_change']['beta_derivative']).applyfunc(s.sympify);z,u,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(z[k])for(i,j),k in COORDS]);M=vec.jacobian(u);rows=list(M.T.to_DM().rref()[1]);iv=M[rows,:].inv();sub={b[0]:0,b[1]:0,b[2]:0,b[3]:-b[12]};base=H*s.Matrix(b);base=base.subs(sub);BB=B.subs(dict(zip(u,base)));flows={}
    xx=s.Matrix(X[1:]);Q=s.zeros(4);Q[1,3]=-xx[0]**2;Q-=Q.T
    for j in[0,2,3]:
     A=sum((stab[i,j]*sp_basis()[i]for i in range(10)),s.zeros(4));g=s.eye(4)+t*A;assert A*A==s.zeros(4);assert g*J[1:,1:]*g.T==J[1:,1:];assert(g.inv()*Q.subs(dict(zip(xx,g*xx)),simultaneous=True)*g.inv().T-Q).applyfunc(s.expand)==s.zeros(4);C=g.T*BB.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g;nv=(iv*s.Matrix([C[i,j].diff(z[k])for(i,j),k in COORDS])[rows,:]).applyfunc(s.expand);flows[j]=nv[1:,:]
    y=s.Matrix(b).subs(sub);par=[]
    for j,f in [(3,lambda y:2*y[5]/(3*y[12])),(2,lambda y:y[13]/(6*y[12])),(0,lambda y:y[10]/y[12])]:
     tt=s.factor(f(y));par.append(tt);y=flows[j].subs(dict(zip(b,y)),simultaneous=True).subs(t,tt).applyfunc(s.factor)
    assert y[5]==y[10]==y[13]==0 and y[8]==b[8] and y[14]==b[14]+b[5] and y[12]==b[12]
    V=b[8];Z=b[14]+b[5];normal=s.zeros(15,1);normal[3]=-1;normal[8]=normal[12]=normal[14]=1
    for i,f in [(4,b[12]**2/(V*Z**2)),(6,b[12]/Z**2),(7,b[12]**2/(V**2*Z)),(9,b[12]/V**2),(11,b[12]/(V*Z))]:normal[i]=s.factor(y[i]*f)
    Ds=D.subs(sub);assert Ds[8,[0,1,2,3,5]]==s.Matrix([[0,0,0,0,b[8]]]);assert(Ds[14,:]+Ds[5,:])[:,[0,1,2,3,5]]==s.Matrix([[0,0,0,0,-Z]])
    a=s.symbols('a0:5');bn=s.Matrix([0,0,0,0,-1,a[0],0,a[1],a[2],1,a[3],0,a[4],1,0,1]);u,act=beta_actions(m);S=bn.jacobian(a).row_join(act.subs(dict(zip(u,bn)))).row_join(bn);wit=minor_witness(S.subs(dict(zip(a,range(1,6)))));assert wit['rank']==16

    # Complete the finite normalization before restricting b2.
    prior=json.loads((ROOT/'reports/S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101.json').read_text());flow4=s.Matrix(prior['finite_beta_flows']['4']).applyfunc(s.sympify);h12={b[0]:0,b[1]:0,b[3]:-b[12]};first=flow4.subs(h12).subs(t,-b[2]/(3*b[12])).applyfunc(s.factor);assert first[0]==first[1]==first[2]==0 and first[3]==-first[12] and first[12]==b[12]
    # Verify the two diagonal group weights and the scalar normalizer exactly.
    rho,sigma,scale=s.symbols('rho sigma scale');g=s.diag(rho,1/rho,sigma,1/sigma);assert(g.inv()*Q.subs(dict(zip(xx,g*xx)),simultaneous=True)*g.inv().T-rho**3*sigma*Q).applyfunc(s.factor)==s.zeros(4);C=g.T*BB.subs(dict(zip(z,g*s.Matrix(z))),simultaneous=True)*g;tb=(iv*s.Matrix([C[i,j].diff(z[k])for(i,j),k in COORDS])[rows,:]).applyfunc(s.factor)
    weights=[]
    for i in range(15):
     if y[i]==0:continue
     q1=s.cancel(D[i,1]/b[i]);q5=s.cancel(D[i,5]/b[i]);assert q1.is_Integer and q5.is_Integer
     assert s.factor(tb[i+1]-base[i+1]*rho**q1*sigma**q5)==0
     moved=s.factor(scale*y[i]*(scale*b[12])**q1*(1/(scale*V))**q5);num=s.fraction(s.cancel(moved-normal[i]))[0];assert s.factor(s.rem(num,V*Z*scale**2-1,scale))==0;weights.append([i,str(q1),str(q5)])
    # The three excluded source divisors have proper beta images, including closures.
    u,act=beta_actions(m);bad=[]
    for label,badsub in [('b12',dict(h12,**{str(b[12]):0})),('b8',dict(sub,**{str(b[8]):0})),('b14_plus_b5',dict(sub,**{str(b[14]):-b[5]}))]:
     badsub={s.Symbol(k)if isinstance(k,str)else k:v for k,v in badsub.items()};bb=(H*s.Matrix(b)).subs(badsub,simultaneous=True).subs(badsub,simultaneous=True);ps=sorted(bb.free_symbols,key=str);DD=bb.jacobian(ps).row_join(act.subs(dict(zip(u,bb))));NN=DD.T.to_DM().nullspace().to_Matrix().T;assert NN.cols==1 and(DD.T.to_DM()*NN.to_DM()).is_zero_matrix;sample=dict(zip(ps,range(1,len(ps)+1)));ww=minor_witness(DD.subs(sample));assert ww['rank']==15;bad.append(dict(factor=label,source_beta=list(map(str,bb)),parameters=list(map(str,ps)),beta_differential=[[str(v)for v in row]for row in DD.tolist()],left_kernel=[[str(v)for v in row]for row in NN.tolist()],rank_witness=ww,image_dimension=15))
    deps=['reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_MINIMAL_FLAG_ENDPOINT_INVARIANTS_20261101.json','reports/S61_MINIMAL_FLAG_ENDPOINT_TABLE_20261101.json','code/s61_model.py','code/primitive_jet_elliptic_curve_20261101.py','code/weight1_structure.py','code/s61_minimal_rankone_A1_removal_20261101.py'];paths=[Path(__file__)]+[ROOT/f for f in deps]
    return dict(id='S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101',all_passed=True,finite_flows={str(j):list(map(str,v))for j,v in flows.items()},first_flow4_parameter=str(-b[2]/(3*b[12])),first_flow4_result=list(map(str,first)),unipotent_parameters=list(map(str,par)),unipotent_result=list(map(str,y)),normal_beta15=list(map(str,normal)),diagonal_weights=weights,scalar_relation='b8*(b14+b5)*scale^2=1',normalized_beta16=list(map(str,bn)),slice_parameter_order=list(map(str,a)),dominance_minor=wit,proper_images=bad,scope=['Five-variable normal form for the entire free-M0 flag endpoint stratum over sufficiently general beta.','Uses Sp4 and the existing weighted scaling symmetry; it is not a quotient classification or a saturation-removal theorem.','The other three remaining minimal flag strata are separate.'],proof_note='research/notes/S61_MINIMAL_FREE_FLAG_NORMAL_FORM_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},seconds=time.time()-start)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
