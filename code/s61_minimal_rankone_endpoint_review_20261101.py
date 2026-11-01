#!/usr/bin/env python3
"""Independent raw135 review of the rank-one A1 binary endpoint."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda x:s.Matrix(x).applyfunc(s.sympify)


def run():
    start=time.time();rp=ROOT/'reports/S61_MINIMAL_RANKONE_ENDPOINT_20261101.json';r=json.loads(rp.read_text());assert r['all_passed']
    for p,h in r['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==2]
    qi,ui={i:j for j,i in enumerate(two)},{i:j for j,i in enumerate(one)}
    eqs=model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    encoded=[[(tuple(qi[i]for i in mon if i in qi),tuple(ui[i]for i in mon if i in ui),s.Rational(a,b))for mon,a,b in eq]for eq in eqs]
    U,G,local=r['universal_double_line'],r['source_geometry'],r['point_local_data'];b=s.symbols('b0:9');z=s.Symbol('s')
    beta,q0,delta,n=M(U['beta']),M(U['q0']),M(U['delta']),M(U['n']);qb=q0+z*delta
    assert (qb-M(U['qbase'])).applyfunc(s.cancel)==s.zeros(45,1)
    denb,denq,denn=map(s.sympify,[U['beta_den'],U['q_den'],U['n_den']]);R,*_=ring(list(b)+[z],QQ)
    BP=[R.from_expr(s.cancel(v*denb))for v in beta];QP=[R.from_expr(s.cancel(v*denq))for v in qb];NP=[R.from_expr(s.cancel(v*denn))for v in n];DB=R.from_expr(denb)
    for row,eq in enumerate(encoded):
        values=[R.zero,R.zero]
        for qm,bm,c in eq:
            if len(qm)!=2:continue
            a,d=qm;C=c*DB**(2-len(bm))
            for j in bm:C*=BP[j]
            values[0]+=C*QP[a]*QP[d];values[1]+=C*(QP[a]*NP[d]+NP[a]*QP[d])
        assert all(v==0 for v in values),row
    fixed=next(i for i,v in enumerate(q0)if v);assert q0[fixed]==15 and n[fixed]==delta[fixed]==0
    n0=n.subs(z,0).applyfunc(s.cancel);assert(n0-M(U['n0'])).applyfunc(s.cancel)==s.zeros(45,1)
    minor=s.factor(delta.row_join(n0)[[0,34],:].det());assert minor==s.sympify(U['embedding_minor'])
    env={s.Symbol(k):s.Rational(v)for k,v in r['source_sample'].items()}
    point=M(local['point61']);assert beta.subs(env).applyfunc(s.cancel)==point[one,:] and q0==point[two,:]
    assert int(s.Rational(minor.subs(env)).p)%32003
    for value in list(n)+list(qb):
        num,den=s.fraction(s.cancel(value.subs(env)));assert int(s.Rational(den.subs(z,0)).p)%32003
        assert all(int(s.Rational(c).q)%32003 for c in s.Poly(num,z).all_coeffs()+s.Poly(den,z).all_coeffs())
    chart=[i for i in range(45)if i!=fixed];D=s.zeros(135,45);tc=s.zeros(135,1);polys=[];cache={}
    for row,eq in enumerate(encoded):
        f={}
        for qm,bm,c in eq:
            if bm not in cache:cache[bm]=s.prod(point[one[k]]for k in bm)
            co=c*cache[bm];f[qm]=f.get(qm,0)+co
            if len(qm)==2:
                i,j=qm;D[row,i]+=co*q0[j];D[row,j]+=co*q0[i]
            elif len(qm)==1:tc[row]+=co*q0[qm[0]]
        assert sum(c*s.prod(q0[j]for j in qm)for qm,c in f.items()if len(qm)==2)==0
        polys.append(f)
    A=D[:,chart].row_join(tc);T,L=M(local['tangent_basis']),M(local['cokernel_rows']);assert A==M(local['jacobian'])and A*T==s.zeros(135,10)and L*A==s.zeros(L.rows,45)
    assert T[-1,:]==s.zeros(1,10);free=[i for i in range(45)if i not in A.rref()[1]];assert T[free,:]==s.eye(10)
    forms={i:([0]*10 if i==fixed else list(T[chart.index(i),:]))for i in range(45)};monos=list(combinations_with_replacement(range(10),2));indices={m:i for i,m in enumerate(monos)};Q=s.zeros(135,55)
    for row,f in enumerate(polys):
        for qm,c in f.items():
            ff=[forms[i]for i in qm]+[list(T[-1,:])]*(2-len(qm))
            for i in range(10):
                for j in range(10):Q[row,indices[tuple(sorted((i,j)))]]+=c*ff[0][i]*ff[1][j]
    obstruction=L*Q;assert obstruction==M(local['quadratic_obstruction_basis']);hs=[]
    for coeff in obstruction.tolist():
        H=s.zeros(10)
        for c,(i,j)in zip(coeff,monos):H[i,j]+=c;H[j,i]+=c
        hs+=H.tolist()
    ranks=(D.rank(),A.rank(),obstruction.rank(),s.Matrix(hs).rank());assert ranks==(34,35,5,7)
    tangents=delta.row_join(n0).subs(env)[chart,:].col_join(s.zeros(1,2));coords=tangents[free,:]
    assert T*coords==tangents and coords==M(r['actual_double_line_tangent'])and coords[6,0]==0 and coords[6,1]==s.Rational(593,972)
    ann=M(r['actual_double_line_initial_annihilator']);assert ann.rank()==8 and ann*coords==s.zeros(8,2)
    zz,u,B=primitive_beta(model);aa,cc,dd=s.symbols('aa cc dd');v1=s.Matrix([aa,1,0,cc]);v2=s.Matrix([-cc,0,1,dd]);BB=B.subs(dict(zip(u,point[one,:])))
    section=s.Matrix([(v1.T*BB.subs(dict(zip(zz,v)))*v2)[0]for v in [v1,v2]]);CJ=section.jacobian([aa,cc,dd]).subs({aa:0,cc:0,dd:0});assert CJ==M(r['curve_jacobian'])and CJ.rank()==2
    family=M(r['actual_binary_family_tangent']);assert s.Matrix(hs)*family==s.zeros(50,3)and family[r['binary_parameters'],:].det()!=0
    vec=s.Matrix([B[i,j].diff(zz[k])for(i,j),k in COORDS]);basis=vec.jacobian(u);rows=list(basis.T.rref()[1]);iv=basis[rows,:].inv()
    old=json.loads((ROOT/'reports/S61_MINIMAL_RANKONE_JUNCTION_20261101.json').read_text())['A1_universal'];src=M(old['beta_kernel'])*s.Matrix(b);assert src==M(U['source_beta'])
    P=B.subs(dict(zip(u,src)));g=M(U['endpoint_normalization']);assert(g*J[1:,1:]*g.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
    pulled=g.T*P.subs(dict(zip(zz,g*s.Matrix(zz))),simultaneous=True)*g
    assert(iv*s.Matrix([pulled[i,j].diff(zz[k])for(i,j),k in COORDS])[rows,:]-beta).applyfunc(s.cancel)==s.zeros(16,1)
    x,y,z0,w=X[1:];xv=s.Matrix([x,y,z0,w]);h=s.Symbol('h');dh=1+b[7]*(h-1)/(2*b[6]);shape=2*b[7]**2*(h-1)/(9*b[6]**2)
    P=s.zeros(4);P[0,1]=x*x;P[1,2]=-shape*x*x-2*dh*x*w-3*h*w*w;P-=P.T
    Q0=s.zeros(4);Q0[1,2]=w*w+15*x*x;Q0-=Q0.T
    assert(g.inv()*P.diff(h).subs(dict(zip(xv,g*xv)),simultaneous=True)*g.inv().T+54*b[6]/b[7]*Q0).applyfunc(s.cancel)==s.zeros(4)
    parameters=list(map(s.Symbol,G['parameter_order']));target=M(G['beta']);inverse={s.Symbol(k):s.sympify(v)for k,v in G['source_inverse'].items()};assert(beta.subs(inverse)-target).applyfunc(s.cancel)==s.zeros(16,1)
    H=target.jacobian(parameters);act=s.Matrix.hstack(*[iv*infinitesimal_action(vec,Bb)[rows,:]for Bb in sp_basis()]).subs(dict(zip(u,target)));O=M(G['projective_orbit']);st=M(G['projective_stabilizer']);an=s.Matrix.hstack(*H.T.nullspace()).T
    assert H.rank()==9 and O.rank()==5 and O*st==s.zeros(45,5)
    normal=(an*act*st).applyfunc(s.cancel);assert(normal-M(G['stabilizer_normal_matrix'])).applyfunc(s.cancel)==s.zeros(7,5)and normal.to_DM().rank()==2
    sample={parameters[0]:2,parameters[1]:3,**dict(zip(parameters[2:],[2,3,5,7,11,13,5]))};diff=H.row_join(act);assert diff.subs(sample).rank()==16 and diff.col_join(s.zeros(45,9).row_join(O)).subs(sample).rank()==16
    paths=[Path(__file__),rp,ROOT/'code/s61_minimal_rankone_endpoint_20261101.py',ROOT/'research/notes/S61_MINIMAL_RANKONE_ENDPOINT_20261101.md',ROOT/'research/notes/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.md',ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json',ROOT/'reports/S61_MINIMAL_RANKONE_JUNCTION_20261101.json',ROOT/'reports/S61_MINIMAL_RANKONE_A1_REMOVAL_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_RANKONE_ENDPOINT_REVIEW_20261101',all_passed=True,raw135_all270_coefficients_verified=True,fixed_chart_value=15,raw_point_ranks=list(ranks),flat_double_line_embedding=True,same_binary_smooth_chart=True,finite_endpoint_normalization=True,source_beta_rank=16,source_joint_rank=16,conceptual_review=dict(full_A1_curve_endpoint_coverage_accepted=True,even_Laurent_descent_accepted=True,strict_filtered_kernel_DVR_accepted=True,relative_ambient_t2_and_generic_Z16_removal_accepted=True),scope=['Independent raw135 all-source double-line coefficients, local derivatives/quadratics and finite endpoint normalization.','Primary sextic closure and Macaulay2 decomposition reproduced separately; no independent implementation of that elimination.','Even-Laurent construction is source-audited and replayed; final rational formulas independently satisfy all270 identities.','No free66, remaining boundary exhaustion, or affine exhaustion assertion.'],seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
