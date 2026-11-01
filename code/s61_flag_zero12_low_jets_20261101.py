#!/usr/bin/env python3
"""Universal true cubic initial relations on the seven-variable flag chart."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)
enc=lambda A:[[str(v)for v in row]for row in A.tolist()]


def run():
    start=time.time();model=S61();source=ROOT/'reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json';nf=json.loads(source.read_text());assert nf['all_passed']
    for p,h in nf['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=json.loads((ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json').read_text());H=M(old['flag']['beta_kernel']);b=s.symbols('b0:15');parameters=[b[i]for i in [4,7,8,9,10,11,14]]
    sub={b[0]:1,b[1]:0,b[2]:1,b[3]:0,b[5]:0,b[6]:0,b[12]:0,b[13]:0};one=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==2]
    beta=(H*s.Matrix(b)).subs(sub);q0=M(old['flag']['point61'])[two,:];fixed=two.index(37);chart=[i for i in range(45)if i!=fixed];_,F=specialize(model,dict(zip(one,beta)),True);J=s.zeros(135,45)
    for row,f in enumerate(F):
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon
                if i!=fixed:J[row,chart.index(i)]+=c*q0[j]
                if j!=fixed:J[row,chart.index(j)]+=c*q0[i]
            elif len(mon)==1:J[row,44]+=c*q0[mon[0]]
    rr,C=J.to_DM().rref();rr=rr.to_Matrix();_,rows=J[:,list(C)].T.to_DM().rref();C,rows=list(C),list(rows);block=J[rows,C].to_DM();det=s.factor(block.domain.to_sympy(block.det()));assert det==-s.Rational(1400832,25)
    assert len(C)==36 and J[:,:44].to_DM().rank()==35
    free=[i for i in range(45)if i not in C];T=s.zeros(45,9)
    for j,f in enumerate(free):
        T[f,j]=1
        for i,c in enumerate(C):T[c,j]=-rr[i,f]
    assert J*T==s.zeros(135,9)and T[-1,:]==s.zeros(1,9)
    inv=block.to_field().inv().to_Matrix();PR,*_=ring(parameters,QQ);convert=lambda x:PR.from_expr(s.cancel(x));JP=[[convert(v)for v in row]for row in J.tolist()];UP=[[PR.zero]*36 for _ in range(45)]
    for i,c in enumerate(C):UP[c]=[convert(v)for v in inv[i,:]]
    W={1:[{(j,):convert(v)for j,v in enumerate(T[i,:])if v}for i in range(45)]};terms=[]
    for f in F:
        out=[]
        for mon,c in f.items():
            if fixed in mon:continue
            i,j=[chart.index(k)for k in mon]+[44]*(2-len(mon));cc=convert(c)
            if cc:out.append((i,j,cc))
        terms.append(out)
    def clean(f):return {m:c for m,c in f.items()if c}
    def add_product(out,f,g,c):
        for a,ca in f.items():
            for d,cd in g.items():
                mon=tuple(sorted(a+d));out[mon]=out.get(mon,PR.zero)+c*ca*cd
    def combine(matrix,polys):
        result=[]
        for row in matrix:
            out={}
            for c,f in zip(row,polys):
                if c:
                    for mon,v in f.items():out[mon]=out.get(mon,PR.zero)+c*v
            result.append(clean(out))
        return result
    residual={}
    for degree in [2,3]:
        Q=[]
        for ff in terms:
            out={}
            for da in range(1,degree//2+1):
                db=degree-da
                for i,j,c in ff:
                    add_product(out,W[da][i],W[db][j],c)
                    if da!=db:add_product(out,W[db][i],W[da][j],c)
            Q.append(clean(out))
        W[degree]=combine([[-v for v in row]for row in UP],[Q[i]for i in rows]);correction=combine(JP,W[degree]);residual[degree]=[]
        for a,c in zip(Q,correction):
            out=dict(a)
            for mon,v in c.items():out[mon]=out.get(mon,PR.zero)+v
            residual[degree].append(clean(out))
        assert all(not residual[degree][i]for i in rows)
    assert not W[2][-1]
    m2=list(combinations_with_replacement(range(9),2));m3=list(combinations_with_replacement(range(9),3));ix3={m:i for i,m in enumerate(m3)}
    A2=s.Matrix([[f.get(mon,PR.zero).as_expr()for mon in m2]for f in residual[2]]);_,qrows=A2.T.to_DM().rref();qrows=list(qrows);B=A2[qrows,:];assert len(qrows)==8 and all(v.is_Rational for v in B)
    _,qcols=B.to_DM().rref();qcols=list(qcols);qinv=B[:,qcols].inv();mult=(A2[:,qcols]*qinv).applyfunc(s.expand);assert (A2-mult*B).applyfunc(s.expand)==s.zeros(135,45)
    active={4,5,7,8};assert all(c==0 or set(mon)<=active for row in B.tolist()for mon,c in zip(m2,row));hs=[]
    for row in B.tolist():
        HH=s.zeros(9)
        for c,(i,j)in zip(row,m2):HH[i,j]+=c;HH[j,i]+=c
        hs+=HH.tolist()
    assert s.Matrix(hs).rank()==4
    multiples=[]
    for row in B.tolist():
        for j in range(9):
            v=[s.S.Zero]*165
            for mon,c in zip(m2,row):v[ix3[tuple(sorted(mon+(j,)))]]+=c
            multiples.append(v)
    mr,piv=s.Matrix(multiples).rref();assert len(piv)==60;nonp=[i for i in range(165)if i not in piv]
    def cubic_class(f):
        coeff=[f.get(mon,PR.zero).as_expr()for mon in m3]
        for i,p in enumerate(piv):
            c=coeff[p]
            for j in nonp:
                if mr[i,j]:coeff[j]-=c*mr[i,j]
        return [s.expand(coeff[j])for j in nonp]
    A3=s.Matrix([cubic_class(f)for f in residual[3]])
    # Remove the quadratic leading rows BEFORE treating any cubic as an initial form.
    P=(A3-mult*A3[qrows,:]).applyfunc(s.expand);assert P.to_DM().rank()==7
    prows=[8,9,27,28,67,69,134]
    column_sets=[[40,63,64,83,84,96,97],[63,83,84,96,97,103,104],[67,83,84,96,97,103,104]]
    determinants=[]
    for cs in column_sets:
        dm=P[prows,cs].to_DM();determinants.append(s.factor(dm.domain.to_sympy(dm.det())))
    assert s.cancel(determinants[0]/b[8]**4).is_Rational and determinants[0]!=0
    assert s.cancel(determinants[1].subs(b[8],0)/b[11]).is_Rational and determinants[1].subs(b[8],0)!=0
    assert determinants[2].subs({b[8]:0,b[11]:0}).is_Rational and determinants[2].subs({b[8]:0,b[11]:0})!=0
    t3=s.Matrix([cubic_class(W[3][-1])]);PB=P[prows,:];cs=column_sets[0];coeff=(t3[:,cs].to_DM().to_field()*PB[:,cs].to_DM().to_field().inv()).to_Matrix().applyfunc(s.cancel)
    assert all(c.is_Rational for c in coeff)and(t3-coeff*PB).applyfunc(s.cancel)==s.zeros(1,105)

    # The other chart's candidate b2=2,b8=0 has a proper beta image.
    I=b[0]*b[2]-b[3]*b[12]+b[3]**2;K=b[0]*b[8]-(b[3]+b[12])*b[5]
    before=s.sympify(nf['regular_chart']['before_torus'][8]);relation={b[2]:(2*b[12]**2+b[3]*b[12]-b[3]**2)/b[0]}
    assert s.cancel((before-K/b[0]).subs(relation))==0
    zz,uu,bb=primitive_beta(model);vec=s.Matrix([bb[i,j].diff(zz[k])for(i,j),k in COORDS]);basis=vec.jacobian(uu);br=list(basis.T.rref()[1]);iv=basis[br,:].inv();action=s.Matrix.hstack(*[iv*infinitesimal_action(vec,S)[br,:]for S in sp_basis()]).subs(dict(zip(uu,H*s.Matrix(b))))
    ann=s.zeros(4,16);ann[0,0]=ann[1,2]=1
    for j,variable in enumerate(b):ann[2,j+1]=s.diff(I-2*b[12]**2,variable);ann[3,j+1]=s.diff(K,variable)
    case={b[0]:1,b[1]:0,b[2]:2,b[3]:0,b[5]:0,b[6]:0,b[8]:0,b[12]:1,b[14]:0};N=(ann*action).subs(case).applyfunc(s.expand)
    assert ann.subs(case).rank()==4 and N[0,:]==N[1,:]and N.to_DM().rank()==3
    sparse=lambda mats:{str(d):[[[list(mon),str(c.as_expr())]for mon,c in f.items()]for f in mat]for d,mat in mats.items()}
    paths=[Path(__file__),source,ROOT/'reports/S61_FLAG_UNIFORM_OPEN_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_FLAG_ZERO12_LOW_JETS_20261101',all_passed=True,parameters=list(map(str,parameters)),substitutions={str(k):str(v)for k,v in sub.items()},full_rank=36,leading_rank=35,linear_minor='-1400832/25',pivot_rows=rows,pivot_columns=C,free_coordinates=free,tangent=enc(T),quadratic_rank=8,quadratic_basis_rows=qrows,quadratic_basis=enc(B),quadratic_common_radical=[0,1,2,3,6],quadratic_multiples_rank=60,cubic_rank=7,cubic_quotient_monomials=[list(m3[j])for j in nonp],true_cubic_classes=enc(P),cubic_minor_rows=prows,cubic_minor_columns=column_sets,cubic_minors=list(map(str,determinants)),cubic_minor_cover='b8 !=0; b8=0,b11 !=0; b8=b11=0',t_order_at_least=4,t3_relation=list(map(str,coeff)),lift=sparse(W),residual=sparse(residual),regular_chart_rankdrop_geometry=dict(source_equations=['u4=0','b1=0','I-2*b12^2=0','b0*b8-(b3+b12)*b5=0'],open='b0*b12 !=0',source_dimension=12,normal_matrix=enc(N),normal_rank_bound=3,beta_image_dimension_bound=15),limitations=['Exact only through cubic order; no t-nilpotence or boundary saturation conclusion.','The b12=0 chart is reached through the proved finite normalization after excluding b0*I0=0 as stated there.','The three cubic minors cover all seven parameters, including b8=0.'],proof_note='research/notes/S61_FLAG_ZERO12_LOW_JETS_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
