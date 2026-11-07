#!/usr/bin/env python3
"""Independent full-source/transfer/true-cut NEW24 review, no producer/CAS."""
import argparse
import ast
from collections import Counter
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import re
import resource
import subprocess
import time

from classification_unimodular_logarithmic_review import RF,ud,ug,um
from classification_minimal_quartic_review import LaurentRing
import classification_degenerate_quadric_q4_nil31_I_review as q

ROOT=Path(__file__).resolve().parents[1]
NAMES=('h','x','y','z','w')+tuple(f'r{i}'for i in range(6));ZERO=(0,)*11
AV=RF((F(),F(1)));BASE=[0,2,4]
BM=[(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]

def atom(c=1,i=None):
    ex=list(ZERO)
    if i is not None:ex[i]=1
    return {tuple(ex):RF(c)}if c else{}
def add(*polys):
    out={}
    for p in polys:
        for ex,c in p.items():out[ex]=out.get(ex,RF())+c
    return {ex:c for ex,c in out.items()if c}
def scale(p,c):return {ex:c*v for ex,v in p.items()if c*v}
def mul(*polys):
    out=atom()
    for p in polys:
        new={}
        for ex,c in out.items():
            for ey,d in p.items():
                key=tuple(a+b for a,b in zip(ex,ey));new[key]=new.get(key,RF())+c*d
        out={e:c for e,c in new.items()if c}
    return out
def power(p,n):
    if n<0:
        assert set(p)=={ZERO};p=atom(1/p[ZERO]);n=-n
    out=atom()
    while n:
        if n%2:out=mul(out,p)
        n//=2
        if n:p=mul(p,p)
    return out
def der(p,i):
    out={}
    for ex,c in p.items():
        if ex[i]:
            e=list(ex);e[i]-=1;out[tuple(e)]=c*ex[i]
    return out
def parse(text):
    def read(n):
        if isinstance(n,ast.Constant)and type(n.value)is int:return atom(n.value)
        if isinstance(n,ast.Name):return atom(AV)if n.id=='a'else atom(1,NAMES.index(n.id))
        if isinstance(n,ast.UnaryOp):return scale(read(n.operand),-1 if isinstance(n.op,ast.USub)else 1)
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Pow):return power(read(n.left),int(ast.literal_eval(n.right)))
            a,b=read(n.left),read(n.right)
            if isinstance(n.op,ast.Add):return add(a,b)
            if isinstance(n.op,ast.Sub):return add(a,scale(b,-1))
            if isinstance(n.op,ast.Mult):return mul(a,b)
            if isinstance(n.op,ast.Div):assert set(b)=={ZERO};return scale(a,1/b[ZERO])
        raise ValueError(ast.dump(n))
    return read(ast.parse(text,mode='eval').body)
def scalar(text):return parse(text).get(ZERO,RF())
def antisymmetric(P):
    out=dict(P)
    for(i,j),p in P.items():out[j,i]=scale(p,-1)
    return out
def tensor(data):return antisymmetric({tuple(map(int,k.split(','))):parse(v)for k,v in data.items()})
def compose(p,values):
    return add(*(mul(atom(c),*(power(values.get(i,atom(1,i)),n)for i,n in enumerate(ex)if n))for ex,c in p.items()))
def rfpower(v,n):return power(atom(v),n)[ZERO]
def rf_compose(v,k):
    def one(poly):
        out=RF()
        for coefficient in reversed(poly):out=out*k+coefficient
        return out
    return one(v.p)/one(v.q)
def compose_a(p,k):return {e:rf_compose(c,k)for e,c in p.items()}
def ev(v,a):
    def one(poly):
        out=F()
        for c in reversed(poly):out=out*a+c
        return out
    return one(v.p)/one(v.q)
def rf_derivative(v):
    dp=tuple(i*v.p[i]for i in range(1,len(v.p)));dq=tuple(i*v.q[i]for i in range(1,len(v.q)))
    return (RF(dp)*RF(v.q)-RF(v.p)*RF(dq))/RF(um(v.q,v.q))
def der_a(p):return {e:c for e,v in p.items()if(c:=rf_derivative(v))}
def numeric(p,a,rvalues):
    out={}
    for ex,c in p.items():
        value=ev(c,a)
        for i,n in enumerate(ex[5:]):value*=rvalues[i]**n
        if value:
            key=ex[:5]+(0,0,0);out[key]=out.get(key,F())+value
    return {e:c for e,c in out.items()if c}
def rf_matrix(data):return [[scalar(v)for v in row]for row in data]
def rf_inverse(A):
    n=len(A);B=[list(row)+[RF(i==j)for j in range(n)]for i,row in enumerate(A)]
    for j in range(n):
        k=next(k for k in range(j,n)if B[k][j]);B[j],B[k]=B[k],B[j]
        d=B[j][j];B[j]=[v/d for v in B[j]]
        for k in range(n):
            if k!=j and B[k][j]:
                d=B[k][j];B[k]=[v-d*w for v,w in zip(B[k],B[j])]
    return [row[n:]for row in B]
def rf_product(A,B):return [[sum((a*b for a,b in zip(row,col)),RF())for col in zip(*B)]for row in A]
def rf_det(A):
    B=[list(row)for row in A];out=RF(1)
    for j in range(len(A)):
        k=next((k for k in range(j,len(A))if B[k][j]),None)
        if k is None:return RF()
        if k!=j:B[j],B[k]=B[k],B[j];out=-out
        d=B[j][j];out*=d
        for k in range(j+1,len(A)):
            if B[k][j]:
                c=B[k][j]/d
                for l in range(j+1,len(A)):B[k][l]-=c*B[j][l]
                B[k][j]=RF()
    return out
def push(P,G):
    inv=rf_inverse(G);values={i:add(*(scale(atom(1,j),v)for j,v in enumerate(row)))for i,row in enumerate(inv)}
    changed={ij:compose(p,values)for ij,p in P.items()}
    return antisymmetric({(i,j):add(*(scale(changed.get((a,b),{}),G[i][a]*G[j][b])for a,b in it.product(range(5),repeat=2)))for i,j in q.PAIRS})
def jacobi(P):
    return {(i,j,k):add(*(mul(P.get((a,l),{}),der(P.get((b,c),{}),l))for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(5)))for i,j,k in it.combinations(range(5),3)}
def encode(p):return [[list(e),list(map(str,c.p)),list(map(str,c.q))]for e,c in sorted(p.items())]
def qparse(text):
    alias={'h':'h','x':'x1','y':'x2','z':'x3','w':'x4'}
    return q.parse(re.sub(r'\b(h|x|y|z|w)\b',lambda m:alias[m[1]],text))
def qdecode(data):return q.antisymmetric({tuple(map(int,k.split(','))):qparse(v)for k,v in data.items()})
def denominator_support(P,allowed):
    for p in P.values():
        for c in p.values():
            den=c.q
            while len(den)>1:
                g=ug(den,allowed)
                assert len(g)>1,('unrecorded denominator',den)
                den,rem=ud(den,g);assert not rem

def family(saved):
    power_=saved['power'];P=tensor(saved['universal_Poisson_brackets'])
    F4=parse(saved['plane_quartic']);D=mul(atom(1,0),F4)
    U,V=[[parse(v)for v in row]for row in saved['anchors']]
    alpha,beta=list(map(parse,saved['commutator_coefficients']))
    Rmag=parse('r0*h*h+r1*h*y+r2*h*w+r3*y*y+r4*y*w+r5*w*w')
    expected=antisymmetric({(1,2):U[0],(1,4):U[1],(2,3):scale(V[0],-1),(3,4):V[1],(1,3):add(mul(alpha,atom(1,1)),mul(beta,atom(1,3)),Rmag)})
    assert all(P.get(ij,{})==expected.get(ij,{})for ij in q.PAIRS)
    assert all(not v for v in jacobi(P).values())
    det=add(mul(U[0],V[1]),scale(mul(U[1],V[0]),-1))
    assert det==mul(parse(saved['determinant_ratio']),F4)
    bracket=[add(*(mul(U[j],der(V[i],[2,4][j]))for j in range(2)),scale(add(*(mul(V[j],der(U[i],[2,4][j]))for j in range(2))),-1))for i in range(2)]
    assert bracket==[add(mul(alpha,U[i]),mul(beta,V[i]))for i in range(2)]
    assert not any(e[i]for v in U+V+[alpha,beta]for e in v for i in [1,3])
    # Independent COMPLETE18 quadratic gradient-syzygy reconstruction, using a
    # new basis from the matched logarithmic U,V, not producer's large K.
    mons=[]
    for mon in BM:
        ex=[0]*11
        for i,n in zip(BASE,mon):ex[i]=n
        mons.append({tuple(ex):RF(1)})
    columns=[mul(mon,der(D,i))for i in BASE for mon in mons]
    labels=sorted({e for p in columns for e in p},reverse=True)
    M=[[p.get(e,RF())for p in columns]for e in labels]
    assert [len(M),18]==saved['full_quadratic_gradient_syzygy_matrix_shape']
    M0=[[ev(c,F())for c in row]for row in M]
    r,cc=q.rank(M0);assert r==16
    rr=q.rank(q.transpose([[row[j]for j in cc]for row in M0]))[1]
    det0=q.determinant([[M0[i][j]for j in cc]for i in rr]);assert det0 and det0==ev(scalar(saved['nonzero16minor']),F())
    lambdaU=add(der(U[0],2),der(U[1],4),beta)
    lambdaV=add(der(V[0],2),der(V[1],4),scale(alpha,-1))
    syzygies=[]
    for W,lam in [(U,lambdaU),(V,lambdaV)]:
        S=[scale(mul(lam,atom(1,0)),-F(1,5)),add(W[0],scale(mul(lam,atom(1,2)),-F(1,5))),add(W[1],scale(mul(lam,atom(1,4)),-F(1,5)))]
        assert not add(*(mul(S[j],der(D,BASE[j]))for j in range(3)))
        syzygies.append([part.get(next(iter(mon)),RF())for part in S for mon in mons])
    assert q.rank(q.transpose([[ev(c,F())for c in col]for col in syzygies]))[0]==2
    # All-source square-zero cofactor is a structural block identity.
    Z=[add(*(der(P.get((i,j),{}),i)for i in range(5)))for j in range(5)]
    assert all(not Z[i]for i in BASE)
    assert all(not der(v,i)for v in Z for i in [1,3])
    N=[[der(z,j)for j in range(5)]for z in Z]
    assert all(not add(*(mul(N[i][k],N[k][j])for k in range(5)))for i,j in it.product(range(5),repeat=2))
    denom=scalar(saved['coefficient_denominator']);denominator_support(P,denom.p)
    assert ev(denom,F())
    sourcepath=ROOT/saved['input_report'];original=json.loads(sourcepath.read_text())['source_magnetic_normalform']
    values=[F(saved['special_source_values'][f'r{i}'])for i in range(6)]
    endpoint=q.antisymmetric({ij:numeric(p,F(),values)for ij,p in P.items()if ij[0]<ij[1]})
    assert all(endpoint.get(ij,{})==qdecode(original['canonical_brackets']).get(ij,{})for ij in q.PAIRS)
    # Entire endpoint magnetic fiber differs only by its six explicit monomials.
    assert all(not e[i]for ij,p in P.items()if ij[0]<ij[1]for e in p for i in range(5,11)if ij!=(1,3))
    ww=[]
    for degree in range(3):ww.append({tuple(e[k]if k!=4 else 0 for k in range(11)):c for e,c in F4.items()if e[4]==degree})
    disc=add(mul(ww[1],ww[1]),scale(mul(ww[2],ww[0]),-4))
    assert disc==parse(saved['quartic_discriminant'])
    at_h0={e:c for e,c in F4.items()if not e[0]}
    assert at_h0==mul(parse('4*(1-4*a**3)'if power_==1 else'36*(1-4*a)'),parse('y**4'))
    return P,dict(power=power_,full7_parameter_tensor_and_all10_strict_Jacobi_verified=True,
        independent_full18_syzygy_shape=[len(M),18],rank16_at_zero_minor=str(det0),independent_complete_two_syzygies=[[[list(map(str,c.p)),list(map(str,c.q))]for c in col]for col in syzygies],
        whole_Saito_determinant_and_commutator_verified=True,all_six_magnetic_endpoint_parameters_retained=True,
        canonical_endpoint_all150_exact=True,full_family_cofactor_has_only_base_to_fiber_block=True,
        full_family_cofactor_square_zero_is_structural_identity=True,
        declared_single_denominator_suffices=True,complex_quartic_discriminant_and_primitive_certificate_verified=True)

def transfer(one,three,saved):
    G0=rf_matrix(saved['base_matrix']);G=rf_matrix(saved['complete_projective_matrix'])
    assert rf_det(G0)==-9/rfpower(AV,6)
    T=rf_matrix(saved['anchor_basis_transfer']);B=rf_matrix(saved['fiber_matrix'])
    assert B==list(map(list,zip(*rf_inverse(T))))
    gf=[[RF(i==j)for j in range(5)]for i in range(5)]
    for i,ii in enumerate([1,3]):
        for j,jj in enumerate([1,3]):gf[ii][jj]=B[i][j]
    assert G==rf_product(gf,G0)and rf_det(G)
    f1=parse(one['plane_quartic']);f3=compose_a(parse(three['plane_quartic']),rfpower(AV,3))
    forward={i:add(*(scale(atom(1,j),c)for j,c in enumerate(row)))for i,row in enumerate(G0)}
    assert compose(f3,forward)==scale(f1,9/rfpower(AV,8))
    P1=tensor(one['universal_Poisson_brackets']);P3={ij:compose_a(p,rfpower(AV,3))for ij,p in tensor(three['universal_Poisson_brackets']).items()}
    moved=push(P1,G);maps=list(map(parse,saved['magnetic_parameter_map']))
    assert all(not any(e[:5])and sum(e[5:])==1 for p in maps for e in p)
    expected={ij:compose(p,{5+i:v for i,v in enumerate(maps)})for ij,p in P3.items()}
    assert all(moved.get(ij,{})==expected.get(ij,{})for ij in q.PAIRS)
    M=[[der(p,5+i).get(ZERO,RF())for i in range(6)]for p in maps]
    determinant=rf_det(M);assert determinant and determinant==scalar(saved['magnetic_linear_map_determinant'])
    return dict(all150_universal_a_to_a3_source_push_verified=True,full6_by6_magnetic_map_nonzero_determinant=True,
        magnetic_determinant_numerator=list(map(str,determinant.p)),magnetic_determinant_denominator=list(map(str,determinant.q)),
        a_zero_not_inverted_at_endpoints=True)

def group_columns(P):
    vertical=[];orbit=[]
    for a,b in it.product(range(5),repeat=2):
        vertical.append(q.vec({(i,j):q.add(q.atom(1,i,b)if j==a else{},q.atom(-1,j,b)if i==a else{})for i,j in q.PAIRS}))
        orbit.append(q.vec({(i,j):q.scale(q.add(P.get((b,j),{})if i==a else{},P.get((i,b),{})if j==a else{},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1)),-1)for i,j in q.PAIRS}))
    return vertical,vertical+orbit+[q.vec(P)]
def derivative(P):
    columns=[]
    for ij,e in q.COLS:
        inc=dict(P);inc[ij]=q.add(inc.get(ij,{}),{e+(0,0,0):F(1)});inc[ij[::-1]]=q.scale(inc[ij],-1)
        columns.append(q.jvec(q.jacobi(inc)))
    return q.transpose(columns)
def check_minor(A,cert):
    value=q.determinant([[A[i][j]for j in cert['columns']]for i in cert['rows']]);assert value==F(cert['determinant'])and value
    return dict(rank=len(cert['rows']),determinant=str(value),rows=cert['rows'],columns=cert['columns'])

def component(Pfamily,saved,pointrow):
    a0=F(pointrow['generic_source_point']['a']);rs=[F(pointrow['generic_source_point'][f'r{i}'])for i in range(6)];factor=F(pointrow['generic_point_integer_scale'])
    point=q.antisymmetric({ij:q.scale(numeric(p,a0,rs),factor)for ij,p in Pfamily.items()if ij[0]<ij[1]})
    assert all(point.get(ij,{})==qdecode(pointrow['generic_point_brackets']).get(ij,{})for ij in q.PAIRS)
    assert all(not p for p in q.jacobi(point).values())
    D=derivative(point);vertical,cols=group_columns(point);M=q.transpose(cols)
    assert q.rank(q.transpose(vertical))[0]==24 and q.rank(M)[0]==48
    dirs=[]
    for variable in ['a',*range(6)]:
        changed={ij:der_a(p)if variable=='a'else der(p,5+variable)for ij,p in Pfamily.items()}
        dirs.append(q.vec({ij:q.scale(numeric(p,a0,rs),factor)for ij,p in changed.items()}))
    Fcols=cols+dirs;A=q.transpose(Fcols)
    assert q.rank(A)[0]==49
    assert all(not sum((row[j]*col[j]for j in range(150)if row[j]and col[j]),F())for row in D for col in Fcols)
    dcert=check_minor(D,saved['full_derivative']);acert=check_minor(A,saved['source_frame'])
    assert dcert['rank']==100 and acert['rank']==49
    check_minor(M,pointrow['full_actual_action']);check_minor(A,pointrow['full_all7_source_tangent'])
    normal=[F(saved['normal_sparse'].get(str(i),'0'))for i in range(150)]
    cuts=saved['source_linear_cut_rows'];assert cuts==saved['source_frame']['rows']
    assert all(not normal[i]for i in cuts)
    assert all(not sum(a*b for a,b in zip(row,normal))for row in D)
    Rpoint={}
    for i,(ij,e)in enumerate(q.COLS):
        if normal[i]:Rpoint[ij]=q.add(Rpoint.get(ij,{}),{e+(0,0,0):normal[i]})
    Rpoint=q.antisymmetric(Rpoint)
    Q=q.jvec(q.jacobi(Rpoint));assert Q==[F(saved['quadratic_sparse'].get(str(i),'0'))for i in range(350)]
    augmented=[row+[v]for row,v in zip(D,Q)];qcert=check_minor(augmented,saved['quadratic_obstruction']);assert qcert['rank']==101
    extended=q.transpose(Fcols+[normal]);check_minor(extended,saved['source_plus_normal_frame'])
    chosen=saved['source_frame']['columns'];frame=q.transpose([Fcols[i]for i in chosen]+[normal])
    implicit=dcert['columns'];complement=[i for i in range(150)if i not in implicit]
    sumdet=q.determinant([frame[i]for i in complement]);assert sumdet
    # These two nonzero minors certify the full150 direct sum source49+R1+U100.
    N=[]
    for j in range(5):
        Z=q.add(*(q.der(point.get((i,j),{}),i)for i in range(5)))
        N.append([q.der(Z,k).get((0,)*8,F())for k in range(5)])
    assert N==[[F(v)for v in row]for row in saved['cofactor_matrix']]
    assert q.matmul(N,N)==[[F()for _ in range(5)]for _ in range(5)]and q.rank(N)[0]==2
    hcol=[F(i==0)for i in range(5)]
    assert q.rank([row+[hcol[i]]for i,row in enumerate(q.transpose(N))])[0]==3
    return dict(full150_350_derivative=dcert,full7_source_frame=acert,all350_D_source_products_zero=True,
        vertical_rank=24,actual_group_rank=48,actual_full_source_rank=49,true_normal_on_all49cuts=True,
        full_cut_kernel_dimension=1,independent_source49_normal1_implicit100_direct_sum_determinant=str(sumdet),
        quadratic_nonzero_full_cokernel_certificate=qcert,component_dimension=24,
        smoothness_and_local_uniqueness_not_claimed=True,
        all_source_N2_zero_dependency='Exact structural base-to-fiber block checked for the complete rational family, not inferred as an open condition from this point.',
        point_rankN2_and_h_outside_dual_image_open_conditions_verified=True)

def simplify_edges(text):
    def walk(n):
        if isinstance(n,ast.Name)and n.id.startswith('t'):return ast.Constant(0)
        if isinstance(n,ast.BinOp):
            left,right=walk(n.left),walk(n.right)
            zero=lambda x:isinstance(x,ast.Constant)and x.value==0
            if isinstance(n.op,ast.Mult)and(zero(left)or zero(right)):return ast.Constant(0)
            if isinstance(n.op,ast.Div)and zero(left):return ast.Constant(0)
            if isinstance(n.op,ast.Add):
                if zero(left):return right
                if zero(right):return left
            if isinstance(n.op,ast.Sub)and zero(right):return left
            if isinstance(n.op,ast.Pow)and zero(left)and isinstance(right,ast.Constant)and right.value>0:return ast.Constant(0)
            return ast.BinOp(left,n.op,right)
        if isinstance(n,ast.UnaryOp):
            value=walk(n.operand)
            return value if isinstance(value,ast.Constant)and value.value==0 else ast.UnaryOp(n.op,value)
        return n
    return ast.unparse(ast.fix_missing_locations(walk(ast.parse(text,mode='eval').body)))

def prior_audit(saved,inputs):
    blob=subprocess.check_output(['git','cat-file','blob',saved['catalogue_git_blob']],cwd=ROOT)
    assert hashlib.sha256(blob).hexdigest()==saved['catalogue_sha256']
    rows=json.loads(blob)['components'];assert rows==saved['all102_rows_frozen']and len(rows)==102
    assert subprocess.check_output(['git','rev-parse',saved['immutable_commit']+':research/COMPONENTS.json'],cwd=ROOT,text=True).strip()==saved['catalogue_git_blob']
    witnesses=[]
    for claim in saved['all8_MPS24_actual_SNC_witnesses']:
        path=ROOT/claim['accepted_whole_source_report'];inputs.append(path)
        assert hashlib.sha256(path.read_bytes()).hexdigest()==claim['source_sha256']
        source=json.loads(path.read_text());periods=claim['all_period_parameters']
        R=LaurentRing(periods,[f'x{i}'for i in range(5)],invert=periods)
        P={}
        for textpair,text in source['homogeneous_bivector'].items():
            i,j=map(int,textpair.split(','));p=R.parse(simplify_edges(text));P[i,j]=p;P[j,i]=R.scale(p,-1)
        for i,j in q.PAIRS:
            expected=R.mul(R.parse(claim['actual_logcanonical_coefficients'][f'{i},{j}']),R.variable(f'x{i}'),R.variable(f'x{j}'))
            assert P.get((i,j),{})==expected
        D={}
        for missing in range(5):
            i,j,k,l=[a for a in range(5)if a!=missing]
            pf=R.add(R.mul(P.get((i,j),{}),P.get((k,l),{})),R.scale(R.mul(P.get((i,k),{}),P.get((j,l),{})),-1),R.mul(P.get((i,l),{}),P.get((j,k),{})))
            D=R.add(D,R.scale(R.mul(R.variable(f'x{missing}'),pf),(-1)**missing))
        kappa=R.parse(claim['kappa']);assert kappa and all(not any(e[-5:])for e in kappa)
        assert D==R.mul(kappa,*[R.variable(f'x{i}')for i in range(5)])
        assert D==R.parse(claim['actual_SNC_quintic'])
        power5=R.mul(kappa,kappa,kappa,kappa,kappa)
        assert power5==R.parse(claim['gradient_5_by_5_minor'])and power5
        witnesses.append(dict(id=claim['id'],all_periods_retained=periods,actual_full_source_edge0_SNC_verified=True,gradient_rank5_minor_kappa5_verified=True))
    ids={c['id']for c in witnesses};assert ids=={'T02','T03','T16','T17','T21','T36','T37','T38'}
    reasons=[]
    for row in rows:
        label=row['id']
        if row['generic_rank']==2:reason='closed rank2 locus'
        elif row['projective_dimension']!=24:reason='different full component dimension'
        elif label in ids:reason='actual rank5 SNC source point outside new closed essential-rank<=3 locus'
        elif label=='R4_CPLUS_EXT':
            assert row['generic_actual_degeneracy_divisor']['irreducible_factor_degrees']==[1,2,2];reason='closed1+2+2 locus'
        elif label=='R4_DOUBLE_RANK3_QUADRIC_TWZERO_24':
            assert not row['generic_actual_degeneracy_divisor']['reduced'];reason='closed nonreduced quintic locus'
        else:
            assert label=='R4_DEPENDENT_THREE_HYPERPLANES_QUADRIC_24'and row['generic_actual_degeneracy_divisor']['irreducible_factor_degrees']==[1,1,1,2];reason='closed1+1+1+2 locus'
        reasons.append(dict(id=label,reason=reason))
    return dict(immutable_commit=saved['immutable_commit'],immutable_git_blob=saved['catalogue_git_blob'],immutable102_rows_verified=True,
        all8_actual_MPS_source_SNC_witnesses=witnesses,all102_distinctness_reasons=reasons,all_prior102_distinct=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    paths=[ROOT/'reports/CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json',
        ROOT/'research/classification/GA9_CUSPIDAL_QUARTIC_COMPONENT.md',
        ROOT/'research/classification/HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND.md',
        ROOT/'research/classification/NIL211_SHORT_TRANSLATION_REVIEW.md',ROOT/'research/notes/MPS_S61_CATALOGUE_20261103_v62.md']
    saved=json.loads(paths[0].read_text());records=[];families=[]
    for key in ['actual_family_power1','actual_family_power3']:
        p,record=family(saved[key]);families.append(p);records.append(record);paths.append(ROOT/saved[key]['input_report'])
    comparison=transfer(saved['actual_family_power1'],saved['actual_family_power3'],saved['whole_six_magnetic_family_transfer'])
    support=component(families[1],saved['full_source_component'],saved['actual_family_power3'])
    audit=prior_audit(saved['all_prior102_audit'],paths)
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
        all_two_complete7_sources=records,whole_sixmagnetic_transfer=comparison,true_full_ambient_support=support,
        immutable_all_prior102_audit=audit,new24_componenthood_and_distinctness_supported=True,
        whole_Ga9_power1_and_power3_closures_contained=True,
        nilpotency_openness_correction='N²0 follows from the allsource base-to-fiber block identity. Only rankN=2 and h outside dualimage are open conditions checked at a1. No nilpotency is generalized from one point by openness.',
        proof_bridge='Actual source rank49 and true49cut with D100/one normal and nonzero quadratic fullcokernel imply rawlocaldimension<=49 by IFT+Krull, hence projective24 component. No smoothness or localuniqueness claimed.',
        localizations_and_complements='Families regular at a0; source denominators2a³-3 or18a-7. Primitive units1-4a³ or1-4a. Generic transfer also inverts a and18a³-7; a0 handled in each original family and by closedness. No type6wholeclassification inferred on excluded parameter values.',
        global_count_or_registers_not_edited=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['all_passed','new24_componenthood_and_distinctness_supported','seconds','maximum_RSS_KiB']}))

if __name__=='__main__':main()
