#!/usr/bin/env python3
"""Independent HF8 whole-source, finite-slice and T36/T38 review.

Standard-library Fraction sparse polynomials only. No producer or CAS import.
The existing independent H arithmetic supplies rational matrices and five-variable
polynomials; the Laurent parser supplies arbitrary parameter polynomial rings.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from pathlib import Path
import ast,hashlib,json,resource,signal,time
import classification_plane_m4_special_review as H
from classification_minimal_quartic_review import LaurentRing

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_m2_eight_hf_separate_review.json'
OUT='reports/CLASSIFICATION_PLANE_M2_EIGHT_HF_SEPARATE_REVIEW_20261108.json'
BASE='reports/CLASSIFICATION_PLANE_DOUBLE_LINE_M2_CONSTANT_DOUBLE_THIRTYFOUR_20261108.json'
SUP='reports/CLASSIFICATION_PLANE_M2_DOUBLE_LINE_EIGHT_HF_SUPPORT_20261108.json'
AMB='reports/CLASSIFICATION_PLANE_M2_HF_TWO_FULL_FAMILY_AMBIENT_20261108.json'
SLICE='reports/CLASSIFICATION_PLANE_M2_HF_COMPLETE_TORUS_SLICES_V2_20261108.json'
PAIRS=H.PAIRS;TRIPLES=H.TRIPLES;QUADS=H.QUADS;COORDS=('a','b','c','u','v')

def load(p):return json.loads((ROOT/p).read_text())
def ck(a,b,label):
    if a!=b:raise AssertionError(label)
def skew(R,d):
    out={(i,j):{} for i in range(5)for j in range(5)}
    for (i,j),p in d.items():out[i,j]=p;out[j,i]=R.scale(p,-1)
    return out
def tensor(R,d):return skew(R,{ast.literal_eval(k):R.parse(v)for k,v in d.items()})
def plus(R,*Ps):return skew(R,{ij:R.add(*(P[ij]for P in Ps))for ij in PAIRS})
def basis(R,v):
    d={ij:{}for ij in PAIRS}
    for (ij,e),c in zip(H.COLS,v):
        if c:d[ij]=R.add(d[ij],{(0,)*R.n+e:F(c)})
    return skew(R,d)
def vec(R,P):return [{e[:R.n]+(0,)*5:c for e,c in P[ij].items()if e[R.n:]==ex}for ij,ex in H.COLS]
def strict(R,P):return {ijk:R.add(*(R.mul(P[i,l],R.derivative(P[j,k],l))for i,j,k in (ijk,(ijk[1],ijk[2],ijk[0]),(ijk[2],ijk[0],ijk[1]))for l in range(5)))for ijk in TRIPLES}
def jac(R,P):
    J=strict(R,P);x=[R.variable(n)for n in R.names[R.n:]]
    return {q:R.add(*(R.scale(R.mul(x[q[k]],J[q[:k]+q[k+1:]]),(-1)**k)for k in range(4)))for q in QUADS}
def jvec(R,J):return [{e[:R.n]+(0,)*5:c for e,c in J[q].items()if e[R.n:]==ex}for q,ex in H.ROWS]
def pfs(R,P):
    out=[]
    for at in range(5):
        i,j,k,l=[v for v in range(5)if v!=at]
        out.append(R.add(R.mul(P[i,j],P[k,l]),R.scale(R.mul(P[i,k],P[j,l]),-1),R.mul(P[i,l],P[j,k])))
    return out
def actual(R,P):return R.add(*(R.scale(R.mul(R.variable(R.names[R.n+i]),v),(-1)**i)for i,v in enumerate(pfs(R,P))))
def zero_params(R,P,names):return {e:c for e,c in P.items()if not any(e[R.names.index(n)]for n in names)}
def evaluate(R,P,point):
    out={}
    for e,c in P.items():
        for n,k in zip(R.names[:R.n],e[:R.n]):c*=F(point.get(n,0))**k
        ex=e[R.n:];out[ex]=out.get(ex,F(0))+c
    return {e:c for e,c in out.items()if c}
def divide(R,P,Q):
    assert Q;out={};rem=P.copy();qe=max(Q);qc=Q[qe]
    while rem:
        re=max(rem);assert all(a>=b for a,b in zip(re,qe)),'nonexact polynomial factor'
        ex=tuple(a-b for a,b in zip(re,qe));co=rem[re]/qc
        out=R.add(out,{ex:co});rem=R.add(rem,R.scale(R.mul({ex:co},Q),-1))
    return out

def source_review(packet,old):
    B1=old['full30_weight1_basis'];B2=old['full9_weight2_basis'];out=[];full=[]
    for rec in packet['cases']:
        name=rec['name'];parent=next(c for c in old['cases']if name.startswith(c['name']+'_WALL_'))
        prior=next(w for w in parent['all_allowed_parameter_walls']if w['name']==name)
        params=tuple(rec['all_P1_parameters'])+tuple(rec['all_homogeneous_P2_parameters_retained'])
        R=LaurentRing(params,COORDS);p=R.parse
        P0=tensor(R,rec['complete_P0']);P1=tensor(R,rec['complete_P1'])
        ck(P0,basis(R,prior['full150_P0']),name+' P0')
        K=H.mat(prior['weight1_certificate']['complete_kernel'])
        cols=[sum(F(B1[i][k])*K[j][i]for i in range(30))for j in range(len(K))for k in range(150)]
        want=[R.add(*(R.scale(p(params[j]),cols[j*150+k])for j in range(len(K))))for k in range(150)]
        ck(vec(R,P1),want,name+' full K1 source')
        qvals=list(map(p,rec['all9_inhomogeneous_P2_plus_full_homogeneous_solution']))
        P2=basis(R,[0]*150)
        P2=skew(R,{ij:R.add(*(R.scale(R.mul(qvals[j],{(0,)*R.n+e:F(B2[j][k])}),1)for j in range(9)for k,(pair,e)in enumerate(H.COLS)if pair==ij and F(B2[j][k])))for ij in PAIRS})
        # The complete affine inverse retains each free coordinate and every
        # homogeneous kernel. The old rank certificate is an accepted input.
        cert=prior['weight2_certificate'];rr=rec['same_existing_M2_minor_rows'];cc=rec['same_existing_M2_minor_columns']
        ck(rr,cert['rows'],name+' P2 rows');ck(cc,cert['columns'],name+' P2 columns')
        ck(F(rec['same_existing_nonzero_M2_minor_determinant']),F(cert['minor_determinant']),name+' P2 det')
        free=[j for j in range(9)if j not in cc]
        ck([qvals[j]for j in free],[p(n)for n in rec['all_homogeneous_P2_parameters_retained']],name+' all P2 free parameters')
        P=plus(R,P0,P1,P2);J=jvec(R,jac(R,P))
        ck(J,list(map(p,rec['all350_original_Jacobi_after_elimination'])),name+' original350')
        # Degree-two compatibility includes P0/P2 and P1 squared. Higher
        # equations and all homogeneous P2 coefficients remain in J.
        weight2=[]
        for (q,e),poly in zip(H.ROWS,J):
            w=e[3]+e[4]-sum(i>=3 for i in q)
            weight2.append(poly if w==2 else {})
        ck(weight2,list(map(p,rec['all350_weight2_compatibilities'])),name+' all weight2 compatibility')
        D=actual(R,P);ck(D,p(rec['full_actual_D']),name+' entire D')
        m=rec['genuine_bad_normal_functional_m'];Dzero=zero_params(R,D,[m])
        assert all(e[R.n+3]>0 for e in Dzero),name+' whole m0 u factor'
        # Every raw P2 variable is a base-pair normal quadratic. No normal
        # row is present even before solving the affine P2 equations.
        assert all(not F(v)for row in B2 for(ij,e),v in zip(H.COLS,row)if 3 in ij or 4 in ij)
        if '_E1_'in name:
            for powerrow in rec['all_exact_original_powers_of_m_rows']:
                ck(J[powerrow['original_row']],p(powerrow['nonzero_constant']+'*'+m+'**'+str(powerrow['power'])),name+' exact square')
            assert rec['all_exact_original_powers_of_m_rows'];result='empty: original square forces m=0, whole u factor'
        elif name.startswith('HF_K0_')or name.startswith('HF_K4_'):
            factor=p('-2*b+p3*u+p6*v'if name.startswith('HF_K0_')else'-6*b+p4*u+p7*v')
            quotient=divide(R,D,factor);ck(R.mul(factor,quotient),D,name+' full linear factor');result='empty: proper linear factor with constant nonzero b coefficient'
        else:
            assert all(not v for v in J);assert len(params)==10
            result='whole actual Poisson A10; parent identification below';full.append((rec,R,P))
        out.append(dict(source=name,full_K1_and_FULL9_affine_P2_reconstructed=True,all_homogeneous_P2_and_350_compatibilities_retained=True,all350_original_and_fullPf_verified=True,conclusion=result))
    return out,full

def sparsemat(d):
    m,n=d['shape'];out=[[F(0)]*n for _ in range(m)]
    for key,c in d['entries'].items():i,j=map(int,key.split(','));out[i][j]=F(c)
    return out
def mm(A,B):
    cols=list(zip(*B));return [[sum(a*b for a,b in zip(row,col)if a and b)for col in cols]for row in A]
def cert(M,C):
    rr,cc=C['minor_rows'],C['minor_columns'];rank=C['rank'];assert rank==len(rr)==len(cc)
    det=H.det([[M[i][j]for j in cc]for i in rr]);ck(det,F(C['minor_determinant']),'existing full minor')
    assert det;K=sparsemat(C['complete_kernel_sparse']);free=[i for i in range(len(M[0]))if i not in cc]
    ck(free,C['free_columns'],'free coordinates');assert len(free)==C['nullity']
    ck(mm(M,K),[[F(0)]*len(free)for _ in M],'complete kernel all original rows')
    ck([K[i]for i in free],[[F(i==j)for j in range(len(free))]for i in range(len(free))],'complete kernel independent identity')
    return K,str(det)
def constant_tensor(R,P):return H.skew({ij:{e[R.n:]:c for e,c in P[ij].items()}for ij in PAIRS})
def polynomial_tensor(R,P):return skew(R,{ij:{(0,)*R.n+e:c for e,c in P.get(ij,{}).items()}for ij in PAIRS})
def ambient_review(rec,slice_rec,source):
    R=LaurentRing((),COORDS);P=tensor(R,rec['complete_actual_integer_brackets']);Ph=constant_tensor(R,P)
    sr,SR,SP=source
    anchors=load('reports/CLASSIFICATION_PLANE_M2_HF_TWO_FULL_FAMILY_ANCHORS_20261108.json')['anchors']
    anchor=next(a for a in anchors if a['source']==rec['source'])
    expected=H.skew({ij:H.scale(evaluate(SR,SP[ij],anchor['chosen_fixed_unit_parameters']),F(rec['integer_scale']))for ij in PAIRS})
    ck(H.vector(Ph),H.vector(expected),'same whole-family anchor')
    assert not any(H.jac(Ph).values());J0=H.jac(Ph)
    Dcols=[]
    for ij,e in H.COLS:
        B=H.skew({ij:{e:F(1)}});J=H.jac(H.plus(Ph,B));JB=H.jac(B)
        Dcols.append(H.jvector({q:H.add(J[q],H.scale(J0[q],-1),H.scale(JB[q],-1))for q in QUADS}))
    D=H.transpose(Dcols);X=[H.a(1,i)for i in range(5)];Vcols=[];Gcols=[]
    for i in range(5):
        for j in range(5):
            Y=[{}for _ in range(5)];Y[i]=X[j];Vcols.append(H.vector(H.wedge(X,Y)))
            G=H.skew({(a,b):H.add(Ph.get((j,b),{})if a==i else{},Ph.get((a,j),{})if b==i else{},H.scale(H.mul(X[j],H.der(Ph.get((a,b),{}),i)),-1))for a,b in PAIRS})
            Gcols.append(H.vector(G))
    V=H.transpose(Vcols);M=H.transpose(Vcols+Gcols+[H.vector(Ph)])
    Ks={};dets={}
    for name,mat in [('D350x150',D),('Euler150x25',V),('Euler_GL_scalar150x51',M)]:
        ck(mat,sparsemat(rec['complete_sparse_matrices'][name]),name+' complete independent reconstruction')
        Ks[name],dets[name]=cert(mat,rec['certificates'][name])
    ck(mm(D,M),[[F(0)]*51 for _ in D],'all D times actual group zero')
    for k,v in enumerate(pfs(R,P)):ck(v,R.parse(rec['all_five_four_by_four_Pfaffians'][str(k)]),'full five Pfaffians')
    ck(actual(R,P),R.parse(rec['actual_signed_Euler_Pfaffian']),'signed Euler Pfaffian')
    w=rec['normalized_actual_point_torus']['coordinate_weights'];assert rec['normalized_actual_point_torus']['tensor_scale_exponent']==-1
    wc=[w[i]+w[j]-sum(a*b for a,b in zip(e,w))-1 for(i,j),e in H.COLS]
    wr=[sum(w[i]for i in q)-sum(a*b for a,b in zip(e,w))-2 for q,e in H.ROWS]
    ck(wc,rec['all150_input_weights'],'all150 weights');ck(wr,rec['all350_equation_weights'],'all350 weights')
    assert all(not v or wc[i]==0 for i,v in enumerate(H.vector(Ph)))
    assert all(not v or wr[i]==wc[j]for i,row in enumerate(D)for j,v in enumerate(row))
    gw=[]
    for col in zip(*M):
        chars={wc[i]for i,v in enumerate(col)if v};assert len(chars)<=1;gw.append(next(iter(chars))if chars else None)
    ck(gw,rec['all51_group_column_weights'],'all group characters')
    KD=Ks['D350x150'];free=rec['certificates']['D350x150']['free_columns']
    assert all(not v or wc[i]==wc[free[j]]for i,row in enumerate(KD)for j,v in enumerate(row))
    ng=Counter(wc[j]for j in free);mg=Counter(gw[j]for j in rec['certificates']['Euler_GL_scalar150x51']['minor_columns'])
    ck({str(k):ng[k]-mg[k]for k in set(ng)|set(mg)if ng[k]!=mg[k]},{'-2':1},'whole normal quotient')
    idx=slice_rec['normal_from_existing_kernel_basis_index'];N=[row[idx]for row in KD]
    ck(N,[F(slice_rec['full150_normal_sparse'].get(str(i),0))for i in range(150)],'actual full normal')
    assert all(not v or wc[i]==-2 for i,v in enumerate(N))
    ell=[F(slice_rec['annihilator_of_all51_actual_columns'].get(str(i),0))for i in range(150)]
    ck([sum(a*b for a,b in zip(ell,c))for c in zip(*M)],[F(0)]*51,'whole group annihilator')
    nonzero=sum(a*b for a,b in zip(ell,N));ck(nonzero,F(slice_rec['normal_non_group_exact_functional_value']),'normal non-group value');assert nonzero
    ck(slice_rec['selected101_equation_rows'],rec['certificates']['D350x150']['minor_rows'],'same implicit101 rows')
    U=rec['certificates']['D350x150']['minor_columns'];ck(U,slice_rec['selected101_coordinate_columns'],'same implicit101 columns')
    ck(F(slice_rec['existing101_block_determinant_reverified']),F(dets['D350x150']),'same101det')
    ck(F(slice_rec['existing48_group_minor_determinant_reverified']),F(dets['Euler_GL_scalar150x51']),'same48det')
    return dict(source=rec['source'],all150_350_D_Euler_and_actual51_reconstructed=True,ranks=[101,24,48],complete_kernels_and_all_nonzero_minor_determinants=True,whole_normal_character=-2,normal_non_group_value=str(nonzero),true_G48_U101_N1_direct_sum=True),Ph,N,U,wc

def graph_review(rec,Ph,N,U,wc):
    R=LaurentRing(('T',),COORDS);bound=rec['complete_graph_degree_bound'];assert bound in(3,4)
    coeffs={(i,0):v for i,v in enumerate(H.vector(Ph))if v}
    coeffs.update({(i,1):v for i,v in enumerate(N)if v})
    assert bound==-min(wc)//2
    for layer in rec['all_four_or_three_fixed_graph_layers']:
        degree=layer['degree'];assert degree in range(2,bound+1)
        entries={int(i):F(c)for i,c in layer['entire150_layer_sparse'].items()}
        assert set(entries)<=set(U)and all(wc[i]==-2*degree for i in entries)
        ck({str(U.index(i)):str(c)for i,c in entries.items()},layer['implicit101_coefficients_sparse'],'whole implicit layer')
        coeffs.update({(i,degree):v for i,v in entries.items()})
    ck({str(i)+','+str(j):str(v)for(i,j),v in coeffs.items()},rec['graph_coefficients_sparse'],'entire graph coefficient array')
    d={ij:{}for ij in PAIRS}
    for(i,k),v in coeffs.items():ij,e=H.COLS[i];d[ij][(k,)+e]=v
    P=skew(R,d);J=jac(R,P);assert not any(J.values()),'all350 full polynomial graph'
    ck(rec['complete_all350_equation_degree_bound'],2*bound,'complete finite degree bound')
    ck(rec['all350_full_equation_coefficients_sparse'],{},'entire equation array zero')
    assert all(k>=0 and k<=bound for(k,*_)in next(iter(d.values())))if next(iter(d.values()))else True
    at=F(rec['integer_parameter_rescale_t_equals_scale_times_e']);endpoint=H.skew({ij:evaluate(R,P[ij],{'T':at})for ij in PAIRS})
    ck(H.vector(endpoint),H.vector(constant_tensor(LaurentRing((),COORDS),tensor(LaurentRing((),COORDS),rec['all_integer_e1_graph_brackets']))),'full integer e1 endpoint')
    return dict(source=rec['source'],whole_IFT_graph_degree=bound,all350_full_equation_degrees_checked=2*bound,all_original_coefficients_zero=True,all_graph_layers_in_same_U_and_forced_weights=True,formal_completeness_by_unique_equivariant_IFT=True,integer_e1_parameter=str(at)),endpoint

class QuadraticFieldRing(LaurentRing):
    """QQ[theta]/(theta²-6) with polynomial coordinate coefficients."""
    def __init__(self,coords):super().__init__(('theta',),coords,invert=('theta',))
    def reduce(self,P):
        out={}
        for e,c in P.items():
            q,r=divmod(e[0],2);ee=(r,)+e[1:];out[ee]=out.get(ee,F(0))+c*F(6)**q
        return {e:c for e,c in out.items()if c}
    def add(self,*Ps):return self.reduce(super().add(*Ps))
    def mul(self,*Ps):return self.reduce(super().mul(*Ps))
    def parse(self,s):return self.reduce(super().parse(str(s).replace('sqrt(6)','theta')))

def powpoly(R,P,k):
    out=R.constant(1)
    for _ in range(k):out=R.mul(out,P)
    return out
def substitute(R,P,g):
    out={}
    for e,c in P.items():
        v={e[:R.n]+(0,)*5:c}
        for k,n in enumerate(e[R.n:]):v=R.mul(v,powpoly(R,g[k],n))
        out=R.add(out,v)
    return out
def matmul(R,A,B):return [[R.add(*(R.mul(a,b)for a,b in zip(row,col)))for col in zip(*B)]for row in A]
def detpoly(R,A):
    if len(A)==1:return A[0][0]
    return R.add(*(R.scale(R.mul(A[0][j],detpoly(R,[row[:j]+row[j+1:]for row in A[1:]])),(-1)**j)for j in range(len(A))if A[0][j]))
def push(R,P,M,Mi):
    X=[R.variable(s)for s in R.names[R.n:]]
    old=[R.add(*(R.mul(Mi[i][j],X[j])for j in range(5)))for i in range(5)]
    return skew(R,{(i,j):substitute(R,R.add(*(R.mul(M[i][a],M[j][b],P[a,b])for a in range(5)for b in range(5))),old)for i,j in PAIRS})
def wedge(R,X,Y):return skew(R,{(i,j):R.add(R.mul(X[i],Y[j]),R.scale(R.mul(X[j],Y[i]),-1))for i,j in PAIRS})
def projective_three(R,P):
    X=[R.variable(s)for s in R.names[R.n:]]
    return {q:R.add(R.mul(X[q[0]],P[q[1],q[2]]),R.scale(R.mul(X[q[1]],P[q[0],q[2]]),-1),R.mul(X[q[2]],P[q[0],q[1]]))for q in TRIPLES}
def identity_review(kind,rec,endpoint,anchor):
    R=QuadraticFieldRing(tuple('x'+str(i)for i in range(5)));p=R.parse
    M=[[p(c)for c in row]for row in rec['actual_PGL_matrix']];Mi=[[p(c)for c in row]for row in rec['actual_PGL_inverse']]
    I=[[R.constant(i==j)for j in range(5)]for i in range(5)]
    ck(matmul(R,M,Mi),I,'PGL inverse right');ck(matmul(R,Mi,M),I,'PGL inverse left')
    determinant=detpoly(R,M);ck(determinant,p(rec['actual_PGL_determinant']),'full PGL determinant');assert determinant
    P=polynomial_tensor(R,endpoint);pushed=push(R,P,M,Mi)
    rename=lambda s:str(s).replace('sqrt(6)','theta').replace('a','x0').replace('b','x1').replace('c','x2').replace('u','x3').replace('v','x4')
    # The source polynomial is already in the same exponent slots. No textual
    # coordinate substitution is used on field constants or completed targets.
    original=skew(R,{ast.literal_eval(k):p(rename(v))for k,v in rec['source_integer_e1_brackets'].items()})
    ck(P,original,'source actual e1 tensor')
    targetpacket=load('reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T'+str(kind)+'.json')
    params=rec['completed_T'+str(kind)+'_parameters'];S=LaurentRing(tuple(params),tuple('x'+str(i)for i in range(5)),invert=('p0',))
    # Evaluate the literal full completed family, including every period-
    # denominator correction; p0 is a fixed nonzero field constant here.
    target={}
    for key,value in targetpacket['homogeneous_bivector'].items():
        poly=S.parse(value);value={}
        for e,c in poly.items():
            term={(0,)+e[S.n:]:c}
            for par,k in zip(params,e[:S.n]):
                if k<0:
                    inv=R.inverse_monomial(p(params[par]));term=R.mul(term,powpoly(R,inv,-k))
                else:term=R.mul(term,powpoly(R,p(params[par]),k))
            value=R.add(value,term)
        target[tuple(map(int,key.split(',')))]=value
    T=skew(R,target)
    ck(T,tensor(R,rec['completed_T'+str(kind)+'_brackets']),'full completed MPS corrections')
    ck(pushed,tensor(R,rec['pushed_brackets']),'all150 PGL push')
    # Vector transform Y=M R(M^-1 X).
    X=[p('x'+str(i))for i in range(5)];old=[R.add(*(R.mul(Mi[i][j],X[j])for j in range(5)))for i in range(5)]
    Eold=[p(rename(v))for v in rec['Euler_field_in_original_coordinates']]
    EY=[R.add(*(R.mul(M[i][j],substitute(R,Eold[j],old))for j in range(5)))for i in range(5)]
    ck(EY,[p(s)for s in rec['Euler_field_in_T'+str(kind)+'_coordinates']],'actual Euler push')
    ck(pushed,plus(R,T,wedge(R,X,EY)),'full150 actual MPS plus literal Euler identity')
    diff=plus(R,pushed,skew(R,{ij:R.scale(T[ij],-1)for ij in PAIRS}))
    assert not any(projective_three(R,diff).values());assert not any(jac(R,T).values())
    Ds=actual(R,P);Dt=actual(R,T)
    ck(Ds,p(rename(rec['source_actual_Pfaffian_quintic'])),'source full signed Pf')
    ck(Dt,p(rec['target_actual_Pfaffian_quintic']),'target full signed Pf')
    ck(Dt,R.mul(determinant,substitute(R,Ds,old)),'signed determinant Pf transform')
    # Explicit inverse projective tensor identity, retaining its Euler field.
    back=push(R,T,Mi,M);XX=[R.add(*(R.mul(M[i][j],X[j])for j in range(5)))for i in range(5)]
    backE=[substitute(R,v,XX)for v in EY];backE=[R.add(*(R.mul(Mi[i][j],backE[j])for j in range(5)))for i in range(5)]
    ck(P,plus(R,back,wedge(R,X,backE)),'full150 inverse actual identity')
    if kind==38:
        # Intrinsic correction: at the original anchor the entire E wedge P,
        # not a raw bivector entry, vanishes on the cone.
        Z=LaurentRing(('s','t'),COORDS);pa=polynomial_tensor(Z,anchor)
        g=[Z.parse(v)for v in('-2*s*t**2','b','t**3','3*s**3','3*s**2*t')]
        assert not any(substitute(Z,v,g)for v in projective_three(Z,pa).values())
    return dict(completed_MPS='T'+str(kind),coefficient_field='QQ(sqrt6)'if kind==36 else'QQ',fixed_PGL_determinant=rec['actual_PGL_determinant'],both_PGL_matrix_inverses=True,full150_forward_inverse_actual_Euler_identities=True,all10_intrinsic_projective_identities=True,full_completed_source_all_corrections_reconstructed=True,signed_full_Pfaffian_determinant_identity=True,source_parameter_localizations='none; only fixed nonzero field constants',intrinsic_B_cone_rankzero_correction_checked=(kind==38))

def main():
    started=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(3072*1024**2,)*2)
    contract=load(JOB)
    for b in contract['exact_inputs_and_hashes']:ck(hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest(),b['sha256'],b['path'])
    assert not(ROOT/OUT).exists()
    support,full=source_review(load(SUP),load(BASE));print('whole source8 PASS',flush=True)
    ambient=load(AMB);sl=load(SLICE);results=[]
    ck(ambient['full_original_coefficient_order'],[dict(pair=list(ij),exponent=list(e))for ij,e in H.COLS],'full original150 order')
    ck(ambient['full_original_equation_order'],[dict(quad=list(q),exponent=list(e))for q,e in H.ROWS],'full original350 order')
    for i,(a,b,s)in enumerate(zip(ambient['cases'],sl['cases'],full)):
        ar,Ph,N,U,wc=ambient_review(a,b,s);gr,endpoint=graph_review(b,Ph,N,U,wc)
        kind=36 if i==0 else 38;ident=load('reports/CLASSIFICATION_PLANE_M2_HF_'+('A_T36'if i==0 else'B_T38')+'_ACTUAL_IDENTITY_20261108.json')
        ir=identity_review(kind,ident,endpoint,Ph);results.append(dict(ambient=ar,whole_graph=gr,actual_identification=ir));print('full HF',kind,'PASS',flush=True)
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,scope='All8 saved whole HF sources: six irreducible quintic empty, full A10 A/B contained in known T36/T38 by true complete finite equivariant slices and actual identities. Accepted earlier34/51 whole symbol completeness and G source upper23 are explicit geometric inputs; no new profile or global component count.',support8=support,two_actual_parent_proofs=results,no_producer_or_CAS_import=True,new_M1_or_point_or_profile_count=0,all_parameter_complements_preserved=True,reviewer_source_bytes=Path(__file__).read_text(),reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input_bindings=contract['exact_inputs_and_hashes'],seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/OUT).open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k]for k in('all_passed','seconds','maximum_RSS_KiB')}))

if __name__=='__main__':main()
