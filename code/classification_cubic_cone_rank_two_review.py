#!/usr/bin/env python3
"""Independent full150 cone rank-two G replay; no producer imports or GB.

All geometric polynomial arithmetic is independent sparse Fraction. SymPy is
used only as an exact one-variable rational-function coefficient field for
the three whole lambda modules, with supplied full kernels plus nonzero
minors proving both rank bounds. No sampled lambda values replace a family.
"""
import hashlib,itertools as it,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
import sympy as s
import classification_cubic_scroll_abelian_review as a
import classification_quadric_cone_zero_review as z

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/cubic_cone_rank_two_review.json'
OUT='reports/CLASSIFICATION_CUBIC_CONE_RANK_TWO_REVIEW_20261107.json'
SRC='reports/CLASSIFICATION_CUBIC_SURFACE_TANGENT_TWISTED_CUBIC_CONE_03_20261107.json'
BV='reports/CLASSIFICATION_CUBIC_CONE_RANK_TWO_BASE_VERTEX_20261107.json'
SUP='reports/CLASSIFICATION_CUBIC_CONE_RANK_TWO_BASE_SUPPORT_20261107.json'
MOD='reports/CLASSIFICATION_CUBIC_CONE_RANK_TWO_FIVE_MODULES_20261107.json'
P=a.PAIRS;M2=a.monomials(5,2);COLS=list(it.product(P,M2));NAMES=['x'+str(i)for i in range(5)]

def rows(v):return [{int(j):Q(c)for j,c in r.items()if Q(c)}for r in v]
def mat(v):return [[Q(c)for c in r]for r in v]
def dot(r,v):return sum(c*v[j]for j,c in r.items())
def mm(A,B):return [[sum(x*y for x,y in zip(r,c))for c in zip(*B)]for r in A]
def ident(n):return [[Q(i==j)for j in range(n)]for i in range(n)]
def tensor(v):
    R=[[{}for j in range(5)]for i in range(5)]
    for ((i,j),e),c in zip(COLS,v):
        if c:R[i][j][e]=c;R[j][i][e]=-c
    return R
def vector(R):return [R[i][j].get(e,Q())for(i,j),e in COLS]
def bracket(R,f,g):return a.add(*(a.mul(a.mul(a.diff(f,i),R[i][j]),a.diff(g,j))for i,j in it.product(range(5),repeat=2)))

def complete_source(packet):
    saved=packet['linear'];xs=[z.var(5,i)for i in range(5)]
    qs=[a.parse(t,NAMES)for t in ['x0*x2-x1**2','x0*x3-x1*x2','x1*x3-x2**2']]
    powers=[(3,0,0),(2,1,0),(1,2,0),(0,3,0),(0,0,1)]
    def pull(poly):
        out={}
        for e,c in poly.items():
            f=tuple(sum(e[i]*powers[i][j]for i in range(5))for j in range(3));out[f]=out.get(f,Q())+c
        return {e:c for e,c in out.items()if c}
    assert all(not pull(q)for q in qs)
    target={};conormal={}
    def ins(dest,label,col,poly):
        for e,c in pull(poly).items():dest.setdefault(label+(e,),{})[col]=c
    for col,((i,j),e)in enumerate(COLS):
        R=tensor([Q(k==col)for k in range(150)])
        for g,q in enumerate(qs):
            V=[a.add(*(a.mul(R[k][l],a.diff(q,l))for l in range(5)))for k in range(5)]
            for u,v in P:ins(target,(g,u,v),col,a.add(a.mul(xs[u],V[v]),a.scale(a.mul(xs[v],V[u]),-1)))
        for g,h in it.combinations(range(3),2):
            H=a.add(a.mul(xs[0],bracket(R,qs[g],qs[h])),a.scale(a.mul(qs[h],bracket(R,qs[g],xs[0])),-2),a.scale(a.mul(qs[g],bracket(R,qs[h],xs[0])),2))
            for k in range(5):ins(conormal,(g,h,k),col,a.diff(H,k))
    tl=sorted(target);cl=sorted(conormal);T=[target[k]for k in tl];C=[conormal[k]for k in cl];L=T+C
    assert [[*v[:-1],list(v[-1])]for v in tl]==saved['full_tangent_row_labels']
    assert [[*v[:-1],list(v[-1])]for v in cl]==saved['full_intrinsic_conormal_firstpartial_labels']
    assert T==rows(saved['full_tangent_rows'])and C==rows(saved['full_intrinsic_conormal_rows'])
    assert len(L)==975 and len(a.rref(T))==87
    K=mat(saved['full150_kernel_basis']);KT=list(map(list,zip(*K)))
    assert len(KT)==51 and len(a.rref(KT))==51
    assert all(not dot(r,v)for r in L for v in KT)
    cert=saved['rank_certificate'];minor=a.determinant([[L[i].get(j,Q())for j in cert['columns']]for i in cert['rows']])
    assert len(cert['rows'])==99 and minor==Q(cert['determinant'])==-192
    E=[]
    for i,j in it.product(range(5),repeat=2):
        Y=[{}for k in range(5)];Y[i]=xs[j];E.append(vector(z.wedge(xs,Y)))
    assert len(a.rref(E))==24 and all(not dot(r,v)for r in L for v in E)
    B=mat(saved['full150_by_source_basis']);BT=list(map(list,zip(*B)));family=[tensor(v)for v in BT]
    assert len(BT)==27 and len(a.rref(E+BT))==51 and all(not dot(r,v)for r in L for v in BT)
    for R,st in zip(family,saved['complete_source_basis']):
        assert all(R[i][j]==a.parse(st.get(str((i,j)),'0'),NAMES)for i,j in P)
        assert all(not a.add(*(a.diff(R[i][j],i)for i in range(5)))for j in range(5))
    return B,family,dict(full975x150_reconstructed=True,rank99_minor=str(minor),all51_kernel_products_zero=True,complete51_kernel_rank=True,all_Euler24=True,canonical27_complete=True)

def source_structure(B,family,packet):
    base=[i for i,R in enumerate(family)if any(R[u][v]for u,v in P if v<4)]
    vertex=[i for i in range(27)if i not in base]
    assert base==packet['base_source_indices']==[0,1,2,3,4,5,6,7,16]and vertex==packet['vertex_source_indices']
    for j,R in enumerate(family):
        if j in base:
            assert all(not e[4]for u,v in P if v<4 for e in R[u][v])
            for i in range(4):assert R[i][4]==a.shift(a.add(*(a.diff(R[k][i],k)for k in range(4))),4)
        else:assert all(not e[4]for u,v in P for e in R[u][v])
    J,D=a.jacobi_and_divisor(family);pairs=list(it.combinations_with_replacement(range(27),2));loc={v:i for i,v in enumerate(pairs)}
    pure_pairs=list(it.combinations_with_replacement(range(9),2));mixed_pairs=list(it.product(range(9),range(18)))
    purecols=[loc[tuple(sorted((base[i],base[j])))]for i,j in pure_pairs]
    mixedcols=[loc[tuple(sorted((base[i],vertex[j])))]for i,j in mixed_pairs]
    vv=[loc[tuple(sorted((i,j)))]for i,j in it.combinations_with_replacement(vertex,2)]
    assert all(not r.get(c)for r in J for c in vv)
    take=lambda A,cs:[{j:r[c]for j,c in enumerate(cs)if r.get(c)}for r in A]
    BJ,MJ=take(J,purecols),take(J,mixedcols);BD,MD=take(D,purecols),take(D,mixedcols)
    assert BJ==rows(packet['base_all350_rows'])and MJ==rows(packet['mixed_all350_rows'])
    assert BD==rows(packet['actualD_pure_base_rows'])and MD==rows(packet['actualD_mixed_rows'])
    assert all(not(r and t)for r,t in zip(BJ,MJ))
    RR=a.rref(BJ);assert len(RR)==12
    assert RR=={int(p):{int(j):Q(c)for j,c in r.items()}for p,r in packet['base_exact_row_basis'].items()}
    for p,phi in packet['base_exact_original350_functionals'].items():assert a.add(*(a.scale(BJ[int(i)],Q(c))for i,c in phi.items()))==RR[int(p)]
    cdep=[];all_cdep=True
    for row,e in zip(BD,a.monomials(5,5)):
        if e[4]:
            rem=row.copy()
            for p in sorted(RR):
                if p in rem:rem=a.add(rem,a.scale(RR[p],-rem[p]))
            all_cdep &= not rem
            if row:cdep.append(dict(exponent=list(e),residual={str(j):str(v)for j,v in rem.items()}))
    assert cdep==packet['cdependent_actualD_row_remainders']and all_cdep==packet['all_actualD_cdependent_rows_in_base_Jacobi_quadric_rowspace']
    return base,vertex,BJ,MJ,BD,MD,[RR[int(p)]for p in packet['base_exact_row_basis']],dict(full350_and126_reconstructed=True,whole_base9_vertex18_identity=True,all_vertex_vertex_Jacobi_zero=True,pure_and_mixed_spatial_rows_disjoint=True,base_equation_rank=12,old_cdependent_D_nonmembership_preserved=not all_cdep)

def representation(B,family,base,eq,packet):
    B9=[[r[i]for i in base]for r in B];F=[family[i]for i in base]
    for name,entries in [('h',{(0,0):3,(1,1):1,(2,2):-1,(3,3):-3}),('f',{(0,1):3,(1,2):2,(2,3):1}),('e',{(1,0):1,(2,1):2,(3,2):3})]:
        A=[[Q(entries.get((i,j),0))for j in range(5)]for i in range(5)]
        fields=[a.add(*(a.scale(z.var(5,k),c)for k,c in enumerate(r)))for r in A];acts=[]
        for R in F:
            T=[[{}for j in range(5)]for i in range(5)]
            for i,j in P:T[i][j]=a.add(*(a.mul(fields[k],a.diff(R[i][j],k))for k in range(5)),*(a.scale(R[k][j],-A[i][k])for k in range(5)),*(a.scale(R[i][k],-A[j][k])for k in range(5)))
            acts.append(vector(T))
        M=mat(packet['three_Sym2_action_matrices'][name]);assert mm(B9,M)==list(map(list,zip(*acts)))
    T=mat(packet['complete_coordinate_change']);Ti=mat(packet['complete_inverse']);assert mm(T,Ti)==ident(9)
    texts=['3*t0-t5','3*t4-2*t6','3*t2+3*t3-t7','2*t0-t5-t8','-4*(t1-2*t4+t6)/3','(3*t2+4*t3-2*t7)/3','t0','t1','t2+t3']
    names=['t'+str(i)for i in range(9)]
    assert T==[[a.parse(text,names).get(tuple(int(k==j)for k in range(9)),Q())for j in range(9)]for text in texts]
    std={'h':[[2,0,0],[0,0,0],[0,0,-2]],'f':[[0,0,0],[2,0,0],[0,1,0]],'e':[[0,1,0],[0,0,2],[0,0,0]]}
    for name,M in packet['three_Sym2_action_matrices'].items():assert mm(mm(T,mat(M)),Ti)==[[Q(std[name][i%3][j%3])if i//3==j//3 else Q()for j in range(9)]for i in range(9)]
    Z=[z.var(9,i)for i in range(9)];vals=[a.add(*(a.scale(Z[j],c)for j,c in enumerate(r)))for r in Ti]
    pairs=list(it.combinations_with_replacement(range(9),2))
    pol=[a.add(*(a.scale(a.mul(vals[pairs[j][0]],vals[pairs[j][1]]),c)for j,c in r.items()))for r in eq]
    zn=['u0','u1','u2','v0','v1','v2','w0','w1','w2'];assert pol==[a.parse(t,zn)for t in packet['transformed_all12_quadrics']]
    products=[a.mul(Z[0],Z[3]),a.add(a.mul(Z[0],Z[4]),a.mul(Z[1],Z[3])),a.add(a.mul(Z[0],Z[5]),a.mul(Z[1],Z[4]),a.mul(Z[2],Z[3])),a.add(a.mul(Z[1],Z[5]),a.mul(Z[2],Z[4])),a.mul(Z[2],Z[5])]
    for p,phi in zip(products,packet['all5_cartan_product_original12_functionals']):assert a.add(*(a.scale(pol[int(i)],Q(c))for i,c in phi.items()))==p
    for branch in packet['branches']:
        vanish=range(3)if branch['name']=='Uzero'else range(3,6);keep=[i for i in range(9)if i not in vanish]
        special=[{tuple(e[i]for i in keep):c for e,c in f.items()if not any(e[i]for i in vanish)}for f in pol]
        assert special==[a.parse(t,branch['parameters'])for t in branch['all12_specialized_equations']]
        simple=[a.parse(t,branch['parameters'])for t in branch['equivalent_factored_row_basis']]
        keys=sorted(set().union(*(f.keys()for f in special+simple)));to_rows=lambda L:[{j:f[e]for j,e in enumerate(keys)if e in f}for f in L]
        assert a.rref(to_rows(special))==a.rref(to_rows(simple))
        minors=[a.add(a.mul(z.var(6,i),z.var(6,j+3)),a.scale(a.mul(z.var(6,j),z.var(6,i+3)),-1))for i,j in it.combinations(range(3),2)]
        if branch['name']=='Uzero':minors.append(a.add(a.mul(z.var(6,1),z.var(6,1)),a.scale(a.mul(z.var(6,0),z.var(6,2)),-4)))
        keys=sorted(set().union(*(f.keys()for f in special+minors)));assert a.rref(to_rows(special))==a.rref(to_rows(minors))
    return Ti,dict(actual_SL2_full150_verified=True,all_three_Sym2_blocks=True,all_five_Cartan_original12_functionals=True,both_whole_branches_rowspace_equivalent=True)

def fielddet(M,K):
    A=[r[:]for r in M];ans=K.one
    for i in range(len(A)):
        p=next((p for p in range(i,len(A))if A[p][i]),None)
        if p is None:return K.zero
        if p!=i:A[p],A[i]=A[i],A[p];ans=-ans
        q=A[i][i];ans*=q
        for j in range(i+1,len(A)):
            if A[j][i]:
                r=A[j][i]/q
                for k in range(i+1,len(A)):A[j][k]-=r*A[i][k]
                A[j][i]=K.zero
    return ans

def five_modules(Ti,MJ,BD,MD,packet):
    lam=s.Symbol('lam');K=s.QQ.frac_field(lam);parse=lambda x:K.from_sympy(s.sympify(str(x),locals={'lam':lam}))
    fam=[('A_ss',[0,1,0,0,0,0,0,lam,0],4),('A_nil',[1,0,0,0,0,0,lam,0,0],3),('B_nil',[0,0,0,1,0,0,lam,0,0],3),('W_ss',[0,0,0,0,0,0,0,1,0],3),('W_nil',[0,0,0,0,0,0,1,0,0],2)]
    pair9=list(it.combinations_with_replacement(range(9),2));results=[];walltotal=0
    def check(M,rec,K0,read):
        B=[[read(v)for v in r]for r in rec['complete_kernel']];rank=rec['rank'];free=rec['kernel_free_rows']
        assert len(B)==18 and len(free)==18-rank
        assert [[B[i][j]for j in range(len(free))]for i in free]==[[K0.one if i==j else K0.zero for j in range(len(free))]for i in range(len(free))]
        assert all(not sum((r[j]*B[j][k]for j in range(18)),K0.zero)for r in M for k in range(len(free)))
        minor=fielddet([[M[i][j]for j in rec['pivot_columns']]for i in rec['selected_original350_rows']],K0)
        assert minor and minor==read(rec['maximal_minor'])
        upper=16+rec['base_stratum_dimension']+len(free)-1
        assert upper==rec['whole_projective_image_upper']and upper<=24 and rec['whole_G_closed_by_dimension']
        return B,minor,upper
    for (name,uvw,dim),rec in zip(fam,packet['records']):
        assert name==rec['case']and dim==rec['base_stratum_dimension']
        t=[sum((parse(c)*parse(v)for c,v in zip(r,uvw)),K.zero)for r in Ti]
        assert t==[parse(v)for v in rec['canonical_original_t']]
        M=[[K.zero for j in range(18)]for i in range(350)]
        for i,r in enumerate(MJ):
            for col,c in r.items():p,j=divmod(col,18);M[i][j]+=parse(c)*t[p]
        assert M==[[parse(v)for v in r]for r in rec['all350_module_matrix']]
        D0=[sum((parse(c)*t[pair9[j][0]]*t[pair9[j][1]]for j,c in r.items()),K.zero)for r in BD]
        DV=[[K.zero for j in range(18)]for i in range(126)]
        for i,r in enumerate(MD):
            for col,c in r.items():p,j=divmod(col,18);DV[i][j]+=parse(c)*t[p]
        assert D0==[parse(v)for v in rec['actualD_particular_coefficients']]and DV==[[parse(v)for v in r]for r in rec['actualD_vertex_matrix']]
        B,minor,upper=check(M,rec,K,parse);walls=[]
        if name.startswith(('A_','B_')):
            denom=s.Integer(1)
            for r in B:
                for v in r:denom=s.lcm(denom,s.fraction(K.to_sympy(v))[1])
            exception=s.cancel(K.to_sympy(minor)*denom)
            num,den=s.fraction(exception);own={str(s.Poly(f,lam,domain=s.QQ).monic().as_expr())for f,p in s.factor_list(num*den,lam)[1]}
            saved={str(s.Poly(s.sympify(w['wall']['irreducible_polynomial'],locals={'lam':lam}),lam).monic().as_expr())for w in rec['walls']}
            assert own==saved
            for wall in rec['walls']:
                f=s.Poly(s.sympify(wall['wall']['irreducible_polynomial'],locals={'lam':lam}),lam);assert f.degree()==1
                root=-f.nth(0)/f.nth(1);rr=Q(root)
                evaluate=lambda v:Q(K.to_sympy(v).subs(lam,root))
                WM=[[s.QQ(evaluate(v))for v in row]for row in M]
                wb,wm,wu=check(WM,wall,s.QQ,lambda v:s.QQ(Q(v)))
                assert wall['base_stratum_dimension']==dim-1
                walls.append(dict(root=str(rr),rank=wall['rank'],kernel=18-wall['rank'],image_upper=wu,minor=str(wm)));walltotal+=1
        else:assert not rec['walls']
        results.append(dict(case=name,generic_rank=rec['rank'],generic_kernel=18-rec['rank'],generic_image_upper=upper,all350x18_and_actual126_reconstructed=True,all_generic_kernel_products_and_free_block=True,generic_minor=str(K.to_sympy(minor)),all_exceptional_factors_and_walls=walls))
        print(name,'whole module and',len(walls),'walls PASS',flush=True)
    assert walltotal==19
    return results

def main():
    begun=time.monotonic();signal.alarm(300);resource.setrlimit(resource.RLIMIT_CPU,(290,295));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for r in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256']
    packets=[json.loads((ROOT/p).read_text())for p in [SRC,BV,SUP,MOD]]
    for p in packets:assert hashlib.sha256(p['producer_bytes'].encode()).hexdigest()==p['producer_sha256']
    B,fam,full=complete_source(packets[0]);print('full975/150 and 51/27 complete sources PASS',flush=True)
    base,vertex,BJ,MJ,BD,MD,eq,structure=source_structure(B,fam,packets[1]);print('all350/126 and base9/vertex18 PASS',flush=True)
    Ti,rep=representation(B,fam,base,eq,packets[2]);print('actualSL2/all150 and whole base support PASS',flush=True)
    modules=five_modules(Ti,MJ,BD,MD,packets[3])
    result=dict(task_id='G',all_passed=True,whole_cubic_cone_rank_two_G24=True,producer_imports=False,Groebner_or_jet_run=False,exact_arithmetic='Independent sparse Fraction; QQ(lam) coefficient-field operations only for whole univariate modules.',complete_source=full,base_vertex_structure=structure,whole_representation_and_support=rep,all_five_modules=modules,total_exact_linear_walls=19,whole_geometry_review='Five actual SL2/tensor-scaled orbit families cover every nonzero base solution, including U=V=0 and all nilpotent/zero complements; basezero yields rank at most two.',all_nonreduced_R_and_vertex_retained=True,normalization_flatness_assumed=False,E_membership_or_global_classification_claimed=False,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_bytes=Path(__file__).read_text(),seconds=time.monotonic()-begun,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=result['seconds'],maximum_RSS_KiB=result['maximum_RSS_KiB'])),flush=True)

if __name__=='__main__':main()
