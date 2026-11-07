#!/usr/bin/env python3
"""Independent Fraction reconstruction of the whole E3 principal-orbit exclusion."""
import hashlib,json,resource,time,sys
from pathlib import Path
import classification_unimodular_logarithmic_review as r

def nilpotent_incidence():
    resource.setrlimit(resource.RLIMIT_CPU,(110,111));resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_BOUNDARY_20261107.json';data=json.loads(path.read_text())
    xx=[r.atom(1,i) for i in range(4)];x,y,z,w=xx
    def deriv(p):return [r.der(p,i) for i in range(4)]
    def contraction(v,form):return r.add(*(r.mul(a,b) for a,b in zip(v,form)))
    def lie(v,form):return [r.add(*(r.add(r.mul(v[i],r.der(form[j],i)),r.mul(form[i],r.der(v[i],j))) for i in range(4))) for j in range(4)]
    def wedge(a,b):return [r.add(r.mul(a[i],b[j]),r.scale(r.mul(a[j],b[i]),-1)) for i,j in r.PAIRS]
    def theta_of_Q(Q):
        gam=[Q[5],r.scale(Q[4],-1),Q[3],Q[2],r.scale(Q[1],-1),Q[0]]
        a={}
        for (i,j),v in zip(r.PAIRS,gam):a[i,j]=v;a[j,i]=r.scale(v,-1)
        return [r.scale(r.add(*(r.mul(xx[i],a.get((i,j),{})) for i in range(4))),r.F(1,4)) for j in range(4)]
    def weighted(rr,ss,P):
        dr,ds,dp=deriv(rr),deriv(ss),deriv(P);out=[];nil={}
        for l1,l2,H2 in [(rr,nil,nil),(ss,nil,nil),(nil,rr,nil),(nil,ss,nil),(nil,nil,r.mul(rr,rr)),(nil,nil,r.mul(rr,ss)),(nil,nil,r.mul(ss,ss))]:
            out.append([r.add(r.mul(P,r.add(r.mul(l1,a),r.mul(l2,b))),r.scale(r.mul(r.add(r.mul(rr,l1),r.mul(ss,l2)),c),r.F(-1,2)),r.mul(H2,r.add(r.mul(ss,a),r.scale(r.mul(rr,b),-1)))) for a,b,c in zip(dr,ds,dp)])
        return out
    def signature(form):return tuple(tuple(sorted(p.items())) for p in form)
    def monomials(degree,indices):
        out=[]
        for inds in r.combinations_with_replacement(indices,degree):
            p=r.atom()
            for i in inds:p=r.mul(p,xx[i])
            out.append(p)
        return out
    records=[];expected_dims={'31':7,'22':7,'4':3,'211':15}
    assert {c['jordan_type'] for c in data['cases']}==set(expected_dims)
    for c in data['cases']:
        begin=time.monotonic();name=c['jordan_type']
        Y={'31':[{},x,y,{}],'22':[{}, {},x,y],'4':[{},x,y,z],'211':[{}, {},{},x]}[name]
        assert Y==list(map(r.parse,c['nilpotent_vector_field']))
        m3=list(map(r.parse,c['cubic_monomial_basis']));m4=list(map(r.parse,c['quartic_monomial_basis']))
        assert {signature([p]) for p in m3}=={signature([p]) for p in monomials(3,range(4))} and len(m3)==20
        assert {signature([p]) for p in m4}=={signature([p]) for p in monomials(4,range(4))} and len(m4)==35
        e3=[next(iter(p)) for p in m3];e4=[next(iter(p)) for p in m4]
        columns=[]
        for i in range(4):
            for mon in m3:
                theta=[mon if j==i else {} for j in range(4)]
                features=[contraction(xx,theta),contraction(Y,theta)]
                columns.append([p.get(e,r.RF()) for p in features for e in e4]+[p.get(e,r.RF()) for p in lie(Y,theta) for e in e3])
        M=r.transpose(columns);assert M==r.matrix(c['constraint_matrix'])
        if name=='31':expected=weighted(x,w,r.add(r.mul(y,y),r.scale(r.mul(x,z),-2)))
        elif name=='22':expected=weighted(x,y,r.add(r.mul(x,w),r.scale(r.mul(y,z),-1)))
        elif name=='4':expected=[theta_of_Q(wedge(Y,V)) for V in [[r.scale(x,-5),r.scale(y,-1),r.scale(z,3),r.scale(w,7)],[{},{},x,y],[{},{},{},x]]]
        else:
            expected=[]
            for i,j in [(0,1),(0,2),(1,2)]:
                for mon in monomials(2,range(3)):
                    form=[{} for _ in range(4)];form[j]=r.mul(xx[i],mon);form[i]=r.scale(r.mul(xx[j],mon),-1);expected.append(form)
        forms=[list(map(r.parse,f)) for f in c['expected_basic_forms']]
        assert {signature(f) for f in forms}=={signature(f) for f in expected}
        K=r.transpose([[p.get(e,r.RF()) for p in f for e in e3] for f in forms])
        assert K==r.matrix(c['expected_basis_matrix']) and r.annihilate(M,K)
        dim=expected_dims[name];assert c['kernel_dimension']==dim
        md=r.det(r.submatrix(M,c['pivot_rows'],c['pivot_columns']))
        kd=r.det(r.submatrix(K,c['kernel_pivot_rows'],c['kernel_pivot_columns']))
        assert md and md==r.scalar(c['pivot_determinant']) and len(c['pivot_rows'])==80-dim
        assert kd and kd==r.scalar(c['kernel_pivot_determinant']) and len(c['kernel_pivot_rows'])==dim
        records.append({'jordan_type':name,'all80_coefficients_retained':True,'constraint_matrix_reconstructed':True,
            'complete_kernel_dimension':dim,'mathematical_expected_basis_independently_reconstructed':True,
            'two_nonzero_complementary_rank_minors':True,'all_passed':True,'seconds':time.monotonic()-begin})
        print(name,'PASS',dim,round(time.monotonic()-begin,3),flush=True)
    # Exact polynomial identity in formal b=n0,k=n1, r=x,s=y.
    b,k=r.atom(1,4),r.atom(1,5);ss=y;rr=x;omega=[ss,r.scale(rr,-1)]
    alpha=[r.mul(b,ss),r.mul(r.add(r.scale(k,-2),r.scale(b,-1)),rr)]
    for i,Fm in enumerate([r.mul(ss,ss),r.mul(rr,ss),r.mul(rr,rr)]):
        value=[r.add(r.mul(Fm,a),r.mul(r.mul(k,r.mul(rr,ss)),r.der(Fm,j))) for j,a in enumerate(alpha)]
        assert value==[r.mul(r.mul(r.add(b,r.scale(k,i)),Fm),v) for v in omega]
    assert data['forcing_operator_eigenvalues']==['b','b + k','b + 2*k']
    out={'task_id':'G','all_passed':True,'cases':records,'forcing_operator_identity_independently_reconstructed':True,
        'method':'Standard-library Fraction sparse-polynomial reconstruction of all80 basis columns, full constraint matrices and explicit mathematical kernel bases; complementary nonzero minors prove completeness. No producer/CAS imported.',
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,Path(r.__file__)]},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'Arithmetic completeness for all four nonzero nilpotent Jordan types and weighted forcing identity. Projective incidence coverage, dense L112 inclusion and regular centralizer conjugation require the separate geometric proof.'}
    (root/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_BOUNDARY_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'cases':4,'seconds':out['seconds']}))

def main():
    if '--nilpotent-incidence' in sys.argv:return nilpotent_incidence()
    resource.setrlimit(resource.RLIMIT_CPU,(25,26));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_20261107.json';c=json.loads(path.read_text())
    X=[r.scale(r.atom(1,i),r.F(a,4)) for i,a in enumerate([-5,-1,3,7])]
    Y=[{},r.atom(4,0),r.atom(4,1),r.atom(4,2)]
    Q=[r.add(r.mul(Y[i],X[j]),r.scale(r.mul(X[i],Y[j]),-1)) for i,j in r.PAIRS]
    assert Q==list(map(r.parse,c['highest_bivector'])) and not r.pf(Q) and not any(r.divergence(Q))
    cols=[];mon3=set()
    for i in range(6):
        for j in range(4):
            L=[r.atom(1,j) if k==i else {} for k in range(6)];cross=r.cross(Q,L);mon3.update(cross)
            cols.append(([v.get(r.ZERO,r.RF()) for v in r.divergence(L)],cross))
    M=r.transpose([d+[p.get(e,r.RF()) for e in sorted(mon3,reverse=True)] for d,p in cols])
    assert M==r.matrix(c['full_linear_matrix'])
    md=r.det(r.submatrix(M,c['full_linear_pivot_rows'],c['full_linear_pivot_columns']))
    assert len(c['full_linear_pivot_rows'])==20 and md and md==r.scalar(c['full_linear_pivot_determinant'])
    T=r.transpose([r.linearvec([r.der(p,j) for p in Q]) for j in range(4)])
    assert T==r.matrix(c['translation_basis']) and r.annihilate(M,T)
    td=r.det(r.submatrix(T,c['translation_cut_rows'],c['translation_pivot_columns']))
    assert len(c['translation_cut_rows'])==4 and td and td==r.scalar(c['translation_pivot_determinant'])
    assert c['normal_dimension']==0 and not any(map(r.parse,c['normal_linear_bivector']))
    mon2=[]
    for i,j in r.combinations_with_replacement(range(4),2):
        e=[0]*r.N;e[i]+=1;e[j]+=1;mon2.append(tuple(e))
    G=r.transpose([[r.cross(Q,[r.atom() if k==i else {} for k in range(6)]).get(e,r.RF()) for e in mon2] for i in range(6)])
    assert G==r.matrix(c['constant_equation_matrix'])
    gd=r.det(r.submatrix(G,c['constant_pivot_rows'],c['constant_pivot_columns']))
    assert len(c['constant_pivot_rows'])==6 and gd and gd==r.scalar(c['constant_pivot_determinant'])
    assert not any(map(r.parse,c['constant_source'])) and not r.parse(c['constant_Pfaffian']) and c['all_remaining_positive_degree_Pfaffian_coefficients']==[]
    helper=Path(r.__file__)
    out={'task_id':'G','all_passed':True,'method':'Independent standard-library Fraction reconstruction from the two linear vector fields in Pym2014 Section3.6; no producer or CAS arithmetic imported.',
      'highest_Q_is_minus_X_wedge_Y':True,'highest_Pfaffian_and_divergence_zero':True,'full24_linear_matrix_rank':20,'full_kernel_equals_translation_image_dimension':4,'true_translation_normal_dimension':0,'full6_constant_matrix_rank':6,
      'linear_rank_minor':c['full_linear_pivot_determinant'],'translation_rank_minor':c['translation_pivot_determinant'],'constant_rank_minor':c['constant_pivot_determinant'],
      'conclusion':'Every affine quadratic unimodular extension over the principal E3 highest orbit with constant Pfaffian has Pfaffian zero. The prescribed nonzero-constant-Pfaffian support is empty.',
      'boundary_scope':'No orbit-closure extension is inferred; closed L112, pullback and common-factor boundary obligations remain separate.',
      'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,helper]},'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    (root/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds'],'ranks':[20,4,6]}))
if __name__=='__main__':main()
