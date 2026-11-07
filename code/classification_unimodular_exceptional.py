#!/usr/bin/env python3
"""All thirty low-degree coefficients over the unique E3 highest-part orbit."""
import hashlib,itertools,json,resource,signal,time,sys
from pathlib import Path
import sympy as s

def nilpotent_incidence():
    """Full horizontal/basic cubic-form spaces for every nonzero nilpotent type."""
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    xx=s.symbols('x y z w');x,y,z,w=xx
    def mons(degree,variables=xx):
        return sorted({s.prod(a) for a in itertools.combinations_with_replacement(variables,degree)},key=str,reverse=True)
    mon3,mon4=mons(3),mons(4)
    theta_basis=[s.Matrix([mon if k==i else 0 for k in range(4)]) for i in range(4) for mon in mon3]
    def features(theta,Y):
        pol=[sum(a*b for a,b in zip(xx,theta)),sum(a*b for a,b in zip(Y,theta))]
        out=[s.Poly(p,*xx).coeff_monomial(mon) for p in pol for mon in mon4]
        lie=[sum(Y[i]*s.diff(theta[j],xx[i])+theta[i]*s.diff(Y[i],xx[j]) for i in range(4)) for j in range(4)]
        return out+[s.Poly(p,*xx).coeff_monomial(mon) for p in lie for mon in mon3]
    def wedge(a,b):return [s.expand(a[i]*b[j]-a[j]*b[i]) for i,j in itertools.combinations(range(4),2)]
    def theta_of_Q(Q):
        gam=[Q[5],-Q[4],Q[3],Q[2],-Q[1],Q[0]];g=s.zeros(4)
        for (i,j),v in zip(itertools.combinations(range(4),2),gam):g[i,j]=v;g[j,i]=-v
        return s.Matrix([s.expand(sum(xx[i]*g[i,j] for i in range(4))/4) for j in range(4)])
    def weighted(r,t,P):
        dr=s.Matrix([s.diff(r,a) for a in xx]);dt=s.Matrix([s.diff(t,a) for a in xx]);dp=s.Matrix([s.diff(P,a) for a in xx]);out=[]
        for l1,l2,H2 in [(r,0,0),(t,0,0),(0,r,0),(0,t,0),(0,0,r*r),(0,0,r*t),(0,0,t*t)]:
            out.append((P*(l1*dr+l2*dt)-(r*l1+t*l2)*dp/2+H2*(t*dr-r*dt)).applyfunc(s.expand))
        return out
    records=[]
    for name,Y in [('31',(0,x,y,0)),('22',(0,0,x,y)),('4',(0,x,y,z)),('211',(0,0,0,x))]:
        Y=s.Matrix(Y);M=s.Matrix.hstack(*(s.Matrix(features(theta,Y)) for theta in theta_basis))
        _,pc=M.rref();_,pr=M.T.rref()
        if name=='31':forms=weighted(x,w,y*y-2*x*z)
        elif name=='22':forms=weighted(x,y,x*w-y*z)
        elif name=='4':
            forms=[theta_of_Q(wedge(Y,V)) for V in [s.Matrix([-5*x,-y,3*z,7*w]),s.Matrix([0,0,x,y]),s.Matrix([0,0,0,x])]]
        else:
            forms=[]
            for i,j in itertools.combinations(range(3),2):
                for mon in mons(2,xx[:3]):
                    f=s.zeros(4,1);f[j]=xx[i]*mon;f[i]=-xx[j]*mon;forms.append(f)
        K=s.Matrix.hstack(*(s.Matrix([s.Poly(v,*xx).coeff_monomial(mon) for v in f for mon in mon3]) for f in forms))
        assert M*K==s.zeros(M.rows,K.cols)
        _,kc=K.rref();_,kr=K.T.rref()
        assert len(pc)+len(kc)==80
        records.append({'jordan_type':name,'nilpotent_vector_field':list(map(str,Y)),
            'cubic_monomial_basis':list(map(str,mon3)),'quartic_monomial_basis':list(map(str,mon4)),
            'constraint_matrix':M.tolist(),'constraint_rank':len(pc),'pivot_rows':list(pr),'pivot_columns':list(pc),
            'pivot_determinant':str(M.extract(pr,pc).det()),'expected_basic_forms':[list(map(str,f)) for f in forms],
            'expected_basis_matrix':K.tolist(),'kernel_dimension':len(kc),'kernel_pivot_rows':list(kr),'kernel_pivot_columns':list(kc),
            'kernel_pivot_determinant':str(K.extract(kr,kc).det()),'all_basis_forms_are_horizontal_tangent_and_invariant':True})
    b,k=s.symbols('b k')
    # Three monomials r^i s^(2-i), i=0,1,2; forcing operator is diagonal.
    rr,ss=s.symbols('r s');Fmons=[ss*ss,rr*ss,rr*rr];omega=s.Matrix([ss,-rr]);alpha=s.Matrix([b*ss,(-2*k-b)*rr])
    forcing=[]
    for Fm in Fmons:
        v=(Fm*alpha+k*rr*ss*s.Matrix([s.diff(Fm,rr),s.diff(Fm,ss)])).applyfunc(s.expand)
        i=len(forcing);assert (v-(b+i*k)*Fm*omega).applyfunc(s.expand)==s.zeros(2,1);forcing.append(str(b+i*k))
    here=Path(__file__);out={'task_id':'G','cases':records,'forcing_operator_eigenvalues':forcing,
        'forcing_operator_determinant':str(s.factor(s.prod(b+i*k for i in range(3)))),
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),
        'scope':'Complete 80-coefficient horizontal, nilpotent-tangent and nilpotent-invariant cubic one-form spaces for every nonzero nilpotent Jordan type; projective closed incidence and geometric classification of these spaces are separate proof obligations.',
        'seconds':time.monotonic()-start}
    (here.resolve().parents[1]/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_BOUNDARY_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps({'ranks':[(c['jordan_type'],c['constraint_rank'],c['kernel_dimension']) for c in records],'forcing':forcing,'seconds':out['seconds']}))

def main():
    if '--nilpotent-incidence' in sys.argv:return nilpotent_incidence()
    started=time.monotonic();signal.alarm(30);resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    x,y,z,w=s.symbols('x y z w');xx=(x,y,z,w);pairs=list(itertools.combinations(range(4),2))
    Q=s.Matrix([5*x*x,5*x*y,5*x*z,y*y+3*x*z,y*z+7*x*w,7*y*w-3*z*z])
    def pf(p):return s.expand(p[0]*p[5]-p[1]*p[4]+p[2]*p[3])
    assert pf(Q)==0
    def mat(p):
        a=s.zeros(4)
        for (i,j),v in zip(pairs,p):a[i,j]=v;a[j,i]=-v
        return a
    def div(p):
        a=mat(p);return [sum(s.diff(a[i,j],xx[i]) for i in range(4)) for j in range(4)]
    assert not any(map(s.expand,div(Q)))
    ll=s.symbols('l0:24');L=s.Matrix([sum(ll[4*i+j]*xx[j] for j in range(4)) for i in range(6)])
    eq=div(L)+s.Poly(pf(Q+L)-pf(Q)-pf(L),*xx).coeffs();M,_=s.linear_eq_to_matrix(eq,ll)
    _,cols=M.rref();_,rows=M.T.rref();K=s.Matrix.hstack(*M.nullspace())
    T=s.Matrix([s.diff(p,v) for p in Q for v in xx]).jacobian(xx)
    assert M*K==s.zeros(M.rows,K.cols) and M*T==s.zeros(M.rows,4)
    _,trrows=T.T.rref();_,trcols=T.rref();cuts=s.zeros(len(trrows),24)
    for i,j in enumerate(trrows):cuts[i,j]=1
    normals=(M.col_join(cuts)).nullspace();NK=s.Matrix.hstack(*normals) if normals else s.zeros(24,0)
    nn=s.symbols('n0:'+str(NK.cols));coeff=NK*s.Matrix(nn) if nn else s.zeros(24,1)
    Ln=s.Matrix([sum(coeff[4*i+j]*xx[j] for j in range(4)) for i in range(6)])
    cc=s.symbols('c0:6');C=s.Matrix(cc);mon2=[xx[i]*xx[j] for i in range(4) for j in range(i,4)]
    G,_=s.linear_eq_to_matrix([s.Poly(pf(Q+C)-pf(Q)-pf(C),*xx).coeff_monomial(v) for v in mon2],cc)
    _,gc=G.rref();_,gr=G.T.rref();free=[i for i in range(6) if i not in gc];aa=s.symbols('a0:'+str(len(free)))
    Cs=s.zeros(6,1)
    for j,a in zip(free,aa):Cs[j]=a
    rhs=-s.Matrix([s.Poly(pf(Ln),*xx).coeff_monomial(mon2[i]) for i in gr])
    if aa:rhs-=G.extract(gr,free)*s.Matrix(aa)
    sol=G.extract(gr,gc).inv()*rhs
    for j,a in zip(gc,sol):Cs[j]=s.factor(a)
    allpf=s.Poly(pf(Q+Ln+Cs),*xx)
    residual=[{'monomial':list(mon),'coefficient':str(s.factor(v))} for mon,v in allpf.terms() if sum(mon) and v]
    here=Path(__file__);out={'task_id':'G','highest_source':'Pym arXiv:1403.6444v1, Section3.6 equation13, unique principal E3 orbit','highest_bivector':list(map(str,Q)),
      'full_linear_matrix':M.tolist(),'full_linear_rank':len(cols),'full_linear_pivot_rows':list(rows),'full_linear_pivot_columns':list(cols),'full_linear_pivot_determinant':str(M.extract(rows,cols).det()),'full_linear_kernel':K.tolist(),
      'translation_basis':T.tolist(),'translation_rank':len(trrows),'translation_cut_rows':list(trrows),'translation_pivot_columns':list(trcols),'translation_pivot_determinant':str(T.extract(trrows,trcols).det()),
      'normal_dimension':NK.cols,'normal_linear_bivector':list(map(str,Ln)),
      'constant_equation_matrix':G.tolist(),'constant_pivot_rows':list(gr),'constant_pivot_columns':list(gc),'constant_pivot_determinant':str(G.extract(gr,gc).det()),'constant_source':list(map(str,Cs)),
      'all_remaining_positive_degree_Pfaffian_coefficients':residual,'constant_Pfaffian':str(s.factor(allpf.coeff_monomial(1))),
      'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-started,
      'scope':'All24 linear and6 constant coefficients retained, modulo only the complete translation image with certified minor. No statement about E3 boundary follows from this orbit calculation alone.'}
    root=here.resolve().parents[1];(root/'reports/CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps({k:out[k] for k in ['full_linear_rank','translation_rank','normal_dimension','normal_linear_bivector','constant_source','all_remaining_positive_degree_Pfaffian_coefficients','constant_Pfaffian','seconds']}))
if __name__=='__main__':main()
