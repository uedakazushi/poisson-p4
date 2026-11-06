#!/usr/bin/env python3
"""All thirty low-degree coefficients over the unique E3 highest-part orbit."""
import hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s

def main():
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
