#!/usr/bin/env python3
"""Two whole-source orbit limits through the covered regular rank3-quadric locus."""
import hashlib,itertools,json,resource,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
x=s.symbols('h x1 x2 x3 x4');h,b,c,r,z=x;ts=s.symbols('tau0:2');t=s.Symbol('t');E=s.Matrix(x)
ENV={str(v):v for v in x+ts+(t,)}
def decode(data):
    P=s.zeros(5)
    for key,v in data.items():i,j=map(int,key.split(','));P[i,j]=s.sympify(v,locals=ENV)
    return P-P.T

def encode(P):return {f'{i},{j}':str(s.factor(P[i,j])) for i,j in itertools.combinations(range(5),2) if P[i,j]!=0}
def rawjac(P):
    return {ij:s.expand(sum(P[a,l]*s.diff(P[b,c],x[l]) for a,b,c in [ij,(ij[1],ij[2],ij[0]),(ij[2],ij[0],ij[1])] for l in range(5))) for ij in itertools.combinations(range(5),3)}

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(60,61));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();records={};inputs=[]
    specs=[('CZERO_NONZERO_I',s.Rational(5,2),144*(2*h*c+b*b)-144*t*h*r+48*t*b*c+8*t*t*c*c-t**4*r*r,2),
           ('CZERO_NONZERO_II',s.Rational(1,48),186624*(2*h*c+b*b)+31104*t*h*r-10368*t*b*c+288*t*t*c*c-t**4*r*r,1)]
    for label,k,q,lam in specs:
        source=ROOT/('reports/CLASSIFICATION_Q3_ANIL3_'+label+'_WEIGHTED_20261107.json');inputs.append(source);saved=json.loads(source.read_text())
        P=decode(saved['full_curve_brackets']).subs({ts[0]:t,ts[1]:k*t*t},simultaneous=True).applyfunc(s.expand)
        assert (P.subs(t,0)-decode(saved['canonical_brackets'])).applyfunc(s.expand)==s.zeros(5)
        J=rawjac(P)
        assert all(s.expand(sum((-1)**i*x[a]*J[inds[:i]+inds[i+1:]] for i,a in enumerate(inds)))==0 for inds in itertools.combinations(range(5),4))
        div=s.Matrix([sum(s.diff(P[i,j],x[i]) for i in range(5)) for j in range(5)])
        Pi=(P-(E*div.T-div*E.T)/5).applyfunc(s.expand)
        assert all(s.expand(sum(s.diff(Pi[i,j],x[i]) for i in range(5)))==0 for j in range(5))
        assert all(v==0 for v in rawjac(Pi).values())
        cof=[]
        for j in range(5):
            val=s.expand(sum(s.diff(q,x[i])*Pi[i,j] for i in range(5)));quo,rem=s.div(val,q,*x)
            assert rem==0;cof.append(s.expand(quo))
        A=s.Matrix(cof).jacobian(x);assert s.Matrix(cof)==A*E
        L=s.Symbol('L');char=s.factor(A.charpoly(L).as_expr())
        assert s.expand(char-L*(L-lam*t)*(L+lam*t)*(L-2*lam*t)*(L+2*lam*t))==0
        H=s.hessian(q,x);idx=list(itertools.combinations(range(5),4))
        assert all(s.expand(H.extract(I,J).det())==0 for I in idx for J in idx)
        minor=s.factor(H[:3,:3].det());assert minor and not minor.has(t)
        pf=[]
        for i in range(5):
            a,b1,c1,d=[j for j in range(5) if j!=i]
            pf.append(s.expand((-1)**i*(Pi[a,b1]*Pi[c1,d]-Pi[a,c1]*Pi[b1,d]+Pi[a,d]*Pi[b1,c1])))
        D=s.expand(sum(v*p for v,p in zip(x,pf))/5);C,rem=s.div(D,q,*x);assert rem==0;C=s.expand(C)
        assert all(s.expand(v-s.diff(D,xi))==0 for v,xi in zip(pf,x))
        assert s.Poly(C,*x,t).total_degree()>=3 and all(sum(e[:5])==3 for e,v in s.Poly(C,*x,t).terms())
        records[label]={'source_weighted_report':str(source.relative_to(ROOT)),'source_orbit':'Czero_nonzero_'+label[-1] if label.endswith('_I') else 'Czero_nonzero_II',
          'normal_parameter_curve':{'tau0':'t','tau1':str(k*t*t)},'literal_actual_curve':encode(P),
          'all350_projective_Jacobi_coefficients_zero':True,'specialization_is_full_source_orbit_representative':True,
          'standard_unimodular_curve':encode(Pi),'all10_strict_Jacobi_zero_after_gauge':True,
          'quadric':str(q),'quadric_rank_all_t':3,'quadric_nonzero_3_minor':str(minor),'all_4_minors_zero':True,
          'canonical_cofactor':list(map(str,cof)),'all5_q_normal_equations_zero':True,'cofactor_characteristic_polynomial':str(char),
          'cofactor_eigenvalues':[str(a*t) for a in [0,lam,-lam,2*lam,-2*lam]],'actual_cubic':str(s.factor(C)),
          'actual_divisor_factorization':'D=quadric*actual_cubic','all5_actual_Pfaffian_gradient_identities':True,
          'specialized_cubic':str(s.factor(C.subs(t,0))),
          'generic_cubic_irreducible_reason':'At t=0 this is a nonzero scalar multiple of the universally proved irreducible canonical source cubic. Reducible homogeneous cubics form a closed subset, so a Zariski-open neighborhood of0 retains irreducibility.',
          'generic_membership':'For nonzero t in this open set, cofactor is regular semisimple with eigenvalue ratio2 and quadric rank3. The covered rank3 regular support theorem gives MPS membership; closedness gives t=0 and the whole orbit.'}
    result={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
      'input_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
      'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_source':Path(__file__).read_text(),
      'recognition_dependency':'research/classification/DEGENERATE_QUADRIC_CUBIC_REGULAR_SUPPORT.md',
      'whole_source_dependency':'research/classification/DEGENERATE_QUADRIC_Q3_ANIL3_SUPPORT.md',
      'seconds':time.monotonic()-start}
    out=ROOT/'reports/CLASSIFICATION_Q3_ANIL3_CURVES_20261107.json'
    with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':result['seconds'],'whole_orbits_covered':2}))
if __name__=='__main__':main()
