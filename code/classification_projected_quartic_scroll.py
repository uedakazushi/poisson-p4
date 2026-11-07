#!/usr/bin/env python3
"""Whole projective rankzero sources on projected quartic scroll models.

The generic S22 nodal projection is one actual source-preserving group orbit.
All150 homogeneous bivector coefficients and full Euler quotient are retained;
the conormal bracket uses its correct projective dehomogenization numerator.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
from fractions import Fraction as Q
import argparse,hashlib,json,resource,signal,time
import sympy as S
from sympy.polys.matrices import DomainMatrix
from classification_projected_veronese_zero_review import monomials,plus,scale,sparse_rref
ROOT=Path(__file__).resolve().parents[1];PAIRS=list(combinations(range(5),2));TRIPLES=list(combinations(range(5),3));QUADS=list(combinations(range(5),4))


def bimon(first,second):return [(i,first-i,j,second-j)for i in range(first,-1,-1)for j in range(second,-1,-1)]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--surface',choices=['s22_node','s22_cusp','s13_node','s13_cusp'],default='s22_node');ap.add_argument('--factors',action='store_true');args=ap.parse_args()
    if args.factors:return factors()
    destination=ROOT/('reports/CLASSIFICATION_PROJECTED_QUARTIC_'+args.surface.upper()+'_20261107.json')
    if destination.exists():raise FileExistsError(destination)
    start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    x=S.Matrix(S.symbols('x0:5'));ss,tt,u,v,rr=S.symbols('ss tt u v rr');param=(ss,tt,u,v)
    if args.surface=='s22_node':
        q=S.Matrix([ss*u*u-tt*v*v,ss*u*v,ss*v*v,tt*u*u,tt*u*v]);F2=x[1]*x[4]-x[2]*x[3];F3=x[0]*x[2]*x[3]-x[1]**2*x[3]+x[4]**2*x[2]
    elif args.surface=='s22_cusp':
        q=S.Matrix([ss*u*u-tt*u*v,ss*u*v,ss*v*v,tt*u*u,tt*v*v]);F2=x[0]*x[2]-x[1]**2+x[1]*x[4];F3=x[1]**2*x[4]-x[2]**2*x[3]
    elif args.surface=='s13_node':
        param=(ss,tt,rr);q=S.Matrix([rr*ss,rr*tt,ss**3-tt**3,ss*ss*tt,ss*tt*tt]);F2=x[0]*x[4]-x[1]*x[3];F3=x[2]*x[3]*x[4]-x[3]**3+x[4]**3
    else:
        param=(ss,tt,rr);q=S.Matrix([rr*ss,rr*tt,ss**3,ss*tt*tt,tt**3]);F2=x[0]*x[4]-x[1]*x[3];F3=x[2]*x[4]**2-x[3]**3
    def parameter_monomials(coordinate_degree):
        if args.surface.startswith('s22'):return bimon(coordinate_degree,2*coordinate_degree)
        degree=3*coordinate_degree
        return sorted([(a,degree-2*c-a,c)for c in range(degree//2+1)for a in range(degree-2*c+1)],reverse=True)
    sub=dict(zip(x,q));columns=[(ij,e)for ij in PAIRS for e in monomials(5,2)]
    labels=[(tri,e)for tri in TRIPLES for e in parameter_monomials(3)];index={label:i for i,label in enumerate(labels)};M=S.zeros(len(labels),150)
    for j,(ij,e)in enumerate(columns):
        coefficient=S.prod(q[i]**n for i,n in enumerate(e))
        for tri in TRIPLES:
            if set(ij).issubset(tri):
                k=next(k for k in tri if k not in ij)
                for exponent,c in S.Poly(S.expand((-1)**tri.index(k)*q[k]*coefficient),*param).terms():M[index[tri,exponent],j]=c
    dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(S.QQ)
    rank=lambda A:len(dm(A).rref()[1]);kernel=dm(M).nullspace().to_Matrix().T;assert M*kernel==S.zeros(M.rows,kernel.cols)
    def decode(coefficients):
        P=S.zeros(5)
        for c,(ij,e)in zip(coefficients,columns):P[ij]+=c*S.prod(xx**n for xx,n in zip(x,e))
        for i,j in PAIRS:P[j,i]=-P[i,j]
        return P.applyfunc(S.expand)
    def vector(P):return S.Matrix([S.Poly(S.expand(P[ij]),*x).coeff_monomial(e)for ij,e in columns])
    eulers=[]
    for a in range(5):
        for b in range(5):
            field=S.zeros(5,1);field[a]=x[b];eulers.append(vector(x*field.T-field*x.T))
    Euler=S.Matrix.hstack(*eulers);assert rank(Euler)==24 and M*Euler==S.zeros(M.rows,25)
    canonical=[]
    for j in range(kernel.cols):
        P=decode(kernel[:,j]);div=S.Matrix([sum(S.diff(P[i,k],x[i])for i in range(5))for k in range(5)])
        canonical.append(vector(P-(x*div.T-div*x.T)/5))
    C=S.Matrix.hstack(*canonical);red,piv=dm(C.T).rref();basis=red.to_Matrix()[:len(piv),:].T
    assert rank(Euler.row_join(basis))==kernel.cols and basis.cols+24==kernel.cols
    assert S.expand(F2.subs(sub,simultaneous=True))==0 and S.expand(F3.subs(sub,simultaneous=True))==0
    d2=S.Matrix([S.diff(F2,xx)for xx in x]);d3=S.Matrix([S.diff(F3,xx)for xx in x]);d0=S.Matrix([1,0,0,0,0]);connection_labels=[(i,e)for i in range(5)for e in parameter_monomials(5)];conn=S.zeros(len(connection_labels),basis.cols);ci={label:i for i,label in enumerate(connection_labels)}
    numerators=[]
    for j in range(basis.cols):
        P=decode(basis[:,j]);H=(d2.T*P*d3)[0];H20=(d2.T*P*d0)[0];H03=(d0.T*P*d3)[0];numerator=S.expand(x[0]*H-3*F3*H20-2*F2*H03);numerators.append(str(numerator))
        assert S.expand(numerator.subs(sub,simultaneous=True))==0
        for i in range(5):
            for e,c in S.Poly(S.expand(S.diff(numerator,x[i]).subs(sub,simultaneous=True)),*param).terms():
                if c:conn[ci[i,e],j]=c
    needed=dm(conn).nullspace().to_Matrix().T;assert conn*needed==S.zeros(conn.rows,needed.cols);final=basis*needed;d=final.cols;z=S.symbols('z0:'+str(d));P=decode(final*S.Matrix(z))
    print('full150 zero rank',rank(M),'kernel',kernel.cols,'projective',basis.cols,'abelian rank',rank(conn),'final source',d,flush=True)
    strict={tri:S.expand(sum(P[a,l]*S.diff(P[b,c],x[l])for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}
    jac={quad:S.expand(sum((-1)**i*x[quad[i]]*strict[quad[:i]+quad[i+1:]]for i in range(4)))for quad in QUADS};J=[S.Poly(jac[quad],*x).coeff_monomial(e)for quad in QUADS for e in monomials(5,4)]
    zcols=list(combinations_with_replacement(range(d),2));JJ=[]
    for p in J:JJ.append({k:Q(c)for k,(i,j)in enumerate(zcols)if(c:=S.Poly(p,*z).coeff_monomial(z[i]*z[j]))})
    actual=S.expand(sum((-1)**i*x[i]*(P[j,k]*P[l,m]-P[j,l]*P[k,m]+P[j,m]*P[k,l])for i in range(5)for j,k,l,m in [[t for t in range(5)if t!=i]]))
    DD=[{k:Q(c)for k,(i,j)in enumerate(zcols)if(c:=S.Poly(S.Poly(actual,*x).coeff_monomial(e),*z).coeff_monomial(z[i]*z[j]))}for e in monomials(5,5)]
    rb,proof,chosen=sparse_rref(JJ,True);allmember=True;targetcert=[]
    for i,target in enumerate(DD):
        residual=dict(target);functional={}
        for p in sorted(rb):
            if p in residual:
                a=residual[p];residual=plus(residual,scale(rb[p],-a));functional=plus(functional,scale(proof[p],a))
        allmember &=not residual
        if target:targetcert.append({'actual_D_row':i,'original350_functional':functional,'residual':residual})
    chosenrows=list(dm(M.T).rref()[1]);chosencols=list(dm(M).rref()[1]);cr=list(dm(conn.T).rref()[1]);cc=list(dm(conn).rref()[1])
    group_dimension={'s22_node':1,'s22_cusp':2,'s13_node':2,'s13_cusp':3}[args.surface]
    out={'task_id':'G','model':args.surface,'surface_parametrization':list(map(str,q)),'parametrization_variables':list(map(str,param)),
        'full150_coefficient_columns':[{'pair':ij,'exponents':e}for ij,e in columns],'all280_zero_on_surface_rows':[{'triple':tri,'exponents':e}for tri,e in labels],
        'complete_zero_matrix':{str(i):{str(j):str(M[i,j])for j in range(150)if M[i,j]}for i in range(M.rows)},'zero_matrix_shape':list(M.shape),'zero_matrix_rank':rank(M),
        'zero_rank_minor':{'rows':chosenrows,'columns':chosencols,'determinant':str(M[chosenrows,chosencols].det())},'full_zero_kernel_basis':kernel.tolist(),'Euler_kernel_dimension':24,
        'complete150_by_projective_source_matrix':basis.tolist(),'projective_zero_source_dimension':basis.cols,'generic_ideal_generators':[str(F2),str(F3)],
        'dense_chart_inversions':['x0','x2','x3'],'excluded_chart_scope':'Only a necessary generic conormal condition is imposed; if it holds at genericS, polynomial partials vanishonwholeintegralS. No rank4candidate islost by these dense-chart choices.',
        'projective_conormal_bracket_numerators':numerators,'all330_conormal_firstpartial_labels':[{'ambient_partial':i,'exponents':e}for i,e in connection_labels],
        'complete_abelian_matrix':conn.tolist(),'abelian_matrix_shape':list(conn.shape),'abelian_matrix_rank':rank(conn),'abelian_rank_minor':{'rows':cr,'columns':cc,'determinant':str(conn[cr,cc].det())},
        'complete_projective_by_abelian_source_matrix':needed.tolist(),'complete150_by_abelian_source_matrix':final.tolist(),'complete_abelian_source_dimension':d,
        'full_reduced_source':{str(i)+','+str(j):str(P[i,j])for i,j in PAIRS if P[i,j]!=0},'all350_Jacobi':list(map(str,J)),'all350_parameter_quadric_rows':JJ,'all126_actual_D_parameter_quadric_rows':DD,
        'Jacobi_row_rank':len(rb),'all_actual_D_rows_in_full_Jacobi_rowspace':allmember,'actual_D':str(actual),'actual_D_membership_certificates':targetcert,
        'whole_fixed_source_dimension_upper_before_Jacobi':d,'source_preserving_group_dimension_lower':group_dimension,'whole_projective_PGL_image_dimension_upper_before_Jacobi':24+d-1-group_dimension,
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-start,'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':f'Complete projective rankzero and necessaryabelianconormal source for thisentire{args.surface}projectionorbit. Fullsource/Euler/allJacobi/actualD retained. Globalclassificationofotherprojectioncentreorbits andtheconeone-parameterfamily is separate.'}
    destination.write_text(json.dumps(out,indent=2,default=str)+'\n');print('final J rank',len(rb),'actualD membership',allmember,'preJac PGLbound',23+d-group_dimension,'seconds',out['seconds'],flush=True)


def factors():
    destination=ROOT/'reports/CLASSIFICATION_PROJECTED_QUARTIC_SCROLL_FACTORS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    start=time.monotonic();x=S.symbols('x0:5');records=[]
    for model in ['s22_node','s22_cusp','s13_node','s13_cusp']:
        path=ROOT/('reports/CLASSIFICATION_PROJECTED_QUARTIC_'+model.upper()+'_20261107.json');saved=json.loads(path.read_text());n=saved['complete_abelian_source_dimension'];z=S.symbols('z0:'+str(n));env={str(v):v for v in (*x,*z)}
        D=S.Poly(S.sympify(saved['actual_D'],locals=env),*x,*z,domain=S.QQ);Qpoly=S.Poly(S.sympify(saved['generic_ideal_generators'][0],locals=env),*x,*z,domain=S.QQ);quotient,remainder=S.div(D,Qpoly);assert remainder.is_zero
        rows=[{int(j):Q(c)for j,c in row.items()}for row in saved['all350_parameter_quadric_rows']];rb,proof,_=sparse_rref(rows,True);cols=list(combinations_with_replacement(range(n),2));cert=[]
        if model.startswith('s13'):
            for i in range(1,n):
                target={cols.index((i,i)):Q(1)};residual=dict(target);w={}
                for p in sorted(rb):
                    if p in residual:
                        a=residual[p];residual=plus(residual,scale(rb[p],-a));w=plus(w,scale(proof[p],a))
                assert not residual and plus(*(scale(rows[k],c)for k,c in w.items()))==target
                cert.append({'squared_parameter':str(z[i]),'original350_Jacobi_functional':{str(k):str(c)for k,c in w.items()}})
            assert D.as_expr().subs({v:0 for v in z[1:]})==0
        records.append({'model':model,'input_report':str(path.relative_to(ROOT)),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'unique_quadric_factor':str(Qpoly.as_expr()),'exact_cubic_quotient_of_actual_D':str(quotient.as_expr()),'exact_polynomial_division_remainder_zero':True,
            'irreducible_quintic_branch_empty_before_Jacobi':True,'S13_square_parameter_functionals':cert,'S13_actual_D_zero_on_all_reduced_Poisson_source':bool(cert)})
    out={'task_id':'G','all_four_model_branches_empty_for_irreducible_quintic':True,'cases':records,'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-start,
        'scope':'Exactfactorcertificatesconsumeeachimmutablefull150/11/conormal source report. Nooldproducerreportisrewritten. Completecentreorbitcoverageistheseparategeometricalnote; S04coneandotherdegree surfacesarenotincluded.'}
    destination.write_text(json.dumps(out,indent=2)+'\n');print('four exactquadric factors andS13squareJacfunctional PASS',out['seconds'])


if __name__=='__main__':main()
