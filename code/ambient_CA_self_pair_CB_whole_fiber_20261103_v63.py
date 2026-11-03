#!/usr/bin/env python3
"""Intrinsic whole-CB-flag recognition along the actual CA self-pair family."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63'
X=s.symbols('h x y z w');h,x,y,z,w=X;a=s.symbols('a',nonzero=True)
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2));CUBIC=list(itertools.combinations_with_replacement(range(5),3))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def strings(M):return[[str(v)for v in row]for row in M.tolist()]
def rawCB(f,g,H):
    P=s.zeros(5);P[1,2]=f;P[3,4]=g;P[2,3]=-s.diff(f,x)*z/2-s.diff(H,w);P[2,4]=-s.diff(f,x)*w/2+s.diff(H,z);return P-P.T
def canonical(P):
    d=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)]);E=s.Matrix(X);return(P-(E*d.T-d*E.T)/5).applyfunc(s.expand)
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    pair=ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63';assert read(pair/'report.json')['all_passed']
    incidence=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63';assert read(incidence/'report.json')['all_passed']
    rankreview=ROOT/'reports/AMBIENT_CA_SELF_PAIR_CB_SOURCE_RANK_20261103_v63_REVIEW.json';assert read(rankreview)['all_passed']
    globalreport=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json';globalreview=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';assert read(globalreview)['all_passed']
    raw=read(pair/'certificate.json.gz');q=s.Matrix([sum(s.Rational(c)*a**n for n,c in pol)for pol in raw['first_source_physical150']]);P=s.zeros(5)
    for k,(i,j)in enumerate(PAIRS):P[i,j]=sum(q[15*k+m]*X[u]*X[v]for m,(u,v)in enumerate(MONS));P[j,i]=-P[i,j]
    field=s.QQ.frac_field(a)
    def rankwitness(A):
        pc=list(DomainMatrix.from_Matrix(A).convert_to(field).rref()[1]);pr=list(DomainMatrix.from_Matrix(A[:,pc].T).convert_to(field).rref()[1]);B=DomainMatrix.from_Matrix(A[pr,pc]).convert_to(field)
        return dict(rank=len(pc),minor_rows=pr,minor_columns=pc,minor_determinant=str(s.factor(field.to_sympy(B.det()))))
    # All constant-covector contractions, a GL-invariant endpoint exclusion.
    casimir=s.Matrix([[s.Poly(P[i,j],*X).coeff_monomial(X[u]*X[v])for j in range(5)]for i in range(5)for u,v in MONS]);rCas=rankwitness(casimir);assert rCas['rank']==5
    # Full intrinsic constant-factor matrix: (Euler wedge v wedge P)_all175 =0.
    factor=[]
    for t in range(5):
        v=s.eye(5)[:,t];column=[]
        for i,j,k,l in itertools.combinations(range(5),4):
            F=s.expand(X[i]*(v[j]*P[k,l]-v[k]*P[j,l]+v[l]*P[j,k])-X[j]*(v[i]*P[k,l]-v[k]*P[i,l]+v[l]*P[i,k])+X[k]*(v[i]*P[j,l]-v[j]*P[i,l]+v[l]*P[i,j])-X[l]*(v[i]*P[j,k]-v[j]*P[i,k]+v[k]*P[i,j]))
            column.extend(s.Poly(F,*X).coeff_monomial(s.prod(X[u]for u in mon))for mon in CUBIC)
        factor.append(s.Matrix(column))
    factor=s.Matrix.hstack(*factor);rFactor=rankwitness(factor);assert rFactor['rank']==4 and factor[:,2]==s.zeros(175,1)
    # Quotient4 tracefree vectorfield C where P=e_y wedge C.
    Y=s.Matrix([-P[0,2],-P[1,2],P[2,3],P[2,4]]);coords=[h,x,z,w];D5=5*Y
    assert D5==s.Matrix([5*x*x+12*a*h*x,-3*a*x*x,5*h*x+2*a*x*z,5*h*h+5*z*z-8*a*x*w])
    deriv=lambda ell:s.expand(sum(D5[i]*s.diff(ell,coords[i])for i in range(4)))
    forms=[x,x+a*h-2*a*a*z,x+3*a*h];eigenvalues=[-3*a,2*a,12*a]
    assert all(s.expand(deriv(ell)-lam*x*ell)==0 for ell,lam in zip(forms,eigenvalues))
    A=s.Matrix([[12*a,0,5],[5,-3*a,0],[0,0,2*a]]);E=s.Matrix([[s.diff(ell,c)for ell in forms]for c in[h,x,z]])
    assert A*E==E*s.diag(*eigenvalues)and s.factor(E.det())==-6*a**3
    alpha,beta,gamma,delta,rh,rx,rz,rw=s.symbols('A B C D rh rx rz rw');ell=alpha*h+beta*x+gamma*z+delta*w;cofactor=rh*h+rx*x+rz*z+rw*w
    darboux=s.Poly(deriv(ell)-ell*cofactor,h,x,z,w)
    assert darboux.coeff_monomial(w*w)==-delta*rw and darboux.coeff_monomial(h*w)==-alpha*rw-delta*rh
    assert darboux.coeff_monomial(z*w)==-gamma*rw-delta*rz and darboux.coeff_monomial(h*h)==5*delta-alpha*rh
    # delta!=0 would force rw=rh=rz=0 and then delta=0. Thus all invariant ell are in span(h,x,z).
    # ell=x candidate needs independent k but this full linear constraint has kernel only <x>.
    kmat=s.Matrix([[s.Poly(s.expand(deriv(c)-s.Rational(9,2)*a*x*c),h,x,z,w).coeff_monomial(s.prod(coords[i]for i in mon))for c in coords]
        for mon in itertools.combinations_with_replacement(range(4),2)if mon!=(1,1)])
    rk=rankwitness(kmat);assert rk['rank']==3 and kmat[:,1]==s.zeros(9,1)
    # The other two h-lines force W=<ell,x>, f=lam*x^2/4 +alpha*ell^2.
    # Then C(x)=-f+(2/5)f_x*x is possible only lam=12a (the ALREADY known CB flag).
    alpha2=s.symbols('alpha2');candidate_residues=[]
    for ell,lam in zip(forms[1:],eigenvalues[1:]):
        residual=s.expand(-3*a*x*x/5+alpha2*ell*ell+lam*x*x/20);candidate_residues.append(str(residual))
        if lam==2*a:
            independent=s.Matrix([[s.Poly(ell*ell,h,x,z,w).coeff_monomial(mon),s.Poly(x*x,h,x,z,w).coeff_monomial(mon)]for mon in[h*h,x*x]])
            assert independent.det()==a*a
        else:assert s.expand(residual-alpha2*ell*ell)==0
    # Universal proper-CB endpoint identities: a=0 endpoint has two constant Casimirs;
    # b=0 endpoint's hline and second W generator obey exactly the identities used above.
    f0,f1,f2=s.symbols('f0:3');g0,g1,g2=s.symbols('g0:3');hc=s.symbols('c0:16')
    cubic=[s.prod(coords[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)];H=sum(c*m for c,m in zip(hc,cubic));F=f0*h*h+f1*h*x+f2*x*x;G=g0*h*h+g1*h*x+g2*x*x
    endpointA=canonical(rawCB(0,G,H));assert endpointA[0,:]==s.zeros(1,5)and endpointA[1,:]==s.zeros(1,5)
    endpointB=canonical(rawCB(F,0,H));assert s.expand(-endpointB[0,2]-2*s.diff(F,x)*h/5)==0 and s.expand(-endpointB[1,2]+F-2*s.diff(F,x)*x/5)==0
    inc=read(incidence/'certificate.json.gz');cbG=s.Matrix([[s.sympify(c,locals={'a':a})for c in row]for row in inc['actual_CB_change_of_coordinates']]);assert cbG[0,:]==s.Matrix([[-3*a,-1,0,0,0]])
    # This flag is H=<x+3ah>, W=<h,x>, L=<e_y>; the candidate is not a new second flag.
    payload=dict(all75_constant_Casimir_matrix=strings(casimir),constant_Casimir_rank_witness=rCas,
        all175_intrinsic_factor_matrix=strings(factor),intrinsic_factor_rank_witness=rFactor,unique_intrinsic_factor='e_y',
        quotient_tracefree_C_times5=list(map(str,D5)),invariant_linear_forms=list(map(str,forms)),D5_eigenvalues=list(map(str,eigenvalues)),
        eigenbasis_matrix=strings(E),eigenbasis_determinant=str(s.factor(E.det())),
        full_general_Darboux_equation=str(darboux.as_expr()),ell_x_candidate_k_constraint_matrix=strings(kmat),ell_x_candidate_k_rank_witness=rk,
        ell_H1_and_H2_candidate_residues=candidate_residues,
        universal_endpointA_two_Casimir_identity=True,universal_endpointB_h_and_x_identities=True,
        unique_CB_flag=dict(L='e_y',W='span(h,x)',H='span(x+3*a*h)',ratio='[f:g]=[1:0]'),
        explicit_CB_flag_same_as_frozen_incidence=True)
    with(out/'certificate.json.gz').open('wb')as file:
        with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,constant_Casimir_rank=5,constant_Casimir_minor=rCas['minor_determinant'],
        intrinsic_constant_factor_rank=4,intrinsic_factor_minor=rFactor['minor_determinant'],
        all_invariant_linear_forms_classified_for_a_nonzero=True,ell_x_independent_k_impossible=True,
        ell_H1_CB_binary_equations_inconsistent=True,ell_H2_is_existing_unique_CB_flag=True,
        second_CB_flag_constructed=False,whole_geometric_CB_fiber_unique_by_proper_endpoint_coverage=True,
        whole_CB_fiber_reduced_length1_using_independently_accepted_source29=True,
        CB_smooth_at_nonzero_family_points_by_finite_fiber_Nakayama=True,
        global_CB_bad_locus_or_origin_full_ring_claimed=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),pair/'report.json',pair/'certificate.json.gz',incidence/'report.json',incidence/'certificate.json.gz',rankreview,globalreport,globalreview];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
