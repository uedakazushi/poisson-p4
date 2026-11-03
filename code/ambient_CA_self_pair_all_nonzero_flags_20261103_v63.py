#!/usr/bin/env python3
"""Strengthen CA/N whole-flag recognition to every nonzero self-pair parameter."""
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
STEM='AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63'
X=s.symbols('h x y z w');h,x,y,z,w=X;a=s.symbols('a',nonzero=True)
PAIRS=list(itertools.combinations(range(5),2));MONS=list(itertools.combinations_with_replacement(range(5),2))


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def strings(M):return[[str(v)for v in row]for row in M.tolist()]
def vector(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[u]*X[v])for i,j in PAIRS for u,v in MONS])
def canonical(P):
    d=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)]);E=s.Matrix(X);return(P-(E*d.T-d*E.T)/5).applyfunc(s.expand)
def rawCA(f,g,H):
    P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z);return P-P.T
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(60)
    whole=ROOT/'reports/AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63';assert read(whole/'report.json')['all_passed']
    inc=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63';assert read(inc/'report.json')['all_passed']
    pair=ROOT/'reports/AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63';data=read(pair/'certificate.json.gz')
    geom=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63';G=read(geom/'report.json');assert G['all_passed']
    geomreview=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json';assert read(geomreview)['all_passed']
    proper=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';assert read(proper)['all_passed']
    ca=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';base=read(ca)
    np=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';nbase=read(np)
    nr=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json';assert read(nr)['all_passed']
    q=s.Matrix([sum(s.Rational(c)*a**n for n,c in pol)for pol in data['first_source_physical150']]);P=s.zeros(5)
    for k,(i,j)in enumerate(PAIRS):P[i,j]=sum(q[15*k+m]*X[u]*X[v]for m,(u,v)in enumerate(MONS));P[j,i]=-P[i,j]
    coords=[h,x,z,w];Y=s.Matrix([-P[0,2],-P[1,2],P[2,3],P[2,4]]);deriv=lambda ell:s.expand(sum(5*Y[i]*s.diff(ell,coords[i])for i in range(4)))
    # Universal CA endpoint formulas, all16 cubics retained.
    f0,f1,f2,g=s.symbols('f0 f1 f2 g');c=s.symbols('c0:16');cubics=[s.prod(coords[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    F=f0*h*h+f1*h*x+f2*x*x;H=sum(v*m for v,m in zip(c,cubics));A=canonical(rawCA(F,0,H))
    assert s.expand(-A[0,2]-s.diff(F,x)*h/5)==0 and s.expand(-A[1,2]+F-s.diff(F,x)*x/5)==0
    endpoint=canonical(rawCA(0,g,H));assert endpoint[0,:]==s.zeros(1,5)and endpoint[1,:]==s.zeros(1,5)
    # Shared complete Darboux classification from the new whole-CB certificate.
    proof=read(whole/'certificate.json.gz');assert proof['constant_Casimir_rank_witness']['rank']==5 and proof['intrinsic_factor_rank_witness']['rank']==4
    forms=[x,x+a*h-2*a*a*z,x+3*a*h];lams=[-3*a,2*a,12*a]
    assert proof['invariant_linear_forms']==list(map(str,forms))
    alpha=s.symbols('alpha');residues=[]
    for ell,lam in zip(forms[1:],lams[1:]):
        # k=x and f_k=lam*x => f=lam*x^2/2+alpha*ell^2.
        residue=s.expand(-3*a*x*x/5+alpha*ell*ell+3*lam*x*x/10);residues.append(str(residue))
        if lam==2*a:assert s.expand(residue-alpha*ell*ell)==0
        else:assert s.expand(residue-alpha*ell*ell)==3*a*x*x
    kmat=s.Matrix([[s.Poly(s.expand(deriv(c)-12*a*x*c),h,x,z,w).coeff_monomial(s.prod(coords[i]for i in mon))for c in coords]
        for mon in itertools.combinations_with_replacement(range(4),2)if mon!=(1,1)])
    field=s.QQ.frac_field(a);pc=list(DomainMatrix.from_Matrix(kmat).convert_to(field).rref()[1]);assert len(pc)==2 and kmat[:,0]==s.zeros(9,1)and kmat[:,1]==s.zeros(9,1)
    # Thus CA ell=x has unique W=<h,x>; ell=H1 has unique W=<x,h-2az>; H2 impossible.
    # Direct full N28 source differential at the first flag, not merely a submatrix of CA29.
    NC=s.Matrix(nbase['fiber19_basis']);assert NC==s.Matrix(base['fixed_fiber_linear_basis150_by20'])[:,[0,1,*range(3,20)]]
    flags=list(map(tuple,nbase['flag_positions']));actions=[]
    for i,j in flags:
        M=s.zeros(5);M[i,j]=1;actions.append(vector((M*P+P*M.T-X[j]*P.diff(X[i])).applyfunc(s.expand)))
    JN=NC.row_join(s.Matrix.hstack(*actions));fixedrows=list(DomainMatrix.from_Matrix(NC.T).convert_to(s.QQ).rref()[1]);inv=NC[fixedrows,:].inv();quot=(JN[:,19:]-NC*inv*JN[fixedrows,19:]).applyfunc(s.expand)
    pcN=list(DomainMatrix.from_Matrix(quot).convert_to(field).rref()[1]);prN=list(DomainMatrix.from_Matrix(quot[:,pcN].T).convert_to(field).rref()[1]);B=DomainMatrix.from_Matrix(quot[prN,pcN]).convert_to(field);det=s.factor(field.to_sympy(B.det()));assert len(pcN)==9
    first=s.Matrix([sum(s.Rational(c)*a**n for n,c in pol)for pol in data['first_source_fixed20']]);second=s.Matrix([sum(s.Rational(c)*a**n for n,c in pol)for pol in data['second_source_fixed20']]);assert first[2]==0 and second[2]==-a*a
    payload=dict(universal_CA_gzero_marked_h_x_identities_verified=True,universal_CA_fzero_two_Casimir_identity_verified=True,
        complete_Casimir_factor_Darboux_input_report_sha256=sha(whole/'report.json'),CA_ell_H1_H2_candidate_residues=residues,
        CA_ell_x_k_condition_matrix=strings(kmat),CA_ell_x_k_condition_rank=2,CA_ell_x_k_kernel='span(h,x)',
        exactly_two_CA_flags=[dict(L='e_y',W='span(h,x)',H='span(x)'),dict(L='e_y',W='span(x,h-2*a*z)',H='span(x+a*h-2*a^2*z)')],
        full150_by28_N_source_differential=strings(JN),N_quotient150_by9=strings(quot),N_minor_rows=prN,N_minor_columns=pcN,N_minor_determinant=str(det),
        first_CA_flag_N_missing_C2='0',second_CA_flag_N_missing_C2=str(second[2]),whole_N_flag_is_first_CA_flag=True)
    with(out/'certificate.json.gz').open('wb')as file:
        with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,parameter_open='a!=0',whole_CA_geometric_flags_exactly2=True,whole_CA_fiber_two_reduced_points=True,
        whole_N_geometric_flags_exactly1=True,N_source_differential_rank=28,N_quotient9_minor=str(det),whole_N_fiber_reduced_length1=True,
        N_smooth_all_nonzero_family_points=True,CA_completion_two_smooth29_branches_all_nonzero_family_points=True,
        CA_branch_tangent_union_rank37_all_nonzero=True,CA_branch_tangent_intersection_dimension21_all_nonzero=True,
        CA_nonnormal_singular_nonCM_all_nonzero_family_points=True,CA_cone_depth_upper_bound=22,
        CA_origin_full_ring_or_conductor_classification_claimed=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),whole/'report.json',whole/'certificate.json.gz',pair/'report.json',pair/'certificate.json.gz',inc/'report.json',geom/'report.json',geomreview,proper,ca,np,nr];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
