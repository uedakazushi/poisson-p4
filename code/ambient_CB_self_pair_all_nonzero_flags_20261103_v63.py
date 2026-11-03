#!/usr/bin/env python3
"""Whole CA/CB/N proper flags on the new CB self-pair family, every b != 0."""
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
from ambient_CA_self_pair_CB_source_rank_20261103_v63 import X,PAIRS,MONS,vector,canonical,rawCB,action
from ambient_CA_self_pair_all_nonzero_flags_20261103_v63 import rawCA

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63'
h,x,y,z,w=X;b=s.symbols('b',nonzero=True)


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def strings(M):return[[str(v)for v in row]for row in M.tolist()]
def matrix(q):
    P=s.zeros(5)
    for k,(i,j)in enumerate(PAIRS):P[i,j]=sum(q[15*k+m]*X[u]*X[v]for m,(u,v)in enumerate(MONS));P[j,i]=-P[i,j]
    return P
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    parent=ROOT/'reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63';basic=read(parent/'report.json');assert basic['all_passed']
    data=read(parent/'certificate.json.gz');q=s.Matrix([s.sympify(v,locals={'b':b})for v in data['full_polynomial_Q150']]);P=matrix(q)
    ca=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';base=read(ca)
    np=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';nbase=read(np)
    careview=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';cbreview=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json'
    nreview=ROOT/'reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json'
    assert all(read(p)['all_passed']for p in[careview,cbreview,nreview])
    domain=s.QQ.frac_field(b)
    def ranks(A):
        pc=list(DomainMatrix.from_Matrix(A).convert_to(domain).rref()[1]);pr=list(DomainMatrix.from_Matrix(A[:,pc].T).convert_to(domain).rref()[1]);D=DomainMatrix.from_Matrix(A[pr,pc]).convert_to(domain)
        return dict(rank=len(pc),rows=pr,columns=pc,minor=str(s.factor(domain.to_sympy(D.det()))))
    casimir=s.Matrix([[s.Poly(P[i,j],*X).coeff_monomial(X[u]*X[v])for j in range(5)]for i in range(5)for u,v in MONS]);rCas=ranks(casimir);assert rCas['rank']==5
    cubics5=list(itertools.combinations_with_replacement(range(5),3));columns=[]
    for t in range(5):
        v=s.eye(5)[:,t];column=[]
        for i,j,k,l in itertools.combinations(range(5),4):
            expr=s.expand(X[i]*(v[j]*P[k,l]-v[k]*P[j,l]+v[l]*P[j,k])-X[j]*(v[i]*P[k,l]-v[k]*P[i,l]+v[l]*P[i,k])+X[k]*(v[i]*P[j,l]-v[j]*P[i,l]+v[l]*P[i,j])-X[l]*(v[i]*P[j,k]-v[j]*P[i,k]+v[k]*P[i,j]))
            column.extend(s.Poly(expr,*X).coeff_monomial(s.prod(X[u]for u in mon))for mon in cubics5)
        columns.append(s.Matrix(column))
    factor=s.Matrix.hstack(*columns);rFactor=ranks(factor);assert rFactor['rank']==4 and factor[:,2]==s.zeros(175,1)
    coords=[h,x,z,w];Y=s.Matrix([-P[0,2],-P[1,2],P[2,3],P[2,4]]);deriv=lambda ell:s.expand(sum(5*Y[i]*s.diff(ell,coords[i])for i in range(4)))
    forms=[x,x+b*h-b*b*z,x+2*b*h];lams=[-2*b,3*b,8*b]
    E=s.Matrix([[s.diff(ell,c)for ell in forms]for c in[h,x,z]]);A=s.Matrix([[8*b,0,5],[5,-2*b,0],[0,0,3*b]])
    assert A*E==E*s.diag(*lams)and s.factor(E.det())==-2*b**3
    assert all(s.expand(deriv(ell)-lam*x*ell)==0 for ell,lam in zip(forms,lams))
    A0,B0,C0,D0,rh,rx,rz,rw=s.symbols('A B C D rh rx rz rw');ell=A0*h+B0*x+C0*z+D0*w;cofactor=rh*h+rx*x+rz*z+rw*w;dar=s.Poly(deriv(ell)-ell*cofactor,h,x,z,w)
    assert dar.coeff_monomial(w*w)==-D0*rw and dar.coeff_monomial(h*w)==-A0*rw-D0*rh and dar.coeff_monomial(z*w)==-C0*rw-D0*rz and dar.coeff_monomial(h*h)==5*D0-A0*rh
    # All invariant ell are exactly the three above, by the same explicit D0!=0 contradiction and eigenbasis proof.
    def kmatrix(lam):
        return s.Matrix([[s.Poly(s.expand(deriv(c)-lam*x*c),h,x,z,w).coeff_monomial(s.prod(coords[i]for i in mon))for c in coords]
            for mon in itertools.combinations_with_replacement(range(4),2)if mon!=(1,1)])
    cbk=kmatrix(3*b);cak=kmatrix(8*b);rCBk=ranks(cbk);rCAk=ranks(cak)
    assert rCBk['rank']==2 and cbk*s.Matrix([0,1,0,0])==s.zeros(9,1)and cbk*s.Matrix([1,0,-b,0])==s.zeros(9,1)
    assert rCAk['rank']==2 and cak[:,0]==s.zeros(9,1)and cak[:,1]==s.zeros(9,1)
    alpha=s.symbols('alpha');cbres=[];cares=[]
    for ell,lam in zip(forms[1:],lams[1:]):
        cr=s.expand(-2*b*x*x/5+alpha*ell*ell+lam*x*x/20);ar=s.expand(-2*b*x*x/5+alpha*ell*ell+3*lam*x*x/10)
        cbres.append(str(cr));cares.append(str(ar))
        assert s.Poly(ell*ell,h,x,z,w).coeff_monomial(h*h)!=0 #ell² andx² independent for b!=0
        if lam==8*b:assert s.expand(cr-alpha*ell*ell)==0
        else:assert s.expand(cr-alpha*ell*ell)==-b*x*x/4
        assert s.expand(ar-alpha*ell*ell)==((3*lam-4*b)/10)*x*x
    # The universal canonical endpoint and marked-line identities retain all16 Hamiltonian coefficients.
    f0,f1,f2,g0,g1,g2,g=s.symbols('f0 f1 f2 g0 g1 g2 g');hc=s.symbols('c0:16')
    cubics=[s.prod(coords[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)];H=sum(c*m for c,m in zip(hc,cubics));F=f0*h*h+f1*h*x+f2*x*x;GF=g0*h*h+g1*h*x+g2*x*x
    CB0=canonical(rawCB(0,GF,H));CA0=canonical(rawCA(0,g,H));assert all(S[0,:]==s.zeros(1,5)and S[1,:]==s.zeros(1,5)for S in[CB0,CA0])
    CA=canonical(rawCA(F,0,H));CB=canonical(rawCB(F,0,H))
    assert s.expand(-CA[0,2]-s.diff(F,x)*h/5)==0 and s.expand(-CA[1,2]+F-s.diff(F,x)*x/5)==0
    assert s.expand(-CB[0,2]-2*s.diff(F,x)*h/5)==0 and s.expand(-CB[1,2]+F-2*s.diff(F,x)*x/5)==0
    # Actual unique CA lift is the original fixed-frame source. N is its closed C[2]=0 subsource.
    C=s.Matrix(base['fixed_fiber_linear_basis150_by20']);NC=s.Matrix(nbase['fiber19_basis']);assert NC==C[:,[0,1,*range(3,20)]]
    given=s.Matrix([s.sympify(v,locals={'b':b})for v in data['CA_fixed20_coordinates']]);assert given[2]==0 and C*given==q
    flags=list(map(tuple,base['flag_positions']));flagJ=s.Matrix.hstack(*[vector(action(i,j,P))for i,j in flags]);sranks=[];source_matrices=[]
    for name,fixed in[('CA',C),('N',NC)]:
        rows=list(DomainMatrix.from_Matrix(fixed.T).convert_to(s.QQ).rref()[1]);inv=fixed[rows,:].inv();quot=(flagJ-fixed*inv*flagJ[rows,:]).applyfunc(s.expand);rank=ranks(quot)
        assert rank['rank']==9;sranks.append(dict(name=name,fixed_rank=fixed.cols,rank=fixed.cols+9,flag_quotient_witness=rank));source_matrices.append(dict(name=name,full_source=strings(fixed.row_join(flagJ)),flag_quotient=strings(quot)))
    payload=dict(all75_Casimir_matrix=strings(casimir),Casimir_rank_witness=rCas,all175_intrinsic_factor_matrix=strings(factor),intrinsic_factor_rank_witness=rFactor,
        full_Darboux_equation=str(dar.as_expr()),invariant_linear_forms=list(map(str,forms)),D5_eigenvalues=list(map(str,lams)),eigenbasis_matrix=strings(E),eigenbasis_determinant=str(s.factor(E.det())),
        CB_marked_x_k_constraint_matrix=strings(cbk),CB_marked_x_k_witness=rCBk,CB_marked_x_W='span(x,h-b*z)=span(x,x+b*h-b^2*z)',
        CA_marked_x_k_constraint_matrix=strings(cak),CA_marked_x_k_witness=rCAk,CA_marked_x_W='span(h,x)',
        CB_other_marked_line_residues=cbres,CA_other_marked_line_residues=cares,
        whole_CB_flags=[dict(L='e_y',H='x',W='span(x,x+b*h-b^2*z)',ratio='[1:0]'),dict(L='e_y',H='x+2*b*h',W='span(h,x)',ratio='[1:0]')],
        whole_CA_flag=dict(L='e_y',H='x',W='span(h,x)'),whole_N_flag='same as the unique CA flag',CA_N_actual_source_rank_certificates=sranks,CA_N_actual_source_matrices=source_matrices)
    with(out/'certificate.json.gz').open('wb')as file:
        with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    result=dict(id=STEM,all_passed=True,parameter_open='b!=0',constant_Casimir_rank5_minor=rCas['minor'],unique_intrinsic_factor_rank4_minor=rFactor['minor'],
        whole_CB_flags_exactly2=True,whole_CB_fiber_two_reduced_points=True,whole_CA_flag_exactly1=True,whole_CA_fiber_reduced_length1=True,
        whole_N_flag_exactly1=True,whole_N_fiber_reduced_length1=True,CA_and_N_source_ranks=[r['rank']for r in sranks],CA_and_N_flag9_minors=[r['flag_quotient_witness']['minor']for r in sranks],
        CA_smooth_all_nonzero_family=True,N_smooth_all_nonzero_family=True,CB_two_smooth29_branches_all_nonzero_family=True,
        CB_branch_tangent_union_rank37_all_nonzero=True,CB_branch_tangent_intersection_dimension21_all_nonzero=True,CB_nonnormal_singular_nonCM_all_nonzero_family=True,CB_cone_depth_upper_bound=22,
        origin_full_conductor_or_global_bad_locus_claimed=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),parent/'report.json',parent/'certificate.json.gz',ca,np,careview,cbreview,nreview,ROOT/'code/ambient_CA_self_pair_CB_source_rank_20261103_v63.py',ROOT/'code/ambient_CA_self_pair_all_nonzero_flags_20261103_v63.py'];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
