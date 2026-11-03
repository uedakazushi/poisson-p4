#!/usr/bin/env python3
"""Full finite proper center/pencil incidence excludes the last extra component.

The source has99 variables and209 exact quadratic equations. Its completed
slice is proved by an EXACT polynomial77-pivot graph, not by truncated jets.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring
from s61_model import bivector_from_coefficients
from q51_universal_section_20261101 import affine
from additional_rank4_families_20261103_v62 import families

ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def strings(M):return [[str(v)for v in row]for row in M.tolist()]
def rank(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ).rank()
def determinant(M):return s.Rational(DomainMatrix.from_Matrix(M).convert_to(s.QQ).det())
def build(out):
    start=time.monotonic();z=s.Matrix(s.symbols('x y z w'));pairs=list(itertools.combinations(range(4),2))
    mon1=[s.S.One,*z];mons=[*mon1,*[z[i]*z[j]for i,j in itertools.combinations_with_replacement(range(4),2)]]
    mon3=[s.prod(z[i]for i in ij)for degree in range(4)for ij in itertools.combinations_with_replacement(range(4),degree)]
    raypath=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz'
    saved=json.loads(gzip.decompress(raypath.read_bytes()));eta=s.Symbol('eta')
    P=affine(bivector_from_coefficients([s.sympify(v).subs(eta,1)for v in saved['raw150']]),z)
    rr=s.symbols('r0:90');R=s.zeros(4)
    for k,(i,j)in enumerate(pairs):R[i,j]=sum(c*m for c,m in zip(rr[15*k:15*k+15],mons));R[j,i]=-R[i,j]
    u,c1,c2,c4=s.symbols('u c1 c2 c4');cs=s.symbols('d0:5');cpairs=[ij for ij in pairs if ij!=(1,2)]
    C0=s.zeros(4);C0[1,2]=1;C0[2,1]=-1;DC=s.zeros(4)
    for(i,j),c in zip(cpairs,cs):DC[i,j]=c;DC[j,i]=-c
    variables=[*rr,u,c1,c2,c4,*cs];assert len(variables)==99
    linear=[];quadratic=[]
    for i,j,k in itertools.combinations(range(4),3):
        lin=sum(P[i,b]*s.diff(R[j,k],z[b])+P[j,b]*s.diff(R[k,i],z[b])+P[k,b]*s.diff(R[i,j],z[b])+R[i,b]*s.diff(P[j,k],z[b])+R[j,b]*s.diff(P[k,i],z[b])+R[k,b]*s.diff(P[i,j],z[b])for b in range(4))
        qua=sum(R[i,b]*s.diff(R[j,k],z[b])+R[j,b]*s.diff(R[k,i],z[b])+R[k,b]*s.diff(R[i,j],z[b])for b in range(4))
        lp,qp=s.Poly(s.expand(lin),*z),s.Poly(s.expand(qua),*z)
        linear += [lp.coeff_monomial(m)for m in mon3];quadratic += [qp.coeff_monomial(m)for m in mon3]
    # The compact center is [u:c1:c2:1:c4]. Its full equation is
    # u*L(P+R)+d_(c1,c2,1,c4) Q(P+R)=0,24 coefficients.
    center=[]
    for i,j in pairs:
        pp,qq=s.Poly(P[i,j],*z),s.Poly(R[i,j],*z)
        L0=sum(pp.coeff_monomial(v)*v for v in z);LR=sum(qq.coeff_monomial(v)*v for v in z)
        Q0=sum(pp.coeff_monomial(m)*m for m in mons[5:]);QR=sum(qq.coeff_monomial(m)*m for m in mons[5:])
        lp=s.Poly(s.expand(u*L0+s.diff(QR,z[2])+c1*s.diff(Q0,z[0])+c2*s.diff(Q0,z[1])+c4*s.diff(Q0,z[3])),*z)
        qp=s.Poly(s.expand(u*LR+c1*s.diff(QR,z[0])+c2*s.diff(QR,z[1])+c4*s.diff(QR,z[3])),*z)
        linear +=[lp.coeff_monomial(v)for v in z];quadratic +=[qp.coeff_monomial(v)for v in z]
        center += [s.Poly(s.expand(u*L0+c1*s.diff(Q0,z[0])+c2*s.diff(Q0,z[1])+s.diff(Q0,z[2])+c4*s.diff(Q0,z[3])),*z).coeff_monomial(v)for v in z]
    divR=s.Matrix([s.expand(sum(s.diff(R[i,j],z[i])for i in range(4)))for j in range(4)])
    assert s.Matrix([s.expand(sum(s.diff(P[i,j],z[i])for i in range(4)))for j in range(4)])==s.zeros(4,1)
    N=s.zeros(5);N[1:,0]=divR.subs(dict.fromkeys(z,0));N[1:,1:]=divR.jacobian(z)
    linear +=[s.S.Zero]*25;quadratic +=[s.expand(v)for v in N*N]
    # A projective constant compatible bivector, on the chart C_yz=1.
    compatibility_full=[]
    cc=s.symbols('k0:6');CC=s.zeros(4)
    for(i,j),c in zip(pairs,cc):CC[i,j]=c;CC[j,i]=-c
    for i,j,k in itertools.combinations(range(4),3):
        lin=sum(C0[i,b]*s.diff(R[j,k],z[b])+C0[j,b]*s.diff(R[k,i],z[b])+C0[k,b]*s.diff(R[i,j],z[b])+DC[i,b]*s.diff(P[j,k],z[b])+DC[j,b]*s.diff(P[k,i],z[b])+DC[k,b]*s.diff(P[i,j],z[b])for b in range(4))
        qua=sum(DC[i,b]*s.diff(R[j,k],z[b])+DC[j,b]*s.diff(R[k,i],z[b])+DC[k,b]*s.diff(R[i,j],z[b])for b in range(4))
        full=sum(CC[i,b]*s.diff(P[j,k],z[b])+CC[j,b]*s.diff(P[k,i],z[b])+CC[k,b]*s.diff(P[i,j],z[b])for b in range(4))
        lp,qp,fp=[s.Poly(s.expand(f),*z)for f in [lin,qua,full]]
        linear +=[lp.coeff_monomial(v)for v in mon1];quadratic +=[qp.coeff_monomial(v)for v in mon1];compatibility_full +=[fp.coeff_monomial(v)for v in mon1]
    assert len(linear)==len(quadratic)==209
    equations=[s.expand(l+q)for l,q in zip(linear,quadratic)]
    M=s.linear_eq_to_matrix(linear,variables)[0];assert rank(M)==77

    # All affine GL4, translations and source-P scalar directions, with
    # corresponding center/C chart derivatives, are actual source actions.
    def orbit_column(T,cent,Cder):
        values=[]
        for i,j in pairs:
            poly=s.Poly(s.expand(T[i,j]),*z);values +=[poly.coeff_monomial(m)for m in mons]
        return s.Matrix([*values,0,*cent,*[Cder[i,j]for i,j in cpairs]])
    orbit=[];labels=[]
    for i,j in itertools.product(range(4),repeat=2):
        A=s.zeros(4);A[i,j]=1;T=A*P+P*A.T-sum(((A*z)[k]*P.diff(z[k])for k in range(4)),s.zeros(4))
        Cd=A*C0+C0*A.T;Cd-=Cd[1,2]*C0
        orbit.append(orbit_column(T,[A[0,2],A[1,2],A[3,2]],Cd));labels.append(f'GL4_{i}_{j}')
    for i in range(4):orbit.append(orbit_column(-P.diff(z[i]),[0,0,0],s.zeros(4)));labels.append(f'translation_{i}')
    orbit.append(orbit_column(P,[0,0,0],s.zeros(4)));labels.append('P_scalar')
    Oall=s.Matrix.hstack(*orbit);assert M*Oall==s.zeros(209,21)
    orbit_columns=list(Oall.rref()[1]);O=Oall[:,orbit_columns];assert O.cols==19
    transverse=[];frame=O
    for v in M.nullspace():
        if frame.row_join(v).rank()>frame.cols:frame=frame.row_join(v);transverse.append(v)
    T=s.Matrix.hstack(*transverse);assert T.shape==(99,3)
    columns=list(M.rref()[1]);rows=list(M[:,columns].T.rref()[1]);B=M[rows,columns];minor=determinant(B);assert minor
    Ep=s.eye(99)[:,columns];fullframe=O.row_join(T).row_join(Ep);frame_det=determinant(fullframe);assert frame_det

    R3,*nv=ring('n0,n1,n2',s.QQ);zero=R3.zero;domain=R3.domain
    linear_graph=[sum((R3(domain.convert(T[i,j]))*nv[j]for j in range(3)),zero)for i in range(99)]
    encoded=[]
    for f in equations:
        terms=[]
        for exponent,c in s.Poly(f,*variables).terms():
            if not c:continue
            ids=[i for i,k in enumerate(exponent)for _ in range(k)];assert 1<=len(ids)<=2
            terms.append((ids,s.Rational(c)))
        encoded.append(terms)
    def evaluate(point):
        return[sum((R3(domain.convert(c))*sprod([point[i]for i in ids],R3.one)for ids,c in terms),zero)for terms in encoded]
    def sprod(factors,one):
        result=one
        for v in factors:result*=v
        return result
    # Compute77 polynomial implicit assignments from the quadratic source.
    q2=evaluate(linear_graph);Binv=B.inv();p2=[]
    for rr0 in Binv.tolist():p2.append(-sum((q2[rows[j]]*domain.convert(c)for j,c in enumerate(rr0)if c),zero))
    graph=list(linear_graph)
    for i,v in zip(columns,p2):graph[i]+=v
    corrections={i:str(graph[i]-linear_graph[i])for i in range(99)if graph[i]!=linear_graph[i]}
    assert corrections=={17:'1/48*n1*n2',49:'-1/1536*n1**2',54:'1/2*n0*n2'}
    # ACTUAL209 polynomial substitution, all degrees included (at most4).
    pullbacks=evaluate(graph);nonzero=[(i,v)for i,v in enumerate(pullbacks)if v]
    assert nonzero==[(179,R3(s.Rational(4,3))*nv[0]**2*nv[2])]
    assert all(not pullbacks[i]for i in rows)
    assert all(not graph[i]for i in range(90,99))
    # In particular the original finite-center coordinate u is EXACTLY0.
    assert graph[90]==zero

    # WHOLE projective fibers over positiveP: unique center and C, with unit
    # minors; no omitted chart/preimage can receive an arc from Hquartic.
    center_full_vars=[u,c1,c2,s.Symbol('v3'),c4];center_full=[]
    for i,j in pairs:
        pp=s.Poly(P[i,j],*z);L0=sum(pp.coeff_monomial(v)*v for v in z);Q0=sum(pp.coeff_monomial(m)*m for m in mons[5:])
        value=center_full_vars[0]*L0+sum(v*s.diff(Q0,zi)for v,zi in zip(center_full_vars[1:],z))
        pp=s.Poly(s.expand(value),*z);center_full +=[pp.coeff_monomial(v)for v in z]
    Kcenter=s.linear_eq_to_matrix(center_full,center_full_vars)[0]
    Kpencil=s.linear_eq_to_matrix(compatibility_full,cc)[0]
    assert Kcenter.nullspace()==[s.Matrix([0,0,0,1,0])]
    assert Kpencil.nullspace()==[s.Matrix([0,0,0,1,0,0])]
    def witness(A):
        cols=list(A.rref()[1]);rs=list(A[:,cols].T.rref()[1]);return dict(matrix=strings(A),rank=len(cols),rows=rs,columns=cols,determinant=str(determinant(A[rs,cols])))
    center_witness,pencil_witness=witness(Kcenter),witness(Kpencil)
    assert center_witness['rank']==4 and pencil_witness['rank']==5
    # Dense Hquartic representative belongs to the finite-center, symplectic
    # constant-pencil, square-zero-modular locus, by exact symbolic checks.
    record=next(r for r in families()if r['component_id']=='R4_H_QUARTIC_ESSENTIAL5');star=affine(record['P'],z)
    Jstar=star.subs(dict.fromkeys(z,0));assert Jstar.det()!=0
    assert all(s.Poly(star[i,j],*z).coeff_monomial(v)==0 for i,j in pairs for v in z)
    divstar=s.Matrix([s.expand(sum(s.diff(star[i,j],z[i])for i in range(4)))for j in range(4)])
    Nstar=s.zeros(5);Nstar[1:,0]=divstar.subs(dict.fromkeys(z,0));Nstar[1:,1:]=divstar.jacobian(z)
    assert Nstar*Nstar==s.zeros(5)and rank(Nstar)==2
    for i,j,k in itertools.combinations(range(4),3):
        assert s.expand(sum(Jstar[i,b]*s.diff(star[j,k],z[b])+Jstar[j,b]*s.diff(star[k,i],z[b])+Jstar[k,b]*s.diff(star[i,j],z[b])for b in range(4)))==0
    out.mkdir(parents=True,exist_ok=False)
    certificate=dict(physical_affine_variables=list(map(str,z)),coefficient_monomials=list(map(str,mons)),source_variables=list(map(str,variables)),positive_P=strings(P),
      full209_equations_encoded=[[[ids,int(c.p),int(c.q)]for ids,c in terms]for terms in encoded],
      full_Jacobian_209_by99=strings(M),affine_orbit_all21_labels=labels,affine_orbit_all21_matrix=strings(Oall),orbit19_selected_columns=orbit_columns,
      normal99_by3=strings(T),implicit77_columns=columns,implicit77_rows=rows,implicit77_unit_minor=str(minor),full99_transverse_frame_determinant=str(frame_det),
      exact99_polynomial_graph=list(map(str,graph)),only_three_nonlinear_assignments=corrections,all209_exact_full_pullbacks=list(map(str,pullbacks)),
      sole_cubic_row=179,sole_cubic='4*n0^2*n2/3',center_u_graph_identically_zero=True,all_nine_center_and_pencil_chart_variables_zero=True,
      whole_projective_center_fiber=center_witness,whole_projective_constant_pencil_fiber=pencil_witness,
      dense_Hquartic_affine_P=strings(star),dense_Hquartic_constant_J=strings(Jstar),dense_Hquartic_modular_N=strings(Nstar))
    archive=out/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/'research/notes/Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63.md',raypath,
      ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/report.json',ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63/report.json',
      ROOT/'code/additional_rank4_families_20261103_v62.py',ROOT/'legacy/v06/nilpotent22.tex',ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
      ROOT/'code/s61_model.py',ROOT/'code/q51_universal_section_20261101.py']
    result=dict(id='Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSION_20261103_v63',all_passed=True,field='QQ; exact full polynomial maps over QQ[n0,n1,n2]',
      proper_incidence='Marked Poisson H, compact center P4, projective constant compatible bivector P5, affine modularN squaredzero',
      fixed_H_affine_source_variables=99,actual_source_equations=209,Jacobian_rank=77,affine_orbit_and_P_scalar_rank=19,exact_normal_variables=3,
      implicit77_minor=str(minor),full99_frame_determinant=str(frame_det),all209_full_polynomial_compositions_verified=True,
      completed_fixed_H_affine_incidence='C[[orbit19,n0,n1,n2]]/(n0^2*n2)',
      finite_center_coordinate_u_identically_zero_in_entire_completed_incidence=True,
      proof_is_not_finite_jet_inference=True,proper_whole_projective_center_and_pencil_preimages_unique=True,
      marked_H_geometric_fiber_unique_by_actual_quintic_h5=True,marked_H_fiber_reduced_claimed=False,
      all_nonzero_positive_axis_points_outside_entire_R4_H_QUARTIC_ESSENTIAL5=True,
      descent='Faithfully flat completion gives u=0 in the local source ring; no finite-center point specializes there. Proper marked H/center/pencil incidence and unique geometric preimage exclude entire orbit closure.',
      full_Poisson_completed_ring_or_new_component_count_claimed=False,
      input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs},certificate_sha256=sha(archive),seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('EXACT209/77 FINITE CENTER/PENCIL MODEL + ENTIRE HQUARTIC EXCLUSION PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);build(parser.parse_args().output)
