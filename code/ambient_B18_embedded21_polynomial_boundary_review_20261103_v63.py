#!/usr/bin/env python3
"""Independent polynomial E(W,lambda) and old embedded21 boundary review.

No author module, new source jet, Groebner basis, or full-origin D computation.
Uses independent coefficient Lie derivatives and direct immutable S61 reconstruction.
"""
import argparse
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

from q_firstjet_spin4_mixed_discriminant_review_20261103_v63 import det

ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_EMBEDDED21_POLYNOMIAL_BOUNDARY_20261103_v63_REVIEW'


def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main(output):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/AMBIENT_B18_EMBEDDED21_POLYNOMIAL_BOUNDARY_20261103_v63/report.json'
    certpath=author.parent/'certificate.json.gz';ar=read(author);cert=read(certpath);assert ar['all_passed']
    recognition=ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63/report.json'
    recognition_review=ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63_REVIEW.json'
    recognition_evidence=recognition.parent/'recognition.json.gz'
    nsource=ROOT/'reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'
    oldpath=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    bundlepath=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    builderpath=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json'
    checked={}
    for path in [author,recognition,recognition_review]:
        data=read(path);assert data['all_passed']
        for key in ['input_sha256','evidence_sha256']:
            for name,want in data.get(key,{}).items():
                found=next((p for p in [ROOT/name,path.parent/name]if p.is_file()and sha(p)==want),None)
                assert found is not None,(path,name);checked[str(found.relative_to(ROOT))]=want
    assert read(recognition)['nilradical_module']=='C[[u1,...,u21]]*F'
    h,x,y,z,w,W,lam,r=s.symbols('h x y z w W lam r');v=[h,x,y,z,w]
    pairs=list(itertools.combinations(range(5),2));monpairs=list(itertools.combinations_with_replacement(range(5),2))
    mons=[s.prod(v[i]for i in mon)for mon in monpairs]
    def vector(P):return [s.expand(s.Poly(P[i,j],*v).coeff_monomial(mon))for i,j in pairs for mon in mons]
    def canonical(P):
        divergence=s.Matrix([sum(s.diff(P[i,j],v[i])for i in range(5))for j in range(5)])
        V=s.Matrix(v);return (P-(V*divergence.T-divergence*V.T)/5).applyfunc(s.expand)
    E=s.zeros(5);E[0,2]=-x*x;E[2,3]=h*x-lam*h*h;E[2,4]=h*h+z*z;E[3,4]=W*x*x;E-=E.T
    assert E==s.Matrix(cert['polynomial_family_matrix'])
    assert all(q==s.sympify(saved)for q,saved in zip(vector(E),cert['all150_polynomial_family_coefficients']))
    # Check all25 divergence and350 raw degree3 Jacobi coefficients explicitly.
    linear_mons=v;cubic_mons=[s.prod(v[i]for i in mon)for mon in itertools.combinations_with_replacement(range(5),3)]
    for j in range(5):
        div=s.Poly(s.expand(sum(s.diff(E[i,j],v[i])for i in range(5))),*v)
        assert all(div.coeff_monomial(mon)==0 for mon in linear_mons)
    for i,j,k in itertools.combinations(range(5),3):
        jac=s.expand(sum(E[i,l]*s.diff(E[j,k],v[l])+E[j,l]*s.diff(E[k,i],v[l])+E[k,l]*s.diff(E[i,j],v[l])for l in range(5)))
        polynomial=s.Poly(jac,*v);assert all(polynomial.coeff_monomial(mon)==0 for mon in cubic_mons)
    assert canonical(E)==E
    # Rebuild the actual fixed N19 Hamiltonian family and compare every basis coefficient.
    swap=s.Matrix([[0,1,0,0,0],[1,0,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]])
    def raw(f,g,H):
        P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z);return P-P.T
    def transport(P):return (swap*P.subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
    cubicvars=[h,x,z,w]
    potentials=[s.prod(cubicvars[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if 2 in mon or 3 in mon]
    fixed=[transport(canonical(raw(f,0,0)))for f in [h*h,h*x]]+[transport(canonical(raw(0,1,0)))]+[transport(canonical(raw(0,0,H)))for H in potentials]
    CN=s.Matrix.hstack(*(s.Matrix(vector(P))for P in fixed));assert CN.shape==(150,19)
    assert CN==s.Matrix(read(nsource)['fiber19_basis'])
    alpha=s.Matrix(list(map(s.sympify,cert['actual_fixed_N19_source_coefficients'])))
    assert CN*alpha==s.Matrix(vector(E))
    # Whole-family transport is direct diagonal entry scaling, with no group search.
    A=s.diag(1,1,r,1,1);V=s.Matrix(v)
    transformed=(A*E.subs(dict(zip(v,A.inv()*V)),simultaneous=True)*A.T/r).applyfunc(s.cancel)
    assert transformed==E.subs(W,W/r)
    at=E.subs({h:0,x:1,y:0,z:0,w:0})
    assert s.expand(at.extract([0,2,3,4],[0,2,3,4]).det()-W*W)==0
    # At this point Euler is e_x and the nonzero block has no e_x entries,
    # so this is also the actual projective rank, not merely raw matrix rank.
    assert at[1,:]==s.zeros(1,5)and at[:,1]==s.zeros(5,1)
    assert at.subs(W,0).rank()==2
    k0,k1,k2,k5=s.symbols('k0 k1 k2 k5')
    H=h*(x*w+z*z)+k0*x*x*z+k1*x*x*w+k2*x*z*z+k5*z**3
    oldP=raw(h*h,1,H)
    old_variables=s.symbols('x0:5')
    stored_old=s.Matrix(cert['old_free4_embedded_source_matrix']).subs(dict(zip(old_variables,v)),simultaneous=True)
    assert oldP==stored_old
    frame=s.Matrix([[0,1,0,0,0],[-1,0,0,0,0],[0,0,1,0,0],[0,-1,0,1,0],[-1,0,0,0,1]])
    inverse=frame.inv();back=(inverse*oldP.subs(dict(zip(v,frame*V)),simultaneous=True)*inverse.T).applyfunc(s.expand)
    assert back==s.Matrix(cert['old_free4_back_transported_physical_family'])
    assert back.subs({k0:1,k1:lam,k2:0,k5:s.Rational(1,3)})==E.subs(W,1)
    # Direct reconstruction through immutable S61 JSON, without normalize().
    old=read(oldpath);full_u=list(map(s.sympify,old['families']['N']['full_61_coordinates']))
    allowed={k0,k1,k2,k5};symbols=set().union(*(q.free_symbols for q in full_u))
    substitutions={q:0 for q in symbols-allowed}
    u=[s.expand(q.subs(substitutions,simultaneous=True))for q in full_u]
    assert len(u)==61 and set().union(*(q.free_symbols for q in u))<=allowed
    source=ROOT/'legacy/v02/data';slice_cert=read(source/'symplectic_slice_certificate.json');slice_data=read(source/'symplectic_slice_equations.json')
    reduced=read(source/'symplectic_slice_reduced_equations.json')
    base=raw(h*h,1,0);raw_values=vector(base)
    for index,q in zip(slice_cert['raw_free_indices'],u):raw_values[index]=q
    for row in slice_data['eliminated_raw_variables']:
        raw_values[row['raw_index']]=s.expand(sum(s.Rational(n,d)*s.prod(u[i]for i in mon)for mon,n,d in row['polynomial']))
    recon=s.zeros(5)
    for (i,j),offset in zip(pairs,range(0,150,15)):
        recon[i,j]=s.expand(sum(raw_values[offset+k]*mon for k,mon in enumerate(mons)));recon[j,i]=-recon[i,j]
    assert canonical(recon)==oldP
    assert all(raw_values[i]==0 for i in slice_cert['vertical_pivots'])
    assert all(raw_values[i]==0 for i in slice_cert['perspective_orbit_pivots_raw'])
    assert all(raw_values[i]==vector(base)[i]for i in slice_cert['constant_raw_indices'])
    assert all(not q.free_symbols for q,weight in zip(u,reduced['parameter_weights'])if weight==1)
    recognize=read(recognition_evidence);at_old={k0:1,k1:0,k2:0,k5:s.Rational(1,3)}
    assert [q.subs(at_old)for q in u]==list(map(s.sympify,recognize['exact_normalized61']))
    bundle,builder=read(bundlepath),read(builderpath)
    q=[u[i]for i in bundle['q_indices_zero_based']]
    K=[q[i]for i in builder['common_q_selected_rows']];remaining=[q[i]for i in builder['remaining_q_rows']]
    expected=[0,-k1,0,0,k0,2*k2,0,3*k5,0,0]
    assert K==list(map(s.sympify,cert['old_target10_embedded_support_coordinates']))==expected
    assert all(value==0 for value in remaining)
    # Coefficient-index Lie action: E_ab P+P E_ab^T−x_b*d_aP.
    P0=E.subs({W:0,lam:0});base_rows={}
    for ii,jj in pairs:
        base_rows[ii,jj]={mon:Q(str(c))for mon,c in s.Poly(P0[ii,jj],*v).terms()if c}
    def rowpoly(i,j):
        if i<j:return base_rows[i,j]
        if i>j:return {mon:-c for mon,c in base_rows[j,i].items()}
        return {}
    columns=[]
    for generator in cert['origin_GL19_generator_indices']:
        aa,bb=divmod(generator,5);column=[]
        for ii,jj in pairs:
            poly={}
            for sourcepoly in ([rowpoly(bb,jj)]if ii==aa else[])+([rowpoly(ii,bb)]if jj==aa else[]):
                for mon,c in sourcepoly.items():poly[mon]=poly.get(mon,Q(0))+c
            for mon,c in rowpoly(ii,jj).items():
                if mon[aa]:
                    target=list(mon);target[aa]-=1;target[bb]+=1;target=tuple(target)
                    poly[target]=poly.get(target,Q(0))-mon[aa]*c
            column += [poly.get(tuple(mon.count(k)for k in range(5)),Q(0))for mon in monpairs]
        columns.append(column)
    columns += [list(map(lambda q:Q(str(q)),vector(E.diff(W)))),list(map(lambda q:Q(str(q)),vector(E.diff(lam))))]
    Jac=list(map(list,zip(*columns)))
    assert Jac==[[Q(x)for x in row]for row in cert['full150_by21_origin_Jacobian']]
    selected=cert['origin21_minor_rows'];minor=det([Jac[i]for i in selected]);assert minor==32==Q(cert['origin21_minor_determinant'])
    evidence={
        'all25_divergence_all350_raw_Jacobi_coefficients_zero':True,
        'full_actual_N19_basis_and_E_membership_all150':True,
        'projective_rank_witness_Euler_excluded':True,
        'whole_W_lambda_group_identity':True,
        'independent_direct_S61_JSON_reconstruction_not_author_normalize':True,
        'all61_free4_family_same_physical_class_and_perspective_zero':True,
        'old_target_free4_values':[str(x)for x in K],
        'all35_remaining_q_zero':True,
        'own_coefficient_Lie_all150x21_entries_match':True,
        'origin21_minor_rows':selected,'origin21_minor':str(minor),
        'connected_group_and_closure_proof_in_review_note':True,
        'dominant_rational_map_not_birationality':True,
        'formal_regular_image_branch_only':True,
    }
    (output/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,author,certpath,recognition,recognition_review,recognition_evidence,nsource,oldpath,bundlepath,builderpath,
           source/'symplectic_slice_certificate.json',source/'symplectic_slice_equations.json',source/'symplectic_slice_reduced_equations.json',ROOT/'sources/v02.zip',
           ROOT/'code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py',
           ROOT/'research/notes/AMBIENT_B18_EMBEDDED21_POLYNOMIAL_BOUNDARY_20261103_v63.md']
    report={
        'id':ID,'all_passed':True,
        'all_polynomial_physical_divergence_Jacobi_N19_group_and_rank_identities_accepted':True,
        'old_free4_embedded_prime_family_full61_literal_binding_accepted':True,
        'connected_global_Ass_action_fixes_embedded21_support':True,
        'global_embedded21_support_contains_entire_E_W_lambda_and_P0':True,
        'dominant_rational21_parameter_map_accepted':True,
        'one_regular21_dimensional_formal_branch_and_transverse2_image_at_P0':True,
        'whole_germ_smooth_irreducible_rationality_birationality_or_full_D_type_claim':False,
        'new_Poisson_component_or_origin_Ass_exhaustion_claim':False,
        'no_author_code_new_source_jet_GB_or_CAS_rerun':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},
        'checked_dependency_sha256':checked,
        'evidence_sha256':{'independent-certificate.json':sha(output/'independent-certificate.json')},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'seconds':time.monotonic()-start,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'origin_minor':str(minor),'seconds':report['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
