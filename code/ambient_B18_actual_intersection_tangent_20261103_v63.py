#!/usr/bin/env python3
"""Actual B18 cotangent embeddings from even kernel-square and I/mI classes."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ACTUAL_INTERSECTION_TANGENT_20261103_v63'
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for q in iter(lambda:f.read(1<<20),b''):h.update(q)
    return h.hexdigest()
def same(A,B):assert A.shape==B.shape and(A-B).applyfunc(s.expand)==s.zeros(*A.shape)
def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120);out.mkdir(parents=True,exist_ok=False)
    h,x,y,z,w,p=s.symbols('h x y z w p');X=s.Matrix([h,x,y,z,w]);pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[u]*X[v])for i,j in pairs for u,v in mons])
    def canonical(P):
        d=s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)])
        return(P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(f,g,H,cb=False):
        P=s.zeros(5);P[1,2]=f;P[3,4]=g if cb else g*h*h
        P[2,3]=(-s.diff(f,x)*z/2 if cb else 0)-s.diff(H,w)
        P[2,4]=(-s.diff(f,x)*w/2 if cb else 0)+s.diff(H,z)
        return P-P.T
    swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
    change=lambda P:(swap*P.subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
    cubics=[s.prod((h,x,z,w)[i]for i in m)for m in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in m)]
    ey=s.eye(5)[:,2];D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T
    def lie(i,j,P):
        E=s.zeros(5);E[i,j]=1;V=E*X
        return(E*P+P*E.T-sum((V[k]*P.diff(X[k])for k in range(5)),s.zeros(5))).applyfunc(s.expand)
    O=s.Matrix.hstack(*(vec(lie(i,j,P0))for i,j in itertools.product(range(5),repeat=2)))
    reduction_path=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json';reduction=read(reduction_path)
    first_path=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63/certificate.json.gz';first=read(first_path)
    quad_path=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/certificate.json.gz';quad=read(quad_path)
    cubic_path=ROOT/'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63/certificate.json.gz';cubic=read(cubic_path)
    model_path=ROOT/'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63/report.json';modelreview=ROOT/'reports/AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW/report.json'
    assert read(modelreview)['actual_CA_CB_separate_full_formal_models_accepted']and read(model_path)['cone_embedding_dimension']==39
    Rrows=reduction['common28_linear_target_rows']+[reduction['common29th_target_row']];assert len(Rrows)==29
    records={};planes=[];coordinates=None
    for name in['CA','CB']:
        cb=name=='CB';r=reduction['families'][name]
        basis=[change(canonical(raw(f,0,0,cb)))for f in[h*h,h*x,x*x]]
        basis+=[change(canonical(raw(0,h*h if cb else 1,0,cb)))];basis+=[change(canonical(raw(0,0,H,cb)))for H in cubics]
        C=s.Matrix.hstack(*map(vec,basis));oi=reduction['orbit_GL_generator_indices'];extras=r['source_fiber_basis_indices']
        J=O[:,oi].row_join(C[:,extras]);same(J,s.Matrix(r['linear_J']))
        # Pure kernel curve is reconstructed by a FULL physical finite matrix
        # substitution. Its p² coefficient is even because R(p)=O(p²),t=O(p).
        M=s.Matrix(r['source_kernel_Lie_matrix']);assert M**3==s.zeros(5)
        co=s.Matrix(r['source_velocity']);alpha=s.zeros(20,1);alpha[0]=-1;alpha[7]=-1;alpha[11]=1;alpha[16]=s.Rational(1,3);alpha-=p*co
        f=alpha[0]*h*h+alpha[1]*h*x+alpha[2]*x*x;H=sum(q*m for q,m in zip(alpha[4:,0],cubics))
        src=change(canonical(raw(f,-alpha[3]*f if cb else alpha[3],H,cb)))
        G=s.eye(5)+p*M+p*p*M*M/2;Gi=s.eye(5)-p*M+p*p*M*M/2;same(G*Gi,s.eye(5))
        U=vec((G*src.subs(dict(zip(X,Gi*X)),simultaneous=True)*G.T).applyfunc(s.expand))
        same(U.subs(p,0),vec(P0));same(U.diff(p).subs(p,0),s.zeros(150,1))
        W=U.applyfunc(lambda q:s.Poly(q,p).nth(2));same(W,s.Matrix(r['raw_second_coefficient']))
        even_source=J.row_join(W);M29=even_source[Rrows,:];assert M29.det()!=0
        E=even_source*M29.inv();same(E[Rrows,:],s.eye(29))
        # Actual I/mI basis, not arbitrary source Taylor monomials.
        A=s.Matrix(first['CA_CB_full_firstjet_data'][name]['actual']);q=quad[name];rr=q['actual_seven_odd_basis_rows'];cc=q['selected_actual_base_columns']
        assert rr==[18,19,28,29,60,67,72]
        L=A[:,cc]*A[rr,cc].inv();same(L*A[rr,:],A);same(L,s.Matrix(q['all_actual_linear_combinations']))
        Qraw=s.Matrix(q['quadratic_coefficient_matrix'])[:,:2];qr=q['independent_quadratic_rows'];assert qr==[46,95]
        Q=Qraw*Qraw[qr,:].inv();same(Q*Qraw[qr,:],Qraw)
        actual_cubic=s.Matrix(cubic[name]['actual_all150_odd_cubic_coefficients']);remaining=actual_cubic-L*actual_cubic[rr,:]
        remaining-=Q*remaining[qr,:];same(remaining,s.Matrix(cubic[name]['actual_zero_linear_and_quad_cubic_classes']))
        assert[i for i,q in enumerate(remaining)if q]==[63]
        kappa=remaining[63];expected=s.Rational(12,125)if name=='CA'else-s.Rational(4,243);assert kappa==expected
        C3=remaining/kappa;B=L.row_join(Q).row_join(C3);same(B[Rrows,:],s.zeros(29,10))
        T=E.row_join(B);oddrows=rr+qr+[63];selected=Rrows+oddrows;assert len(set(selected))==39
        if coordinates is None:coordinates=selected
        assert coordinates==selected
        selected_matrix=T[selected,:];det=s.factor(selected_matrix.det());assert det!=0
        normalized=T*selected_matrix.inv();same(normalized[selected,:],s.eye(39));planes.append(normalized)
        records[name]=dict(full150_by28_physical_even_linear_J=[[str(q)for q in row]for row in J.tolist()],physical_pure_kernel_second_W=list(map(str,W)),even29_minor=str(s.factor(M29.det())),
            full150_by29_even_cotangent_embedding=[[str(q)for q in row]for row in E.tolist()],actual_odd_linear_basis_rows=rr,actual_odd_quadratic_basis_rows=qr,actual_odd_cubic_basis_row=63,actual_odd_cubic_normalization=str(kappa),
            full150_by10_actual_odd_cotangent_embedding=[[str(q)for q in row]for row in B.tolist()],full150_by39_actual_tangent=[[str(q)for q in row]for row in T.tolist()],selected39_minor=str(det))
    same(*planes);assert planes[0].rank()==39
    assert s.Matrix.hstack(s.Matrix(records['CA']['full150_by39_actual_tangent']),s.Matrix(records['CB']['full150_by39_actual_tangent'])).rank()==39
    # The radial line lies in this actual tangent; q106=1 is a scalar chart.
    assert vec(P0)[106]==1 and planes[0].row_join(vec(P0)).rank()==39
    paths=[ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
           ROOT/'reports/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_CA_CB_MINIMAL_LENGTHS_AND_R_CM_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
           ROOT/'reports/AMBIENT_R19_BOUNDARY_SPECTRAL_MODEL_20261103_v63_REVIEW.json']
    assert all(read(p)['all_passed']for p in paths)
    assert read(paths[4])['Hilbert_Samuel_multiplicity']==2 and read(paths[5])['generic_N_length1_by_common_reduced_Z_base_and_flat_functionfield_A10_intersection_A9']
    evidence=dict(common29_even_target_rows=Rrows,common10_odd_target_rows=oddrows,common39_physical_coordinate_rows=coordinates,
        common150_by39_normalized_tangent=[[str(q)for q in row]for row in planes[0].tolist()],families=records)
    with(out/'certificate.json.gz').open('wb')as f:
        with gzip.GzipFile(fileobj=f,mode='wb',filename='',mtime=0)as g:g.write((json.dumps(evidence,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/f'research/notes/{STEM}.md',reduction_path,first_path,quad_path,cubic_path,model_path,modelreview,*paths]
    result=dict(id=STEM,all_passed=True,point='generic B18 P0=e_y wedge (x²,0,0,hx,h²+z²)',
        independent_physical_raw_CA_CB_bases_all25_GL_derivatives_and_pure_kernel_curve_second_coefficient=True,
        even_R29_cotangent_from_J28_and_true_kernel_square_not_arbitrary_source_jet_span=True,
        old_common_R29_odd_minimal_number10_from_intrinsic_edim39_not_new_base_alignment=True,
        actual_I_mod_mI_basis_7linear_2quadratic_1cubic_and_all150_coefficients=True,
        full_actual_CA_and_CB_tangent_dimensions=[39,39],union_tangent_dimension=39,
        actual_full_CA_and_CB_tangent_planes_equal=True,actual_D_tangent_dimension=39,
        common39_physical_coordinate_rows=coordinates,
        actual_D_cone_dimension=28,actual_D_cone_embedding_dimension=39,actual_D_Hilbert_Samuel_multiplicity=2,actual_D_cone_depth_upper_bound=20,
        actual_D_projective_dimension=27,actual_D_projective_embedding_dimension=38,actual_D_projective_depth_upper_bound=19,
        multiplicity_from_only_topdim_N_generic_length1_times_N_multiplicity2=True,
        depth_bound_from_incomparable_minimal_R20_component=True,
        full_D_origin_ring_exact_depth_or_nilradical_linear_visibility_claimed=False,
        T_N_dimension37_not_used_to_infer_reduced_union_tangent_or_nilvisibility=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},evidence_sha256={'certificate.json.gz':sha(out/'certificate.json.gz')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:result[k]for k in['all_passed','full_actual_CA_and_CB_tangent_dimensions','actual_full_CA_and_CB_tangent_planes_equal','actual_D_cone_embedding_dimension','actual_D_Hilbert_Samuel_multiplicity','seconds']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
