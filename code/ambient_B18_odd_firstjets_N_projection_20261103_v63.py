#!/usr/bin/env python3
"""Actual odd firstjets in common target coordinates and N pinch projection.

Only coefficients of the full source needed for the exact firstjets are
truncated. The rank2 algebra then proves the whole projected N formal type.
"""
import argparse,gzip,hashlib,itertools,json,time
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector,action
ROOT=Path(__file__).resolve().parents[1];h,x,y,z,w=X
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(M):return [[str(v)for v in row]for row in M.tolist()]
def coefficient(matrix,p,degree,zero):return s.Matrix([s.Poly(v,p).nth(degree).subs(zero)for v in matrix])
def oddjet(U,J,p,K,linearrows,branchrow):
    zero=dict.fromkeys([p,*K],0);inverse=J[linearrows,:].inv()
    assert U.diff(p).subs(zero)==s.zeros(150,1)
    W=coefficient(U,p,2,dict.fromkeys(K,0));Z=coefficient(U,p,3,dict.fromkeys(K,0))
    mixed=s.zeros(150,19).row_join(U.diff(p).jacobian(K).subs(zero));assert mixed.cols==J.cols
    rem=lambda V:(V-J*inverse*V[linearrows,:]).applyfunc(s.expand)
    WW=rem(W);s20=-inverse*W[linearrows,:];ZZ=rem(Z+mixed*s20);EE=rem(mixed*inverse)
    weight=WW[branchrow];assert weight!=0
    normal=lambda V:(V-WW/weight*V[branchrow,:]).applyfunc(s.expand)
    tangent=normal(EE);normalcol=normal(ZZ)/weight
    L=J[branchrow,:]*inverse
    actual=tangent.row_join(normalcol)
    # eta=r_branch−L*Y; the actualraw coordinate basis is (Y,r_branch).
    actual[:,:J.cols]=tangent-normalcol*L
    discriminator=s.Matrix([list(-L/weight)+[1/weight]])
    return dict(W=W,correctedW=WW,rawZ=Z,correctedZ=ZZ,mixed=mixed,tangent=tangent,normal=normalcol,actual=actual,discriminant=discriminator,L=L,weight=weight)
def main(out):
    started=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json';data=json.loads(source.read_text())
    dependencies=[ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json']
    for path in dependencies:
        check=json.loads(path.read_text());assert check['all_passed'] is True
        for name,digest in check['input_sha256'].items():assert sha(ROOT/name)==digest,name
    p=s.Symbol('p');kk=s.symbols('k0:9');Y=data['common28_linear_target_rows'];extra=data['common29th_target_row'];coordinate126=data['common_target126_raw_coordinate_indices']
    ey=s.Matrix([0,0,1,0,0]);D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D.T-D*ey.T;v0=vector(P0)
    allO=s.Matrix.hstack(*[vector(action(s.eye(5)[:,i]*s.eye(5)[j,:],P0))for i,j in itertools.product(range(5),repeat=2)])
    O=allO[:,data['orbit_GL_generator_indices']]
    records={};disc=[]
    for name,entry in data['families'].items():
        U=s.Matrix(list(map(s.sympify,entry['core_full150_coordinates'])));zero=dict.fromkeys([p,*kk],0)
        assert U.subs(zero)==v0
        J=O.row_join(U.jacobian(kk).subs(zero));assert J==s.Matrix(entry['linear_J'])
        result=oddjet(U,J,p,list(kk),Y,extra)
        A=result['actual'];assert result['tangent'].rank()==6 and A.rank()==7
        assert A[coordinate126[29:],:].rank()==7
        assert A[Y+[extra],:]==s.zeros(29,29)
        records[name]={key:mat(value)if isinstance(value,s.MatrixBase)else str(value)for key,value in result.items()}
        disc.append(result['discriminant'])
    assert s.Matrix.vstack(*disc).rank()==2
    # Exact N source is the CA source's k1=0 subbundle, then the sameH19.
    ca=data['families']['CA'];assert s.sympify(ca['source_velocity'][2])==0
    U=s.Matrix([s.sympify(v).subs(kk[1],0)for v in ca['core_full150_coordinates']]);Kn=[k for i,k in enumerate(kk)if i!=1]
    Jn=O.row_join(U.jacobian(Kn).subs(dict.fromkeys([p,*Kn],0)))
    augmented=s.Matrix(ca['entire_source_tangent'])[:,[i for i in range(29)if i!=20]]
    assert augmented.rank()==28
    # Check N's fixed-fiber constraint before the moving kernel flag action.
    M=s.Matrix(ca['source_kernel_Lie_matrix']);g=s.eye(5)+p*M+p*p*M*M/2;gi=s.eye(5)-p*M+p*p*M*M/2
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    PU=s.zeros(5)
    for offset,(i,j)in enumerate(pairs):
        PU[i,j]=sum(U[15*offset+k]*X[a]*X[b]for k,(a,b)in enumerate(mons));PU[j,i]=-PU[i,j]
    fixed=(gi*PU.subs(dict(zip(X,g*s.Matrix(X))),simultaneous=True)*gi.T).applyfunc(s.expand)
    cubics=list(map(s.sympify,data['cubic_monomials']));sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    change=lambda Q:(sw*Q.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    basis=[change(canonical(raw(f,0,0)))for f in[h*h,h*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H)))for H in cubics]
    C=s.Matrix.hstack(*map(vector,basis));cr=list(C.T.rref()[1]);left=C[cr,:].inv()
    assert C*left*vector(fixed)[cr,:]==vector(fixed)
    nr=[i for i in Y if i!=27];assert len(nr)==27
    N=oddjet(U,Jn,p,Kn,nr,27)
    assert N['weight']==-3
    assert set(nr+[27])==set(Y)
    # Reorder derivative coordinates to the actual common first28 rowsY.
    oldorder=nr+[27];permutation=[oldorder.index(i)for i in Y]
    bN=N['actual'][61,permutation];nuN=N['discriminant'][:,permutation]
    expected={16:-1,65:-6,105:4,112:1,117:-4}
    assert list(bN)==[expected.get(row,0)for row in Y]
    assert bN[0,Y.index(27)]==0
    assert s.Matrix.vstack(nuN,bN).rank()==2
    scalar=s.zeros(1,28);scalar[0,Y.index(106)]=1
    assert s.Matrix.vstack(nuN,bN,scalar).rank()==3
    certificate=dict(common_Y28=Y,common29th_row=extra,CA_CB_full_firstjet_data=records,
        CA_CB_discriminant_differentials=mat(s.Matrix.vstack(*disc)),
        N_source_exact_k1zero=True,N_augmented_source28_rank=True,N_own_linear27_rows=nr,N_own_branch_row=27,
        N_full_firstjet_data={key:mat(value)if isinstance(value,s.MatrixBase)else str(value)for key,value in N.items()},
        N_odd_extra61_differential_in_actual_Y_order=mat(bN),N_discriminant_differential_in_actual_Y_order=mat(nuN),
        N_regular_scalar_coordinate=106,N_independent_discriminant_odd_and_scalar_rank=3)
    (out/'certificate.json.gz').write_bytes(gzip.compress((json.dumps(certificate,separators=(',',':'))+'\n').encode(),mtime=0))
    note=ROOT/'research/notes/AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63.md'
    paths=[Path(__file__),note,source,*dependencies,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py']
    result=dict(id='AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63',all_passed=True,
        all_finite_source150_maps_and28_derivatives_bound_to_physical_Lie=True,
        exact_mixed_p_source_and_cubic_implicit_corrections=True,
        actual_odd_coefficient_linear_ranks_CA_CB={'tangent28':6,'whole29':7},
        branch_discriminant_differentials_in_common_R29_independent=True,
        N_source_restriction_k1zero_verified_all_orders_and_etale28=True,
        N_own_R28_is_common_first28_actual_target_coordinates=True,
        N_extra29th_projection_odd_coefficient_nonzero_regular_parameter=True,
        N_odd_linear_form='-Y16-6*Y65+4*Y105+Y112-4*Y117',
        N_branch_discriminant_odd_coefficient_and_scalar_independent=True,
        exact_projected_N_cone_image_completion='C[[g1,...,g26,nu,b,z]]/(z^2-nu*b^2)',
        exact_projected_N_projective_image_completion='C[[g1,...,g25,nu,b,z]]/(z^2-nu*b^2)',
        N_projection_kernel_in_common_R29='(z61-a(Y))^2-nu_N(Y)*b(Y)^2; a,b,nu_N finite formal coefficient recipes',
        N_projection_is_not_smooth_graph=True,
        independent_formal_normal_coordinates_not_full_common_equation_expansion=True,
        full_odd_ideals_or_actual_D_completed_ring_not_claimed=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path)for path in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ACTUAL ODD LINEAR7 + COMMON R29 N PINCH PROJECTION PASS',result['seconds'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);main(ap.parse_args().output.resolve())
