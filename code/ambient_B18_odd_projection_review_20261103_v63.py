#!/usr/bin/env python3
"""Independent physical source and augmented-block odd-jet projection audit."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)


def mat(rows):return s.Matrix([[s.sympify(c) for c in row] for row in rows])
def same(A,B):assert A.shape==B.shape and (A-B).applyfunc(s.expand)==s.zeros(*A.shape)


def main(output,replay):
    started=time.monotonic();folder=ROOT/'reports'/STEM
    report=read(folder/'report.json');assert report['all_passed']
    for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in report['evidence_sha256'].items():assert sha(folder/name)==digest,name
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/ambient_B18_odd_firstjets_N_projection_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=40)
    assert (folder/'certificate.json.gz').read_bytes()==(replay/'certificate.json.gz').read_bytes()
    cert=read(folder/'certificate.json.gz')
    scopepath=ROOT/'reports/AMBIENT_B18_ODD_FIRSTJET_COORDINATE_SCOPE_20261103_v63.json'
    scope=read(scopepath);assert scope['all_passed']
    for name,digest in scope['input_sha256'].items():assert sha(ROOT/name)==digest,name
    oldpath=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    old=read(oldpath);Y=old['common28_linear_target_rows'];branch=old['common29th_target_row']
    h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);p=s.Symbol('p');K=s.symbols('k0:9')
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[a]*X[b]) for i,j in pairs for a,b in mons])
    def canonical(P):
        d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(name,f,g,H):
        P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h if name=='CA' else g
        P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
        if name=='CB':P[2,3]-=s.diff(f,x)*z/2;P[2,4]-=s.diff(f,x)*w/2
        return P-P.T
    swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
    def change(P):return (swap*P.subs({h:x,x:h},simultaneous=True)*swap.T).applyfunc(s.expand)
    def lie(M,P):return (M*P+P*M.T-sum(((M*X)[i]*P.diff(X[i]) for i in range(5)),s.zeros(5))).applyfunc(s.expand)
    def flag_derivative(M):
        L=M[:,2];L[2]=0
        ell=-M[1,:].T;ell[1]=0
        wedge=[]
        for i,j in pairs:
            v=-(M[0,i]*int(j==1)-M[0,j]*int(i==1))-(int(i==0)*M[1,j]-int(j==0)*M[1,i])
            wedge.append(0 if (i,j)==(0,1) else v)
        return s.Matrix(list(L)+list(ell)+wedge)
    D=s.Matrix([x*x,0,0,h*x,h*h+z*z]);ey=s.Matrix([0,0,1,0,0]);P0=ey*D.T-D*ey.T;v0=vec(P0)
    all_generators=[s.eye(5)[:,i]*s.eye(5)[j,:] for i,j in itertools.product(range(5),repeat=2)]
    selected=[all_generators[i] for i in old['orbit_GL_generator_indices']]
    O=s.Matrix.hstack(*[vec(lie(M,P0)) for M in selected]);HF=s.Matrix.hstack(*map(flag_derivative,selected))
    assert O.rank()==19
    cubics=[s.prod((h,x,z,w)[i] for i in q) for q in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in q)]
    assert list(map(str,cubics))==old['cubic_monomials']
    def source(name):
        basis=[change(canonical(raw(name,f,0,0))) for f in [h*h,h*x,x*x]]
        basis += [change(canonical(raw(name,0,1 if name=='CA' else h*h,0)))]
        basis += [change(canonical(raw(name,0,0,H))) for H in cubics]
        C=s.Matrix.hstack(*map(vec,basis));rows=list(C.T.rref()[1])
        M=s.zeros(5);M[1,0]=1;M[0,3]=-2 if name=='CA' else s.Rational(1,2)
        velocity=C[rows,:].inv()*vec(lie(M,P0))[rows,:];same(C*velocity,vec(lie(M,P0)))
        alpha=s.zeros(20,1);alpha[0]=-1;alpha[7]=-1;alpha[11]=1;alpha[16]=s.Rational(1,3)
        alpha-=p*velocity
        for k,j in zip(K,[1,2,3,9,10,14,15,18,19]):alpha[j]+=k
        same(alpha,s.Matrix(list(map(s.sympify,old['families'][name]['source_alpha_coefficients']))))
        f=alpha[0]*h*h+alpha[1]*h*x+alpha[2]*x*x
        H=sum(a*mon for a,mon in zip(alpha[4:,0],cubics));g=alpha[3] if name=='CA' else -alpha[3]*f
        fixed=change(canonical(raw(name,f,g,H)))
        G=s.eye(5)+p*M+p*p*M*M/2;Gi=s.eye(5)-p*M+p*p*M*M/2;same(G*Gi,s.eye(5))
        U=vec((G*fixed.subs(dict(zip(X,Gi*X)),simultaneous=True)*G.T).applyfunc(s.expand))
        same(U,s.Matrix(list(map(s.sympify,old['families'][name]['core_full150_coordinates']))))
        return U,basis,M,G,Gi,velocity
    def odd(U,parameters,coordinates):
        zero=dict.fromkeys([p,*parameters],0);same(U.subs(zero),v0);same(U.diff(p).subs(zero),s.zeros(150,1))
        J=O.row_join(U.jacobian(parameters).subs(zero))
        W=s.Matrix([s.Poly(f,p).nth(2).subs(dict.fromkeys(parameters,0)) for f in U])
        Z=s.Matrix([s.Poly(f,p).nth(3).subs(dict.fromkeys(parameters,0)) for f in U])
        mixed=s.zeros(150,19).row_join(U.diff(p).jacobian(parameters).subs(zero))
        # A single augmented inverse implements both even-coordinate changes.
        # It is independent of the producer's successive rem/q operations.
        augmented=J.row_join(W);inverse=augmented[coordinates,:].inv()
        linear_odd=mixed.row_join(Z)*inverse
        linear_odd=(linear_odd-augmented*inverse*linear_odd[coordinates,:]).applyfunc(s.expand)
        discriminator=inverse[-1,:]
        return J,W,Z,mixed,linear_odd,discriminator
    records={};discs=[];CA=None
    for name in ['CA','CB']:
        U,basis,M,G,Gi,velocity=source(name)
        J,W,Z,mixed,linear,disc=odd(U,list(K),Y+[branch])
        evidence=cert['CA_CB_full_firstjet_data'][name]
        for computed,key in [(W,'W'),(Z,'rawZ'),(mixed,'mixed'),(linear,'actual'),(disc,'discriminant')]:same(computed,mat(evidence[key]))
        same(J,mat(old['families'][name]['linear_J']))
        assert linear.rank()==7 and linear[:,0:28].rank()==7
        # Rank6 refers to Y at fixed ADJUSTED eta=r61-LY, not fixed raw r61.
        L=J[branch,:]*J[Y,:].inv()
        adjusted_Y=linear[:,:28]+linear[:,-1]*L
        same(adjusted_Y,mat(evidence['tangent']))
        assert adjusted_Y.rank()==6
        assert (linear*s.Matrix.hstack(*disc.nullspace())).rank()==6
        same(linear[Y+[branch],:],s.zeros(29,29))
        assert linear[old['common_target126_raw_coordinate_indices'][29:],:].rank()==7
        augmented=J.row_join(s.zeros(150,1)).col_join(HF.row_join(s.zeros(20,9)).row_join(flag_derivative(M)))
        same(augmented,mat(old['families'][name]['entire_source_tangent']));assert augmented.rank()==29
        discs.append(disc);records[name]=dict(whole29_odd_rank=7,
            raw_Y28_at_fixed_raw_row61_rank=7,
            adjusted_Y28_at_fixed_eta_rank=6,source_etale_augmented_rank=29)
        if name=='CA':CA=(U,basis,M,G,Gi,velocity,augmented)
    assert s.Matrix.vstack(*discs).rank()==2
    U,basis,M,G,Gi,velocity,augmented=CA;assert velocity[2]==0
    UN=U.subs(K[1],0);KN=[v for i,v in enumerate(K) if i!=1]
    sourceN=augmented[:,[i for i in range(29) if i!=20]];assert sourceN.rank()==28
    # Independently remove the moving flag and check the full N19 subspace.
    PU=s.zeros(5)
    for off,(i,j) in enumerate(pairs):
        PU[i,j]=sum(UN[15*off+k]*X[a]*X[b] for k,(a,b) in enumerate(mons));PU[j,i]=-PU[i,j]
    fixed=vec((Gi*PU.subs(dict(zip(X,G*X)),simultaneous=True)*Gi.T).applyfunc(s.expand))
    Nbasis=[basis[0],basis[1],*basis[3:]];NC=s.Matrix.hstack(*map(vec,Nbasis));nr=list(NC.T.rref()[1])
    same(NC*NC[nr,:].inv()*fixed[nr,:],fixed)
    J,W,Z,mixed,linear,disc=odd(UN,KN,Y)
    Ncert=cert['N_full_firstjet_data'];own_order=[i for i in Y if i!=27]+[27]
    perm=[own_order.index(i) for i in Y]
    same(linear,mat(Ncert['actual'])[:,perm]);same(disc,mat(Ncert['discriminant'])[:,perm])
    b=linear[61,:];same(b,mat(cert['N_odd_extra61_differential_in_actual_Y_order']))
    assert linear.rank()==7 and (linear*s.Matrix.hstack(*disc.nullspace())).rank()==6
    expected={16:-1,65:-6,105:4,112:1,117:-4}
    assert list(b)==[expected.get(row,0) for row in Y]
    assert b[Y.index(27)]==0 and disc[Y.index(27)]==-s.Rational(1,3)
    scalar=s.zeros(1,28);scalar[Y.index(106)]=1
    assert s.Matrix.vstack(disc,b).rank()==2 and s.Matrix.vstack(disc,b,scalar).rank()==3
    assert v0[106]==1
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'certificate.json.gz',oldpath,
           ROOT/'research/notes'/f'{STEM}.md',
           ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
           ROOT/'reports/AMBIENT_B18_ODD_DISCRIMINANT_REDUCTION_20261103_v63_REVIEW.json',scopepath,
           ROOT/'code/ambient_B18_odd_firstjet_coordinate_scope_20261103_v63.py',
           ROOT/'research/notes/AMBIENT_B18_ODD_FIRSTJET_COORDINATE_SCOPE_20261103_v63.md']
    result=dict(id=STEM+'_REVIEW',all_passed=True,
        fresh_producer_certificate_byte_identical=True,
        independent_full150_CA_CB_source_and_physical_GL19_reconstruction=True,
        both_CB_mu_f_cross_terms_and_CA_N_subbundle_all_orders_bound=True,
        independent_augmented_block_inverse_odd_firstjet_formula=True,
        full_mixed_p_source_and_cubic_implicit_corrections_checked=True,
        common29_odd_linear_ranks=records,independent_CA_CB_discriminant_differentials=True,
        N_exact_k1zero_source_etale28_and_full_fixed_flag_N19_membership=True,
        N_actual_base_is_commonY28_and_extra61_odd_differential_checked=True,
        N_discriminant_odd_and_scalar_differentials_independent=True,
        coordinate_scope_followup_bound_and_intrinsic_discriminant_tangent_ranks6_checked=True,
        whole_projected_N_completion_pinched_hypersurface_proof_accepted=True,
        cone_projected_image='C[[g26,nu,b,z]]/(z^2-nu*b^2)',
        projective_projected_image='C[[g25,nu,b,z]]/(z^2-nu*b^2)',
        partial_actual_projection_only_not_full_N_target_or_D_type=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','whole_projected_N_completion_pinched_hypersurface_proof_accepted','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
