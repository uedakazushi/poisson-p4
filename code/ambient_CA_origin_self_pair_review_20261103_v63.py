#!/usr/bin/env python3
"""Independent physical polynomial tensors/flags for the actual CA selfpair."""
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
STEM='AMBIENT_CA_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63'
def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:return json.load(f)
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()
def same(A,B):assert A.shape==B.shape and (A-B).applyfunc(s.expand)==s.zeros(*A.shape)


def main(output,replay):
    started=time.monotonic();folder=ROOT/'reports'/STEM;report=read(folder/'report.json');assert report['all_passed'] and report['actual_polynomial_self_pair_family_proved']
    for name,digest in report['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in report['evidence_sha256'].items():assert sha(folder/name)==digest,name
    subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'code/ambient_CA_origin_polynomial_self_pair_20261103_v63.py'),
                    '--output',str(replay)],check=True,timeout=30)
    assert (folder/'certificate.json.gz').read_bytes()==(replay/'certificate.json.gz').read_bytes()
    cert=read(folder/'certificate.json.gz');h,x,y,z,w,a=s.symbols('h x y z w a');X=s.Matrix([h,x,y,z,w])
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    def vec(P):return s.Matrix([s.Poly(P[i,j],*X).coeff_monomial(X[i0]*X[i1]) for i,j in pairs for i0,i1 in mons])
    def canonical(P):
        d=s.Matrix([sum(s.diff(P[i,j],X[i]) for i in range(5)) for j in range(5)])
        return (P-(X*d.T-d*X.T)/5).applyfunc(s.expand)
    def raw(f,g,H):
        P=s.zeros(5);P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
        return P-P.T
    sw=s.eye(5);sw[0,0]=sw[1,1]=0;sw[0,1]=sw[1,0]=1
    def change(P):return (sw*P.subs({h:x,x:h},simultaneous=True)*sw.T).applyfunc(s.expand)
    cubics=[s.prod((h,x,z,w)[i] for i in mon) for mon in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in mon)]
    basis=[change(canonical(raw(f,0,0))) for f in [h*h,h*x,x*x]]+[change(canonical(raw(0,1,0)))]+[change(canonical(raw(0,0,H))) for H in cubics]
    C=s.Matrix.hstack(*map(vec,basis));criticalpath=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';critical=read(criticalpath)
    same(C,s.Matrix(critical['fixed_fiber_linear_basis150_by20']))
    ey=s.Matrix([0,0,1,0,0]);D0=s.Matrix([x*x,0,0,h*x,h*h+z*z]);P0=ey*D0.T-D0*ey.T
    D=s.Matrix([x*x+12*a*h*x/5,-3*a*x*x/5,0,h*x+2*a*x*z/5,h*h+z*z-8*a*x*w/5])
    P=ey*D.T-D*ey.T
    same(vec(P),vec(P0)+a*(-3*C[:,1]-C[:,9]))
    decode=lambda terms:sum(s.Rational(c)*a**degree for degree,c in terms)
    same(vec(P),s.Matrix(list(map(decode,cert['first_source_physical150']))))
    M=s.zeros(5);M[1,0]=1;M[0,3]=-2;assert M**3==s.zeros(5) and M**2!=s.zeros(5)
    G=s.eye(5)+a*M+a*a*M*M;Gi=s.eye(5)-a*M
    same(G*Gi,s.eye(5));same(Gi*G,s.eye(5));assert G.det()==1
    savedG=s.Matrix([[decode(q) for q in row] for row in cert['group_g5_by5']]);same(G,savedG)
    savedGi=s.Matrix([[decode(q) for q in row] for row in cert['group_inverse5_by5']]);same(Gi,savedGi)
    second=(G*P.subs(dict(zip(X,Gi*X)),simultaneous=True)*G.T).applyfunc(s.expand)
    reverse=(Gi*second.subs(dict(zip(X,G*X)),simultaneous=True)*Gi.T).applyfunc(s.expand)
    same(reverse,P)
    same(vec(second),s.Matrix(list(map(decode,cert['second_source_physical150']))))
    rows=list(C.T.rref()[1]);inverse=C[rows,:].inv();second_coeff=inverse*vec(second)[rows,:]
    same(C*second_coeff,vec(second));same(second_coeff,s.Matrix(list(map(decode,cert['second_source_fixed20']))))
    assert all(s.Poly(q,a).degree()<=2 for q in second_coeff)
    # Direct symbolic physical derivative identities, independently of the
    # producer's sparse coefficient-ring differentiation.
    for tensor in [P,second]:
        divergence=[s.expand(sum(s.diff(tensor[i,j],X[i]) for i in range(5))) for j in range(5)]
        assert divergence==[0]*5
        for i,j,k in itertools.combinations(range(5),3):
            jac=s.expand(sum(tensor[i,l]*s.diff(tensor[j,k],X[l])+tensor[j,l]*s.diff(tensor[k,i],X[l])+tensor[k,l]*s.diff(tensor[i,j],X[l]) for l in range(5)))
            assert jac==0
    first_flagH=s.Matrix([[0,1,0,0,0]]);secondH=first_flagH*G;secondW=G[:2,:];secondL=Gi*ey
    same(secondL,ey);same(secondH,s.Matrix([[a,1,0,-2*a*a,0]]));same(secondW[0,:],s.Matrix([[1,0,0,-2*a,0]]))
    assert secondW[:,[0,1]].det()==1 and secondH[0]==a
    same(secondW*secondL,s.zeros(2,1))
    # H' wedge H0 has nonzero coefficient a whenever a is nonzero.
    assert s.det(s.Matrix([[secondH[0],secondH[1]],[0,1]]))==a
    line=s.Matrix.hstack(vec(P0),vec(P).diff(a));assert line[[16,20],:].det()==-s.Rational(12,5)
    assert vec(P)[20]==-1
    # Projective Poisson rank is exactly2, not merely affine-rank<=2.
    assert s.expand(h*P[1,2]-x*P[0,2]+y*P[0,1]-x*x*(x+3*a*h))==0
    for i,j,k,l in itertools.combinations(range(5),4):assert s.expand(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])==0
    accepted=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';assert read(accepted)['all_passed']
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}_REVIEW.md',folder/'report.json',folder/'certificate.json.gz',criticalpath,accepted,
           ROOT/'research/notes'/f'{STEM}.md']
    result=dict(id=STEM+'_REVIEW',all_passed=True,fresh_author_certificate_byte_identical=True,
        independent_raw_CA_canonical20basis_and_explicit_physical_tensor_reconstructed=True,
        actual_full_QQa_finite5matrix_push_and_reverse_push_checked=True,
        all150_exact_second_CA_source_membership_and_coefficients_checked=True,
        both_sources_all25_divergence_and350_strict_Jacobi_coefficients_zero=True,
        actual_two_flags_same_tensor_and_distinct_for_every_a_nonzero=True,
        flag_Hprime='x+a*h-2*a^2*z',group_determinant1_and_no_parameter_denominator=True,
        first_source_nonconstant_projective_line_and_nonzero_every_a=True,
        entire_family_generic_projective_Poisson_rank2_not4=True,
        actual_polynomial_self_pair_family_proved_without_jet_extrapolation=True,
        source_normalization_immersions_or_full_fiber_count_not_inferred=True,
        no_entire_conductor_full_origin_I_or_branch_classification_claim=True,
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['all_passed','actual_polynomial_self_pair_family_proved_without_jet_extrapolation','seconds']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True);args=parser.parse_args();main(args.output.resolve(),args.replay.resolve())
