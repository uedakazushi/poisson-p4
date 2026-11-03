#!/usr/bin/env python3
"""A monic-quintic cyclic open and explicit norm localization of knownfive.

The generic coefficients are finite inverse recipes, not serialized enormous
polynomials. Saved finite-field data certify the inverses; exact algebra proves
the generic statements. The complementary support remains unresolved.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from s61_affine_normal_projection_20261102 import P, pivot_minor
from s61_global_annihilator_review_20261103_v62 import resultant, trim
from s61_nonsquare_cubic_row_basis_review_20261103_v62 import determinant_mod

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def divide(f,g):
    f,g=trim(f),trim(g);q=[0]*max(0,len(f)-len(g)+1)
    inv=pow(g[-1],-1,P)
    for i in range(len(f)-len(g),-1,-1):
        c=f[i+len(g)-1]*inv%P;q[i]=c
        for j,v in enumerate(g):f[i+j]=(f[i+j]-c*v)%P
    return trim(q),trim(f)


def inverse(A):
    A=np.asarray(A,dtype=np.int64).copy()%P;n=A.shape[0];U=np.eye(n,dtype=np.int64)
    for j in range(n):
        k=j+int(np.flatnonzero(A[j:,j])[0]);A[[j,k]]=A[[k,j]];U[[j,k]]=U[[k,j]]
        c=pow(int(A[j,j]),-1,P);A[j]=A[j]*c%P;U[j]=U[j]*c%P
        factors=A[:,j].copy();factors[j]=0
        A=(A-factors[:,None]*A[j])%P;U=(U-factors[:,None]*U[j])%P
    assert np.array_equal(A,np.eye(n,dtype=np.int64));return U


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    cyclic=ROOT/'reports/S61_AFFINE_CYCLIC_OPEN_20261103_v62';saved=json.loads((cyclic/'report.json').read_text())
    assert saved['all_passed']
    for name,value in saved['input_sha256'].items():assert sha(ROOT/name)==value
    for name,value in saved['evidence_sha256'].items():assert sha(cyclic/name)==value
    lines=(cyclic/'native-certificate.txt').read_text().splitlines()
    chi=list(map(int,next(line for line in lines if line.startswith('CHI ')).split()[1:]))
    coordinates={int(line.split()[1]):list(map(int,line.split()[2:]))for line in lines if line.startswith('COORDINATE ')}
    values=np.fromstring((cyclic/'quadrics.in').read_text(),sep=' ',dtype=np.int64)
    assert values[:4].tolist()==[29,134,P,17];o=4
    high=values[o:o+134*435].reshape(134,435);o+=134*435
    cross=values[o:o+134*29*17].reshape(134,29,17);o+=134*29*17
    base=values[o:o+134*153].reshape(134,153);o+=134*153;v=values[o:]
    assert len(v)==17 and int(v[0])==P-3
    pairs=list(itertools.combinations_with_replacement(range(29),2));bpairs=list(itertools.combinations_with_replacement(range(17),2))
    yp=np.array([int(v[i])*int(v[j])%P for i,j in bpairs],dtype=np.int64)
    linear=cross@v%P;constant=base@yp%P
    native_input=args.output/'probe.in'
    native_input.write_text('1031\n'+' '.join(map(str,chi))+'\n'+'\n'.join(' '.join(map(str,coordinates[i]))for i in range(29))+'\n'+
        '\n'.join(' '.join(map(str,a.ravel()))for a in[high,linear,constant])+'\n')
    source=ROOT/'code/s61_affine_quintic_open_probe_20261103_v62.cpp';binary=args.output/'probe'
    flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
    compiled=subprocess.run(['c++','-O2','-std=c++17',str(source),*flags,'-o',str(binary)],capture_output=True,text=True,check=True,timeout=60)
    (args.output/'compile.log').write_text(compiled.stdout+compiled.stderr)
    process=subprocess.run([str(binary.resolve()),str(native_input.resolve()),str((args.output/'native.txt').resolve())],capture_output=True,text=True,check=True,timeout=30)
    (args.output/'native.log').write_text(process.stdout+process.stderr);assert 'MONIC_QUINTIC_HIGH_BLOCK_AND_BEZOUT_PASS 0 16461' in process.stdout
    output=(args.output/'native.txt').read_text().splitlines()
    def polynomial(name):return list(map(int,next(line for line in output if line.startswith(name+' ')).split()[1:]))
    r,g,u,w=[polynomial(key)for key in['RESIDUAL','GCD','BEZOUT_CHI','BEZOUT_RESIDUAL']]
    assert len(r)==1031 and len(g)==6 and g[-1]==1 and len(u)<=1025 and len(w)<=1026
    combination=np.convolve(np.array(u,dtype=np.int64),np.array(chi,dtype=np.int64))
    second=np.convolve(np.array(w,dtype=np.int64),np.array(r,dtype=np.int64))
    total=np.zeros(max(len(combination),len(second)),dtype=np.int64)
    total[:len(combination)]+=combination;total[:len(second)]+=second
    assert trim(total%P)==g
    q1,rem1=divide(chi,g);q2,rem2=divide(r,g);assert not rem1 and not rem2
    assert resultant(q1,q2)==16461
    derivative=[i*g[i]%P for i in range(1,6)];assert resultant(g,derivative)==19651
    assert divide(coordinates[0],g)[1]==[1]
    gcd_degrees=polynomial('INDIVIDUAL_GCD_DEGREES');assert gcd_degrees==[5]*134
    # The same original full Jacobian gives a cofactor inverse for all16
    # remainder normal directions after transport through the cyclic map.
    known_path=ROOT/'reports/S61_AFFINE_ENLARGED_SECTION_20261103_v62/report.json'
    known=json.loads(known_path.read_text())['known_five_mod_32003']
    scale=pow(int(v[0]),-1,P);Y=v*scale%P;Js=[];Zs=[]
    ymon=np.array([int(Y[i])*int(Y[j])%P for i,j in bpairs],dtype=np.int64)
    for point in known:
        Z=np.array(point+[0]*25,dtype=np.int64)*scale%P;Zs.append(Z.tolist())
        zm=np.array([int(Z[i])*int(Z[j])%P for i,j in pairs],dtype=np.int64)
        assert np.all((high@zm+np.einsum('rkb,k,b->r',cross,Z,Y,optimize=True)+base@ymon)%P==0)
        JZ=cross@Y%P;JY=np.einsum('rkb,k->rb',cross,Z)%P
        for j,(a,b)in enumerate(pairs):
            JZ[:,a]=(JZ[:,a]+high[:,j]*Z[b])%P;JZ[:,b]=(JZ[:,b]+high[:,j]*Z[a])%P
        for j,(a,b)in enumerate(bpairs):
            JY[:,a]=(JY[:,a]+base[:,j]*Y[b])%P;JY[:,b]=(JY[:,b]+base[:,j]*Y[a])%P
        J=np.concatenate([JZ,JY[:,1:]],axis=1);assert len(pivot_minor(J)[0])==45;Js.append(J)
    rng=np.random.default_rng(621605)
    for trial in range(20):
        C=rng.integers(0,17,size=(45,134),dtype=np.int64)
        determinants=[determinant_mod((C@J%P).tolist(),P)for J in Js]
        if all(determinants):break
    else:raise AssertionError('no common Jacobian combination')
    assert trial==0 and determinants==[20321,16946,11970,18492,6923]
    norm=1
    for J,d in zip(Js,determinants):
        L=inverse(C@J%P)@C%P;assert np.array_equal(L@J%P,np.eye(45,dtype=np.int64))
        B=L[29:];assert not np.any(B@J[:,:29]%P)and np.array_equal(B@J[:,29:]%P,np.eye(16,dtype=np.int64))
        norm=norm*d%P
    assert norm==25750
    (args.output/'jacobian-certificate.json').write_text(json.dumps(dict(prime=P,chart='Y1=1',scale=scale,
        shifted_center=(v*scale%P).tolist(),known_Z_scaled=Zs,common_combination=C.tolist(),
        determinants=determinants,norm=norm,all_left_inverse_identities_checked=True),indent=2)+'\n')
    (args.output/'finite-recipe.json').write_text(json.dumps(dict(
        chosen_original_residual=0,high_coefficient_block=dict(rows='coefficients s^5,...,s^2055',
            first_columns='s^j*chi, 0<=j<1025',second_columns='s^j*r0,0<=j<1026',size=2051,
            nonzero_determinant_witness_cofactor_resultant=16461),
        monic_relation='Solve highblock*(u,v)=e_s5; g=u*chi+v*r0 is monic degree5',
        exact_module_presentation_shape=[5,675],retained_residuals='rem_g(chi), rem_g(F1(P)),...,rem_g(F134(P))',
        normal_recognition_recipe=['Y1=1, yk=Y_(k+1)-v_(k+1)/v1',
            'F=R[s]/g, F0=Rbeta[s]/g_at_y0; shrink beta by Disc(g0) and Norm(det(C*Jac_original))',
            'B=last16rows((C*Jac_original)^(-1)*C) overF0; lift its degree<5 coefficients toF',
            'rj=sum_i Bji*rem_g(Fi(P)), divide coefficientwise by y using telescoping after clearing unit denominators',
            'rj=sum_k Mjk*yk, M_at_y0=I16; invert Norm_F/R(det M)',
            'Then residualideal=(y), C=F0, t=1/v1; affine chart is the knownfive algebra'],
        global_open='Product of setup minors, Y1, cyclicDelta, highblockdet, Disc(g), Norm(t), Jacobian norm, denominator norms, Norm(detM)',
        all_open_factors_are_units_at_known_five=True,complement_resolved=False,
        generic_coefficients_representation='Exact finite inverses and norm recipes; no expanded generic coefficient table'),indent=2)+'\n')
    binary.unlink()
    inputs=[Path(__file__),source,cyclic/'report.json',cyclic/'native-certificate.txt',cyclic/'quadrics.in',known_path,
        ROOT/'code/s61_global_annihilator_review_20261103_v62.py',ROOT/'code/s61_nonsquare_cubic_row_basis_review_20261103_v62.py',
        ROOT/'research/notes/S61_AFFINE_GLOBAL_ACTION_SPREAD_20261103_v62.md',
        ROOT/'research/notes/S61_AFFINE_QUINTIC_OPEN_20261103_v62.md']
    result=dict(id='S61_AFFINE_QUINTIC_OPEN_20261103_v62',all_passed=True,prime=P,
        first_original_residual_gcd_degree=5,monic_quintic_relation_proved_on_explicit_cyclic_open=True,
        high_block_size=2051,high_block_nonzero_cofactor_resultant=16461,
        exact_small_module_presentation=[5,675],all_original_equations_retained=134,
        common_Jacobian_norm=25750,squarefree_quintic_discriminant=19651,
        explicit_norm_localized_normal_recognition_proved=True,
        p_integral_finite5_direction_module_on_specified_open=True,
        knownfive_finite_free_over_beta_after_further_norm_localization=True,
        generic_polynomial_coefficients_serialized=False,generic_polynomial_recipe_defined=True,
        global_full_affine_exhaustion_claimed=False,complement_of_explicit_open_resolved=False,
        limitations=['All reduced or nonreduced extra support is confined to the complementary denominator/norm locus; that locus is not excluded.',
            'Generic quintic and normal-annihilator coefficients are exact finite inverse recipes, not expanded stored polynomials.',
            'The proof uses the previously established universal knownfive lower bound and homogeneous corrected cyclic presentation.'],
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
        evidence_sha256={p.name:sha(p)for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS explicit monicquintic open;5x675;norm-localizedknownfive;complementOPEN;',result['seconds'])


if __name__=='__main__':main()
