#!/usr/bin/env python3
"""Independent whole first nil31-long quartic22 component and prior96 audit.

No producer, producer helpers, CAS, NumPy or legacy imports. All90source,
Laurent normalizer and full150/350 matrices use prior independent arithmetic.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

from classification_minimal_quartic_review import LaurentRing,determinant
from classification_quadratic_cofactor_review import verify_source
from classification_nil22_image_components_review import substitute_sparse,derivative_parameter,evaluate
import classification_translated_minimal_quartic_review as exact

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(it.combinations(range(5),2))
PARAMETERS=('a','n3','n4','n8','n12')


def matrix(R,entries):
    P=[[{}for _ in range(5)]for _ in range(5)]
    for i,j in PAIRS:P[i][j],P[j][i]=entries.get((i,j),{}),R.scale(entries.get((i,j),{}),-1)
    return P


def strict(R,P):
    J=[R.add(*(R.mul(P[a][k],R.derivative(P[b][c],k))for a,b,c in [(i,j,l),(j,l,i),(l,i,j)]for k in range(5)))for i,j,l in it.combinations(range(5),3)]
    assert J==[{}]*10


def modular(R,P):
    return [R.add(*(R.derivative(P[i][j],i)for i in range(5)))for j in range(5)]


def point(R,P,values,scale=Q(1)):
    return exact.arithmetic.antisymmetric({ij:value for ij in PAIRS if(value:=evaluate(R,P[ij[0]][ij[1]],values,scale))})


def whole_source(original,support,normal):
    original_check=verify_source(original,5)
    assert original_check['rank']==75 and len(original_check['free_indices'])==15
    R=LaurentRing(PARAMETERS,('h','x','y','z','w'),invert=('a',))
    oldR=LaurentRing(tuple(f'n{i}'for i in range(15)),('h','x','y','z','w'))
    parameter_texts=['54*a*n12*n4**2-18*a*n3*n4','-27*a*n4**2','0','n3','n4','36*a*n8**2','0','-12*a*n8','n8','0','18*a*n12*n4-6*a*n3','-9*a*n4','n12','-1/6','a']
    bindings=list(map(R.parse,parameter_texts))
    saved=support['all_two_complete_families'][0]
    for i,value in enumerate(bindings):assert value==R.parse(saved['complete_source_substitutions'][f'n{i}'])
    spatial=[R.variable(name)for name in ('h','x','y','z','w')]
    entries={ij:substitute_sparse(oldR,oldR.parse(text),R,bindings,spatial,True)for ij,text in zip(it.combinations(range(1,5),2),original['brackets'])}
    P=matrix(R,entries);strict(R,P)
    N=modular(R,P)
    assert N==[{},R.variable('y'),R.variable('z'),R.variable('h'),{}]
    canonical_texts=['x*z/2-y**2/3','h*x/2-y*z/6','-w*y/6-3*x*z+2*y**2','h*y/3-z**2/6','-3*h*x-w*z/6+y*z','-h*w/6-2*h*y+z**2']
    canonical=matrix(R,dict(zip(it.combinations(range(1,5),2),map(R.parse,canonical_texts))))
    assert all(canonical[i][j]==R.parse(normal['canonical_homogeneous_brackets'].get(f'{i},{j}','0'))for i,j in PAIRS)
    Gtexts=[['1','0','0','0','0'],['-108*a*n4*n8-12*n12*n4-2*n3','1','0','-3*n8','-3*a*n4'],['-3*n8','0','1','0','0'],['0','0','0','1','0'],['36*a*n8+6*n12','0','0','0','a']]
    invtexts=[['1','0','0','0','0'],['-6*n12*n4+2*n3','1','0','3*n8','3*n4'],['3*n8','0','1','0','0'],['0','0','0','1','0'],['-36*n8-6*n12/a','0','0','0','1/a']]
    G=[[R.parse(v)for v in row]for row in Gtexts];inverse=[[R.parse(v)for v in row]for row in invtexts]
    for i,j in it.product(range(5),repeat=2):
        assert G[i][j]==R.parse(normal['canonical_to_universal_affine_matrix'][i][j])
        assert inverse[i][j]==R.parse(normal['universal_to_canonical_affine_matrix'][i][j])
        assert R.add(*(R.mul(G[i][k],inverse[k][j])for k in range(5)))==R.constant(i==j)
    assert determinant(R,G)==R.variable('a')
    inverse_coordinates=[R.add(*(R.mul(inverse[i][j],spatial[j])for j in range(5)))for i in range(5)]
    substituted=[[substitute_sparse(R,canonical[i][j],R,[R.variable(a)for a in PARAMETERS],inverse_coordinates)for j in range(5)]for i in range(5)]
    actual=[[R.add(*(R.mul(R.mul(G[i][k],substituted[k][l]),G[j][l])for k,l in it.product(range(5),repeat=2)))for j in range(5)]for i in range(5)]
    assert actual==P
    point_values=dict(zip(PARAMETERS,map(Q,[1,0,0,0,0])))
    witness=point(R,P,point_values,Q(6))
    derivatives=[point(R,[[derivative_parameter(P[i][j],k)for j in range(5)]for i in range(5)],point_values,Q(6))for k in range(5)]
    pf=R.add(R.mul(canonical[1][2],canonical[3][4]),R.scale(R.mul(canonical[1][3],canonical[2][4]),-1),R.mul(canonical[1][4],canonical[2][3]))
    F=R.parse('9*h**2*x**2-18*h*x*y*z+6*x*z**3+8*h*y**3-3*y**2*z**2')
    assert pf==R.scale(F,Q(1,6))
    assert R.parse('9*x**2-18*x*y*z+6*x*z**3+8*y**3-3*y**2*z**2')==R.parse('9*(x-y*z+z**3/3)**2+8*(y-z**2/2)**3')
    assert any(ex[R.n]==0 for ex in F)
    # Optional geometric identification; this is not used for distinctness.
    assert R.scale(F,-108)==R.parse('(3*z)**2*(6*y)**2-4*h*(6*y)**3-4*(3*z)**3*(6*x)-27*h**2*(6*x)**2+18*h*(3*z)*(6*y)*(6*x)')
    # Independently recover the old nil4 cofactor from its actual tensor.
    oldtexts=['36*h*w+6*x*z-8*y**2','6*w*x-10*y*z','-10*w*y-z**2','6*w*y-9*z**2','-12*w*z','-12*w**2']
    oldP=matrix(R,dict(zip(it.combinations(range(1,5),2),map(R.parse,oldtexts))))
    assert modular(R,oldP)==[{},R.parse('36*y'),R.parse('36*z'),R.parse('36*w'),{}]
    return witness,derivatives,dict(original_full90_source=original_check,all_five_original_source_parameters=True,
        universal_strict10_Jacobi_zero=True,universal_h_Casimir_modular=['0','y','z','h','0'],
        all150_actual_PGL_normalizer_coefficients_reconstructed=True,normalizer_determinant='a',
        only_parameter_inverted='a',a_zero_complement='Already excluded completely for H+irreduciblequartic by the accepted fullsupport review.',
        whole_first_five_parameter_source_is_one_actual_PGL_orbit=True,canonical_integer_scale=6,
        canonical_affine_Pfaffian='(3/2)*(x-yz+z^3/3)^2+(4/3)*(y-z^2/2)^3',
        irreducibility='Polynomial automorphism (u,v,z,w)=(x-yz+z³/3,y-z²/2,z,w) gives quadratic in u with nonsquare v³ discriminant; homogeneous F has no h factor.',
        actual_quartic_factor='9h²x²-18hxyz+6xz³+8hy³-3y²z²',
        actual_quintic_reduced=True,unique_simple_hyperplane='h=0',full_source_including_all_other_parameter_zero_values=True)


def certify(M):
    rank,cols=exact.arithmetic.rank(M)
    rows=exact.arithmetic.rank([[row[j]for row in M]for j in cols])[1]
    determinant=exact.minor(M,rows,cols);assert determinant
    return dict(shape=[len(M),len(M[0])],rank=rank,rows=rows,columns=cols,determinant=str(determinant),
        full_matrix_sha256=hashlib.sha256(json.dumps(M,default=str,separators=(',',':')).encode()).hexdigest())


def ambient(witness,derivatives,saved):
    # All coefficients of the actual frozen integer point are checked first.
    ring=LaurentRing((),('h','x','y','z','w'))
    expected=point(ring,matrix(ring,{ij:ring.parse(saved['witness_brackets'].get(f'{ij[0]},{ij[1]}','0'))for ij in PAIRS}),{})
    assert witness==expected
    D,G,rows=exact.matrices(witness)
    columns=[exact.vector(p)for p in derivatives]
    F=[row+[col[i]for col in columns]for i,row in enumerate(G)]
    assert all(sum(v*col[j]for j,v in row)==0 for row in rows for col in columns)
    certificates={name:certify(M)for name,M in [('D',D),('vertical',[row[:25]for row in G]),('GL5_vertical_scalar',G),('whole_actual_five_source',F)]}
    assert [a['rank']for a in certificates.values()]==[103,24,47,47]
    return dict(independent_full_matrices=certificates,all350_Jacobirows_and150variables=True,
        full_D_times_all_actual_GL5_vertical_scalar_and_source_directions_zero=True,
        full_kernel_dimension=47,projective_orbit_dimension=22,projective_tangent_dimension=22,
        regular_smooth_point_of_full_Jacobi_scheme=True,local_unique_full_component=True,
        entire_PGL_orbit_closure_is_a_full_component=True,projective_component_dimension=22,
        whole_first_five_source_covered_by_this_component=True,dense_open_actual_PGL_orbit_dimension=22,rational_quotient_dimension=0)


def cofactor_flag():
    new=[[Q(0)for j in range(5)]for i in range(5)];old=[[Q(0)for j in range(5)]for i in range(5)]
    new[1][2]=new[2][3]=new[3][0]=Q(1)
    old[1][2]=old[2][3]=old[3][4]=Q(1)
    def mul(A,B):return [[sum(A[i][k]*B[k][j]for k in range(5))for j in range(5)]for i in range(5)]
    h=[Q(i==0)for i in range(5)]
    records=[]
    for label,N in [('new',new),('old_nil4',old)]:
        square=mul(N,N);cube=mul(square,N);fourth=mul(cube,N)
        assert [exact.arithmetic.rank(M)[0]for M in [N,square,cube]]==[3,2,1]
        assert all(not x for row in fourth for x in row)
        transpose=[list(row)for row in zip(*N)];cubeT=[list(row)for row in zip(*cube)]
        rank1=exact.arithmetic.rank([row+[v]for row,v in zip(transpose,h)])[0]
        rank3=exact.arithmetic.rank([row+[v]for row,v in zip(cubeT,h)])[0]
        assert (rank1,rank3)==((3,1)if label=='new'else(4,2))
        records.append(dict(label=label,N=[[str(v)for v in row]for row in N],Jordan_type=[4,1],
            h_in_image_N_transpose=rank1==3,h_in_image_N_transpose_cubed=rank3==1))
    return dict(whole_source_unique_H_intrinsic_cofactor_flags=records,
        new_h_in_third_dual_image_old_h_outside_first_dual_image=True,
        invariance='The unique simple H determines the unique H-Casimir homogeneous representative. Its linear divergence transforms by conjugacy under PGL and by scalar under tensor scaling, so the covector-image incidence is intrinsic.',
        same_dimension_component_distinction='Both component closures have dense actual22 PGL orbits. Equality would force these two dense locallyclosed orbits to intersect, contradicting their intrinsic H-N flags. No global closedness of the image flag is assumed.')


def prior_audit(snapshot):
    rows=snapshot['components'];assert len(rows)==96 and len({a['id']for a in rows})==96
    exclusions={};same=[];nonreduced=[]
    for row in rows:
        name=row['id'];assert type(row['projective_dimension'])is int
        if row['generic_rank']==2:exclusions[name]='closed global rank<=2 versus new nonzero rank4 quintic';continue
        if row['projective_dimension']!=22:exclusions[name]='different full-component projective dimension after componenthood';continue
        same.append(name);divisor=row['generic_actual_degeneracy_divisor'];assert divisor
        if divisor['reduced']is False:
            assert max(divisor['corresponding_multiplicities'])>=2
            exclusions[name]='old whole closure satisfies closed nonreducedquintic condition; new hF4 is reduced'
            nonreduced.append(name)
        elif name=='R4_TWO_HYPERPLANES_CUBIC_22':
            assert divisor['irreducible_factor_degrees']==[1,1,3]
            exclusions[name]='old whole closure in properclosed 1+1+3 factor image; new has irreduciblequartic factor'
        else:
            assert name=='R4_NIL4_HYPERPLANE_QUARTIC_22'
            assert divisor['irreducible_factor_degrees']==[1,4]
            exclusions[name]='old dense22orbit uniqueH covector outside imNdual, new inside im(Ndual)^3'
    assert len(same)==15 and len(nonreduced)==13
    return dict(prior_component_count=96,all_prior_dimensions_present=True,same22_rank4_components=same,
        thirteen_old_nonreduced22_components=nonreduced,exclusions_for_every_prior96_component=exclusions,
        distinct_from_every_prior96_component=True,
        closed_conditions='Nonreducedquintic is union of properprojective images (G,R)->G²R for degG1or2; 1+1+3 is a properprojective multiplicationimage. Include zero via affinecones.',
        adopted_dependency='Prior componenthood/dimensions/genericfactor data are the audited immutable prior96; only the newlyneeded oldnil4 H-cofactor flag is recomputed here.',
        publication_novelty_claimed=False)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(170,175));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    original_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    support_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_SUPPORT_20261107.json'
    normal_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_NORMALFORM_20261107.json'
    ambient_path=ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_ORBIT1_AMBIENT_20261107.json'
    catalogue_path=ROOT/'reports/CLASSIFICATION_QUARTIC_NIL31_LONG_PRIOR96_20261107.json'
    normal=json.loads(normal_path.read_text())['nil31_long_normalform']['all_two_actual_orbit_normalforms'][0]
    witness,derivatives,source=whole_source(json.loads(original_path.read_text())['full_nilpotent_sources'][4],json.loads(support_path.read_text())['nil31_long_support'],normal)
    full_ambient=ambient(witness,derivatives,json.loads(ambient_path.read_text())['source_point'])
    paths=[Path(__file__),ROOT/'code/classification_quadratic_cofactor_review.py',ROOT/'code/classification_minimal_quartic_review.py',ROOT/'code/classification_nil22_image_components_review.py',ROOT/'code/classification_translated_minimal_quartic_review.py',ROOT/'code/classification_nil22_weighted_slice_review.py',original_path,support_path,normal_path,ambient_path,catalogue_path,
        ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_LONG_REVIEW_20261107.json',ROOT/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL_REGULAR_COMPONENT_REVIEW_20261107.json']
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_or_producer_helper_import=True,no_CAS_or_NumPy_import=True,
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
        complete_first_source=source,full_ambient_component_proof=full_ambient,
        intrinsic_H_cofactor_distinction=cofactor_flag(),all_prior96_audit=prior_audit(json.loads(catalogue_path.read_text())),
        accepted_new_full_component_relative_to_prior96=True,projective_dimension=22,generic_rank=4,
        second_source_component_membership_claimed=False,reviewer_source_bytes=Path(__file__).read_text(),
        seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as handle:json.dump(result,handle,indent=2);handle.write('\n')
    print(json.dumps({k:result[k]for k in ['all_passed','projective_dimension','seconds','maximum_RSS_KiB']}))


if __name__=='__main__':main()
