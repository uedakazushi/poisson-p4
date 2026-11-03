#!/usr/bin/env python3
"""One fixed accepted NEW25 point, rebuilt via DAG and independent QQ(s) pairs.

No author point module is imported. Exactly sample[0] is evaluated; no new
point is selected and no invariant evaluation matrix or degree probe is made.
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

from q_firstjet_spin4_mixed_discriminant_review_20261103_v63 import (
    mat, mm, tp, det, inv,
)

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63_REVIEW'


def read(path):return json.loads(gzip.decompress(path.read_bytes())if path.suffix=='.gz'else path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(output):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63/report.json'
    acpath=author.parent/'certificate.json.gz';ar=read(author);ac=read(acpath);assert ar['all_passed']
    samplepath=ROOT/'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63/samples.json'
    samplereview=ROOT/'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63_REVIEW/report.json'
    mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
    mapreview=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.json'
    primitivepath=ROOT/'reports/PRIMITIVE_JET_ELLIPTIC_CURVE_20261101.json'
    origpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    construction=ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63_REVIEW/report.json'
    slice_cert=ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63/certificate.json'
    circuitpath=slice_cert.parent/'mixed-discriminant-circuit.json.gz'
    checked={}
    for p in [author,samplereview,mapreview,construction]:
        report=read(p);assert report['all_passed']
        for name,want in report['input_sha256'].items():assert sha(ROOT/name)==want,name;checked[name]=want
    samples=read(samplepath);sample=samples['actual_QQ_source_points'][0]
    assert ac['sample_index']==0 and ac['accepted_source_point']==sample
    point=sample['point'];payload=read(mappath);values=[]
    for i,node in enumerate(payload['nodes']):
        op,value=node
        if op=='q':q=Q(value)
        elif op=='v':q=Q(point[payload['input_names'][value]])
        elif op=='+':assert all(0<=j<i for j in value);q=sum((values[j]for j in value),Q(0))
        elif op=='*':
            assert all(0<=j<i for j in value);q=Q(1)
            for j in value:q*=values[j]
        elif op=='inv':assert 0<=value<i and values[value];q=1/values[value]
        else:raise AssertionError(op)
        values.append(q)
    coordinates=[values[i]for i in payload['output_u1_through_u61']];assert len(coordinates)==61
    vp=values[payload['normalization']['value_pivot_node']]
    pf=values[payload['normalization']['value_pfaffian_node']]
    assert vp==Q(sample['value_pivot'])==Q(2051,12)and pf==Q(sample['value_pfaffian'])==Q(117727,2592)
    # A different beta reconstruction from the S61 output, rather than simply
    # accepting the stored16-vector used by the author point program.
    primitive=read(primitivepath);assert primitive['all_passed']
    for name,want in primitive['source_sha256'].items():assert sha(ROOT/name)==want,name;checked[name]=want
    z=s.symbols(' '.join(primitive['affine_variables']))
    substitutions={s.Symbol('u'+str(i+1)):s.Rational(q.numerator,q.denominator)for i,q in enumerate(coordinates)}
    pairs=list(itertools.combinations(range(4),2));coords=[(ij,k)for ij in pairs for k in range(4)]
    beta=[]
    for (i,j),k in coords:
        expression=s.sympify(primitive['primitive_beta'][i][j]).coeff(z[k]).subs(substitutions)
        beta.append(Q(str(expression)))
    B=mat(read(origpath)['primitive24_by16_basis'])
    assert beta==[sum(c*Q(q)for c,q in zip(row,sample['beta16']))for row in B]
    assert beta==list(map(Q,ac['beta24']))
    psi=[[beta[coords.index((ij,k))]for k in range(4)]for ij in pairs[:-1]]
    assert psi==mat(ac['psi5x4'])
    chi=[(-1)**i*det([row for j,row in enumerate(psi)if j!=i])for i in range(5)]
    spin=read(spinpath);metric=mat(spin['quadratic_form5']);metric_inv=inv(metric)
    normal=[sum(c*x for c,x in zip(row,chi))for row in metric_inv]
    q=sum(normal[i]*metric[i][j]*normal[j]for i in range(5)for j in range(5))
    J=mat([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]);Ji=inv(J)
    gram=mm(mm(tp(psi),metric),psi);p0=det(mm(gram,Ji))
    assert q==Q(-10025632819168382287872,1019784129324739375)<0 and q==16*p0
    assert chi==list(map(Q,ac['normal_cofactor_chi']))and normal==list(map(Q,ac['normal_vector5']))
    assert p0==Q(ac['p0_detC'])
    # The quadratic algebra is represented by plain Fraction pairs, not the
    # author's overloaded K class or any symbolic algebraic-number package.
    Z=(Q(0),Q(0));ONE=(Q(1),Q(0));ROOTS=(Q(0),Q(1))
    def lift(a):return (Q(a),Q(0))
    def ka(x,y):return (x[0]+y[0],x[1]+y[1])
    def kn(x):return (-x[0],-x[1])
    def ks(x,y):return ka(x,kn(y))
    def km(x,y):return (x[0]*y[0]+q*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
    def ki(x):
        norm=x[0]**2-q*x[1]**2;assert norm
        return (x[0]/norm,-x[1]/norm)
    def kd(x,y):return km(x,ki(y))
    def kc(x):return (x[0],-x[1])
    def ksum(xs):
        out=Z
        for x in xs:out=ka(out,x)
        return out
    def identity(n):return [[ONE if i==j else Z for j in range(n)]for i in range(n)]
    def kmm(A,B):return [[ksum(km(x,y)for x,y in zip(row,col))for col in zip(*B)]for row in A]
    def kmat(A):return [[lift(x)for x in row]for row in A]
    def kparse(A):return [[tuple(map(Q,x))for x in row]for row in A]
    def kinv(A):
        n=len(A);M=[row[:]+e for row,e in zip(A,identity(n))]
        for j in range(n):
            i=next(i for i in range(j,n)if M[i][j]!=Z);M[i],M[j]=M[j],M[i]
            p=ki(M[j][j]);M[j]=[km(x,p)for x in M[j]]
            for i in range(n):
                if i!=j:
                    c=M[i][j];M[i]=[ks(x,km(c,y))for x,y in zip(M[i],M[j])]
        assert [row[:n]for row in M]==identity(n)
        return [row[n:]for row in M]
    def kdet(A):
        A=[row[:]for row in A];d=ONE
        for j in range(len(A)):
            i=next((i for i in range(j,len(A))if A[i][j]!=Z),None)
            if i is None:return Z
            if i!=j:A[i],A[j]=A[j],A[i];d=kn(d)
            p=A[j][j];d=km(d,p)
            for i in range(j+1,len(A)):
                c=kd(A[i][j],p)
                for k in range(j+1,len(A)):A[i][k]=ks(A[i][k],km(c,A[j][k]))
        return d
    gamma=list(map(mat,spin['gamma_matrices']))
    N=[[sum(normal[e]*gamma[e][i][j]for e in range(5))for j in range(4)]for i in range(4)]
    assert N==mat(ac['Clifford_normal4x4'])
    assert mm(N,N)==[[q if i==j else Q(0)for j in range(4)]for i in range(4)]
    I=identity(4);NK=kmat(N);omega=kmat(Ji)
    projectors=[[[kd(ka(I[i][j],km(lift(sign),kd(NK[i][j],ROOTS))),lift(2))for j in range(4)]for i in range(4)]for sign in [1,-1]]
    assert projectors==list(map(kparse,ac['projectors']))
    def pairing(a,b):return ksum(km(km(a[i],omega[i][j]),b[j])for i in range(4)for j in range(4))
    vectors=[];column_pairs=[];pairings=[]
    for projector in projectors:
        columns=tp(projector);i=next(i for i,v in enumerate(columns)if any(x!=Z for x in v))
        j=next(j for j,v in enumerate(columns)if pairing(columns[i],v)!=Z)
        p=pairing(columns[i],columns[j]);vectors.extend([columns[i],[kn(kd(x,p))for x in columns[j]]])
        column_pairs.append([i,j]);pairings.append(p)
    frame=tp(vectors);frame_inverse=kinv(frame);g=tp(frame_inverse)
    assert column_pairs==ac['eigenvector_projector_column_pairs']
    assert pairings==[tuple(map(Q,x))for x in ac['symplectic_pairings_before_normalization']]
    assert frame==kparse(ac['eigenframe4x4'])and g==kparse(ac['physical_pullback_g4x4'])
    assert kmm(kmm(tp(frame),omega),frame)==omega
    assert kmm(kmm(tp(g),kmat(J)),g)==kmat(J)and kdet(g)==ONE
    diagonal=[[km(ROOTS,lift((1 if i<2 else -1)if i==j else 0))for j in range(4)]for i in range(4)]
    assert kmm(kmm(frame_inverse,NK),frame)==diagonal
    # Direct tensor contraction, independently of the author's form-matrix loop.
    transformed=[]
    for (i,j),k in coords:
        transformed.append(ksum(km(km(g[l][k],lift(beta[coords.index(((a,b),l))])),
                                   ks(km(g[a][i],g[b][j]),km(g[b][i],g[a][j])))
                                  for l in range(4)for a,b in pairs))
    assert transformed==[tuple(map(Q,x))for x in ac['transformed_beta24']]
    lookup={coord:x for coord,x in zip(coords,transformed)}
    for k in range(4):assert ka(lookup[((0,1),k)],lookup[((2,3),k)])==Z
    for i,j,k in itertools.combinations(range(4),3):
        assert ka(ks(lookup[((j,k),i)],lookup[((i,k),j)]),lookup[((i,j),k)])==Z
    newpsi=[[lookup[(ij,k)]for k in range(4)]for ij in pairs[:-1]]
    assert newpsi==kparse(ac['transformed_psi5x4'])and newpsi[0]==[Z]*4
    assert kdet(newpsi[1:])==ROOTS
    # Extract the12 biform coefficients through a left inverse of the independently
    # accepted physical24x12 embedding, rather than copying the author's formulas.
    E=mat(read(slice_cert)['physical_embedding24x12'])
    augmented=[row[:]for row in tp(E)];chosen=[];r=0
    for j in range(24):
        i=next((i for i in range(r,12)if augmented[i][j]),None)
        if i is None:continue
        augmented[i],augmented[r]=augmented[r],augmented[i];p=augmented[r][j]
        augmented[r]=[x/p for x in augmented[r]]
        for i in range(r+1,12):
            c=augmented[i][j];augmented[i]=[x-c*y for x,y in zip(augmented[i],augmented[r])]
        chosen.append(j);r+=1
        if r==12:break
    assert len(chosen)==12
    left=inv([E[i]for i in chosen])
    uv=[ksum(km(lift(c),transformed[j])for c,j in zip(row,chosen))for row in left]
    assert all(ksum(km(lift(c),x)for c,x in zip(row,uv))==target for row,target in zip(E,transformed))
    assert uv==[tuple(map(Q,x))for x in ac['biform12']]
    def conv(A,B):
        C=[Z]*(len(A)+len(B)-1)
        for i,a in enumerate(A):
            for j,b in enumerate(B):C[i+j]=ka(C[i+j],km(a,b))
        return C
    v0,v1=uv[6:9],uv[9:12]
    F=conv(uv[:2],conv(v1,v1));middle=conv(uv[2:4],conv(v0,v1));last=conv(uv[4:6],conv(v0,v0))
    F=[ka(ks(x,y),z)for x,y,z in zip(F,middle,last)]
    assert F==[tuple(map(Q,x))for x in ac['binary_quintic_coefficients']]
    syl=[[Z]*4 for i in range(4)]
    for row,poly,shift in [(0,v0,0),(1,v0,1),(2,v1,0),(3,v1,1)]:
        for i,x in enumerate(poly):syl[row][i+shift]=x
    resultant=kdet(syl);assert resultant!=Z and resultant==tuple(map(Q,ac['V_quadratic_resultant']))
    syl5=[[Z]*9 for i in range(9)]
    for shift in range(4):
        for i,x in enumerate(F):syl5[shift][i+shift]=x
    for shift in range(5):
        for i in range(5):syl5[4+shift][i+shift]=km(lift(i+1),F[i+1])
    assert F[5]!=Z
    disc=kd(kdet(syl5),F[5]);D=kd(disc,km(resultant,resultant))
    assert disc==tuple(map(Q,ac['binary_quintic_discriminant']))and D==tuple(map(Q,ac['mixed_discriminant']))
    assert D[1]==0 and D[0]!=0
    circuit=read(circuitpath)
    def evaluate(inputs):
        out=[]
        for i,(op,val)in enumerate(circuit['nodes']):
            if op=='q':x=lift(Q(val))
            elif op=='v':x=inputs[val]
            elif op=='+':x=ksum(out[j]for j in val)
            elif op=='*':
                x=ONE
                for j in val:x=km(x,out[j])
            elif op=='inv':x=ki(out[val])
            else:raise AssertionError(op)
            out.append(x)
        return out
    nodes=evaluate(uv);nodes_conjugate=evaluate([kc(x)for x in uv]);assert len(nodes)==170
    assert nodes==[tuple(map(Q,x))for x in ac['all170_circuit_nodes']]
    assert nodes_conjugate==list(map(kc,nodes))
    o=circuit['outputs'];assert [nodes[i]for i in o['binary_quintic_coefficients']]==F
    assert nodes[o['resultant']]==resultant and nodes[o['quintic_discriminant']]==disc and nodes[o['mixed_discriminant']]==D
    assert D[0].numerator.bit_length()==420 and D[0].denominator.bit_length()==319
    evidence={
        'fixed_sample_index':0,'fixed_source_point':point,'fresh_DAG_node_count':len(payload['nodes']),
        'fresh_all61_coordinates':[str(x)for x in coordinates],
        'fresh_primitive_beta24':[str(x)for x in beta],
        'normal_square_q':str(q),'p0_detC':str(p0),
        'eigenframe_projector_columns':column_pairs,'all_frame_and_full24_pullback_coefficients_match':True,
        'biform_extraction_left_inverse_rows':chosen,
        'independent_quadratic_field_pair_arithmetic':True,
        'independent_quadratic_and_quintic_Sylvester_determinants':True,
        'mixed_D_QQ_numerator':str(D[0].numerator),'mixed_D_QQ_denominator':str(D[0].denominator),
        'mixed_D_nonzero_QQ_not_only_finite_field':True,
        'all170_nodes_and_involution_verified':True,
        'normalizer_frame_independence_uses_prior_accepted_evenness':True,
        'no_additional_source_point_tested':True,
    }
    (output/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,author,acpath,samplepath,samplereview,mappath,mapreview,primitivepath,origpath,spinpath,construction,slice_cert,circuitpath,
           ROOT/'code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py',
           ROOT/'research/notes/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63.md']
    report={
        'id':ID,'all_passed':True,'fixed_existing_NEW25_sample_index':0,'number_source_points_tested':1,
        'fresh_DAG_and_independent_S61_to_beta_reconstruction':True,
        'own_QQ_s_pair_arithmetic_projectors_frame_and_all24_pullback':True,
        'quadratic_resultant_and_quintic_discriminant_by_independent_Sylvester_determinants':True,
        'canonical_mixed_D_exact_QQ_nonzero':True,
        'mixed_D_numerator_bits':420,'mixed_D_denominator_bits':319,
        'all170_circuit_nodes_and_field_conjugation_match':True,
        'mixed_D_and_any_p0_power_cleared_numerator_rejected_as_NEW25_equation':True,
        'mixed_discriminant_construction_unchanged':True,
        'pole_order_k_new_degree_lower_bound_or_NEW25_equation_claim':False,
        'no_new_point_or_high_degree_probe':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},
        'checked_dependency_sha256':checked,
        'evidence_sha256':{'independent-certificate.json':sha(output/'independent-certificate.json')},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({'all_passed':True,'D_exact_QQ_nonzero':True,'seconds':report['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
