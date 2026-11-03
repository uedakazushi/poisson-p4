#!/usr/bin/env python3
"""Universal parallel-Lagrangian first-jet singularity and six actual exclusions.

No generic sample vanishing is promoted to an identity. The universal jet has
eight free primitive coefficients, and its Pfaffian curve has an explicit
singular point. Every full symbolic source family is checked against the
geometric hypothesis before closure/dominance is applied.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from s61_model import X,S61
from q51_universal_section_20261101 import affine
from additional_rank4_families_20261103_v62 import families
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
from primitive_jet_elliptic_curve_20261101 import primitive_beta

ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(itertools.combinations(range(4),2))
NAMES=['R4_ZERO_LOG_5H','R4_MINIMAL_M','R4_A0','R4_DPM','R4_CPLUS_EXT','R4_CMINUS_MAG']
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def matrix(M):return [[str(v)for v in row]for row in M.tolist()]
def pf4(M):return s.expand(M[0,1]*M[2,3]-M[0,2]*M[1,3]+M[0,3]*M[1,2])
def skew(entries):
    M=s.zeros(4)
    for(i,j),v in entries.items():M[i,j]=v;M[j,i]=-v
    return M
def curve(beta,z):
    pp=s.symbols('p0:5');aa,bb,cc,dd,ee=pp
    Pl=skew({(0,1):aa,(0,2):bb,(0,3):cc,(1,2):dd,(1,3):ee,(2,3):-aa})
    Star=skew({(0,1):-aa,(0,2):-ee,(0,3):dd,(1,2):cc,(1,3):-bb,(2,3):aa})
    cov=s.Matrix([sum(s.diff(beta[i,j],v)*Pl[i,j]for i,j in PAIRS)for v in z])
    Phi=Star.row_join(cov).col_join((-cov.T).row_join(s.zeros(1)))
    F=s.Matrix([(-1)**i*pf4(Phi.minor_submatrix(i,i))for i in range(5)])
    assert (F-(Pl*cov).col_join(s.Matrix([pf4(Pl)]))).applyfunc(s.expand)==s.zeros(5,1)
    return pp,Pl,cov,Phi,F

def build(out):
    started=time.monotonic();z=s.Matrix(s.symbols('z0:4'));xx,yy,zz,ww=z
    J=skew({(0,1):s.Integer(1),(2,3):s.Integer(1)})
    # Full compatible linear jet with zero base-base component and no fiber
    # variables in its four base-fiber components: ten coefficients remain.
    A,B,C,D,E,F,G,H,K,L=s.symbols('A B C D E F G H K L')
    raw=skew({(0,1):A*xx+B*zz,(0,3):C*xx+D*zz,(1,2):E*xx+F*zz,
       (1,3):G*xx+H*zz+(B-C)*yy+(-K-F)*ww,(2,3):K*xx+L*zz})
    compatibility=[]
    for i,j,k in itertools.combinations(range(4),3):
        value=s.expand(sum(J[i,b]*s.diff(raw[j,k],z[b])+J[j,b]*s.diff(raw[k,i],z[b])+J[k,b]*s.diff(raw[i,j],z[b])for b in range(4)))
        assert value==0;compatibility.append(str(value))
    alpha=s.Matrix([s.diff(raw[0,1]+raw[2,3],v)/5 for v in z]);assert alpha[1]==alpha[3]==0
    v=J*alpha
    primitive=(raw-2*(alpha.dot(z))*J+z*v.T-v*z.T).applyfunc(s.expand)
    assert s.expand(primitive[0,1]+primitive[2,3])==0
    a,b,c,d,e,f,g,hh=s.symbols('a b c d e f g h')
    universal=skew({(0,1):a*xx+b*zz,(0,3):c*xx+d*zz,(1,2):e*xx+f*zz,
       (1,3):g*xx+hh*zz+(b-c)*yy+(a-f)*ww,(2,3):-a*xx-b*zz})
    substitutions={a:A-3*alpha[0],b:B-2*alpha[2],c:C-alpha[2],d:D,e:E,f:F+alpha[0],g:G,hh:H}
    assert (primitive-universal.subs(substitutions,simultaneous=True)).applyfunc(s.expand)==s.zeros(4)
    beta=(-J*universal*J).applyfunc(s.expand)
    assert beta[1,3]==0
    for i,j,k in itertools.combinations(range(4),3):
        assert s.expand(s.diff(beta[j,k],z[i])-s.diff(beta[i,k],z[j])+s.diff(beta[i,j],z[k]))==0
    pp,Pl,cov,Phi,eqs=curve(beta,z);singular=dict(zip(pp,[0,0,0,0,1]))
    assert cov.subs(singular)==s.zeros(4,1)and eqs.subs(singular)==s.zeros(5,1)
    jac=eqs.jacobian(pp).subs(singular).applyfunc(s.expand)
    assert all(jac[:,i]==s.zeros(5,1)for i in [0,2,3,4])and jac[4,1]==-1
    # Thus rank is exactly one over QQ[a,...,h], including every parameter
    # specialization. On the smooth Pluecker Q3 all curve differentials vanish.

    # Every rank-four constant value with zero base-base block admits a value
    # normalization preserving the parallel vertical plane; no point sample.
    t=s.symbols('t0:5');A2=s.Matrix([[t[0],t[1]],[t[2],t[3]]]);B2=s.Matrix([[0,t[4]],[-t[4],0]])
    C0=s.zeros(2).row_join(A2).col_join((-A2.T).row_join(B2));inv=A2.inv().T
    shear=-inv*B2*inv.T/2
    normal=s.eye(2).row_join(s.zeros(2)).col_join(shear.row_join(inv))
    standard=s.zeros(2).row_join(s.eye(2)).col_join((-s.eye(2)).row_join(s.zeros(2)))
    assert (normal*C0*normal.T-standard).applyfunc(s.cancel)==s.zeros(4)
    assert s.factor(C0.det())==s.factor(A2.det()**2)
    assert s.cancel(normal.det()-1/A2.det())==0

    records={}
    for record in families():
        name=record['component_id']
        if name not in NAMES:continue
        P=affine(record['P'],z)
        assert P[0,2]==0,name
        derivatives=[]
        for i,j in [(0,1),(0,3),(2,1),(2,3)]:
            for fiber in [yy,ww]:
                value=s.expand(s.diff(P[i,j],fiber));assert value==0,(name,i,j,fiber);derivatives.append(str(value))
        records[name]=dict(full_symbolic_affine_P=matrix(P),full_source_parameters=list(map(str,record['parameters'])),
          coefficient_relation=str(record['field_relation'])if record['field_relation']is not None else None,
          zero_base_base_component=True,all_eight_base_fiber_derivatives_zero=derivatives,
          actual_component_identification=record['source'],anchor=record['source_anchor'])
    assert set(records)==set(NAMES)

    ray_path=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz'
    ray=json.loads(gzip.decompress(ray_path.read_bytes()));eta=s.Symbol('eta')
    # Independently write the primitive beta at eta=1, already reconstructed
    # from physical150 and all61 coordinates in the prerequisite packet.
    Bplus=skew({(0,2):s.Rational(27,8)*zz,(0,3):s.Rational(3,2)*xx,(1,3):-s.Rational(2048,27)*yy})
    zold,uold,betold=primitive_beta(S61());us=list(map(s.sympify,ray['S61']))
    actual=betold.subs({v:us[int(str(v)[1:])-1].subs(eta,1)for v in uold})
    assert actual.subs(dict(zip(zold,z)),simultaneous=True)==Bplus
    mons=list(itertools.combinations_with_replacement(range(5),2));physical_pairs=list(itertools.combinations(range(5),2))
    for value,((i,j),(k,l))in zip(map(s.sympify,ray['raw150']),itertools.product(physical_pairs,mons)):
        if value:
            assert i+j-k-l==( -1 if value.has(eta)else 5)
    beta24=[s.diff(Bplus[i,j],v)for i,j in PAIRS for v in z]
    c4,c6,delta,j=curve_invariants_from_beta(beta24)
    assert c4==0 and c6==-949978046398464 and delta==-522255954073519803473068032 and j==0
    assert delta!=0
    out.mkdir(parents=True,exist_ok=False)
    certificate=dict(universal_primitive_linear_jet=matrix(universal),universal_primitive_beta=matrix(beta),
      raw_compatible_ten_parameter_jet=matrix(raw),primitive_projection_alpha=list(map(str,alpha)),eight_parameter_substitution={str(k):str(v)for k,v in substitutions.items()},
      plucker_variables=list(map(str,pp)),plucker_matrix=matrix(Pl),canonical_covector=list(map(str,cov)),canonical_Phi=matrix(Phi),five_Pfaffians=list(map(str,eqs)),
      universal_singular_point=[0,0,0,0,1],entire_Jacobian_at_singular_point=matrix(jac),
      general_constant_value=matrix(C0),vertical_plane_preserving_normalizer=matrix(normal),families=records,
      positive_ray_primitive_beta_at_eta1=matrix(Bplus),positive_ray_c4=str(c4),positive_ray_c6=str(c6),positive_ray_Delta=str(delta),positive_ray_j=str(j))
    archive=out/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/'research/notes/Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63.md',ray_path,
       ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/report.json',
       ROOT/'code/additional_rank4_families_20261103_v62.py',ROOT/'code/additional_rank4_s61_circuits_20261103_v62.py',
       ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json',ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
       ROOT/'code/elliptic_curve_invariants_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',
       ROOT/'code/q51_universal_section_20261101.py',ROOT/'code/s61_model.py']
    inputs += [ROOT/record['actual_component_identification']for record in records.values()]
    result=dict(id='Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63',all_passed=True,field='QQ[ten raw jet coefficients], QQ[eight primitive coefficients], full source parameter rings',
      full_universal_compatible_to_primitive_jet_identity=True,
      full_universal_five_Pfaffian_singular_point_certificate=True,Jacobian_rank_at_singular_point=1,
      rational_value_normalizer_preserves_parallel_lagrangian_fibers=True,normalizer_open='det(A)!=0 iff full value rank4',
      all_six_full_symbolic_source_builders_satisfy_hypotheses=True,
      six_actual_component_firstjet_image_closures_satisfy_Delta_identically_zero=NAMES,
      discriminant_identity_proof='Universal explicit singular point plus characteristic-zero Pfaffian genus-one discriminant criterion; not interpolation, finite samples, or an expanded Delta composition.',
      positive_axis_eta1_c4=str(c4),positive_axis_eta1_c6=str(c6),positive_axis_eta1_Delta=str(delta),
      all_nonzero_positive_axis_points_outside_six_actual_components=True,
      transport_to_entire_nonzero_axis='The global37 positive torus gives one Cstar orbit; intrinsic singularity is invariant under projective frame and nonzero scalar.',
      MPS_T02_through_T38_universal_Delta_zero_claimed=False,new_component_or_known_count_increase_claimed=False,
      input_sha256={str(path.relative_to(ROOT)):sha(path)for path in dict.fromkeys(inputs)},certificate_sha256=sha(archive),seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('UNIVERSAL PARALLEL LAGRANGIAN JET SINGULARITY + SIX ACTUAL EXCLUSIONS PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args().output)
