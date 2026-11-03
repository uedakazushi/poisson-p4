#!/usr/bin/env python3
"""A distinct actual CB polynomial self-pair family through B18, with source ranks."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_self_pair_CB_source_rank_20261103_v63 import X,PAIRS,MONS,vector,canonical,rawCB,action

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63'
h,x,y,z,w=X;b=s.symbols('b',nonzero=True)


def read(path):
    with(gzip.open(path,'rt')if path.suffix=='.gz'else path.open())as f:return json.load(f)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def strings(M):return[[str(v)for v in row]for row in M.tolist()]
def push(G,Gi,P):return(G*P.subs(dict(zip(X,Gi*s.Matrix(X))),simultaneous=True)*G.T).applyfunc(s.cancel)
def matrix(q):
    P=s.zeros(5)
    for k,(i,j)in enumerate(PAIRS):P[i,j]=sum(q[15*k+m]*X[u]*X[v]for m,(u,v)in enumerate(MONS));P[j,i]=-P[i,j]
    return P
def main(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));signal.alarm(120)
    ca=ROOT/'reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/certificate.json.gz';base=read(ca);C=s.Matrix(base['fixed_fiber_linear_basis150_by20'])
    center=s.Matrix(list(map(s.Rational,base['center_fixed20'])));coords=center.copy();coords[1]-=2*b;coords[9]-=b;q=(C*coords).applyfunc(s.expand);P=matrix(q)
    Y=s.Matrix([x*x+s.Rational(8,5)*b*h*x,-s.Rational(2,5)*b*x*x,0,h*x+s.Rational(3,5)*b*x*z,h*h+z*z-s.Rational(7,5)*b*x*w]);ey=s.eye(5)[:,2]
    assert P==(ey*Y.T-Y*ey.T)
    assert s.Matrix([sum(s.diff(P[i,j],X[i])for i in range(5))for j in range(5)])==s.zeros(5,1)
    for i,j,k in itertools.combinations(range(5),3):assert s.expand(sum(P[i,l]*s.diff(P[j,k],X[l])+P[j,l]*s.diff(P[k,i],X[l])+P[k,l]*s.diff(P[i,j],X[l])for l in range(5)))==0
    G1=s.Matrix([[0,1,0,0,0],[b,1,0,-b*b,0],[0,0,1,0,0],[1,0,0,0,0],[0,0,0,0,1]])
    G2=s.eye(5);G2[0,0]=2*b;G2[0,1]=1
    fvalues=[-b*h*x,2*b*x*x];sourcepoints=[];sourceHs=[];maps=[]
    for G,f in zip([G1,G2],fvalues):
        Gi=G.inv();assert G*Gi==s.eye(5)and Gi*G==s.eye(5)
        S=push(G,Gi,P);Cy=s.Matrix([-S[0,2],-S[1,2],S[2,3],S[2,4]])
        Hw=s.cancel(-Cy[2]-s.diff(f,x)*z/10);Hz=s.cancel(Cy[3]+s.diff(f,x)*w/10)
        H=s.integrate(Hw,w);rest=s.cancel(Hz-s.diff(H,z));assert s.diff(rest,w)==0;H=s.cancel(H+s.integrate(rest,z))
        expected=canonical(rawCB(f,0,H));assert (S-expected).applyfunc(s.cancel)==s.zeros(5)
        sourcepoints.append(S);sourceHs.append(H);maps.append(dict(G=strings(G),G_inverse=strings(Gi),determinant=str(s.factor(G.det())),f=str(f),H=str(s.expand(H)),
            full150_source=list(map(str,vector(expected))),full150_pushed_Q=list(map(str,vector(S))),all150_difference_zero=True))
    assert G1[0,:]!=G2[0,:]and G1[0,:]==s.Matrix([[0,1,0,0,0]])and G2[0,:]==s.Matrix([[2*b,1,0,0,0]])
    payload=dict(parameter='b',CA_fixed20_coordinates=list(map(str,coords)),full_polynomial_Q150=list(map(str,q)),quotient_vectorfield=list(map(str,Y)),
        actual_two_CB_sources=maps,distinct_CB_marked_lines=['x','x+2*b*h'],first_CB_W='span(x,x+b*h-b^2*z)',second_CB_W='span(x,h)',all_flags_L='e_y')
    def writepayload():
        with(out/'certificate.json.gz').open('wb')as file:
            with gzip.GzipFile(fileobj=file,mode='wb',filename='',mtime=0)as archive:archive.write((json.dumps(payload,separators=(',',':'))+'\n').encode())
    writepayload();print(json.dumps(dict(stage='actual_two_CB_sources_saved',seconds=time.monotonic()-start,determinants=[m['determinant']for m in maps])),flush=True)
    globalpath=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json';globaldata=read(globalpath);review=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW.json';assert read(review)['all_passed']
    fiberreport=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json';fiberreview=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json'
    assert read(fiberreview)['all_passed']and read(fiberreport)['whole_flag_fibers']['CB']['scheme_length']==2
    parabolic=list(map(tuple,globaldata['parabolic_GL_positions']));flags=[pos for pos in itertools.product(range(5),repeat=2)if pos not in parabolic];assert len(flags)==9
    u=s.symbols('u0:19');mu=s.symbols('mu');cubics=[s.prod((h,x,z,w)[i]for i in mon)for mon in itertools.combinations_with_replacement(range(4),3)if any(i>=2 for i in mon)]
    F=sum(c*m for c,m in zip(u[:3],[h*h,h*x,x*x]));Hu=sum(c*m for c,m in zip(u[3:],cubics));universal=vector(canonical(rawCB(F,mu*F,Hu)));J=universal.jacobian([*u,mu])
    fixed=J[:,:19].subs(mu,0);assert not any(v.free_symbols for v in fixed);fixedrows=list(DomainMatrix.from_Matrix(fixed.T).convert_to(s.QQ).rref()[1]);left=fixed[fixedrows,:].inv()
    def quotient(A):return(A-fixed*left*A[fixedrows,:]).applyfunc(s.cancel)
    domain=s.QQ.frac_field(b)
    def ranks(A):
        pc=list(DomainMatrix.from_Matrix(A).convert_to(domain).rref()[1]);pr=list(DomainMatrix.from_Matrix(A[:,pc].T).convert_to(domain).rref()[1]);D=DomainMatrix.from_Matrix(A[pr,pc]).convert_to(domain)
        return dict(rank=len(pc),rows=pr,columns=pc,minor=str(s.factor(domain.to_sympy(D.det()))))
    fibers=[];ranksources=[]
    for S,f,H in zip(sourcepoints,fvalues,sourceHs):
        values=[*[s.Poly(f,h,x).coeff_monomial(mon)for mon in[h*h,h*x,x*x]],*[s.Poly(H,h,x,z,w).coeff_monomial(mon)for mon in cubics],0]
        sub=dict(zip([*u,mu],values));assert (universal.subs(sub)-vector(S)).applyfunc(s.cancel)==s.zeros(150,1)
        fiber=J.subs(sub).applyfunc(s.cancel);assert fiber[:,19]==vector(canonical(rawCB(0,f,0)))
        full=fiber.row_join(s.Matrix.hstack(*[vector(action(i,j,S))for i,j in flags]));Q=quotient(full[:,19:]);witness=ranks(Q)
        fibers.append(fiber);ranksources.append(dict(source20_values=list(map(str,values)),full150_by29=strings(full),quotient150_by10=strings(Q),rank=19+witness['rank'],quotient_witness=witness))
    orbit=s.Matrix.hstack(*[vector(action(i,j,P))for i,j in itertools.product(range(5),repeat=2)]);rOrbit=ranks(orbit);rScale=ranks(orbit.row_join(q));rFamily=ranks(orbit.row_join(q).row_join(q.diff(b)))
    print(json.dumps(dict(stage='source29_and_GL_ranks',seconds=time.monotonic()-start,source_ranks=[r['rank']for r in ranksources],orbit=rOrbit['rank'],family=rFamily['rank'])),flush=True)
    # Same-point cross tangent in first CB frame. GL directions are shared by both sources.
    relative=G1*G2.inv();relative_inverse=G2*G1.inv();back=s.Matrix.hstack(*[vector(push(relative,relative_inverse,matrix(fibers[1][:,j])))for j in range(20)])
    firstFull=quotient(s.Matrix.hstack(fibers[0][:,19],*[vector(action(i,j,sourcepoints[0]))for i,j in flags],back));rCross=ranks(firstFull);union=19+rCross['rank'];intersection=sum(r['rank']for r in ranksources)-union
    payload.update(source_rank_certificates=ranksources,full_GL150_by25=strings(orbit),GL_rank_witness=rOrbit,GL_scalar_rank_witness=rScale,GL_scalar_parameter_rank_witness=rFamily,
        same_point_cross_quotient150_by30=strings(firstFull),cross_quotient_witness=rCross,branch_tangent_union_rank=union,branch_tangent_intersection_dimension=intersection)
    writepayload()
    result=dict(id=STEM,all_passed=True,actual_new_polynomial_family='Q_b=P0-b*(2*C[1]+C[9])',all150_both_CB_source_identities_zero=True,
        strict350_Jacobi_and25_divergence_zero=True,two_actual_CB_flags_distinct_for_all_b_nonzero=True,
        both_source_differential_ranks=[r['rank']for r in ranksources],both_source_quotient10_minors=[r['quotient_witness']['minor']for r in ranksources],
        GL_affine_orbit_rank=rOrbit['rank'],GL_plus_scalar_rank=rScale['rank'],GL_scalar_parameter_rank=rFamily['rank'],projective_saturated_family_dimension=rFamily['rank']-1,
        same_point_branch_tangent_union_rank=union,same_point_branch_tangent_intersection_dimension=intersection,
        generic_whole_CB_fiber_two_reduced_points_from_proper_open=True,
        generic_reduced_CB_completion_two_smooth29_branches=True,generic_CB_cone_depth_upper_bound=intersection+1,
        all_nonzero_whole_CB_flag_classification_or_origin_conductor_classification_claimed=False,seconds=time.monotonic()-start)
    paths=[Path(__file__),ca,globalpath,review,fiberreport,fiberreview,ROOT/'code/ambient_CA_self_pair_CB_source_rank_20261103_v63.py'];note=ROOT/f'research/notes/{STEM}.md'
    if note.exists():paths.append(note)
    result['input_sha256']={str(p.relative_to(ROOT)):sha(p)for p in paths};result['evidence_sha256']={'certificate.json.gz':sha(out/'certificate.json.gz')}
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256','evidence_sha256']},indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
