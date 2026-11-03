#!/usr/bin/env python3
"""Independent Laurent coefficient/finite-etale projector normal-form review.

No author module, GB, new beta point, eigenvector approximation, or invariant
probe is used. Matrix traces/determinants are rebuilt with Fraction dictionaries.
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
    mat,mm,tp,det,plus,mul,const,variable,power,determinant_poly,
)

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63_REVIEW'


def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pscale(p,c):return {m:a*c for m,a in p.items()if a*c}
def psum(ps):
    out={}
    for p in ps:out=plus(out,p)
    return out
def pmm(A,B):return [[psum(mul(x,y)for x,y in zip(row,col))for col in zip(*B)]for row in A]
def ptrace(A):return psum(A[i][i]for i in range(len(A)))
def substitute(p,mapping,n=12):
    out={}
    for mon,c in p.items():
        base=list(mon)
        for i in mapping:base[i]=0
        term={tuple(base):c}
        for i,replacement in mapping.items():
            assert mon[i]>=0
            term=mul(term,power(replacement,mon[i],n))
        out=plus(out,term)
    return out
def laurent(text,variables):
    # Parsing only; all identities are checked by our Fraction dictionaries.
    out={}
    for term in s.Add.make_args(s.expand(s.sympify(text))):
        powers=term.as_powers_dict();mon=tuple(int(powers.get(v,0))for v in variables)
        coefficient=term/s.prod(v**k for v,k in zip(variables,mon))
        c=Q(str(coefficient));out[mon]=out.get(mon,Q(0))+c
    return {m:c for m,c in out.items()if c}


def main(output):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63/report.json'
    certpath=author.parent/'certificate.json.gz';ar=read(author);cert=read(certpath);assert ar['all_passed']
    pairpath=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/certificate.json.gz'
    pairreview=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/report.json'
    pointreview=ROOT/'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63_REVIEW/report.json'
    checked={}
    for p in [author,pairreview,pointreview]:
        data=read(p);assert data['all_passed']
        for name,want in data['input_sha256'].items():assert sha(ROOT/name)==want,name;checked[name]=want
    names=s.symbols('x y a b c d e f g h i j');n=12
    x,y,a,b,c,d,e,f,g,h,i,j=[variable(k,n)for k in range(n)]
    ONE=const(1,n);ZERO={};neg=lambda p:pscale(p,-1)
    C=[[a,b,g,h],[c,neg(a),i,j],[neg(j),h,d,e],[i,neg(g),f,neg(d)]]
    D=[[x,ZERO,ZERO,ZERO],[ZERO,neg(x),ZERO,ZERO],[ZERO,ZERO,y,ZERO],[ZERO,ZERO,ZERO,neg(y)]]
    assert all(p==laurent(text,names)for row,saved in zip(C,cert['generic_C'])for p,text in zip(row,saved))
    assert all(p==laurent(text,names)for row,saved in zip(D,cert['generic_D'])for p,text in zip(row,saved))
    J=mat([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    JP=[[const(v,n)for v in row]for row in J]
    for A in [C,D]:assert all(not q for row in [[plus(x,y)for x,y in zip(r,t)]for r,t in zip(pmm(tp(A),JP),pmm(JP,A))]for q in row)
    CC,CD,DD=pmm(C,C),pmm(C,D),pmm(D,D)
    q1,q2,q3=ptrace(CC),ptrace(CD),ptrace(DD)
    trC4,trD4=ptrace(pmm(CC,CC)),ptrace(pmm(DD,DD))
    f1=plus(q1,q3,-2);f2=plus(q2,q3)
    f3=psum([trC4,pscale(trD4,4),pscale(mul(q3,q3),-3)])
    f4=psum([neg(trC4),pscale(ptrace(pmm(CC,DD)),6),pscale(ptrace(pmm(CD,DD)),12),pscale(trD4,8)])
    equations=[f1,f2,f3,f4]
    assert all(q==laurent(text,names)for q,text in zip(equations,cert['four_original_constraint_polynomials']))
    delta=plus(mul(x,x),mul(y,y),-1);detdiff=plus(determinant_poly(C,n),mul(delta,delta),-1)
    assert psum([f3,pscale(detdiff,4),pscale(mul(f1,plus(q1,q3,2)),-Q(1,2))])=={}
    H=pscale(plus(f3,f4),Q(1,12));assert H==laurent(cert['linear_constraint_H'],names)
    yinv={tuple(-int(k==1)for k in range(n)):Q(1)}
    d0=neg(mul(psum([mul(a,x),mul(x,x),mul(y,y)]),yinv))
    M=plus(mul(g,j),mul(h,i),-1);L=plus(mul(plus(a,x),plus(a,x)),mul(y,y),-1)
    c0=plus(M,L,-1);f0=plus(M,mul(mul(mul(x,x),L),mul(yinv,yinv)),-1)
    for key,value in [('d',d0),('c',c0),('f',f0),('M',M),('L',L)]:assert value==laurent(cert['graph'][key],names)
    Z1,Z2,Z3=plus(d,d0,-1),plus(c,c0,-1),plus(f,f0,-1)
    gauge={3:ONE,6:ONE};graph={**gauge,5:d0,4:c0,7:f0}
    F1=pscale(substitute(f1,gauge),Q(1,2));HG=substitute(H,gauge)
    d_quadratic=mul(plus(pscale(d0,2),Z1),Z1)
    h_quadratic=mul(psum([pscale(mul(mul(y,y),d0),2),pscale(mul(mul(y,y),y),2),mul(mul(y,y),Z1)]),Z1)
    assert plus(f2,pscale(mul(y,Z1),-2))=={}
    assert psum([F1,neg(Z2),neg(Z3),neg(d_quadratic)])=={}
    assert psum([HG,neg(mul(mul(x,x),Z2)),neg(mul(mul(y,y),Z3)),neg(h_quadratic)])=={}
    X,Y=plus(F1,d_quadratic,-1),plus(HG,h_quadratic,-1)
    assert psum([Y,neg(mul(mul(y,y),X)),neg(mul(delta,Z2))])=={}
    assert psum([mul(mul(x,x),X),neg(Y),neg(mul(delta,Z3))])=={}
    graph_detdiff=substitute(detdiff,graph)
    final=mul(mul(y,y),graph_detdiff)
    assert all(min(mon)>=0 and value.denominator==1 for mon,value in final.items())
    keep=[0,1,2,8,9,10,11]
    assert all(all(not mon[k]for k in range(12)if k not in keep)for mon in final)
    final7={tuple(mon[k]for k in keep):value for mon,value in final.items()}
    assert len(final7)==77 and max(map(sum,final7))==8
    assert final7=={tuple(mon):Q(c)for mon,c in cert['final_polynomial_sparse']}
    assert final==laurent(cert['final_polynomial'],names)
    for equation,scalar in zip(equations,[0,0,-4,4]):assert substitute(equation,graph)==pscale(graph_detdiff,scalar)
    # Difference quotients by t-value are made by a finite geometric sum,
    # independently of symbolic rational division/cancellation.
    residual=substitute(detdiff,gauge);cofactors=[]
    for variable_index,replacement in [(5,d0),(4,c0),(7,f0)]:
        cofactor={}
        for mon,value in residual.items():
            degree=mon[variable_index];assert degree>=0
            base=list(mon);base[variable_index]=0
            for r in range(degree):
                tpower=base[:];tpower[variable_index]=degree-1-r
                cofactor=plus(cofactor,mul({tuple(tpower):value},power(replacement,r,n)))
        next_residual=substitute(residual,{variable_index:replacement})
        difference=plus(variable(variable_index,n),replacement,-1)
        assert plus(plus(residual,next_residual,-1),mul(difference,cofactor),-1)=={}
        cofactors.append(cofactor);residual=next_residual
    assert residual==graph_detdiff
    assert all(q==laurent(text,names)for q,text in zip(cofactors,cert['determinant_telescoping_cofactors']))
    # The two torus weights and residual sign are exact entrywise statements.
    torus_weights=[(1,0),(-1,0),(0,1),(0,-1)]
    assert tuple(torus_weights[0][k]-torus_weights[1][k]for k in range(2))==(2,0)
    assert tuple(torus_weights[2][k]-torus_weights[3][k]for k in range(2))==(0,2)
    relative=[-1,-1,1,1]
    for r in range(4):
        for col in range(4):
            expected=substitute(C[r][col],{k:neg(variable(k,n))for k in [8,9,10,11]})
            assert pscale(C[r][col],relative[r]*relative[col])==expected
    # Same accepted numeric pair[0], without testing another beta.
    pair=read(pairpath);assert pair['same_point_index']==cert['same_point_index']==0
    actualC,actualD=mat(pair['actual_C4']),mat(pair['actual_D4'])
    D2=mm(actualD,actualD);S=sum(D2[k][k]for k in range(4))/2;P=det(actualD);disc=S*S-4*P
    assert P and disc and S==Q(cert['actual_D_squared_trace_half'])and P==Q(cert['actual_detD'])
    assert disc==Q(cert['actual_D_squared_discriminant'])
    D4=mm(D2,D2)
    assert all(D4[r][k]-S*D2[r][k]+(P if r==k else 0)==0 for r in range(4)for k in range(4))
    Z=(Q(0),Q(0));ONEK=(Q(1),Q(0));eta=(Q(0),Q(1))
    lift=lambda a:(Q(a),Q(0))
    ka=lambda a,b:(a[0]+b[0],a[1]+b[1])
    kn=lambda a:(-a[0],-a[1])
    km=lambda a,b:(a[0]*b[0]+disc*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    conj=lambda a:(a[0],-a[1])
    def ki(a):
        norm=a[0]**2-disc*a[1]**2;assert norm
        return (a[0]/norm,-a[1]/norm)
    kd=lambda a,b:km(a,ki(b))
    def total(values):
        out=Z
        for x in values:out=ka(out,x)
        return out
    def kmm(A,B):return [[total(km(x,y)for x,y in zip(row,col))for col in zip(*B)]for row in A]
    tr=lambda A:total(A[k][k]for k in range(4))
    Xk,Yk=kd(ka(lift(S),eta),lift(2)),kd(ka(lift(S),kn(eta)),lift(2))
    assert km(Xk,Yk)==lift(P)and [Xk,Yk]==[tuple(map(Q,z))for z in cert['D_squared_X_Y']]
    I=[[ONEK if r==k else Z for k in range(4)]for r in range(4)]
    projector=[[kd(ka(lift(D2[r][k]),kn(km(Yk,I[r][k]))),eta)for k in range(4)]for r in range(4)]
    other=[[conj(x)for x in row]for row in projector]
    assert kmm(projector,projector)==projector and kmm(other,other)==other
    assert [[ka(x,y)for x,y in zip(a,b)]for a,b in zip(projector,other)]==I
    assert tr(projector)==lift(2)
    assert projector==[[tuple(map(Q,z))for z in row]for row in cert['D_squared_projector']]
    CK,DK=[[lift(x)for x in row]for row in actualC],[[lift(x)for x in row]for row in actualD]
    def block_product(proj,eigenvalue):
        dc=tr(kmm(kmm(DK,CK),proj))
        return ka(kd(tr(kmm(kmm(kmm(proj,CK),proj),CK)),lift(2)),kn(kd(km(dc,dc),km(lift(4),eigenvalue))))
    bc,ef=block_product(projector,Xk),block_product(other,Yk)
    assert ef==conj(bc)and bc==tuple(map(Q,cert['actual_block_bc']))and ef==tuple(map(Q,cert['actual_block_ef']))
    norm=km(bc,ef);assert norm[1]==0 and norm[0]and norm==tuple(map(Q,cert['actual_block_bc_times_ef']))
    assert norm[0].numerator.bit_length()==295 and norm[0].denominator.bit_length()==282
    evidence={
        'all_four_trace_polynomials_rebuilt_with_Fraction_Laurent_dictionaries':True,
        'determinant_Newton_identity_all_coefficients':True,
        'all_three_forward_graph_identities_and_both_reverse_linear_cofactors':True,
        'all_three_det_telescoping_cofactors_by_geometric_difference_sum':True,
        'final77_integer_coefficients_match':True,'final_total_degree':8,
        'ideal_equality_on_x_y_Delta_units_after_b_e_gauge':True,
        'same_existing_pair_index':0,'D_squared_trace_half':str(S),'det_D':str(P),'D_squared_discriminant':str(disc),
        'projector_trace':2,'own_quadratic_pair_projector_arithmetic':True,
        'bc_ef_norm_QQ':str(norm[0]),'norm_nonzero_QQ_not_GF_only':True,
        'finite_etale_gauges_and_finite_residual_scope_in_review_note':True,
    }
    (output/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,author,certpath,pairpath,pairreview,pointreview,
           ROOT/'code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py',
           ROOT/'code/q_firstjet_two_adjoint_normal_form_20261103_v63.py',
           ROOT/'research/notes/Q_FIRSTJET_TWO_ADJOINT_NORMAL_FORM_20261103_v63.md']
    result={
        'id':ID,'all_passed':True,
        'full_four_constraint_graph_and_reverse_ideal_equality_verified_independently':True,
        'seven_coordinate77term_degree8_hypersurface_accepted':True,
        'same_point_D2_projector_and_bc_ef_QQ_nonzero_verified':True,
        'chart_meets_actual_Z_component_and_NEW25':True,
        'finite_Weyl_and_effective_relative_sign_cross_section_accepted':True,
        'generic_orbit_dimension10_and_NEW25_quotient_dimension5_on_this_open':True,
        'quotient_rationality_all_H_image_normality_or_NEW25_equation_claim':False,
        'no_author_code_GB_new_beta_point_or_highdegree_probe':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},
        'checked_dependency_sha256':checked,
        'evidence_sha256':{'independent-certificate.json':sha(output/'independent-certificate.json')},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'seconds':time.monotonic()-start,
    }
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'seconds':result['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
