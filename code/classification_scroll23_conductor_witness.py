#!/usr/bin/env python3
"""One exact S(2,3)->P3 projection with an absolute conductor certificate.

No GB/point search. A closed-rational-one-form matrix proves absolute
irreducibility: each absolute factor would contribute an independent dlog.
"""
import argparse,hashlib,itertools as it,json,resource,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix

ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    ss,tt,uu,vv,T,r,x,y,z=s.symbols('s t u v T r x y z');A=[ss**2,ss*tt,tt**2,0];B=[ss**3,tt**3,ss**2*tt,ss**3+tt**3]
    minors=[s.expand(A[i]*B[j]-A[j]*B[i])for i,j in it.combinations(range(4),2)];gcd=s.Poly(minors[0],ss,tt,domain=s.QQ)
    for value in minors[1:]:gcd=s.gcd(gcd,s.Poly(value,ss,tt,domain=s.QQ))
    assert gcd.total_degree()==0
    monA=[ss**2,ss*tt,tt**2];monB=[ss**3,ss**2*tt,ss*tt**2,tt**3]
    matrix=s.Matrix([[s.Poly(value,ss,tt).coeff_monomial(mon)for value in A]for mon in monA]+[[s.Poly(value,ss,tt).coeff_monomial(mon)for value in B]for mon in monB]);assert len(DomainMatrix.from_Matrix(matrix).convert_to(s.QQ).rref()[1])==4
    d=1+T**3;Z=[r+1,r*T+T**3,r*T**2+T,d]
    first=s.expand(d*(y-T*x)+T-T**3);second=s.expand(d*(z-T*T*x)-T+T*T);res=s.resultant(first,second,T);factors=s.factor_list(res,x,y,z)[1];candidates=[]
    def homogeneous_pullback(F,derivative=False):
        out=0
        for exponent,coefficient in s.Poly(F,x,y,z).terms():
            if derivative:
                if exponent[2]==0:continue
                coefficient*=exponent[2];exponent=(exponent[0],exponent[1],exponent[2]-1);power=4-sum(exponent)
            else:power=5-sum(exponent)
            out+=coefficient*Z[0]**exponent[0]*Z[1]**exponent[1]*Z[2]**exponent[2]*d**power
        return s.expand(out)
    for factor,multiplicity in factors:
        if s.Poly(factor,x,y,z).total_degree()==5 and homogeneous_pullback(factor)==0:candidates.append(factor)
    assert len(candidates)==1
    F=s.Poly(candidates[0],x,y,z,domain=s.QQ).primitive()[1].as_expr();assert homogeneous_pullback(F)==0
    mapx=Z[0]/d;mapy=Z[1]/d;J=s.factor(s.det(s.Matrix([mapx,mapy]).jacobian((T,r))));Jnum=s.cancel(d**3*J);assert s.denom(Jnum)==1
    numerator=homogeneous_pullback(F,True);quotient,remainder=s.div(s.Poly(numerator,r,T,domain=s.QQ),s.Poly(Jnum,r,T,domain=s.QQ));assert remainder.is_zero
    G=quotient.primitive()[1].as_expr();factor=s.cancel(quotient.as_expr()/G);assert factor.is_Rational and factor!=0
    assert s.expand(numerator-factor*Jnum*G)==0
    polynomial=s.Poly(G,r,T,domain=s.QQ);dr=polynomial.degree(r);dt=polynomial.degree(T);assert dr==3
    assert all(j<=6-i for (i,j),coefficient in polynomial.terms())
    homogenized=s.expand(sum(coefficient*uu**i*vv**(3-i)*ss**(6-i-j)*tt**j for(i,j),coefficient in polynomial.terms()))
    assert homogenized.subs(ss,0)!=0 and homogenized.subs(vv,0)!=0
    squarefree=s.gcd(s.gcd(polynomial,s.Poly(s.diff(G,r),r,T)),s.Poly(s.diff(G,T),r,T));assert squarefree.total_degree()==0
    Path('/tmp/scroll23_conductor_witness_progress.json').write_text(json.dumps(dict(status='projection_and_global_conductor_complete',F=str(F),Jacobian=str(J),G=str(G),global_section=str(homogenized)),indent=2)+'\n')
    # A/f dr+B/f dT closed; rectangular degree bounds include every dlog factor.
    unknowns=[('A',i,j)for i in range(dr)for j in range(dt+1)]+[('B',i,j)for i in range(dr+1)for j in range(dt)]
    columns=[]
    for kind,i,j in unknowns:
        mon=r**i*T**j
        value=G*s.diff(mon,T)-mon*s.diff(G,T)if kind=='A'else-G*s.diff(mon,r)+mon*s.diff(G,r)
        columns.append(dict(s.Poly(s.expand(value),r,T,domain=s.QQ).terms()))
    exponents=sorted(set().union(*(set(value)for value in columns)),reverse=True);M=s.Matrix([[col.get(e,0)for col in columns]for e in exponents]);DM=DomainMatrix.from_Matrix(M).convert_to(s.QQ);RR,pivots=DM.rref();pivots=list(pivots);rank=len(pivots);assert rank==len(unknowns)-1
    free=next(j for j in range(len(unknowns))if j not in pivots);K=s.zeros(len(unknowns),1);K[free]=1;R=RR.to_Matrix()
    for row,pivot in enumerate(pivots):K[pivot]=-R[row,free]
    assert M*K==s.zeros(M.rows,1)
    dlog=s.Matrix([s.Poly(s.diff(G,r)if kind=='A'else s.diff(G,T),r,T).coeff_monomial(r**i*T**j)for kind,i,j in unknowns]);assert dlog!=s.zeros(len(unknowns),1)and M*dlog==s.zeros(M.rows,1)
    rows=list(DomainMatrix.from_Matrix(M.T).convert_to(s.QQ).rref()[1]);det=M[rows,pivots].det(method='domain-ge');assert det!=0
    inputs=[Path(__file__),ROOT/'research/classification/RANK_TWO_SEXTIC_SURFACE_NONCONICAL_DIVISOR_EXCLUSION.md',ROOT/'research/classification/RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md']
    result=dict(task_id='G',all_passed=True,one_explicit_projection_only=True,section_A2=list(map(str,A)),section_B3=list(map(str,B)),all_six_fiber_independence_minors=list(map(str,minors)),binary_gcd=str(gcd.as_expr()),all_four_projection_sections_independent=True,
        finite_birational_normalization_justification='Basepointfree subsystem of the ample scrollH is finite. Image is not a plane because4sections are independent. Prime degreeH^2=5 therefore gives imagequintic and degree1 map, so the smoothscroll is its normalization.',
        affine_parameter_map_Z=list(map(str,Z)),implicit_resultant_factors=[dict(factor=str(a),multiplicity=b)for a,b in factors],exact_quintic_equation=str(F),quintic_pullback_zero=True,
        Jacobian_xy=str(J),Jacobian_numerator=str(Jnum),conductor_affine_polynomial=str(G),conductor_global_section_3H_minus3f=str(homogenized),canonical_frame_scalar=str(factor),global_conductor_identity='G=(1+T^3)*F_z(map)/Jac_xy up to the recorded nonzero scalar.',full_homogeneous_canonical_identity_verified=True,no_s_infinity_or_v_zero_component=True,affine_conductor_squarefree=True,
        closed_form_absolute_factor_certificate=dict(bidegree=[dr,dt],unknown_order=unknowns,all_spatial_row_exponents=exponents,full_matrix=[[str(v)for v in row]for row in M.tolist()],rank=rank,kernel_dimension=1,full_kernel=list(map(str,K)),dlog_f_vector=list(map(str,dlog)),nonzero_rank_minor=dict(rows=rows,columns=pivots,determinant=str(det))),
        absolute_irreducibility_logic='OverC every distinct absolute factor f_i contributes (f/f_i)*df_i in these rectangular degree bounds. Their rationalforms dlog(f_i) are independent by residues along f_i. The exactQQmatrix still has kernel dimension1 after extension toC, so the squarefree f has one absolute factor. No missing global boundary component then gives an integral whole conductor.',
        whole_scroll_conductor_absolutely_irreducible=True,scope='One finite projection witness for properness of the split-conductor locus. Full rank2S6 conical G dimension claim additionally needs the separate geometric conductor-family/residue/dimension proof; no Poisson tensor is sampled or counted.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in inputs},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as f:json.dump(result,f,indent=2,default=str);f.write('\n')
    print(json.dumps(dict(all_passed=True,conductor_bidegree=[dr,dt],closed_form_rank=rank,unknowns=len(unknowns),seconds=result['seconds'])))


if __name__=='__main__':main()
