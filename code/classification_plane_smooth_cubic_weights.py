#!/usr/bin/env python3
"""Exact full normal-weight arithmetic for a plane with smooth cubic restriction.

All ten cubic coefficients, all30 weight1 and all9 weight2 coefficients are
symbolic. Geometry (connection rigidity and smooth-cubic injectivity) is
separate; this checks every Schouten/Euler coefficient used in the proof.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import hashlib,json,resource,signal,time
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
PAIRS=list(combinations(range(5),2));TRIPLES=list(combinations(range(5),3));QUADS=list(combinations(range(5),4))


def monomials(n,d):
    result=[]
    for terms in combinations_with_replacement(range(n),d):
        e=[0]*n
        for i in terms:e[i]+=1
        result.append(tuple(e))
    return sorted(result,reverse=True)


def wedge(first,second):return first*second.T-second*first.T


def bilinear(P,Q,x):
    return {tri:s.expand(sum(P[a,l]*s.diff(Q[b,c],x[l])+Q[a,l]*s.diff(P[b,c],x[l])for a,b,c in [tri,(tri[1],tri[2],tri[0]),(tri[2],tri[0],tri[1])]for l in range(5)))for tri in TRIPLES}


def projective(strict,x):
    return {quad:s.expand(sum((-1)**i*x[quad[i]]*strict[quad[:i]+quad[i+1:]]for i in range(4)))for quad in QUADS}


def all350(jac,x):
    return [str(s.Poly(jac[quad],*x).coeff_monomial(e))for quad in QUADS for e in monomials(5,4)]


def main():
    destination=ROOT/'reports/CLASSIFICATION_PLANE_SMOOTH_CUBIC_NORMAL_WEIGHTS_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    x=s.Matrix(s.symbols('x0 x1 x2 z w'));base=s.Matrix([*x[:3],0,0]);z,w=x[3:];a=s.Symbol('a');coeff=s.symbols('f0:10')
    f=sum(c*s.prod(x[i]**e[i]for i in range(3))for c,e in zip(coeff,monomials(3,3)))
    def gradient(g):
        P=s.zeros(5)
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:P[i,j]=s.diff(g,x[k]);P[j,i]=-P[i,j]
        return P
    Pf=gradient(f)
    def linear(name):
        symbols=s.symbols(name+'0:9');P=s.zeros(5)
        for k,(i,j)in enumerate(combinations(range(3),2)):
            P[i,j]=sum(symbols[3*k+l]*x[l]for l in range(3));P[j,i]=-P[i,j]
        return P,symbols
    L,lp=linear('l');M,mp=linear('m');vectors={name:s.Matrix([*s.symbols(name.lower()+'0:3'),0,0])for name in ['B','C','D','E']}
    B,C,D,E=[vectors[name]for name in ['B','C','D','E']];ez=s.Matrix([0,0,0,1,0]);ew=s.Matrix([0,0,0,0,1])
    P1=z*L+w*M+z*z*wedge(D,ew)+z*w*(wedge(B,ez)+wedge(E,ew))+w*w*wedge(C,ez)
    def U(V):return Pf.applyfunc(lambda value:-sum(V[i]*s.diff(value,x[i])for i in range(3)))
    Uv={name:U(V)for name,V in vectors.items()}
    # Smoothness supplies injectivity of V->V(f). The identity is polynomial for every cubic.
    for name,V in vectors.items():
        assert s.expand(x[0]*Uv[name][1,2]-x[1]*Uv[name][0,2]+x[2]*Uv[name][0,1]+2*sum(V[i]*s.diff(f,x[i])for i in range(3)))==0
    base_strict_L=bilinear(Pf,L,x)[0,1,2];base_strict_M=bilinear(Pf,M,x)[0,1,2]
    # Curl criterion has no numerical cubic specialization.
    ell=s.Matrix([L[1,2],-L[0,2],L[0,1]]);curl=s.Matrix([s.diff(ell[2],x[1])-s.diff(ell[1],x[2]),s.diff(ell[0],x[2])-s.diff(ell[2],x[0]),s.diff(ell[1],x[0])-s.diff(ell[0],x[1])])
    curl_dot=s.expand(sum(curl[i]*s.diff(f,x[i])for i in range(3)))
    assert s.expand(base_strict_L+curl_dot)==0
    q2=s.symbols('q0:6');g=sum(c*s.prod(x[i]**e[i]for i in range(3))for c,e in zip(q2,monomials(3,2)))
    h2=s.symbols('h0:6');h=sum(c*s.prod(x[i]**e[i]for i in range(3))for c,e in zip(h2,monomials(3,2)))
    records=[]
    for name,K in [('two_roots',a*z*w),('double_root',a*z*z),('zero',s.Integer(0))]:
        P0=Pf+K*wedge(ez,ew);raw=bilinear(P0,P1,x);jac=projective(raw,x)
        if name=='two_roots':normal=a*z**3*D+a*w**3*C
        elif name=='double_root':normal=-a*z**3*E+a*z*z*w*B+2*a*z*w*w*C
        else:normal=s.zeros(5,1)
        Tz=K*M+z*w*Uv['B']+w*w*Uv['C'];Tw=-K*L+z*z*Uv['D']+z*w*Uv['E']
        expected={tri:s.Integer(0)for tri in TRIPLES};expected[0,1,2]=z*base_strict_L+w*base_strict_M
        for i,j in combinations(range(3),2):expected[i,j,3]=Tz[i,j];expected[i,j,4]=Tw[i,j]
        for i in range(3):expected[i,3,4]=normal[i]
        assert all(s.expand(raw[tri]-expected[tri])==0 for tri in TRIPLES)
        forcing=wedge(base,normal)+z*Tw-w*Tz
        for i,j in combinations(range(3),2):assert s.expand(jac[i,j,3,4]-forcing[i,j])==0
        if name=='two_roots':closed_L=Uv['E']/a;closed_M=-Uv['B']/a;closed_B=B;closed_C=s.zeros(5,1);closed_D=s.zeros(5,1);closed_E=E
        elif name=='double_root':closed_L=Uv['D']/a-wedge(base,E);closed_M=Uv['E']/a;closed_B=s.zeros(5,1);closed_C=s.zeros(5,1);closed_D=D;closed_E=E
        else:closed_L=gradient(g);closed_M=gradient(h);closed_B=closed_C=closed_D=closed_E=s.zeros(5,1)
        kernel=z*closed_L+w*closed_M+z*z*wedge(closed_D,ew)+z*w*(wedge(closed_B,ez)+wedge(closed_E,ew))+w*w*wedge(closed_C,ez)
        assert all(s.expand(value)==0 for value in projective(bilinear(P0,kernel,x),x).values())
        if K!=0:
            Y=s.Matrix([*s.symbols('y0:3'),0,0]);Z=s.Matrix([*s.symbols('v0:3'),0,0]);field=z*Y+w*Z
            orbit=s.zeros(5)
            for i,j in PAIRS:
                orbit[i,j]=s.expand(-sum(field[k]*s.diff(P0[i,j],x[k])for k in range(5))+sum(P0[k,j]*s.diff(field[i],x[k])+P0[i,k]*s.diff(field[j],x[k])for k in range(5)));orbit[j,i]=-orbit[i,j]
            first=s.Matrix([*[s.Poly(orbit[i,3],z,w).coeff_monomial(z*z)for i in range(3)],0,0]);second=s.Matrix([*[s.Poly(orbit[i,4],z,w).coeff_monomial(w*w)for i in range(3)],0,0])
            canonical_orbit=orbit+wedge(x,z*first+w*second)
            replacement={}
            if name=='two_roots':replacement.update({B[i]:-a*Z[i]for i in range(3)});replacement.update({E[i]:a*Y[i]for i in range(3)})
            else:replacement.update({D[i]:a*Y[i]for i in range(3)});replacement.update({E[i]:a*Z[i]for i in range(3)})
            assert all(s.expand(value)==0 for value in canonical_orbit-kernel.subs(replacement,simultaneous=True))
        t2=s.symbols('r0:9');P2=s.zeros(5)
        for k,(i,j)in enumerate(combinations(range(3),2)):
            P2[i,j]=t2[3*k]*z*z+t2[3*k+1]*z*w+t2[3*k+2]*w*w;P2[j,i]=-P2[i,j]
        weight2=projective(bilinear(P0,P2,x),x)
        for i,j in combinations(range(3),2):assert s.expand(weight2[i,j,3,4]+2*K*P2[i,j])==0
        actual=s.expand(sum((-1)**i*x[i]*(P0[j,k]*P0[l,m]-P0[j,l]*P0[k,m]+P0[j,m]*P0[k,l])for i in range(5)for j,k,l,m in [[t for t in range(5)if t!=i]]))
        assert s.expand(actual-3*f*K)==0
        records.append({'normal_quadratic':name,'all350_full_weight1_Jacobi':all350(jac,x),'complete30_source':{str(i)+','+str(j):str(s.expand(P1[i,j]))for i,j in PAIRS if P1[i,j]!=0},
            'strict_full_firstweight_formula_verified':True,'all_base2_normal2_forcing_coefficients':{str(i)+','+str(j):str(s.expand(forcing[i,j]))for i,j in combinations(range(3),2)},
            'all350_kernel_substitution_zero':True,'complete_kernel_representative':{str(i)+','+str(j):str(s.expand(kernel[i,j]))for i,j in PAIRS if kernel[i,j]!=0},
            'actual_shear_orbit_equals_full_kernel':K!=0,'all9_weight2_injectivity_formula_verified':True,'all350_weight2_Jacobi':all350(weight2,x),'actual_D_equals_three_f_K':True})
        print(name,'all350/full30/kernel/shear/weight2 exact PASS',flush=True)
    result={'task_id':'G','all_passed':True,'all10_cubic_coefficients_symbolic':True,'smoothness_used_only_in_geometric_constant_direction_injectivity':True,
        'nonzero_a_localization_only_for_two_nonzero_binary_quadratic_types':True,'cubic':str(f),'normal_filtration_complete30_and9':True,'cases':records,
        'producer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_bytes':Path(__file__).read_text(),'seconds':time.monotonic()-start,
        'maximum_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'scope':'Exact full Schouten/Euler arithmetic with arbitrary ternary cubic coefficients. All-smooth injectivity, conormal scalar rigidity, source completeness and finite shear geometry are separate proofs in the note; no numerical Hesse specialization or generic rank argument.'}
    destination.write_text(json.dumps(result,indent=2)+'\n')





def reduced_m3_nine_main():
    """One whole finite batch after the complete root-cover geometry."""
    import argparse
    from fractions import Fraction as F
    import classification_degenerate_quadric_q4_nil31_I_review as q
    ap=argparse.ArgumentParser();ap.add_argument('--reduced-m3-nine',action='store_true');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    def addp(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in PAIRS})
    def act(P):
        terms=[]
        for missing in range(5):
            i,j,k,l=[a for a in range(5)if a!=missing]
            pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
            terms.append(q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
        return q.add(*terms)
    def cert(M,n):
        A=[[F(v)for v in row]for row in M];r=0;cols=[]
        for j in range(n):
            k=next((i for i in range(r,len(A))if A[i][j]),None)
            if k is None:continue
            A[r],A[k]=A[k],A[r];d=A[r][j];A[r]=[v/d for v in A[r]]
            for i in range(len(A)):
                if i!=r and A[i][j]:
                    d=A[i][j];A[i]=[a-d*b for a,b in zip(A[i],A[r])]
            cols.append(j);r+=1
        K=[]
        for j in range(n):
            if j in cols:continue
            v=[F()for _ in range(n)];v[j]=F(1)
            for row,pivot in zip(A,cols):v[pivot]=-row[j]
            assert all(not sum(a*b for a,b in zip(row,v))for row in M)
            K.append(v)
        rows=q.rank(q.transpose(M))[1] if M and n else []
        det=q.determinant([[M[i][j]for j in cols]for i in rows]) if r else F(1)
        assert det and len(rows)==r
        return dict(rank=r,nullity=n-r,rows=rows,columns=cols,minor_determinant=str(det),complete_kernel=[[str(v)for v in row]for row in K],full_kernel_substitution_zero=True),K
    def pdecode(text,names):
        poly=s.Poly(s.sympify(text,locals=dict(zip(names,s.symbols(' '.join(names))))),*s.symbols(' '.join(names)),domain=s.QQ)
        return {e+(0,0,0):F(int(c.p),int(c.q))for e,c in poly.terms()if c}
    def mat(M):return [[str(v)for v in row]for row in M]
    first=[]
    for normal in [3,4]:
        for ij in combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,mon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    assert len(first)==30 and len(second)==9
    d3mon=[e for e in q.monomials(5)if e[3]+e[4]==3];assert len(d3mon)==24
    euler=[]
    for a in range(5):
        for b in range(5):
            euler.append(q.antisymmetric({ij:q.add(q.atom(1,ij[0],b)if ij[1]==a else {},q.atom(-1,ij[1],b)if ij[0]==a else {})for ij in PAIRS}))
    EC,EK=cert(q.transpose([q.vec(P)for P in euler]),25);assert EC['rank']==24
    files=['research/classification/inputs/plane_tangent_conic_m3_four_first_symbols_20261108.json','research/classification/inputs/plane_triangle_m3_five_first_symbols_20261108.json']
    records=[]
    for fp in files:
        data=json.loads((ROOT/fp).read_text());names=data['variables']
        for profile in data['profiles']:
            P0=q.antisymmetric({tuple(names.index(v)for v in ij.split(',')):pdecode(text,names)for ij,text in profile['raw_brackets'].items()})
            assert not any(q.jvec(q.jacobi(P0))) and not act(P0)
            def op(P):
                v=q.jvec(q.jacobi(addp(P0,P)));own=q.jvec(q.jacobi(P));return [a-b for a,b in zip(v,own)]
            J1=q.transpose([op(P)for P in first]);J2=q.transpose([op(P)for P in second])
            R1,Kraw=cert(J1,30);R2,K2=cert(J2,9)
            D3=[[act(addp(P0,P)).get(e+(0,0,0),F())for P in first]for e in d3mon]
            target=[pdecode(t,names)for t in profile.get('whole_parallel_G3_basis',profile.get('whole_leading_target_basis',[]))]
            targetrows=[[t.get(e+(0,0,0),F())for e in d3mon]for t in target]
            TC,ann=cert(targetrows,24)
            constrained=J1+[[sum(a*b for a,b in zip(v,col))for col in zip(*D3)]for v in ann]
            RC,K1=cert(constrained,30)
            # Verify the entire P0 orbit in the actual Euler quotient, not a selected stabilizer.
            orbit=[]
            for a,b in [(a,b)for block in [range(3),range(3,5)]for a in block for b in block]:
                out={}
                for i,j in PAIRS:
                    out[i,j]=q.add(q.scale(q.mul(q.atom(1,b),q.der(P0.get((i,j),{}),a)),-1),P0.get((b,j),{})if i==a else {},P0.get((i,b),{})if j==a else {})
                orbit.append(q.antisymmetric(out))
            GM=q.transpose([q.vec(P)for P in euler+orbit+[P0]]);GC,GK=cert(GM,39)
            B=GC['rank']-24;assert B<=profile['first_symbol_source_dimension_upper']
            dim=len(K1)+len(K2);badzero=all(not v[k]for v in K1 for k in range(21,24))
            # Every P0 here preserves u=0. If every full firstweight bad row is zero,
            # all higher u-rows have factor u, forcing a proper factor of actual D.
            assert all(not e[3]==0 for (i,j),p in P0.items()if j==3 for e in p)
            bound=B+dim+5
            record=dict(input=fp,name=profile['name'],full150_P0=list(map(str,q.vec(P0))),full350_P0_Jacobi_zero=True,actual_P0_zero=True,
                complete350x30=mat(J1),weight1_unconstrained_certificate=R1,complete350x9=mat(J2),weight2_certificate=R2,
                full24x30_actual_D3=mat(D3),full_parallel_target_rows=mat(targetrows),target_certificate=TC,complete_target_annihilator=mat(ann),
                constrained_weight1_matrix=mat(constrained),constrained_weight1_certificate=RC,
                full150x39_Euler_Levi_scale=mat(GM),filtered_group_certificate=GC,first_symbol_dimension=B,proof_upper=profile['first_symbol_source_dimension_upper'],
                full_higher_dimension_upper=dim,whole_projective_image_upper=bound,entire_bad_u_row_zero=badzero,actual_D_has_proper_u_factor=badzero,G24_closed=badzero or bound<=24)
            records.append(record);print(profile['name'],B,len(K1),len(K2),bound,'factor',badzero,flush=True)
            result=dict(task_id='G',one_complete_nine_profile_batch=True,cases=records,completed_profiles=len(records),all_passed=len(records)==9 and all(r['G24_closed']for r in records),
                complete30_and9=True,no_quadratic_jet_or_GB=True,source_input_files=[dict(path=p,sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest())for p in files],
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_source=Path(__file__).read_text(),
                seconds=time.monotonic()-begun,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                scope='Whole higher necessary supersource for the nine finite first-symbol orbits. Root-cover geometric exhaustion and all-point E assignment are separate obligations. No nonclosing profile is certified as G closed.')
            args.output.write_text(json.dumps(result,indent=2)+'\n')




def triple_line_m2_double_main():
    """One whole eight-symbol batch for triple-line normal-double m2."""
    import argparse
    from fractions import Fraction as F
    import classification_degenerate_quadric_q4_nil31_I_review as q
    ap=argparse.ArgumentParser();ap.add_argument('--triple-line-m2-double',action='store_true');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    def addp(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in PAIRS})
    def act(P):
        terms=[]
        for missing in range(5):
            i,j,k,l=[a for a in range(5)if a!=missing]
            pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
            terms.append(q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
        return q.add(*terms)
    def cert(M,n):
        A=[[F(v)for v in row]for row in M];r=0;cols=[]
        for j in range(n):
            k=next((i for i in range(r,len(A))if A[i][j]),None)
            if k is None:continue
            A[r],A[k]=A[k],A[r];d=A[r][j];A[r]=[v/d for v in A[r]]
            for i in range(len(A)):
                if i!=r and A[i][j]:
                    d=A[i][j];A[i]=[a-d*b for a,b in zip(A[i],A[r])]
            cols.append(j);r+=1
        K=[]
        for j in range(n):
            if j in cols:continue
            v=[F()for _ in range(n)];v[j]=F(1)
            for row,pivot in zip(A,cols):v[pivot]=-row[j]
            assert all(not sum(a*b for a,b in zip(row,v))for row in M)
            K.append(v)
        rows=q.rank(q.transpose(M))[1] if M and n else []
        det=q.determinant([[M[i][j]for j in cols]for i in rows]) if r else F(1)
        assert det and len(rows)==r
        return dict(rank=r,nullity=n-r,rows=rows,columns=cols,minor_determinant=str(det),complete_kernel=[[str(v)for v in row]for row in K],full_kernel_substitution_zero=True),K
    def pdecode(text,names):
        poly=s.Poly(s.sympify(text,locals=dict(zip(names,s.symbols(' '.join(names))))),*s.symbols(' '.join(names)),domain=s.QQ)
        return {e+(0,0,0):F(int(c.p),int(c.q))for e,c in poly.terms()if c}
    def mat(M):return [[str(v)for v in row]for row in M]
    first=[]
    for normal in [3,4]:
        for ij in combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,mon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    assert len(first)==30 and len(second)==9
    d3mon=[e for e in q.monomials(5)if e[3]+e[4]==3];assert len(d3mon)==24
    euler=[]
    for a in range(5):
        for b in range(5):
            euler.append(q.antisymmetric({ij:q.add(q.atom(1,ij[0],b)if ij[1]==a else {},q.atom(-1,ij[1],b)if ij[0]==a else {})for ij in PAIRS}))
    EC,EK=cert(q.transpose([q.vec(P)for P in euler]),25);assert EC['rank']==24
    files=['research/classification/inputs/plane_triple_line_m2_constant_double_eight_first_symbols_20261108.json']
    records=[]
    for fp in files:
        data=json.loads((ROOT/fp).read_text());names=data['variables']
        for profile in data['profiles']:
            P0=q.antisymmetric({tuple(names.index(v)for v in ij.split(',')):pdecode(text,names)for ij,text in profile['raw_brackets'].items()})
            assert not any(q.jvec(q.jacobi(P0))) and act(P0)==pdecode(profile['G2'],names)
            def op(P):
                v=q.jvec(q.jacobi(addp(P0,P)));own=q.jvec(q.jacobi(P));return [a-b for a,b in zip(v,own)]
            J1=q.transpose([op(P)for P in first]);J2=q.transpose([op(P)for P in second])
            R1,Kraw=cert(J1,30);R2,K2=cert(J2,9)
            RC,K1=R1,Kraw
            # Verify the entire P0 orbit in the actual Euler quotient, not a selected stabilizer.
            orbit=[]
            for a,b in [(a,b)for block in [range(3),range(3,5)]for a in block for b in block]:
                out={}
                for i,j in PAIRS:
                    out[i,j]=q.add(q.scale(q.mul(q.atom(1,b),q.der(P0.get((i,j),{}),a)),-1),P0.get((b,j),{})if i==a else {},P0.get((i,b),{})if j==a else {})
                orbit.append(q.antisymmetric(out))
            GM=q.transpose([q.vec(P)for P in euler+orbit+[P0]]);GC,GK=cert(GM,39)
            B=GC['rank']-24
            dim=len(K1)+len(K2);badzero=all(not v[k]for v in K1 for k in range(21,24))
            def combination(bs,vs):return q.antisymmetric({ij:q.add(*(q.scale(P.get(ij,{}),c)for P,c in zip(bs,vs)if c))for ij in PAIRS})
            families=[P0]+[combination(first,v)for v in K1]+second
            afactor=all(e[0]>0 for P in families for i in range(1,5)for e,c in P.get((0,i),{}).items()if c)
            ufactor=all(e[3]>0 for P in families for i in [0,1,2,4]for e,c in P.get((i,3),{}).items()if c)
            assert ufactor==badzero
            bound=B+dim+5
            record=dict(input=fp,name=profile['name'],full150_P0=list(map(str,q.vec(P0))),full350_P0_Jacobi_zero=True,actual_P0_zero=False,
                complete350x30=mat(J1),weight1_unconstrained_certificate=R1,complete350x9=mat(J2),weight2_certificate=R2,
                actual_G2=profile['G2'],constrained_weight1_certificate=RC,
                full150x39_Euler_Levi_scale=mat(GM),filtered_group_certificate=GC,first_symbol_dimension=B,
                full_higher_dimension_upper=dim,whole_projective_image_upper=bound,entire_bad_u_row_zero=badzero,actual_D_has_proper_u_factor=ufactor,actual_D_has_proper_a_factor=afactor,G24_closed=ufactor or afactor or bound<=24)
            records.append(record);print(profile['name'],B,len(K1),len(K2),bound,'u_factor',ufactor,'a_factor',afactor,flush=True)
            result=dict(task_id='G',one_complete_eight_profile_batch=True,cases=records,completed_profiles=len(records),all_passed=len(records)==8 and all(r['G24_closed']for r in records),
                complete30_and9=True,no_quadratic_jet_or_GB=True,source_input_files=[dict(path=p,sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest())for p in files],
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_source=Path(__file__).read_text(),
                seconds=time.monotonic()-begun,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                scope='Whole higher necessary supersource for all triple-line constant-double m2 symbols. Whole geometric exhaustion is independently reviewed; all-point E assignment is separate. No nonclosing profile is certified as G closed.')
            args.output.write_text(json.dumps(result,indent=2)+'\n')



def triple_line_m2_sf01_main():
    """One whole five-symbol batch for triple-line squarefree flag01 m2."""
    import argparse
    from fractions import Fraction as F
    import classification_degenerate_quadric_q4_nil31_I_review as q
    ap=argparse.ArgumentParser();ap.add_argument('--triple-line-m2-sf01',action='store_true');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    def addp(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{})for P in Ps))for ij in PAIRS})
    def act(P):
        terms=[]
        for missing in range(5):
            i,j,k,l=[a for a in range(5)if a!=missing]
            pf=q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{})))
            terms.append(q.scale(q.mul(q.atom(1,missing),pf),(-1)**missing))
        return q.add(*terms)
    def cert(M,n):
        A=[[F(v)for v in row]for row in M];r=0;cols=[]
        for j in range(n):
            k=next((i for i in range(r,len(A))if A[i][j]),None)
            if k is None:continue
            A[r],A[k]=A[k],A[r];d=A[r][j];A[r]=[v/d for v in A[r]]
            for i in range(len(A)):
                if i!=r and A[i][j]:
                    d=A[i][j];A[i]=[a-d*b for a,b in zip(A[i],A[r])]
            cols.append(j);r+=1
        K=[]
        for j in range(n):
            if j in cols:continue
            v=[F()for _ in range(n)];v[j]=F(1)
            for row,pivot in zip(A,cols):v[pivot]=-row[j]
            assert all(not sum(a*b for a,b in zip(row,v))for row in M)
            K.append(v)
        rows=q.rank(q.transpose(M))[1] if M and n else []
        det=q.determinant([[M[i][j]for j in cols]for i in rows]) if r else F(1)
        assert det and len(rows)==r
        return dict(rank=r,nullity=n-r,rows=rows,columns=cols,minor_determinant=str(det),complete_kernel=[[str(v)for v in row]for row in K],full_kernel_substitution_zero=True),K
    def pdecode(text,names):
        poly=s.Poly(s.sympify(text,locals=dict(zip(names,s.symbols(' '.join(names))))),*s.symbols(' '.join(names)),domain=s.QQ)
        return {e+(0,0,0):F(int(c.p),int(c.q))for e,c in poly.terms()if c}
    def mat(M):return [[str(v)for v in row]for row in M]
    first=[]
    for normal in [3,4]:
        for ij in combinations(range(3),2):
            for k in range(3):first.append(q.antisymmetric({ij:q.atom(1,normal,k)}))
    for normal,mon in [(3,[3,4]),(3,[4,4]),(4,[3,3]),(4,[3,4])]:
        for i in range(3):first.append(q.antisymmetric({(i,normal):q.atom(1,*mon)}))
    second=[q.antisymmetric({ij:q.atom(1,*mon)})for ij in combinations(range(3),2)for mon in [[3,3],[3,4],[4,4]]]
    assert len(first)==30 and len(second)==9
    d3mon=[e for e in q.monomials(5)if e[3]+e[4]==3];assert len(d3mon)==24
    euler=[]
    for a in range(5):
        for b in range(5):
            euler.append(q.antisymmetric({ij:q.add(q.atom(1,ij[0],b)if ij[1]==a else {},q.atom(-1,ij[1],b)if ij[0]==a else {})for ij in PAIRS}))
    EC,EK=cert(q.transpose([q.vec(P)for P in euler]),25);assert EC['rank']==24
    files=['research/classification/inputs/plane_triple_line_m2_sf01_five_first_symbols_20261108.json']
    records=[]
    for fp in files:
        data=json.loads((ROOT/fp).read_text());names=data['variables']
        for profile in data['profiles']:
            P0=q.antisymmetric({tuple(names.index(v)for v in ij.split(',')):pdecode(text,names)for ij,text in profile['raw_brackets'].items()})
            assert not any(q.jvec(q.jacobi(P0))) and act(P0)==pdecode(profile['G2'],names)
            def op(P):
                v=q.jvec(q.jacobi(addp(P0,P)));own=q.jvec(q.jacobi(P));return [a-b for a,b in zip(v,own)]
            J1=q.transpose([op(P)for P in first]);J2=q.transpose([op(P)for P in second])
            R1,Kraw=cert(J1,30);R2,K2=cert(J2,9)
            RC,K1=R1,Kraw
            # Verify the entire P0 orbit in the actual Euler quotient, not a selected stabilizer.
            orbit=[]
            for a,b in [(a,b)for block in [range(3),range(3,5)]for a in block for b in block]:
                out={}
                for i,j in PAIRS:
                    out[i,j]=q.add(q.scale(q.mul(q.atom(1,b),q.der(P0.get((i,j),{}),a)),-1),P0.get((b,j),{})if i==a else {},P0.get((i,b),{})if j==a else {})
                orbit.append(q.antisymmetric(out))
            GM=q.transpose([q.vec(P)for P in euler+orbit+[P0]]);GC,GK=cert(GM,39)
            B=GC['rank']-24
            dim=len(K1)+len(K2);badzero=all(not v[k]for v in K1 for k in range(21,24))
            def combination(bs,vs):return q.antisymmetric({ij:q.add(*(q.scale(P.get(ij,{}),c)for P,c in zip(bs,vs)if c))for ij in PAIRS})
            families=[P0]+[combination(first,v)for v in K1]+second
            afactor=all(e[0]>0 for P in families for i in range(1,5)for e,c in P.get((0,i),{}).items()if c)
            ufactor=all(e[3]>0 for P in families for i in [0,1,2,4]for e,c in P.get((i,3),{}).items()if c)
            assert ufactor==badzero
            bound=B+dim+5
            record=dict(input=fp,name=profile['name'],full150_P0=list(map(str,q.vec(P0))),full350_P0_Jacobi_zero=True,actual_P0_zero=False,
                complete350x30=mat(J1),weight1_unconstrained_certificate=R1,complete350x9=mat(J2),weight2_certificate=R2,
                actual_G2=profile['G2'],constrained_weight1_certificate=RC,
                full150x39_Euler_Levi_scale=mat(GM),filtered_group_certificate=GC,first_symbol_dimension=B,
                full_higher_dimension_upper=dim,whole_projective_image_upper=bound,entire_bad_u_row_zero=badzero,actual_D_has_proper_u_factor=ufactor,actual_D_has_proper_a_factor=afactor,G24_closed=ufactor or afactor or bound<=24)
            records.append(record);print(profile['name'],B,len(K1),len(K2),bound,'u_factor',ufactor,'a_factor',afactor,flush=True)
            result=dict(task_id='G',one_complete_five_profile_batch=True,cases=records,completed_profiles=len(records),all_passed=len(records)==5 and all(r['G24_closed']for r in records),
                complete30_and9=True,no_quadratic_jet_or_GB=True,source_input_files=[dict(path=p,sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest())for p in files],
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),producer_source=Path(__file__).read_text(),
                seconds=time.monotonic()-begun,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                scope='Whole higher necessary supersource for all triple-line squarefree flag01 m2 symbols. Whole geometric exhaustion is independently reviewed; all-point E assignment is separate. No nonclosing profile is certified as G closed.')
            args.output.write_text(json.dumps(result,indent=2)+'\n')


def double_line_m2_double_main():
    """One bounded whole34 R=2L+M normal-double source batch."""
    import argparse,ast
    from fractions import Fraction as F
    import classification_plane_m2_concurrent_double_review as v
    from classification_unimodular_logarithmic_review import RF
    H=v.H;compact=v.compact;AF=v.Quotient
    ap=argparse.ArgumentParser(description=double_line_m2_double_main.__doc__)
    ap.add_argument('--double-line-m2-double',action='store_true');ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    assert not args.output.exists();checkpoint=args.output.with_suffix('.partial.json');assert not checkpoint.exists()
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(590,595));resource.setrlimit(resource.RLIMIT_AS,(3072*1024**2,)*2)
    job=json.loads((ROOT/'research/classification/jobs/plane_double_line_m2_constant_double_thirtyfour.json').read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    PPAR=RF((F(),F(1)));coords=('a','b','c','u','v');symbol=s.Symbol('p')
    def parse(text):
        def go(n):
            if isinstance(n,ast.Constant):return H.a(F(n.value))
            if isinstance(n,ast.Name):return H.scale(H.a(1),PPAR)if n.id=='p'else H.a(1,coords.index(n.id))
            if isinstance(n,ast.UnaryOp):return H.scale(go(n.operand),-1 if isinstance(n.op,ast.USub)else 1)
            if isinstance(n,ast.BinOp):
                if isinstance(n.op,ast.Add):return H.add(go(n.left),go(n.right))
                if isinstance(n.op,ast.Sub):return H.add(go(n.left),H.scale(go(n.right),-1))
                if isinstance(n.op,ast.Mult):return H.mul(go(n.left),go(n.right))
                if isinstance(n.op,ast.Div):
                    den=go(n.right);assert set(den)=={(0,)*5};return H.scale(go(n.left),1/den[(0,)*5])
                if isinstance(n.op,ast.Pow):
                    out=H.a(1)
                    for _ in range(n.right.value):out=H.mul(out,go(n.left))
                    return out
            raise ValueError(ast.dump(n))
        return go(ast.parse(text,mode='eval').body)
    def enc(x):
        x=compact(x)
        if isinstance(x,RF):return dict(numerator=list(map(str,x.p)),denominator=list(map(str,x.q)))
        if isinstance(x,AF):return dict(residue=list(map(str,x.p)),modulus=list(map(str,AF.modulus)))
        return str(x)
    def mat(M):return[[enc(x)for x in row]for row in M]
    def pol(P):return[dict(exponent=list(e),coefficient=enc(c))for e,c in sorted(P.items())]
    def certificate(M,n):
        B=[[compact(x)for x in row]for row in M];ids=list(range(len(B)));rows=[];piv=[];det=F(1)
        for j in range(n):
            r=len(piv);k=next((k for k in range(r,len(B))if B[k][j]),None)
            if k is None:continue
            B[r],B[k]=B[k],B[r];ids[r],ids[k]=ids[k],ids[r];d=B[r][j];rows.append(ids[r]);piv.append(j);det=compact(det*d)
            B[r]=[compact(x/d)if x else F()for x in B[r]];active=[l for l in range(j,n)if B[r][l]]
            for k in range(r+1,len(B)):
                c=B[k][j]
                if c:
                    for l in active:B[k][l]=compact(B[k][l]-c*B[r][l])
        for r in range(len(piv)-1,-1,-1):
            j=piv[r];active=[l for l in range(j+1,n)if B[r][l]]
            for k in range(r):
                c=B[k][j]
                if c:
                    for l in active:B[k][l]=compact(B[k][l]-c*B[r][l])
                    B[k][j]=F()
        K=[]
        for j in range(n):
            if j not in piv:
                w=[F()for _ in range(n)];w[j]=F(1)
                for row,i in zip(B,piv):w[i]=-row[j]
                K.append(w)
        for row in M:
            nz=[(i,c)for i,c in enumerate(row)if c]
            assert all(not sum((c*w[i]for i,c in nz if w[i]),F())for w in K)
        if isinstance(det,AF):assert det*det.inv()==1
        return dict(rank=len(piv),nullity=n-len(piv),rows=rows,columns=piv,minor_determinant=enc(det),complete_kernel=mat(K),all_original_rows_annihilate_kernel=True),K,det
    B1,B2=H.sources();bad=[[Q.get((i,3),{}).get((0,0,0,0,2),F())for Q in B1]for i in range(3)]
    empty_group=v.D.group(H.skew({}));assert len(H.rref([[F(x)for x in row[:25]]for row in empty_group])[1])==24
    assert all(not Q.get((i,3),{})for Q in B2 for i in range(5)if i!=3)
    def operators(P):return v.D.op(P,B1),v.D.op(P,B2),v.D.group(P)
    def process(name,P,target,pdim,field,prepared=None):
        assert not any(H.jac(P).values()),(name,'P0 full350 Jacobi')
        assert H.actual(P)==target and target,(name,'P0 actualG2')
        J1,J2,G=operators(P)if prepared is None else prepared
        C1,K1,d1=certificate(J1,30);C2,K2,d2=certificate(J2,9);CG,KG,dg=certificate(G,39)
        badmaps=[[sum((a*b for a,b in zip(row,w)),F())for w in K1]for row in bad]
        assert all(e[3]for i in range(5)if i!=3 for e in P.get((i,3),{}))
        ufactor=not any(any(r)for r in badmaps);stabilizer=37-CG['rank'];fiber=stabilizer+6
        upper=24-fiber+pdim+len(K1)+len(K2)
        row=dict(name=name,field=field,source_parameter_dimension=pdim,full150_P0=list(map(enc,H.vector(P))),full350_P0_Jacobi_zero=True,
            actualG2=pol(target),complete350x30=mat(J1),complete350x9=mat(J2),weight1_certificate=C1,weight2_certificate=C2,
            full150x39_Euler_Levi_scale=mat(G),normalizer_certificate=CG,Euler_rank=24,graded_projective_Levi_stabilizer=stabilizer,actual_filtered_group_fiber_dimension=fiber,
            first_symbol_affine_dimension_upper=CG['rank']-24+pdim,higher_fiber_upper=len(K1)+len(K2),whole_projective_image_upper=upper,
            all_bad_urow_maps=mat(badmaps),actual_u_factor_for_every_FULL9_P2_affine_solution=ufactor,G24_closed=ufactor or upper<=24)
        print(json.dumps(dict(name=name,field=field,K1=len(K1),K2=len(K2),group=fiber,pdim=pdim,upper=upper,u_factor=ufactor,closed=row['G24_closed'])),flush=True)
        return row,(d1,d2,dg),K1+K2+KG
    def factors_from(dets,K):
        found={}
        for value in list(dets)+[x for row in K for x in row]:
            z=RF(value)
            for pp in(z.p,z.q):
                if len(pp)>1:
                    expr=sum(s.Rational(c.numerator,c.denominator)*symbol**j for j,c in enumerate(pp))
                    for g,multiplicity in s.factor_list(expr,symbol)[1]:
                        coeff=tuple(F(int(c.p),int(c.q))for c in reversed(s.Poly(g,symbol).monic().all_coeffs()));found[coeff]=str(g)
        return found
    ip='research/classification/inputs/plane_double_line_m2_constant_double_thirtyfour_symbols_20261108.json'
    data=json.loads((ROOT/ip).read_text());assert len(data['profiles'])==34
    records=[];generic_count=wall_count=algebraic_count=0;sourcebytes=Path(__file__).read_text();sourcesha=hashlib.sha256(sourcebytes.encode()).hexdigest()
    for profile in data['profiles']:
        P=H.skew({tuple(coords.index(x)for x in ij.split(',')):parse(text)for ij,text in profile['raw_brackets'].items()})
        target=parse(profile['G2']);pdim=profile['source_parameter_dimension'];field='QQ(p)'if pdim else'QQ'
        ops=operators(P);record,dets,K=process(profile['name'],P,target,pdim,field,ops)
        record['source_profile']=profile
        if pdim:
            generic_count+=1;factors=factors_from(dets,K);excluded=set(map(F,profile['domain_excluded_p_values']));walls=[];skip=[]
            for g,label in sorted(factors.items()):
                if len(g)==2 and -g[0]in excluded:
                    skip.append(dict(p=str(-g[0]),minimal_polynomial=list(map(str,g)),reason='actual leadingG2 zero; normal order/partition outside this source'));continue
                if len(g)==2:
                    point=-g[0];description=dict(p=str(point),minimal_polynomial=list(map(str,g)));wallfield='QQ'
                else:
                    AF.modulus=g;point=AF((0,1));description=dict(minimal_polynomial=list(map(str,g)),all_complex_roots_covered_by_quotient_unit_minors=True);wallfield='QQ[p]/('+label+')';algebraic_count+=1
                WP=v.specialize_tensor(P,point);WT=v.specialize_poly(target,point)
                prepared=tuple([[v.hom(x,point)for x in row]for row in O]for O in ops)
                wall,_,_=process(profile['name']+'_WALL_'+label,WP,WT,0,wallfield,prepared);wall.update(description);walls.append(wall);wall_count+=1
            record.update(all_generic_minor_and_complete_kernel_factors=[list(map(str,g))for g in factors],all_allowed_parameter_walls=walls,
                excluded_leading_points=skip,source_domain_excluded_p_values=profile['domain_excluded_p_values'],whole_parameter_line_G24_closed=record['G24_closed']and all(w['G24_closed']for w in walls))
        else:record['whole_parameter_line_G24_closed']=record['G24_closed']
        records.append(record)
        checkpoint.write_text(json.dumps(dict(task_id='G',completed=records,all34_complete=False,producer_sha256=sourcesha,producer_source=sourcebytes,inputs=job['exact_inputs_and_hashes']),indent=2)+'\n')
    assert generic_count==11
    result=dict(task_id='G',all_arithmetic_passed=True,all34_complete=True,all34_G24_closed=all(r['whole_parameter_line_G24_closed']for r in records),
        generic_parameter_lines=generic_count,fixed_sources=23,allowed_wall_fields=wall_count,nonrational_wall_fields=algebraic_count,cases=records,
        full30_weight1_basis=[list(map(str,H.vector(Q)))for Q in B1],full9_weight2_basis=[list(map(str,H.vector(Q)))for Q in B2],
        full350_row_labels=[dict(quad=list(q),monomial=list(e))for q,e in H.ROWS],inputs=job['exact_inputs_and_hashes'],producer_sha256=sourcesha,producer_source=sourcebytes,
        seconds=time.monotonic()-begun,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Exactly full34 necessary R=2L+M m2constant-double sources, eleven whole QQ(p) lines and23fixed inputs; full9 P2 retained. All unclosed linear supersources preserved. Independent arithmetic review pending.',
        next_if_not_closed='Read only the saved entire unclosed source in Ore/canonical geometry; no quadratic jet/minor/profile or GB series.')
    args.output.write_text(json.dumps(result,indent=2)+'\n');checkpoint.unlink()
    print(json.dumps({k:result[k]for k in('all_arithmetic_passed','all34_G24_closed','allowed_wall_fields','seconds','maximum_RSS_KiB')}))


if __name__=='__main__':
    import sys
    if '--double-line-m2-double' in sys.argv:double_line_m2_double_main()
    elif '--reduced-m3-nine' in sys.argv:reduced_m3_nine_main()
    elif '--triple-line-m2-double' in sys.argv:triple_line_m2_double_main()
    elif '--triple-line-m2-sf01' in sys.argv:triple_line_m2_sf01_main()
    else:main()
