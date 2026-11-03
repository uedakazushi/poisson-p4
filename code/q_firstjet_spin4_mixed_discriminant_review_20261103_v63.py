#!/usr/bin/env python3
"""Independent Fraction/coefficients review of the Spin4 mixed discriminant.

No author module is imported. No composed12-variable discriminant, new first-jet
sample, high-degree invariant probe, or generic pole-order calculation is made.
"""
import argparse
from collections import defaultdict
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

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63_REVIEW'


def read(path):
    return json.loads(gzip.decompress(path.read_bytes())if path.suffix=='.gz' else path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def mat(A):return [[Q(x)for x in row]for row in A]
def zero(n,m):return [[Q(0)for j in range(m)]for i in range(n)]
def tp(A):return list(map(list,zip(*A)))
def mm(A,B):return [[sum((x*y for x,y in zip(row,col)),Q(0))for col in zip(*B)]for row in A]
def add(A,B,c=1):return [[x+c*y for x,y in zip(a,b)]for a,b in zip(A,B)]
def det(A):
    A=[row[:]for row in A];answer=Q(1)
    for j in range(len(A)):
        i=next((i for i in range(j,len(A))if A[i][j]),None)
        if i is None:return Q(0)
        if i!=j:A[i],A[j]=A[j],A[i];answer=-answer
        p=A[j][j];answer*=p
        for i in range(j+1,len(A)):
            c=A[i][j]/p
            for k in range(j+1,len(A)):A[i][k]-=c*A[j][k]
    return answer
def rank(A):
    A=[row[:]for row in A];r=0
    for j in range(len(A[0])):
        i=next((i for i in range(r,len(A))if A[i][j]),None)
        if i is None:continue
        A[i],A[r]=A[r],A[i];p=A[r][j]
        A[r]=[x/p for x in A[r]]
        for i in range(r+1,len(A)):
            c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r
def inv(A):
    n=len(A);B=[row[:]+[Q(i==j)for j in range(n)]for i,row in enumerate(A)]
    for j in range(n):
        i=next(i for i in range(j,n)if B[i][j]);B[i],B[j]=B[j],B[i]
        p=B[j][j];B[j]=[x/p for x in B[j]]
        for i in range(n):
            if i!=j:
                c=B[i][j];B[i]=[x-c*y for x,y in zip(B[i],B[j])]
    return [row[n:]for row in B]
def plus(A,B,c=1):
    C=dict(A)
    for mon,a in B.items():
        C[mon]=C.get(mon,Q(0))+c*a
        if not C[mon]:del C[mon]
    return C
def mul(A,B):
    C=defaultdict(Q)
    for m,a in A.items():
        for n,b in B.items():C[tuple(x+y for x,y in zip(m,n))]+=a*b
    return {m:a for m,a in C.items()if a}
def power(A,k,n):
    C={(0,)*n:Q(1)}
    for j in range(k):C=mul(C,A)
    return C
def variable(i,n):return {tuple(int(i==j)for j in range(n)):Q(1)}
def const(a,n):return {(0,)*n:Q(a)}if a else {}
def evalp(A,values):return sum((a*s.prod(Q(x)**k for x,k in zip(values,m))for m,a in A.items()),Q(0))
def parsep(text,variables):
    return {tuple(m):Q(str(c))for m,c in s.Poly(s.sympify(text),*variables).terms()if c}
def determinant_poly(A,nvariables):
    # Subset DP Laplace expansion, independent of resultant/discriminant APIs.
    states={0:const(1,nvariables)};size=len(A)
    for row in range(size):
        new={}
        for mask,value in states.items():
            for j in range(size):
                if mask>>j&1 or not A[row][j]:continue
                sign=(-1)**sum((mask>>k)&1 for k in range(j+1,size))
                target=mask|(1<<j)
                new[target]=plus(new.get(target,{}),mul(value,A[row][j]),sign)
        states=new
    return states[(1<<size)-1]
def convolution(A,B):
    C=[{}for _ in range(len(A)+len(B)-1)]
    for i,a in enumerate(A):
        for j,b in enumerate(B):C[i+j]=plus(C[i+j],mul(a,b))
    return C
def compose(A,values,n):
    C={}
    for m,c in A.items():
        term=const(c,n)
        for v,k in zip(values,m):term=mul(term,power(v,k,n))
        C=plus(C,term)
    return C


def main(output):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63/report.json'
    certificate=author.parent/'certificate.json';circuitpath=author.parent/'mixed-discriminant-circuit.json.gz'
    spinreview=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json'
    spincert=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    original=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    ar,cert,circuit,spin,orig=map(read,[author,certificate,circuitpath,spincert,original])
    assert ar['all_passed']and read(spinreview)['all_passed']
    checked={}
    for p in [author,spinreview]:
        for name,want in read(p)['input_sha256'].items():assert sha(ROOT/name)==want,name;checked[name]=want
    pairs=list(itertools.combinations(range(4),2));coords=[(ij,k)for ij in pairs for k in range(4)]
    J=mat([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    Ji=inv(J);gamma=[]
    for i,j in pairs[:-1]:
        alpha=zero(4,4);alpha[i][j]=1;alpha[j][i]=-1
        if (i,j)==(0,1):alpha[2][3]=-1;alpha[3][2]=1
        gamma.append(mm(alpha,Ji))
    assert gamma==list(map(mat,spin['gamma_matrices']))
    # Build each tensor coefficient from its symmetric spinor entries.
    P=[zero(5,4)for _ in range(12)]
    for i in range(3):
        for j in range(2):
            col=2*i+j
            if i==0:P[col][1+j][0]=1
            elif i==1:P[col][1+j][1]=Q(1,2);P[col][3+j][0]=Q(1,2)
            else:P[col][3+j][1]=1
    for i in range(2):
        for j in range(3):
            col=6+3*i+j
            if j==0:P[col][1+2*i][2]=1
            elif j==1:P[col][1+2*i][3]=Q(1,2);P[col][2+2*i][2]=Q(1,2)
            else:P[col][2+2*i][3]=1
    E=[[P[l][pairs.index(ij)][k]if ij!=(2,3)else Q(0)for l in range(12)]for ij,k in coords]
    assert E==mat(cert['physical_embedding24x12'])and rank(E)==12
    B=mat(orig['primitive24_by16_basis']);E16=mat(cert['accepted_embedding16x12'])
    assert rank(B)==16 and mm(B,E16)==E
    trace_cross=[]
    for out in range(4):
        trace_cross.append([gamma[e][out][k]for e in range(1,5)for k in range(4)])
    assert rank(trace_cross)==4
    for psi in P:
        assert all(sum(gamma[e][i][k]*psi[e][k]for e in range(5)for k in range(4))==0 for i in range(4))
        beta=[psi[pairs.index(ij)][k]if ij!=(2,3)else Q(0)for ij,k in coords]
        def entry(i,j,k):return beta[coords.index(((i,j),k))]if i<j else -beta[coords.index(((j,i),k))]
        assert all(entry(j,k,i)-entry(i,k,j)+entry(i,j,k)==0 for i,j,k in itertools.combinations(range(4),3))
    def forms(beta):
        out=[zero(4,4)for k in range(4)]
        for value,((i,j),k)in zip(beta,coords):out[k][i][j]=value;out[k][j][i]=-value
        return out
    generators=[]
    for i in range(4):
        for j in range(i,4):
            symmetric=zero(4,4);symmetric[i][j]=symmetric[j][i]=1
            generators.append(mm(J,symmetric))
    indices=[i for i,A in enumerate(generators)if all(A[j][k]==A[k][j]==0 for j in range(2)for k in range(2,4))]
    assert indices==cert['spin4_sp4_generator_indices']and len(indices)==6
    monU=[(2-i,i,1-j,j)for i in range(3)for j in range(2)]
    monV=[(1-i,i,2-j,j)for i in range(2)for j in range(3)]
    actions=[]
    for index in indices:
        A=generators[index];R=zero(12,12)
        for offset,mons in [(0,monU),(6,monV)]:
            for col,mon in enumerate(mons):
                for k in range(4):
                    if not mon[k]:continue
                    for l in range(4):
                        if A[k][l]:
                            target=list(mon);target[k]-=1;target[l]+=1
                            R[offset+mons.index(tuple(target))][offset+col]+=mon[k]*A[k][l]
        physical=[]
        for col in range(12):
            b=forms([row[col]for row in E]);bs=[]
            for k in range(4):
                value=add(mm(tp(A),b[k]),mm(b[k],A))
                for i in range(4):value=add(value,b[i],A[i][k])
                bs.append(value)
            physical.append([bs[k][i][j]for(i,j),k in coords])
        assert mm(E,R)==tp(physical)
        actions.append(R)
    assert actions==list(map(mat,cert['six_biform_action_matrices']))
    g=mat([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    assert g==mat(cert['normalizer4'])and mm(mm(tp(g),J),g)==J
    swap=zero(12,12)
    for i in range(2):
        for j in range(3):swap[2*j+i][6+3*i+j]=-1
    for i in range(3):
        for j in range(2):swap[6+3*j+i][2*i+j]=1
    assert swap==mat(cert['normalizer_biform12'])
    for col in range(12):
        b=forms([row[col]for row in E]);bs=[]
        for k in range(4):
            combination=zero(4,4)
            for i in range(4):combination=add(combination,b[i],g[i][k])
            bs.append(mm(mm(tp(g),combination),g))
        assert [bs[k][i][j]for(i,j),k in coords]==[sum(E[r][i]*swap[i][col]for i in range(12))for r in range(24)]
    # All quartic delta coefficients, independently by subset determinant.
    linearpsi=[[{}for k in range(4)]for e in range(4)]
    for e in range(4):
        for k in range(4):linearpsi[e][k]={tuple(int(j==l)for j in range(12)):P[l][e+1][k]for l in range(12)if P[l][e+1][k]}
    delta=determinant_poly(linearpsi,12);vars12=s.symbols(' '.join(cert['variables12']))
    assert delta==parsep(cert['normal_component_delta'],vars12)
    delta_swap={}
    for mon,value in delta.items():
        target=[0]*12;sign=1
        for i,k in enumerate(mon):
            j=next(j for j,x in enumerate(swap[i])if x);target[j]+=k;sign*=swap[i][j]**k
        delta_swap[tuple(target)]=delta_swap.get(tuple(target),Q(0))+value*sign
    assert delta_swap=={mon:-value for mon,value in delta.items()}
    Q5=mat(spin['quadratic_form5']);Qinverse=inv(Q5);assert det(Q5)==Q(1,16)
    for i in range(5):
        for j in range(5):
            minor=[[Q5[r][c]for c in range(5)if c!=j]for r in range(5)if r!=i]
            assert (-1)**(i+j)*det(minor)==det(Q5)*Qinverse[j][i]
    assert det([row[1:]for row in Q5[1:]])==Q(1,16)and det(Ji)==1
    # Independent quadratic Sylvester resultant and6x6 restriction determinant.
    vv=[variable(i,6)for i in range(6)];v0,v1=vv[:3],vv[3:]
    sylvester=[[{}for _ in range(4)]for _ in range(4)]
    for row,values,shift in [(0,v0,0),(1,v0,1),(2,v1,0),(3,v1,1)]:
        for i,q in enumerate(values):sylvester[row][i+shift]=q
    resultant=determinant_poly(sylvester,6)
    vnames=vars12[6:];assert resultant==parsep(cert['resultant'],vnames)
    quadratics=[convolution(v1,v1),[{m:-c for m,c in q.items()}for q in convolution(v0,v1)],convolution(v0,v0)]
    restriction=[[{}for _ in range(6)]for _ in range(6)]
    for col in range(6):
        shift=col%2
        for i,q in enumerate(quadratics[col//2]):restriction[i+shift][col]=q
    assert all(restriction[i][j]==parsep(cert['restriction_matrix'][i][j],vnames)for i in range(6)for j in range(6))
    assert determinant_poly(restriction,6)==power(resultant,3,6)
    quintic=[]
    for i in range(6):
        q={}
        for j in range(6):
            for vm,c in restriction[i][j].items():q[tuple(int(j==k)for k in range(6))+vm]=c
        quintic.append(q)
        assert q==parsep(cert['quintic_coefficients'][i],vars12)
    # Universal quintic Disc:9x9 Sylvester determinant / leading coefficient.
    ff=[variable(i,6)for i in range(6)];syl=[[{}for j in range(9)]for i in range(9)]
    for shift in range(4):
        for i,q in enumerate(ff):syl[shift][i+shift]=q
    for shift in range(5):
        for i in range(5):syl[4+shift][i+shift]={m:(i+1)*c for m,c in ff[i+1].items()}
    res=determinant_poly(syl,6);disc={}
    for mon,c in res.items():assert mon[5]>=1;disc[mon[:5]+(mon[5]-1,)]=c
    saved={tuple(mon):Q(c)for mon,c in cert['generic_quintic_discriminant59_terms']}
    assert disc==saved and len(disc)==59 and all(sum(mon)==8 for mon in disc)
    # Prove the unexpanded circuit universally by cutting at its six quintic
    # outputs and resultant, rather than extrapolating from fixture values.
    assert circuit['input_names']==cert['variables12']
    nodes=circuit['nodes'];outputs=circuit['outputs']
    assert len(nodes)==170
    for i,(op,value)in enumerate(nodes):
        if op in ['+','*']:assert all(0<=j<i for j in value)
        elif op=='inv':assert 0<=value<i
        elif op=='v':assert 0<=value<12
        else:assert op=='q';Q(value)
    lowcache={}
    def lower(i):
        if i in lowcache:return lowcache[i]
        op,value=nodes[i]
        if op=='q':out=const(Q(value),12)
        elif op=='v':out=variable(value,12)
        elif op=='+':
            out={}
            for j in value:out=plus(out,lower(j))
        elif op=='*':
            out=const(1,12)
            for j in value:out=mul(out,lower(j))
        else:raise AssertionError('An inverse before the quintic/resultant cut')
        lowcache[i]=out;return out
    fcuts=outputs['binary_quintic_coefficients'];rcut=outputs['resultant']
    assert [lower(i)for i in fcuts]==quintic
    assert lower(rcut)=={(0,)*6+m:c for m,c in resultant.items()}
    cuts={node:variable(i,7)for i,node in enumerate(fcuts)};cuts[rcut]=variable(6,7)
    topcache={}
    def upper(i):
        if i in cuts:return cuts[i]
        if i in topcache:return topcache[i]
        op,value=nodes[i]
        if op=='q':out=const(Q(value),7)
        elif op=='+':
            out={}
            for j in value:out=plus(out,upper(j))
        elif op=='*':
            out=const(1,7)
            for j in value:out=mul(out,upper(j))
        elif op=='inv':
            p=upper(value);assert len(p)==1
            mon,c=next(iter(p.items()));out={tuple(-k for k in mon):1/c}
        else:raise AssertionError('An original input bypasses the quintic/resultant cut')
        topcache[i]=out;return out
    assert upper(outputs['quintic_discriminant'])=={m+(0,):c for m,c in disc.items()}
    assert upper(outputs['mixed_discriminant'])=={m+(-2,):c for m,c in disc.items()}
    fixed=list(map(Q,cert['fixed_biform_values']));swapped=[sum(x*y for x,y in zip(row,fixed))for row in swap]
    assert swapped==list(map(Q,cert['fixed_swapped_values']))
    fixture_values=[]
    for values in [fixed,swapped]:
        f=[evalp(q,values)for q in quintic];r=evalp(resultant,values[6:]);d=evalp(disc,f)
        fixture_values.append((r,d/r**2))
        computed=[];degrees=[]
        for op,value in circuit['nodes']:
            if op=='q':x=Q(value);degree=(0,0)
            elif op=='v':x=values[value];degree=(1,0)if value<6 else(0,1)
            elif op=='+':
                x=sum((computed[j]for j in value),Q(0));nonzero={degrees[j]for j in value if circuit['nodes'][j]!=['q','0']};assert len(nonzero)<=1;degree=next(iter(nonzero),(0,0))
            elif op=='*':
                x=Q(1)
                for j in value:x*=computed[j]
                degree=tuple(sum(degrees[j][i]for j in value)for i in range(2))
            elif op=='inv':x=1/computed[value];degree=tuple(-i for i in degrees[value])
            else:raise AssertionError(op)
            computed.append(x);degrees.append(degree)
        o=circuit['outputs']
        assert [computed[j]for j in o['binary_quintic_coefficients']]==f
        assert computed[o['resultant']]==r and computed[o['quintic_discriminant']]==d
        assert computed[o['mixed_discriminant']]==d/r**2
        assert all(degrees[j]==(1,2)for j in o['binary_quintic_coefficients'])
        assert degrees[o['resultant']]==(0,4)and degrees[o['quintic_discriminant']]==(8,16)
        assert degrees[o['mixed_discriminant']]==(8,8)
    assert fixture_values==[(Q(50),Q(53268630696)),(Q(37),Q(53268630696))]
    # The boundary valuation fixture uses the same universal59 coefficients,
    # now with one epsilon variable only, never the12-variable composition.
    eps=variable(0,1);vc=[{},const(1,1),{},eps,{},const(1,1)]
    uc=[1,1,3,2,7,5];f_eps=[]
    for i in range(6):
        q={}
        for j in range(6):q=plus(q,compose(restriction[i][j],vc,1),uc[j])
        f_eps.append(q)
    fixture_disc=compose(disc,f_eps,1)
    assert compose(resultant,vc,1)==eps and min(m[0]for m in fixture_disc)==2
    quotient={(m[0]-2,):c for m,c in fixture_disc.items()}
    es=s.Symbol('e');assert quotient==parsep(cert['fixture_discriminant_quotient_by_e2'],[es])
    assert quotient[(0,)]==1466325
    evidence={
        'physical24x12_embedding':[[str(q)for q in row]for row in E],
        'full_cross_gamma_rank':4,'slice_dimension':12,
        'six_generator_indices':indices,'all_six_physical_biform_intertwinings':True,
        'full_finite_normalizer24_identities':True,'all_delta_coefficients_and_odd_parity':True,
        'all25_compound_Q_identities_for_Cauchy_Binet':True,
        'slice_p0_equals_delta_squared_over16':True,
        'independent_resultant_terms':len(resultant),
        'independent_restriction_determinant_R_cubed':True,
        'independent_quintic_discriminant_Sylvester_DP59_terms':[[list(m),str(c)]for m,c in sorted(disc.items())],
        'boundary_fixture_disc_R_valuation':2,'boundary_fixture_constant':1466325,
        'normalizer_scaling_fixture':[[str(a),str(b)]for a,b in fixture_values],
        'circuit_all170_nodes_checked':len(circuit['nodes'])==170,
        'universal_circuit_lower12_coefficient_and_upper7_Laurent_cut_proof':True,
        'mixed_bidegree':[8,8],'generic_DVR_and_irreducibility_proof_in_review_note':True,
        'associated_bundle_descent_regular_on_p0_open_only':True,
    }
    (output/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    result={
        'id':ID,'all_passed':True,
        'independent_physical_embedding_gamma_six_actions_and_normalizer':True,
        'restriction_determinant_and_all59_discriminant_coefficients_verified_independently':True,
        'full_unexpanded170_node_circuit_verified_as_coefficient_identity':True,
        'exact_extraneous_resultant_power2_accepted':True,
        'irreducible_mixed_D_bidegree8_8_total16_accepted':True,
        'normalizer_evenness_accepted':True,
        'associated_bundle_isomorphism_and_regular_invariant_descent_on_p0_open_accepted':True,
        'descended_rational_invariant_degree':16,
        'pole_order_k_bound_or_degree24_polynomial_claim':False,
        'NEW25_relation_or_vanishing_claim':False,
        'no_new_firstjet_point_or_highdegree_expansion':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,author,certificate,circuitpath,spinreview,spincert,original,
            ROOT/'code/q_firstjet_spin4_mixed_discriminant_20261103_v63.py',
            ROOT/'research/notes/Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63.md',
            ROOT/'research/notes/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW.md']},
        'checked_dependency_sha256':checked,
        'evidence_sha256':{'independent-certificate.json':sha(output/'independent-certificate.json')},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({'all_passed':True,'seconds':result['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
