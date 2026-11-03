#!/usr/bin/env python3
"""Independent full octic coefficients and same-point two-adjoint recovery.

Uses a different four-bit exponent encoding, own Fraction coefficient algebra,
and exact minors. No author module, new beta point, GB, or degree22+ probe.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import math
from pathlib import Path
import json
import resource
import signal
import time

from q_firstjet_spin4_mixed_discriminant_review_20261103_v63 import mat,mm,tp,det,inv,add

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW'


def read(p):return json.loads(gzip.decompress(p.read_bytes())if p.suffix=='.gz'else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def zero(n,m):return [[Q(0)for j in range(m)]for i in range(n)]
def pa(A,B,c=1):
    C=dict(A)
    for m,a in B.items():
        C[m]=C.get(m,0)+c*a
        if not C[m]:del C[m]
    return C
def pm(A,B):
    C=defaultdict(int)
    for m,a in A.items():
        for n,b in B.items():C[m+n]+=a*b
    return {m:a for m,a in C.items()if a}
def psum(values):
    out={}
    for p in values:out=pa(out,p)
    return out
def pt(A):return list(map(list,zip(*A)))
def pmm(A,B):
    return [[psum(pm(x,y)for x,y in zip(row,col))for col in zip(*B)]for row in A]
def pconstant(A):return [[{0:q}if q else{}for q in row]for row in A]
def ptrace(A):return psum(A[i][i]for i in range(len(A)))
def tracepair(A,B):return psum(pm(A[i][j],B[j][i])for i in range(4)for j in range(4))
def trace(A):return sum(A[i][i]for i in range(len(A)))
def from9(m):
    out=0;degree=0
    for i in range(16):
        e=m%9;m//=9;degree+=e;out+=e<<(4*i)
    assert m==0 and degree==8
    return out
def evaluate(poly,beta):
    answer=Q(0)
    for m,c in poly.items():
        term=Q(c)
        for i,b in enumerate(beta):term*=b**((m>>(4*i))&15)
        answer+=term
    return answer
def pivot_columns(A):
    A=[row[:]for row in A];r=0;out=[]
    for j in range(len(A[0])):
        i=next((i for i in range(r,len(A))if A[i][j]),None)
        if i is None:continue
        A[i],A[r]=A[r],A[i];p=A[r][j];A[r]=[x/p for x in A[r]]
        for i in range(r+1,len(A)):
            c=A[i][j];A[i]=[x-c*y for x,y in zip(A[i],A[r])]
        out.append(j);r+=1
        if r==len(A):break
    return out
def witness(M,record):
    sub=[[M[i][j]for j in record['columns']]for i in record['rows']]
    assert sub==mat(record['minor'])
    d=det(sub);assert d==Q(record['determinant'])and d
    return str(d)


def main(output):
    started=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    output.mkdir(parents=True,exist_ok=False)
    author=ROOT/'reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63/report.json'
    certpath=author.parent/'certificate.json.gz';ar=read(author);cert=read(certpath);assert ar['all_passed']
    spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
    spinreview=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63_REVIEW/report.json'
    pointpath=ROOT/'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63/certificate.json.gz'
    pointreview=ROOT/'reports/Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63_REVIEW/report.json'
    origpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
    spin,point,orig=map(read,[spinpath,pointpath,origpath]);assert read(spinreview)['all_passed']and read(pointreview)['all_passed']
    checked={}
    for path in [author,spinreview,pointreview]:
        for name,want in read(path)['input_sha256'].items():assert sha(ROOT/name)==want,name;checked[name]=want
    assert point['sample_index']==cert['same_point_index']==0
    beta=list(map(Q,point['accepted_source_point']['beta16']))
    J=mat([[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]);Ji=inv(J)
    metric=mat(spin['quadratic_form5']);B=mat(orig['primitive24_by16_basis'])
    pairs=list(itertools.combinations(range(4),2));coords=[(ij,k)for ij in pairs for k in range(4)]
    psi_poly=[[{1<<(4*l):B[coords.index((ij,k))][l]for l in range(16)if B[coords.index((ij,k))][l]}
               for k in range(4)]for ij in pairs[:-1]]
    # Reconstruct both quadratic covariants from the original physical embedding.
    Cfresh=pmm(pmm(pmm(pt(psi_poly),pconstant(metric)),psi_poly),pconstant(Ji))
    Vfresh=pmm(pmm(pmm(psi_poly,pconstant(Ji)),pt(psi_poly)),pconstant(metric))
    generators=[]
    for i in range(4):
        for j in range(i,4):
            symmetric=zero(4,4);symmetric[i][j]=symmetric[j][i]=1
            generators.append(mm(J,symmetric))
    At=list(map(tp,generators));assert At==list(map(mat,cert['pair_coordinate_basis4']))
    alpha=[]
    for i,j in pairs[:-1]:
        A=zero(4,4);A[i][j]=1;A[j][i]=-1
        if(i,j)==(0,1):A[2][3]=-1;A[3][2]=1
        alpha.append(A)
    vector_reps=[]
    for A in generators:
        columns=[]
        for E in alpha:
            value=add(mm(tp(A),E),mm(E,A))
            assert value[0][1]+value[2][3]==0
            columns.append([value[i][j]for i,j in pairs[:-1]])
        vector_reps.append(tp(columns))
    assert vector_reps==list(map(mat,spin['action5_matrices']))
    R5=[[vector_reps[k][i][j]for k in range(10)]for i in range(5)for j in range(5)]
    r5rows=pivot_columns(tp(R5));assert len(r5rows)==10
    R5left=inv([R5[i]for i in r5rows])
    Vflat=[p for row in Vfresh for p in row]
    Dcoordinates=[psum({m:c*q for m,q in Vflat[i].items()}for c,i in zip(row,r5rows)if c)for row in R5left]
    Dfresh=[[psum({m:At[k][i][j]*q for m,q in Dcoordinates[k].items()}for k in range(10)if At[k][i][j])for j in range(4)]for i in range(4)]
    assert all(psum({m:R5[r][k]*q for m,q in Dcoordinates[k].items()}for k in range(10)if R5[r][k])==Vflat[r]for r in range(25))
    quadratic=spin['quadratic_monomials'];encoded=[(1<<(4*i))+(1<<(4*j))for i,j in quadratic]
    stored=[]
    for key in ['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136']:
        stored.append([{m:Q(c)for m,c in zip(encoded,row)if Q(c)}for row in spin[key]])
    assert [p for row in Cfresh for p in row]==stored[0]
    assert [p for row in Dfresh for p in row]==stored[1]
    gamma=list(map(mat,spin['gamma_matrices']))
    assert all(psum({m:gamma[e][r][k]*c for m,c in psi_poly[e][k].items()}for e in range(5)for k in range(4)if gamma[e][r][k])=={}for r in range(4))
    def cleared(A):
        out=[]
        for row in A:
            rr=[]
            for poly in row:
                q={m:4*c for m,c in poly.items()};assert all(c.denominator==1 for c in q.values())
                rr.append({m:c.numerator for m,c in q.items()})
            out.append(rr)
        return out
    Cpol,Dpol=cleared(Cfresh),cleared(Dfresh)
    CC,CD,DD=pmm(Cpol,Cpol),pmm(Cpol,Dpol),pmm(Dpol,Dpol)
    traces=[tracepair(A,B)for A,B in [(CC,CC),(CC,CD),(CC,DD),(CD,CD),(CD,DD),(DD,DD)]]
    qpol=[ptrace(A)for A in [CC,CD,DD]]
    products=[pm(qpol[i],qpol[j])for i in range(3)for j in range(i,3)]
    polynomials=traces+products
    for actual,saved in zip(polynomials,cert['trace_polynomials']+cert['product_polynomials']):
        assert actual=={from9(m):c for m,c in saved}
    relations=cert['seven_complete_relation_rows'];assert len(relations)==7
    assert relations==[
        [0]*6+[1,0,0,0,0,-4], [0]*6+[0,1,0,0,0,2], [0]*6+[0,0,1,0,0,-2],
        [0]*6+[0,0,0,1,0,-1], [0]*6+[0,0,0,0,1,1],
        [1,0,0,0,0,4,0,0,0,0,0,-3],[-1,0,6,0,12,8,0,0,0,0,0,0]]
    for relation in relations:assert psum({m:c*a for m,a in poly.items()}for c,poly in zip(relation,polynomials)if c)=={}
    relation_det=det([[Q(row[j])for j in [6,7,8,9,10,0,2]]for row in relations]);assert relation_det==6
    selected=[from9(m)for m in cert['coefficient_selected_monomials']]
    coefficient_matrix=[[Q(p.get(m,0))for p in polynomials]for m in selected]
    coefficient_det=witness(coefficient_matrix,cert['coefficient_rank5_witness']);assert Q(coefficient_det)==8388608
    # All following numeric operations use the SAME accepted QQ beta[0].
    psi=[[evaluate(p,beta)for p in row]for row in psi_poly]
    C=[[evaluate(p,beta)for p in row]for row in Cfresh]
    D=[[evaluate(p,beta)for p in row]for row in Dfresh]
    V=[[evaluate(p,beta)for p in row]for row in Vfresh]
    assert psi==mat(cert['actual_psi5x4'])==mat(point['psi5x4'])
    assert C==mat(cert['actual_C4'])and D==mat(cert['actual_D4'])and V==mat(cert['actual_V5'])
    Cprime=mm(mm(Ji,C),J)
    assert mm(V,psi)==mm(psi,Cprime)
    L=zero(24,20)
    for col in range(20):
        E=zero(5,4);E[col//4][col%4]=1
        difference=add(mm(V,E),mm(E,Cprime),-1)
        values=[x for row in difference for x in row]
        values += [gamma[col//4][r][col%4]for r in range(4)]
        for r,value in enumerate(values):L[r][col]=value
    assert L==mat(cert['linear_recovery24x20'])
    psivec=[x for row in psi for x in row]
    assert all(sum(c*x for c,x in zip(row,psivec))==0 for row in L)
    recovery_det=witness(L,cert['recovery_rank19_witness'])
    free=cert['recovery_free_column'];assert psivec[free]
    normalized=[x/psivec[free]for x in psivec]
    assert normalized==list(map(Q,cert['recovery_normalized_kernel_at_same_point']))
    # Recover all20 pair coordinates and differentiate their full136 coefficients.
    G=[[A[i][j]for A in At]for i in range(4)for j in range(4)]
    grows=pivot_columns(tp(G));Gleft=inv([G[i]for i in grows]);assert len(grows)==10
    pairpoly=[]
    for fresh in [Cfresh,Dfresh]:
        flat=[p for row in fresh for p in row]
        pairpoly += [psum({m:c*q for m,q in flat[i].items()}for c,i in zip(row,grows)if c)for row in Gleft]
    paircoeff=[[poly.get(m,Q(0))for m in encoded]for poly in pairpoly]
    assert paircoeff==mat(cert['full_pair_coordinate_coefficients20x136'])
    Jac=zero(20,16)
    for k,(i,j)in enumerate(quadratic):
        for row in range(20):
            Jac[row][i]+=paircoeff[row][k]*beta[j]
            Jac[row][j]+=paircoeff[row][k]*beta[i]
    assert Jac==mat(cert['jacobian20x16'])
    jacobian_det=witness(Jac,cert['jacobian_rank16_witness'])
    def trword(word):
        value=[[Q(i==j)for j in range(4)]for i in range(4)]
        for c in word:value=mm(value,C if c=='C'else D)
        return trace(value)
    def derivative(word,letter,E):
        value=[[Q(i==j)for j in range(4)]for i in range(4)];jet=zero(4,4)
        for c in word:
            factor=C if c=='C'else D
            jet=add(mm(jet,factor),mm(value,E)if c==letter else zero(4,4))
            value=mm(value,factor)
        return trace(jet)
    q3=trword('DD');gradient=zero(4,20)
    assert trword('CC')-2*q3==trword('CD')+q3==0
    assert trword('CCCC')+4*trword('DDDD')-3*q3**2==0
    assert -trword('CCCC')+6*trword('CCDD')+12*trword('CDDD')+8*trword('DDDD')==0
    for i in range(20):
        letter='C'if i<10 else'D';E=At[i%10]
        dcc,dcd,ddd,dc4,dd4=[derivative(w,letter,E)for w in ['CC','CD','DD','CCCC','DDDD']]
        values=[dcc-2*ddd,dcd+ddd,dc4+4*dd4-6*q3*ddd,
                -dc4+6*derivative('CCDD',letter,E)+12*derivative('CDDD',letter,E)+8*dd4]
        for r,x in enumerate(values):gradient[r][i]=x
    assert gradient==mat(cert['four_constraints_gradient4x20'])and mm(gradient,Jac)==zero(4,16)
    gradient_det=witness(gradient,cert['constraint_gradient_rank4_witness'])
    nonzero_index=cert['scale_fixing_entry_in_C_then_D_row_major']
    assert ([x for row in C for x in row]+[x for row in D for x in row])[nonzero_index]==Q(cert['scale_fixing_nonzero_entry'])!=0
    # Independent complete symmetric-power DP, with all weight multiplicities.
    weights=list(map(tuple,orig['Cartan_weights']));assert len(weights)==16
    dp=[defaultdict(int)for d in range(25)];dp[0][(0,0)]=1
    for w in weights:
        for degree in range(1,25):
            for (x,y),count in dp[degree-1].items():dp[degree][(x+w[0],y+w[1])]+=count
    weyl=[]
    for perm in [(0,1),(1,0)]:
        for signs in itertools.product([-1,1],repeat=2):
            rho=(2,1);wrho=(signs[0]*rho[perm[0]],signs[1]*rho[perm[1]])
            weyl.append(((2-wrho[0],1-wrho[1]),(-1 if perm==(1,0)else 1)*signs[0]*signs[1]))
    chars=[]
    for degree in range(25):
        total=sum(dp[degree].values());assert total==math.comb(15+degree,degree)
        terms={w:(sign,dp[degree].get(w,0))for w,sign in weyl}
        value=sum(sign*mult for sign,mult in terms.values())
        saved=cert['character_records_through24'][degree]
        assert saved['degree']==degree and value==saved['invariant_dimension']
        assert terms=={tuple(t['weight']):(t['sign'],t['multiplicity'])for t in saved['Weyl_terms']}
        chars.append({'degree':degree,'total_weight_monomials':total,'invariant_dimension':value})
    assert chars[22]['invariant_dimension']==21 and chars[24]['invariant_dimension']==67
    evidence={
        'fresh_C_and_D_from_original20_by16_Psi_embedding':True,
        'all25_V_D_quadratic_coefficients_and_full_gamma_trace_verified':True,
        'independent_exponent_encoding':'Four bits per variable, no carry for totaldegree8',
        'all12_octic_coefficients_match':True,
        'trace_support_sizes':[len(p)for p in traces],
        'seven_relation_minor':str(relation_det),'coefficient_rank5_minor':coefficient_det,
        'same_point_index':0,'recovery_rank19_minor':recovery_det,'pair_jacobian_rank16_minor':jacobian_det,
        'constraint_gradient_rank4_minor':gradient_det,
        'all24_by20_recovery_all20_by16_Jac_all4_by20_gradient_entries_match':True,
        'universal_recovery_identity_from_associativity':True,
        'characters_through24':chars,
        'generic_sign_recovery_and_local_component_proof_in_review_note':True,
    }
    (output/'independent-certificate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    note=ROOT/'research/notes'/(ID+'.md')
    files=[Path(__file__),note,author,certpath,spinpath,spinreview,pointpath,pointreview,origpath,
           ROOT/'code/q_firstjet_two_adjoint_relations_20261103_v63.py',
           ROOT/'code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py',
           ROOT/'research/notes/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63.md']
    result={
        'id':ID,'all_passed':True,
        'fresh_physical_quadratic_covariants_and_all12_octic_coefficients_verified':True,
        'all12_octic_rank_exact5_and_complete_seven_relations':True,
        'universal_V_D_Psi_recovery_identity_accepted':True,
        'same_existing_point_only_rank19_recovery_rank16_Jac_and_rank4_gradient_verified':True,
        'generic_pair_map_degree_exact2_signs_only':True,
        'image_function_field_equals_sign_invariants_and_Sp4_invariant_fields_agree':True,
        'image_one_generically_smooth_codim4_component_of_four_equations':True,
        'all_locus_exhaustion_normality_regular_ring_generation_claim':False,
        'character_only_Inv22':21,'character_only_Inv24':67,
        'new_beta_points_or_degree22plus_evaluations':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},
        'checked_dependency_sha256':checked,
        'evidence_sha256':{'independent-certificate.json':sha(output/'independent-certificate.json')},
        'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-started,
    }
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    signal.alarm(0);print(json.dumps({'all_passed':True,'seconds':result['seconds']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'reports'/ID)
    main(parser.parse_args().output.resolve())
