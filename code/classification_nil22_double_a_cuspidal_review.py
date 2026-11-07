#!/usr/bin/env python3
"""Independent Fraction(epsilon) full-six-magnetic actual-source review.

No producer/CAS imports. Reconstructs the constant anchor frame, all-six map,
full150 rational family, strict10/projective350 and actual Euler Pfaffian.
"""
import ast
from fractions import Fraction as F
from itertools import combinations
import hashlib,json,resource,signal,time
from pathlib import Path
import classification_plane_rank_zero_rank_two_review as r
import classification_plane_high_m_rank_at_least_two_review as h

ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/nil22_double_a_cuspidal_review.json'
OUT='reports/CLASSIFICATION_NIL22_DOUBLE_A_CUSPIDAL_REVIEW_20261107.json'
N=11;ZERO=(0,)*N;EPS=r.RHO
PARENT=('h','x','y','z','w')+tuple(f'r{i}'for i in range(6))
TARGET=('h','q1','q2','p1','p2')+tuple(f'm{i}'for i in range(6))
OLDTARGET=('h','x1','x2','x3','x4')+tuple(f'm{i}'for i in range(6))


def atom(c=1,k=None):
    e=[0]*N
    if k is not None:e[k]=1
    return {tuple(e):r.RF(c)}if c else{}
def parse(text,names):
    def visit(node):
        if isinstance(node,ast.Constant)and type(node.value)is int:return atom(node.value)
        if isinstance(node,ast.Name):
            if node.id=='a':return atom(EPS+F(1,4))
            if node.id=='epsilon':return atom(EPS)
            return atom(1,names.index(node.id))
        if isinstance(node,ast.UnaryOp):return r.q.scale(visit(node.operand),-1 if isinstance(node.op,ast.USub)else 1)
        if isinstance(node,ast.BinOp):
            a,b=visit(node.left),visit(node.right)
            if isinstance(node.op,ast.Add):return r.q.add(a,b)
            if isinstance(node.op,ast.Sub):return r.q.add(a,r.q.scale(b,-1))
            if isinstance(node.op,ast.Mult):return r.q.mul(a,b)
            if isinstance(node.op,ast.Div):assert set(b)=={ZERO};return r.q.scale(a,r.RF(1)/b[ZERO])
            if isinstance(node.op,ast.Pow):
                assert set(b)<={ZERO};k=h.constant(b.get(ZERO,r.RF()));assert k.denominator==1 and k>=0
                out=atom()
                for _ in range(int(k)):out=r.q.mul(out,a)
                return out
        raise ValueError(ast.dump(node))
    return visit(ast.parse(str(text),mode='eval').body)
def subst(P,images):
    out={}
    for ex,c in P.items():
        term=atom(c)
        for k,power in enumerate(ex):
            for _ in range(power):term=r.q.mul(term,images[k])
        out=r.q.add(out,term)
    return out
def endpoint(P):return {e:r.evaluate(c,F(0))for e,c in P.items()if r.evaluate(c,F(0))}
def mulmat(A,B):return [[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]
def identity(n):return [[F(i==j)for j in range(n)]for i in range(n)]
def solve(A,b):
    n=len(A[0]);B=[list(row)+[v]for row,v in zip(A,b)];j0=0
    for j in range(n):
        k=next(k for k in range(j0,len(B))if B[k][j]);B[j0],B[k]=B[k],B[j0];v=B[j0][j];B[j0]=[x/v for x in B[j0]]
        for k in range(len(B)):
            if k!=j0:
                v=B[k][j];B[k]=[x-v*y for x,y in zip(B[k],B[j0])]
        j0+=1
    assert all(not any(row[:n])and not row[n]for row in B[n:])
    return [row[n]for row in B[:n]]


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2);signal.alarm(90)
    start=time.monotonic();job=json.loads((ROOT/JOB).read_text())
    for row in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    saved=json.loads((ROOT/'reports/CLASSIFICATION_NIL22_DOUBLE_A_CUSPIDAL_LIMIT_20261107.json').read_text())
    parent=json.loads((ROOT/'reports/CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json').read_text())['actual_family_power3']
    targetrow=json.loads((ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json').read_text())['cases'][0]
    assert targetrow['name']=='nil22_image_double_A'
    P=r.q.antisymmetric({tuple(map(int,ij.split(','))):parse(v,PARENT)for ij,v in parent['universal_Poisson_brackets'].items()})
    target=r.q.antisymmetric({tuple(map(int,ij.split(','))):parse(v,OLDTARGET)for ij,v in targetrow['full_six_magnetic_source'].items()})
    x=[atom(1,i)for i in range(N)];base=[r.q.add(r.q.scale(x[3],-4),r.q.scale(x[0],-3)),r.q.add(r.q.scale(x[3],-48),r.q.scale(x[4],-F(9,2)),r.q.scale(x[0],-54))]
    B=[[F(-4),F(0)],[F(-48),-F(9,2)]];Bi=r.q.inverse(B);ct=[F(-3),F(-54)]
    assert r.q.determinant(B)==18
    baseimages=[x[0],x[1],base[0],x[2],base[1]]+x[5:]
    Mp=[[endpoint(subst(P.get((i,j),{}),baseimages))for j in (2,4)]for i in (1,3)]
    Nt=[[r.q.add(*(r.q.scale(Mp[i][k],Bi[j][k])for k in range(2)))for j in range(2)]for i in range(2)]
    A=[]
    for i in range(2):
        equations=[];rhs=[]
        for j in range(2):
            T=target.get((i+1,j+3),{});mons=set(Nt[0][j])|set(Nt[1][j])|set(T)
            for e in sorted(mons):
                equations.append([h.constant(Nt[k][j].get(e,r.RF()))for k in range(2)]);rhs.append(h.constant(T.get(e,r.RF())))
        A.append(solve(equations,rhs))
    assert A==[[F(v)for v in row]for row in saved['actual_fiber_parent_to_target']]
    detA=r.q.determinant(A);assert detA==F(5,24192)
    G=[[F()for j in range(5)]for i in range(5)];G[0][0]=1
    for i in range(2):
        for j,index in enumerate((1,3)):G[i+1][index]=A[i][j]
        G[i+3][0]=-sum(a*b for a,b in zip(Bi[i],ct))
        for j,index in enumerate((2,4)):G[i+3][index]=Bi[i][j]
    Gi=r.q.inverse(G);assert mulmat(G,Gi)==identity(5)
    detG=r.q.determinant(G);assert detG==F(-5,435456)
    assert G==[[F(v)for v in row]for row in saved['whole_projective_parent_to_target']]
    mons=[(0,0),(0,3),(0,4),(3,3),(3,4),(4,4)]
    magmons=[(0,0),(0,2),(0,4),(2,2),(2,4),(4,4)]
    magcolumns=[]
    for i,j in magmons:
        value=r.q.scale(r.q.mul(baseimages[i],baseimages[j]),detA)
        column=[]
        for u,v in mons:
            e=[0]*N;e[u]+=1;e[v]+=1;column.append(h.constant(value.get(tuple(e),r.RF())))
        magcolumns.append(column)
    M=list(map(list,zip(*magcolumns)));Mi=r.q.inverse(M)
    assert M==[[F(v)for v in row]for row in saved['endpoint_magnetic_matrix']]
    detM=r.q.determinant(M);assert detM==F(15625,1909592158002053382144)and mulmat(M,Mi)==identity(6)
    rmap=[r.q.add(*(r.q.scale(x[5+j],c)for j,c in enumerate(row)))for row in Mi]
    assert rmap==[parse(v,TARGET)for v in saved['actual_parent_r_binding']]
    images=[r.q.add(*(r.q.scale(x[j],v)for j,v in enumerate(row)))for row in Gi]+rmap
    S={ij:subst(v,images)for ij,v in P.items()}
    curve=r.q.antisymmetric({ij:r.q.add(*(r.q.scale(S.get((a,b),{}),G[ij[0]][a]*G[ij[1]][b])for a in range(5)for b in range(5)))for ij in r.PAIRS})
    published=r.q.antisymmetric({tuple(map(int,ij.split(','))):parse(v,TARGET)for ij,v in saved['complete_regular_family'].items()})
    assert all(curve.get(ij,{})==published.get(ij,{})for ij in r.PAIRS)
    assert all(endpoint(curve.get(ij,{}))==target.get(ij,{})for ij in r.PAIRS)
    for i,j,k in combinations(range(5),3):
        assert not r.q.add(*(r.q.mul(curve.get((a,l),{}),r.q.der(curve.get((b,c),{}),l))for a,b,c in [(i,j,k),(j,k,i),(k,i,j)]for l in range(5)))
    fullJ,D=h.jacobi_and_divisor(curve,6);assert all(not j for j in fullJ)
    assert D==parse(saved['actual_projective_degeneracy_quintic'],TARGET)
    C=parse('(384*p1**2*p2+512*h*p1**2+384*h*p1*p2-9*h*p2**2-432*h**3)/384',TARGET)
    explicit=parse('5*h/(36*epsilon-5)',TARGET)
    expected=r.q.mul(explicit,r.q.add(r.q.mul(x[0],C),parse('epsilon*(4*p1+3*h)**4/6',TARGET)))
    assert D==expected and endpoint(D)==r.q.scale(r.q.mul(r.q.mul(x[0],x[0]),C),-1)
    Fparent=parse(parent['plane_quartic'],PARENT)
    C3=parse('h*w**2+12*w*y**2-144*y**3',PARENT)
    assert endpoint(Fparent)==r.q.mul(x[0],C3)
    assert subst(C3,baseimages)==r.q.scale(C,-864)
    denominator=r.RF((F(-960),F(6912)))
    assert denominator==parse(saved['family_denominator'],TARGET)[ZERO]
    for p in curve.values():
        for c in p.values():assert c.q[0]and not r.ud(denominator.p,c.q)[1]
    assert r.evaluate(parse(parent['coefficient_denominator'],PARENT)[ZERO],F(0))==F(-5,2)
    out=dict(task_id='E_multiple_hyperplane',all_passed=True,producer_or_CAS_imports=False,
             exact_input_bindings=job['exact_inputs_and_hashes'],known_component='R4_GA9_CUSPIDAL_PLANE_QUARTIC_24',
             constant_anchor_frame_independently_solved=True,actual_PGL_G=[[str(v)for v in row]for row in G],
             det_fiber=str(detA),det_PGL=str(detG),det_all6_magnetic_map=str(detM),
             full_all6_magnetic_inverse_checked=True,all150_endpoint_and_rational_family_checked=True,
             all10_strict_and350_projective_Jacobi_zero=True,actual_Euler_Pfaffian_full_family_checked=True,
             full_denominator_unit_at_endpoint=True,whole_source_PGL_closure_in_known_component=True,
             nil22_image_unrecognized_sources_remaining=['stable_line_conic','stable_three_lines'],
             new_component_or_global_exhaustion_claimed=False,
             elapsed_seconds=time.monotonic()-start,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ['all_passed','elapsed_seconds','peak_rss_kib']}))


if __name__=='__main__':main()
