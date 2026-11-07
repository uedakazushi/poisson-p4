#!/usr/bin/env python3
"""Independent all eight fixed-Q ENTIRE affine W modules, Fraction only.
Normal forms are whole actual-orbit representatives, not chosen higherpoints.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import hashlib,json,resource,signal,time
import classification_degenerate_quadric_q4_nil31_I_review as q
import classification_plane_rank_zero_nonkernel_review as H
from classification_minimal_quartic_review import LaurentRing
ROOT=Path(__file__).resolve().parents[1]
JOB='research/classification/jobs/plane_ore_module_review.json'
OUT='reports/CLASSIFICATION_PLANE_ORE_MODULE_REVIEW_20261107.json'
R=LaurentRing((),('a','b','c','z','w'))
def p(text):return {e+(0,0,0):c for e,c in R.parse(text).items()}
def encode(v):return {str(i):str(c)for i,c in enumerate(v)if c}
def Wtensor(W):return q.antisymmetric({(0,i):v for i,v in enumerate(W,1)})
def bases():
    out=[]
    for i in range(2):
        for mon in ('b*z','c*z','b*w','c*w'):
            W=[{}for _ in range(4)];W[i]=p(mon);out.append(Wtensor(W))
    for i in range(2):
        for mon in ('z**2','z*w','w**2'):
            W=[{}for _ in range(4)];W[i]=p(mon);out.append(Wtensor(W))
    for i,mon in ((2,'z*w'),(2,'w**2'),(3,'z**2'),(3,'z*w')):
        W=[{}for _ in range(4)];W[i]=p(mon);out.append(Wtensor(W))
    assert len(out)==18
    assert all(not any(q.jacobi(H.Padd(a,b)).values())for a,b in combinations(out,2))
    return out
def model(name,delta):
    if name=='A':
        bc='-z*b'+(' + w**2'if delta else'');xb=['0','z**2'];xc=['2*z**2','3*z*w'];W0=['0','0','0','b*z']
    elif name=='B':
        bc='z*b-2*w*c'+(' + z*w'if delta else'');xb=['-2*z*w','-3*w**2'];xc=['0','-z*w'];W0=['0','0','0','b*z']
    elif name=='S3':
        bc='w*c'+(' + z**2'if delta else'');xb=['3*z*w','2*w**2'];xc=['w**2','0'];W0=['0','0','4*b*z','2*b*w']
    else:
        assert name=='S1';bc='-z*c'+(' + z*w'if delta else'');xb=['-2*z**2','-z*w'];xc=['0','w**2'];W0=['0','0','2*b*z','0']
    Q=q.antisymmetric({(1,2):p(bc),(1,3):p(xb[0]),(1,4):p(xb[1]),(2,3):p(xc[0]),(2,4):p(xc[1])})
    return Q,H.Padd(Q,Wtensor([p(t)for t in W0]))
def case(name,delta,B):
    Q,P=model(name,delta);assert not any(q.jacobi(Q).values())
    const=q.jvec(q.jacobi(P));columns=[]
    for T in B:columns.append([a-b for a,b in zip(q.jvec(q.jacobi(H.Padd(P,T))),const)])
    M=q.transpose(columns);indices=[i for i,(quad,e)in enumerate(q.ROWS)if 0 in quad and e[0]==0];assert len(indices)==140
    assert all(not v for i,v in enumerate(const)if i not in indices)
    assert all(not v for i,row in enumerate(M)if i not in indices for v in row)
    rank,cols=q.rank(M);_,rows=q.rank(q.transpose(M));assert len(rows)==len(cols)==rank
    A=[[M[i][j]for j in cols]for i in rows];det=q.determinant(A);assert det
    inv=q.inverse(A);x=[F(0)for _ in B]
    for j,row in zip(cols,inv):x[j]=-sum((v*const[i]for i,v in zip(rows,row)),F(0))
    residual=[sum((v*w for v,w in zip(row,x)),F(0))+c for row,c in zip(M,const)]
    augmented=[row+[v]for row,v in zip(M,const)];ar,ac=q.rank(augmented)
    consistent=not any(residual);assert consistent==(ar==rank)
    kernel=[]
    for free in range(18):
        if free in cols:continue
        v=[F(0)for _ in B];v[free]=1
        for j,row in zip(cols,inv):v[j]=-sum((a*M[i][free]for i,a in zip(rows,row)),F(0))
        assert all(not sum((a*b for a,b in zip(row,v)),F(0))for row in M);kernel.append(v)
    base=5 if name in ('A','B') else 6
    orbit={'A':5,'B':6,'S3':5,'S1':5}[name]+delta
    k=18-rank;bound=base+orbit+k+5
    result=dict(anchor=name,central_delta=delta,Q_full150=encode(q.vec(Q)),baseline_full150=encode(q.vec(P)),
                full140_rows_in_full350=indices,complete140_by18=[[str(v)for v in M[i]]for i in indices],
                full350_constant=encode(const),rank=rank,rank_minor_rows=rows,rank_minor_columns=cols,rank_minor_determinant=str(det),
                full18_kernel=[[str(v)for v in row]for row in q.transpose(kernel)],affine_consistent=consistent,
                affine_particular=list(map(str,x))if consistent else None,full350_particular_residual=encode(residual),
                affine_fiber_dimension=k if consistent else None,affine_symbol_B=base,Q_orbit_dimension_upper=orbit,
                restored_projective_image_upper=bound if consistent else None,G24_closed=not consistent or bound<=24)
    if consistent:
        point=H.Padd(P,q.antisymmetric({ij:q.add(*(q.scale(T.get(ij,{}),v)for T,v in zip(B,x)))for ij in q.PAIRS}))
        D0=H.divisor(point);Ds=[]
        for vector in kernel:
            direction=q.antisymmetric({ij:q.add(*(q.scale(T.get(ij,{}),v)for T,v in zip(B,vector)))for ij in q.PAIRS})
            Ds.append(q.add(H.divisor(H.Padd(point,direction)),q.scale(D0,-1)))
        polynomials=[poly for poly in [D0]+Ds if poly]
        common=tuple(min(e[i]for poly in polynomials for e in poly)for i in range(5))if polynomials else None
        result['actual_D_affine_particular']={','.join(map(str,e[:5])):str(v)for e,v in D0.items()}
        result['actual_D_all_affine_kernel_coefficients']=[{','.join(map(str,e[:5])):str(v)for e,v in poly.items()}for poly in Ds]
        result['universal_monomial_divisor']=list(common)if common is not None else None
        result['whole_family_rank_le2_or_coordinate_factor']=not polynomials or any(common)
        if result['whole_family_rank_le2_or_coordinate_factor']:result['G24_closed']=True
    if not consistent:
        _,rr=q.rank(q.transpose(augmented));minor=q.determinant([[augmented[i][j]for j in ac]for i in rr]);assert minor
        result['inconsistency_minor']=dict(rows=rr,columns=ac,determinant=str(minor),rank=ar)
    return result
def main():
    started=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    assert not(ROOT/OUT).exists();job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    B=bases();records=[]
    for name in ('A','B','S3','S1'):
        for delta in (0,1):
            r=case(name,delta,B);records.append(r);print(name,delta,'rank',r['rank'],'k',r['affine_fiber_dimension'],'bound',r['restored_projective_image_upper'],'closed',r['G24_closed'],flush=True)
    result=dict(task_id='G',all_arithmetic_passed=True,all_passed=all(r['G24_closed']for r in records),
                all_eight_entire_affine_modules=records,full18_source_basis=[encode(q.vec(P))for P in B],
                all350_Jacobi_equivalent_to_true140_affine_rows=True,no_CAS_or_producer_imports=True,
                scope='Whole canonical Q orbit modules, with complete W18, no higher unknown localized. Fullsource/Ore/normalizer proof remains an explicit mathematical dependency; no proof inferred from an anchor alone.',
                inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    (ROOT/OUT).write_text(json.dumps(result,indent=2)+'\n');print('PASS arithmetic',OUT,round(result['seconds'],3),flush=True)
if __name__=='__main__':main()
