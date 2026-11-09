#!/usr/bin/env python3
"""Independent full150 strict-Casimir review of two-line m4 quintics.

The whole Jordan parameter is treated in Fraction(tau), its only minor wall
is treated in Fraction. No producer, producer helper or CAS is imported.
"""
import ast,functools,hashlib,json,resource,time
from fractions import Fraction as F
from pathlib import Path
import classification_plane_m4_special_review as q
from classification_unimodular_logarithmic_review import RF

ROOT=Path(__file__).resolve().parents[1]
INPUT='reports/CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_20261107.json'
JOB='research/classification/jobs/plane_m4_two_line_casimir_review.json'
OUTPUT='reports/CLASSIFICATION_PLANE_M4_TWO_LINE_CASIMIR_REVIEW_20261107.json'
T=RF((F(0),F(1)))


def compact(v):
    if isinstance(v,RF)and v.q==(F(1),)and len(v.p)<=1:return v.p[0]if v.p else F()
    return v


@functools.lru_cache(None)
def scalar(text):
    def walk(n):
        if isinstance(n,ast.Constant):return F(n.value)
        if isinstance(n,ast.Name):assert n.id=='tau';return T
        if isinstance(n,ast.UnaryOp):return -walk(n.operand)if isinstance(n.op,ast.USub)else walk(n.operand)
        a=walk(n.left)
        if isinstance(n.op,ast.Pow):
            out=F(1)
            for _ in range(n.right.value):out=compact(out*a)
            return out
        b=walk(n.right)
        if isinstance(n.op,ast.Add):return compact(a+b)
        if isinstance(n.op,ast.Sub):return compact(a-b)
        if isinstance(n.op,ast.Mult):return compact(a*b)
        if isinstance(n.op,ast.Div):return compact(a/b)
        raise ValueError(ast.dump(n))
    return walk(ast.parse(str(text),mode='eval').body)


def determinant(M):
    B=[list(r)for r in M];out=F(1)
    for j in range(len(B)):
        k=next((k for k in range(j,len(B))if B[k][j]),None)
        if k is None:return F()
        if k!=j:B[j],B[k]=B[k],B[j];out=-out
        pivot=B[j][j];out=compact(out*pivot)
        active=[l for l in range(j+1,len(B))if B[j][l]]
        for k in range(j+1,len(B)):
            if B[k][j]:
                c=compact(B[k][j]/pivot)
                for l in active:B[k][l]=compact(B[k][l]-c*B[j][l])
                B[k][j]=F()
    return out


def source_matrix(f,with_divergence=True):
    gradient=[q.der(f,i)for i in range(5)];mon6=q.monomials(6);mon1=q.monomials(1)
    columns=[]
    for(i,j),e in q.COLS:
        m={e:F(1)};cas=[{}for _ in range(5)];dv=[{}for _ in range(5)]
        cas[i]=q.mul(m,gradient[j]);cas[j]=q.scale(q.mul(m,gradient[i]),-1)
        if with_divergence:dv[i]=q.scale(q.der(m,j),-1);dv[j]=q.der(m,i)
        columns.append([cas[a].get(t,F())for a in range(5)for t in mon6]+[dv[a].get(t,F())for a in range(5)for t in mon1])
    return q.transpose(columns)


def exact_matrix(saved):return [[scalar(v)for v in row]for row in saved]


def verify_complete_source(M,c):
    rank=c['rank'];rows,cols=c['minor_rows'],c['minor_columns'];assert len(rows)==len(cols)==rank
    assert len(set(rows))==len(set(cols))==rank
    value=determinant([[M[i][j]for j in cols]for i in rows]);assert value==scalar(c['minor_determinant'])and value
    K=exact_matrix(c['complete_kernel']);free=[j for j in range(150)if j not in cols]
    assert len(free)==c['nullity']and rank+len(free)==150
    assert [K[j]for j in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
    for row in M:
        entries=[(i,v)for i,v in enumerate(row)if v]
        for k in range(len(free)):
            assert not sum((v*K[i][k]for i,v in entries if K[i][k]),F())
    base_columns=[at for at,(ij,e)in enumerate(q.COLS)if 4 not in ij];assert len(base_columns)==90
    assert all(not v for i in base_columns for v in K[i])
    return dict(rank=rank,nullity=len(free),minor_determinant=repr(value),kernel_free_coordinate_indices=free,
        full150_kernel=[[repr(v)for v in row]for row in K],all_original_rows_annihilate_kernel=True,all90_baseblock_rows_zero=True,common_vertex_exterior_factor=True,rank4_source_empty=True),value,K


def main():
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    data=json.loads((ROOT/INPUT).read_text());assert data['all_exact_arithmetic_passed']
    assert hashlib.sha256(data['producer_bytes'].encode()).hexdigest()==data['producer_sha256']
    assert data['all150_quadratic_basis']==[[i,j,list(e)]for(i,j),e in q.COLS]
    assert data['all210_degree6_monomials']==list(map(list,q.monomials(6)))
    sem=q.add(q.a(1,0,2,2,2,2),q.a(1,1,2,2,3,3),q.a(1,3,3,3,3,3))
    jordan=q.add(q.a(1,0,2,2,3,3),q.a(1,1,2,2,2,3),q.a(1,2,2,2,2,2),q.a(1,3,3,3,3,3))
    perturb=q.a(1,2,3,3,3,3);records=[]
    for row,base,variable in zip(data['cases'],[sem,q.add(sem,perturb),jordan],[False,False,True]):
        M0=source_matrix(base);Mt=source_matrix(perturb,False)if variable else None
        M=[[compact(a+T*b)if b else a for a,b in zip(r,s)]for r,s in zip(M0,Mt)]if variable else M0
        assert M==exact_matrix(row['full1075by150_Casimir_and_divergence_matrix'])
        c,value,K=verify_complete_source(M,row['full_canonical_source_certificate']);c['type']=row['type']
        if variable:
            val=RF(value);assert val.q==(F(1),)and len(val.p)==7 and all(not v for v in val.p[:6])and val.p[6]
            assert row['generic_minor_all_factors']==[['tau',6]]
            assert len(row['all_exceptional_fields'])==1
            wall=row['all_exceptional_fields'][0];assert wall['minimal_polynomial']=='tau'and F(wall['tau'])==0
            assert wall['all_complex_embeddings_covered']
            wc,_,_=verify_complete_source(M0,wall['source_certificate']);c['all_parameter_walls']=[dict(tau='0',**wc)]
            # Every coefficient denominator is a power of tau, whose full wall
            # was just checked independently. No isolated complex root omitted.
            for line in K:
                for v in line:
                    den=RF(v).q;assert sum(bool(a)for a in den)==1
            c['generic_minor_only_zero_tau0']=True
            c['whole_complex_parameter_covered']=True
        else:assert not row['generic_minor_all_factors']and not row['all_exceptional_fields']
        records.append(c);print(row['type'],'full1075/150/kernel/vertex90 PASS',flush=True)
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all_highrank_generic211_two_line_sources_empty=True,cases=records,
        proof='All actualD normalforms and full Jordan line are necessary supersources for every support<=2line highrank211 m4 stratum. Canonical strict Casimir/div0 kernels have only vertex-paired brackets, so every such tensor is partial_w wedge V and rank<=2, contradicting rank4.',
        scope='Only m4 highrank generic211 with zero-cubic support atmost2lines. Primitive whole actualD normalforms/irreducibility/shear/nonzero complements read separately. m2/m3 other sources unaffected.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),no_CAS_or_producer_imports=True,seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/OUTPUT).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,whole_source_empty=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


def _complete_kernel_without_vertex_assumption(M,c):
    rows,cols=c['minor_rows'],c['minor_columns'];rank=c['rank']
    assert len(rows)==len(cols)==rank and len(set(rows))==len(set(cols))==rank
    value=determinant([[M[i][j]for j in cols]for i in rows]);assert value==scalar(c['minor_determinant'])and value
    K=exact_matrix(c['complete_kernel']);free=[j for j in range(150)if j not in cols]
    assert len(free)==c['nullity']and rank+len(free)==150
    assert [K[i]for i in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
    for row in M:
        entries=[(i,v)for i,v in enumerate(row)if v]
        for k in range(len(free)):assert not sum((v*K[i][k]for i,v in entries if K[i][k]),F())
    return dict(rank=rank,nullity=len(free),minor_determinant=repr(value),kernel_free_coordinate_indices=free,all_original_rows_annihilate_kernel=True),K


def _all_five_pfaffians(P):
    out=[]
    for missing in range(5):
        i,j,k,l=[v for v in range(5)if v!=missing]
        out.append(q.scale(q.add(q.mul(P.get((i,j),{}),P.get((k,l),{})),q.scale(q.mul(P.get((i,k),{}),P.get((j,l),{})),-1),q.mul(P.get((i,l),{}),P.get((j,k),{}))),(-1)**missing))
    return out


def squarefree_main():
    """Complete fixed six-source review; no common-vertex shortcut assumed."""
    import itertools
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(230,235));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    jobpath='research/classification/jobs/plane_m4_nonreduced_squarefree_review.json'
    inputpath='reports/CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_CASIMIR_20261107.json'
    outputpath='reports/CLASSIFICATION_PLANE_M4_NONREDUCED_SQUAREFREE_REVIEW_20261107.json'
    job=json.loads((ROOT/jobpath).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    data=json.loads((ROOT/inputpath).read_text());assert data['all_exact_arithmetic_passed']
    assert hashlib.sha256(data['producer_bytes'].encode()).hexdigest()==data['producer_sha256']
    assert data['all150_quadratic_basis']==[[i,j,list(e)]for(i,j),e in q.COLS]
    assert data['all210_degree6_monomials']==list(map(list,q.monomials(6)))
    secant=q.add(q.a(1,0,2,2,2,2),q.a(1,1,3,3,3,3))
    p=q.a(1,2,2,2,3,3);r=q.a(1,2,2,3,3,3)
    compression=q.add(q.a(1,0,2,2,2,3),q.a(1,1,2,3,3,3),q.a(1,2,2,2,2,2),q.a(1,3,3,3,3,3))
    mixed=q.add(q.a(1,0,2,2,2,2),q.a(1,1,2,3,3,3),q.a(1,3,3,3,3,3))
    Fs=[secant,q.add(secant,p),q.add(secant,p,r),compression,mixed,q.add(mixed,p)]
    records=[]
    assert len(data['cases'])==len(Fs)==6
    for row,f,dim in zip(data['cases'],Fs,[5,3,2,1,1,1]):
        M=source_matrix(f);assert M==exact_matrix(row['full1075by150_Casimir_and_divergence_matrix'])
        cert,K=_complete_kernel_without_vertex_assumption(M,row['full_canonical_source_certificate']);assert cert['nullity']==dim
        basis=[]
        for v in q.transpose(K):
            P={ij:{}for ij in q.PAIRS}
            for(ij,e),c in zip(q.COLS,v):
                if c:P[ij][e]=c
            basis.append(q.skew(P))
        pairs=list(itertools.combinations_with_replacement(range(dim),2))
        for i,j in pairs:
            current=basis[i]if i==j else q.plus(basis[i],basis[j])
            assert not any(_all_five_pfaffians(current))
        assert row['entire_necessary_source_rank_at_most2']
        records.append(dict(type=row['type'],F_sparse=[dict(exponent=e,coefficient=str(c))for e,c in sorted(f.items())],complete_source_certificate=cert,full150_kernel=[[str(v)for v in line]for line in K],all_quadratic_parameter_pairs=pairs,all_five_principal_Pfaffian_components_zero_for_every_pair=True,entire_necessary_source_rank_at_most2=True,rank4_source_empty=True))
        print(row['type'],'full1075/150/kernel/all5Pf-polarizations PASS',flush=True)
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all_six_highrank_special_squarefree_sources_empty=True,whole_nonreduced_m4_squarefree_source_covered_by_essential3_eight=True,cases=records,
        geometric_review='Nonexceptional J2/I3 invariance yields Delta=cL6. Root cover is etale on all P2\\L=A2 regardless of additional original connection poles; PicA2/unit determinants/base-degree<=1 force constant four roots. Special stable Witt/Hankel pencils and support<=2line give exactly six actualD normalforms after all shear/scaling zero-nonzero splits. Full kernels have all5Pf polarizations zero, so highrank sources empty; rank1 actualD essential3 uses accepted eight-source theorem.',
        scope='Only nonreducedR, m4 generic squarefree leading quartic: all points in eight essential3 components, G<=23. Other repeated-leading packets and m2/m3/E remain distinct.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),no_CAS_or_producer_imports=True,seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/outputpath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,highrank_empty=True,whole_SF_essential3_recognition=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


if __name__=='__main__':
    import sys
    if '--squarefree'in sys.argv:squarefree_main()
    else:main()
