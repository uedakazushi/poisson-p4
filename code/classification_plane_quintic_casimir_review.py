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


def reduced_m3_main():
    """Independent full350/24/150 replay of all nine reduced-m3 whole sources."""
    from itertools import product
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jobpath='research/classification/jobs/plane_reduced_m3_nine_separate_review.json'
    inputpath='reports/CLASSIFICATION_PLANE_REDUCED_M3_NINE_FULLHIGHER_20261108.json'
    outputpath='reports/CLASSIFICATION_PLANE_REDUCED_M3_NINE_SEPARATE_REVIEW_20261108.json'
    assert not(ROOT/outputpath).exists();job=json.loads((ROOT/jobpath).read_text())
    for row in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    data=json.loads((ROOT/inputpath).read_text());assert data['all_passed']and len(data['cases'])==9
    assert hashlib.sha256(data['producer_source'].encode()).hexdigest()==data['producer_sha256']
    def parse(text,names):
        def walk(n):
            if isinstance(n,ast.Constant):return q.a(F(n.value))
            if isinstance(n,ast.Name):return q.a(1,names.index(n.id))
            if isinstance(n,ast.UnaryOp):return q.scale(walk(n.operand),-1 if isinstance(n.op,ast.USub)else 1)
            a=walk(n.left)
            if isinstance(n.op,ast.Pow):
                out=q.a(1)
                for _ in range(n.right.value):out=q.mul(out,a)
                return out
            b=walk(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            if isinstance(n.op,ast.Div):
                zero=(0,)*5;assert set(b)=={zero};return q.scale(a,1/b[zero])
            raise ValueError(ast.dump(n))
        return walk(ast.parse(str(text),mode='eval').body)
    def exact(M):return[[F(v)for v in row]for row in M]
    def verify(M,C,n):
        r=C['rank'];rr,cc=C['rows'],C['columns'];assert len(M[0])==n and len(rr)==len(cc)==r
        assert len(set(rr))==len(set(cc))==r
        minor=q.det([[M[i][j]for j in cc]for i in rr]);assert minor==F(C['minor_determinant'])and minor
        K=exact(C['complete_kernel']);assert len(K)==n-r==C['nullity']
        free=[j for j in range(n)if j not in cc]
        assert [[row[j]for j in free]for row in K]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
        for row in M:
            nz=[(i,v)for i,v in enumerate(row)if v]
            assert all(not sum((c*v[i]for i,c in nz),F())for v in K)
        return K,str(minor)
    first,second=q.sources();assert len(first)==30 and len(second)==9
    mon3=[e for e in q.monomials(5)if e[3]+e[4]==3];assert len(mon3)==24
    X=[q.a(1,i)for i in range(5)]
    eulers=[q.vector(q.skew({(i,j):q.add(q.a(1,i,b)if j==a else{},q.a(-1,j,b)if i==a else{})for i,j in q.PAIRS}))for a,b in product(range(5),repeat=2)]
    assert len(q.rref(q.transpose(eulers))[1])==24
    levi=[(a,b)for a,b in product(range(5),repeat=2)if(a<3)==(b<3)]
    inputs={}
    for item in data['source_input_files']:
        name=item['path'];assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==item['sha256']
        packet=json.loads((ROOT/name).read_text());inputs[name]=packet
    records=[]
    for c in data['cases']:
        packet=inputs[c['input']];names=packet['variables'];profile=next(v for v in packet['profiles']if v['name']==c['name'])
        P=q.skew({tuple(names.index(s)for s in key.split(',')):parse(text,names)for key,text in profile['raw_brackets'].items()})
        assert q.vector(P)==list(map(F,c['full150_P0']))and not any(q.jac(P).values())and not q.actual(P)
        J1=q.transpose([q.jvector(q.jac(q.plus(P,V)))for V in first]);J2=q.transpose([q.jvector(q.jac(q.plus(P,V)))for V in second])
        assert J1==exact(c['complete350x30'])and J2==exact(c['complete350x9'])
        K0,minor0=verify(J1,c['weight1_unconstrained_certificate'],30);K2,minor2=verify(J2,c['weight2_certificate'],9)
        cross=[q.actual(q.plus(P,V))for V in first]
        assert all(sum(e[3:])==3 for p in cross for e in p)
        D3=[[p.get(e,F())for p in cross]for e in mon3];assert D3==exact(c['full24x30_actual_D3'])
        texts=profile.get('whole_parallel_G3_basis',profile.get('whole_leading_target_basis'))
        assert texts is not None
        T=[[parse(text,names).get(e,F())for e in mon3]for text in texts];assert T==exact(c['full_parallel_target_rows'])
        ann,tm=verify(T,c['target_certificate'],24)
        assert ann==exact(c['complete_target_annihilator'])and ann==q.rref(T)[2]
        constrained=J1+[[sum(a*b for a,b in zip(row,col))for col in zip(*D3)]for row in ann]
        assert constrained==exact(c['constrained_weight1_matrix'])
        K1,minor1=verify(constrained,c['constrained_weight1_certificate'],30)
        orbit=[]
        for a,b in levi:
            O=q.skew({(i,j):q.add(P.get((b,j),{})if i==a else{},P.get((i,b),{})if j==a else{},q.scale(q.mul(X[b],q.der(P.get((i,j),{}),a)),-1))for i,j in q.PAIRS})
            orbit.append(q.vector(O))
        G=q.transpose(eulers+orbit+[q.vector(P)]);assert G==exact(c['full150x39_Euler_Levi_scale'])
        KG,mg=verify(G,c['filtered_group_certificate'],39)
        B=c['filtered_group_certificate']['rank']-24;assert B==c['first_symbol_dimension']==c['proof_upper']==profile['first_symbol_source_dimension_upper']
        badmon=(0,0,0,0,2)
        ufactor=all(e[3]for ij,p in P.items()if 3 in ij for e in p)
        for row in K1:
            Q=q.skew({ij:q.add(*(q.scale(V.get(ij,{}),a)for V,a in zip(first,row)if a))for ij in q.PAIRS})
            ufactor=ufactor and all(not Q.get((i,3),{}).get(badmon,F())for i in range(3))
        assert ufactor==c['entire_bad_u_row_zero']==c['actual_D_has_proper_u_factor']
        higher=len(K1)+len(K2);upper=B+higher+5
        assert higher==c['full_higher_dimension_upper']and upper==c['whole_projective_image_upper']
        assert ufactor or upper<=24;assert c['G24_closed']
        records.append(dict(name=c['name'],full150_P0_and350_Jacobi_and_actualPf_zero=True,
            full350x30_and350x9_reconstructed=True,full24_D3_map_and_all_targets_reconstructed=True,
            target_rank=c['target_certificate']['rank'],target_complete_annihilator_dimension=len(ann),target_minor=tm,
            unconstrained_K1=len(K0),constrained_K1=len(K1),K2=len(K2),J1_minor=minor1,J2_minor=minor2,
            full150x39_filtered_group_reconstructed=True,filtered_group_rank=c['filtered_group_certificate']['rank'],filtered_group_minor=mg,
            complete_original_kernels_checked=True,first_symbol_B=B,higher_necessary_dimension_upper=higher,projective_image_upper=upper,
            proper_actual_u_factor_for_entire_source=ufactor,G24_closed=True))
        print(c['name'],'independent350/24/150/kernel/factor PASS',flush=True)
    assert [c['name']for c in records]==['S','N','D','M','SF_C3_UNUSED','TRIPLE_HX2','TRIPLE_HY2','TRIPLE_HYZ','TRIPLE_HXY']
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all_nine_whole_necessary_sources_G24_closed=True,cases=records,
        mathematical_bridge='AllfullP1 lies in verified target-constrained kerJ1. Over each P1, fullP2 is empty or a kerJ2 affine-torsor subset. Euler24/Levi13/scale gives B=rankG-24 and wholeprojective image<=B+k1+k2+5. A fullbad-u-row zero gives actualu factor for every solution, independently of quadratic compatibility.',
        scope='Only arithmetic and whole necessary-source dimension/factor implication. Exhaustive root-cover geometry and allpoint E recognition are separate reviews. No producer module was imported or executed; no CAS or numerical arithmetic.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),
        immutable_producer_snapshot=inputpath+':producer_source',seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/outputpath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,all9_G24_closed=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


def moving_double_m3_main():
    """Independent whole eleven-source and three Ore-module reconstruction."""
    import itertools,signal
    from classification_minimal_quartic_review import LaurentRing
    import classification_plane_double_line_m3_ore_review as L
    Qmake=_moving_double_Qmake;Wmake=_moving_double_Wmake
    begun=time.monotonic();signal.alarm(240)
    resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jobpath='research/classification/jobs/plane_nonreduced_m3_moving_double_separate_review.json'
    outputpath='reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_SEPARATE_REVIEW_20261108.json'
    job=json.loads((ROOT/jobpath).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    paths=['reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_20261107.json','reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_20261107.json']
    first,ore=[json.loads((ROOT/p).read_text())for p in paths]
    for d in[first,ore]:assert d['all_exact_arithmetic_passed']and hashlib.sha256(d['producer_bytes'].encode()).hexdigest()==d['producer_sha256']
    names=('a','b','c','z','w');Y=(0,1,3,4);zero=(0,)*5;X=[q.a(1,i)for i in range(5)]
    @functools.lru_cache(None)
    def parse(text,at=None):
        def walk(n):
            if isinstance(n,ast.Constant):return{zero:F(n.value)}if n.value else{}
            if isinstance(n,ast.Name):
                if n.id=='tau':return{zero:T if at is None else at}
                return X[names.index(n.id)]
            if isinstance(n,ast.UnaryOp):return q.scale(walk(n.operand),-1)if isinstance(n.op,ast.USub)else walk(n.operand)
            a=walk(n.left)
            if isinstance(n.op,ast.Pow):
                out={zero:F(1)}
                for _ in range(n.right.value):out=q.mul(out,a)
                return out
            b=walk(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            assert isinstance(n.op,ast.Div)and set(b)=={zero};return q.scale(a,1/b[zero])
        return walk(ast.parse(str(text),mode='eval').body)
    def sc(v,at=None):
        if at is None:return scalar(str(v))
        u=RF(scalar(str(v)));ev=lambda pp:sum((c*at**i for i,c in enumerate(pp)),F());return ev(u.p)/ev(u.q)
    def mat(v,at=None):return[[sc(t,at)for t in row]for row in v]
    def vec(v,at=None):return[sc(t[0]if isinstance(t,list)else t,at)for t in v]
    def savedP(v,at=None):return q.skew({ast.literal_eval(ij):parse(text,at)for ij,text in v.items()})
    def op(P,V):
        base=q.jvector(q.jac(P));selfq=q.jvector(q.jac(V));t=q.jvector(q.jac(q.plus(P,V)))
        return[compact(a-b-c)for a,b,c in zip(t,base,selfq)]
    def cert(M,C,at=None):
        rr,cc=C['minor_rows'],C['minor_columns'];r=C['rank'];n=len(M[0]);assert len(rr)==len(cc)==r and len(set(rr))==len(set(cc))==r
        d=determinant([[M[i][j]for j in cc]for i in rr]);assert d and d==sc(C['minor_determinant'],at)
        K=mat(C['complete_kernel'],at);free=[j for j in range(n)if j not in cc];assert len(free)==C['nullity']and len(K)==n
        assert[K[i]for i in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
        for row in M:
            nz=[(i,a)for i,a in enumerate(row)if a]
            for k in range(len(free)):assert not sum((a*K[i][k]for i,a in nz if K[i][k]),F())
        return K,d
    def rftext(c):
        if not isinstance(c,RF):return str(c)
        poly=lambda vv:'+'.join('('+str(a)+')*tau**'+str(i)for i,a in enumerate(vv)if a)or'0'
        return '('+poly(c.p)+')/('+poly(c.q)+')'
    def evalP(P,t):
        return{e:sc(rftext(c),t)for e,c in P.items()if sc(rftext(c),t)}
    def shrink(p,roots):
        a=list(RF(p).p)
        for r in roots:
            while len(a)>1 and not sum((c*r**i for i,c in enumerate(a)),F()):
                b=[F()]* (len(a)-1);b[-1]=a[-1]
                for i in range(len(a)-2,0,-1):b[i-1]=a[i]+r*b[i]
                assert a[0]+r*b[0]==0;a=b
        assert len(a)<=1 and a,('uncovered polynomial roots',a,roots)
    def cover_denominators(objects,roots):
        for v in objects:
            u=RF(v);shrink(RF(u.q),roots)
    def dividepoly(P,D):
        assert D;R=P.copy();lead=max(D);coef=D[lead]
        while R:
            e=max(R);assert all(a>=b for a,b in zip(e,lead)),('factor mismatch',P,D)
            exp=tuple(a-b for a,b in zip(e,lead));c=R[e]/coef
            R=q.add(R,q.scale(q.mul({exp:c},D),-1))
        return True
    B1,B2=q.sources();B1=B1[:18]+B1[24:27]+B1[18:21]+B1[27:30]+B1[21:24]
    base2=[e[:3]for e in q.monomials(2)if not any(e[3:])];base3=[e[:3]for e in q.monomials(3)if not any(e[3:])]
    mon3=[e+(j,3-j)for e in base2 for j in range(4)];mon2=[e+(j,2-j)for e in base3 for j in range(3)]
    assert list(map(list,mon3))==first['all24_D3_monomials']and list(map(list,mon2))==first['all30_K_monomials']
    R=LaurentRing(('tau',),names,invert=('a','b'));rp=R.parse
    def reconstruct_forms(c,at=None):
        fs='a**3'if c['type'].startswith('TRIPLE')else'a**2*b';f=parse(fs,at)
        P={(0,1):q.der(f,2),(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)}
        for form in c['full_connection_inverse_forms']:
            i,j=form['i'],form['j'];beta=[rp(t)for t in form['beta']];om=[rp(t)for t in form['horizontal_polynomial_form']]
            assert om==[R.add(R.scale(R.mul(rp(fs),v),3),R.scale(R.derivative(rp(fs),k),-int(i==j)))for k,v in enumerate(beta)]
            assert not R.add(*(R.mul(rp(names[k]),om[k])for k in range(3)))
            om=[parse(t,at)for t in form['horizontal_polynomial_form']]
            V=[q.scale(q.add(q.der(om[2],1),q.scale(q.der(om[1],2),-1)),F(1,3)),q.scale(q.add(q.der(om[0],2),q.scale(q.der(om[2],0),-1)),F(1,3)),q.scale(q.add(q.der(om[1],0),q.scale(q.der(om[0],1),-1)),F(1,3))]
            assert V==[parse(t,at)for t in form['inverse_tracefree_field']]
            inv=[q.add(q.mul(X[2],V[1]),q.scale(q.mul(X[1],V[2]),-1)),q.add(q.mul(X[0],V[2]),q.scale(q.mul(X[2],V[0]),-1)),q.add(q.mul(X[1],V[0]),q.scale(q.mul(X[0],V[1]),-1))]
            assert inv==om
            for k in range(3):P[k,3+i]=q.add(P.get((k,3+i),{}),q.scale(q.mul(V[k],X[3+j]),-1))
        return q.skew(P)
    def shear(P):
        cols=[]
        for normal in(3,4):
            for k in range(3):
                G=[{}for _ in range(5)];G[k]=X[normal]
                O=q.skew({(i,j):q.add(*(q.scale(q.mul(G[l],q.der(P.get((i,j),{}),l)),-1)for l in range(5)),*(q.add(q.mul(P.get((l,j),{}),q.der(G[i],l)),q.mul(P.get((i,l),{}),q.der(G[j],l)))for l in range(5)))for i,j in q.PAIRS})
                S=[q.add(q.scale(X[3],O.get((i,3),{}).get((0,0,0,2,0),F())),q.scale(X[4],O.get((i,4),{}).get((0,0,0,0,2),F())))for i in range(3)]+[{},{}]
                O=q.plus(O,q.wedge(X,S));raw=q.vector(O)
                out=[]
                for V in B1:
                    loc=next((ij,e)for ij,p in V.items()if ij in q.PAIRS for e in p);out.append(O.get(loc[0],{}).get(loc[1],F()))
                assert raw==q.vector(q.skew({ij:q.add(*(q.scale(V.get(ij,{}),a)for V,a in zip(B1,out)if a))for ij in q.PAIRS}));cols.append(out)
        return q.transpose(cols)
    first_results=[];first_kernel={};whole_minor_polys=[]
    for parent in first['cases']:
        genericP=reconstruct_forms(parent)
        variants=[(parent,None)]+[(dict(w,type=parent['type']),F(w['alpha']))for w in parent.get('all_exceptional_sources',[])]
        roots=[F(w['alpha'])for w in parent.get('all_exceptional_sources',[])]
        if roots:shrink(sc(parent['all_exception_polynomial']),roots)
        for c,t in variants:
            Pbase=genericP if t is None else q.skew({ij:evalP(p,t)for ij,p in genericP.items()})
            P=savedP(c['full150_P0'],t);Kp=parse(c['unique_global_K'],t)
            assert all(P.get(ij,{})==Pbase.get(ij,{})for ij in q.PAIRS if ij!=(3,4));assert P[3,4]==Kp
            D=q.actual(Pbase);KM=[]
            for nm in[(3,3),(3,4),(4,4)]:
                DD=q.actual(q.plus(Pbase,q.skew({(3,4):q.a(1,*nm)})))
                KM.append([compact(DD.get(e,F())-D.get(e,F()))for e in mon2])
            KM=q.transpose(KM);constant=[D.get(e,F())for e in mon2]
            assert KM==mat(c['full30_global_K_matrix'],t)and constant==vec(c['full30_global_K_constant'],t)
            _,kd=cert(KM,c['K_certificate'],t);_,kad=cert([r+[-a]for r,a in zip(KM,constant)],c['K_augmented_certificate'],t)
            assert not q.actual(P)and not any(q.jac(P).values())
            J1=q.transpose([op(P,V)for V in B1]);J2=q.transpose([op(P,V)for V in B2])
            assert J1==mat(c['full350_weight1_operator'],t)and J2==mat(c['full350_weight2_operator'],t)
            D3=q.transpose([[q.actual(q.plus(P,V)).get(e,F())for e in mon3]for V in B1]);assert D3==mat(c['full24by30_actual_D3_map'],t)
            target=parse(parent['exact_nonzero_leading_target'],t);tv=[target.get(e,F())for e in mon3];index=next(i for i,v in enumerate(tv)if v)
            Ann=[]
            for i in range(24):
                if i!=index:row=[F()]*24;row[i]=1;row[index]=-tv[i]/tv[index];Ann.append(row)
            assert Ann==mat(parent['full_line_annihilator'],t)
            M=J1+[[sum((a*b for a,b in zip(row,col)),F())for col in zip(*D3)]for row in Ann];assert M==mat(c['full_weight1_plus_line_matrix'],t)
            K1,m1=cert(M,c['weight1_certificate'],t);K2,m2=cert(J2,c['weight2_certificate'],t);Sh=shear(P);assert Sh==mat(c['complete_actual_shear'],t);KS,ms=cert(Sh,c['actual_shear_certificate'],t)
            assert all(not sum((a*b for a,b in zip(row,col)),F())for row in M for col in zip(*Sh))
            image=[[sum((a*b for a,b in zip(row,col)),F())for col in zip(*K1)]for row in D3];assert image==mat(c['whole_D3_kernel_image'],t)
            empty=not any(any(row)for row in image);assert empty==c['whole_m3_source_empty']
            upper=parent['whole_fixed_plane_symbol_dimension_upper']+5-int(t is not None and parent['whole_alpha_line'])+len(K1[0])+len(K2[0]);assert upper==c['whole_projective_image_upper']
            if t is None and roots:
                for d in[kd,kad,m1,m2,ms]:shrink(RF(d),roots);cover_denominators([d],roots)
                cover_denominators([v for k in[K1,K2,KS]for line in k for v in line],roots)
            first_kernel[parent['type'],t]=(K1,P)
            first_results.append(dict(type=parent['type'],alpha=str(t)if t is not None else'whole',empty_by_fullD3=empty,K1=len(K1[0]),K2=len(K2[0]),image_upper=upper,complete_original_rows=True))
        print(parent['type'],'full30/9/24/allwalls PASS',flush=True)
    # Full22 source actual Euler gauges and exact joint inverse, not a chosen Q slice.
    coverage=[]
    for row in ore['full_two_source_coverage']:
        kind=row['type'];K,Pbase=first_kernel[kind,None];assert len(K[0])==13
        params=tuple('p'+str(i)for i in range(13))+tuple('r'+str(i)for i in range(9))+('tau',)
        R=LaurentRing(params,names);p=R.parse;P=L.skew(R,{ij:p(next(c for c in first['cases']if c['type']==kind)['full150_P0'][str(ij)])for ij in q.PAIRS})
        for j in range(22):
            for ij in q.PAIRS:
                raw=q.add(*(q.scale(V.get(ij,{}),K[i][j])for i,V in enumerate(B1)if K[i][j]))if j<13 else B2[j-13].get(ij,{})
                for ex,coef in raw.items():
                    cc=p(rftext(coef));term=R.mul(cc,p(params[j]));mono=R.constant(1)
                    for k,n in enumerate(ex):mono=R.mul(mono,L.power(R,p(names[k]),n))
                    P[ij]=R.add(P[ij],R.mul(term,mono));P[ij[::-1]]=R.scale(P[ij],-1)
        gauge=[R.scale(R.derivative(P[0,2],2),-1),R.scale(R.derivative(P[1,2],2),-1),{},{},{}];assert gauge==[p(v)for v in row['full_source_first_Euler_gauge']]
        xx=list(map(p,names));G=L.skew(R,{(i,j):R.add(P[i,j],R.scale(R.add(R.mul(xx[i],gauge[j]),R.scale(R.mul(xx[j],gauge[i]),-1)),-1))for i,j in q.PAIRS});assert all(not R.derivative(G[ij],2)for ij in q.PAIRS)
        Q=L.skew(R,{ij:G[ij]for ij in itertools.combinations(Y,2)});W=[G[2,i]for i in Y]
        def cf(poly,mon):return L.coef(R,poly,next(iter(p(mon)))[R.n:])
        ell=R.add(R.scale(R.mul(cf(W[2],'z**2'),p('z')),-1),R.scale(R.mul(cf(W[3],'w**2'),p('w')),-1));assert ell==p(row['second_gauge_c_component'])
        W=[R.add(v,R.mul(ell,xx[i]))for v,i in zip(W,Y)];u=cf(Q[0,1],'b*z'if kind.startswith('REPEATED')else'a*z');central=R.add(Q[0,1],R.scale(R.mul(u,p('b*z'if kind.startswith('REPEATED')else'a*z')),-1))
        expected=Qmake(R,u,central,kind.startswith('REPEATED'))
        assert all(Q[ij]==expected[ij]for ij in q.PAIRS)
        Wzero=Wmake(R,kind.startswith('REPEATED'));WB=L.Wbasis(R);coords=[u]+[cf(central,s)for s in('z**2','z*w','w**2')]
        for j in range(18):
            at=next(i for i in range(4)if WB[j][i]);coords.append(L.coef(R,R.add(W[at],R.scale(Wzero[at],-1)),next(iter(WB[j][at]))[R.n:]))
        assert W==[R.add(w,*[R.mul(coords[4+j],WB[j][i])for j in range(18)])for i,w in enumerate(Wzero)]
        J=[]
        for cp in coords:
            assert all(sum(ex[:22])==1 and not any(ex[R.n:])for ex in cp)
            line=[]
            for j in range(22):
                coeff={ex:c for ex,c in cp.items()if ex[j]};expr={ex[22:23]+(0,)*5:c for ex,c in coeff.items()}
                # Joint map is constant, with tau cancellations checked directly.
                assert all(not ex[0]for ex in expr);line.append(sum(expr.values(),F()))
            J.append(line)
        assert J==mat(row['joint22by22_source_coordinate_matrix']);det=determinant(J);assert det==sc(row['joint_nonzero_constant_determinant'])and det
        coverage.append(dict(type=kind,full13plus9=True,joint22_determinant=str(det),all_Euler_and_c_independence=True))
    # Universal strictQ, actual central/shear/diagonal identities and all source inverse charts.
    R=LaurentRing(('tau','u','f0','f1','f2','sa','sb','sc','sd','s','r','k','l'),names,invert=('s','r','k','l'));p=R.parse;xx=list(map(p,names));Fpoly=p('f0*z**2+f1*z*w+f2*w**2')
    for double in[True,False]:
        Q=Qmake(R,p('u'),Fpoly,double);assert all(not a for a in L.strict(R,Q).values())
        P=L.pad(R,Q,Wmake(R,double));D=L.actual(R,P);lead={e:c for e,c in D.items()if e[-2]+e[-1]==3}
        expected=p('3*u*(2/3-2*tau)*(-1/3-2*tau)*z*(b*z+a*w)**2/2'if double else'6*u*z*(b*z+a*w)**2');assert lead==expected
        g=list(map(p,('a+sa*z+sb*w','b+sc*z+sd*w','c','z','w')))
        delta=p('u*((-1/3-2*tau)*sc*z**2+(2/3-2*tau)*(sa+sd)*z*w+(5/3-2*tau)*sb*w**2)'if double else'u*((2*sc-sd-sa)*z**2+(2*sd+2*sa-2*sb)*z*w+2*sb*w**2)')
        L.transform_identity(R,Q,g,Qmake(R,p('u'),R.add(Fpoly,delta),double))
        # Diagonals and tensor character; both the whole Q and W0 are preserved.
        g=list(map(p,('s*a','s*l*b/k'if double else's*b','r*c','k*z','l*w'if double else'k*w')))
        expectedQ=Qmake(R,p('s**2*u/(r*k)'),{},double)
        L.transform_identity(R,Qmake(R,p('u'),{},double),g,expectedQ,p('s/r'))
        L.transform_identity(R,L.pad(R,L.skew(R,{}),Wmake(R,double)),g,L.pad(R,L.skew(R,{}),Wmake(R,double)),p('s/r'))
    assert q.det([[0,0,2],[2,-2,0],[-1,2,0]])==4
    module_results=[]
    for parent in ore['cases']:
        variants=[(parent,None)]+[(w,F(w['alpha']))for w in parent.get('all_exceptional_sources',[])if not w.get('source_empty_by_original_full_D3_image')]
        roots=[F(w['alpha'])for w in parent.get('all_exceptional_sources',[])]
        if roots:shrink(sc(parent['all_exception_polynomial']),roots)
        for c,t in variants:
            Q=q.skew({(Y[i],Y[j]):parse(v,t)for key,v in parent['full_Q'].items()for i,j in[ast.literal_eval(key)]});W=[parse(v,t)for v in parent['W0']]
            def pad(Wv):return q.plus(Q,q.skew({(2,i):v for i,v in zip(Y,Wv)}))
            WB=[]
            for at,terms in[(0,['a*z','b*z','a*w','b*w','z**2','z*w','w**2']),(1,['a*z','b*z','a*w','b*w','z**2','z*w','w**2']),(2,['z*w','w**2']),(3,['z**2','z*w'])]:
                for s in terms:v=[{}for _ in range(4)];v[at]=parse(s,t);WB.append(v)
            P=pad(W);C=[compact(2*v)for v in q.jvector(q.jac(P))];M=q.transpose([[compact(2*(a-b))for a,b in zip(q.jvector(q.jac(pad([q.add(a,b)for a,b in zip(W,V)]))),q.jvector(q.jac(P)))]for V in WB])
            assert C==vec(c['full350_constant'],t)and M==mat(c['full350_matrix'],t)
            K,mi=cert(M,c['linear_certificate'],t);KA,ma=cert([r+[-v]for r,v in zip(M,C)],c['augmented_certificate'],t)
            if c['status']=='inconsistent':assert c['augmented_certificate']['rank']>c['linear_certificate']['rank'];factor='inconsistent'
            else:
                assert c['status']=='complete_affine_module'and c['augmented_certificate']['rank']==c['linear_certificate']['rank']
                part=[sc(v,t)for v in c['particular_solution']];assert all(not(a+sum((b*d for b,d in zip(row,part)),F()))for row,a in zip(M,C))
                D=q.actual(P);Dl=[q.add(q.actual(pad([q.add(a,b)for a,b in zip(W,V)])),q.scale(D,-1))for V in WB]
                Dp=q.add(D,*[q.scale(v,a)for v,a in zip(Dl,part)]);Dk=[q.add(*[q.scale(Dl[i],K[i][j])for i in range(18)])for j in range(len(K[0]))]
                assert Dp==parse(c['full_D_particular'],t)and Dk==[parse(v,t)for v in c['full_D_kernel_coefficients']]
                f=parse(c['proper_D_factor'],t);assert f
                assert 0<max(map(sum,f))<5
                for v in[Dp]+Dk:dividepoly(v,f)
                factor=c['proper_D_factor']
                if t is None and roots:cover_denominators(part+[v for line in K for v in line]+[v for poly in[Dp]+Dk+[f]for v in poly.values()],roots)
            if t is None and roots:
                for d in[mi,ma]:shrink(RF(d),roots);cover_denominators([d],roots)
            module_results.append(dict(type=parent['type'],alpha=str(t)if t is not None else'whole',rank=c['linear_certificate']['rank'],kernel=c['linear_certificate']['nullity'],full350_affine=True,whole_actual_D_factor_or_inconsistent=factor,G_irreducible_source_empty=True))
        print(parent['type'],'all350/18/augmented/fullD/allwalls PASS',flush=True)
    result=dict(task_id='G',all_independent_arithmetic_passed=True,whole_degree_1_0_G_image_upper23=True,
        all11_first_connections_and_walls=first_results,full_two_source_bijections=coverage,all_three_Ore_modules_and_walls=module_results,
        universal_strict_Q_Jacobi_and_corrected_leading=True,all_actual_central_and_diagonal_inverse_identities=True,
        whole_geometry_independently_read=True,all_parameter_minors_denominators_and_zero_complements=True,
        no_producer_or_CAS_imports=True,proof_scope='Whole nonreducedR/m3 generic21 primitive repeated-degree1/simple-degree0 only. Nine first connections empty or image<=23, remaining two entire sources bijectively transfer to three Ore modules all irreducible-D empty. E membership of two weak-bound strata remains open. No global completion or newcomponent.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/outputpath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


def _moving_double_Qmake(R,u,central,double):
    p=R.parse;mul=R.mul;add=R.add
    values=[add(mul(u,p('b*z'if double else'a*z')),central),mul(u,p('(2/3-2*tau)*z**2'if double else'2*z**2')),mul(u,p('(5/3-2*tau)*z*w'if double else'-z**2+2*z*w')),mul(u,p('(2*tau-2/3)*z*w'if double else'-2*z*w')),mul(u,p('(2*tau-5/3)*w**2'if double else'z*w-2*w**2'))]
    import classification_plane_double_line_m3_ore_review as L
    return L.skew(R,{ij:v for ij,v in zip([(0,1),(0,3),(0,4),(1,3),(1,4)],values)})


def _moving_double_Wmake(R,double):
    return list(map(R.parse,('a**2','-2*a*b','3*tau*a*z','-2*a*w+(1-3*tau)*b*z')if double else('0','-3*a**2','-3*a*z','3*(a+b)*z')))


def moving_simple_m3_main():
    """Independent whole eleven-source and three Ore-module reconstruction."""
    import itertools,signal
    from classification_minimal_quartic_review import LaurentRing
    import classification_plane_double_line_m3_ore_review as L
    Qmake=_moving_simple_Qmake;Wmake=_moving_simple_Wmake
    begun=time.monotonic();signal.alarm(360)
    resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jobpath='research/classification/jobs/plane_nonreduced_m3_moving_simple_separate_review.json'
    outputpath='reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_SEPARATE_REVIEW_20261108.json'
    job=json.loads((ROOT/jobpath).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    paths=['reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_20261107.json','reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_ORE_20261108.json']
    paths.append('reports/CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_MODULES_20261108.json')
    first,ore,modules=[json.loads((ROOT/p).read_text())for p in paths]
    for d in[first,ore,modules]:assert d['all_exact_arithmetic_passed']and hashlib.sha256(d['producer_bytes'].encode()).hexdigest()==d['producer_sha256']
    names=('a','b','c','z','w');Y=(0,1,3,4);zero=(0,)*5;X=[q.a(1,i)for i in range(5)]
    @functools.lru_cache(None)
    def parse(text,at=None):
        def walk(n):
            if isinstance(n,ast.Constant):return{zero:F(n.value)}if n.value else{}
            if isinstance(n,ast.Name):
                if n.id=='tau':return{zero:T if at is None else at}
                return X[names.index(n.id)]
            if isinstance(n,ast.UnaryOp):return q.scale(walk(n.operand),-1)if isinstance(n.op,ast.USub)else walk(n.operand)
            a=walk(n.left)
            if isinstance(n.op,ast.Pow):
                out={zero:F(1)}
                for _ in range(n.right.value):out=q.mul(out,a)
                return out
            b=walk(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            assert isinstance(n.op,ast.Div)and set(b)=={zero};return q.scale(a,1/b[zero])
        return walk(ast.parse(str(text),mode='eval').body)
    def sc(v,at=None):
        if at is None:return scalar(str(v))
        u=RF(scalar(str(v)));ev=lambda pp:sum((c*at**i for i,c in enumerate(pp)),F());return ev(u.p)/ev(u.q)
    def mat(v,at=None):return[[sc(t,at)for t in row]for row in v]
    def vec(v,at=None):return[sc(t[0]if isinstance(t,list)else t,at)for t in v]
    def savedP(v,at=None):return q.skew({ast.literal_eval(ij):parse(text,at)for ij,text in v.items()})
    def op(P,V):
        base=q.jvector(q.jac(P));selfq=q.jvector(q.jac(V));t=q.jvector(q.jac(q.plus(P,V)))
        return[compact(a-b-c)for a,b,c in zip(t,base,selfq)]
    def cert(M,C,at=None):
        rr,cc=C['minor_rows'],C['minor_columns'];r=C['rank'];n=len(M[0]);assert len(rr)==len(cc)==r and len(set(rr))==len(set(cc))==r
        d=determinant([[M[i][j]for j in cc]for i in rr]);assert d and d==sc(C['minor_determinant'],at)
        K=mat(C['complete_kernel'],at);free=[j for j in range(n)if j not in cc];assert len(free)==C['nullity']and len(K)==n
        assert[K[i]for i in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
        for row in M:
            nz=[(i,a)for i,a in enumerate(row)if a]
            for k in range(len(free)):assert not sum((a*K[i][k]for i,a in nz if K[i][k]),F())
        return K,d
    def rftext(c):
        if not isinstance(c,RF):return str(c)
        poly=lambda vv:'+'.join('('+str(a)+')*tau**'+str(i)for i,a in enumerate(vv)if a)or'0'
        return '('+poly(c.p)+')/('+poly(c.q)+')'
    def evalP(P,t):
        return{e:sc(rftext(c),t)for e,c in P.items()if sc(rftext(c),t)}
    def shrink(p,roots):
        a=list(RF(p).p)
        for r in roots:
            while len(a)>1 and not sum((c*r**i for i,c in enumerate(a)),F()):
                b=[F()]* (len(a)-1);b[-1]=a[-1]
                for i in range(len(a)-2,0,-1):b[i-1]=a[i]+r*b[i]
                assert a[0]+r*b[0]==0;a=b
        assert len(a)<=1 and a,('uncovered polynomial roots',a,roots)
    def cover_denominators(objects,roots):
        for v in objects:
            u=RF(v);shrink(RF(u.q),roots)
    def dividepoly(P,D):
        assert D;R=P.copy();lead=max(D);coef=D[lead]
        while R:
            e=max(R);assert all(a>=b for a,b in zip(e,lead)),('factor mismatch',P,D)
            exp=tuple(a-b for a,b in zip(e,lead));c=R[e]/coef
            R=q.add(R,q.scale(q.mul({exp:c},D),-1))
        return True
    B1,B2=q.sources();B1=B1[:18]+B1[24:27]+B1[18:21]+B1[27:30]+B1[21:24]
    base2=[e[:3]for e in q.monomials(2)if not any(e[3:])];base3=[e[:3]for e in q.monomials(3)if not any(e[3:])]
    mon3=[e+(j,3-j)for e in base2 for j in range(4)];mon2=[e+(j,2-j)for e in base3 for j in range(3)]
    assert list(map(list,mon3))==first['all24_D3_monomials']and list(map(list,mon2))==first['all30_K_monomials']
    R=LaurentRing(('tau',),names,invert=('a','b'));rp=R.parse
    def reconstruct_forms(c,at=None):
        fs='a**3'if 'TRIPLE'in c['type']else'a**2*b';f=parse(fs,at)
        P={(0,1):q.der(f,2),(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)}
        for form in c['full_connection_inverse_forms']:
            i,j=form['i'],form['j'];beta=[rp(t)for t in form['beta']];om=[rp(t)for t in form['horizontal_polynomial_form']]
            assert om==[R.add(R.scale(R.mul(rp(fs),v),3),R.scale(R.derivative(rp(fs),k),-int(i==j)))for k,v in enumerate(beta)]
            assert not R.add(*(R.mul(rp(names[k]),om[k])for k in range(3)))
            om=[parse(t,at)for t in form['horizontal_polynomial_form']]
            V=[q.scale(q.add(q.der(om[2],1),q.scale(q.der(om[1],2),-1)),F(1,3)),q.scale(q.add(q.der(om[0],2),q.scale(q.der(om[2],0),-1)),F(1,3)),q.scale(q.add(q.der(om[1],0),q.scale(q.der(om[0],1),-1)),F(1,3))]
            assert V==[parse(t,at)for t in form['inverse_tracefree_field']]
            inv=[q.add(q.mul(X[2],V[1]),q.scale(q.mul(X[1],V[2]),-1)),q.add(q.mul(X[0],V[2]),q.scale(q.mul(X[2],V[0]),-1)),q.add(q.mul(X[1],V[0]),q.scale(q.mul(X[0],V[1]),-1))]
            assert inv==om
            for k in range(3):P[k,3+i]=q.add(P.get((k,3+i),{}),q.scale(q.mul(V[k],X[3+j]),-1))
        return q.skew(P)
    def shear(P):
        cols=[]
        for normal in(3,4):
            for k in range(3):
                G=[{}for _ in range(5)];G[k]=X[normal]
                O=q.skew({(i,j):q.add(*(q.scale(q.mul(G[l],q.der(P.get((i,j),{}),l)),-1)for l in range(5)),*(q.add(q.mul(P.get((l,j),{}),q.der(G[i],l)),q.mul(P.get((i,l),{}),q.der(G[j],l)))for l in range(5)))for i,j in q.PAIRS})
                S=[q.add(q.scale(X[3],O.get((i,3),{}).get((0,0,0,2,0),F())),q.scale(X[4],O.get((i,4),{}).get((0,0,0,0,2),F())))for i in range(3)]+[{},{}]
                O=q.plus(O,q.wedge(X,S));raw=q.vector(O)
                out=[]
                for V in B1:
                    loc=next((ij,e)for ij,p in V.items()if ij in q.PAIRS for e in p);out.append(O.get(loc[0],{}).get(loc[1],F()))
                assert raw==q.vector(q.skew({ij:q.add(*(q.scale(V.get(ij,{}),a)for V,a in zip(B1,out)if a))for ij in q.PAIRS}));cols.append(out)
        return q.transpose(cols)
    first_results=[];first_kernel={};whole_minor_polys=[]
    for parent in first['cases']:
        genericP=reconstruct_forms(parent)
        variants=[(parent,None)]+[(dict(w,type=parent['type']),F(w['alpha']))for w in parent.get('all_exceptional_sources',[])]
        roots=[F(w['alpha'])for w in parent.get('all_exceptional_sources',[])]
        if roots:shrink(sc(parent['all_exception_polynomial']),roots)
        for c,t in variants:
            Pbase=genericP if t is None else q.skew({ij:evalP(p,t)for ij,p in genericP.items()})
            D=q.actual(Pbase);KM=[]
            for nm in[(3,3),(3,4),(4,4)]:
                DD=q.actual(q.plus(Pbase,q.skew({(3,4):q.a(1,*nm)})))
                KM.append([compact(DD.get(e,F())-D.get(e,F()))for e in mon2])
            KM=q.transpose(KM);constant=[D.get(e,F())for e in mon2]
            assert KM==mat(c['full30_global_K_matrix'],t)and constant==vec(c['full30_global_K_constant'],t)
            _,kd=cert(KM,c['K_certificate'],t);_,kad=cert([r+[-a]for r,a in zip(KM,constant)],c['K_augmented_certificate'],t)
            if c.get('global_K_inconsistent'):
                assert c['K_augmented_certificate']['rank']>c['K_certificate']['rank'];first_results.append(dict(type=parent['type'],alpha=str(t),global_K_inconsistent=True));continue
            P=savedP(c['full150_P0'],t);Kp=parse(c['unique_global_K'],t)
            assert all(P.get(ij,{})==Pbase.get(ij,{})for ij in q.PAIRS if ij!=(3,4));assert P[3,4]==Kp
            assert not q.actual(P)and not any(q.jac(P).values())
            J1=q.transpose([op(P,V)for V in B1]);J2=q.transpose([op(P,V)for V in B2])
            assert J1==mat(c['full350_weight1_operator'],t)and J2==mat(c['full350_weight2_operator'],t)
            D3=q.transpose([[q.actual(q.plus(P,V)).get(e,F())for e in mon3]for V in B1]);assert D3==mat(c['full24by30_actual_D3_map'],t)
            target=parse(parent['exact_nonzero_leading_target'],t);tv=[target.get(e,F())for e in mon3];index=next(i for i,v in enumerate(tv)if v)
            Ann=[]
            for i in range(24):
                if i!=index:row=[F()]*24;row[i]=1;row[index]=-tv[i]/tv[index];Ann.append(row)
            assert Ann==mat(parent['full_line_annihilator'],t)
            M=J1+[[sum((a*b for a,b in zip(row,col)),F())for col in zip(*D3)]for row in Ann];assert M==mat(c['full_weight1_plus_line_matrix'],t)
            K1,m1=cert(M,c['weight1_certificate'],t);K2,m2=cert(J2,c['weight2_certificate'],t);Sh=shear(P);assert Sh==mat(c['complete_actual_shear'],t);KS,ms=cert(Sh,c['actual_shear_certificate'],t)
            assert all(not sum((a*b for a,b in zip(row,col)),F())for row in M for col in zip(*Sh))
            image=[[sum((a*b for a,b in zip(row,col)),F())for col in zip(*K1)]for row in D3];assert image==mat(c['whole_D3_kernel_image'],t)
            om=next(v for v in parent['full_connection_inverse_forms']if v['i']==v['j']==1)
            f=parse('a**3'if'TRIPLE'in parent['type']else'a**2*b',t)
            qform=[q.add(parse(v,t),q.der(f,i))for i,v in enumerate(om['horizontal_polynomial_form'])]
            QM=[[F()]*30 for _ in range(6)]
            for i,ex in enumerate([e for e in q.monomials(2)if not any(e[3:])]):
                for j in range(3):QM[i][-3+j]=qform[j].get(ex,F())
            assert QM==mat(c['full6by30_constant_root_character_matrix'],t)
            assert all(not sum((a*b for a,b in zip(row,col)),F())for row in QM for col in zip(*K1))
            empty=not any(any(row)for row in image);assert empty==c['whole_m3_source_empty']
            upper=parent['whole_fixed_plane_symbol_dimension_upper']+5-int(t is not None and parent['whole_alpha_line'])+len(K1[0])+len(K2[0]);assert upper==c['whole_projective_image_upper']
            if t is None and roots:
                for d in[kd,kad,m1,m2,ms]:shrink(RF(d),roots);cover_denominators([d],roots)
                cover_denominators([v for k in[K1,K2,KS]for line in k for v in line],roots)
            first_kernel[parent['type'],t]=(K1,P)
            first_results.append(dict(type=parent['type'],alpha=str(t)if t is not None else'whole',empty_by_fullD3=empty,K1=len(K1[0]),K2=len(K2[0]),image_upper=upper,complete_original_rows=True))
        print(parent['type'],'full30/9/24/allwalls PASS',flush=True)
    # Three full bad-column exclusions, independently from the complete original kernels.
    bad_results=[]
    for ex in modules['three_full_bad_column_zero_sources']:
        t=F(ex['alpha'])if ex['alpha']is not None else None;K,P=first_kernel[ex['type'],t]
        assert K==mat(ex['complete_P1_kernel'])and all(not v for line in K[-3:]for v in line)
        assert all(e[3]for ij,p in P.items()if 3 in ij for e in p)
        assert all(not V.get((i,3),{}).get((0,0,0,0,2),F())for V in B2 for i in range(3))
        bad_results.append(dict(type=ex['type'],alpha=ex['alpha'],entire_z_row_divisible=True,whole_actualD_z_factor=True))
    class FieldRing(LaurentRing):
        def constant(self,c):return{self.zero:RF(c)}if c else{}
        def parse(self,text):
            def walk(n):
                if isinstance(n,ast.Constant):return self.constant(F(n.value))
                if isinstance(n,ast.Name):return self.constant(T)if n.id=='tau'else self.variable(n.id)
                if isinstance(n,ast.UnaryOp):return self.scale(walk(n.operand),-1)if isinstance(n.op,ast.USub)else walk(n.operand)
                a=walk(n.left)
                if isinstance(n.op,ast.Pow):
                    out=self.constant(1)
                    for _ in range(n.right.value):out=self.mul(out,a)
                    return out
                b=walk(n.right)
                if isinstance(n.op,ast.Add):return self.add(a,b)
                if isinstance(n.op,ast.Sub):return self.add(a,self.scale(b,-1))
                if isinstance(n.op,ast.Mult):return self.mul(a,b)
                assert isinstance(n.op,ast.Div);return self.mul(a,self.inverse_monomial(b))
            return walk(ast.parse(str(text),mode='eval').body)
    def get_original(label,at):
        if label=='DEG1_HA_CA_OTHER_WALL':return first_kernel['DEG1_HA_CA_OTHER',F(-8,3)]
        if label=='DEG1_HB_CA_OTHER_WALL':return first_kernel['DEG1_HB_CA_OTHER',F(4,3)]
        return first_kernel[label,at]
    def source_gauge(label,at,vertex,zero_gauge,claimed_gauge):
        K,Pbase=get_original(label,at);n=len(K[0]);params=tuple('p'+str(i)for i in range(n))+tuple('r'+str(i)for i in range(9))
        R=FieldRing(params,names);p=R.parse;xx=list(map(p,names));P=L.skew(R,{ij:p(polytext(Pbase.get(ij,{})))for ij in q.PAIRS})
        Z=list(map(p,zero_gauge));P=L.skew(R,{(i,j):R.add(P[i,j],R.scale(R.add(R.mul(xx[i],Z[j]),R.scale(R.mul(xx[j],Z[i]),-1)),-1))for i,j in q.PAIRS})
        for j in range(n+9):
            for ij in q.PAIRS:
                raw=q.add(*(q.scale(V.get(ij,{}),K[i][j])for i,V in enumerate(B1)if K[i][j]))if j<n else B2[j-n].get(ij,{})
                pp=p(polytext(raw));pp=R.mul(pp,p(params[j]));P[ij]=R.add(P[ij],pp);P[ij[::-1]]=R.scale(P[ij],-1)
        gauge=[R.scale(R.derivative(P.get((i,vertex),{}),vertex),-1)if i!=vertex else{}for i in range(5)]
        assert gauge==[p(t)for t in claimed_gauge]
        G=L.skew(R,{(i,j):R.add(P[i,j],R.scale(R.add(R.mul(xx[i],gauge[j]),R.scale(R.mul(xx[j],gauge[i]),-1)),-1))for i,j in q.PAIRS})
        return R,G,params
    def polytext(P):
        return '+'.join('('+rftext(c)+')'+''.join('*'+names[i]+'**'+str(k)for i,k in enumerate(e)if k)for e,c in P.items())or'0'
    def coefficient(R,P,e):return L.coef(R,P,e)
    def constant_part(R,P):return{e:c for e,c in P.items()if not any(e[:R.n])}
    def source_coords(R,polys,params):
        out=[]
        for pp in polys:
            assert all(sum(e[:R.n])==1 and not any(e[R.n:])for e in pp)
            out.append([pp.get(tuple(int(i==j)for i in range(R.n))+(0,)*5,F())for j in range(len(params))])
        return out
    # Retain the old seven-source diagnosis, including its genuine b-dependent failure.
    ore_diagnostics=[]
    for parent in ore['cases']:
        variants=[(parent,None)]+[(w,F(w['alpha']))for w in parent.get('all_exceptional_sources',[])if not w.get('whole_m3_source_empty_by_original_certificate')]
        for c,at in variants:
            R,G,params=source_gauge(parent['type'],at,c['vertex_index'],c['weight_zero_Euler_gauge'],c['full_positive_weight_Euler_gauge']);p=R.parse;vertex=c['vertex_index'];indices=tuple(i for i in range(5)if i!=vertex)
            bad=any(R.derivative(G[ij],vertex)for ij in q.PAIRS);assert bad==(not c['source_covered'])
            if bad:ore_diagnostics.append(dict(type=parent['type'],alpha=str(at),genuine_vertex_dependence=True,source_covered=False));continue
            Q=L.skew(R,{ij:G[ij]for ij in itertools.combinations(indices,2)});W=[G[vertex,i]for i in indices]
            mono=lambda text:next(iter(p(text)))[R.n:]
            cf=lambda pp,text:coefficient(R,pp,mono(text))
            ell=R.add(R.scale(R.mul(cf(W[2],'z**2'),p('z')),-1),R.scale(R.mul(cf(W[3],'w**2'),p('w')),-1));assert ell==p(c['remaining_Euler_c_component'])
            W=[R.add(v,R.mul(ell,p(names[i])))for v,i in zip(W,indices)];Q0={ij:constant_part(R,Q[ij])for ij in Q}
            assert all(Q0[indices[i],indices[j]]==p(t)for key,t in c['full_Q0'].items()for i,j in[ast.literal_eval(key)])
            W0=[constant_part(R,v)for v in W];assert W0==[p(t)for t in c['full_W0']]
            localmons=[tuple(sum(k==i for k in pair)for i in indices)for pair in itertools.combinations_with_replacement(indices,2)]
            mons=[tuple(e[indices.index(i)]if i in indices else 0 for i in range(5))for e in localmons]
            qcoords=[coefficient(R,R.add(Q[ij],R.scale(Q0[ij],-1)),e)for ij in itertools.combinations(indices,2)for e in mons]
            J=source_coords(R,qcoords,params);assert J==mat(c['full60_Q_source_matrix'],at)
            qc=dict(c['Q_space_certificate']);K,md=cert(J,qc,at)
            wcoords=[]
            for i in range(2):wcoords.extend(cf(R.add(W[i],R.scale(W0[i],-1)),term)for term in[names[indices[0]]+'*z',names[indices[1]]+'*z',names[indices[0]]+'*w',names[indices[1]]+'*w','z**2','z*w','w**2'])
            wcoords.extend([cf(R.add(W[2],R.scale(W0[2],-1)),'z*w'),cf(R.add(W[2],R.scale(W0[2],-1)),'w**2'),cf(R.add(W[3],R.scale(W0[3],-1)),'z**2'),cf(R.add(W[3],R.scale(W0[3],-1)),'z*w')])
            assert source_coords(R,wcoords,params)==mat(c['full_W18_source_matrix'],at)
            ore_diagnostics.append(dict(type=parent['type'],alpha=str(at),full_original_source_vertex_gauge=True,Q_rank=qc['rank'],W18_map_reconstructed=True,source_covered=True))
    # Four complete joint22 source bijections to the announced canonical Q spaces.
    coverage=[]
    for row in modules['full_four_source_coverage']:
        label=row['type'];kind='HA'if'_HA_'in label else'HB'if'_HB_'in label else'TR'if'TRIPLE'in label else'CB'
        R,G,params=source_gauge(label,None,2,['0']*5,row['full_Euler_gauge']);p=R.parse;assert len(params)==22 and all(not R.derivative(G[ij],2)for ij in q.PAIRS)
        Q=L.skew(R,{ij:G[ij]for ij in itertools.combinations(Y,2)});W=[G[2,i]for i in Y];cf=lambda pp,text:coefficient(R,pp,next(iter(p(text)))[R.n:])
        ell=R.add(R.scale(R.mul(cf(W[2],'z**2'),p('z')),-1),R.scale(R.mul(cf(W[3],'w**2'),p('w')),-1));W=[R.add(w,R.mul(ell,p(names[i])))for w,i in zip(W,Y)]
        u=cf(Q[1,3],'z**2')if kind in['HA','TR']else cf(Q[0,4],'z*w')if kind=='HB'else cf(Q[0,3],'z**2')
        central=Q[0,1]
        if kind=='HB':central=R.add(central,R.scale(R.mul(u,p('(3*tau+2)*a*w/(3*tau-4)+b*z')),-1))
        expected=Qmake(R,u,central,kind);assert all(Q[ij]==expected[ij]for ij in q.PAIRS)
        W0=Wmake(R,kind);WB=L.Wbasis(R);coords=[u]+[cf(central,t)for t in('z**2','z*w','w**2')]
        for j in range(18):i=next(i for i in range(4)if WB[j][i]);coords.append(coefficient(R,R.add(W[i],R.scale(W0[i],-1)),next(iter(WB[j][i]))[R.n:]))
        assert W==[R.add(w,*[R.mul(coords[4+j],WB[j][i])for j in range(18)])for i,w in enumerate(W0)]
        J=source_coords(R,coords,params);assert J==mat(row['full_joint22_coordinate_matrix']);det=determinant(J);assert det==sc(row['nonzero_constant_joint_determinant'])and det==F(-1)
        coverage.append(dict(type=label,joint22_determinant=str(det),full13plus_all9_source=True,all_c_independence_and_QW_identity=True))
    # Actual transformations; rational tau coefficients stay in QQ(tau).
    R=FieldRing(('u','f0','f1','f2','sa','sb','sc','sd','s','r','k','l'),names,invert=('s','r','k','l'));p=R.parse;xx=list(map(p,names));central=p('f0*z**2+f1*z*w+f2*w**2');geometries=[]
    for kind in['HA','HB','TR','CB']:
        Q=Qmake(R,p('u'),central,kind);assert all(not v for v in L.strict(R,Q).values())
        P=L.pad(R,Q,Wmake(R,kind));D=L.actual(R,P);lead={e:c for e,c in D.items()if e[-2]+e[-1]==3}
        row=next(r for r in modules['full_four_source_coverage']if('HA'if'_HA_'in r['type']else'HB'if'_HB_'in r['type']else'TR'if'TRIPLE'in r['type']else'CB')==kind)
        original=next(c for c in first['cases']if c['type']==row['type']);amp=p(row['exact_D3_target_amplitude']);assert lead==R.mul(amp,p(original['exact_nonzero_leading_target']))
        g=list(map(p,('a+sa*z+sb*w','b+sc*z+sd*w','c','z','w')))
        delta=p('u*((-(3*tau+5)*sd/(3*tau+8)-sa)*z**2-2*sb*z*w)'if kind=='HA'else'u*(-sc*z**2-((1-3*tau)/(3*tau-4))*sb*w**2)'if kind=='HB'else'u*(-(sd+sa+sb)*z**2-2*sb*z*w)'if kind=='TR'else'u*((sc+sb)*z**2+2*sd*z*w)')
        L.transform_identity(R,Q,g,Qmake(R,p('u'),R.add(central,delta),kind))
        g=list(map(p,('s*a','s*l*b/k'if kind in['HA','HB']else's*b','r*c','k*z','l*w'if kind in['HA','HB']else'k*w')))
        up=p('s**2*l*u/(r*k**2)'if kind=='HA'else's**2*u/(r*k)')
        L.transform_identity(R,Qmake(R,p('u'),{},kind),g,Qmake(R,up,{},kind),p('s/r'))
        L.transform_identity(R,L.pad(R,L.skew(R,{}),Wmake(R,kind)),g,L.pad(R,L.skew(R,{}),Wmake(R,kind)),p('s/r'))
        geometries.append(dict(type=kind,strict_Q_Jacobi=True,full_leading_amplitude=True,actual_central_and_diagonal_identities=True))
    module_results=[]
    for parent in modules['cases']:
        variants=[(parent,None)]+[(w,F(w['alpha']))for w in parent.get('all_exceptional_sources',[])if not w.get('source_empty_by_original_full_D3_image')]
        roots=[F(w['alpha'])for w in parent.get('all_exceptional_sources',[])]
        if roots:shrink(sc(parent['all_exception_polynomial']),roots)
        for c,t in variants:
            Q=q.skew({(Y[i],Y[j]):parse(v,t)for key,v in parent['full_Q'].items()for i,j in[ast.literal_eval(key)]});W=[parse(v,t)for v in parent['W0']]
            def pad(Wv):return q.plus(Q,q.skew({(2,i):v for i,v in zip(Y,Wv)}))
            WB=[]
            for at,terms in[(0,['a*z','b*z','a*w','b*w','z**2','z*w','w**2']),(1,['a*z','b*z','a*w','b*w','z**2','z*w','w**2']),(2,['z*w','w**2']),(3,['z**2','z*w'])]:
                for s in terms:v=[{}for _ in range(4)];v[at]=parse(s,t);WB.append(v)
            P=pad(W);C=[compact(2*v)for v in q.jvector(q.jac(P))];M=q.transpose([[compact(2*(a-b))for a,b in zip(q.jvector(q.jac(pad([q.add(a,b)for a,b in zip(W,V)]))),q.jvector(q.jac(P)))]for V in WB])
            assert C==vec(c['full350_constant'],t)and M==mat(c['full350_matrix'],t)
            K,mi=cert(M,c['linear_certificate'],t);KA,ma=cert([r+[-v]for r,v in zip(M,C)],c['augmented_certificate'],t)
            if c['status']=='inconsistent':assert c['augmented_certificate']['rank']>c['linear_certificate']['rank'];factor='inconsistent'
            else:
                assert c['status']=='complete_affine_module'and c['augmented_certificate']['rank']==c['linear_certificate']['rank']
                part=[sc(v,t)for v in c['particular_solution']];assert all(not(a+sum((b*d for b,d in zip(row,part)),F()))for row,a in zip(M,C))
                D=q.actual(P);Dl=[q.add(q.actual(pad([q.add(a,b)for a,b in zip(W,V)])),q.scale(D,-1))for V in WB]
                Dp=q.add(D,*[q.scale(v,a)for v,a in zip(Dl,part)]);Dk=[q.add(*[q.scale(Dl[i],K[i][j])for i in range(18)])for j in range(len(K[0]))]
                assert Dp==parse(c['full_D_particular'],t)and Dk==[parse(v,t)for v in c['full_D_kernel_coefficients']]
                if c['whole_D_zero']:
                    assert not Dp and not any(Dk);f={};factor='zeroD'
                elif c['proper_D_factor']is not None:
                    f=parse(c['proper_D_factor'],t);assert f and 0<max(map(sum,f))<5
                    for v in[Dp]+Dk:dividepoly(v,f)
                    factor=c['proper_D_factor']
                else:
                    f={};factor='weak_image_bound_'+str(c['whole_projective_image_upper'])
                    assert c['whole_projective_image_upper']==c['image_offset']+len(K[0])<=24
                if t is None and roots:cover_denominators(part+[v for line in K for v in line]+[v for poly in[Dp]+Dk+[f]for v in poly.values()],roots)
            if t is None and roots:
                for d in[mi,ma]:shrink(RF(d),roots);cover_denominators([d],roots)
            expected=dict(HA_0=17,HB_0=17,HB_EXTRA_0=16,TRIPLE_0=16,CB2_0=17,HA_1=18,HB_1=18,HB_EXTRA_1=17,TRIPLE_1=17,CB2_1=18)[parent['type']]-int(t is not None)
            assert c['image_offset']==expected
            module_results.append(dict(type=parent['type'],alpha=str(t)if t is not None else'whole',rank=c['linear_certificate']['rank'],kernel=c['linear_certificate']['nullity'],full350_affine=True,whole_actual_D_factor_or_inconsistent=factor,G_irreducible_source_empty=c.get('G_empty',True),whole_image_upper=c.get('whole_projective_image_upper')))
        print(parent['type'],'all350/18/augmented/fullD/allwalls PASS',flush=True)
    result=dict(task_id='G',all_independent_arithmetic_passed=True,whole_degrees_0_1_and_0_2_G_image_upper24=True,
        all37_first_connections_and_walls=first_results,three_whole_bad_column_exclusions=bad_results,all7_vertex_diagnostics=ore_diagnostics,full_four_source_bijections=coverage,universal_actual_transformation_identities=geometries,all_ten_modules_and_walls=module_results,
        four_EXACT_C_sources_empty_by_full_constant_root_contraction=True,all_global_normalizer_inverse_and_primitive_quadratic_case_splits_read=True,
        producer_summary_correction='HB_0 at tau0 and HA_1 at2/15 are consistent modules without a commonfactor; each whole image<=23 by source bound. They are not declared empty. CA2 vertex failure remains a real diagnosis but fullbadcolumn excludes its irreducible source.',
        no_producer_or_CAS_imports=True,scope='Whole nonredR/m3 generic21 moving-simple primitive degrees(0,1)/(0,2) G<=24. Weak-bound strata E unresolved. Otherm/roots/reducedR outside.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/outputpath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


def _moving_simple_Qmake(R,u,central,kind):
    p=R.parse;mul=R.mul;add=R.add
    if kind=='HA':values=[central,{},mul(u,p('-(3*tau+5)*z**2/(3*tau+8)')),mul(u,p('z**2')),mul(u,p('2*z*w'))]
    elif kind=='HB':values=[add(central,mul(u,p('(3*tau+2)*a*w/(3*tau-4)+b*z'))),{},mul(u,p('z*w')),mul(u,p('-(3*tau+2)*z*w/(3*tau-4)')),mul(u,p('-(6*tau+1)*w**2/(3*tau-4)'))]
    elif kind=='TR':values=[central,{},mul(u,p('-z**2')),mul(u,p('z**2')),mul(u,p('2*z*w+z**2'))]
    else:assert kind=='CB';values=[central,mul(u,p('z**2')),mul(u,p('2*z*w')),{},mul(u,p('-z**2'))]
    import classification_plane_double_line_m3_ore_review as L
    return L.skew(R,{ij:v for ij,v in zip([(0,1),(0,3),(0,4),(1,3),(1,4)],values)})


def _moving_simple_Wmake(R,kind):
    return list(map(R.parse,dict(HA=('a**2','-2*a*b','-2*a*z','3*tau*a*w+(3*tau+5)*b*z'),HB=('a**2','-2*a*b','a*z','3*tau*a*w+(3*tau+2)*b*z'),TR=('0','-3*a**2','0','-3*a*w+3*(a-b)*z'),CB=('a**2','-2*a*b','4*a*z','10*a*w-3*a*z-3*b*w'))[kind]))


def triple_line_sf01_main():
    """Independent full350/24/150 replay of all nine reduced-m3 whole sources."""
    from itertools import product
    begun=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jobpath='research/classification/jobs/plane_triple_line_m2_sf01_five_separate_review.json'
    inputpath='reports/CLASSIFICATION_PLANE_TRIPLE_LINE_M2_SF01_FIVE_FULLHIGHER_20261108.json'
    outputpath='reports/CLASSIFICATION_PLANE_TRIPLE_LINE_M2_SF01_FIVE_SEPARATE_REVIEW_20261108.json'
    assert not(ROOT/outputpath).exists();job=json.loads((ROOT/jobpath).read_text())
    for row in job['exact_inputs_and_hashes']:
        assert hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
    data=json.loads((ROOT/inputpath).read_text());assert data['all_passed']and len(data['cases'])==5
    assert hashlib.sha256(data['producer_source'].encode()).hexdigest()==data['producer_sha256']
    def parse(text,names):
        def walk(n):
            if isinstance(n,ast.Constant):return q.a(F(n.value))
            if isinstance(n,ast.Name):return q.a(1,names.index(n.id))
            if isinstance(n,ast.UnaryOp):return q.scale(walk(n.operand),-1 if isinstance(n.op,ast.USub)else 1)
            a=walk(n.left)
            if isinstance(n.op,ast.Pow):
                out=q.a(1)
                for _ in range(n.right.value):out=q.mul(out,a)
                return out
            b=walk(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            if isinstance(n.op,ast.Div):
                zero=(0,)*5;assert set(b)=={zero};return q.scale(a,1/b[zero])
            raise ValueError(ast.dump(n))
        return walk(ast.parse(str(text),mode='eval').body)
    def exact(M):return[[F(v)for v in row]for row in M]
    def verify(M,C,n):
        r=C['rank'];rr,cc=C['rows'],C['columns'];assert len(M[0])==n and len(rr)==len(cc)==r
        assert len(set(rr))==len(set(cc))==r
        minor=q.det([[M[i][j]for j in cc]for i in rr]);assert minor==F(C['minor_determinant'])and minor
        K=exact(C['complete_kernel']);assert len(K)==n-r==C['nullity']
        free=[j for j in range(n)if j not in cc]
        assert [[row[j]for j in free]for row in K]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
        for row in M:
            nz=[(i,v)for i,v in enumerate(row)if v]
            assert all(not sum((c*v[i]for i,c in nz),F())for v in K)
        return K,str(minor)
    first,second=q.sources();assert len(first)==30 and len(second)==9
    mon3=[e for e in q.monomials(5)if e[3]+e[4]==3];assert len(mon3)==24
    X=[q.a(1,i)for i in range(5)]
    eulers=[q.vector(q.skew({(i,j):q.add(q.a(1,i,b)if j==a else{},q.a(-1,j,b)if i==a else{})for i,j in q.PAIRS}))for a,b in product(range(5),repeat=2)]
    assert len(q.rref(q.transpose(eulers))[1])==24
    levi=[(a,b)for a,b in product(range(5),repeat=2)if(a<3)==(b<3)]
    inputs={}
    for item in data['source_input_files']:
        name=item['path'];assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==item['sha256']
        packet=json.loads((ROOT/name).read_text());inputs[name]=packet
    records=[]
    for c in data['cases']:
        packet=inputs[c['input']];names=packet['variables'];profile=next(v for v in packet['profiles']if v['name']==c['name'])
        P=q.skew({tuple(names.index(s)for s in key.split(',')):parse(text,names)for key,text in profile['raw_brackets'].items()})
        assert q.vector(P)==list(map(F,c['full150_P0']))and not any(q.jac(P).values())
        assert q.actual(P)==parse(profile['G2'],names)==parse(c['actual_G2'],names)and q.actual(P)
        J1=q.transpose([q.jvector(q.jac(q.plus(P,V)))for V in first]);J2=q.transpose([q.jvector(q.jac(q.plus(P,V)))for V in second])
        assert J1==exact(c['complete350x30'])and J2==exact(c['complete350x9'])
        K0,minor0=verify(J1,c['weight1_unconstrained_certificate'],30);K2,minor2=verify(J2,c['weight2_certificate'],9)
        K1,minor1=verify(J1,c['constrained_weight1_certificate'],30)
        assert K0==K1
        orbit=[]
        for a,b in levi:
            O=q.skew({(i,j):q.add(P.get((b,j),{})if i==a else{},P.get((i,b),{})if j==a else{},q.scale(q.mul(X[b],q.der(P.get((i,j),{}),a)),-1))for i,j in q.PAIRS})
            orbit.append(q.vector(O))
        G=q.transpose(eulers+orbit+[q.vector(P)]);assert G==exact(c['full150x39_Euler_Levi_scale'])
        KG,mg=verify(G,c['filtered_group_certificate'],39)
        B=c['filtered_group_certificate']['rank']-24;assert B==c['first_symbol_dimension']==10
        badmon=(0,0,0,0,2)
        ufactor=all(e[3]for ij,p in P.items()if 3 in ij for e in p)
        for row in K1:
            Q=q.skew({ij:q.add(*(q.scale(V.get(ij,{}),a)for V,a in zip(first,row)if a))for ij in q.PAIRS})
            ufactor=ufactor and all(not Q.get((i,3),{}).get(badmon,F())for i in range(3))
        assert ufactor==c['entire_bad_u_row_zero']==c['actual_D_has_proper_u_factor']
        higher=len(K1)+len(K2);upper=B+higher+5
        assert higher==c['full_higher_dimension_upper']and upper==c['whole_projective_image_upper']
        assert ufactor and all(not p for V in second for ij,p in V.items()if 3 in ij);assert c['G24_closed']
        records.append(dict(name=c['name'],full150_P0_and350_Jacobi_and_nonzero_G2=True,
            full350x30_and350x9_reconstructed=True,
            unconstrained_K1=len(K0),constrained_K1=len(K1),K2=len(K2),J1_minor=minor1,J2_minor=minor2,
            full150x39_filtered_group_reconstructed=True,filtered_group_rank=c['filtered_group_certificate']['rank'],filtered_group_minor=mg,
            complete_original_kernels_checked=True,first_symbol_B=B,higher_necessary_dimension_upper=higher,projective_image_upper=upper,
            proper_actual_u_factor_for_entire_source=ufactor,G24_closed=True))
        print(c['name'],'independent350/150/kernel/fullu-factor PASS',flush=True)
    assert [c['name']for c in records]==['SF01_p_-3','SF01_p_0','SF01_p_3','SF01_p_-3_2','SF01_p_-1']
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all_five_SF01_whole_irreducible_sources_empty=True,cases=records,
        mathematical_bridge='All fullP1 belongs to complete originalkerJ1. Its wholebad-u-row is zero and originalP0 u-row divisible byu. AllnineP2 has no normalcolumn, so inhomogeneous P2particular and every fullsolution preserve u-row divisor. EulerPf gives properufactor, contradicting irreduciblequintic. Group/bounds are diagnostic only.',
        scope='Whole arithmetic and fullu-factor exclusion for allfive R3L/m2 SF01 fibers. Exhaustive inverse geometry is independently accepted input205a. SF00nil/otherR/orders remain open. No producer import or CAS.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),
        immutable_producer_snapshot=inputpath+':producer_source',seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/outputpath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,all5_SF01_empty=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))

def concurrent_nonore_m3_main():
    """Independent 56 whole-source Fraction/RF/quotient-algebra replay."""
    import itertools,signal
    import classification_plane_m2_concurrent_double_review as C
    begun=time.monotonic();signal.alarm(540);resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jp='research/classification/jobs/plane_concurrent_m3_nonore_separate_review.json'
    opath='reports/CLASSIFICATION_PLANE_CONCURRENT_M3_NONORE_SEPARATE_REVIEW_20261108.json'
    job=json.loads((ROOT/jp).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    ipath='reports/CLASSIFICATION_PLANE_CONCURRENT_M3_NONORE_HIGHER_20261108.json';data=json.loads((ROOT/ipath).read_text())
    assert data['all_exact_arithmetic_passed']and hashlib.sha256(data['producer_source_bytes'].encode()).hexdigest()==data['producer_sha256']
    B1,B2=q.sources();read=lambda M:C.matrix(M)
    def tensor(v):
        P={ij:{}for ij in q.PAIRS}
        for(ij,e),a in zip(q.COLS,v):
            if a:P[ij][e]=a
        return q.skew(P)
    # Canonical complete raw36/Euler6 sources can differ only in recorded order.
    oldB1=read(data['full30_weight1_basis']);oldB2=read(data['full9_weight2_basis'])
    assert sorted(map(tuple,oldB1))==sorted(tuple(q.vector(P))for P in B1)
    assert sorted(map(tuple,oldB2))==sorted(tuple(q.vector(P))for P in B2)
    B1=list(map(tensor,oldB1));B2=list(map(tensor,oldB2));basis=B1+B2
    assert data['full350_row_labels']==[dict(quad=list(quad),monomial=list(e))for quad,e in q.ROWS]
    quadratic={tuple(row['pair']):row for row in data['universal_full39_quadratic_Jacobi_and_actualD']}
    assert len(quadratic)==276
    for i,j in itertools.combinations_with_replacement(range(39),2):
        P=basis[i]if i==j else q.plus(basis[i],basis[j]);J=q.jvector(q.jac(P));D=q.actual(P)
        if any(J)or D:
            row=quadratic.pop((i,j));assert {str(k):str(v)for k,v in enumerate(J)if v}==row['all_nonzero_Jacobi_rows'];assert D==C.actualpoly(row['actualD'])
    assert not quadratic
    a,b,c,u,v=[q.a(1,i)for i in range(5)];f=q.mul(q.mul(a,b),q.add(a,q.scale(b,-1)));Pf=q.skew({(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)})
    Xall=[q.a(1,i)for i in range(5)]
    eulers=[q.vector(q.skew({(i,j):q.add(q.a(1,i,l)if j==k else{},q.a(-1,j,l)if i==k else{})for i,j in q.PAIRS}))for k,l in itertools.product(range(5),repeat=2)]
    levi=[(k,l)for k,l in itertools.product(range(5),repeat=2)if(k<3)==(l<3)]
    def models(label,t):
        z=q.add(q.mul(c,u),q.mul(a,v));M=b;N=q.add(a,q.scale(b,-1));targets=[]
        if label in['J','J_ell_M','J_ell_M_eps1','J_ell_N_eps1']:
            A=q.mul(c,q.add(a,q.scale(b,t)))if label=='J'else q.mul(c,b)if label in['J_ell_M','J_ell_M_eps1']else q.mul(c,N)
            if label=='J_ell_M_eps1':A=q.add(A,q.mul(N,N))
            if label=='J_ell_N_eps1':A=q.add(A,q.mul(M,M))
            Z=[{},q.scale(q.der(A,2),-3),q.scale(q.der(A,1),3)]
            X=[{},{},q.add(a,q.scale(b,-2))];Y=[{},{},q.add(q.scale(a,4),q.scale(b,-8))];k=F()
            zz=q.add(q.mul(A,u),q.mul(q.mul(M,N),v));targets=[q.a(1,0,0,3,3,3),q.mul(q.mul(u,u),zz)]
        else:
            jordan='_Jordan_'in label
            if label.startswith('Tj_'):j=int(label.split('_')[1]);main='Tj'
            else:j=None;main=label
            Z=[q.add(q.scale(a,1+t),q.scale(b,-3)),q.scale(b,t-2),q.scale(c,1-2*t)]
            if jordan:Z[2]=q.add(Z[2],N if t==0 else M)
            if main=='F01':
                X=[{},{},q.add(a,b)];Y=[{},{},q.add(q.scale(a,1-3*t),q.scale(b,4))];k=t+1
                targets=[q.mul(q.mul(q.mul(M,u),u),z)]
                if t==-1:targets.append(q.mul(q.mul(N,N),q.a(1,3,3,3)))
                if t==0:targets.append(q.mul(q.mul(M,N),q.a(1,3,3,3)))
                if t==1:targets +=[q.mul(q.mul(M,M),q.a(1,3,3,3)),q.mul(u,q.mul(z,z))]
            elif main=='F10':
                X=[{},{},q.add(q.scale(a,3*t-2),b)];Y=[{},{},q.add(q.scale(a,-2),q.scale(b,4))];k=4*t-2
                targets=[q.mul(u,q.mul(z,z))]
                if t==0:targets +=[q.mul(q.mul(N,N),q.a(1,3,3,3)),q.mul(q.mul(N,z),q.a(1,3,3))]
                if t==1:targets +=[q.mul(q.mul(M,M),q.a(1,3,3,3)),q.mul(q.mul(M,z),q.a(1,3,3))]
                if t==F(1,2):targets.append(q.mul(q.mul(M,N),q.a(1,3,3,3)))
            else:
                assert main=='Tj';X=[{},{},q.add(q.scale(a,3*j-3*t-2),b)];Y=[{},{},q.add(q.scale(a,3*j-6*t-2),q.scale(b,4))];k=3*j-2*t-2
                h=q.a(1)
                for _ in range(j):h=q.mul(h,M)
                for _ in range(2-j):h=q.mul(h,N)
                targets=[q.mul(h,q.a(1,3,3,3))]
                if not jordan:
                    if t==j:targets.append(q.mul(N,q.mul(q.mul(u,u),z)))
                    if t==j-1:targets.append(q.mul(M,q.mul(q.mul(u,u),z)))
                    if t==F(j,2):targets.append(q.mul(u,q.mul(z,z)))
        P=q.plus(Pf,q.skew({(i,3):q.mul(u,X[i])for i in range(3)}),q.skew({(i,4):q.add(q.mul(u,Z[i]),q.mul(v,Y[i]))for i in range(3)}),q.skew({(3,4):q.scale(q.mul(u,u),k)}))
        return P,targets
    def poly(v):return C.actualpoly(v)
    def verify_case(row,label,t):
        P,targets=models(label,t);assert q.vector(P)==[C.scalar(v)for v in row['full150_P0']]
        assert not q.actual(P)and not any(q.jac(P).values())
        mons=list(map(tuple,row['all24_D3_monomials']));assert len(mons)==24 and len(set(mons))==24 and set(mons)=={e for e in q.monomials(5)if e[3]+e[4]==3}
        target=[poly(v)for v in row['whole_parallel_targets']]
        # Saved ordering is immaterial, but reconstruct the whole row space both ways.
        T=[[p.get(e,F())for e in mons]for p in target];Ann,mt=C.certificate(T,row['target_span_certificate']);assert Ann==read(row['full_line_or_space_annihilator'])
        expectedT=[[p.get(e,F())for e in mons]for p in targets]
        for test in expectedT:assert all(not sum((x*y for x,y in zip(test,h)),F())for h in Ann)
        # Distinct leading constructions have exactly the independently enumerated dimensions.
        assert len(targets)==row['target_span_certificate']['rank']
        assert C.determinant([[expectedT[i][j]for j in row['target_span_certificate']['columns']]for i in range(len(targets))])
        J1=q.transpose([q.jvector(q.jac(q.plus(P,V)))for V in B1]);J2=q.transpose([q.jvector(q.jac(q.plus(P,V)))for V in B2])
        assert J1==read(row['full350by30_J1'])and J2==read(row['full350by9_J2'])
        Dlin=[q.actual(q.plus(P,V))for V in basis];assert Dlin==[poly(v)for v in row['all_full39_actualD_linear_terms']]
        D3=[[p.get(e,F())for p in Dlin[:30]]for e in mons];assert D3==read(row['full24by30_D3_map'])
        M=J1+[[sum((a*b for a,b in zip(h,col)),F())for col in zip(*D3)]for h in Ann];assert M==read(row['full_weight1_plus_target_matrix'])
        K1,m1=C.certificate(M,row['weight1_certificate']);K2,m2=C.certificate(J2,row['weight2_certificate'])
        orbit=[]
        for k,l in levi:
            O=q.skew({(i,j):q.add(P.get((l,j),{})if i==k else{},P.get((i,l),{})if j==k else{},q.scale(q.mul(Xall[l],q.der(P.get((i,j),{}),k)),-1))for i,j in q.PAIRS});orbit.append(q.vector(O))
        G=q.transpose(eulers+orbit+[q.vector(P)]);assert G==read(row['full150_Euler_Levi_scale_matrix']);KG,mg=C.certificate(G,row['normalizer_certificate'])
        assert row['graded_Levi_projective_stabilizer_dimension']==37-row['normalizer_certificate']['rank']
        assert row['actual_filtered_group_fiber_dimension']==43-row['normalizer_certificate']['rank']
        bad=[];image=[]
        for k in K1:
            V=q.skew({ij:q.add(*(q.scale(B.get(ij,{}),a)for B,a in zip(B1,k)if a))for ij in q.PAIRS})
            bad.append([V.get((i,3),{}).get((0,0,0,0,2),F())for i in range(3)]);image.append([sum((a*b for a,b in zip(r,k)),F())for r in D3])
        assert bad==read(row['full_bad_C_maps_on_complete_kernel'])and not any(any(r)for r in bad)
        assert image==read(row['full_D3_image_on_complete_kernel'])
        assert all(e[3]for ij,p in P.items()if 3 in ij for e in p)
        assert not any(p for V in B2 for ij,p in V.items()if 3 in ij)
        assert row['all_bad_C_zero_for_every_solution']and row['G24_closed']
        dim=len(K1)+len(K2);upper=24+row['source_parameter_dimension']+dim-row['actual_filtered_group_fiber_dimension']
        assert dim==row['fixed_P0_higher_necessary_dimension_upper']and upper==row['whole_projective_image_dimension_upper']
        return dict(type=row['type'],field=row['field'],whole_target_rank=len(targets),K1=len(K1),K2=len(K2),J1_minor=repr(m1),J2_minor=repr(m2),group_minor=repr(mg),all_original350_24_150_rows_reconstructed=True,complete_kernels_and_units=True,whole_bad_u_row_zero=True,entire_actualD_u_factor=True,whole_irreducible_source_empty=True,diagnostic_dimension_upper=upper),[mt,m1,m2,mg]+[x for K in[Ann,K1,K2,KG]for line in K for x in line]
    records=[];allwallcount=0
    for parent in data['cases']:
        label=parent['type'];variable=parent['source_parameter_dimension']==1
        if variable:t=T
        elif '_Jordan_'in label:t=F(label.split('_')[-1])
        else:t=F()
        r,coeffs=verify_case(parent,label,t);records.append(r)
        if variable:
            factors=[C.monic(v)for v in parent['all_generic_minor_and_kernel_factors']]
            savedwall=[C.monic(w['minimal_polynomial'])for w in parent['all_parameter_walls']];assert set(savedwall)==set(factors)and len(savedwall)==len(set(savedwall))
            for scalar_value in coeffs:
                if isinstance(scalar_value,RF):C.supported(scalar_value.q,factors)
            for scalar_value in coeffs[:4]:
                if scalar_value:C.supported(RF(scalar_value).p,factors)
            jumps=dict(F01=[-1,0,1],F10=[0,F(1,2),1],Tj_0=[0,1,-1],Tj_1=[1,0,F(1,2)],Tj_2=[0,2,1],J=[])[label]
            for x in jumps:assert (-F(x),F(1))in factors
            for wall in parent['all_parameter_walls']:
                mod=C.monic(wall['minimal_polynomial'])
                if len(mod)==2:
                    t=-mod[0];assert t==F(wall['parameter'])
                else:
                    assert mod==(F(3),F(3),F(1))and wall['all_complex_embeddings_covered'];C.Quotient.modulus=mod;t=C.Quotient((0,1))
                rr,_=verify_case(wall,label,t);records.append(rr);allwallcount+=1
        print(label,'full350/24/150/parallel/allwalls/wholeufactor PASS',flush=True)
    assert len(data['cases'])==15 and len(records)==56 and allwallcount==41
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all56_nonore_whole_irreducible_sources_empty=True,cases=records,
        all_original_full39_universal_quadratic_pairs_reconstructed=True,whole_geometry_and_all_actual_normalizer_inverses_read=True,
        all_parameter_roots_and_quadratic_field_embeddings_covered=True,corrected_free_r_convention_checked_from_M=True,
        scope='Only allconcurrent m3 nonOre firstsymbol sources; entire higher irreducibleD empty by fullu factor. Ore/c-independent SF remainingconstruction and globalG-E remain separate.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),no_producer_or_CAS_imports=True,seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/opath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,all56_empty=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


def corrected_triple90_main():
    """Independent whole eleven-source and three Ore-module reconstruction."""
    import itertools,signal
    from classification_minimal_quartic_review import LaurentRing
    import classification_plane_double_line_m3_ore_review as L
    begun=time.monotonic();signal.alarm(540)
    resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jobpath='research/classification/jobs/plane_nonreduced_m3_triple_corrected90_separate_review.json'
    outputpath='reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED90_SEPARATE_REVIEW_20261108.json'
    job=json.loads((ROOT/jobpath).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    paths=['reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json','reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_20261108.json']
    first,old_literal=[json.loads((ROOT/p).read_text())for p in paths]
    for d in[first,old_literal]:assert d['all_exact_arithmetic_passed']and hashlib.sha256(d['producer_bytes'].encode()).hexdigest()==d['producer_sha256']
    names=('a','b','c','z','w');Y=(0,1,3,4);zero=(0,)*5;X=[q.a(1,i)for i in range(5)]
    @functools.lru_cache(None)
    def parse(text,at=None):
        def walk(n):
            if isinstance(n,ast.Constant):return{zero:F(n.value)}if n.value else{}
            if isinstance(n,ast.Name):
                if n.id=='tau':return{zero:T if at is None else at}
                return X[names.index(n.id)]
            if isinstance(n,ast.UnaryOp):return q.scale(walk(n.operand),-1)if isinstance(n.op,ast.USub)else walk(n.operand)
            a=walk(n.left)
            if isinstance(n.op,ast.Pow):
                out={zero:F(1)}
                for _ in range(n.right.value):out=q.mul(out,a)
                return out
            b=walk(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            assert isinstance(n.op,ast.Div)and set(b)=={zero};return q.scale(a,1/b[zero])
        return walk(ast.parse(str(text),mode='eval').body)
    def sc(v,at=None):
        if at is None:return scalar(str(v))
        u=RF(scalar(str(v)));ev=lambda pp:sum((c*at**i for i,c in enumerate(pp)),F());return ev(u.p)/ev(u.q)
    def mat(v,at=None):return[[sc(t,at)for t in row]for row in v]
    def vec(v,at=None):return[sc(t[0]if isinstance(t,list)else t,at)for t in v]
    def savedP(v,at=None):return q.skew({ast.literal_eval(ij):parse(text,at)for ij,text in v.items()})
    def op(P,V):
        base=q.jvector(q.jac(P));selfq=q.jvector(q.jac(V));t=q.jvector(q.jac(q.plus(P,V)))
        return[compact(a-b-c)for a,b,c in zip(t,base,selfq)]
    def cert(M,C,at=None):
        rr,cc=C['minor_rows'],C['minor_columns'];r=C['rank'];n=len(M[0]);assert len(rr)==len(cc)==r and len(set(rr))==len(set(cc))==r
        d=determinant([[M[i][j]for j in cc]for i in rr]);assert d and d==sc(C['minor_determinant'],at)
        K=mat(C['complete_kernel'],at);free=[j for j in range(n)if j not in cc];assert len(free)==C['nullity']and len(K)==n
        assert[K[i]for i in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
        for row in M:
            nz=[(i,a)for i,a in enumerate(row)if a]
            for k in range(len(free)):assert not sum((a*K[i][k]for i,a in nz if K[i][k]),F())
        return K,d
    def rftext(c):
        if not isinstance(c,RF):return str(c)
        poly=lambda vv:'+'.join('('+str(a)+')*tau**'+str(i)for i,a in enumerate(vv)if a)or'0'
        return '('+poly(c.p)+')/('+poly(c.q)+')'
    def evalP(P,t):
        return{e:sc(rftext(c),t)for e,c in P.items()if sc(rftext(c),t)}
    def shrink(p,roots):
        a=list(RF(p).p)
        for r in roots:
            while len(a)>1 and not sum((c*r**i for i,c in enumerate(a)),F()):
                b=[F()]* (len(a)-1);b[-1]=a[-1]
                for i in range(len(a)-2,0,-1):b[i-1]=a[i]+r*b[i]
                assert a[0]+r*b[0]==0;a=b
        assert len(a)<=1 and a,('uncovered polynomial roots',a,roots)
    def cover_denominators(objects,roots):
        for v in objects:
            u=RF(v);shrink(RF(u.q),roots)
    def dividepoly(P,D):
        assert D;R=P.copy();lead=max(D);coef=D[lead]
        while R:
            e=max(R);assert all(a>=b for a,b in zip(e,lead)),('factor mismatch',P,D)
            exp=tuple(a-b for a,b in zip(e,lead));c=R[e]/coef
            R=q.add(R,q.scale(q.mul({exp:c},D),-1))
        return True
    B1,B2=q.sources();B1=B1[:18]+B1[24:27]+B1[18:21]+B1[27:30]+B1[21:24]
    base2=[e[:3]for e in q.monomials(2)if not any(e[3:])];base3=[e[:3]for e in q.monomials(3)if not any(e[3:])]
    mon3=[e+(j,3-j)for e in base2 for j in range(4)];mon2=[e+(j,2-j)for e in base3 for j in range(3)]
    assert list(map(list,mon3))==first['all24_D3_monomials']and list(map(list,mon2))==first['all30_K_monomials']
    R=LaurentRing(('tau',),names,invert=('a','b'));rp=R.parse
    def reconstruct_forms(c,at=None):
        fs='a**3'if 'TRIPLE'in c['type']else'a**2*b';f=parse(fs,at)
        P={(0,1):q.der(f,2),(0,2):q.scale(q.der(f,1),-1),(1,2):q.der(f,0)}
        for form in c['full_connection_inverse_forms']:
            i,j=form['i'],form['j'];beta=[rp(t)for t in form['beta']];om=[rp(t)for t in form['horizontal_polynomial_form']]
            assert om==[R.add(R.scale(R.mul(rp(fs),v),3),R.scale(R.derivative(rp(fs),k),-int(i==j)))for k,v in enumerate(beta)]
            assert not R.add(*(R.mul(rp(names[k]),om[k])for k in range(3)))
            om=[parse(t,at)for t in form['horizontal_polynomial_form']]
            V=[q.scale(q.add(q.der(om[2],1),q.scale(q.der(om[1],2),-1)),F(1,3)),q.scale(q.add(q.der(om[0],2),q.scale(q.der(om[2],0),-1)),F(1,3)),q.scale(q.add(q.der(om[1],0),q.scale(q.der(om[0],1),-1)),F(1,3))]
            assert V==[parse(t,at)for t in form['inverse_tracefree_field']]
            inv=[q.add(q.mul(X[2],V[1]),q.scale(q.mul(X[1],V[2]),-1)),q.add(q.mul(X[0],V[2]),q.scale(q.mul(X[2],V[0]),-1)),q.add(q.mul(X[1],V[0]),q.scale(q.mul(X[0],V[1]),-1))]
            assert inv==om
            for k in range(3):P[k,3+i]=q.add(P.get((k,3+i),{}),q.scale(q.mul(V[k],X[3+j]),-1))
        pp=next(v for v in c['full_connection_inverse_forms']if v['i']==v['j']==0)
        qq=next(v for v in c['full_connection_inverse_forms']if v['i']==v['j']==1)
        rr=next(v for v in c['full_connection_inverse_forms']if v['i']==1 and v['j']==0)
        bp,bq,br=([rp(t)for t in v['beta']]for v in[pp,qq,rr])
        for i,j in itertools.combinations(range(3),2):
            for bv in[bp,bq]:assert not R.add(R.derivative(bv[j],i),R.scale(R.derivative(bv[i],j),-1))
            assert not R.add(R.derivative(br[j],i),R.scale(R.derivative(br[i],j),-1),R.scale(R.mul(R.add(bq[i],R.scale(bp[i],-1)),br[j]),-1),R.mul(R.add(bq[j],R.scale(bp[j],-1)),br[i]))
        return q.skew(P)
    def shear(P):
        cols=[]
        for normal in(3,4):
            for k in range(3):
                G=[{}for _ in range(5)];G[k]=X[normal]
                O=q.skew({(i,j):q.add(*(q.scale(q.mul(G[l],q.der(P.get((i,j),{}),l)),-1)for l in range(5)),*(q.add(q.mul(P.get((l,j),{}),q.der(G[i],l)),q.mul(P.get((i,l),{}),q.der(G[j],l)))for l in range(5)))for i,j in q.PAIRS})
                S=[q.add(q.scale(X[3],O.get((i,3),{}).get((0,0,0,2,0),F())),q.scale(X[4],O.get((i,4),{}).get((0,0,0,0,2),F())))for i in range(3)]+[{},{}]
                O=q.plus(O,q.wedge(X,S));raw=q.vector(O)
                out=[]
                for V in B1:
                    loc=next((ij,e)for ij,p in V.items()if ij in q.PAIRS for e in p);out.append(O.get(loc[0],{}).get(loc[1],F()))
                assert raw==q.vector(q.skew({ij:q.add(*(q.scale(V.get(ij,{}),a)for V,a in zip(B1,out)if a))for ij in q.PAIRS}));cols.append(out)
        return q.transpose(cols)
    first_results=[];first_kernel={};whole_minor_polys=[];reused_main=0;reused_total=0
    old_by_name={c['type']:c for c in old_literal['cases']};new_by_name={c['type']:c for c in first['cases']}
    for parent in first['cases']:
        label=parent['type']
        if label.startswith('E'):
            expectedB=10 if 'EXACT_B'in label or 'LOG_GENERIC'in label else 9 if 'ETA0_'in label or 'ETA1_'in label else 11
        else:
            index=int(label.rsplit('_',1)[1]);expectedB=([4,6,7,8]if'ZERO'in label else[7,8,9]if'LINEAR'in label else[8,9])[index]
        assert parent['whole_fixed_plane_symbol_dimension_upper']==expectedB
        genericP=reconstruct_forms(parent)
        variants=[(parent,None)]+[(dict(w,type=parent['type']),F(w['alpha']))for w in parent.get('all_exceptional_sources',[])]
        roots=[F(w['alpha'])for w in parent.get('all_exceptional_sources',[])]
        if roots:shrink(sc(parent['all_exception_polynomial']),roots)
        for c,t in variants:
            Pbase=genericP if t is None else q.skew({ij:evalP(p,t)for ij,p in genericP.items()})
            if c.get('full_special_source_elsewhere'):
                e=int(parent['type'][1]);eta=F(4,3)-e-t;assert eta in[-1,0,1,2]and'LOG_GENERIC'in parent['type']
                kind=parent['type'].rsplit('_',1)[1];index={'ZERO':{-1:2,0:0,1:0,2:3},'DX':{-1:3,0:1,1:1,2:4},'Y':{-1:0,0:4,1:2,2:2}}[kind][eta]
                fixed=new_by_name['E'+str(e)+'_LOG_ETA'+str(eta)+'_'+str(index)]
                scale=F(1,2)if kind=='DX'and eta==-1 else F(-1)if kind=='DX'and eta==2 else F(1)
                FP=reconstruct_forms(fixed);trans=q.skew({ij:{ex:v*scale**(int(ij[0]==4)+int(ij[1]==4)-ex[4])for ex,v in pp.items()}for ij,pp in FP.items()if ij in q.PAIRS})
                assert q.vector(Pbase)==q.vector(trans)
                first_results.append(dict(type=parent['type'],alpha=str(t),redirected_to_full_special_fixed_source=fixed['type'],actual_normal_w_scale=str(scale),full_higher_fiber_bijection=True))
                continue
            D=q.actual(Pbase);KM=[]
            for nm in[(3,3),(3,4),(4,4)]:
                DD=q.actual(q.plus(Pbase,q.skew({(3,4):q.a(1,*nm)})))
                KM.append([compact(DD.get(e,F())-D.get(e,F()))for e in mon2])
            KM=q.transpose(KM);constant=[D.get(e,F())for e in mon2]
            assert KM==mat(c['full30_global_K_matrix'],t)and constant==vec(c['full30_global_K_constant'],t)
            _,kd=cert(KM,c['K_certificate'],t);_,kad=cert([r+[-a]for r,a in zip(KM,constant)],c['K_augmented_certificate'],t)
            if c.get('global_K_inconsistent'):
                assert c['K_augmented_certificate']['rank']>c['K_certificate']['rank'];first_results.append(dict(type=parent['type'],alpha=str(t),global_K_inconsistent=True));continue
            P=savedP(c['full150_P0'],t);Kp=parse(c['unique_global_K'],t)
            assert all(P.get(ij,{})==Pbase.get(ij,{})for ij in q.PAIRS if ij!=(3,4));assert P[3,4]==Kp
            assert not q.actual(P)and not any(q.jac(P).values())
            J1=q.transpose([op(P,V)for V in B1]);J2=q.transpose([op(P,V)for V in B2])
            assert J1==mat(c['full350_weight1_operator'],t)and J2==mat(c['full350_weight2_operator'],t)
            D3=q.transpose([[q.actual(q.plus(P,V)).get(e,F())for e in mon3]for V in B1]);assert D3==mat(c['full24by30_actual_D3_map'],t)
            target=parse(parent['exact_nonzero_leading_target'],t);tv=[target.get(e,F())for e in mon3];index=next(i for i,v in enumerate(tv)if v)
            Ann=[]
            for i in range(24):
                if i!=index:row=[F()]*24;row[i]=1;row[index]=-tv[i]/tv[index];Ann.append(row)
            assert Ann==mat(parent['full_line_annihilator'],t)
            M=J1+[[sum((a*b for a,b in zip(row,col)),F())for col in zip(*D3)]for row in Ann];assert M==mat(c['full_weight1_plus_line_matrix'],t)
            K1,m1=cert(M,c['weight1_certificate'],t);K2,m2=cert(J2,c['weight2_certificate'],t)
            image=[[sum((a*b for a,b in zip(row,col)),F())for col in zip(*K1)]for row in D3];assert image==mat(c['whole_D3_kernel_image'],t)
            om=next(v for v in parent['full_connection_inverse_forms']if v['i']==v['j']==1)
            f=parse('a**3'if'TRIPLE'in parent['type']else'a**2*b',t)
            qform=[q.add(parse(v,t),q.der(f,i))for i,v in enumerate(om['horizontal_polynomial_form'])]
            QM=[[F()]*30 for _ in range(6)]
            for i,ex in enumerate([e for e in q.monomials(2)if not any(e[3:])]):
                for j in range(3):QM[i][-3+j]=qform[j].get(ex,F())
            assert QM==mat(c['full6by30_constant_root_character_matrix'],t)
            assert all(not sum((a*b for a,b in zip(row,col)),F())for row in QM for col in zip(*K1))
            empty=not any(any(row)for row in image);assert empty==c['whole_m3_source_empty']
            upper=parent['whole_fixed_plane_symbol_dimension_upper']+5-int(t is not None and parent['whole_alpha_line'])+len(K1[0])+len(K2[0]);assert upper==c['whole_projective_image_upper']
            if t is None and roots:
                for d in[kd,kad,m1,m2]:shrink(RF(d),roots);cover_denominators([d],roots)
                cover_denominators([v for k in[K1,K2]for line in k for v in line],roots)
            badrows=all(not v for line in K1[-3:]for v in line)
            assert badrows==c['all_three_bad_normal_kernel_rows_zero']
            proper=badrows and all(e[3]for ij,pp in P.items()if 3 in ij for e in pp)
            assert proper==bool(c['proper_D_factor_by_bad_normal_rows'])
            assert c['G_closed']==(empty or proper or upper<=24)
            assert c['whole_irreducible_source_empty']==(empty or proper)
            meta=c.get('reused_exact_P0_operator_certificate')
            if meta:
                assert meta['path']==paths[1];prior=old_by_name[meta['source_type']]
                old_at=F(meta['alpha'])if meta['alpha']is not None else None
                if old_at is not None:prior=next(w for w in prior['all_exceptional_sources']if F(w['alpha'])==old_at)
                assert P==savedP(prior['full150_P0'],old_at)
                assert J1==mat(prior['full350_weight1_operator'],old_at)and J2==mat(prior['full350_weight2_operator'],old_at)
                assert D3==mat(prior['full24by30_actual_D3_map'],old_at)and M==mat(prior['full_weight1_plus_line_matrix'],old_at)
                reused_total+=1;reused_main+=int(t is None)
            first_kernel[parent['type'],t]=(K1,P)
            first_results.append(dict(type=parent['type'],alpha=str(t)if t is not None else'whole',empty_by_fullD3=empty,proper_actual_z_factor=proper,G_closed=c['G_closed'],K1=len(K1[0]),K2=len(K2[0]),image_upper=upper,complete_original_rows=True))
        print(parent['type'],'full30/9/24/allwalls PASS',flush=True)
    assert len(first['cases'])==90 and sum(c['whole_alpha_line']for c in first['cases'])==15
    assert reused_main==first['exact_main_analyses_reused']==45
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all90_corrected_sources_and110_wall_records_checked=True,full_original_wall_sources=74,actual_special_fiber_redirects=36,all_source_flatness_correct_minus_checked=True,
        cases=first_results,exact_literal_reuses_main=reused_main,exact_literal_reuses_total=reused_total,
        old_plus_wholecoverage_not_used=True,whole_G_closed_by_this_phase=False,
        pending='Correctedtable residual28 labels need validliteral14 fullsource bridge/strictsupport18modules+34walls and newfixed originalJacobi ideal factor. No oldcoverage reinstatement.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),no_producer_or_CAS_imports=True,seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/outputpath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,sources90=True,literal_reuses45=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


def corrected_triple_modules_factor_main():
    """Independent eighteen affine modules and newfixed original ideal lift."""
    import itertools,signal
    from classification_minimal_quartic_review import LaurentRing
    import classification_plane_double_line_m3_ore_review as L
    begun=time.monotonic();signal.alarm(360);resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jp='research/classification/jobs/plane_nonreduced_m3_triple_modules_factor_separate_review.json';opath='reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_MODULES_FACTOR_SEPARATE_REVIEW_20261108.json'
    job=json.loads((ROOT/jp).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    modulepath='reports/CLASSIFICATION_PLANE_M3_TRIPLE_VALID_ORE_MODULES_20261108.json';sourcepath='reports/CLASSIFICATION_PLANE_M3_TRIPLE_NEW_FIXED_ORIGINAL_SUPPORT_20261108.json';factorpath='reports/CLASSIFICATION_PLANE_M3_TRIPLE_NEW_FIXED_ACTUAL_FACTOR_20261108.json'
    data,support,factor=[json.loads((ROOT/p).read_text())for p in[modulepath,sourcepath,factorpath]]
    for d in[data,support,factor]:assert d['all_exact_arithmetic_passed']and hashlib.sha256(d['producer_bytes'].encode()).hexdigest()==d['producer_sha256']
    names=('a','b','c','z','w');Y=(0,1,3,4);zero=(0,)*5;X=[q.a(1,i)for i in range(5)]
    def sc(v,at=None):
        u=scalar(str(v))
        if at is None:return u
        u=RF(u);ev=lambda p:sum((a*at**i for i,a in enumerate(p)),F());return ev(u.p)/ev(u.q)
    def mat(m,at=None):return[[sc(v,at)for v in row]for row in m]
    def vec(m,at=None):return[sc(v[0]if isinstance(v,list)else v,at)for v in m]
    @functools.lru_cache(None)
    def parse(text,at=None):
        def walk(n):
            if isinstance(n,ast.Constant):return{zero:F(n.value)}if n.value else{}
            if isinstance(n,ast.Name):return{zero:T if at is None else at}if n.id=='tau'else X[names.index(n.id)]
            if isinstance(n,ast.UnaryOp):return q.scale(walk(n.operand),-1)if isinstance(n.op,ast.USub)else walk(n.operand)
            a=walk(n.left)
            if isinstance(n.op,ast.Pow):
                p={zero:F(1)}
                for _ in range(n.right.value):p=q.mul(p,a)
                return p
            b=walk(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            assert isinstance(n.op,ast.Div)and set(b)=={zero};return q.scale(a,1/b[zero])
        return walk(ast.parse(str(text),mode='eval').body)
    def cert(M,C,at=None):
        rr,cc=C['minor_rows'],C['minor_columns'];rank=C['rank'];n=len(M[0]);assert len(rr)==len(cc)==rank
        det=determinant([[M[i][j]for j in cc]for i in rr]);assert det and det==sc(C['minor_determinant'],at)
        K=mat(C['complete_kernel'],at);free=[j for j in range(n)if j not in cc];assert len(K)==n and len(free)==C['nullity']
        assert[K[i]for i in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
        for row in M:
            nz=[(i,a)for i,a in enumerate(row)if a]
            for j in range(len(free)):assert not sum((a*K[i][j]for i,a in nz if K[i][j]),F())
        return K,det
    def supported(poly,roots):
        aa=list(poly)
        for r in roots:
            while len(aa)>1 and not sum((a*r**i for i,a in enumerate(aa)),F()):
                b=[F()]*(len(aa)-1);b[-1]=aa[-1]
                for i in range(len(aa)-2,0,-1):b[i-1]=aa[i]+r*b[i]
                assert not aa[0]+r*b[0];aa=b
        assert len(aa)<=1 and aa
    def division(P,D):
        assert D;R=P.copy();e=max(D)
        while R:
            f=max(R);assert all(a>=b for a,b in zip(f,e));k=tuple(a-b for a,b in zip(f,e));v=R[f]/D[e];R=q.add(R,q.scale(q.mul({k:v},D),-1))
    def strict(Q):return[q.add(*(q.mul(Q.get((a,l),{}),q.der(Q.get((b,c),{}),l))for a,b,c in[t,(t[1],t[2],t[0]),(t[2],t[0],t[1])]for l in Y))for t in itertools.combinations(Y,3)]
    records=[];wall_count=0;redirects=[]
    for parent in data['cases']:
        roots=[F(w['alpha'])for w in parent.get('all_exceptional_fields',[])]
        for w in parent.get('all_exceptional_fields',[]):
            if w.get('source_excluded_here'):
                r=F(w['alpha']);kind=parent['type']
                assert(kind.startswith('E1DX')and r==F(-1,6))or(kind.startswith('E2R')and r in[F(1,3),F(-2,3)])or(kind.startswith('PHI0_RHO')and r==0)
                if kind.startswith('E2R')and r==F(-2,3):assert all(any(c['type']=='E2_D0_'+str(i)for c in data['cases'])for i in range(3))
                redirects.append(dict(type=kind,alpha=str(r),reason=w['reason']));wall_count+=1
        variants=[(parent,None)]+[(w,F(w['alpha']))for w in parent.get('all_exceptional_fields',[])if not w.get('source_excluded_here')]
        if roots:supported(RF(sc(parent['all_exception_polynomial'])).p,roots)
        for row,at in variants:
            Q=q.skew({(Y[i],Y[j]):parse(v,at)for key,v in parent['full_Q'].items()for i,j in[ast.literal_eval(key)]});assert not any(strict(Q))
            W=[parse(v,at)for v in parent['W0']];WB=[]
            for i,terms in[(0,['a*z','b*z','a*w','b*w','z**2','z*w','w**2']),(1,['a*z','b*z','a*w','b*w','z**2','z*w','w**2']),(2,['z*w','w**2']),(3,['z**2','z*w'])]:
                for t in terms:p=[{}for _ in range(4)];p[i]=parse(t,at);WB.append(p)
            pad=lambda vv:q.plus(Q,q.skew({(2,i):v for i,v in zip(Y,vv)}))
            P=pad(W);C=[compact(2*v)for v in q.jvector(q.jac(P))]
            M=q.transpose([[compact(2*(a-b))for a,b in zip(q.jvector(q.jac(pad([q.add(a,b)for a,b in zip(W,V)]))),q.jvector(q.jac(P)))]for V in WB])
            assert M==mat(row['full350_affine_matrix'],at)and C==vec(row['full350_affine_constant'],at)
            K,mi=cert(M,row['linear_certificate'],at);KA,ma=cert([r+[-c]for r,c in zip(M,C)],row['augmented_certificate'],at)
            if row['status']=='inconsistent':assert row['augmented_certificate']['rank']>row['linear_certificate']['rank'];conclusion='inconsistent'
            else:
                assert row['status']=='whole_affine_module'and row['augmented_certificate']['rank']==row['linear_certificate']['rank']
                part=vec(row['particular_solution'],at);assert all(not(a+sum((b*c for b,c in zip(r,part)),F()))for r,a in zip(M,C))
                D=q.actual(P);Dl=[q.add(q.actual(pad([q.add(a,b)for a,b in zip(W,V)])),q.scale(D,-1))for V in WB]
                Dp=q.add(D,*[q.scale(p,c)for p,c in zip(Dl,part)]);Dk=[q.add(*[q.scale(Dl[i],K[i][j])for i in range(18)])for j in range(len(K[0]))]
                assert Dp==parse(row['full_actual_D_particular'],at)and Dk==[parse(v,at)for v in row['full_actual_D_kernel_coefficients']]
                if row['whole_actual_D_zero']:assert not Dp and not any(Dk);factorpoly={};conclusion='zeroD'
                elif row['proper_D_factor']is not None:
                    factorpoly=parse(row['proper_D_factor'],at);assert factorpoly and 0<max(map(sum,factorpoly))<5
                    for p in[Dp]+Dk:division(p,factorpoly)
                    conclusion=row['proper_D_factor'];assert row['G_empty']
                else:
                    assert parent['type']=='E1DX_1'and at in[F(4,21),F(1,12)]and len(K[0])==6
                    assert row['image_dimension_offset']==17 and row['whole_projective_image_upper']==23 and not row['G_empty']
                    factorpoly={};conclusion='weak_full_source_image_upper23'
                assert row['G_closed']
                if at is None and roots:
                    for v in part+[c for line in K for c in line]+[c for p in[Dp]+Dk+[factorpoly]for c in p.values()]:supported(RF(v).q,roots)
            if at is None and roots:
                for v in[mi,ma]:supported(RF(v).p,roots);supported(RF(v).q,roots)
            if at is not None:wall_count+=1
            records.append(dict(type=parent['type'],alpha=str(at)if at is not None else'whole',full350_W18_reconstructed=True,rank=row['linear_certificate']['rank'],nullity=row['linear_certificate']['nullity'],complete_augmented_and_original_kernel=True,entire_actualD_factor_or_contradiction=conclusion,whole_irreducible_source_empty=row['G_empty']))
        print(parent['type'],'full350/W18/fullD/allwalls PASS',flush=True)
    assert len(data['cases'])==18 and wall_count==34 and len(records)==44 and len(redirects)==8
    # Newfixed source: reconstruct all15 firstkernel plus FULL9 secondsource.
    v2=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json').read_text());base=next(c for c in v2['cases']if c['type']=='E2_LOG_ETA-1_4');K=mat(base['weight1_certificate']['complete_kernel']);assert len(K)==30 and len(K[0])==15
    sourcevars=support['source_variables'];assert sourcevars==factor['all_source_variables']and len(sourcevars)==24
    R=LaurentRing(tuple(sourcevars),names);p=R.parse
    B1,B2=q.sources();B1=B1[:18]+B1[24:27]+B1[18:21]+B1[27:30]+B1[21:24]
    P=L.skew(R,{ij:p(v)for key,v in base['full150_P0'].items()for ij in[ast.literal_eval(key)]})
    for j in range(24):
        for ij in q.PAIRS:
            raw=q.add(*(q.scale(V.get(ij,{}),K[i][j])for i,V in enumerate(B1)if K[i][j]))if j<15 else B2[j-15].get(ij,{})
            terms={}
            for ex,c in raw.items():e=[0]*R.n+list(ex);e[j]=1;terms[tuple(e)]=F(c)
            P[ij]=R.add(P[ij],terms);P[ij[::-1]]=R.scale(P[ij],-1)
    assert all(P[ast.literal_eval(ij)]==p(v)for ij,v in support['full150_original_Pi'].items())
    JJ=L.projective(R,P);J=[L.coef(R,JJ[quad],ex)for quad,ex in q.ROWS]
    assert J==[p(t)for t in support['full350_original_Jacobi_equations']]
    D=L.actual(R,P);assert D==p(support['full_actual_D'])
    restricted={e:c for e,c in D.items()if not e[R.n+3]};assert restricted==p(factor['full_actual_D_restriction'])
    for row in factor['three_equations_and_original_row_lifts']:
        target=p(row['target']);rhs=R.add(*(R.mul(p(v['coefficient']),J[v['row']])for v in row['exact_original350_row_lift']));assert target==rhs
    lift=R.add(*(R.mul(p(v['coefficient']),J[v['row']])for v in factor['exact_full_original350_D_restriction_lift']));assert lift==restricted
    assert not factor['inverted_source_coordinates']and factor['whole_fixed_source_irreducible_quintic_empty']
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all18_modules_and34_wall_records_closed=True,module_cases=records,explicit_excluded_chart_redirects=redirects,
        newfixed=dict(full15_kernel_plus_all9=True,full150_raw_source_reconstructed=True,all350_original_polynomial_Jacobi_rows=True,full_actualD_reconstructed=True,all_three_original_QQ_row_lifts=True,exact_four_original_row_polynomial_D_restriction_lift=True,no_source_coordinate_inverted=True,whole_irreducible_source_empty=True),
        whole_G_closed_by_this_phase=False,pending='Validliteral14 sourcejointmaps/strictsupport and correctedtable residualsource inclusion are separate finalbridge; oldpluscoverage remainswithdrawn.',
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),no_producer_or_CAS_imports=True,seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/opath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,modules18=True,focused_factor=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


def corrected_triple_literal_bridge_main():
    """Independent literal14 Euler/joint maps and corrected-table actual inclusion."""
    import itertools,signal
    from classification_minimal_quartic_review import LaurentRing
    import classification_plane_double_line_m3_ore_review as L
    begun=time.monotonic();signal.alarm(540);resource.setrlimit(resource.RLIMIT_AS,(2048*1024**2,)*2)
    jp='research/classification/jobs/plane_nonreduced_m3_triple_literal_bridge_separate_review.json';opath='reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_LITERAL_BRIDGE_SEPARATE_REVIEW_20261108.json'
    job=json.loads((ROOT/jp).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    old90=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_20261108.json').read_text());ore=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_ORE_20261108.json').read_text());v2=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_V2_20261108.json').read_text())
    for d in[old90,ore,v2]:assert d['all_exact_arithmetic_passed']and hashlib.sha256(d['producer_bytes'].encode()).hexdigest()==d['producer_sha256']
    names=('a','b','c','z','w');Y=(0,1,3,4);zero=(0,)*5;X=[q.a(1,i)for i in range(5)]
    def sc(v,at=None):
        s=scalar(str(v))
        if at is None:return s
        u=RF(s);ev=lambda p:sum((c*at**i for i,c in enumerate(p)),F());return ev(u.p)/ev(u.q)
    def mat(M,at=None):return[[sc(v,at)for v in r]for r in M]
    @functools.lru_cache(None)
    def parse(text,at=None):
        def walk(n):
            if isinstance(n,ast.Constant):return{zero:F(n.value)}if n.value else{}
            if isinstance(n,ast.Name):return{zero:T if at is None else at}if n.id=='tau'else X[names.index(n.id)]
            if isinstance(n,ast.UnaryOp):return q.scale(walk(n.operand),-1)if isinstance(n.op,ast.USub)else walk(n.operand)
            a=walk(n.left)
            if isinstance(n.op,ast.Pow):
                p={zero:F(1)}
                for _ in range(n.right.value):p=q.mul(p,a)
                return p
            b=walk(n.right)
            if isinstance(n.op,ast.Add):return q.add(a,b)
            if isinstance(n.op,ast.Sub):return q.add(a,q.scale(b,-1))
            if isinstance(n.op,ast.Mult):return q.mul(a,b)
            assert isinstance(n.op,ast.Div)and set(b)=={zero};return q.scale(a,1/b[zero])
        return walk(ast.parse(str(text),mode='eval').body)
    def savedP(row,at=None):return q.skew({ast.literal_eval(ij):parse(v,at)for ij,v in row['full150_P0'].items()})
    def rftext(v):
        if not isinstance(v,RF):return str(v)
        poly=lambda vv:'+'.join('('+str(c)+')*tau**'+str(i)for i,c in enumerate(vv)if c)or'0'
        return '('+poly(v.p)+')/('+poly(v.q)+')'
    def polytext(P):return '+'.join('('+rftext(c)+')'+''.join('*'+names[i]+'**'+str(k)for i,k in enumerate(e)if k)for e,c in P.items())or'0'
    class FR(LaurentRing):
        def __init__(self,params,at=None):super().__init__(params,names);self.at=at
        def constant(self,c):return{self.zero:RF(c)}if c else{}
        def parse(self,text):
            def walk(n):
                if isinstance(n,ast.Constant):return self.constant(F(n.value))
                if isinstance(n,ast.Name):return self.constant(T if self.at is None else self.at)if n.id=='tau'else self.variable(n.id)
                if isinstance(n,ast.UnaryOp):return self.scale(walk(n.operand),-1)if isinstance(n.op,ast.USub)else walk(n.operand)
                a=walk(n.left)
                if isinstance(n.op,ast.Pow):
                    p=self.constant(1)
                    for _ in range(n.right.value):p=self.mul(p,a)
                    return p
                b=walk(n.right)
                if isinstance(n.op,ast.Add):return self.add(a,b)
                if isinstance(n.op,ast.Sub):return self.add(a,self.scale(b,-1))
                if isinstance(n.op,ast.Mult):return self.mul(a,b)
                assert isinstance(n.op,ast.Div);return self.mul(a,self.inverse_monomial(b))
            return walk(ast.parse(str(text),mode='eval').body)
    def cert(M,C,at=None):
        rr,cc=C['minor_rows'],C['minor_columns'];rank=C['rank'];n=len(M[0]);d=determinant([[M[i][j]for j in cc]for i in rr]);assert d and d==sc(C['minor_determinant'],at)
        K=mat(C['complete_kernel'],at);free=[j for j in range(n)if j not in cc];assert len(K)==n and len(free)==C['nullity']
        assert[K[i]for i in free]==[[F(i==j)for j in range(len(free))]for i in range(len(free))]
        for row in M:
            nz=[(i,c)for i,c in enumerate(row)if c]
            for j in range(len(free)):assert not sum((c*K[i][j]for i,c in nz if K[i][j]),F())
        return K,d
    def supported(poly,roots):
        pp=list(poly)
        for r in roots:
            while len(pp)>1 and not sum((c*r**i for i,c in enumerate(pp)),F()):
                qq=[F()]*(len(pp)-1);qq[-1]=pp[-1]
                for i in range(len(pp)-2,0,-1):qq[i-1]=pp[i]+r*qq[i]
                assert not pp[0]+r*qq[0];pp=qq
        assert len(pp)<=1 and pp
    B1,B2=q.sources();B1=B1[:18]+B1[24:27]+B1[18:21]+B1[27:30]+B1[21:24]
    mon3=list(map(tuple,old90['all24_D3_monomials']));old_by={r['type']:r for r in old90['cases']}
    selected={}
    for c in ore['cases']:
        label=c['type']
        if label.startswith('EXACT_B_FIXED_'):
            outer=old_by['E0_EXACT_B_'+label.rsplit('_',1)[1]];src=next(w for w in outer['all_exceptional_sources']if F(w['alpha'])==F(-2,3));src=dict(src,exact_nonzero_leading_target=outer['exact_nonzero_leading_target'])
        else:src=old_by[label]
        selected[label]=src
    def source_coords(R,polys):
        out=[]
        for pp in polys:
            assert all(sum(e[:R.n])==1 and not any(e[R.n:])for e in pp)
            out.append([pp.get(tuple(int(i==j)for i in range(R.n))+(0,)*5,F())for j in range(R.n)])
        return out
    def constpart(R,P):return{e:c for e,c in P.items()if not any(e[:R.n])}
    results=[];source_descriptors={}
    for parent in ore['cases']:
        label=parent['type'];src=selected[label];roots=[F(w['alpha'])for w in parent.get('all_exceptional_sources',[])];variants=[(parent,None)]+[(w,F(w['alpha']))for w in parent.get('all_exceptional_sources',[])]
        if roots:supported(RF(sc(parent['all_exception_polynomial'])).p,roots)
        for row,at in variants:
            P0=savedP(src,at);assert not q.actual(P0)and not any(q.jac(P0).values())
            J1=q.transpose([q.jvector(q.jac(q.plus(P0,V)))for V in B1]);D3=q.transpose([[q.actual(q.plus(P0,V)).get(e,F())for e in mon3]for V in B1])
            target=parse(src['exact_nonzero_leading_target'],at);tv=[target.get(e,F())for e in mon3];idx=next(i for i,v in enumerate(tv)if v);Ann=[]
            for i in range(24):
                if i!=idx:h=[F()]*24;h[i]=1;h[idx]=-tv[i]/tv[idx];Ann.append(h)
            M=J1+[[sum((a*b for a,b in zip(h,col)),F())for col in zip(*D3)]for h in Ann]
            if row.get('whole_m3_source_empty_by_original_certificate'):
                assert at is not None;rr,piv,Krows=q.rref(M);assert all(not sum((a*b for a,b in zip(r,k)),F())for r in D3 for k in Krows)
                results.append(dict(type=label,alpha=str(at),entire_original_D3_kernel_image_zero=True));continue
            K,mi=cert(M,row['full_original_M1_certificate'],at)
            if all(not c for r in K[-3:]for c in r):
                assert all(e[3]for ij,p in P0.items()if 3 in ij for e in p)
                results.append(dict(type=label,alpha=str(at),full_original_bad_z_rows_zero=True,entire_actualD_z_factor=True));continue
            n=len(K[0]);params=tuple('p'+str(i)for i in range(n))+tuple('r'+str(i)for i in range(9));R=FR(params,at);p=R.parse;xx=list(map(p,names));P=L.skew(R,{ij:p(polytext(P0.get(ij,{})))for ij in q.PAIRS})
            for j in range(n+9):
                for ij in q.PAIRS:
                    raw=q.add(*(q.scale(V.get(ij,{}),K[i][j])for i,V in enumerate(B1)if K[i][j]))if j<n else B2[j-n].get(ij,{})
                    poly=R.mul(p(polytext(raw)),p(params[j]));P[ij]=R.add(P[ij],poly);P[ij[::-1]]=R.scale(P[ij],-1)
            gauge=[R.scale(R.derivative(P.get((i,2),{}),2),-1)if i!=2 else{}for i in range(5)];assert gauge==[p(t)for t in row['full_positive_weight_Euler_gauge']]
            G=L.skew(R,{(i,j):R.add(P[i,j],R.scale(R.add(R.mul(xx[i],gauge[j]),R.scale(R.mul(xx[j],gauge[i]),-1)),-1))for i,j in q.PAIRS});assert all(not R.derivative(G[ij],2)for ij in q.PAIRS)
            Q=L.skew(R,{ij:G[ij]for ij in itertools.combinations(Y,2)});W=[G[2,i]for i in Y];cf=lambda pp,text:L.coef(R,pp,next(iter(p(text)))[R.n:])
            ell=R.add(R.scale(R.mul(cf(W[2],'z**2'),p('z')),-1),R.scale(R.mul(cf(W[3],'w**2'),p('w')),-1));assert ell==p(row['remaining_Euler_c_component'])
            W=[R.add(v,R.mul(ell,xx[i]))for v,i in zip(W,Y)];W0=[constpart(R,v)for v in W];assert W0==[p(t)for t in row['full_W0']]
            Q0={ij:constpart(R,Q[ij])for ij in Q};assert all(Q0[Y[i],Y[j]]==p(t)for key,t in row['full_Q0'].items()for i,j in[ast.literal_eval(key)])
            mon2=[e for e in q.monomials(2)if not e[2]]
            coords=[L.coef(R,R.add(Q[ij],R.scale(Q0[ij],-1)),e)for ij in itertools.combinations(Y,2)for e in mon2]
            Qmap=source_coords(R,coords);assert Qmap==mat(row['full60_Q_source_matrix'],at);QK,qmi=cert(Qmap,row['Q_space_certificate'],at);nq=row['Q_space_certificate']['rank']
            WB=L.Wbasis(R);wc=[]
            for j in range(18):i=next(i for i in range(4)if WB[j][i]);wc.append(L.coef(R,R.add(W[i],R.scale(W0[i],-1)),next(iter(WB[j][i]))[R.n:]))
            assert W==[R.add(v,*[R.mul(wc[j],WB[j][i])for j in range(18)])for i,v in enumerate(W0)]
            Wmap=source_coords(R,wc);assert Wmap==mat(row['full_W18_source_matrix'],at)
            joint=[Qmap[i]for i in row['Q_space_certificate']['minor_rows']]+Wmap;assert joint==mat(row['whole_QW_source_joint_matrix'],at)
            JK,jmi=cert(joint,row['joint_QW_certificate'],at);assert len(joint)==len(params)==row['joint_QW_certificate']['rank']and not row['joint_QW_certificate']['nullity']
            # Reconstruct the canonical whole Q-space basis from independent original coefficients.
            RC=FR(tuple('Q'+str(i)for i in range(nq))+tuple('W'+str(i)for i in range(18)),at);pc=RC.parse;cols=row['Q_space_certificate']['minor_columns'];expected={}
            for k,ij in enumerate(itertools.combinations(Y,2)):
                poly=pc('0')
                for j,col in enumerate(cols):
                    for z,e in enumerate(mon2):
                        coefficient=Qmap[10*k+z][col]
                        if coefficient:
                            term=RC.mul(pc('Q'+str(j)),pc(rftext(coefficient)))
                            for i,power in enumerate(e):term=RC.mul(term,L.power(RC,pc(names[i]),power))
                            poly=RC.add(poly,term)
                expected[ij]=poly
            QC=L.skew(RC,expected);assert all(QC[Y[i],Y[j]]==pc(t)for key,t in row['full_affine_Q'].items()for i,j in[ast.literal_eval(key)])
            strict=L.strict(RC,QC);assert all(strict[tuple(Y[i]for i in ast.literal_eval(key))]==pc(t)for key,t in row['full_strict_Q_Jacobi'].items())
            WC=[pc(t)for t in row['full_W0']];WBC=L.Wbasis(RC);WC=[RC.add(w,*[RC.mul(pc('W'+str(j)),WBC[j][i])for j in range(18)])for i,w in enumerate(WC)]
            DD=L.actual(RC,L.pad(RC,QC,WC));assert DD==pc(row['full_actual_D_polynomial'])
            common=row.get('whole_Q_common_constant_v');common_ok=False
            if common is not None:
                vv=[RC.constant(sc(t,at))for t in common];assert any(vv)
                assert all(not RC.add(RC.mul(vv[i],QC[Y[j],Y[k]]),RC.scale(RC.mul(vv[j],QC[Y[i],Y[k]]),-1),RC.mul(vv[k],QC[Y[i],Y[j]]))for i,j,k in itertools.combinations(range(4),3))
                assert all(not RC.add(*(RC.mul(vv[k],RC.derivative(QC[ij],Y[k]))for k in range(4)))for ij in itertools.combinations(Y,2));common_ok=True
            results.append(dict(type=label,alpha=str(at)if at is not None else'whole',full30_plus_all9_Euler_source=True,Q_rank=nq,full_joint_bijection_dimension=len(params),joint_minor=repr(jmi),full_strict_Q_and_actualD_reconstructed=True,whole_common_constant_vector=common_ok))
            if at is None and roots:
                for value in[mi,qmi,jmi]:supported(RF(value).p,roots);supported(RF(value).q,roots)
                for value in[c for KK in[K,QK,JK]for line in KK for c in line]+[sc(c,at)for c in common or[]]:supported(RF(value).q,roots)
        print(label,'fullrawEuler/QWjoint/strict/common/allwalls PASS',flush=True)
    # Actual normal diagonal maps of every unclosed corrected point to valid literal families.
    def transformed(P,s):return q.skew({ij:{e:c*s**(int(ij[0]==4)+int(ij[1]==4)-e[4])for e,c in pp.items()}for ij,pp in P.items()if ij in q.PAIRS})
    inclusions=[];focused=0
    for parent in v2['cases']:
        variants=[(parent,None)]+[(w,F(w['alpha']))for w in parent.get('all_exceptional_sources',[])if not w.get('full_special_source_elsewhere')]
        for row,at in variants:
            if row['G_closed']:continue
            label=parent['type'];scale=F(1);old_at=at
            if label=='E2_LOG_ETA-1_4':focused+=1;inclusions.append(dict(type=label,alpha=str(at),target='newfixed_four_original_rows'));continue
            if label in selected:target=label
            elif 'EXACT_B'in label:target='EXACT_B_FIXED_'+label.rsplit('_',1)[1];old_at=None
            else:
                e=int(label[1]);eta=int(label.split('_ETA')[1].split('_')[0]);index=int(label.rsplit('_',1)[1]);old_at=F(4,3)-e-eta
                if (eta in[0,1]and index==0)or(eta==2 and index==3):target='E2_LOG_GENERIC_ZERO';assert e==2
                else:
                    assert(eta in[0,1]and index==1)or(eta==-1 and index==3)or(eta==2 and index==4)
                    target='E'+str(e)+'_LOG_GENERIC_DX';scale=F(1,2)if eta==-1 else F(-1)if eta==2 else F(1)
            assert target in selected;newP=savedP(row,at);oldP=savedP(selected[target],old_at)
            assert q.vector(transformed(newP,scale))==q.vector(oldP)
            newtarget=parse(parent['exact_nonzero_leading_target'],at);oldtarget=parse(selected[target]['exact_nonzero_leading_target'],old_at);assert newtarget==oldtarget
            inclusions.append(dict(type=label,alpha=str(at)if at is not None else'whole',target=target,target_alpha=str(old_at)if old_at is not None else'fixed_or_whole',actual_normal_w_scale=str(scale),whole_fullhigher_fiber_bijection=True))
    assert focused==1
    phase1=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_CORRECTED90_SEPARATE_REVIEW_20261108.json').read_text());phase2=json.loads((ROOT/'reports/CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_MODULES_FACTOR_SEPARATE_REVIEW_20261108.json').read_text());assert phase1['all_independent_arithmetic_passed']and phase2['all_independent_arithmetic_passed']
    result=dict(task_id='G',all_independent_arithmetic_passed=True,all_validliteral14_source_maps_and_walls_reconstructed=True,literal_sources=results,all_corrected_unclosed_actual_source_inclusions=inclusions,
        strict_reduced_support_and_actual_normalizer_bridge='Independent handproof in final separate note; allliteral strictQ polynomials/fullQWsource were reconstructed here. No oldpluscoverage used.',
        all_corrected_nonreduced_m3_generictriple_G_image_upper24=True,two_E1DX1_weakwalls_E_unresolved=True,whole_pointwise_E_covered=False,
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),no_producer_or_CAS_imports=True,seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with(ROOT/opath).open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_independent_arithmetic_passed=True,whole_G24_closed=True,seconds=result['seconds'],RSS=result['maximum_RSS_KiB'])))


if __name__=='__main__':
    import sys
    if '--corrected-triple-literal-bridge'in sys.argv:corrected_triple_literal_bridge_main()
    elif '--corrected-triple-modules-factor'in sys.argv:corrected_triple_modules_factor_main()
    elif '--corrected-triple90'in sys.argv:corrected_triple90_main()
    elif '--concurrent-nonore-m3'in sys.argv:concurrent_nonore_m3_main()
    elif '--triple-line-sf01'in sys.argv:triple_line_sf01_main()
    elif '--moving-simple-m3'in sys.argv:moving_simple_m3_main()
    elif '--moving-double-m3'in sys.argv:moving_double_m3_main()
    elif '--reduced-m3'in sys.argv:reduced_m3_main()
    elif '--squarefree'in sys.argv:squarefree_main()
    else:main()
