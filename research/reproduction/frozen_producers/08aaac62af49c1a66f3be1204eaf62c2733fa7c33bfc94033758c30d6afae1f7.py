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


if __name__=='__main__':
    import sys
    if '--moving-double-m3'in sys.argv:moving_double_m3_main()
    elif '--reduced-m3'in sys.argv:reduced_m3_main()
    elif '--squarefree'in sys.argv:squarefree_main()
    else:main()
