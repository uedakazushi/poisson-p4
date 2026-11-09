#!/usr/bin/env python3
"""Independent whole free-quintic source, true-cut support and prior103 audit.

Fraction sparse arithmetic only. Producer/legacy engine/SymPy are not imported.
The nonzero quadratic cokernel coefficient is the weakest support certificate;
no cubic or completed graph is computed after it succeeds.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,product
from math import lcm
from pathlib import Path
import hashlib,json,re,resource,signal,subprocess,time
import classification_degenerate_quadric_q4_nil31_I_review as q
from classification_minimal_quartic_review import LaurentRing,determinant
from classification_ga9_cuspidal_quartic_review import simplify_edges
import classification_ga9_distinct_component_review as G

ROOT=Path(__file__).resolve().parents[1]
PREFIX='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_'
SOURCE=PREFIX+'SAITO_20261107.json';AMBIENT=PREFIX+'AMBIENT_20261107.json';SUPPORT=PREFIX+'COMPONENT_SUPPORT_20261107.json';AUDIT=PREFIX+'PRIOR103_20261107.json'
PROJECTION='research/classification/inputs/essential_three_quintic_prior103_20261107.json'
OUT=PREFIX+'REVIEW_20261107.json';JOB=PREFIX+'REVIEW_JOB_20261107.json'


def sparse_vector(v,n):return [F(v.get(str(i),'0')) for i in range(n)]
def sparse_rows(rows,n):return [sparse_vector(row,n) for row in rows]
def matmul(A,B):return [[sum((a*b for a,b in zip(row,col) if a and b),F(0)) for col in zip(*B)] for row in A]
def plus(*Ps):return q.antisymmetric({ij:q.add(*(P.get(ij,{}) for P in Ps)) for ij in q.PAIRS})
def tensor(v):return q.antisymmetric({ij:{ex+(0,0,0):c for (pair,ex),c in zip(q.COLS,v) if pair==ij and c} for ij in q.PAIRS})
def minor(M,C):
    rr,cc=C['rows'],C['columns'];assert len(rr)==len(cc)==C['rank']
    value=q.determinant([[M[i][j] for j in cc] for i in rr]);assert value and value==F(C.get('determinant',C.get('minor_determinant')))
    return str(value)


def whole_source(row):
    params=('g00','g01','g10','g11')+tuple('r'+str(i) for i in range(6))
    R=LaurentRing(params,('x','y','z','u','v'));p=R.parse;x=[R.variable(n) for n in ('x','y','z','u','v')]
    f=p('x**5+x**2*y**3+y**4*z');assert f==p(row['polynomial'])
    mon=[p(t) for t in ('x**2','x*y','x*z','y**2','y*z','z**2')]
    target=sorted((e for e in product(range(7),repeat=3) if sum(e)==6),reverse=True)
    grad=[R.derivative(f,i) for i in range(3)]
    cols=[R.mul(m,g) for g in grad for m in mon]
    M=[[P.get((0,)*10+e+(0,0),F(0)) for P in cols] for e in target]
    saved=row['complete_quadratic_gradient_kernel'];assert M==[[F(v) for v in r] for r in saved['all_28_by_18_rows']]
    assert target==[tuple(e) for e in saved['target_exponents']]
    K=[[F(v) for v in k] for k in saved['full_18coefficient_basis']]
    assert len(K)==2 and q.rank(K)[0]==2 and not any(any(r) for r in matmul(M,q.transpose(K)))
    rank,cc=q.rank(M);_,rr=q.rank(q.transpose(M));assert rank==16
    det=q.determinant([[M[i][j] for j in cc] for i in rr]);assert det
    fields=[[R.add(*(R.scale(mon[j],k[6*i+j]) for j in range(6))) for i in range(3)] for k in K]
    assert fields==[[p(t) for t in v] for v in saved['vector_fields']]
    saito=determinant(R,[[x[i],fields[0][i],fields[1][i]] for i in range(3)])
    coefficient=F(saved['determinant_over_F']);assert coefficient==F(-9,64)
    assert saito==p(saved['saito_determinant'])==R.scale(f,coefficient)
    U=[R.add(R.mul(p('g00'),u),R.mul(p('g01'),v)) for u,v in zip(*fields)]
    V=[R.add(R.mul(p('g10'),u),R.mul(p('g11'),v)) for u,v in zip(*fields)]
    div=lambda W:R.add(*(R.derivative(W[i],i) for i in range(3)))
    magnetic=R.add(*(R.mul(p('r'+str(i)),m) for i,m in enumerate(mon)))
    data={(i,3):U[i] for i in range(3)};data.update({(i,4):V[i] for i in range(3)})
    data[3,4]=R.add(R.scale(R.mul(div(V),x[3]),-1),R.mul(div(U),x[4]),magnetic)
    P={}
    for ij in q.PAIRS:
        P[ij]=data.get(ij,{});P[ij[::-1]]=R.scale(P[ij],-1)
    source=row['whole_cotangent_source']
    assert all(P[ij]==p(source['all_ten_brackets'][str(ij[0])+','+str(ij[1])]) for ij in q.PAIRS)
    for i,j,k in combinations(range(5),3):
        J=R.add(*(R.mul(P.get((a,l),{}),R.derivative(P.get((b,c),{}),l)) for a,b,c in ((i,j,k),(j,k,i),(k,i,j)) for l in range(5)))
        assert not J
    assert all(not R.add(*(R.derivative(P.get((i,j),{}),j) for j in range(5))) for i in range(5))
    assert all(not R.add(*(R.mul(P.get((i,j),{}),R.derivative(f,j)) for j in range(5))) for i in range(5))
    pf=[]
    for missing in range(5):
        i,j,k,l=[a for a in range(5) if a!=missing]
        pf.append(R.add(R.mul(P[i,j],P[k,l]),R.scale(R.mul(P[i,k],P[j,l]),-1),R.mul(P[i,l],P[j,k])))
    D=R.add(*(R.scale(R.mul(x[i],v),(-1)**i) for i,v in enumerate(pf)))
    assert pf==[p(t) for t in source['all_five_raw_pfaffians']]
    assert D==p(source['actual_quintic'])==R.scale(R.mul(p('g00*g11-g01*g10'),f),-coefficient)
    # Absolute irreducibility: primitive degree-one polynomial in z.
    assert R.derivative(f,2)==p('y**4')
    assert all(not e[10+2] or e[10+2]==1 for e in f)
    constant={e:c for e,c in f.items() if e[12]==0}
    assert {e:c for e,c in constant.items() if e[11]==0}==p('x**5')
    values=(1,0,0,1,0,0,0,0,0,0)
    def evaluate(poly):
        out={}
        for ex,c in poly.items():
            for power,v in zip(ex[:10],values):c*=v**power
            e=ex[10:]+(0,0,0);out[e]=out.get(e,F(0))+c
        return {e:c for e,c in out.items() if c}
    point=q.antisymmetric({ij:evaluate(P[ij]) for ij in q.PAIRS})
    scale=lcm(*(c.denominator for ij in q.PAIRS for c in point.get(ij,{}).values()))
    point=q.antisymmetric({ij:q.scale(point.get(ij,{}),scale) for ij in q.PAIRS})
    return point,dict(full28_by18_kernel_rank=16,complete_two_syzygies=True,independent16minor=str(det),
                      all_frame4_and_magnetic6_parameters_retained=True,all10strictJacobi_divergence_Casimir=True,
                      determinant_over_F=str(coefficient),all_rawPfaffians_actualD=True,absolute_irreducibility='Primitive z-linear with coefficient y4, constant at y0=x5 over C',integer_scale=scale)


def prior103(audit):
    projection=json.loads((ROOT/PROJECTION).read_text());meta=projection['immutable_catalogue']
    blob=subprocess.check_output(['git','cat-file','blob',meta['git_blob']],cwd=ROOT)
    assert hashlib.sha256(blob).hexdigest()==meta['raw_sha256']
    assert subprocess.check_output(['git','rev-parse',meta['commit']+':'+meta['path']],cwd=ROOT,text=True).strip()==meta['git_blob']
    catalogue={r['id']:r for r in json.loads(blob)['components']};rows=projection['all103_comparison_rows'];assert len(catalogue)==len(rows)==103
    assert rows==audit['all103_comparison_rows']
    for row in rows:assert all(catalogue[row['id']][k]==v for k,v in row.items())
    mps=[]
    for c in audit['all29_same23_MPS_actual_SNC_witnesses']:
        path=ROOT/c['accepted_source_report'];assert hashlib.sha256(path.read_bytes()).hexdigest()==c['source_sha256']
        raw=json.loads(path.read_text())
        periods=sorted(set(re.findall(r'\bp[0-9]+\b',' '.join(raw['homogeneous_bivector'].values()))))
        assert set(c['full_period_parameters'])<=set(periods)
        R=LaurentRing(periods,tuple(f'x{i}' for i in range(5)),invert=periods)
        P={}
        for pair,text in raw['homogeneous_bivector'].items():
            i,j=map(int,pair.split(','));P[i,j]=R.parse(simplify_edges(text));P[j,i]=R.scale(P[i,j],-1)
        D={}
        for missing in range(5):
            i,j,k,l=[a for a in range(5) if a!=missing]
            pf=R.add(R.mul(P.get((i,j),{}),P.get((k,l),{})),R.scale(R.mul(P.get((i,k),{}),P.get((j,l),{})),-1),R.mul(P.get((i,l),{}),P.get((j,k),{})))
            D=R.add(D,R.scale(R.mul(R.variable('x'+str(missing)),pf),(-1)**missing))
        kappa=R.parse(c['kappa']);assert kappa and all(not any(e[-5:]) for e in kappa)
        assert D==R.parse(c['actual_SNC_quintic'])==R.mul(kappa,*[R.variable('x'+str(i)) for i in range(5)])
        assert R.mul(kappa,kappa,kappa,kappa,kappa)==R.parse(c['nonzero_gradient_five_minor'])
        mps.append(c['id'])
    assert len(set(mps))==29
    irred=[]
    for c in audit['both_old_irreducible_quintic_actual_rank5_witnesses']:
        path=ROOT/c['accepted_report'];assert hashlib.sha256(path.read_bytes()).hexdigest()==c['source_sha256'];raw=json.loads(path.read_text())
        if c['id']=='R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23':
            tensor_data=raw['canonical_eigen_coordinate_bivector'];P=q.antisymmetric({tuple(map(int,k.split(','))):G.modelparse(v) for k,v in tensor_data.items()})
        else:
            row=next(r for r in raw['nil22_long_quartic_curve_recognition']['all_five_actual_curve_recognitions'] if r['anchor']==5)
            P=q.antisymmetric({tuple(map(int,k.split(','))):G.parse(v) for k,v in row['primitive_eigen_tensor'].items()})
        D=G.quintic(P);R=LaurentRing((),tuple(f'x{i}' for i in range(5)))
        expected=R.parse(c['actual_quintic']);assert {e[:5]:v for e,v in D.items()}==expected
        minor=q.determinant([[q.der(D,i).get(tuple(e)+(0,0,0),F(0)) for i in range(5)] for e in c['minor_rows']])
        assert minor and minor==F(c['nonzero_gradient_five_minor']);irred.append(c['id'])
    witness=set(mps+irred);reasons=[]
    for row in rows:
        if row['generic_rank']!=4:reason='different generic rank'
        elif row['projective_dimension']!=23:reason='different full component dimension'
        elif row['id'] in witness:reason='actual essential5 source point outside closed new essential<=3 locus'
        else:
            divisor=row['generic_actual_degeneracy_divisor'];assert divisor['irreducible_factor_degrees']!=[5] or not divisor.get('reduced',True)
            reason='accepted generic reducible/nonreduced divisor forces closed multiplication locus, excludes new irreducible F'
        reasons.append(dict(id=row['id'],reason=reason))
    return dict(immutable_catalogue=meta,all103_rows_verified_from_Git_blob=True,all29_actual_fullperiod_SNC_sources=mps,
                both_accepted_actual_essential5_quintics=irred,all103_distinctness_reasons=reasons,all_prior103_distinct=True)


def main():
    start=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    destination=ROOT/OUT
    if destination.exists():raise FileExistsError(destination)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    src=json.loads((ROOT/SOURCE).read_text())['candidates'][0];amb=json.loads((ROOT/AMBIENT).read_text())['full_ambient'];sup=json.loads((ROOT/SUPPORT).read_text())['component_support']
    P,source_review=whole_source(src);assert source_review['integer_scale']==int(amb['canonical_point_integer_scale'])==80
    R0=LaurentRing((),('x','y','z','u','v'))
    legacy=q.antisymmetric({tuple(map(int,k.split(','))):{e+(0,0,0):c for e,c in R0.parse(v).items()} for k,v in amb['canonical_point'].items()})
    assert q.vec(P)==q.vec(legacy) and not any(q.jacobi(P).values())
    D=q.transpose([q.jvec(q.jacobi(plus(P,q.antisymmetric({ij:{e+(0,0,0):F(1)}})))) for ij,e in q.COLS])
    assert D==sparse_rows(amb['full_D_sparse_rows'],150)
    # All actual Euler/GL/scalar columns, independently assembled.
    Ecols=[];Ocols=[]
    for i,j in product(range(5),repeat=2):
        E=q.antisymmetric({(a,b):q.add(q.atom(1,a,j) if b==i else {},q.atom(-1,b,j) if a==i else {}) for a,b in q.PAIRS})
        O=q.antisymmetric({(a,b):q.add(q.mul(q.atom(1,j),q.der(P.get((a,b),{}),i)),q.scale(P.get((j,b),{}) if a==i else {},-1),q.scale(P.get((a,j),{}) if b==i else {},-1)) for a,b in q.PAIRS})
        Ecols.append(q.vec(E));Ocols.append(q.vec(O))
    M=q.transpose(Ecols+Ocols+[q.vec(P)])
    assert M==q.transpose([sparse_vector(v,150) for v in amb['all51_actual_direction_columns']])
    assert not any(any(row) for row in matmul(D,M)) and q.rank(M)[0]==48
    mminor=minor(M,amb['actual_Euler_GL_scalar_source'])
    cut=amb['true_cut'];N=sparse_vector(cut['all150_true_normals'][0],150);R=tensor(N)
    assert not any(row[0] for row in matmul(D,[[v] for v in N]))
    sigma=amb['Cstar']['coordinate_weights'];assert sigma==[0,0,0,1,1]
    bw=[sum(sigma[i] for i in ij)-sum(a*b for a,b in zip(e,sigma))-1 for ij,e in q.COLS]
    ow=[sum(sigma[i] for i in quad)-sum(a*b for a,b in zip(e,sigma))-2 for quad,e in q.ROWS]
    assert bw==amb['Cstar']['relative_input_weights'] and ow==amb['Cstar']['relative_output_weights']
    assert all(not v or bw[j]==ow[i] for i,row in enumerate(D) for j,v in enumerate(row))
    assert all(not v or bw[i]==0 for i,v in enumerate(q.vec(P))) and all(not v or bw[i]==-1 for i,v in enumerate(N))
    assert min(bw)==-3 and max(bw)==1 and cut['normal_weights']==[-1]
    Fcols=[sparse_vector(v,150) for v in cut['actual48_basis']];Fmat=q.transpose(Fcols)
    cursor=0
    for block in amb['full_weight_blocks']:
        for j in block['actual_source_column_indices']:
            assert Fcols[cursor]==[M[i][j] if i in set(block['input_columns']) else F(0) for i in range(150)]
            assert Fcols[cursor]==[M[i][j] for i in range(150)];cursor+=1
    assert cursor==48
    cutrows=cut['coordinate_rows'];assert len(set(cutrows))==48 and all(not N[i] for i in cutrows)
    detcut=q.determinant([Fmat[i] for i in cutrows]);assert detcut==F(cut['source_cut48_determinant']) and detcut
    dc=cut['independent_implicit_coordinate_indices'];dr=sup['implicit_rows'];assert dc==sup['implicit_columns'] and len(dc)==len(dr)==101
    A=[[D[i][j] for j in dc] for i in dr];inverse=[[F(v) for v in r] for r in sup['inverse_101_by_101']]
    eye=[[F(i==j) for j in range(101)] for i in range(101)]
    assert matmul(A,inverse)==eye and matmul(inverse,A)==eye
    # Directsum determinant with101 coordinate columns reduces to49x49.
    remaining=[i for i in range(150) if i not in dc];perm=remaining+dc
    sign=(-1)**sum(perm[i]>perm[j] for i in range(150) for j in range(i+1,150))
    split=sign*q.determinant([Fmat[i]+[N[i]] for i in remaining]);assert split==F(cut['full150_direct_sum_determinant']) and split
    assert q.rank(q.transpose(Fcols+[N]))[0]==49
    # True one-variable quadratic obstruction; stop here after nonzero pairing.
    assert len(sup['stages'])==1 and sup['stages'][0]['degree']==2
    stage=sup['stages'][0];Q=q.jvec(q.jacobi(R));assert Q==sparse_vector(stage['raw_quadratic_source'],350)
    T=sparse_vector(stage['unique_selected_IFT_correction'],150)
    assert all(not v or i in dc and bw[i]==-2 for i,v in enumerate(T))
    residual=[sum((a*b for a,b in zip(row,T) if a and b),F(0))+v for row,v in zip(D,Q)]
    assert residual==sparse_vector(stage['full350_residual'],350) and all(not residual[i] for i in dr)
    cert=sup['weakest_support_certificate'];ell=sparse_vector(cert['cokernel_functional'],350)
    assert len([v for v in ell if v])==11 and all(not v for v in matmul([ell],D)[0])
    pairing=sum(a*b for a,b in zip(ell,Q));assert pairing==F(cert['nonzero_initial_coefficient'])==F(-57897,220448)
    assert cert['initial_degree']==2 and pairing==residual[cert['residual_row']]
    grad=G.quintic(P);assert q.rank([[q.der(grad,i).get(ex+(0,0,0),F(0)) for i in range(5)] for ex in q.monomials(4)])[0]==3
    audit=prior103(json.loads((ROOT/AUDIT).read_text())['prior103_audit'])
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,source_review=source_review,
        full350_by150_D_reconstructed=True,D_rank=101,actual_Euler_GL_scalar_rank=48,actual_nonzero48minor=mminor,
        all51_actual_columns_reconstructed=True,full150_directsum_determinant=str(split),true48cut_and_normal=True,
        coordinate_weights=sigma,normalized_tensor_character=1,normal_character=-1,prior_IFT_degree_upper=3,
        whole_GL2_frame_plus6mag_one_orbit='Accepted uniform essential3 theorem gives full GL2+Rsource the same23orbit; no extra param fixed.',
        quadratic_source_all350=True,all350_unique_selected_correction_residual=True,cokernel_all150_D_zero=True,
        exact_quadratic_nonzero_pairing=str(pairing),stopped_after_weakest_nonzero_quadratic=True,
        reduced_truecut_support='origin',projective_orbitclosure_is_full_component=True,projective_component_dimension=23,
        componenthood_scope='FormalIFT solves101 independent Jacobi equations on genuine48cut, leavingone normal. Nonzero initial n² residual makes reducedcutdim0. Krull upperraw48 matchesactualorbit/gauge48, so full23orbitclosure is a component. No extra cubic ornilpotence claim.',
        prior103_audit=audit,new_distinct_component_supported=True,no_global_count_changed=True,
        no_producer_legacy_or_CAS_imports=True,inputs=job['exact_inputs_and_hashes'],
        reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)
def eight_ambient_review():
    """Independent all-seven full ambient proof; old first104 code is untouched."""
    start=time.monotonic();signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    source='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCES_20261107.json'
    ambient='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_AMBIENT_20261107.json'
    jobpath='research/classification/jobs/essential_three_quintic_eight_ambient_review.json'
    destination=ROOT/'reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_AMBIENT_REVIEW_20261107.json'
    if destination.exists():raise FileExistsError(destination)
    job=json.loads((ROOT/jobpath).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    sources={c['id']:c for c in json.loads((ROOT/source).read_text())['eight_curve_sources']['all_eight_candidate_sources']}
    all_data=json.loads((ROOT/ambient).read_text())['eight_full_ambient'];records=[]
    R0=LaurentRing((),('x','y','z','u','v'))
    convert=lambda p:{e+(0,0,0):v for e,v in R0.parse(p).items()}
    for case in all_data['all_eight_cases']:
        if case['id']=='C45':
            for key in ('immutable_full_ambient','immutable_component_support'):
                b=case[key];assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256']
            continue
        a=case['full_ambient'];row=sources[case['id']];fields=[[convert(t) for t in v] for v in row['complete_quadratic_gradient_kernel']['vector_fields']]
        U,V=fields;div=lambda W:q.add(*(q.der(W[i],i) for i in range(3)))
        p={(i,3):U[i] for i in range(3)};p.update({(i,4):V[i] for i in range(3)})
        p[3,4]=q.add(q.scale(q.mul(div(V),q.atom(1,3)),-1),q.mul(div(U),q.atom(1,4)))
        scale=int(a['canonical_point_integer_scale']);P=q.antisymmetric({ij:q.scale(v,scale) for ij,v in p.items()})
        expected=q.antisymmetric({tuple(map(int,k.split(','))):convert(t) for k,t in a['canonical_point'].items()})
        assert q.vec(P)==q.vec(expected) and not any(q.jacobi(P).values())
        D=q.transpose([q.jvec(q.jacobi(plus(P,q.antisymmetric({ij:{e+(0,0,0):F(1)}})))) for ij,e in q.COLS])
        assert D==sparse_rows(a['full_D_sparse_rows'],150)
        Ecols=[];Ocols=[]
        for i,j in product(range(5),repeat=2):
            E=q.antisymmetric({(u,v):q.add(q.atom(1,u,j) if v==i else {},q.atom(-1,v,j) if u==i else {}) for u,v in q.PAIRS})
            O=q.antisymmetric({(u,v):q.add(q.mul(q.atom(1,j),q.der(P.get((u,v),{}),i)),q.scale(P.get((j,v),{}) if u==i else {},-1),q.scale(P.get((u,j),{}) if v==i else {},-1)) for u,v in q.PAIRS})
            Ecols.append(q.vec(E));Ocols.append(q.vec(O))
        M=q.transpose(Ecols+Ocols+[q.vec(P)])
        assert M==q.transpose([sparse_vector(v,150) for v in a['all51_actual_direction_columns']])
        assert q.rank(M)[0]==48 and not any(any(row) for row in matmul(D,M))
        dm=minor(D,a['full_350_by_150_D']);mm=minor(M,a['actual_Euler_GL_scalar_source']);rank=a['full_350_by_150_D']['rank']
        assert rank in (101,102)
        sigma=a['Cstar']['coordinate_weights'];assert sigma==[0,0,0,1,1]
        bw=[sum(sigma[i] for i in ij)-sum(c*d for c,d in zip(e,sigma))-1 for ij,e in q.COLS]
        ow=[sum(sigma[i] for i in quad)-sum(c*d for c,d in zip(e,sigma))-2 for quad,e in q.ROWS]
        assert bw==a['Cstar']['relative_input_weights'] and ow==a['Cstar']['relative_output_weights']
        assert all(not v or bw[j]==ow[i] for i,rowD in enumerate(D) for j,v in enumerate(rowD))
        assert all(not v or bw[i]==0 for i,v in enumerate(q.vec(P)))
        cut=a['true_cut'];normals=[sparse_vector(v,150) for v in cut['all150_true_normals']]
        assert len(normals)==102-rank
        Fcols=[sparse_vector(v,150) for v in cut['actual48_basis']];Fmat=q.transpose(Fcols)
        cursor=0
        for block in a['full_weight_blocks']:
            for j in block['actual_source_column_indices']:
                assert Fcols[cursor]==[M[i][j] for i in range(150)]
                assert all(not v or bw[i]==block['weight'] for i,v in enumerate(Fcols[cursor]));cursor+=1
        assert cursor==48 and len(set(cut['coordinate_rows']))==48
        cm=q.determinant([Fmat[i] for i in cut['coordinate_rows']]);assert cm==F(cut['source_cut48_determinant']) and cm
        for n in normals:
            assert not any(rowD[0] for rowD in matmul(D,[[v] for v in n]))
            assert all(not v or bw[i]==-1 for i,v in enumerate(n)) and all(not n[i] for i in cut['coordinate_rows'])
        dc=cut['independent_implicit_coordinate_indices'];assert len(set(dc))==rank
        remaining=[i for i in range(150) if i not in dc];perm=remaining+dc
        sign=(-1)**sum(perm[i]>perm[j] for i in range(150) for j in range(i+1,150))
        split=sign*q.determinant([Fmat[i]+[n[i] for n in normals] for i in remaining])
        assert split==F(cut['full150_direct_sum_determinant']) and split
        assert q.rank(q.transpose(Fcols+normals))[0]==150-rank
        pairing=None
        if rank==101:
            w=case['weakest_one_normal_quadratic'];assert dc==w['IFT_columns'] and w['true_normal_character']==-1
            dr=w['IFT_rows'];A=[[D[i][j] for j in dc] for i in dr];inverse=[[F(v) for v in r] for r in w['all101inverse_rows']]
            eye=[[F(i==j) for j in range(101)] for i in range(101)]
            assert matmul(A,inverse)==eye and matmul(inverse,A)==eye
            assert q.determinant(A)==F(w['nonzero101minor'])
            Q=q.jvec(q.jacobi(tensor(normals[0])));assert Q==sparse_vector(w['full350_Jacobi_quadratic_source'],350)
            T=sparse_vector(w['unique_cut_quadratic_correction'],150)
            assert all(not v or i in dc and bw[i]==-2 for i,v in enumerate(T))
            residual=[sum((x*y for x,y in zip(rowD,T) if x and y),F(0))+v for rowD,v in zip(D,Q)]
            assert residual==sparse_vector(w['all350_true_quadratic_residual'],350) and all(not residual[i] for i in dr)
            ell=sparse_vector(w['cokernel_functional'],350);assert all(not v for v in matmul([ell],D)[0])
            pairing=sum(x*y for x,y in zip(ell,Q));assert pairing and pairing==F(w['nonzero_pairing'])==residual[w['residual_row']]
        grad=G.quintic(P);assert q.rank([[q.der(grad,i).get(ex+(0,0,0),F(0)) for i in range(5)] for ex in q.monomials(4)])[0]==3
        records.append(dict(id=case['id'],full350_by150_reconstructed=True,full51_actual_columns_reconstructed=True,D_rank=rank,M_rank=48,D_minor=dm,M_minor=mm,true48cut_full150_directsum=str(split),quadratic_pairing=None if pairing is None else str(pairing),projective_component_dimension=23,smooth=rank==102,no_higher_obstruction_run=True))
    assert len(records)==7
    meta=all_data['immutable_prior104_catalogue'];blob=subprocess.check_output(['git','cat-file','blob',meta['git_blob']],cwd=ROOT)
    assert hashlib.sha256(blob).hexdigest()==meta['raw_sha256'] and len(json.loads(blob)['components'])==104
    assert subprocess.check_output(['git','rev-parse',meta['commit']+':'+meta['path']],cwd=ROOT,text=True).strip()==meta['git_blob']
    binding=all_data['accepted_prior103_audit'];assert hashlib.sha256((ROOT/binding['path']).read_bytes()).hexdigest()==binding['sha256']
    projection=all_data['prior104_projection'];assert hashlib.sha256((ROOT/projection['path']).read_bytes()).hexdigest()==projection['sha256']
    oldreview=json.loads((ROOT/OUT).read_text());assert oldreview['all_passed'] and oldreview['prior103_audit']['all_prior103_distinct']
    types=[row['expected_curve_type'] for row in sources.values()];assert len(set(types))==8
    pairs={tuple(sorted((r['first'],r['second']))) for r in all_data['all28_source_closure_distinctness_reasons']}
    assert pairs=={tuple(sorted(p)) for p in combinations(sources,2)}
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,all_seven_full_ambient_proofs=records,immutable_prior104=meta,prior103_accepted_independent_audit_reused='reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_REVIEW_20261107.json',all7_new_vs104='All prior103 full-component reasons apply uniformly to essential3 irreducible23components; previous104 maximum cusp multiplicity4 differs from all seven max<=3.',all28_mutual_distinctness='Intrinsic essential-three plane quotient and distinct verified cuspconfigurations; two dense23orbits cannot have the same23closure.',curve_classification_and_full_source_recognition='Separate author proof/input obligation; this review confirms all ambient facts and the conditional all104/mutual comparison, not the independent Puiseux completeness audit.',no_new_count_registered=True,no_producer_CAS_imports=True,inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.monotonic()-start,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    destination.write_text(json.dumps(result,indent=2)+'\n');print('PASS all7fullambient',round(result['seconds'],3),'seconds',flush=True)


if __name__=='__main__':
    import sys
    if '--eight-ambient' in sys.argv:eight_ambient_review()
    else:main()
