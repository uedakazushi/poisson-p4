#!/usr/bin/env python3
"""Source/support preflight and one exact polynomial constrained83 block build.

No kernel substitution, selector, Gaussian, Gram, source point, or prime.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import resource
import signal
import time

ROOT = Path(__file__).resolve().parents[1]
ID = 'S61_SPARSE_GENERIC5_CONSTRAINED83_RAW_BLOCK_20261104_v63'


def read(p):
    b = p.read_bytes()
    return json.loads(gzip.decompress(b) if p.suffix == '.gz' else b)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def poly(p):
    return {tuple(e): Q(c) for e,c in p if Q(c)}


def plus_support(a,b):
    return {tuple(x+y for x,y in zip(e,f)) for e in a for f in b}


def source():
    folder = ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63'
    known = read(folder/'known-five-input.json.gz')
    raw = read(folder/'pulled134-input.json.gz')
    c = [poly(p) for p in known['q_center']]
    d = [poly(p) for p in known['Q52_minus_center']]
    n = [[poly(p) for p in row] for row in known['N45x20']]
    l = [[poly(p) for p in row] for row in known['L51x20']]
    rows = []
    for terms in raw['equations']:
        eq = {}
        for mon,e,v in terms:
            p = eq.setdefault(tuple(mon),{})
            e = tuple(e)
            p[e] = p.get(e,Q(0)) + Q(v)
        rows.append({mon:{e:v for e,v in p.items() if v} for mon,p in eq.items()})
    assert len(rows)==134 and len(c)==len(d)==len(n)==45 and len(l)==51
    return known,raw,c,d,n,l,rows


def prepare(out,work):
    started = time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    signal.alarm(60)
    out.mkdir(parents=True,exist_ok=False)
    work.mkdir(parents=True,exist_ok=False)
    known,raw,c,d,n,l,rows = source()
    # The old five-column frame is literally the proposed unconstrained frame.
    old_dir = ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_FIVE_SELECTOR_20261104_v63'
    matrix = read(old_dir/'five-column-matrix.json.gz')
    assert matrix['column_order']==['P*F(C)','J(C)*W0','J(C)*W1','polarQ(W0,d)','polarQ(W1,d)']
    old_counts = [sum(len(row[j]['terms']) for row in matrix['matrix134x5']) for j in range(5)]
    assert sum(old_counts)==692953
    del matrix
    pivots = [3,6,13]
    normal = [i for i in range(45) if i not in pivots]
    counts = dict(Fcenter=0,JcenterN=0,polar_d_N=0,raw_polar_d_nonpivot42=0)
    maxima = dict.fromkeys(counts,0)
    support_rows = []
    pairs = dict.fromkeys(counts,0)
    for eq in rows:
        fc=set(); jc=[set() for _ in range(45)]; pd=[set() for _ in range(45)]
        for mon,p in eq.items():
            if not mon:
                fc.update(p)
            elif len(mon)==1:
                i,=mon; fc.update(plus_support(p,c[i])); jc[i].update(p)
                pairs['Fcenter']+=len(p)*len(c[i])
            else:
                assert len(mon)==2
                i,j=mon
                fc.update(plus_support(plus_support(p,c[i]),c[j]))
                jc[i].update(plus_support(p,c[j]));jc[j].update(plus_support(p,c[i]))
                pd[i].update(plus_support(p,d[j]));pd[j].update(plus_support(p,d[i]))
                pairs['Fcenter']+=len(p)*len(c[i])*len(c[j])
        ja=[set() for _ in range(20)]; pb=[set() for _ in range(20)]
        for i in range(45):
            for j in range(20):
                ja[j].update(plus_support(jc[i],n[i][j]))
                pb[j].update(plus_support(pd[i],n[i][j]))
                pairs['JcenterN']+=len(jc[i])*len(n[i][j])
                pairs['polar_d_N']+=len(pd[i])*len(n[i][j])
        blocks={'Fcenter':[fc],'JcenterN':ja,'polar_d_N':pb,
                'raw_polar_d_nonpivot42':[pd[i] for i in normal]}
        support_rows.append({name:sum(map(len,p)) for name,p in blocks.items()})
        for name,p in blocks.items():
            counts[name]+=sum(map(len,p))
            maxima[name]=max(maxima[name],max(map(len,p)))
    lterms = sum(len(p) for row in l for p in row)
    rawterms = sum(len(p) for row in rows for p in row.values())
    input_terms = rawterms+sum(map(len,c))+sum(map(len,d))+sum(len(p) for row in n+l for p in row)
    output_upper = sum(counts.values())+2*lterms
    aggregate = input_terms+output_upper
    record={
        'id':ID,'phase':'source_support_preflight_only',
        'coefficient_field':'QQ(H,L,s0,q0,t0)',
        'actual_matrix_shape':[236,83],
        'column_order':'[c,alpha20,beta20,n42]',
        'row_blocks':[134,51,51],
        'top_formula':'[F(C), J(C)N, polar_Q(d,N), polar_Q(d,e_normal42)]',
        'bottom_formula':'[0,L,0,0] then [0,0,L,0]',
        'normal_q_rows42':normal,
        'source_terms':dict(raw134=rawterms,center=sum(map(len,c)),d=sum(map(len,d)),N=sum(len(p) for row in n for p in row),L=lterms,total=input_terms),
        'top_block_support_upper':counts,'max_single_entry_support_upper':maxima,
        'row_support_upper':support_rows,'partial_scalar_pair_upper':pairs,
        'two_L_copies_terms':2*lterms,'output_total_support_upper':output_upper,
        'input_plus_output_support_upper':aggregate,
        'old_literal_five_frame_columns':old_counts,
        'old_692953_frame_job_aborted_as_exact_duplicate':True,
        'constant_unit_prefix_scope':'Only L51x20; NOT this new constrained236x83 matrix.',
        'no_kernel_T_or_P_or_Delta_substitution':True,
        'no_selector_Gaussian_Gram_point_prime_or_GB':True,
        'hard_seconds':60,'hard_memory_bytes':2_000_000_000,
        'hard_artifact_bytes':32*1024*1024,'hard_live_expanded_terms':150000,
        'arithmetic_authorized_gate_passed':aggregate+10000<=150000,
        'seconds':time.monotonic()-started}
    ipaths=[ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63/known-five-input.json.gz',ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63/pulled134-input.json.gz',old_dir/'five-column-matrix.json.gz',ROOT/'reports/S61_SPARSE_GENERIC5_ALL670_SOURCE_BINDING_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_SPARSE_GENERIC5_CROSS_RANK_COVARIANCE_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_SPARSE_GENERIC5_PLANE_RANK3_CIRCUIT_20261104_v63_REVIEW/report.json',ROOT/'reports/S61_SPARSE_GENERIC5_CONSTANT_UNIT_PREFIX_20261104_v63/report.json']
    for p in ipaths:
        if p.name=='report.json':assert read(p)['all_passed']
    record['input_sha256']={str(p.relative_to(ROOT)):sha(p) for p in ipaths}
    (out/'source-support-preflight.json').write_text(json.dumps(record,indent=2)+'\n')
    if record['arithmetic_authorized_gate_passed']:
        inp=work/'source-input.txt'
        with inp.open('w') as f:
            def put(p):
                f.write(str(len(p))+'\n')
                for e,v in sorted(p.items()):f.write(str(v)+' '+' '.join(map(str,e))+'\n')
            f.write('5 134 45 20 51\n')
            for p in c+d:put(p)
            for matrix in (n,l):
                for row in matrix:
                    for p in row:put(p)
            for eq in rows:
                f.write(str(len(eq))+'\n')
                for mon,p in sorted(eq.items()):
                    f.write(' '.join(map(str,[len(mon),*mon]))+'\n');put(p)
        assert inp.stat().st_size<32*1024*1024
        record['work_input_sha256']=sha(inp)
        record['work_input_bytes']=inp.stat().st_size
        (out/'source-support-preflight.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['phase','source_terms','top_block_support_upper','two_L_copies_terms','output_total_support_upper','input_plus_output_support_upper','arithmetic_authorized_gate_passed','seconds']}))


def finish(out,work):
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    signal.alarm(60)
    assert not (out/'report.json').exists()
    pre=read(out/'source-support-preflight.json')
    assert pre['arithmetic_authorized_gate_passed']
    for p,h in pre['input_sha256'].items():assert sha(ROOT/p)==h,p
    assert sha(work/'source-input.txt')==pre['work_input_sha256']
    native=read(out/'actual-236x83-matrix.json')
    assert native['all_passed'] and native['shape']==[236,83]
    assert native['all_QQ_coefficient_identity_checks']==10988
    assert native['peak_live_expanded_terms']<=150000
    entries={}
    colterms=[0]*83;rowterms=[0]*236
    constant=[];monomial=[]
    for i,j,p in native['sparse_entries']:
        assert (i,j) not in entries and 0<=i<236 and 0<=j<83
        assert Q(p['content']) and p['terms']
        assert len({tuple(t[1:]) for t in p['terms']})==len(p['terms'])
        assert all(int(t[0]) and len(t)==6 and min(t[1:])>=0 for t in p['terms'])
        entries[i,j]=p;colterms[j]+=len(p['terms']);rowterms[i]+=len(p['terms'])
        if len(p['terms'])==1:
            monomial.append([i,j])
            if not any(p['terms'][0][1:]):constant.append([i,j])
    assert len(entries)==native['nonzero_entries']==9431
    assert sum(colterms)==native['total_matrix_terms']==49039
    source_folder=ROOT/'reports/S61_SPARSE_GENERIC5_INPUT_PREFLIGHT_20261104_v63'
    known=read(source_folder/'known-five-input.json.gz')
    for i,row in enumerate(known['L51x20']):
        for j,encoded in enumerate(row):
            expected=poly(encoded)
            for target in [(134+i,1+j),(185+i,21+j)]:
                actual=entries.get(target)
                actual={} if actual is None else {tuple(t[1:]):Q(actual['content'])*int(t[0]) for t in actual['terms']}
                assert actual==expected
    assert all(i<134 or (134<=i<185 and 1<=j<=20) or (185<=i<236 and 21<=j<=40) for i,j in entries)

    def support_matching(edges):
        # A graph only: not Gaussian elimination or a determinant assertion.
        neighbors=[[] for _ in range(83)]
        for row,col in edges:neighbors[col].append(row)
        occupied={}
        def visit(col,seen):
            for row in neighbors[col]:
                if row in seen:continue
                seen.add(row)
                if row not in occupied or visit(occupied[row],seen):
                    occupied[row]=col
                    return True
            return False
        for col in range(83):visit(col,set())
        return sorted([[row,col] for row,col in occupied.items()])
    support={
        'literal_QQ_constant_entries':len(constant),'literal_parameter_monomial_entries':len(monomial),
        'QQ_constant_entries':constant,'parameter_monomial_entries':monomial,
        'constant_support_matching':support_matching(constant),
        'monomial_support_matching':support_matching(monomial),
        'matching_is_NOT_rank_or_an_executable_elimination_sequence':True,
        'column_parameter_terms':colterms,'row_parameter_terms':rowterms,
        'no_constant_or_monomial_elimination_executed':True}
    (out/'actual-support-and-unit-prefix-frontier.json').write_text(json.dumps(support,indent=2)+'\n')
    proof={
        'coefficient_field':'K=QQ(H,L,s0,q0,t0)',
        'raw134_rows_and_all45_q_unknowns_source_preserved':True,
        'old_frame_identity':'Old S=[P*F(C),J(C)*W0,J(C)*W1,polar(W0,d),polar(W1,d)] with W=P*N*T. Thus substituting alpha=T*a,beta=T*b gives the same old five frame divided by parameter P. Recomputing it would repeat692953 terms.',
        'actual236_by83':'Top134 [F(C),J(C)N,polar(d,N),polar(d,e_normal42)], bottom51 [0,L,0,0], bottom51 [0,0,L,0].',
        'column_order':'c; alpha0..19; beta0..19; n0..41 on actual q_nonpivot rows',
        'normal_q_rows42':pre['normal_q_rows42'],
        'q_to_historical_u_indices_zero_based':read(source_folder/'pulled134-input.json.gz')['q_to_u_indices_zero_based'],
        'exact_source_identity_checks':'10988 QQ polynomial comparisons: J(C)N vs direct alpha coefficient, polar(d,N) vs direct beta*sigma coefficient, and raw42 polar vs direct normal*sigma coefficient. Diagonal q_i² contributes twice. Both appended L copies independently compared to source.',
        'sigma_binding':'q=C+N*T*(a,b)+d*sigma, d=Q52-C. Q52 has sigma1, Q51 has sigma-1/2, and the toric triangle has sigma0. This is the original plane basis, not sigma identified with an unscaled E-pivot coordinate.',
        'same_plane_image_five_basis':'The actual five evaluations are the old canonical five functionals in a different parameter basis. Centered trace-zero triangle and sigma values0,-1/2,1 give f0=1+sigma-2sigma²+Q0, f1=a+Q1,f2=b+Q2,f3=a*sigma,f4=b*sigma. The known five-functionals map is injective; no new rank job.',
        'rank_equivalence_proof':[
            'Accepted L51x20 has rank18 and kernelT20x2, hence the two bottom constraints have rank36 and restrict alpha,beta to their two-dimensional kernels.',
            'On that kernel the top matrix is the actual47 columns [the five old plane-image columns/P, raw d-polar42].',
            'Modulo the five-dimensional plane image, raw d-polar is its actual residual d-polar. Accepted RAW134 cross126 injectivity gives injectivity42 for this restriction. The actual E-column contraction d=(H*q0*t0/144)*E0 is nonzero and parameter-only.',
            'Thus the constrained236x83 matrix has column rank36+47=83 overK. This is a corollary of accepted ranks, NOT a newly computed minor or inverse.'],
        'leftinverse_frontier':'A coefficient leftinverse of this small actual matrix would give a parameter-only42-direction inverse modulo the plane image without substituting T. No leftinverse,83Gaussian,Gram,or next pivot is computed.',
        'old_kernel_content_GCD_job':'Not launched: the constrained representation supersedes that preparatory route for this packet.',
        'not_claimed':['126 inverse','remaining84 coupled directions solved','normal42 ideal membership','saturation or residual empty','generic fibre exhaustion','full generic coefficient rank recomputed']}
    (out/'source-and-rank-equivalence-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    note=ROOT/f'research/notes/{ID}.md'
    paths=list(pre['input_sha256'])+[str(Path(__file__).relative_to(ROOT)),f'code/{Path(__file__).stem}.cpp',str(note.relative_to(ROOT))]
    report={
        'id':ID,'all_passed':True,'actual_QQ_236_by83_matrix_constructed':True,
        'all134_raw_sources_and_both_L_constraints_preserved':True,
        'no_kernel_substitution_or_692953_frame_repetition':True,
        'actual_total_coefficient_terms':49039,'actual_nonzero_entries':9431,
        'source_plus_output_support_upper':82343,
        'peak_live_expanded_terms':native['peak_live_expanded_terms'],
        'max_single_polynomial_terms':native['max_single_polynomial_terms'],
        'all_QQ_coefficient_identity_checks':10988,
        'rank83_derived_from_accepted_RAW134_cross_and_plane_ranks_only':True,
        'literal_constant_entries':len(constant),'literal_monomial_entries':len(monomial),
        'no_Gaussian_Gram_new_minor_point_prime_GB_or_exhaustion':True,
        'kernel_column_gcd_or_unit_elimination_executed':False,
        'input_sha256':{p:sha(ROOT/p) for p in paths},
        'work_input_sha256':sha(work/'source-input.txt'),
        'work_executable_sha256':sha(work/'native'),
        'evidence_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()},
        'native_seconds':native['native_seconds'],
        'finish_seconds':time.monotonic()-started,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('input_sha256','evidence_sha256')}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--work',required=True,type=Path)
    p.add_argument('--finish',action='store_true')
    a=p.parse_args();(finish if a.finish else prepare)(a.output.resolve(),a.work.resolve())
