#!/usr/bin/env python3
"""Independent fixed-formal source28 coefficients; never calls the native producer."""
import argparse
import gzip
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import q_new_family_v_zero_degree26_restriction_review_20261104_v63 as q

ROOT=Path(__file__).resolve().parents[1]
BASE='Q_NEW_FAMILY_V_ZERO_DEGREE28_RESTRICTION_20261104_v63'
STEM=BASE+'_REVIEW'
P=q.P
STAGE='inputs'

def decode_records(records):
    return {q.pack(e):int(c)%P for e,c in records if int(c)%P}

def constant(v=1):return {0:v%P} if v%P else {}

def power(a,n):
    out=constant()
    for _ in range(n):out=q.mul(out,a)
    return out

def product(values):
    out=constant()
    for a in values:out=q.mul(out,a)
    return out

def main(out):
    global STAGE
    started=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('independent source28 review60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    folder=ROOT/'reports'/BASE
    old=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63'
    previous=ROOT/'reports/Q_NEW_FAMILY_V_ZERO_DEGREE26_RESTRICTION_20261104_v63_REVIEW'
    universal=ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63'
    gatepaths=[previous/'report.json',ROOT/'reports/Q_NEW_FAMILY_DEGREE28_PRODUCT_BASIS_20261104_v63_REVIEW/report.json']
    primary=q.read(folder/'report.json');native=q.read(folder/'native-certificate.json');recipe=q.read(folder/'recipe.json')
    gates=[q.read(p) for p in gatepaths]
    assert primary['all_passed'] and native['all_passed'] and all(g['all_passed'] for g in gates)
    assert gates[0]['QQ_restriction_injective_and_NEW_defining_degree_lower_bound']==28
    assert gates[1]['complete_QQ_Inv28_basis']
    checked={}
    for reportpath in [folder/'report.json',*gatepaths]:
        d=q.read(reportpath)
        for key in ['input_sha256','checked_dependency_sha256']:
            for n,h in d.get(key,{}).items():
                assert q.sha(ROOT/n)==h,n;checked[n]=h
        for n,h in d.get('evidence_sha256',{}).items():
            p=reportpath.parent/n;assert q.sha(p)==h,str(p);checked[str(p.relative_to(ROOT))]=h
    old_evidence=q.read(previous/'independent-polynomials.json.gz')
    beta=[decode_records(z) for z in old_evidence['beta16']]
    low=[decode_records(z) for z in old_evidence['low_values19']]
    traces=[decode_records(z) for z in old_evidence['trace_values8']]
    grads=[decode_records(z) for z in old_evidence['all128_ambient_trace_gradients']]
    assert len(beta)==16 and sum(map(len,beta))==12 and len(low)==19
    STAGE='bind input prefix and every appended ambient coefficient'
    oldtext=gzip.decompress((old/'native-input.txt.gz').read_bytes()).decode()
    text=gzip.decompress((folder/'native-input.txt.gz').read_bytes()).decode()
    assert text.startswith(oldtext)
    tok=iter(oldtext.split())
    assert (int(next(tok)),int(next(tok)))==(16,136)
    for expected in beta:
        found={}
        for _ in range(int(next(tok))):
            c=q.mod(next(tok));e=[int(next(tok)) for _ in range(4)]
            if c:found[q.pack(e)]=c
        assert found==expected
    mon=[tuple(int(next(tok)) for _ in range(2)) for _ in range(136)]
    cc=[[int(next(tok))%P for _ in range(136)] for _ in range(32)]
    hc=[[int(next(tok))%P for _ in range(136)] for _ in range(512)]
    Oi=[q.mod(next(tok)) for _ in range(256)]
    oldselected=[int(next(tok)) for _ in range(45)];assert next(tok,None) is None
    suffix=iter(text[len(oldtext):].split())
    for count,expected in [(19,low),(128,grads),(8,traces)]:
        assert int(next(suffix))==count
        for want in expected:
            found={}
            for _ in range(int(next(suffix))):
                c=q.mod(next(suffix));e=[int(next(suffix)) for _ in range(4)]
                if c:found[q.pack(e)]=c
            assert found==want
    selected24=[int(next(suffix)) for _ in range(67)]
    selected28=[int(next(suffix)) for _ in range(108)]
    assert next(suffix,None) is None
    full24=ROOT/'reports/Q_NEW_FAMILY_FULL_DEGREE24_HESSIAN_BASIS_20261104_v63/native-certificate.json.gz'
    assert selected24==q.read(full24)['degree24_selected_columns']
    definitions=q.read(universal/'definitions.json')
    assert selected28==definitions['selected_basis_indices']==native['selected_universal_candidates']
    assert recipe['universal_basis']==definitions and recipe['original_scalar_xi0_power']==14
    assert recipe['source_polynomial_variables']==['H','L','s0','q0']
    STAGE='independent C,D pencil and ambient Poisson cross-check'
    quad=[q.mul(beta[i],beta[j]) for i,j in mon]
    cd=[q.total(q.scale(v,c*pow(4,-1,P)%P) for c,v in zip(row,quad) if c) for row in cc]
    C,D=cd[:16],cd[16:]
    moment=[q.add(q.scale(a,-1),q.scale(b,2)) for a,b in zip(C,D)]
    F4=q.tracepair(moment,moment,4);assert F4==low[9]
    pencil=[{} for _ in range(5)]
    for perm in itertools.permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        # Leibniz formula, with the determinant-variable degree indexed separately.
        for letters in itertools.product(range(2),repeat=4):
            value=product((C,D)[letters[i]][4*i+perm[i]] for i in range(4))
            k=sum(letters);pencil[k]=q.add(pencil[k],value,sign)
    assert q.mul(F4,F4)==q.add(q.scale(pencil[0],400),q.scale(pencil[4],1600))
    for a,b,expected in [(0,1,low[0]),(0,2,low[1]),(1,2,low[2])]:
        value=q.total(q.scale(q.mul(grads[16*a+i],grads[16*b+j]),Oi[16*i+j]) for i in range(16) for j in range(16) if Oi[16*i+j])
        assert value==expected
    for i in range(8):assert q.total(q.mul(beta[j],grads[16*i+j]) for j in range(16))==q.scale(traces[i],8 if i<5 else 12)
    STAGE='full16 by16 theta transport through seventh power'
    matrices=[q.total(q.scale(v,c*pow(2,-1,P)%P) for c,v in zip(row,quad) if c) for row in hc]
    theta,Hess=matrices[:256],matrices[256:]
    level=Hess
    H12=None
    for degree in range(1,8):
        level=q.mm(theta,level,16)
        if degree==5:H12=q.total(level[16*i+i] for i in range(16))
    H16=q.total(level[16*i+i] for i in range(16))
    assert H12==low[18]
    primitive12=[low[15],low[16],low[17],H12]
    q14=low[:3]
    generators=[F4,*pencil,*primitive12,H16,*q14]
    assert generators==[q.saved_poly(p) for p in q.read(folder/'source-generators14.json')]
    STAGE='semantic complete24 basis and all108 selected28 source columns'
    spectral=[F4,*pencil[1:]]
    def inv24(z):
        k=z['kind']
        if k=='spectral':return product(power(p,e) for p,e in zip(spectral,z['F4_p1_p2_p3_p4_exponents']))
        if k=='primitive12_times_spectral12':
            j=z['spectral'];v=power(F4,3) if j==0 else q.mul(F4,pencil[j])
            return q.mul(primitive12[z['primitive']],v)
        if k=='primitive12_product':return product(primitive12[i] for i in z['indices'])
        assert k=='primitive16_times_Inv8'
        return q.mul(H16,power(F4,2) if z['index']==0 else pencil[z['index']])
    basis24=[inv24(z) for z in definitions['degree24_basis_definitions']]
    columns=[];descended=[];QQterms=0
    for j,index in enumerate(selected28):
        z=definitions['candidates'][index];k=z['kind']
        if k=='F4_times_Inv24_basis':value=q.mul(F4,basis24[z['Inv24_basis_index']])
        elif k=='primitive12_times_p_p':value=product([primitive12[z['primitive12_index']],*[pencil[i] for i in z['p_indices']]])
        elif k=='primitive16_times_primitive12':value=q.mul(H16,primitive12[z['primitive12_index']])
        else:
            assert k=='Inv14_product';value=product(q14[i] for i in z['indices'])
        saved=q.read(folder/f'column{j}.json')
        assert saved['universal_candidate_index']==index and saved['xi0_power']==14
        assert value==q.saved_poly(saved['HLsq_polynomial'])
        QQterms+=len(saved['HLsq_polynomial']['terms'])
        assert len(saved['HLsq_polynomial']['terms'])<=3696
        for coefficient,a,b,c,d in saved['HLsq_polynomial']['terms']:
            q.mod(coefficient)
            assert b%2==0 and (2*a+3*b)%4==0 and b<=8 and a+b<=28 and c+d<=(28-3*b)//2
        reduced={}
        for m,c in value.items():
            a,b,ss,qq=q.unpack(m)
            assert b%2==0 and (2*a+3*b)%4==0
            e=((2*a+3*b)//4,b//2,ss,qq)
            assert e not in reduced;reduced[e]=c
        columns.append(value);descended.append(reduced)
    rows=[tuple(e) for e in native['coefficient_monomials_r0_B0_s0_q0']]
    matrix=[[col.get(e,0) for col in descended] for e in rows]
    assert matrix==native['selected_coefficient_matrix']
    assert len(rows)==len(set(rows))==native['rank']
    # Verify the actual frozen 107 minor, without choosing another point/prime/minor.
    fixed=q.read(folder/'modular-selected-minor.json')
    pivots=fixed['pivot_columns'];free=fixed['free_columns']
    assert native['rank']==fixed['rank']==107 and free==[106]
    assert pivots==list(range(106))+[107] and fixed['coefficient_rows']==[list(e) for e in rows]
    square=[[row[i] for i in pivots] for row in matrix]
    determinant=q.det(square)
    assert square==fixed['matrix'] and determinant==fixed['minor']==23998
    iv=q.inverse(square);right=[row[106] for row in matrix]
    solution=[-sum(a*b for a,b in zip(row,right))%P for row in iv]
    relation=[0]*108;relation[106]=1
    for i,c in zip(pivots,solution):relation[i]=c
    # One whole-polynomial identity certifies the upper rank too. It is ONLY mod p.
    assert q.total(q.scale(col,c) for col,c in zip(columns,relation))=={}
    assert len(columns)==108 and max(max(q.unpack(m)) for col in columns for m in col)<q.RADIX
    STAGE='save frozen independent receipt'
    out.mkdir(parents=True,exist_ok=False)
    ep=out/'independent-polynomials.json.gz'
    evidence={'prime':P,'variables':['H','L','s0','q0'],'packed_base64_no_carry_bound':28,
        'source_generators14':[q.sparse_records(p) for p in generators],
        'all108_source_columns':[q.sparse_records(p) for p in columns],
        'coefficient_monomials':rows,'independent107_by108_rows':matrix,
        'fixed107_pivot_columns':pivots,'independent107_by107_minor':square,
        'minor_mod32003':determinant,'all_coefficient_modular_kernel_vector':relation}
    with gzip.open(ep,'wt') as f:json.dump(evidence,f);f.write('\n')
    for p in [Path(__file__),ROOT/'code/q_new_family_v_zero_degree26_restriction_review_20261104_v63.py',ROOT/'research/notes'/(STEM+'.md'),ROOT/'research/notes'/(BASE+'.md'),folder/'report.json',folder/'recipe.json',folder/'native-certificate.json',previous/'independent-polynomials.json.gz',universal/'definitions.json',full24,*gatepaths]:
        checked[str(p.relative_to(ROOT))]=q.sha(p)
    report={'id':STEM,'all_passed':True,
        'accepted_actual_dense_v0_and_ambient16_derivative_binding_reused':True,
        'every_input_prefix_and_appended_low19_gradient128_trace8_bound':True,
        'independent_4by4_Leibniz_pencil_and_three_ambient_brackets':True,
        'full16by16_theta5_H12_checked_theta7_H16_reconstructed':True,
        'all14_generator_and108_source_column_coefficients_independently_verified_mod32003':True,
        'all11556_fixed_coefficient_entries_and11449_minor_entries_verified':True,
        'GF32003_exact_coefficient_rank':107,'coefficient_minor_mod32003':determinant,
        'all_coefficient_modular_kernel_identity_verified':True,
        'QQ_rank_lower_bound':107,'QQ_rank107_or_QQ_relation_proved':False,
        'NEW_defining_degree_lower_bound_unchanged':28,
        'degree28_equation_degree30_lower_bound_or_source30_claim':False,
        'saved_QQ_column_term_count':QQterms,
        'author_native_replay_newpoint_newprime_or_new_basis_search':False,
        'input_sha256':checked,'evidence_sha256':{'independent-polynomials.json.gz':q.sha(ep)},
        'seconds':time.monotonic()-started,'polynomial_products':q.PRODUCTS,
        'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','coefficient_minor_mod32003','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:main(a.output)
    except (TimeoutError,MemoryError) as e:
        a.output.mkdir(parents=True,exist_ok=True)
        (a.output/'blocked-review.json').write_text(json.dumps({'all_passed':False,'stage':STAGE,'reason':str(e),'no_automatic_larger_retry':True})+'\n')
        raise
