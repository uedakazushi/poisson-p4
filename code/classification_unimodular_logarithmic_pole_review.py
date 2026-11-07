#!/usr/bin/env python3
"""Independently reconstruct x^2 yz or x^2 q extension certificates.

Only Fraction(u) / sparse polynomial primitives from the independent L112
checker are reused.  Neither the producer nor a computer algebra system is
imported; the highest bivector is reconstructed as a wedge of vector fields.
"""
import classification_unimodular_logarithmic_review as arithmetic
from classification_unimodular_logarithmic_review import (
    F, RF, N, ZERO, PAIRS, atom, add, scale, mul, der, parse, scalar,
    matrix, transpose, rank, det, submatrix, annihilate, pf, cross,
    divergence, xycoeff, linearvec, roots_supported, combinations_with_replacement)
import hashlib, json, resource, time, sys
from pathlib import Path


def main():
    resource.setrlimit(resource.RLIMIT_CPU, (110, 111))
    resource.setrlimit(resource.RLIMIT_AS, (1024*1024**2, 1024*1024**2))
    started = time.monotonic()
    root = Path(__file__).resolve().parents[1]
    double_quadric = '--double-quadric' in sys.argv
    label = 'DOUBLE_QUADRIC' if double_quadric else 'DOUBLE_LINES'
    path = root/('reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_'+label+'_20261107.json')
    data = json.loads(path.read_text())
    records, coverage = [], {}
    expected_quadrics = {
        'x2q_rank2_binary_squarefree_family':'z*w+y*y+u*x*x',
        'x2q_rank2_binary_mixed':'z*w+x*y',
        'x2q_rank2_binary_eigen_square':'z*w+x*x',
        'x2q_rank2_binary_zero':'z*w',
        'x2q_rank1_aligned_zero':'z*z/2+x*w',
        'x2q_rank1_aligned_nonzero':'z*z/2+x*w+y*y',
        'x2q_rank1_transverse_zero':'z*z/2+y*w',
        'x2q_rank1_transverse_nonzero':'z*z/2+y*w+x*x',
        'x2q_rank0_full_mixed':'x*z+y*w'}
    if double_quadric:
        assert {c['name'] for c in data['cases']}==set(expected_quadrics)
        assert [(c['name'],c['residue']) for c in data['cases'] if c['residue']!='generic u']==[('x2q_rank2_binary_squarefree_family','0')]
    else:
        assert {c['name'] for c in data['cases']} == {'x2yz_b_nonzero', 'x2yz_b_zero'}
    for c in data['cases']:
        begin = time.monotonic()
        arithmetic.U = RF((F(0), F(1))) if c['residue']=='generic u' else RF(F(c['residue']))
        x, y, z, w = (atom(1, i) for i in range(4))
        if double_quadric:
            q = parse(expected_quadrics[c['name']]); assert q==parse(c['quadric'])
            dq = [der(q,i) for i in range(4)]
            gamma = [add(scale(q,2*int(i==0 and j==1)),
                mul(add(y,x),add(scale(dq[j],int(i==0)),scale(dq[i],-int(j==0)))),
                scale(mul(x,add(scale(dq[j],int(i==1)),scale(dq[i],-int(j==1)))),-1)) for i,j in PAIRS]
            Q = [gamma[5],scale(gamma[4],-1),gamma[3],gamma[2],scale(gamma[1],-1),gamma[0]]
        else:
            a, b = (arithmetic.U, RF(1)) if c['name']=='x2yz_b_nonzero' else (RF(1), RF(0))
            V = [x, scale(y, -2), {}, add(w, scale(x, a))]
            W = [{}, y, scale(z, -1), scale(x, b)]
            Q = [add(mul(V[i], W[j]), scale(mul(V[j], W[i]), -1)) for i,j in PAIRS]
        assert Q == list(map(parse, c['highest_bivector']))
        assert not pf(Q) and not any(divergence(Q))
        columns, mon3 = [], set()
        for i in range(6):
            for j in range(4):
                L = [atom(1,j) if k==i else {} for k in range(6)]
                cr = cross(Q,L); mon3.update(cr)
                columns.append(([v.get(ZERO,RF()) for v in divergence(L)], cr))
        mon3 = sorted(mon3, reverse=True)
        M = transpose([d+[cr.get(e,RF()) for e in mon3] for d,cr in columns])
        assert M == matrix(c['full_linear_matrix'])
        K = matrix(c['full_kernel']); r = c['rank']
        assert annihilate(M,K) and rank(K)==24-r
        pivot = det(submatrix(M,c['pivot_rows'],c['pivot_columns']))
        assert pivot and pivot==scalar(c['pivot_determinant']) and len(c['pivot_rows'])==r
        T = transpose([linearvec([der(p,j) for p in Q]) for j in range(4)])
        assert T==matrix(c['translation_basis']) and annihilate(M,T)
        tr = c['translation_rank']; assert rank(T)==tr
        trrows, trcols = c['translation_cut_rows'], []
        for j in range(4):
            if rank(submatrix(T,trrows,trcols+[j]))>len(trcols): trcols.append(j)
        trminor = det(submatrix(T,trrows,trcols))
        assert trminor and trminor==scalar(c['translation_minor'])
        Ln = list(map(parse,c['normal_linear_bivector'])); nd = c['normal_dimension']
        NK = transpose([linearvec([der(p,4+i) for p in Ln]) for i in range(nd)]) if nd else [[] for _ in range(24)]
        assert nd+tr==24-r and annihilate(M,NK) and rank(NK)==nd
        assert all(not v for i in trrows for v in NK[i])
        assert rank([v+t for v,t in zip(T,NK)])==24-r
        mon2 = []
        for i,j in combinations_with_replacement(range(4),2):
            e=[0]*N; e[i]+=1; e[j]+=1; mon2.append(tuple(e))
        G = transpose([[cross(Q,[atom() if k==i else {} for k in range(6)]).get(e,RF()) for e in mon2] for i in range(6)])
        assert G==matrix(c['constant_equation_matrix'])
        gr,gc = c['constant_pivot_rows'],c['constant_pivot_columns']
        gm = det(submatrix(G,gr,gc))
        assert gm and gm==scalar(c['constant_minor']) and rank(G)==len(gc)
        C = list(map(parse,c['constant_source']))
        free = [i for i in range(6) if i not in gc]
        assert all(C[i]==atom(1,9+j) for j,i in enumerate(free))
        Pf = pf([add(v,l,d) for v,l,d in zip(Q,Ln,C)])
        coeff = xycoeff(Pf)
        assert all(mon2[i][:4] not in coeff for i in gr)
        actual={tuple(sorted(p.items())) for e,p in coeff.items() if sum(e)>0}
        expected={tuple(sorted(parse(p).items())) for p in c['remaining_full_Pfaffian_equations']}
        assert actual==expected and coeff.get((0,)*4,{})==parse(c['constant_Pfaffian'])
        mixed = c['name']=='x2q_rank0_full_mixed'
        if mixed:
            # Reduced support: an actual residual is -2*n0^2/5. After n0=0,
            # another actual residual is -5*n1^2/2. No radical/GB assumption.
            remaining = list(map(parse,c['remaining_full_Pfaffian_equations']))
            assert parse('-2*n0**2/5') in remaining
            after_n0 = [{e:v for e,v in p.items() if not e[4]} for p in remaining]
            assert parse('-5*n1**2/2') in after_n0
            L0 = [{e:v for e,v in p.items() if not e[4] and not e[5]} for p in Ln]
            assert all(not der(p,j) for p in L0[:5] for j in (2,3))
            assert all(not der(der(p,k),j) for p in Q[:5] for k in range(4) for j in (2,3))
        else:
            assert not Pf and not parse(c['constant_Pfaffian'])
        if c['residue']=='generic u':
            allowed = {F(d['residue']) for d in data['cases'] if d['name']==c['name'] and d['residue']!='generic u'}
            found = set()
            for v in (pivot,trminor,gm): found |= roots_supported(v.p,allowed)|roots_supported(v.q,allowed)
            for row in K+NK:
                for v in row: found |= roots_supported(v.q,allowed)
            for p in C:
                for v in p.values(): found |= roots_supported(v.q,allowed)
            assert found==allowed
            coverage[c['name']] = list(map(str,sorted(found)))
        records.append({'name':c['name'],'full_linear_rank':r,'translation_rank':tr,
            'normal_dimension':nd,'constant_rank':len(gc),'pivot_determinant':c['pivot_determinant'],
            'constant_minor':c['constant_minor'],'all_original_coefficient_matrices_and_complete_kernels':True,
            'true_translation_slice_and_inverse_constant_solve':True,'all_parameter_values_covered':True,
            'entire_reduced_extension_source_has_identically_zero_Pfaffian':not mixed,
            'whole_translation_including_primitive_evaluation_rank_at_most_three':mixed,
            'all_passed':True,'seconds':time.monotonic()-begin})
        print(c['name'],'PASS',round(time.monotonic()-begin,3),flush=True)
    assert len(records)==(10 if double_quadric else 2)
    assert len(coverage)==(9 if double_quadric else 2)
    result={'task_id':'G','all_passed':True,'cases':records,'exceptional_coverage':coverage,
        'method':'Independent standard-library Fraction(u), sparse polynomial coefficient reconstruction and wedge construction; no producer or CAS imported.',
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (path,Path(arithmetic.__file__))},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'All ten x^2 q source certificates, including the unique exceptional parameter0, all residual equations and full-source primitive rank<=3 implication for rank0 mixed; geometric completeness of quadric reduction is separate.' if double_quadric else 'Both complete x^2 yz independent-numerator sources, for every parameter. Dependent linears, zero nilpotent source and E3 boundary are separate geometric proof dependencies.',
        'seconds':time.monotonic()-started}
    (root/('reports/CLASSIFICATION_UNIMODULAR_LOGARITHMIC_'+label+'_REVIEW_20261107.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'cases':len(records),'seconds':result['seconds']}))


if __name__=='__main__': main()
