#!/usr/bin/env python3
"""Exact ideal of the CLOSED weight6 free-coordinate slice of the 35-model.

This retains a and the four weight6 free inputs, sets all other positive free
inputs to zero, and evaluates the full finite DAG. It is not a jet truncation
and makes no assertion about all branches of the full35-variable model.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
ASTAR=s.Rational(2,9)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output,model_path):
    started=time.monotonic();raw=gzip.decompress(model_path.read_bytes());model=json.loads(raw)
    old_work=ROOT/'work/q-weighted-model-v62/model-first/model.json.gz'
    work_comparison=None
    if old_work.exists():
        old=json.loads(gzip.decompress(old_work.read_bytes()))
        assert set(model)-set(old)=={'reduced_a_open'}and not(set(old)-set(model))
        assert all(model[key]==value for key,value in old.items())
        work_comparison='All common fields identical; official payload adds reduced_a_open metadata.'
    validation_path=model_path.parent/'validation/validation.json'
    model_review_path=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63_REVIEW.json'
    resonance_path=ROOT/'reports/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62.json'
    validation,review,resonance=[json.loads(p.read_text())for p in [validation_path,model_review_path,resonance_path]]
    for source in [validation,review,resonance]:
        assert source['all_passed']
        for name,want in source.get('input_sha256',{}).items():assert sha(ROOT/name)==want,name
    assert validation['archive_sha256']==sha(model_path)
    R,a=ring('a',s.QQ);K=R.to_field();FD=K.to_domain()
    P,*ss=ring('s0,s1,s2,s3',K)
    S,w0,w1,w2,eta=ring('w0,w1,w2,eta',K);ww=[w0,w1,w2]
    zero=(0,0,0,0);u=K(a-ASTAR)
    def evaluate_poly(p):return sum(s.Rational(c)*ASTAR**m[0]for m,c in p.items())
    def specialize(c):
        den=evaluate_poly(c.denom);assert den!=0
        return evaluate_poly(c.numer)/den
    def regular(c):return evaluate_poly(c.denom)!=0
    def unit(c):return regular(c)and evaluate_poly(c.numer)!=0
    selected_inputs=[i for i,w in enumerate(model['input_weights'])if w==6]
    assert selected_inputs==[12,13,24,34]
    assert [model['input_names'][i]for i in selected_inputs]==['z11','z12','z23','z33']
    inputs=[P(K(a))]+[ss[selected_inputs.index(i)]if i in selected_inputs else P.zero
                       for i in range(1,len(model['input_names']))]
    values=[];inverse_count=0
    for op,*data in model['nodes']:
        if op=='q':value=P(K(s.Rational(data[0])))
        elif op=='v':value=inputs[data[0]]
        elif op=='+':value=sum((values[i]for i in data[0]),P.zero)
        elif op=='*':
            value=P.one
            for i in data[0]:
                if not values[i]:value=P.zero;break
                value*=values[i]
        elif op=='inv':
            previous=values[data[0]];assert previous and all(not any(m)for m in previous)
            coefficient=previous[zero];assert unit(coefficient)
            value=P(coefficient**(-1));inverse_count+=1
        else:raise ValueError(op)
        values.append(value)
    original=[values[i]for i in model['all350_residual_nodes']];assert len(original)==350
    counts=Counter(sum(next(iter(p)))for p in original if p)
    assert counts=={1:2,2:31,3:2}
    assert all(len({sum(m)for m in p})<=1 for p in original)
    assert max(sum(m)for p in original for m in p)==3

    D=27*a*a-60*a+8;H=36*a*a+13*a-8
    ell_coefficients=[162*a*a-72*a-64,54*a*a+78*a-16,27*a*a+39*a-8,288*a+128]
    ell=sum((P(K(c))*x for c,x in zip(ell_coefficients,ss)),P.zero)
    assert original[296]==-(9*K(a)-2)*ell/(12*K(D))and original[334]==-original[296]
    tau=-s.Rational(3,10)*H;assert evaluate_poly(tau)==1
    # Coordinate order is (w0,w1,w2,eta). eta is eps on the known Rstar curve.
    T=DomainMatrix([[K.one,K.zero,K.zero,K(-4)],
                    [K.zero,K.one,K.zero,K(-6)],
                    [K.zero,K.zero,K.one,K(4)],
                    *[[-K(c)/80 for c in ell_coefficients]]],(4,4),FD)
    assert T.det()==K(tau)
    t=(eta+sum((S(K(ell_coefficients[i])/80)*ww[i]for i in range(3)),S.zero))/K(tau)
    images=[w0+4*t,w1+6*t,w2-4*t,t]
    def transform(p):
        result=S.zero
        for mon,c in p.items():
            term=S(c)
            for degree,image in zip(mon,images):
                if degree:term*=image**degree
            result+=term
        return result
    assert transform(ss[0]-4*ss[3])==w0
    assert transform(ss[1]-6*ss[3])==w1
    assert transform(ss[2]+4*ss[3])==w2
    assert transform(-ell/80)==eta
    residual=[transform(p)for p in original]
    assert residual[296]==S(60/K(D))*u*eta

    qmons=[(2,0,0,0),(1,1,0,0),(1,0,1,0),(0,2,0,0),(0,1,1,0),
           (0,0,2,0),(1,0,0,1),(0,1,0,1),(0,0,1,1)]
    pure_eta2=(0,0,0,2)
    quadratic_indices=[i for i,p in enumerate(residual)if p and all(sum(m)==2 for m in p)]
    coefficients=[[residual[i].get(m,K.zero)for m in qmons]for i in quadratic_indices]
    special_matrix=s.Matrix([[specialize(c)for c in row]for row in coefficients])
    selected=list(special_matrix.T.rref()[1]);assert len(selected)==9
    row_indices=[quadratic_indices[i]for i in selected]
    assert row_indices==[10,46,47,92,96,168,224,233,276]
    N=DomainMatrix([coefficients[i]for i in selected],(9,9),FD)
    determinant=N.det();assert unit(determinant)and specialize(determinant)==-534600
    inverse=N.inv();assert(N*inverse).to_dense()==DomainMatrix.eye(9,FD).to_dense()
    b=[residual[i].get(pure_eta2,K.zero)for i in row_indices]
    assert all(regular(c/u)for c in b)
    J=[u*eta]+[S.from_dict({m:K.one})for m in qmons]

    # All350 containments in the standard ideal, with explicit polynomial lifts.
    forward=[]
    for polynomial in residual:
        row=[S.zero]*10
        for mon,c in polynomial.items():
            if not any(mon[:3]):
                assert mon[3]>0;coefficient=c/u;assert regular(coefficient)
                remainder=list(mon);remainder[3]-=1
                row[0]+=S.from_dict({tuple(remainder):coefficient})
            else:
                index=next(j for j,m in enumerate(qmons,1)if all(x>=y for x,y in zip(mon,m)))
                remainder=tuple(x-y for x,y in zip(mon,qmons[index-1]))
                row[index]+=S.from_dict({remainder:c})
        assert sum((x*y for x,y in zip(row,J)),S.zero)==polynomial
        forward.append(row)
    # Ten reverse containments in the actual350-row ideal.
    reverse=[[S.zero]*350 for _ in range(10)]
    reverse[0][296]=S(K(D)/60)
    inv_rows=inverse.to_list()
    for j in range(9):
        correction=K.zero
        for k,i in enumerate(row_indices):
            coefficient=inv_rows[j][k];reverse[j+1][i]+=S(coefficient)
            correction-=coefficient*(b[k]/u)*K(D)/60
        reverse[j+1][296]+=S(correction)*eta
    for generator,row in zip(J,reverse):
        assert sum((x*y for x,y in zip(row,residual)),S.zero)==generator

    # A finite explicit localization suffices; u is never inverted.
    source_open=R.from_expr(s.sympify(model['reduced_a_open']))
    h=R.from_expr(s.sqf_part((source_open*tau*determinant.numer).as_expr(),R.symbols[0]))
    assert evaluate_poly(h)!=0
    denominator_cache=set()
    def supported_denominator(c):
        assert regular(c);key=str(c.denom)
        if key in denominator_cache:return
        rem=c.denom
        while not rem.is_ground:
            divisor=rem.gcd(h);assert not divisor.is_ground,key
            rem=rem.exquo(divisor)
        denominator_cache.add(key)
    for p in [*original,*residual,*images,*J,
              *(p for row in forward for p in row),*(p for row in reverse for p in row)]:
        for c in p.values():supported_denominator(c)
    for row in inv_rows:
        for c in row:supported_denominator(c)
    for row in T.to_list():
        for c in row:supported_denominator(c)

    # Compare the actual raw150 reconstruction with Q*+eps Rstar-75eps² T.
    ray=[4,6,-4,1]
    def on_axis(p):
        answer=[s.S.Zero]*(1+max((sum(m)for m in p),default=0))
        for mon,c in p.items():answer[sum(mon)]+=specialize(c)*s.prod(r**d for r,d in zip(ray,mon))
        while answer and answer[-1]==0:answer.pop()
        return answer
    actual_axis=[on_axis(values[i])for i in model['raw150_reconstruction_nodes']]
    X=s.symbols('x0:5');pairs=list(itertools.combinations(range(5),2))
    mons=list(itertools.combinations_with_replacement(range(5),2));basis=list(itertools.product(pairs,mons))
    Q={(1,4):ASTAR*X[0]**2,(2,3):-s.Rational(2,9)*X[0]**2,
       (2,4):s.Rational(4,3)*X[0]*X[1],(3,4):s.Rational(2,3)*X[0]*X[2]+X[1]**2}
    Rstar={tuple(map(int,key.split(','))):s.sympify(value)for key,value in resonance['integrable_R'].items()}
    expected_axis=[]
    for pair,mon in basis:
        monomial=s.prod(X[j]for j in mon)
        coeff=[s.Poly(Q.get(pair,0),*X).coeff_monomial(monomial),
               s.Poly(Rstar.get(pair,0),*X).coeff_monomial(monomial),
               s.Integer(-75)if pair==(0,1)and mon==(4,4)else s.S.Zero]
        while coeff and coeff[-1]==0:coeff.pop()
        expected_axis.append(coeff)
    assert actual_axis==expected_axis
    reference=validation['resonance_curve']['all150_reconstruction_coefficients']
    assert actual_axis==[[s.Rational(v)for v in row]for row in reference]
    assert all(on_axis(p)==[]for p in original)

    def scalar(c):
        return [[[m[0],str(co)]for m,co in sorted(c.numer.items())],
                [[m[0],str(co)]for m,co in sorted(c.denom.items())]]
    def poly(p):return [[list(m),scalar(c)]for m,c in sorted(p.items())]
    def sparse_matrix(matrix):return [[[j,poly(p)]for j,p in enumerate(row)if p]for row in matrix]
    certificate=dict(schema='Sparse polynomials in(w0,w1,w2,eta); coefficients are(num,den), each sparse in a.',
                     variables=['w0','w1','w2','eta'],parameter='a',center='2/9',
                     coordinate_matrix=[[scalar(c)for c in row]for row in T.to_list()],
                     inverse_images_of_s0_through_s3=list(map(poly,images)),
                     all350_residuals_after_coordinate_change=list(map(poly,residual)),
                     ten_standard_generators=list(map(poly,J)),
                     forward350_by10=sparse_matrix(forward),reverse10_by350=sparse_matrix(reverse),
                     selected_quadratic_rows=row_indices,
                     quadratic_matrix=[[scalar(c)for c in row]for row in N.to_list()],
                     quadratic_matrix_inverse=[[scalar(c)for c in row]for row in inv_rows],
                     actual_raw150_resonance_curve=[[str(c)for c in row]for row in actual_axis])
    output.mkdir(parents=True,exist_ok=False)
    encoded=(json.dumps(certificate,sort_keys=True,separators=(',',':'))+'\n').encode()
    archive=output/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write(encoded)
    factors=s.factor_list(h.as_expr(),R.symbols[0])[1]
    factor_values=[dict(polynomial=str(f),value_at_center=str(f.subs(R.symbols[0],ASTAR)))for f,power in factors]
    note=ROOT/'research/notes/Q_WEIGHT6_CLOSED_SUBSLICE_20261103_v63.md'
    inputs_paths=[Path(__file__),note,model_path,model_review_path,validation_path,resonance_path]
    report=dict(id='Q_WEIGHT6_CLOSED_SUBSLICE_20261103_V63',all_passed=True,
                scope='EXACT CLOSED slice of the finite35-model, retaining a and only four weight6 free inputs.',
                coefficient_ring='QQ[a,1/h(a)]',center_a='2/9',
                comparison_to_requested_work_payload=work_comparison,
                retained_model_inputs=[model['input_names'][i]for i in selected_inputs],
                closed_coordinate_names=['s0','s1','s2','s3'],retained_input_indices=selected_inputs,
                all_other_positive_free_inputs_set_to_zero=True,full_DAG_evaluated=True,no_jet_truncation=True,
                all350_equations_retained=True,nonzero_equations_by_ordinary_degree={str(k):v for k,v in sorted(counts.items())},
                inverse_nodes_evaluated=inverse_count,all_inverse_nodes_ground_units_at_center=True,
                primitive_ell_coefficients=list(map(str,ell_coefficients)),eta='-ell/80',
                w=['s0-4*s3','s1-6*s3','s2+4*s3'],
                coordinate_matrix=[[str(c)for c in row]for row in T.to_list()],
                coordinate_determinant=str(s.factor(tau.as_expr())),coordinate_determinant_at_center='1',
                exact_ideal='((a-2/9)*eta,(w0,w1,w2)*(w0,w1,w2,eta))',
                linear_generator_source_row=296,linear_generator_multiplier='(27*a^2-60*a+8)/60',
                quadratic_generator_order=[str(p)for p in J[1:]],selected_quadratic_rows=row_indices,
                quadratic_unit_determinant=str(s.factor(determinant.as_expr())),
                quadratic_unit_determinant_at_center=str(specialize(determinant)),
                forward_actual350_into_standard10_verified=True,reverse_standard10_into_actual350_verified=True,
                explicit_localization_unit=str(s.factor(h.as_expr())),localization_factors_and_values=factor_values,
                all_certificate_denominators_supported_on_h=True,distinct_denominators_checked=len(denominator_cache),
                eta_axis_at_center_inputs=[4,6,-4,1],
                actual_raw150_axis_equals_Qstar_eps_Rstar_minus75_eps2_T=True,
                generic_a_fiber_length=4,special_a_fiber_dimension=1,
                special_a_fiber_structure='QQ[eta] plus three square-zero w directions annihilated by eta',
                reduced_closed_slice='(w0,w1,w2,(a-2/9)*eta)',
                completed_closed_slice_primary_decomposition='(eta,(w0,w1,w2)^2) intersect (a-2/9,w0,w1,w2)',
                primary_generic_lengths=[4,1],nilradical_square_zero=True,
                nilradical_free_rank_over_QQ_u=3,regular_parameter='u+eta',
                completed_closed_slice_Cohen_Macaulay=True,Hilbert_Samuel_multiplicity=5,
                completed_closed_slice_embedded_primes=False,
                full35_branch_or_component_classification_claimed=False,new_global_component_claimed=False,
                evidence_sha256={'certificate.json.gz':sha(archive)},decoded_certificate_sha256=hashlib.sha256(encoded).hexdigest(),
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs_paths},
                seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact closed weight6 slice:350 forward,10 reverse,raw150 Rstar;',report['seconds'],'seconds',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve(),args.model.resolve())
