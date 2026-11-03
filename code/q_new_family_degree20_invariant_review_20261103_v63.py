#!/usr/bin/env python3
"""Independent fresh QQ DAG points, integer trace-word minor and Weyl20 audit."""
import argparse,gzip,hashlib,itertools,json,math,resource,signal,sys,time
from collections import Counter,defaultdict
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63_REVIEW'
def read(p):
    with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for q in iter(lambda:f.read(1<<20),b''):h.update(q)
    return h.hexdigest()
def fm(M):return[[F(str(q))for q in row]for row in M.tolist()]
def mm(A,B):
    cols=list(zip(*B))
    return[[sum((a*b for a,b in zip(row,col)if a and b),0)for col in cols]for row in A]
def mv(A,v):return[sum((q*c for q,c in zip(row,v)if q and c),F(0))for row in A]
def bareiss(A):
    A=[list(row)for row in A];n=len(A);previous=1;sign=1
    for k in range(n-1):
        if not A[k][k]:
            pivot=next(i for i in range(k+1,n)if A[i][k]);A[k],A[pivot]=A[pivot],A[k];sign=-sign
        p=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value=p*A[i][j]-A[i][k]*A[k][j];q,r=divmod(value,previous);assert r==0;A[i][j]=q
            A[i][k]=0
        previous=p
    return sign*A[-1][-1]
def grouped_character(weights,maximum=20):
    chars=[defaultdict(int)for _ in range(maximum+1)];chars[0][0,0]=1
    for(a,b),multiplicity in sorted(Counter(weights).items()):
        updated=[defaultdict(int)for _ in range(maximum+1)]
        for d,old in enumerate(chars):
            for k in range(maximum-d+1):
                factor=math.comb(multiplicity+k-1,k);target=updated[d+k]
                for(u,v),n in old.items():target[u+k*a,v+k*b]+=factor*n
        chars=updated
    denominator=Counter({(0,0):1})
    for a,b in[(2,0),(0,2),(1,1),(1,-1)]:
        new=denominator.copy()
        for(u,v),n in denominator.items():new[u+a,v+b]-=n
        denominator=Counter({w:n for w,n in new.items()if n})
    for d,char in enumerate(chars):assert sum(char.values())==math.comb(15+d,d)
    return chars,[sum(n*char[w]for w,n in denominator.items())for char in chars]
def main(out):
    start=time.monotonic();signal.alarm(120);resource.setrlimit(resource.RLIMIT_AS,(2200000000,2200000000));sys.set_int_max_str_digits(0);out.mkdir(parents=True,exist_ok=False)
    parent=ROOT/'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_EXCLUSION_20261103_v63';author=read(parent/'report.json');assert author['all_passed']
    for key in['input_sha256','evidence_sha256']:
        for p,h in author[key].items():assert sha(ROOT/p if key=='input_sha256'else parent/p)==h,p
    samples=read(parent/'samples.json');exact=read(parent/'exact-evaluation.json');probe_path=ROOT/'reports/Q_NEW_FAMILY_DEGREE20_INVARIANT_PROBE_20261103_v63/report.json';probe=read(probe_path)
    info_path=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json';info=read(info_path)
    weights=list(map(tuple,info['Cartan_weights']));chars,dims=grouped_character(weights);assert dims[20]==36 and dims[19]==0
    assert all((a+b)%2 for a,b in weights),'central minus-I acts by minus one'
    for record in probe['character_records']:
        d=record['degree'];assert dims[d]==record['invariant_dimension']
        for term in record['Weyl_terms']:assert chars[d][tuple(term['weight'])]==term['multiplicity']
    from s61_model import S61
    from primitive_jet_elliptic_curve_20261101 import primitive_beta
    from weight1_structure import sp_basis
    model=S61();z,u,beta=primitive_beta(model);ones=[i for i,w in enumerate(model.reduced['parameter_weights'])if w==1]
    coefficients=[beta[i,j].coeff(t)for i,j in itertools.combinations(range(4),2)for t in z]
    betaU=s.Matrix([[s.expand(q).coeff(v)for v in u]for q in coefficients])
    basis=s.Matrix(info['primitive24_by16_basis']);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv()
    BU,BL,Basis=map(fm,[betaU,left,basis])
    map_path=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz';payload=read(map_path)
    betas=[];point_records=[]
    for idx,(sample,row)in enumerate(zip(samples['actual_QQ_source_points'],samples['probe_rows'])):
        point=sample['point'];assert point==probe['actual_source_points'][row]
        values=[]
        for node in payload['nodes']:
            op=node[0]
            if op=='q':v=F(node[1])
            elif op=='v':v=F(point[payload['input_names'][node[1]]])
            elif op=='+':v=sum((values[j]for j in node[1]),F(0))
            elif op=='*':v=math.prod(values[j]for j in node[1])
            elif op=='inv':assert values[node[1]];v=1/values[node[1]]
            else:raise AssertionError(op)
            values.append(v)
        lin=[values[payload['output_u1_through_u61'][i]]for i in ones];b24=mv(BU,lin);b16=mv(BL,[b24[i]for i in rows]);assert mv(Basis,b16)==b24
        assert list(map(str,b16))==sample['beta16'];betas.append(b16)
        vp=values[payload['normalization']['value_pivot_node']];pf=values[payload['normalization']['value_pfaffian_node']];assert vp and pf
        assert str(vp)==sample['value_pivot']and str(pf)==sample['value_pfaffian']
        point_records.append(dict(probe_row=row,fresh_source_point=point,beta16=list(map(str,b16)),value_pivot=str(vp),value_pfaffian=str(pf)))
    assert len(betas)==36 and len(set(samples['probe_rows']))==36
    O=s.Matrix(info['symplectic_form16']);Oi=O.inv();reps=[s.Matrix(M)for M in info['all_ten_representation_matrices']]
    gs=sp_basis();gram=s.Matrix([[s.trace(A*B)for B in gs]for A in gs]);Gi=gram.inv()
    # Bind all native scalar inputs, but never compile or execute its evaluator.
    tokens=iter((parent/'native-input.txt').read_text().split());assert[next(tokens),next(tokens)]==['36','108']
    def bind(M):
        assert(int(next(tokens)),int(next(tokens)))==M.shape
        for q in M:assert F(next(tokens))==F(str(q))
    for M in[O,Oi,Gi,*reps]:bind(M)
    words=[''.join(map(str,q['word']))for q in probe['definitions'][:108]]
    assert[next(tokens)for _ in range(108)]==words
    columns=samples['selected_columns'];selected_words=[tuple(q)for q in samples['selected_words']]
    assert[int(next(tokens))for _ in range(36)]==columns
    assert selected_words==[tuple(probe['definitions'][j]['word'])for j in columns]
    assert len(set(selected_words))==36 and all(len(q)==10 for q in selected_words)
    for b16 in betas:bind(s.Matrix(list(map(str,b16))))
    assert next(tokens,None)is None
    of,oif,gif=map(fm,[O,Oi,Gi]);rf=list(map(fm,reps));integer_minor=[];row_denominators=[];records=[]
    for idx,bb in enumerate(betas):
        vectors=[mv(of,mv(M,bb))for M in rf];V=[list(row)for row in zip(*vectors)]
        mu=[sum((q*c for q,c in zip(bb,v)),F(0))/2 for v in vectors];co=mv(gif,mu)
        theta=[[sum((c*M[i][j]for c,M in zip(co,rf)),F(0))for j in range(16)]for i in range(16)]
        rest=mm(mm(mm(oif,V),gif),list(map(list,zip(*V))));H=[[theta[i][j]+rest[i][j]for j in range(16)]for i in range(16)]
        den=math.lcm(*(q.denominator for M in[theta,H]for row in M for q in row));row_denominators.append(den**10)
        matrices=[[[int(den*q)for q in row]for row in M]for M in[theta,H]]
        cache={():[[int(i==j)for j in range(16)]for i in range(16)]}
        def word(q):
            q=tuple(q)
            if q not in cache:cache[q]=mm(word(q[:-1]),matrices[q[-1]])
            return cache[q]
        numerators=[sum(word(q)[i][i]for i in range(16))for q in selected_words]
        qq=[F(n,den**10)for n in numerators]
        assert qq==[F(exact['evaluation'][idx][j])for j in columns]
        integer_minor.append(numerators);records.append(dict(probe_row=samples['probe_rows'][idx],integer_matrix_denominator=str(den),selected_exact36=list(map(str,qq))))
    determinant=F(bareiss(integer_minor),math.prod(row_denominators));assert determinant and determinant==F(exact['exact_selected36_minor'])
    assert abs(determinant.numerator).bit_length()==8320 and determinant.denominator.bit_length()==5697
    assert determinant.numerator*pow(determinant.denominator,-1,32003)%32003==9474
    accepted=ROOT/'reports/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.json';prior=ROOT/'reports/Q_NEW_FAMILY_DEGREE18_INVARIANT_EXCLUSION_20261103_v63_REVIEW/report.json'
    assert read(accepted)['theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved']and read(prior)['all_passed']
    mapreview=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.json';assert read(mapreview)['all_passed']
    evidence=dict(fresh_QQ_points=point_records,integer_word_records=records,full_character_dimensions0through20=dims,exact36_minor=str(determinant))
    (out/'independent-integer-review.json').write_text(json.dumps(evidence,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',parent/'report.json',parent/'samples.json',parent/'exact-evaluation.json',parent/'native-input.txt',parent/'native.log',probe_path,info_path,map_path,accepted,prior,mapreview,ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    result=dict(id=STEM,all_passed=True,producer_and_native_not_rerun=True,
        all36_selected_actual_source_points_freshly_evaluated_on_exact_QQ_DAG=True,
        selected_fresh_new_points_outside_old20=sum(i>=20 for i in samples['probe_rows']),
        primitive24_and_all16_basis_coordinates_value_pivots_Pfaffians_checked=True,
        every_native_input_scalar_bound_to_fresh_matrices_words_and_betas=True,
        independent_scaled_integer_sequential_word_products_selected_QQ_values=1296,
        independent_integer_Bareiss36_minor_exactly_equals_author=True,
        full_inv20_dimension36_by_distinct_grouped_weight_binomial_and_positive_root_denominator=True,
        all_degree0through20_Weyl_terms_and_total_symmetric_dimensions_checked=True,
        central_minus_I_rules_out_every_odd_degree=True,
        complete_degree20_restriction_injective=True,defining_hypersurface_equation_degree_lower_bound=22,
        all3888_QQ_values_recomputed_not_claimed=True,
        degree22_probe_equation_identification_or_global_fiber_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={'independent-integer-review.json':sha(out/'independent-integer-review.json')},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:result[k]for k in['all_passed','selected_fresh_new_points_outside_old20','independent_scaled_integer_sequential_word_products_selected_QQ_values','defining_hypersurface_equation_degree_lower_bound','seconds']},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
