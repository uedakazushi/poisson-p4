#!/usr/bin/env python3
"""Independent full character, anti-action, actual jets and degree<=14 audit.

The character upper bounds use Newton/Adams recursion and the positive-root
denominator, independently of the producers' coin DP/signed permutations.
"""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61, X, pushforward
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import J, COORDS, sp_basis

ROOT = Path(__file__).resolve().parents[1]
STEMS = ['Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63',
         'Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63',
         'Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63']


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path): return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes())
def matrix(rows): return s.Matrix([[s.Rational(x) for x in row] for row in rows])


def own_action(B, A):
    """Three covariant slots of a linear2-form, constructed entry-by-entry."""
    lookup = {coord:i for i,coord in enumerate(COORDS)}
    def value(i,j,k,col):
        if i == j: return s.S.Zero
        return (1 if i < j else -1) * B[lookup[((min(i,j),max(i,j)),k)],col]
    return s.Matrix(24,16,lambda row,col:sum(
        A[r,COORDS[row][0][0]]*value(r,COORDS[row][0][1],COORDS[row][1],col)
        +A[r,COORDS[row][0][1]]*value(COORDS[row][0][0],r,COORDS[row][1],col)
        +A[r,COORDS[row][1]]*value(*COORDS[row][0],r,col) for r in range(4)))


def character_review(weights):
    multiplicity = Counter(weights)
    # Newton identity d*h_d=sum_{k=1}^d p_k*h_(d-k), p_k=Sum e^(k*lambda).
    chars = [Counter({(0,0):1})]
    for d in range(1,15):
        total = Counter()
        for k in range(1,d+1):
            for (a,b),m in multiplicity.items():
                for (u,v),c in chars[d-k].items(): total[u+k*a,v+k*b] += m*c
        assert all(c % d == 0 for c in total.values())
        chars.append(Counter({weight:c//d for weight,c in total.items() if c}))
        assert sum(chars[-1].values()) == int(s.binomial(15+d,d))
    roots = [(2,0),(0,2),(1,1),(1,-1)]
    denominator = Counter({(0,0):1})
    for a,b in roots:
        updated = denominator.copy()
        for (u,v),c in denominator.items(): updated[u+a,v+b] -= c
        denominator = Counter({weight:c for weight,c in updated.items() if c})
    dims = [sum(coefficient*char[weight] for weight,coefficient in denominator.items()) for char in chars]
    assert dims == [1,0,0,0,1,0,0,0,5,0,0,0,9,0,3]
    return chars,dims


def own_dag(payload, point):
    result = []
    for node in payload['nodes']:
        op = node[0]
        if op == 'q': v = Fraction(node[1])
        elif op == 'v': v = Fraction(point[payload['input_names'][node[1]]])
        elif op == '+': v = sum((result[i] for i in node[1]),Fraction(0))
        elif op == '*':
            v = Fraction(1)
            for i in node[1]: v *= result[i]
        elif op == 'inv': v = 1/result[node[1]]
        else: raise AssertionError(op)
        result.append(v)
    return result


def necklaces(length):
    allwords = set(itertools.product(range(2),repeat=length))
    return sorted({min(word[i:]+word[:i] for i in range(length)) for word in allwords})


def main(output, replay):
    started = time.monotonic()
    records = []
    for stem,short in zip(STEMS,['baseline','degree8','degree12_14']):
        directory = ROOT/'reports'/stem
        report = load(directory/'report.json')
        assert report['all_passed'] is True
        for path,digest in report['input_sha256'].items(): assert sha(ROOT/path)==digest,path
        for path,digest in report['evidence_sha256'].items():
            assert sha(directory/path)==digest,path
            assert (directory/path).read_bytes()==(replay/short/path).read_bytes()
        fresh = load(replay/short/'report.json')
        assert {k:v for k,v in report.items() if k!='seconds'}=={k:v for k,v in fresh.items() if k!='seconds'}
        records.append(report)
    low = load(ROOT/'reports'/STEMS[0]/'invariant-certificates.json')
    d8 = load(ROOT/'reports'/STEMS[1]/'certificate.json')
    d14 = load(ROOT/'reports'/STEMS[2]/'certificate.json')
    B = matrix(low['primitive24_by16_basis'])
    # Independently derive all8 closure/primitive linear constraints.
    triples = list(itertools.combinations(range(4),3))
    K = s.zeros(8,24)
    for col,((i,j),k) in enumerate(COORDS):
        if k not in (i,j):
            order = (k,i,j)
            K[triples.index(tuple(sorted(order))),col]=(-1)**sum(order[a]>order[b] for a in range(3) for b in range(a+1,3))
        if (i,j) in [(0,1),(2,3)]: K[4+k,col]=1
    assert K.rank()==8 and K*B==s.zeros(8,16)
    assert s.Matrix.hstack(*K.nullspace())==B
    rows = list(B.T.rref()[1]);left=B[rows,:].inv()
    generators=sp_basis();rep=[]
    for A in generators:
        full=own_action(B,A);R=left*full[rows,:]
        assert B*R==full
        rep.append(R)
    assert rep==[matrix(R) for R in low['all_ten_representation_matrices']]
    weights=[]
    for col in range(16):
        wt=[]
        for diagonal in [(1,-1,0,0),(0,0,1,-1)]:
            action=left*own_action(B,s.diag(*diagonal))[rows,:]
            assert action.is_diagonal();wt.append(int(action[col,col]))
        weights.append(tuple(wt))
    assert weights==list(map(tuple,low['Cartan_weights']))
    chars,dims=character_review(weights)
    for data,upper in [(d8['invariant_dimensions_through10'],10),(d14['complete_character_dimensions_through14'],14)]:
        assert len(data)==upper+1
        for item in data:
            degree=item['degree'];assert item['invariant_dimension']==dims[degree]
            assert item['total_monomials']==sum(chars[degree].values())
            for term in item['Weyl_terms']:
                assert term.get('multiplicity',term.get('weight_multiplicity'))==chars[degree][tuple(term['weight'])]
    # Independently rebuild the tensor contraction symplectic form.
    Omega=s.Matrix(24,24,lambda a,b:(J[COORDS[a][0][0],COORDS[b][0][0]]*J[COORDS[a][0][1],COORDS[b][0][1]]-
        J[COORDS[a][0][0],COORDS[b][0][1]]*J[COORDS[a][0][1],COORDS[b][0][0]])*J[COORDS[a][1],COORDS[b][1]])
    O=B.T*Omega*B;assert O==matrix(low['symplectic_form16']) and O.T==-O and O.det()!=0
    assert all(R.T*O+O*R==s.zeros(16) for R in rep)
    basis=s.Matrix.hstack(*[s.Matrix(list(A)) for A in generators]);br=list(basis.T.rref()[1]);bl=basis[br,:].inv()
    for i,A in enumerate(generators):
        for j,C in enumerate(generators):
            coordinates=bl*s.Matrix(list(A*C-C*A))[br,:]
            commrep=sum((v*R for v,R in zip(coordinates,rep)),s.zeros(16))
            assert rep[i]*rep[j]-rep[j]*rep[i]==-commrep
            # Universal moment equivariance, with the ANTI-action sign.
            assert rep[i].T*O*rep[j]+O*rep[j]*rep[i]==O*commrep
    G=s.Matrix([[s.trace(A*C) for C in generators] for A in generators]);assert G.det()!=0
    Gi,Oi=G.inv(),O.inv();w4,w6=necklaces(4),necklaces(6)
    assert w4==list(map(tuple,d8['degree8_invariant_definitions']['trace_words']))
    assert w6==list(map(tuple,d14['degree12_trace_words']))
    selected=d8['independent_five_columns']
    def data(b):
        grad=[O*R*b for R in rep]
        mu=s.Matrix([(b.T*v)[0]/2 for v in grad]);co=Gi*mu
        theta=sum((c*R for c,R in zip(co,rep)),s.zeros(16));V=s.Matrix.hstack(*grad)
        H=theta+Oi*V*Gi*V.T
        return mu,co,V,[theta,H]
    def traceword(mats,word):
        value=s.eye(16)
        for letter in word: value=value*mats[letter]
        return s.trace(value)
    def evaluate(b,length):
        mu,co,V,mats=data(b);values=[traceword(mats,word) for word in (w4 if length==4 else w6)]
        if length==4: values.append((mu.T*co)[0]**2)
        return values
    def derivative_columns(b):
        # Forward first-order dual-number propagation through all matrix
        # products; no producer cyclic-gradient implementation is used.
        mu,co,V,mats=data(b);answer=s.zeros(5,16)
        for k in range(16):
            dV=s.Matrix.hstack(*[(O*R)[:,k] for R in rep]);dc=Gi*V.T[:,k]
            dt=sum((c*R for c,R in zip(dc,rep)),s.zeros(16))
            dh=dt+Oi*(dV*Gi*V.T+V*Gi*dV.T);dm=[dt,dh]
            for ii,index in enumerate(selected):
                if index==6:
                    answer[ii,k]=4*(mu.T*co)[0]*(V*co)[k]
                else:
                    product,derivative=s.eye(16),s.zeros(16)
                    for letter in w4[index]:
                        derivative=derivative*mats[letter]+product*dm[letter]
                        product=product*mats[letter]
                    answer[ii,k]=s.trace(derivative)
        return answer
    # Bind actual source jets using an independent Fraction circuit parser
    # AND direct physical homogeneous pushforward+S61 normalization.
    mappath=ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/map.json.gz'
    payload=load(mappath);family=load(ROOT/'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz')
    oldX=s.symbols('h x y z w');P=s.Matrix([[s.sympify(v).subs(dict(zip(oldX,X)),simultaneous=True) for v in row] for row in family['full_canonical_homogeneous_P']])
    model=S61();z,u,beta=primitive_beta(model);ones=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1]
    beta24=[beta[i,j].coeff(t) for i,j in itertools.combinations(range(4),2) for t in z]
    betaU=s.Matrix([[s.expand(c).coeff(v) for v in u] for c in beta24])
    E8,E12,E14=[],[],[];points=[]
    for number,sample in enumerate(d14['all_actual_new_family_samples']):
        point=sample['point'];values=own_dag(payload,point)
        U=s.Matrix([s.Rational(values[i]) for i in payload['output_u1_through_u61']])
        pivot=values[payload['normalization']['value_pivot_node']];pf=values[payload['normalization']['value_pfaffian_node']]
        assert pivot and pf and str(pivot)==sample['value_pivot'] and str(pf)==sample['value_pfaffian']
        linear=s.Matrix([[s.Rational(values[j]) for j in row] for row in payload['normalization']['linear_normalizing_matrix_nodes']])
        translation=s.eye(5)
        for i in range(1,5): translation[i,0]=1-point['base'+str(i)]
        physical=P.subs({s.Symbol(name):point[name] for name in ['n0','n1','n2']})
        normalized=model.normalize(pushforward(pushforward(physical,translation),linear))
        assert list(U)==normalized['parameters']
        b24=betaU*U[ones,:];assert list(b24)==list(map(s.Rational,sample['primitive_beta24']))
        b=left*b24[rows,:];assert B*b==b24
        row8=evaluate(b,4);row12=evaluate(b,6);E8.append(row8);E12.append(row12)
        if number<5: assert row8==list(map(s.Rational,d8['full_evaluation_matrix'][number]))
        assert row12==list(map(s.Rational,d14['degree12_evaluation'][number]))
        if number<3:
            grads=derivative_columns(b)
            assert grads==matrix(d14['degree8_gradients'][number])
            brackets=grads*Oi*grads.T
            row14=[brackets[i,j] for i,j in itertools.combinations(range(5),2)]
            assert row14==list(map(s.Rational,d14['degree14_evaluation'][number]));E14.append(row14)
        points.append(b)
    # The saved quartic is checked against every infinitesimal operator.
    quartic=low['unique_quartic_sparse']
    for R in rep:
        derivative=Counter()
        for mon,text in quartic:
            for i in set(mon):
                rem=list(mon);rem.remove(i)
                for j in range(16):
                    if R[i,j]: derivative[tuple(sorted(rem+[j]))]+=s.Rational(text)*mon.count(i)*R[i,j]
        assert all(v==0 for v in derivative.values())
    qvalue=sum(s.Rational(c)*s.prod(points[0][i] for i in mon) for mon,c in quartic)
    assert qvalue==s.Rational(low['quartic_on_new_beta'])!=0
    minors=[s.Matrix(E8[:5])[:,d8['independent_five_columns']],s.Matrix(E12)[:,d14['degree12_basis_columns']],s.Matrix(E14)[:,d14['degree14_basis_columns']]]
    expected=[d8['exact5_by5_determinant'],d14['degree12_exact9minor'],d14['degree14_first3_exact3minor']]
    determinants=[]
    for M,want in zip(minors,expected):
        det=s.Rational(DomainMatrix.from_Matrix(M).convert_to(s.QQ).det());assert det==s.Rational(want)!=0;determinants.append(str(det))
    note=ROOT/'research/notes/Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW.md'
    paths=[Path(__file__),note,mappath,ROOT/'code/weight1_structure.py',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py']
    for stem in STEMS:
        directory=ROOT/'reports'/stem;paths += [directory/'report.json',directory/('invariant-certificates.json' if stem==STEMS[0] else 'certificate.json')]
    result=dict(id='Q_NEW_FAMILY_INVARIANT_DEGREE16_20261103_v63_REVIEW',all_passed=True,
        all_three_bounded_producers_fresh_replayed_and_all_evidence_byte_equal=True,
        independent_three_covariant_slots_representation_and_all8_constraints=True,
        full_integer_characters_via_independent_Newton_Adams_recursion=True,
        independent_positive_root_denominator_extracts_all_trivial_multiplicities=True,
        invariant_dimensions_through14=dims,
        all100_anti_Lie_brackets_and_universal_moment_covariance_sign_checked=True,
        invariant_symplectic_form_and_trace_Gram_nondegenerate=True,
        theta_and_H_same_conjugation_and_constant_Poisson_bracket_invariance_proved=True,
        all9_actual_QQ_jets_direct_physical_normalization_and_own_DAG_checked=True,
        all_degree8_and12_values_and_first3_full_gradients_brackets_rebuilt=True,
        exact_nonzero_restriction_minors_degree8_12_14=determinants,
        lowdegree_quartic_all10_Lie_derivatives_zero_and_actual_nonzero_value=True,
        defining_hypersurface_degree_lower_bound=16,
        actual_defining_equation_or_general_fiber_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},seconds=time.monotonic()-started)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as stream: json.dump(result,stream,indent=2);stream.write('\n')
    print('INDEPENDENT ALL INVARIANTS THROUGH14 EXCLUDED; DEGREE>=16 PASS',result['seconds'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--replay',type=Path,required=True)
    args=ap.parse_args();main(args.output.resolve(),args.replay.resolve())
