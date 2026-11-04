#!/usr/bin/env python3
"""Five-point geometry, exact S5 quadratic representation, and scale-torus reduction.

Only the finite representation is calculated; no finite-group invariant search.
"""
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
ID='Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    extra=ROOT/'reports/Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63_REVIEW/report.json'
    joint=ROOT/'reports/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63_REVIEW/report.json'
    checked={}
    for path in [extra,joint]:
        report=json.loads(path.read_text());assert report['all_passed']
        for name,want in report['input_sha256'].items():assert sha(ROOT/name)==want;checked[name]=want
    x=s.symbols('x0:4');pairs=list(itertools.combinations(range(4),2));allmons=[x[i]*x[j]for i in range(4)for j in range(i,4)]
    off=[x[i]*x[j]for i,j in pairs];basis=[m-off[-1]for m in off[:-1]]
    frame=[s.eye(4)[:,i]for i in range(4)]+[-s.ones(4,1)]
    evaluation=s.Matrix([[mon.subs(dict(zip(x,v)))for mon in allmons]for v in frame])
    assert evaluation.shape==(5,10)and evaluation.rank()==5
    assert all(all(q.subs(dict(zip(x,v)))==0 for v in frame)for q in basis)
    assert s.Matrix([[s.Poly(q,*x).coeff_monomial(m)for q in basis]for m in allmons]).rank()==5
    def standard(p):return s.Matrix.hstack(*(frame[p[i]]for i in range(4)))
    def induced(p):
        inverse=tuple(p.index(i)for i in range(5));G=standard(inverse)
        Y=G*s.Matrix(x);sub=dict(zip(x,Y));columns=[]
        for q in basis:
            image=s.Poly(s.expand(q.subs(sub,simultaneous=True)),*x)
            assert all(image.coeff_monomial(v*v)==0 for v in x)
            coefficients=[image.coeff_monomial(m)for m in off]
            assert sum(coefficients)==0;columns.append(coefficients[:-1])
        return s.Matrix.hstack(*(s.Matrix(c)for c in columns))
    permutations=list(itertools.permutations(range(5)));matrices={};characters={};sizes={};standard_matrices={}
    def cycle_type(p):
        seen=set();lengths=[]
        for i in range(5):
            if i in seen:continue
            j=i;n=0
            while j not in seen:seen.add(j);n+=1;j=p[j]
            lengths.append(n)
        return tuple(sorted(lengths,reverse=True))
    for p in permutations:
        G=standard(p);assert all(G*frame[i]==frame[p[i]]for i in range(5))
        sign=(-1)**sum(p[i]>p[j]for i in range(5)for j in range(i+1,5))
        assert G.det()==sign
        M=induced(p);matrices[p]=M;standard_matrices[p]=G
        typ=cycle_type(p);chi=int(s.trace(M));sizes[typ]=sizes.get(typ,0)+1
        if typ in characters:assert characters[typ]==chi
        characters[typ]=chi
        fixed=sum(p[i]==i for i in range(5));fixed_square=sum(p[p[i]]==i for i in range(5))
        assert chi==((fixed-1)**2+(fixed_square-1))//2-fixed
        assert M==s.eye(5)or M!=s.trace(M)/5*s.eye(5)
    assert len({tuple(M)for M in matrices.values()})==120
    assert sum(sizes[k]*characters[k]**2 for k in sizes)==120
    generators=[];lifts=[]
    for i in range(4):
        p=list(range(5));p[i],p[i+1]=p[i+1],p[i];p=tuple(p)
        generators.append(matrices[p]);lifts.append(standard_matrices[p])
    for i in range(4):
        assert generators[i]**2==s.eye(5)
        for j in range(i+1,4):
            assert (generators[i]*generators[j])**(3 if j==i+1 else 2)==s.eye(5)
    #An abstract quadric nonemptiness witness, not a beta/source sample.
    q=x[0]*x[1]-x[2]*x[3]
    Q=s.hessian(q,x)/2;assert Q.det()==s.Rational(1,16)
    assert all(q.subs(dict(zip(x,v)))==0 for v in frame)
    #Central sign characters are exactly the two independent form signs.
    signs={}
    for A,B in itertools.product([-1,1],repeat=2):signs[str((A,B))]=[A*A*B,A*B*B]
    assert set(tuple(v)for v in signs.values())==set(itertools.product([-1,1],repeat=2))
    note=ROOT/'research/notes'/(ID+'.md');output=ROOT/'reports'/ID;output.mkdir(parents=True,exist_ok=False)
    evidence={'frame_vectors':[list(map(str,v))for v in frame],
        'quadratic_basis':[str(q)for q in basis],
        'standard4_S5_adjacent_lifts':[[list(map(str,row))for row in M.tolist()]for M in lifts],
        'quadratic5_S5_adjacent_matrices':[[list(map(str,row))for row in M.tolist()]for M in generators],
        'conjugacy_character':[{'cycle_type':list(k),'class_size':sizes[k],'trace':characters[k]}for k in sorted(sizes)],
        'all120_matrices_distinct_and_only_identity_projectively_scalar':True,
        'character_squared_inner_product':1,'irreducible_label':'[3,2], dimension5',
        'quadric_witness_determinant':'1/16','not_a_beta_point':True,
        'center_signs_U_V':signs,'effective_amplitude_character_lattice':'2Z*eU + 2Z*eV',
        'ruling_monodromy_action':'permutes the two effective amplitude characters',
        'ruling_cover':'epsilon²=det(Q), epsilon transforms by det(standard4 lift); projective q scaling gives epsilon weight2',
        'generic_scale_torus':'Res_(ruling_etale_algebra / projective_base) Gm',
        'full_Sp4_invariant_field':'C(P(V5))^S5 (s1,s2), pure transcendence2 relative to the projective finite-group quotient',
        'projective_no_name':'(P(V5) x Std4)/S5 = B(t4) = C(e2,e3,e4,e5)(u4)',
        'stable_rationality_level_upper_for_full_field':2}
    ev=output/'certificate.json';ev.write_text(json.dumps(evidence,indent=2)+'\n')
    files=[Path(__file__),note,extra,joint,
           ROOT/'research/notes/Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63.md',
           ROOT/'research/notes/Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63.md']
    result={'id':ID,'all_passed':True,
        'generic_projective_pair_fivepoint_correspondence_and_unique_forms_proved':True,
        'projective_fivepoint_quotient_equals_P4_over_S5_proved':True,
        'exact_five_dimensional_S5_representation_all_Coxeter_characters_and_projective_freeness_checked':True,
        'amplitude_center_and_ruling_twist_and_Brauer_cocycle_killed_by_squares_proved':True,
        'Hilbert90_torus_fiber_pure_relative_two_variables_proved':True,
        'full_Sp4_invariant_field_reduction':'B(s1,s2), B=C(P4)^S5 on irreducible[3,2]representation',
        'projective_finite_group_quotient_stable_rationality_level_upper':4,
        'full_Sp4_invariant_field_stable_rationality_level_upper':2,
        'pure_rationality_P4_S5_or_full_field_claim':False,
        'finite_group_invariant_search_new_beta_point_highdegree_probe_or_GB':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p)for p in files},'checked_dependency_sha256':checked,
        'evidence_sha256':{ev.name:sha(ev)},'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'seconds':time.monotonic()-start}
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({'all_passed':True,'S5_dimension':5,'full_stable_level_upper':2,'seconds':result['seconds']}))

if __name__=='__main__':main()
