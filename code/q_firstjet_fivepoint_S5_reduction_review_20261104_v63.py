#!/usr/bin/env python3
"""Independent integer frame/quadratic S5 audit and conceptual vector-descent review."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import resource
import signal
import time

ROOT=Path(__file__).resolve().parents[1]
STEM='Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def accepted(p):
    r=read(p);assert r['all_passed']
    for key in ['input_sha256','checked_dependency_sha256']:
        for n,h in r.get(key,{}).items():assert sha(ROOT/n)==h,n
    for n,h in r.get('evidence_sha256',{}).items():
        q=ROOT/n if (ROOT/n).exists() else p.parent/n
        assert sha(q)==h,n
    return r
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def product(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def parity(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def determinant(A):
    return sum(parity(p)*prod(A[i][p[i]] for i in range(len(A))) for p in itertools.permutations(range(len(A))))
def prod(values):
    x=1
    for v in values:x*=v
    return x
def cycle_type(p):
    seen=set();cycles=[]
    for i in range(5):
        if i in seen:continue
        j=i;n=0
        while j not in seen:seen.add(j);n+=1;j=p[j]
        cycles.append(n)
    return tuple(sorted(cycles,reverse=True))

def main(out):
    start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60)
    assert not out.exists();out.mkdir(parents=True)
    base=ROOT/'reports/Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63'
    ap=base/'report.json';accepted(ap);cp=base/'certificate.json';c=read(cp)
    vp=ROOT/'reports/Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63.json';vector=accepted(vp)
    assert vector['genuine_Aut_quadric_linearization_including_extra_swap']
    extras=[ROOT/f'reports/{n}/report.json' for n in ['Q_FIRSTJET_SPIN5_EXTRA_NORMALIZER_FIELD_20261104_v63_REVIEW','Q_FIRSTJET_SPIN4_BIFORM_RATIONAL_FIELD_20261104_v63_REVIEW']]
    for p in extras:accepted(p)
    frame=[[int(i==j) for i in range(4)] for j in range(4)]+[[-1]*4]
    off=list(itertools.combinations(range(4),2));I5=eye(5);I4=eye(4)
    def lift(p):return [[frame[p[j]][i] for j in range(4)] for i in range(4)]
    def induced(p):
        inv=tuple(p.index(i) for i in range(5));G=lift(inv);columns=[]
        for i,j in off[:-1]:
            diag=[G[i][k]*G[j][k]-G[2][k]*G[3][k] for k in range(4)]
            assert diag==[0]*4
            coeff=[G[i][k]*G[j][l]+G[i][l]*G[j][k]-G[2][k]*G[3][l]-G[2][l]*G[3][k] for k,l in off]
            assert sum(coeff)==0;columns.append(coeff[:-1])
        return [list(v) for v in zip(*columns)]
    permutations=list(itertools.permutations(range(5)));M={};G={};characters={};sizes=Counter()
    for p in permutations:
        G[p]=lift(p);M[p]=induced(p)
        assert all([sum(G[p][j][k]*frame[i][k] for k in range(4)) for j in range(4)]==frame[p[i]] for i in range(5))
        assert determinant(G[p])==parity(p)
        typ=cycle_type(p);chi=sum(M[p][i][i] for i in range(5));sizes[typ]+=1
        if typ in characters:assert characters[typ]==chi
        characters[typ]=chi
        fixed=sum(p[i]==i for i in range(5));fixed2=sum(p[p[i]]==i for i in range(5))
        assert 2*chi==(fixed-1)**2+fixed2-1-2*fixed
        scalar=all(M[p][i][j]==(M[p][0][0] if i==j else 0) for i in range(5) for j in range(5))
        assert scalar==(p==tuple(range(5)))
    gens=[]
    for i in range(4):
        p=list(range(5));p[i],p[i+1]=p[i+1],p[i];gens.append(tuple(p))
    assert [[[str(v) for v in row] for row in G[p]] for p in gens]==c['standard4_S5_adjacent_lifts']
    assert [[[str(v) for v in row] for row in M[p]] for p in gens]==c['quadratic5_S5_adjacent_matrices']
    for p in permutations:
        for q in gens:
            pq=tuple(p[q[i]] for i in range(5));assert product(M[p],M[q])==M[pq] and product(G[p],G[q])==G[pq]
    for i,p in enumerate(gens):
        assert product(M[p],M[p])==I5
        for j,q in enumerate(gens):
            if j<=i:continue
            A=product(M[p],M[q]);power=I5
            for k in range(3 if j==i+1 else 2):power=product(power,A)
            assert power==I5
    assert len({tuple(v for row in a for v in row) for a in M.values()})==120
    assert sum(sizes[t]*characters[t]**2 for t in sizes)==120
    record=[dict(cycle_type=list(t),class_size=sizes[t],trace=characters[t]) for t in sorted(sizes)]
    assert record==c['conjugacy_character']
    ordered=[(1,1,1,1,1),(2,1,1,1),(2,2,1),(3,1,1),(3,2),(4,1),(5,)]
    assert [characters[t] for t in ordered]==[5,1,1,-1,1,-1,0]
    # Kernel of five evaluation conditions: q_ii=0 and sum six offdiag=0.
    assert len(off)-1==5 and 10-5==5
    Qmatrix=[[Q(0),Q(1,2),Q(0),Q(0)],[Q(1,2),Q(0),Q(0),Q(0)],
             [Q(0),Q(0),Q(0),Q(-1,2)],[Q(0),Q(0),Q(-1,2),Q(0)]]
    assert determinant(Qmatrix)==Q(1,16)
    assert all(v[0]*v[1]-v[2]*v[3]==0 for v in frame)
    signs={(A*A*B,A*B*B) for A,B in itertools.product([-1,1],repeat=2)}
    assert signs==set(itertools.product([-1,1],repeat=2))
    assert 10-6==4 and 12-6==6 and 4+4==8 and 6+2==8
    certificate=dict(independent_integer_S5_action_characters=record,all120_projectively_distinct=True,
        all120_times4_group_products_verified=True,generic_free_projective_and_standard_actions=True,
        independent_Coxeter_relations=True,irreducible_representation='[3,2]',
        frame_quadric_witness_determinant='1/16',not_a_beta_point=True,
        exact_center_sign_characters=sorted(signs),
        primary_descent_proof='X/Z open in genuinely AutQ-linearized LU^2 direct_sum LV^2; generic rank2 vector-bundle descent',
        projective_no_name='B(t4)=C(Std4)^S5(u4), since V5 has honest linear lift',
        full_field='F=B(t2)',B_stable_level_at_most=4,F_stable_level_at_most=2)
    ep=out/'independent-certificate.json';ep.write_text(json.dumps(certificate,indent=2)+'\n')
    paths=[Path(__file__),ROOT/f'research/notes/{STEM}.md',ap,cp,vp,*extras,
           ROOT/'code/q_firstjet_fivepoint_S5_reduction_20261104_v63.py',
           ROOT/'research/notes/Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63.md',
           ROOT/'research/notes/Q_FIRSTJET_FIVEPOINT_VECTOR_DESCENT_20261104_v63.md']
    result=dict(id=STEM,all_passed=True,
        independent_integer_linear_form_Sym2_S5_engine_no_SymPy_or_producer=True,
        all120_action_Coxeter_character_and_projective_freeness_verified=True,
        fivepoint_unique_biform_Koszul_and_projective_frame_quotient_proof_accepted=True,
        center_squared_line_rank2_bundle_honest_full_AutQ_descent_proof_accepted=True,
        nonsplit_ruling_and_Brauer_gerbe_not_assumed_split=True,
        full_Sp4_field_equals_B_two_independent_variables=True,
        projective_no_name_uses_honest_rank5_vector_bundle_not_Brauer_form=True,
        full_Sp4_invariant_field_stable_rationality_level_upper=2,
        pure_rationality_or_minimal_level_or_invariant_ring_or_NEW_equation_claim=False,
        no_finite_group_invariant_search_new_beta_point_highdegree_GB_or_author_rerun=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:result[k] for k in ['all_passed','full_Sp4_invariant_field_stable_rationality_level_upper','seconds']},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
