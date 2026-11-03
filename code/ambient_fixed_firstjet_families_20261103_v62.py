#!/usr/bin/env python3
"""Exact affine10-plane CA/CB families over their common six-plane of jets."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from additional_rank4_families_20261103_v62 import four,ROOT
from s61_model import S61,X

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic();h,x,y,z,w=X
    t,A,B,C,D,F=s.symbols('t A B C D F');base=[t,A,B,C,D,F];q=s.Symbol('q');k=s.symbols('k0:9')
    quadratic=A*x*z+B*x*w+C*z*z+D*z*w+F*w*w
    cubic_monomials=[s.prod((x,z,w)[i] for i in ids) for ids in itertools.combinations_with_replacement(range(3),3)
                     if ids!=(0,0,0)];assert len(cubic_monomials)==9
    potential=h*quadratic+sum(a*b for a,b in zip(k,cubic_monomials));model=S61();bodies={};vectors={}
    for name,p,gkind in [('CA',5*t/2,'A'),('CB',-5*t,'B')]:
        f=h*h+p*h*x+q*x*x;extra=0 if gkind=='A' else -s.diff(f,x)/2
        P=four([f,0,0,extra*z-s.diff(potential,w),extra*w+s.diff(potential,z),h*h if gkind=='A' else f])
        normalized=model.normalize(P);u=s.Matrix(normalized['parameters']);vectors[name]=u
        constant=u.subs(dict.fromkeys([q,*k],0));matrix=u.jacobian([q,*k])
        assert not set().union(*(v.free_symbols for v in matrix))
        assert u==constant+matrix*s.Matrix([q,*k]);assert matrix.rank()==10
        bodies[name]=dict(coordinates=list(map(str,u)),constant=list(map(str,constant)),direction_matrix=matrix.tolist(),
            perspective=list(map(str,normalized['perspective'])),direction_rank=10)
        print(name,'AFFINE10 PLANE READY',time.monotonic()-start,flush=True)
    weights=[item['weight'] for item in model.certificate['free_variable_basis']]
    indices=[i for i,v in enumerate(weights) if v==1]
    assert all(s.expand(vectors['CA'][i]-vectors['CB'][i])==0 for i in indices)
    MA=s.Matrix(bodies['CA']['direction_matrix']);MB=s.Matrix(bodies['CB']['direction_matrix'])
    assert MA[:,1:]==MB[:,1:];common=MA[:,1:];assert common.rank()==9
    union=MA.row_join(MB);assert union.rank()==11
    # The first nine columns absorb the difference of the independent cubic potentials.
    solve=MA.row_join(-MB[:,:1]);target=s.Matrix(bodies['CB']['constant'])-s.Matrix(bodies['CA']['constant'])
    solution=solve.gauss_jordan_solve(target)[0];assert solve*solution==target
    payload=dict(id='AMBIENT_FIXED_FIRSTJET_FAMILIES_20261103_v62',all_passed=True,
        base_parameter_order=list(map(str,base)),fiber_parameter_order=list(map(str,[q,*k])),
        cubic_potential_monomials=list(map(str,cubic_monomials)),families=bodies,
        common_firstjet_coordinates={f'u{i+1}':str(vectors['CA'][i]) for i in indices},
        common_direction_dimension=9,total_direction_span_dimension=11,
        intersection_dimension=9,intersection_equations=dict(
            CA_q=str(solution[0]),CA_cubic_minus_CB_cubic=list(map(str,solution[1:10])),CB_q=str(solution[10])),
        proof_scope='These are exact affine10-plane subfamilies with the same firstjet; their intersection is the displayed affine9-plane. Exhaustion of whole generic component fibers requires the line-parabolic incidence argument.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/additional_rank4_families_20261103_v62.py']},
        seconds=time.monotonic()-start)
    def convert(v):
        if isinstance(v,s.Rational):return str(v)
        raise TypeError(type(v))
    (args.output/'report.json').write_text(json.dumps(payload,indent=2,default=convert)+'\n')
    print('EXACT INTERSECTION AFFINE9 PASS',payload['intersection_equations'],payload['seconds'],flush=True)
