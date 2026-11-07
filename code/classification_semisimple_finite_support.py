#!/usr/bin/env python3
"""Finite support charts for complete affine semisimple modular sources.

The SS12 mode retains all nine affine coefficients and all seven Jacobi
conditions.  It classifies this fixed spectrum, not arbitrary spectral arcs.
"""
import argparse
import hashlib
import itertools as it
import json
import resource
import time
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ss12(output):
    begun = time.monotonic()
    source_path = ROOT / 'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_SEMISIMPLE_SOURCES_20261107.json'
    row = json.loads(source_path.read_text())['full_semisimple_sources'][11]
    assert row['eigenvalues'] == [-2, -1, 1, 2]
    n = s.symbols('n0:9')
    n0, n1, n2, n3, n4, n5, n6, n7, n8 = n
    coords = s.symbols('x y z w')
    x, y, z, w = coords
    env = {str(a): a for a in n + coords}
    brackets = [s.sympify(p, locals=env) for p in row['brackets']]
    expected = [(n5+n7-1)*x*y, n0*y+(-n5+n8+1)*x*z,
                n1+n2*y*z+(2-n7-n8)*w*x,
                n3+n4*w*x+n5*y*z, n6*z+n7*w*y, n8*w*z]
    assert all(s.expand(a-b) == 0 for a, b in zip(brackets, expected))
    P = s.zeros(4)
    for (i, j), p in zip(it.combinations(range(4), 2), brackets):
        P[i, j], P[j, i] = p, -p
    jac = [s.expand(sum(P[i,l]*s.diff(P[j,k], coords[l])+
                       P[j,l]*s.diff(P[k,i], coords[l])+
                       P[k,l]*s.diff(P[i,j], coords[l]) for l in range(4)))
           for i, j, k in it.combinations(range(4), 3)]
    assert all(s.expand(f-s.sympify(g, locals=env)) == 0
               for f, g in zip(jac, row['all_four_jacobi']))
    generators = [n4*(n7+n8-1), n2*n4, n1*n4-n3*(n7+n8),
                  n0*(n7-1), n6*(n5-n8), n2*(2*n5-1),
                  n0*n6+n1*(1-n5)-n2*n3]
    # Exact constant-span equality proves that these are ALL Jacobi equations.
    coeffs = [c for f in jac for c in s.Poly(f, *coords).coeffs()]
    mons = sorted(set(m for f in coeffs+generators for m, c in s.Poly(f, *n).terms()))
    def vector(f):
        p = s.Poly(f, *n)
        return s.Matrix([[p.coeff_monomial(m) for m in mons]])
    C = s.Matrix.vstack(*(vector(f) for f in coeffs))
    G = s.Matrix.vstack(*(vector(f) for f in generators))
    assert C.rank() == G.rank() == C.col_join(G).rank() == 7
    forward = [list(s.linsolve((G.T, vector(f).T)))[0] for f in coeffs]
    reverse = [list(s.linsolve((C.T, vector(f).T)))[0] for f in generators]
    reverse = [tuple(v.subs({t: 0 for a in v for t in a.free_symbols}) for v in reverse)]
    reverse = reverse[0]
    assert all(s.expand(f-sum(a*b for a,b in zip(v,generators))) == 0
               for f,v in zip(coeffs,forward))
    assert all(s.expand(f-sum(a*b for a,b in zip(v,coeffs))) == 0
               for f,v in zip(generators,reverse))

    charts = []
    def add(name, mapping, nonzero, zero):
        images = [s.sympify(mapping.get(v, v)) for v in n]
        free = sorted(set().union(*(p.free_symbols for p in images)), key=str)
        assert len(free) <= 3
        rename = dict(zip(free, s.symbols('t0:'+str(len(free)))))
        images = [s.cancel(p.subs(rename, simultaneous=True)) for p in images]
        sub = dict(zip(n, images))
        def image(f):
            return s.cancel(s.sympify(f).subs(sub, simultaneous=True))
        assert all(image(f) == 0 for f in generators)
        assert all(image(f) == 0 for f in jac)
        assert all(image(f) == 0 for f in zero)
        assert all(image(f) != 0 for f in nonzero)
        recovered = []
        for t in rename.values():
            hits = [str(n[i]) for i, a in enumerate(images) if a == t]
            if hits:
                recovered.append(dict(parameter=str(t), coordinate=hits[0]))
        charts.append(dict(name=name, dimension=len(free),
                           source_map={str(v): str(a) for v,a in zip(n,images)},
                           nonzero=[str(image(a)) for a in nonzero],
                           source_nonzero=[str(a) for a in nonzero],
                           source_zero=[str(a) for a in zero],
                           denominators=sorted(set(str(s.denom(a)) for a in images)-{'1'}),
                           directly_recovered_parameters=recovered,
                           all_seven_jacobi_and_all_four_raw_jacobi_zero=True))

    # A: both resonant linear coefficients vanish.
    base = {n0: 0, n6: 0}
    add('A41', base|{n2:0,n7:1-n8,n3:n1*n4,n5:1}, [n4,n1], [n0,n6,n2])
    add('A40', base|{n2:0,n7:1-n8,n3:0,n1:0}, [n4], [n0,n6,n2,n1,n3])
    add('A23', base|{n4:0,n5:s.Rational(1,2),n1:2*n2*n3,n7:-n8}, [n2,n3], [n0,n6,n4])
    add('A20', base|{n4:0,n5:s.Rational(1,2),n1:0,n3:0}, [n2], [n0,n6,n4,n1,n3])
    base |= {n2:0,n4:0}
    add('A13', base|{n5:1,n7:-n8}, [n1,n3], [n0,n6,n2,n4])
    add('A10', base|{n5:1,n3:0}, [n1], [n0,n6,n2,n4,n3])
    add('A03', base|{n1:0,n7:-n8}, [n3], [n0,n6,n2,n4,n1])
    add('A00', base|{n1:0,n3:0}, [], [n0,n6,n2,n4,n1,n3])

    # B: n0 != 0,n6=0, hence n7=1. C is obtained by the actual
    # (x,y,z,w)->(w,z,y,x), tensor-sign-minus involution below.
    B = []
    def b(name, mapping, nonzero, zero):
        B.append((name, {n6:0,n7:1}|mapping, [n0]+nonzero, [n6]+zero))
    b('B41',{n2:0,n8:0,n3:n1*n4,n5:1},[n4,n1],[n2,n8])
    b('B40',{n2:0,n8:0,n3:0,n1:0},[n4],[n2,n8,n1,n3])
    b('B23',{n4:0,n5:s.Rational(1,2),n1:2*n2*n3,n8:-1},[n2,n3],[n4])
    b('B20',{n4:0,n5:s.Rational(1,2),n1:0,n3:0},[n2],[n4,n1,n3])
    b('B13',{n2:0,n4:0,n5:1,n8:-1},[n1,n3],[n2,n4])
    b('B10',{n2:0,n4:0,n5:1,n3:0},[n1],[n2,n4,n3])
    b('B03',{n2:0,n4:0,n1:0,n8:-1},[n3],[n2,n4,n1])
    b('B00',{n2:0,n4:0,n1:0,n3:0},[],[n2,n4,n1,n3])
    invol = {n0:n6,n6:n0,n7:1+n8-n5,n8:n5+n7-1}
    for name,mapping,nz,zz in B:
        add(name,mapping,nz,zz)
        images = [v.subs(invol, simultaneous=True).subs(mapping, simultaneous=True) for v in n]
        # Pushforward rather than re-solving a second set of equations.
        add('C'+name[1:],dict(zip(n,images)),
            [a.subs(invol,simultaneous=True) for a in nz],
            [a.subs(invol,simultaneous=True) for a in zz])
    assert all(s.expand(v.subs(invol,simultaneous=True).subs(invol,simultaneous=True)-v)==0 for v in n)
    Q = -s.Matrix([[P[3-i,3-j] for j in range(4)] for i in range(4)])
    csub = {x:w,y:z,z:y,w:x}
    assert all(s.expand(Q[i,j].subs(csub,simultaneous=True)-P[i,j].subs(invol,simultaneous=True))==0
               for i in range(4) for j in range(4))

    # D: neither resonant linear coefficient vanishes.
    base={n7:1,n8:n5}
    add('D4',base|{n5:0,n8:0,n2:0,n1:-n0*n6,n3:-n0*n6*n4},[n0,n6,n4],[n2,n5,n8])
    add('D2',base|{n4:0,n5:s.Rational(1,2),n8:s.Rational(1,2),n3:0,n1:-2*n0*n6},[n0,n6,n2],[n4,n3])
    add('D3',base|{n2:0,n4:0,n5:-1,n8:-1,n1:-n0*n6/2},[n0,n6,n3],[n2,n4])
    add('D0',base|{n2:0,n4:0,n3:0,n1:-n0*n6/(1-n5)},[n0,n6,1-n5],[n2,n4,n3])
    assert len(charts)==28 and max(c['dimension'] for c in charts)==3

    # A simpler intrinsic orbit obstruction than flags/cofactors: every source
    # point has the nontrivial spectral torus in its projective stabilizer.
    weights=(-2,-1,1,2)
    assert all(s.expand(sum(weights[k]*coords[k]*s.diff(P[i,j],coords[k])
                           for k in range(4))-(weights[i]+weights[j])*P[i,j])==0
               for i,j in it.combinations(range(4),2))
    frontier_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_DOUBLE_A_KNOWN_FRONTIER_20261107.json'
    frontier=json.loads(frontier_path.read_text())
    assert frontier['canonical_double_A_actual_PGL_orbit_dimension']==23
    double_path=ROOT/'reports/CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json'
    double=json.loads(double_path.read_text())['cases'][0]
    assert double['ranks_QQ']['orbit']==48 and double['name']=='nil22_image_double_A'
    h,q1,q2,a,b=s.symbols('h x1 x2 x3 x4')
    m=s.symbols('m0:6')
    denv={str(v):v for v in (h,q1,q2,a,b)+m}
    dcoords=(q1,q2,a,b)
    DP=s.zeros(4)
    for key,value in double['full_six_magnetic_source'].items():
        i,j=(int(v)-1 for v in key.split(','))
        v=s.sympify(value,locals=denv).subs(h,1)
        DP[i,j],DP[j,i]=v,-v
    N=s.Matrix([a,b,0,0]);DN=N.jacobian(dcoords)
    LN=DP.applyfunc(lambda f:s.expand(sum(N[k]*s.diff(f,dcoords[k]) for k in range(4))))-DN*DP-DP*DN.T
    assert LN.applyfunc(s.expand)==s.zeros(4)
    assert DN.rank()==2 and DN*DN==s.zeros(4)

    result = dict(task_id='E_multiple_hyperplane',all_passed=True,
        source_input=dict(path=str(source_path.relative_to(ROOT)),sha256=digest(source_path)),
        producer_sha256=digest(Path(__file__)),spectrum=[-2,-1,1,2],
        complete_affine_brackets=[str(f) for f in brackets],
        all_four_jacobi=[str(f) for f in jac],
        seven_jacobi_generators=[str(f) for f in generators],
        all_raw_coefficients=[str(f) for f in coeffs],
        constant_span_both_directions=dict(raw_from_seven=[[str(a) for a in v] for v in forward],
                                          seven_from_raw=[[str(a) for a in v] for v in reverse]),
        involution={str(k):str(v) for k,v in invol.items()},
        actual_involution_pushforward_checked=True,charts=charts,chart_count=28,
        full_fixed_source_dimension_upper=3,centralizer_dimension=4,
        projective_PGL_saturation_dimension_upper=23,
        double_A_exclusion=dict(accepted_orbit_rank=48,projective_orbit_dimension=23,
            exact_whole_double_A_nilpotent_symmetry=['x3','x4','0','0'],
            nilpotent_rank=2,nilpotent_square_zero=True,
            entire_SS12_source_has_nontrivial_actual_spectral_torus=True,
            proof='Each of finitely many irreducible chart-image closures has dimension<=23. If one contains the23dim double_A orbit, they have the same closure; its dense actual SS12 constructible image meets that orbit. Impossible: that orbit has1dim stabilizer generated by the displayed nilpotent field, whereas SS12 has a nonzero semisimple torus.',
            inputs=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p)) for p in [frontier_path,double_path]]),
        scope='Complete fixed SS12 source only. All constants and resonant linear coefficients kept. No exclusion of moving spectral ratios or other sources. Chart coverage and dense-orbit exclusion require the accompanying case proof.',
        seconds=time.monotonic()-begun,maximum_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with output.open('x') as handle:
        json.dump(result,handle,indent=2)
        handle.write('\n')
    print(json.dumps({k:result[k] for k in ('all_passed','chart_count','full_fixed_source_dimension_upper','projective_PGL_saturation_dimension_upper','seconds','maximum_rss_kib')}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ss12',action='store_true',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    ss12(args.output)
