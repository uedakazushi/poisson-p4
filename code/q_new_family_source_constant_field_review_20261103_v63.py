#!/usr/bin/env python3
"""Exact five first integrals and rational constants field of the source kernel."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import signal
import time

import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = 'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW'


def read(p):
    with (gzip.open(p, 'rt') if p.suffix == '.gz' else p.open()) as f:
        return json.load(f)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for a in iter(lambda: f.read(1 << 20), b''):
            h.update(a)
    return h.hexdigest()


def main(out):
    started = time.monotonic()
    signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS, (2200000000, 2200000000))
    out.mkdir(parents=True, exist_ok=False)
    gate = ROOT/'reports/Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW/report.json'
    kernel_path = ROOT/'reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63/kernel-certificate.json.gz'
    assert read(gate)['all_passed']
    kernel = read(kernel_path)
    u, v, c, X, W, B = s.symbols('u v c X W B')
    oldvars = [u, v, X, W, B]
    raw = list(map(s.sympify, kernel['full26_polynomial_kernel'][:6]))
    K = s.Integer(256)
    rexpr = 384*v*X+K+v*(3072*W-1)
    sexpr = 384*v*X+K-v*(3072*W-1)
    def olddelta(f):
        return s.expand(v/s.Integer(3072)*(sum(s.diff(f, z)*a for z, a in zip(oldvars, raw[:5]))+s.diff(f, c)*c*raw[5]))
    r, t0 = s.symbols('r s')
    newvars = [u, r, t0, v, c, B]
    velocity = [s.Integer(0), r*(r*r-K*K), t0*(t0*t0-K*K),
        v*(r*r+t0*t0)+16*u*(r-t0)**2,
        c*(2*K*K-r*r-t0*t0-32*u*(r-t0)**2/v),
        B*(2*K*K-r*r-t0*t0-48*u*(r-t0)**2/v)]
    assert all(s.cancel(olddelta(f)-z.subs({r: rexpr, t0: sexpr}, simultaneous=True)) == 0
        for f, z in zip([u, rexpr, sexpr, v, c, B], velocity))
    Xback = (r+t0-2*K)/(768*v)
    Wback = (1+(r-t0)/(2*v))/3072
    assert s.cancel(rexpr.subs({X: Xback, W: Wback}, simultaneous=True)-r) == 0
    assert s.cancel(sexpr.subs({X: Xback, W: Wback}, simultaneous=True)-t0) == 0
    def delta(f):
        return s.cancel(sum(s.diff(f, z)*a for z, a in zip(newvars, velocity)))
    zr, zs = 1-K*K/r**2, 1-K*K/t0**2
    kval = zr/zs
    vval = (v-8*u*(r-t0)**2/(r*t0))/(r*t0*zr)
    cval = c*v**2/(r*t0*zr**2)
    bval = B*v**3/((r*t0)**2*zr**3)
    integrals = [u, kval, vval, cval, bval]
    assert s.cancel(delta(zr)-2*K*K*zr) == 0 and s.cancel(delta(zs)-2*K*K*zs) == 0
    assert all(delta(f) == 0 for f in integrals)
    k, iv, ic, ib, t = s.symbols('k Iv Ic IB t')
    fieldnames = [u, k, iv, ic, ib, t]
    parameter = (1/r-1/K)/(1/t0-1/K)
    rback = -K*(t*t-k)/(t*t-2*k*t+k)
    sback = K*(t*t-k)/(t*t-2*t+k)
    assert s.cancel(kval.subs({r: rback, t0: sback}, simultaneous=True)-k) == 0
    assert s.cancel(parameter.subs({r: rback, t0: sback}, simultaneous=True)-t) == 0
    assert s.cancel(rback.subs({k: kval, t: parameter}, simultaneous=True)-r) == 0
    assert s.cancel(sback.subs({k: kval, t: parameter}, simultaneous=True)-t0) == 0
    assert s.cancel(K*K*(1/r**2-kval/t0**2)-(1-kval)) == 0
    assert delta(parameter) == s.cancel(K*(r-t0)*parameter)
    assert delta(parameter) != 0
    # Complete rational inverse, checked triangularly to avoid artificial
    # expansion of all six rational functions into one giant expression.
    vback = iv*r*t0*zr+8*u*(r-t0)**2/(r*t0)
    cback = ic*r*t0*zr**2/v**2
    bback = ib*(r*t0)**2*zr**3/v**3
    assert s.cancel(vval.subs(v, vback)-iv) == 0
    assert s.cancel(cval.subs(c, cback)-ic) == 0
    assert s.cancel(bval.subs(B, bback)-ib) == 0
    assert s.cancel(vback.subs(iv, vval)-v) == 0
    assert s.cancel(cback.subs(ic, cval)-c) == 0
    assert s.cancel(bback.subs(ib, bval)-B) == 0
    # Bind the frozen author proof, its different conic parameter, the actual
    # five-source circuit, and the known finite full-jet involution.
    author_folder = ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63'
    ar, ac, five = [read(author_folder/p) for p in ['report.json', 'certificate.json', 'source-five-circuit.json.gz']]
    assert ar['all_passed']
    assert all(s.cancel(s.sympify(a)-b) == 0 for a, b in zip(ac['five_first_integrals'], integrals))
    at = s.sympify(ac['fibre_parameter_t'])
    assert s.cancel(at-(r/t0-1)/(r/K-1)) == 0
    assert s.cancel(delta(at)-s.sympify(ac['delta_t'])) == 0 and delta(at) != 0
    rp = s.sympify(ac['inverse_r']).subs(s.Symbol('K'), k)
    sp = s.sympify(ac['inverse_s']).subs(s.Symbol('K'), k)
    assert s.cancel(kval.subs({r: rp, t0: sp}, simultaneous=True)-k) == 0
    assert s.cancel(at.subs({r: rp, t0: sp}, simultaneous=True)-t) == 0
    assert s.cancel(rp.subs({k: kval, t: at}, simultaneous=True)-r) == 0
    assert s.cancel(sp.subs({k: kval, t: at}, simultaneous=True)-t0) == 0
    vsection = -98304*X/(147456*X**2-(3072*W-1)**2)
    assert s.cancel(s.sympify(ac['section_v'])-vsection) == 0
    rsec = rexpr.subs(v, vsection)
    ssec = sexpr.subs(v, vsection)
    assert s.cancel(at.subs({r: rsec, t0: ssec}, simultaneous=True)-2) == 0
    assert s.cancel(delta(at).subs({r: rsec, t0: ssec}, simultaneous=True)) != 0
    six_folder = ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63'
    six = read(six_folder/'certificate.json')
    old = read(six_folder/'source-six-circuit.json.gz')
    lookup = {}
    for i, (op, arg) in enumerate(five['nodes']):
        if isinstance(arg, list):
            assert all(j < i for j in arg)
            arg = tuple(arg)
        if op == 'inv':
            assert arg < i
        assert (op, arg) not in lookup
        lookup[op, arg] = i
    assert len(five['nodes']) == ar['five_source_circuit_nodes'] == 1032
    zero, one = lookup['q', '0'], lookup['q', '1']
    def opadd(ids):
        ids = sorted(i for i in ids if i != zero)
        return zero if not ids else ids[0] if len(ids) == 1 else lookup['+', tuple(ids)]
    def opmul(ids):
        if zero in ids:
            return zero
        ids = sorted(i for i in ids if i != one)
        return one if not ids else ids[0] if len(ids) == 1 else lookup['*', tuple(ids)]
    def expression(f):
        if f.is_Rational:
            return lookup['q', str(f)]
        if f.is_Symbol:
            return lookup['v', five['input_names'].index(str(f))]
        if f.is_Add:
            return opadd([expression(a) for a in f.args])
        if f.is_Mul:
            return opmul([expression(a) for a in f.args])
        assert f.is_Pow and f.exp.is_Integer
        z = expression(f.base)
        if f.exp < 0:
            z = lookup['inv', z]
        return opmul([z]*abs(int(f.exp)))
    repl = [expression(f) for f in [u, vsection, c, X, W, B]]
    mapping = []
    for op, arg in old['nodes']:
        if op == 'q':
            z = lookup['q', arg]
        elif op == 'v':
            z = repl[arg]
        elif op == 'inv':
            z = lookup['inv', mapping[arg]]
        elif op == '+':
            z = opadd([mapping[j] for j in arg])
        else:
            assert op == '*'
            z = opmul([mapping[j] for j in arg])
        mapping.append(z)
    def remap(a):
        if isinstance(a, int):
            return mapping[a]
        if isinstance(a, list):
            return [remap(v) for v in a]
        return {k: remap(v) for k, v in a.items()}
    assert remap(old['outputs']) == five['outputs']
    q = s.Matrix(s.symbols('q0:4'))
    G = s.diag(-1, 1, -1, 1)
    rawC, rawL = s.Matrix(six['raw_value4']), s.Matrix(six['raw_linear4'])
    invol = {u: -u, v: -v, c: -c, X: -X, W: W, B: B}
    assert (G*(c*rawC)*G.T-(c*rawC).subs(invol, simultaneous=True)).applyfunc(s.cancel) == s.zeros(4)
    assert (G*(c*rawL).subs(dict(zip(q, G*q)), simultaneous=True)*G.T-(c*rawL).subs(invol, simultaneous=True)).applyfunc(s.expand) == s.zeros(4)
    invnew = {u: -u, v: -v, c: -c, r: t0, t0: r, B: B}
    expected = [-u, 1/kval, -kval*vval, -kval**2*cval, -kval**3*bval]
    assert all(s.cancel(f.subs(invnew, simultaneous=True)-e) == 0 for f, e in zip(integrals, expected))
    assert s.cancel(at.subs(invnew, simultaneous=True)-at/(at-1)) == 0
    assert s.cancel(vsection.subs(X, -X)+vsection) == 0
    # On the section, the involution negates only the three variables u,c,X.
    # QQ(u,c,X,W,B)=QQ(X²,u/X,c/X,W,B)(X), X² in the smaller field,
    # and the nontrivial involution X -> -X proves the fixed-field assertion.
    evidence = dict(original_kernel_rescaling='delta=(v/3072)*K_raw',
        r=str(rexpr), s=str(sexpr), X_inverse=str(Xback), W_inverse=str(Wback),
        transformed_variable_order=list(map(str, newvars)), transformed_velocity=list(map(str, velocity)),
        first_integral_order=['u', 'k', 'Iv', 'Ic', 'IB'], first_integrals=list(map(str, integrals)),
        conic_identity='256^2*(1/r^2-k/s^2)=1-k',
        rational_parameter=str(parameter), r_inverse=str(rback), s_inverse=str(sback),
        v_inverse=str(vback), c_inverse=str(cback), B_inverse=str(bback),
        derivative_parameter=str(delta(parameter)),
        rational_field_generators=list(map(str, fieldnames)),
        field_identity='QQ(u,v,c,X,W,B)=QQ(u,k,Iv,Ic,IB)(t)',
        constants_identity='ker(delta)=QQ(u,k,Iv,Ic,IB)',
        symbolic_section_v=str(vsection), source_five_circuit_nodes=1032,
        known_finite_involution_section_fixed_field='QQ(X^2,u/X,c/X,W,B)',
        no_actual_beta_quotient_field_equality_or_index_one_asserted=True)
    (out/'independent-certificate.json').write_text(json.dumps(evidence, indent=2)+'\n')
    paths = [Path(__file__), ROOT/f'research/notes/{STEM}.md', gate, kernel_path,
        *[author_folder/p for p in ['report.json', 'certificate.json', 'source-five-circuit.json.gz']],
        ROOT/'code/q_new_family_source_rational_quotient_20261103_v63.py',
        ROOT/'research/notes/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63.md',
        six_folder/'certificate.json', six_folder/'source-six-circuit.json.gz']
    result = dict(id=STEM, all_passed=True,
        original26_kernel_and_exact_source_coordinate_rescaling_bound=True,
        all6_transformed_derivation_identities_and_rational_X_W_inverse=True,
        all5_rational_first_integrals_have_exact_zero_derivative=True,
        complete_rational_source_field_inverse_and_algebraic_independence=True,
        exact_source_constants_field_by_nonzero_derivative_on_rational_t=True,
        author_different_conic_parameter_and_inverse_exactly_verified=True,
        symbolic_t2_section_and_all1032_nodes_composed_from_reviewed_source6=True,
        known_full_rawjet_order2_equivalence_and_rational_section_fixed_field=True,
        actual_NEW_beta_quotient_field_index_or_birationality_claimed=False,
        no_new_point_native_solver_high_degree_probe_or_GB=True,
        input_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths},
        evidence_sha256={'independent-certificate.json': sha(out/'independent-certificate.json')},
        seconds=time.monotonic()-started)
    (out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ['all_passed', 'exact_source_constants_field_by_nonzero_derivative_on_rational_t', 'seconds']}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True, type=Path)
    main(p.parse_args().output.resolve())
