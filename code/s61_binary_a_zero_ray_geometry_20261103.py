#!/usr/bin/env python3
"""Exact incidence and orbit checks for the seven nonzero a=0 rays.

This does not call the transverse-length producer.  It reconstructs the
primitive beta representation, the highest-equation kernel, and the orbit
motion, and verifies that the SAME incidence dominates the beta base.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time

import sympy as s

from s61_model import S61, X
from s61_nil22_infinity_20261101 import beta_matrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from toric_closed_jet_coordinates_20261102 import coords, pairs, m2
from weight1_structure import COORDS, infinitesimal_action, sp_basis

ROOT = Path(__file__).resolve().parents[1]
PRIME = 32003


def encode(matrix):
    return [[str(z) for z in row] for row in matrix.tolist()]


def unit(z):
    z = s.Rational(z)
    return int(z.p) % PRIME != 0 and int(z.q) % PRIME != 0


def minor(matrix, size=None):
    columns = list(matrix.to_DM().rref()[1])
    if size is not None:
        assert len(columns) == size
    rows = list(matrix[:, columns].T.to_DM().rref()[1])
    value = matrix[rows, columns].det(method='domain-ge')
    assert value
    return dict(rows=rows, columns=columns, determinant=str(value))


def verify():
    started = time.time()
    model = S61()
    one = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i, w in enumerate(model.reduced['parameter_weights']) if w == 2]
    extra_path = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    source = model.reduced['equations'] + [json.loads(extra_path.read_text())['universal_new_equation']]
    point_path = ROOT/'code/certificates/s61_binary_a_zero_points_20261102.json'
    points = json.loads(point_path.read_text())
    data = coords()
    x, y, z, w = X[1:]
    spatial = s.Matrix([x, y, z, w])
    zz, bu, beta_form = primitive_beta(model)
    vector = s.Matrix([beta_form[i,j].diff(zz[k]) for (i,j),k in COORDS])
    basis = vector.jacobian(bu)
    rows = list(basis.T.to_DM().rref()[1])
    left = basis[rows, :].inv()
    acts = s.Matrix.hstack(*[(left*infinitesimal_action(vector, a)[rows, :]).applyfunc(s.expand)
                            for a in sp_basis()])
    assert all(basis*acts[:, j] == infinitesimal_action(vector, a)
               for j, a in enumerate(sp_basis()))
    records = []
    expected = {'xzw': (7,10,2), 'T2': (8,9,2), 'double root': (7,10,3),
                'smooth': (7,10,1), 'T5': (6,12,23), 'T6': (8,9,2)}
    for entry in points:
        q = s.zeros(45,1)
        beta = s.zeros(16,1)
        for j, c in entry['q'].items():
            q[int(j)] = s.Rational(c)
        for j, c in entry['beta'].items():
            beta[int(j)] = s.Rational(c)
        physical = data['E']*q
        tensor = s.zeros(4)
        for k, (i,j) in enumerate(pairs):
            tensor[i,j] = sum(physical[k*10+h]*spatial[v]*spatial[t]
                              for h,(v,t) in enumerate(m2))
        tensor -= tensor.T
        kernel_matrix = beta_matrix(model, tensor)
        H = s.Matrix.hstack(*kernel_matrix.nullspace())
        assert kernel_matrix*beta == s.zeros(kernel_matrix.rows,1)
        assert H*H.gauss_jordan_solve(beta)[0] == beta
        parameters = s.symbols('b:'+str(H.cols))
        relative_beta = H*s.Matrix(parameters)
        relative_values = dict(zip(one, relative_beta)) | dict(zip(two, q))
        for equation in source:
            highest = sum(s.Rational(n,d)*s.prod(relative_values[i] for i in mon)
                          for mon,n,d in equation if sum(i in two for i in mon) == 2)
            assert s.expand(highest) == 0
        orbit = []
        fixed = 21
        assert q[fixed] == 1
        for a in sp_basis():
            moved = -a*tensor-tensor*a.T+sum(((a*spatial)[i]*tensor.diff(spatial[i])
                                             for i in range(4)),s.zeros(4))
            coefficients = s.Matrix([s.Poly(moved[i,j],*spatial).coeff_monomial(spatial[v]*spatial[t])
                                     for i,j in pairs for v,t in m2])
            column = data['Pi']*coefficients
            orbit.append(column-q*column[fixed])
        O = s.Matrix.hstack(*orbit)
        orbit_minor = minor(O)
        hrows = list(H.T.to_DM().rref()[1])
        incidence_minor = H[hrows,:].det()*s.Rational(orbit_minor['determinant'])
        assert unit(incidence_minor)
        d,h,m = expected[entry['type']]
        assert O.to_DM().rank() == d and H.cols == h
        assert all(int(s.denom(c)) % PRIME for c in H)
        actual_beta_motion = acts.subs(dict(zip(bu,beta)))
        sample = H*s.Matrix(range(2,H.cols+2))
        sample_motion = acts.subs(dict(zip(bu,sample)))
        beta_minor = minor(H.row_join(sample_motion),16)
        actual_rank = H.row_join(actual_beta_motion).to_DM().rank()
        point = dict(zip(one,beta)) | dict(zip(two,q)) | {61:s.S.Zero}
        variables = [i for i in range(62) if i != two[fixed]]
        loc = dict(zip(variables,range(61)))
        jacobian = s.zeros(135,61)
        for row,equation in enumerate(source):
            for mon,n,den in equation:
                ids = tuple(mon) + (61,)*(2-sum(i in two for i in mon))
                for variable,exponent in Counter(ids).items():
                    if variable not in loc:
                        continue
                    rest = list(ids)
                    rest.remove(variable)
                    jacobian[row,loc[variable]] += s.Rational(n,den)*exponent*s.prod(point[i] for i in rest)
        ambient_rank = jacobian.to_DM().rank()
        assert ambient_rank == 61-d-h+(0 if m==1 else {'xzw':-1,'T2':-1,'double root':-2,'T5':-3,'T6':-1}[entry['type']])
        if entry['type'] == 'smooth':
            assert ambient_rank == 44 and d+h == 17
            hamiltonian = None
        else:
            hamiltonian = s.integrate(tensor[1,2],w)
            discrepancy = s.expand(-tensor[1,3]-s.diff(hamiltonian,z))
            assert not discrepancy.has(w)
            hamiltonian += s.integrate(discrepancy,z)
            check = s.zeros(4)
            check[1,2] = s.diff(hamiltonian,w)
            check[1,3] = -s.diff(hamiltonian,z)
            check -= check.T
            assert (check-tensor).applyfunc(s.expand) == s.zeros(4)
            if entry['type']=='xzw':
                quadratic = s.cancel(hamiltonian/x)
                assert not quadratic.has(x,y) and s.discriminant(quadratic.subs(w,1),z) != 0
            elif entry['type']=='T2':
                assert s.expand(hamiltonian-(x*w*w/2-z*z*w/44)) == 0
            elif entry['type']=='double root':
                assert s.expand(hamiltonian-w*w*(8*x-z)/16) == 0
                assert s.expand(hamiltonian.subs(z,z+8*x)+z*w*w/16) == 0
            elif entry['type']=='T5':
                assert s.expand(hamiltonian-(x*w*w/2-x*x*z)) == 0
            elif entry['type']=='T6':
                assert s.expand(hamiltonian-(x*w*w/2-z**3/15)) == 0
        records.append(dict(slope=entry['slope'],type=entry['type'],tensor=encode(tensor),
                            hamiltonian=str(s.factor(hamiltonian)) if hamiltonian is not None else None,
                            highest_kernel=encode(H),all135_highest_zero_on_whole_kernel=True,
                            orbit_dimension=d,H_dimension=h,incidence_dimension=d+h,
                            incidence_unit_minor=str(incidence_minor),orbit_minor=orbit_minor,beta_cut_rows=hrows,
                            generic_beta_sample=list(map(str,sample)),beta_dominance_minor=beta_minor,
                            beta_rank_at_given_point=int(actual_rank),ambient_full_jacobian_rank=int(ambient_rank),
                            known_generic_multiplicity=m))
        print(entry['slope'],entry['type'],'beta rank',actual_rank,'generic beta rank16','ambient rank',ambient_rank,flush=True)
    dependencies = [ROOT/'research/notes/S61_NONBINARY_CANONICAL_EXHAUSTION_20261101.md',
                    ROOT/'research/notes/S61_R0_MIXED_CURVE_EXHAUSTION_20261101.md',
                    ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101.md',
                    ROOT/'research/notes/S61_R0_T5_FORMAL_RANK23_20261101.md']
    paths = [Path(__file__),point_path,extra_path,*dependencies,
             ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',
             ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
             ROOT/'code/toric_closed_jet_coordinates_20261102.py']
    return dict(id='S61_BINARY_A_ZERO_RAY_GEOMETRY_20261103',all_passed=True,field='QQ',
                rays=records,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                limitations=['Generic multiplicities are proof dependencies, not re-proved by these linear checks.',
                             'Transverse upper lengths are independently produced elsewhere.',
                             'No full affine-fiber exhaustion follows from these rays.'],seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print(json.dumps({'all_passed':True,'rays':len(report['rays']),'seconds':report['seconds']}))
