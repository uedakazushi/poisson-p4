"""Full positive-weight central slices at a=0, K=0 in the nonsquare chart."""
import sympy as sp
from sympy.polys.rings import ring
from s61_binary_endpoint_full_transverse_20261102 import (
    PRIME, CUTS, B, S61, specialize, load_report)

PIVOT_STEPS = [[71, 46], [75, 8], [73, 45], [70, 34], [0, 9], [3, 12], [9, 4], [62, 25], [65, 30], [72, 42], [74, 7], [117, 21], [69, 18], [87, 36], [1, 10], [2, 11], [6, 15], [8, 17], [49, 2], [64, 29], [90, 38], [99, 6], [122, 28], [92, 39], [4, 13], [5, 14], [7, 16]]

def transversal(b,c,d,e):
    model = S61()
    chart = load_report('S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102')
    old = load_report('S61_BINARY_DEEPER_TRIANGULAR_CHART_20261101')['linear']
    fixed, columns = old['fixed'], old['cols']
    weights = chart['local_coordinate_weights']+[1,1]
    grades = chart['equation_weights']
    normal = sp.symbols('n0 n1')
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    beta = dict.fromkeys(one,sp.S.Zero)
    beta.update({13:b,34:c,35:d,48:e,49:-5*c,11:normal[0],12:normal[1]})
    _, equations = specialize(model,beta,True)
    rows = []
    for grade,equation in zip(grades,equations):
        local = {}
        for monomial,coefficient in equation.items():
            # The chosen mixed coefficient equals one in this projective chart.
            indices = tuple(columns.index(i) for i in monomial if i != fixed)+(44,)*(2-len(monomial))
            if any(i in CUTS for i in indices):
                continue
            for (n0,n1),coefficient in sp.Poly(coefficient,*normal).terms():
                mon = tuple(sorted(indices+(45,)*n0+(46,)*n1))
                local[mon] = local.get(mon,0)+coefficient
        rows.append({mon:coefficient for mon,coefficient in local.items() if coefficient})
        assert all(sum(weights[j] for j in mon) == grade for mon in rows[-1])
    assert [weights[i] for i in CUTS].count(0) == 2
    assert all(w > 0 for i,w in enumerate(weights) if i not in CUTS)
    return rows,weights,grades


def reduce_line(rows,weights,grades):
    """Eliminate 27 coordinates using minors invertible on the entire test line."""
    variables = [i for i in range(47) if i not in CUTS]
    free = [i for i in variables if i not in {j for _,j in PIVOT_STEPS}]
    assert sorted(weights[i] for i in free) == [1]*7+[2]*5+[3]*3
    field = sp.GF(PRIME)
    polynomials,*x = ring('b,'+','.join(f'x{i}' for i in range(len(free))),field)
    coordinates = [polynomials.zero]*47
    for j,v in enumerate(free):
        coordinates[v] = x[j+1]
    coefficient_cache = {}
    def convert(f):
        if f not in coefficient_cache:
            v = polynomials.zero
            for (power,),scalar in sp.Poly(f,B).terms():
                coefficient = int(scalar.p)*pow(int(scalar.q),-1,PRIME) % PRIME
                v += polynomials(coefficient)*x[0]**power
            coefficient_cache[f] = v
        return coefficient_cache[f]
    def evaluate(row):
        value = polynomials.zero
        for mon,c in row.items():
            term = convert(c)
            for i in mon:
                term *= coordinates[i]
            value += term
        return value

    determinants = []
    for grade in range(1,5):
        steps = [(i,j) for i,j in PIVOT_STEPS if weights[j] == grade]
        inds,cols = zip(*steps)
        block = sp.Matrix([[rows[i].get((j,),sp.S.Zero) for j in cols] for i in inds])
        determinant = sp.factor(block.det(method='domain-ge'))
        assert determinant.is_Rational and int(determinant.p) % PRIME and int(determinant.q) % PRIME
        determinants.append(str(determinant))
        inverse = block.to_DM().to_field().inv().to_Matrix()
        values = [evaluate(rows[i]) for i in inds]
        for j,index in enumerate(cols):
            coordinates[index] = -sum((convert(inverse[j,h])*value for h,value in enumerate(values)
                                       if value),polynomials.zero)
    residual = [evaluate(row) if grade <= 4 else polynomials.zero for row,grade in zip(rows,grades)]
    assert all(not residual[i] for i,_ in PIVOT_STEPS)
    return residual, [weights[i] for i in free], determinants
