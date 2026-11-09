"""Exact E8 first-jet necessary equations on v21 fixed-support models.

No parameter localization or radical is taken. Input files are read-only.
Matrices act on column vectors of conormal covectors.
"""
import argparse
import itertools
import json
import time
from pathlib import Path
import sympy as s


def normalized_coefficients(expressions, t, parameters):
    result = {}
    raw_count = 0
    for f in expressions:
        for a in s.Poly(s.expand(f), t).all_coeffs():
            if a == 0:
                continue
            raw_count += 1
            g = s.Poly(a, *parameters, domain=s.QQ).monic().as_expr()
            result[s.srepr(g)] = str(g)
    return [result[k] for k in sorted(result)], raw_count


def build(source, outdir):
    started = time.monotonic()
    full = json.loads(source.read_text())
    model = full['linear_model']
    x = s.symbols('x0:5')
    z = s.symbols('z0:' + str(model['linear_dimension']))
    t = s.Symbol('t')
    locals_ = {str(v): v for v in x + z + (t,)}
    P = s.zeros(5)
    for key, expression in model['P_entries'].items():
        i, j = map(int, key)
        P[i,j] = s.sympify(expression, locals=locals_)
        P[j,i] = -P[i,j]
    residues = [s.sympify(v, locals=locals_) for v in model['residue_expressions']]
    components = []
    all_conditions = []
    for number, component in enumerate(model['components']):
        r = [s.sympify(v, locals=locals_) for v in component['parametrization']]
        pivot = next(i for i in range(5) if r[i] == 1)
        tangent = next(i for i in range(5) if r[i] == t)
        normal = [i for i in range(5) if i not in (pivot, tangent)]
        chart = [i for i in range(5) if i != pivot]
        point = dict(zip(x, r))
        p = s.zeros(5)
        for i in chart:
            for j in chart:
                p[i,j] = s.expand((P[i,j] - x[i]*P[pivot,j] + x[j]*P[pivot,i]).subs(x[pivot],1))
        assert all(s.expand(p[i,j].subs(point)) == 0 for i in chart for j in chart)
        H = s.Matrix(3,3,lambda k,j: s.expand(s.diff(p[tangent,normal[j]],x[normal[k]]).subs(point)))
        B = s.zeros(3)
        for column,(a,b) in enumerate([(1,2),(2,0),(0,1)]):
            i,j = normal[a], normal[b]
            bracket = p[i,j] - s.diff(r[i],t).subs(t,x[tangent])*p[tangent,j] - s.diff(r[j],t).subs(t,x[tangent])*p[i,tangent]
            for k in range(3):
                B[k,column] = s.expand(s.diff(bracket,x[normal[k]]).subs(point))
        R = s.expand(sum(residues[k]*s.sympify(b,locals=locals_) for k,b in zip(component['residue_indices'],component['residue_basis'])))
        assert all(s.expand(v)==0 for v in B-B.T)
        assert s.expand(s.trace(H)+R) == 0
        families = {
            'heisenberg_rank_one': [B.extract(rows,cols).det() for rows in itertools.combinations(range(3),2) for cols in itertools.combinations(range(3),2)],
            'e8_trace_square': [s.trace(H*H)-s.Rational(7,18)*R**2],
            'e8_determinant': [H.det()+s.Rational(1,36)*R**3],
            'e8_center': list((H+R*s.eye(3)/2)*B),
        }
        equations = {}
        counts = {}
        for name, expressions in families.items():
            equations[name], raw_count = normalized_coefficients(expressions,t,z)
            counts[name] = {'raw_nonzero_coefficients':raw_count,'distinct_monic_equations':len(equations[name])}
            all_conditions.extend(expressions)
        components.append({'component':number,'pivot':pivot,'tangent':tangent,'normal':normal,
                           'parametrization':[str(v) for v in r], 'R':str(R),
                           'H':[[str(v) for v in row] for row in H.tolist()],
                           'B':[[str(v) for v in row] for row in B.tolist()],
                           'equations':equations,'counts':counts,
                           'trace_and_symmetry_checked':True})
    combined,raw_count=normalized_coefficients(all_conditions,t,z)
    report={'name':model['name'],'scope':'Necessary E8 generic first-jet equations; not sufficient for Jacobi or normality',
            'linear_dimension':model['linear_dimension'],'components':components,
            'all_e8_equations':combined,'raw_nonzero_e8_coefficients':raw_count,
            'distinct_monic_e8_equations':len(combined),'parameter_division':False,
            'linear_model_relative_path':'../../previous_v21/data/tree_jacobi_models/'+source.name,
            'all_structural_checks_passed':True}
    outdir.mkdir(exist_ok=True,parents=True)
    (outdir/source.name).write_text(json.dumps(report,indent=2)+'\n')
    summary={'name':model['name'],'components':len(components),'e8_equations':len(combined),
             'seconds':round(time.monotonic()-started,2)}
    print(json.dumps(summary),flush=True)
    return summary


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('source_directory',type=Path)
    ap.add_argument('output_directory',type=Path)
    ap.add_argument('--model',action='append')
    args=ap.parse_args()
    reports=[]
    for path in sorted(args.source_directory.glob('*.json')):
        if args.model and path.stem not in args.model: continue
        d=json.loads(path.read_text())
        if 'linear_model' not in d:continue
        reports.append(build(path,args.output_directory))
    (args.output_directory/'summary.json').write_text(json.dumps({'models':reports,'all_structural_checks_passed':True},indent=2)+'\n')
