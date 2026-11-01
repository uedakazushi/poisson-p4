#!/usr/bin/env python3
"""Uniform quadratic rank15 at r=-4/25 with polynomial basis and t2 reduction."""
import argparse,hashlib,json,pickle,subprocess,sys,tempfile,time
from pathlib import Path
from itertools import combinations_with_replacement as cwr
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();start=time.time()
    engines=['code/s61_flag_boundary_eight_frame_engine_20261101.py','code/s61_flag_true_core_frame_engine_20261101.py','code/s61_flag_linear_artifacts_quadratic_engine_20261101.py']
    with tempfile.TemporaryDirectory(prefix='flag-true-core-quad-')as temp:
        wd=Path(temp)
        for i,engine in enumerate(engines):subprocess.run([sys.executable,str(ROOT/engine),str(wd if i<2 else wd/'core-0')],cwd=ROOT,check=True)
        wd=wd/'core-0';data=pickle.loads((wd/'flag-b1-universal-quadratic-cache.pkl').read_bytes());q=json.loads((wd/'flag-b1-universal-quadratic.json').read_text());frame=json.loads((wd/'flag-b1-normalized-frame.json').read_text())
    variables=list(map(s.Symbol,data['variables']));decode=lambda p:sum(s.Rational(co)*s.prod(x**e for x,e in zip(variables,ex))for ex,co in p)
    A=s.Matrix([[decode(p)for p in row]for row in data['A2']]);rows=q['pivot_rows'];B=A[rows,:];_,columns=B.subs({x:0 for x in variables}).to_DM().rref();columns=list(columns);assert len(rows)==len(columns)==15
    M=B[:,columns].to_DM();det=M.domain.to_sympy(M.det());assert det.is_Rational and det
    Q=(M.to_field().inv()*B.to_DM().to_field()).to_Matrix().applyfunc(s.cancel)
    assert all(s.denom(x).is_Number for x in Q)
    assert(A-A[:,columns]*Q).applyfunc(s.expand)==s.zeros(135,55)
    t2=s.Matrix([[decode(p)for p in data['W2'][44]]]);assert any(t2)
    assert(t2-t2[:,columns]*Q).applyfunc(s.expand)==s.zeros(1,55)
    monomials=list(cwr(range(10),2));active=[4,5,7,8];active_columns=[j for j,m in enumerate(monomials)if all(a in active for a in m)]
    rref=s.Matrix(q['quadratic_rref']).applyfunc(s.sympify)
    active_rows=[i for i in range(15)if all(not rref[i,j]for j,m in enumerate(monomials)if not all(a in active for a in m))]
    assert len(active_rows)==6;C=rref[active_rows,:];assert not any(x.free_symbols for x in C)
    assert(C-C[:,columns]*Q).applyfunc(s.expand)==s.zeros(6,55)
    active_cubics=list(cwr(active,3));products=s.zeros(24,20)
    for i in range(6):
        for j,a in enumerate(active):
            for k in active_columns:products[4*i+j,active_cubics.index(tuple(sorted(monomials[k]+(a,))))]+=C[i,k]
    selected=list(products.T.to_DM().rref()[1]);assert len(selected)==20
    inverse=products[selected,:].inv();assert inverse*products[selected,:]==s.eye(20)
    encode=lambda matrix:[[str(x)for x in row]for row in matrix.tolist()]
    paths=[str(Path(__file__).relative_to(ROOT)),*engines,'code/s61_model.py','code/s61_flag_transverse_obstruction_20261101.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json','reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json']
    result=dict(all_passed=True,core='-4/25',variables=data['variables'],jacobian_rank=35,tangent_dimension=10,quadratic_rank=15,hilbert_function_through2=[1,10,40],constant_minor=dict(rows=rows,columns=columns,value=str(det)),polynomial_quadratic_basis=encode(Q),basis_parameter_degree=max(s.Poly(x,*variables).total_degree()for x in Q),raw_t2_nonzero=True,t2_in_quadratic_ideal=True,t_order_at_least3=True,active_variables=active,constant_active_quadrics=encode(C[:,active_columns]),active_monomials=[monomials[k]for k in active_columns],active_cubic_product_selected_rows=selected,active_cubic_product_inverse=encode(inverse),active_quotient_hilbert_function=[1,4,4],active_quotient_length=9,scope='All seven parameter values on the normalized true core, including b8=0.',limitations=['Only true quadratic data and t in m3.','No universal third or higher initial ideal, complete local ring, t-nilpotence or saturation removal.'],proof_note='research/notes/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.md',seconds=time.time()-start,input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in paths})
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('TRUE CORE UNIFORM QUADRATIC PASS',result['seconds'],flush=True)
if __name__=='__main__':main()
