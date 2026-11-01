#!/usr/bin/env python3
"""All four minimal isolated-axis presentations are removed for generic beta."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,infinitesimal_action,sp_basis,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def run():
    start=time.time();m=S61();weights=m.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2];g=s.symbols('g0:9');k,K,L,U,V=s.symbols('k K L U V');b=s.symbols('b0:6')
    exhausted_path=ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json';family_path=ROOT/'reports/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.json';exhausted=json.loads(exhausted_path.read_text());family=json.loads(family_path.read_text());z,u,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(z[h]) for (i,j),h in COORDS]);bas=vec.jacobian(u);rr=list(bas.T.rref()[1]);left=bas[rr,:].inv();actions=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);out=[]
    for block in exhausted['families']:
        kind=block['label'];qexpr=s.Matrix(list(map(s.sympify,(family if kind=='M0' else family['M1_family'])['q_coordinates'])))
        for axis_case in block['isolated_axis_charts']:
            axis=axis_case['axis'];dependent=K if axis==6 else L;free_base=L if axis==6 else K;base_solution=s.solve(s.sympify(axis_case['base_divisor']),dependent)[0];axis_solution=s.solve(s.sympify(axis_case['axis_equation']),g[axis])[0];gs=dict.fromkeys(g,s.S.Zero);gs.update({s.sympify(a):s.sympify(v) for a,v in axis_case['rational_solution'].items()});gs[g[axis]]=g[axis]
            gs={a:s.cancel(v.subs(g[axis],axis_solution).subs(dependent,base_solution)) for a,v in gs.items()}
            assert all(s.cancel(s.sympify(f).subs(gs,simultaneous=True).subs(dependent,base_solution))==0 for f in block['vertex_equations'])
            assert s.cancel(s.sympify(block['G6']).subs(gs,simultaneous=True).subs(dependent,base_solution))==0
            bv=s.Matrix([0,0,0,0,k,b[0],U,b[1],b[2],V,b[3],b[4],b[5],2*k+1,K-5*V,5*U-L]).subs(dependent,base_solution);params=[k,free_base,U,V]+list(b);qfamily=qexpr.subs(gs,simultaneous=True);sample={k:1,free_base:5 if axis==6 else 19,U:2,V:3,**dict(zip(b,[2,3,5,7,11,13]))};beta=bv.subs(sample);q=qfamily.subs(sample);assert all(v.is_Rational for v in q)
            pt=s.zeros(61,1)
            for i,v in zip(one,beta):pt[i]=v
            for i,v in zip(two,q):pt[i]=v
            _,F=specialize(m,dict(zip(one,beta)),True);A=s.zeros(135,45);tc=s.zeros(135,1)
            for row,f in enumerate(F):
                val=0
                for mon,v in f.items():
                    if len(mon)==2:
                        i,j=mon;val+=v*q[i]*q[j];A[row,i]+=v*q[j];A[row,j]+=v*q[i]
                    elif len(mon)==1:tc[row]+=v*q[mon[0]]
                assert val==0
            assert len(A.to_DM().rref()[1])==44;fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];full=A[:,cols].row_join(tc);assert len(full[:134,:].to_DM().rref()[1])==45
            rows=list(full[:134,:].T.to_DM().rref()[1]);minor=full[rows,:];dm=minor.to_DM();det=dm.domain.to_sympy(dm.det());assert det!=0
            D=bv.jacobian(params).subs(sample).row_join(actions.subs(dict(zip(u,beta))));dom=minor_witness(D);assert dom['rank']==16
            out.append(dict(kind=kind,axis=axis,base_equation=axis_case['base_divisor'],source_variables=list(map(str,params)),source_beta=list(map(str,bv)),source_G={str(a):str(v) for a,v in gs.items()},sample_point61=list(map(str,pt)),source_beta_differential=enc(D),dominance_minor=dom,leading_rank=44,chart_fixed_q=fixed,chart_columns=cols,homogenized_jacobian=enc(full),original134_minor_rows=rows,original134_minor_determinant=str(det),homogenized_rank=45,local_ring='C',t_zero=True,whole_generic_beta_axis_closure_removed=True))
            print(json.dumps({'kind':kind,'axis':axis,'full_rank':45}),flush=True)
    paths=[Path(__file__),exhausted_path,family_path,ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_ALL_AXIS_REMOVAL_20261101',all_passed=True,field='QQ',cases=out,conclusion='Every isolated g6/g7 axis presentation in M0 and M1 is absent from the t-saturated closure over sufficiently general beta.',limitations=['Four presentations need not be four distinct irreducible components.','No generic degree of an axis cover is computed.','Does not remove the positive-dimensional tail, free-line, or main-curve families.'],proof_note='research/notes/S61_MINIMAL_ALL_AXIS_REMOVAL_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
if __name__=='__main__':main()
