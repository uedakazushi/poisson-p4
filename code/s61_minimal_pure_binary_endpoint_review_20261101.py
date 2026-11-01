#!/usr/bin/env python3
"""Independent raw135 geometry and local rank audit of pure-cubic binary endpoint."""
import argparse,hashlib,json,time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)


def run():
    start=time.time();rp=ROOT/'reports/S61_MINIMAL_PURE_BINARY_ENDPOINT_20261101.json';report=json.loads(rp.read_text())
    assert report['all_passed']
    for name,digest in report['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2]
    qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)}
    equations=model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    encoded=[[(tuple(qi[i] for i in mon if i in qi),tuple(ui[i] for i in mon if i in ui),s.Rational(n,d)) for mon,n,d in eq] for eq in equations]
    geo=report['source_geometry'];beta=M(geo['beta']);parameters=list(map(s.Symbol,geo['parameters']));deltas=[M(geo['actual_line_delta'])]
    q=M(report['point_local_data']['point61'])[two,:];fixed=next(i for i,v in enumerate(q) if v);products={}
    assert q[fixed]==105 and deltas[0][35]==s.Rational(1,2)
    for delta in deltas:
        for row,eq in enumerate(encoded):
            coefficients=[0,0,0]
            for qm,bm,c in eq:
                if len(qm)!=2:continue
                i,j=qm;values=[q[i]*q[j],q[i]*delta[j]+delta[i]*q[j],delta[i]*delta[j]]
                if not any(values):continue
                if bm not in products:products[bm]=s.prod(beta[k] for k in bm)
                for k,value in enumerate(values):coefficients[k]+=c*products[bm]*value
            assert all(s.cancel(value)==0 for value in coefficients),row

    def local_audit(local,expected):
        point=M(local['point61']);b=point[one,:];qq=point[two,:]
        fx=next(i for i,v in enumerate(qq) if v);chart=[i for i in range(45) if i!=fx]
        derivative=s.zeros(135,45);tcol=s.zeros(135,1);polys=[];cache={}
        for row,eq in enumerate(encoded):
            f={}
            for qm,bm,c in eq:
                if bm not in cache:cache[bm]=s.prod(b[k] for k in bm)
                co=c*cache[bm];f[qm]=f.get(qm,0)+co
                if len(qm)==2:
                    i,j=qm;derivative[row,i]+=co*qq[j];derivative[row,j]+=co*qq[i]
                elif len(qm)==1:tcol[row]+=co*qq[qm[0]]
            assert sum(c*qq[m[0]]*qq[m[1]] for m,c in f.items() if len(m)==2)==0
            polys.append(f)
        A=derivative[:,chart].row_join(tcol);T=M(local['tangent_basis']);L=M(local['cokernel_rows']);n=T.cols
        assert A==M(local['jacobian']) and A*T==s.zeros(135,n) and L*A==s.zeros(L.rows,45)
        assert T[-1,:]==s.zeros(1,n)
        forms={i:([0]*n if i==fx else list(T[chart.index(i),:])) for i in range(45)}
        monos=list(combinations_with_replacement(range(n),2));indices={m:i for i,m in enumerate(monos)};H=s.zeros(135,len(monos))
        for row,f in enumerate(polys):
            for qm,c in f.items():
                ff=[forms[i] for i in qm]+[list(T[-1,:])]*(2-len(qm))
                for i in range(n):
                    for j in range(n):H[row,indices[tuple(sorted((i,j)))]]+=c*ff[0][i]*ff[1][j]
        projected=L*H;assert projected==M(local['quadratic_obstruction_basis'])
        hh=[]
        for coeff in projected.tolist():
            hs=s.zeros(n)
            for c,(i,j) in zip(coeff,monos):hs[i,j]+=c;hs[j,i]+=c
            hh+=hs.tolist()
        ranks=(derivative.rank(),A.rank(),projected.rank(),s.Matrix(hh).rank())
        assert ranks==expected,ranks
        return A,T,chart,point,ranks

    A,T,chart,point,ranks=local_audit(report['point_local_data'],(34,35,5,7))
    free=[i for i in range(45) if i not in A.rref()[1]];assert T[free,:]==s.eye(10)
    sample={s.Symbol(k):s.Rational(v) for k,v in geo['beta_sample'].items()}
    assert beta.subs(sample)==point[one,:]
    for delta,expected in zip(deltas,[report['actual_line_tangent_coordinates']]):
        dt=delta.subs(sample)[chart,:].col_join(s.zeros(1,1));dc=dt[free,:]
        assert T*dc==dt and list(map(str,dc))==expected
    family=M(report['actual_binary_family_tangent']);assert A*T*family==s.zeros(135,3) and family[report['binary_parameters'],:].det()!=0
    z,u,bf=primitive_beta(model);a,c,d=s.symbols('a c d');v1=s.Matrix([a,1,c,0]);v2=s.Matrix([c,0,d,1]);bp=bf.subs(dict(zip(u,point[one,:])))
    section=s.Matrix([(v1.T*bp.subs(dict(zip(z,v)))*v2)[0] for v in [v1,v2]])
    CJ=section.jacobian([a,c,d]).subs({a:0,c:0,d:0});assert CJ==M(report['curve_jacobian']) and CJ.rank()==2
    vec=s.Matrix([bf[i,j].diff(z[k]) for (i,j),k in COORDS]);basis=vec.jacobian(u);br=list(basis.T.rref()[1]);inv=basis[br,:].inv()
    action=s.Matrix.hstack(*[inv*infinitesimal_action(vec,B)[br,:] for B in sp_basis()]);ar=action.subs(dict(zip(u,beta)))
    x,y,zz,w=X[1:];spatial=s.Matrix([x,y,zz,w]);Q=s.zeros(4);Q[1,3]=zz*zz+105*x*x;Q-=Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q
        return s.Matrix(model.normalize(P)['parameters'])[two,:]
    assert coords(Q)==q
    O=[]
    for B in sp_basis():
        dQ=sum(((B*spatial)[i]*Q.diff(spatial[i]) for i in range(4)),s.zeros(4))-B*Q-Q*B.T
        dq=coords(dQ);dq-=q*dq[fixed]/q[fixed];O.append(dq)
    O=s.Matrix.hstack(*O);stab=M(geo['stabilizer']);assert O.rank()==5 and stab.rank()==5 and O*stab==s.zeros(45,5)
    H=beta.jacobian(parameters);ann=M(geo['source_annihilator']);assert ann*H==s.zeros(7,9) and ann.rank()==7
    normal=(ann*ar*stab).applyfunc(s.expand);assert (normal-M(geo['normal_stabilizer_matrix'])).applyfunc(s.expand)==s.zeros(7,5) and normal.to_DM().rank()==2
    D=H.row_join(ar);joint=D.col_join(s.zeros(45,9).row_join(O));assert D.subs(sample).rank()==16 and joint.subs(sample).rank()==16
    assert D.subs(sample)[:,geo['beta_minor_columns']].det()==s.Rational(geo['beta_minor'])
    paths=[Path(__file__),rp,ROOT/'code/s61_minimal_pure_binary_endpoint_20261101.py',ROOT/'research/notes/S61_MINIMAL_PURE_BINARY_ENDPOINT_20261101.md',ROOT/'research/notes/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.md',ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/truncated_local_ideal_20261101.cpp',ROOT/'code/truncated_local_ideal_io_20261101.py']
    return dict(id='S61_MINIMAL_PURE_BINARY_ENDPOINT_REVIEW_20261101',all_passed=True,all_input_hashes_verified=True,raw135_all_H9_actual_line_verified=True,closed_embedding_coordinate='delta[35]=1/2',fixed_chart_value='u38=105',raw135_point_ranks=list(ranks),actual_line_tangent_verified=True,actual_binary_smooth_chart_verified=True,projective_orbit_rank=5,normal_stabilizer_rank=2,beta_rank=16,joint_rank=16,
                conceptual_review=dict(strict_filtered_kernel_torsion_free_accepted=True,DVR_torsion_free_argument_accepted=True,origin_factor_completion_accepted=True,relative_ambient_t2_and_generic16D_removal_accepted=True,flag_degenerations_in_same_proper_closed_locus_accepted=True),
                scope=['Independent QQ raw135 line, point Jacobian and quadratic obstructions, actual smooth binary chart, and source geometry.','Primary five-jet finite-field and Macaulay2 calculations reproduced separately; no independent implementation of that elimination claimed.','DVR argument reuses reviewed semisimple endpoint B with the present hypotheses independently checked.','No full freeness or affine exhaustion assertion.'],seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
