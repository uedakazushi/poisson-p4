#!/usr/bin/env python3
"""Independent rational geometry and raw135 point audit of semisimple endpoint B."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
M=lambda a:s.Matrix(a).applyfunc(s.sympify)

def run():
    start=time.time();rp=ROOT/'reports/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.json';pp=ROOT/'reports/S61_BINARY_PRIME_SPREADING_20261101.json'
    report=json.loads(rp.read_text());prime=json.loads(pp.read_text())
    for r in [report,prime]:
        assert r['all_passed']
        for name,digest in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    model=S61();one=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(model.reduced['parameter_weights']) if w==2];qi={i:j for j,i in enumerate(two)};ui={i:j for j,i in enumerate(one)}
    equations=model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    encoded=[[(tuple(qi[i] for i in mon if i in qi),tuple(ui[i] for i in mon if i in ui),s.Rational(n,d)) for mon,n,d in eq] for eq in equations]
    geo=report['source_geometry'];beta=M(geo['beta']);parameters=list(map(s.Symbol,geo['parameters']));delta=M(geo['actual_line_delta']);fixed=two.index(29);q=s.zeros(45,1);q[fixed]=1
    products={}
    for row,eq in enumerate(encoded):
        coefficients=[0,0,0]
        for qm,bm,c in eq:
            if len(qm)!=2:continue
            i,j=qm;values=[q[i]*q[j],q[i]*delta[j]+delta[i]*q[j],delta[i]*delta[j]]
            if not any(values):continue
            if bm not in products:products[bm]=s.prod(beta[k] for k in bm)
            for k,value in enumerate(values):coefficients[k]+=c*products[bm]*value
        assert all(s.cancel(value)==0 for value in coefficients),row
    sample={s.Symbol(k):s.Rational(v) for k,v in geo['beta_sample'].items()};pbeta=beta.subs(sample);point=M(report['point_local_data']['point61']);assert pbeta==point[one,:]
    chart=[i for i in range(45) if i!=fixed];A=s.zeros(135,45);cache={}
    for row,eq in enumerate(encoded):
        for qm,bm,c in eq:
            if bm not in cache:cache[bm]=s.prod(pbeta[k] for k in bm)
            co=c*cache[bm]
            if len(qm)==2:
                i,j=qm
                if i!=fixed:A[row,chart.index(i)]+=co*q[j]
                if j!=fixed:A[row,chart.index(j)]+=co*q[i]
            elif len(qm)==1:A[row,44]+=co*q[qm[0]]
    local=report['point_local_data'];T=M(local['tangent_basis']);assert A==M(local['jacobian']) and A.rank()==35 and A[:,:44].rank()==34
    assert A*T==s.zeros(135,10) and T[-1,:]==s.zeros(1,10)
    piv=list(A.rref()[1]);free=[i for i in range(45) if i not in piv];assert T[free,:]==s.eye(10)
    tangent=delta.subs(sample)[chart,:].col_join(s.zeros(1,1));dc=tangent[free,:]
    assert T*dc==tangent and list(map(str,dc))==report['actual_line_tangent_coordinates']
    family=M(report['actual_binary_family_tangent']);assert A*T*family==s.zeros(135,3) and family[report['binary_parameters'],:].det()!=0
    z,u,bf=primitive_beta(model);vec=s.Matrix([bf[i,j].diff(z[k]) for (i,j),k in COORDS]);basis=vec.jacobian(u);br=list(basis.T.rref()[1]);inv=basis[br,:].inv()
    action=s.Matrix.hstack(*[inv*infinitesimal_action(vec,B)[br,:] for B in sp_basis()]);ar=action.subs(dict(zip(u,beta)))
    x,y,zz,w=X[1:];spatial=s.Matrix([x,y,zz,w]);Q=s.zeros(4);Q[1,2]=x*w;Q-=Q.T
    def coords(Q):
        P=X[0]**2*J;P[1:,1:]+=Q
        return s.Matrix(model.normalize(P)['parameters'])[two,:]
    assert coords(Q)==q
    O=[]
    for B in sp_basis():
        dQ=sum(((B*spatial)[i]*Q.diff(spatial[i]) for i in range(4)),s.zeros(4))-B*Q-Q*B.T
        dq=coords(dQ);dq-=q*dq[fixed];O.append(dq)
    O=s.Matrix.hstack(*O);stab=M(geo['stabilizer']);assert O.rank()==5 and stab.rank()==5 and O*stab==s.zeros(45,5)
    H=beta.jacobian(parameters);ann=M(geo['source_annihilator']);assert ann*H==s.zeros(6,10) and ann.rank()==6
    normal=(ann*ar*stab).applyfunc(s.expand);assert (normal-M(geo['normal_stabilizer_matrix'])).applyfunc(s.expand)==s.zeros(6,5) and normal.to_DM().rank()==1
    D=H.row_join(ar);joint=D.col_join(s.zeros(45,10).row_join(O));assert D.subs(sample).rank()==16 and joint.subs(sample).rank()==16
    assert D.subs(sample)[:,geo['beta_minor_columns']].det()==s.Rational(geo['beta_minor'])
    paths=[Path(__file__),rp,pp,ROOT/'research/notes/S61_MINIMAL_SEMI_B_ENDPOINT_20261101.md',ROOT/'research/notes/S61_BINARY_PRIME_SPREADING_20261101.md',ROOT/'code/s61_minimal_semi_B_endpoint_20261101.py',ROOT/'code/s61_binary_prime_spreading_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'code/truncated_local_ideal_20261101.cpp',ROOT/'code/truncated_local_ideal_io_20261101.py',ROOT/'research/notes/S61_BINARY_MIXED_STRATUM_REMOVAL_20261101.md']
    return dict(id='S61_MINIMAL_SEMI_B_ENDPOINT_REVIEW_20261101',all_passed=True,all_input_hashes_verified=True,raw135_all_H_actual_line_verified=True,raw135_point_jacobian_verified=True,leading_rank=34,full_rank=35,homogenizer_linear_part_zero=True,actual_line_tangent_verified=True,actual_binary_family_tangent_rank=3,projective_orbit_rank=5,normal_stabilizer_rank=1,beta_rank=16,joint_rank=16,
                conceptual_review=dict(strict_filtered_initial_kernel_argument_accepted=True,C_primary_torsion_free_argument_accepted=True,good_reduction_scope_accepted=True,DVR_equal_rank_specialization_accepted=True,p_flat_line_Tor_identification_accepted=True,generic_Hilbert_and_homogenizer_order_accepted=True,completion_and_QQ_descent_accepted=True,relative_ambient_and_generic16D_removal_accepted=True),
                scope=['Independent raw135 rational line identities, point Jacobian and full source geometry.','The primary finite-field jet producer and M2 primary decomposition were reproduced separately, not replaced by an independent finite-field implementation.','The good-reduction and DVR arguments receive proof/source review, not an inference from a finite-field sample.','No endpoint A or full affine exhaustion assertion.'],seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
