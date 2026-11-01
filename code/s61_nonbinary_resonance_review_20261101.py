#!/usr/bin/env python3
"""Independent local/matrix audit of the two nonbinary resonance nodes."""
import argparse, hashlib, json, time
from itertools import combinations_with_replacement
from pathlib import Path
import sympy as s
from s61_model import S61, X, J
from s61_minimal_formal_rings_review_20261101 import grouped, derivative, mat
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, infinitesimal_action, sp_basis
ROOT=Path(__file__).resolve().parents[1]

def run():
    start=time.time();m=S61();weights=m.reduced['parameter_weights'];one=[i for i,w in enumerate(weights) if w==1];two=[i for i,w in enumerate(weights) if w==2]
    path=ROOT/'reports/S61_NONBINARY_RESONANCE_NODES_20261101.json';source=json.loads(path.read_text());extra=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';eqs=m.reduced['equations']+[json.loads(extra.read_text())['universal_new_equation']]
    x,y,z,w=X[1:];V=s.Matrix([x,y,z,w]);local={}
    def qcoords(G):
        P=X[0]**2*J;P[2,3]+=s.diff(G,w);P[3,2]-=s.diff(G,w);P[2,4]-=s.diff(G,z);P[4,2]+=s.diff(G,z)
        pp=m.normalize(P)['parameters'];return s.Matrix([pp[i] for i in two])
    for label,record in source['cases'].items():
        rec=record['local_data'];point=list(map(s.Rational,rec['point61']));q=s.Matrix([point[i] for i in two]);columns=rec['chart_columns'];_,polys=grouped(eqs,weights,{i:point[i] for i in one});A,tc=derivative(polys,q,columns);F=A.row_join(tc)
        assert F==mat(rec['jacobian']);T=mat(rec['tangent_basis']);C=mat(rec['cokernel_rows']);assert F.rank()==42 and T.rank()==3 and F*T==s.zeros(135,3) and C*F==s.zeros(C.rows,45)
        n=3;mons=list(combinations_with_replacement(range(n),2));mi={a:i for i,a in enumerate(mons)};forms={a:list(T[j,:]) for j,a in enumerate(columns)};forms[next(i for i in range(45) if i not in columns)]=[0]*3;tf=list(T[-1,:]);H=s.zeros(135,len(mons))
        for row,f in enumerate(polys):
            for mon,v in f.items():
                first,second=[forms[i] for i in mon]+[tf]*(2-len(mon))
                for i in range(n):
                    for j in range(n):H[row,mi[tuple(sorted((i,j)))]]+=v*first[i]*second[j]
        ob=C*H;assert ob==mat(rec['quadratic_obstruction_basis']) and ob.rank()==2
        a,b=s.symbols('a b');v=[a,b,-a];qs=[s.expand(sum(r[i]*v[j]*v[k] for i,(j,k) in enumerate(mons))) for r in ob.tolist()];assert s.gcd(*qs)==1
        gb=s.groebner(qs,a,b,order='grevlex');lm=[p.monoms(order=gb.order)[0] for p in gb.polys];std=[(i,j) for i in range(5) for j in range(5) if not any(i>=u and j>=v for u,v in lm)];assert len(std)==4
        vm=mat(record['main_tangent']);vt=mat(record['tail_tangent']);assert s.Matrix.hstack(vm,vt).rank()==2 and vm[0]+vm[2]!=0 and vt[0]+vt[2]!=0
        qt=qcoords(z**3 if label=='ratio3/5' else w**3);assert T*vt==qt[columns,:].col_join(s.zeros(1,1));c=s.Symbol('c');qq=q+c*qt
        assert all(s.expand(sum(v*s.prod(qq[i] for i in mon) for mon,v in f.items() if len(mon)==2))==0 for f in polys)
        local[label]={'homogenized_rank':42,'tangent_dimension':3,'obstruction_rank':2,'transverse_length':4,'independent_branch_tangents':True,'entire_explicit_tail_verified':True}
    # Reconstruct the full universal leading Jacobian, using the stored beta
    # basis only as an input whose dimension and equation annihilation are checked.
    tail=source['tail_universal_rank'];K=mat(tail['beta_kernel']);assert K.rank()==9;bb=s.symbols('b0:9');bv=K*s.Matrix(bb);_,polys=grouped(eqs,weights,dict(zip(one,bv)),minimum=2);q=qcoords(x*z*w+z**3);assert q==mat(tail['q_point']);D,_=derivative(polys,q,list(range(45)));C=mat(tail['constant_kernel']);N=mat(tail['residual_kernel']);assert D[:70,:].rank()==28 and C.rank()==17 and D[:70,:]*C==s.zeros(70,17)
    R=(D[70:,:]*C).applyfunc(s.expand);assert R==mat(tail['residual_matrix']);assert N[tail['free'],:]==s.eye(3);assert (R.to_DM()*N.to_DM()).is_zero_matrix
    # Independently recover beta-action columns and the fixed-Q kernel. This
    # checks that the saved source differential genuinely supplies the main branch.
    zz,bu,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();action=s.Matrix.hstack(*[left*infinitesimal_action(vec,A)[rows,:] for A in sp_basis()]);bvar=s.symbols('p0:16');q=qcoords(x*z*w);_,polys=grouped(eqs[:134],weights,dict(zip(one,bvar)),minimum=2);vals=s.Matrix([s.expand(sum(v*s.prod(q[i] for i in mon) for mon,v in f.items())) for f in polys]);M=vals.jacobian(bvar);HK=s.Matrix.hstack(*M.nullspace());assert HK.cols==10
    Q=s.zeros(4);Q[1,2]=x*z;Q[2,1]=-x*z;Q[1,3]=-x*w;Q[3,1]=x*w;fixed=next(i for i,v in enumerate(q) if v);chart=[i for i in range(45) if i!=fixed];orbit_cols=[]
    for A in sp_basis():
        dQ=sum(((A*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-A*Q-Q*A.T
        P=X[0]**2*J;P[1:,1:]+=dQ;pp=m.normalize(P)['parameters'];dq=s.Matrix([pp[i] for i in two]);dq-=q*dq[fixed]/q[fixed];orbit_cols.append(dq[chart,:])
    orbit=s.Matrix.hstack(*orbit_cols);stab=s.Matrix.hstack(*orbit.nullspace());assert orbit.rank()==7 and stab.cols==3
    for label,record in source['cases'].items():
        rec=record['local_data'];p=list(map(s.Rational,rec['point61']));BD=mat(record['beta_parameter_differential']);assert BD[:,:10]==HK;assert BD[:,10:]==action.subs(dict(zip(bu,[p[i] for i in one])));assert BD.rank()==16
        local[label]['main_source_beta_rank']=16
        HH=mat(record['restricted_beta_kernel']);assert HH.rank()==9 and M*HH==s.zeros(134,9);constraint=s.zeros(1,16);constraint[0,one.index(13)]=1;constraint[0,one.index(48)]=-p[13]/p[48];assert constraint*HH==s.zeros(1,9)
        restricted=HH.row_join(BD[:,10:]);assert restricted==mat(record['restricted_beta_differential']);cols=record['node_beta_minor_columns'];det=restricted[:,cols].det(method='domain-ge');assert det==s.Rational(record['node_beta_minor'])!=0
        assert restricted.col_join(s.zeros(44,9).row_join(orbit)).rank()==16
        bh=HH*s.Matrix(s.symbols('h0:9'));stabaction=(action*stab).subs(dict(zip(bu,bh)),simultaneous=True);assert (constraint*stabaction).applyfunc(s.expand)==s.zeros(1,3) and (M*stabaction).applyfunc(s.expand)==s.zeros(134,3)
        local[label].update(restricted_beta_rank=16,joint_incidence_rank=16,projective_Q_orbit_dimension=7,connected_stabilizer_preserves_Hres=True)
    paths=[Path(__file__),path,extra,ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_formal_rings_review_20261101.py',ROOT/'research/notes/S61_NONBINARY_RESONANCE_NODES_20261101.md',ROOT/'reports/S61_NONBINARY_FORMAL_RING_20261101.json']
    return dict(id='S61_NONBINARY_RESONANCE_REVIEW_20261101',all_passed=True,local=local,universal_tail=dict(first_rank=28,kernel_dimension=17,residual_kernel_dimension=3,leading_rank_upper_bound=42),proof_accepted=True,scope='The two specified local rings and generic-beta removal of the two sixteen-dimensional node-incidence closures; no full infinity exhaustion or affine finiteness claim.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
if __name__=='__main__':main()
