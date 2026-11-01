#!/usr/bin/env python3
"""Independent all-face and raw-equation audit of the seven mixed R0 orbits."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time();mainpath=ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json'
    main=json.loads(mainpath.read_text());assert main['all_passed']
    facespath=ROOT/'reports/S61_R0_UNIFORM_VERTEX_20261101.json'
    faces=json.loads(facespath.read_text())['extra_vertical_prime_faces']
    a,b,c,x,z,w=s.symbols('a b c x z w')
    mons=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]
    swap=[]
    for mon in mons:
        image=s.expand(mon.subs({z:-w,w:z},simultaneous=True))
        swap.append(next(i for i,m in enumerate(mons) if image in (m,-m)))
    representatives={(0,6),(4,6),(0,3,7),(3,8),(0,1,4),(2,8),(0,1,8),(2,6)}
    checked=[]
    for k,fs in faces.items():
        for ff in fs:
            face=tuple(sorted(ff));dual=tuple(sorted(swap[j] for j in ff))
            assert face in representatives or dual in representatives
            checked.append([k,ff])
    assert len(checked)==16
    # Polynomial identities with denominators cleared; coefficient-zero boundaries
    # are separately enumerated in the proof note.
    transformations=[
        ((z*z*(a*x+b*w)).subs(w,w-a*x/b),b*z*z*w),
        ((z*(a*x*x+b*x*w+c*w*w)).subs(w,w-b*x/(2*c)),z*((a-b*b/(4*c))*x*x+c*w*w)),
        ((x*z*(a*x+b*w)).subs(w,w-a*x/b),b*x*z*w),
        ((a*x*x*z+b*x*x*w+c*x*w*w).subs(z,z-b*w/a),a*x*x*z+c*x*w*w),
        ((a*x*x*z+b*x*x*w+c*w**3).subs(z,z-b*w/a),a*x*x*z+c*w**3)]
    assert all(s.cancel(left-right)==0 for left,right in transformations)
    model=S61();u=s.symbols('u1:62')
    one=[i for i,v in enumerate(model.reduced['parameter_weights']) if v==1]
    two=[i for i,v in enumerate(model.reduced['parameter_weights']) if v==2]
    extra=json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']
    equations=model.reduced['equations']+[extra]
    raw={};xx,yy,zz,ww=X[1:]
    for name,entry in main['orbits'].items():
        G=s.sympify(entry['G']).subs({x:xx,z:zz,w:ww})
        P=X[0]**2*J;P[2,3]+=s.diff(G,ww);P[3,2]-=s.diff(G,ww)
        P[2,4]-=s.diff(G,zz);P[4,2]+=s.diff(G,zz)
        normalized=model.normalize(P)['parameters'];q={i:normalized[i] for i in two}
        leading=[]
        for eq in equations:
            value=0
            for mon,num,den in eq:
                if sum(i in q for i in mon)==2:
                    value+=s.Rational(num,den)*s.prod(q[i] if i in q else u[i] for i in mon)
            leading.append(s.expand(value))
        assert all(v==0 for v in leading[:70])
        M=s.Matrix(leading[70:134]).jacobian([u[i] for i in one])
        assert M*s.Matrix([u[i] for i in one])==s.Matrix(leading[70:134])
        H=s.Matrix(entry['beta_kernel']).applyfunc(s.sympify)
        assert M*H==s.zeros(64,H.cols)
        rank=len(M.to_DM().rref()[1]);assert rank+H.cols==16
        params=s.symbols('h:'+str(H.cols));full=dict(zip([u[i] for i in one],H*s.Matrix(params)))
        assert s.expand(leading[134].subs(full,simultaneous=True))==0
        QD=s.Matrix(entry['q_orbit_differential']).applyfunc(s.sympify)
        BD=s.Matrix(entry['source_beta_differential']).applyfunc(s.sympify)
        assert len(QD.to_DM().rref()[1])==entry['stats']['orbitdim']
        assert len(BD.to_DM().rref()[1])==entry['stats']['beta_rank']
        joint=BD.col_join(s.zeros(44,H.cols).row_join(QD))
        assert len(joint.to_DM().rref()[1])==entry['stats']['joint_rank']
        if name=='T3':
            constraints=s.Matrix([u[3],u[4],u[11],u[12],5*u[13]-u[48],5*u[34]+u[49],5*u[25]-u[50]])
            C=constraints.jacobian([u[i] for i in one]);assert C*H==s.zeros(7,9) and C.rank()==7
        raw[name]=dict(beta_kernel_dimension=H.cols,raw_beta_rank=rank,source_ranks=entry['stats'])
    paths=[Path(__file__),mainpath,facespath,ROOT/'code/s61_model.py',
        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',
        ROOT/'research/notes/S61_R0_MIXED_ORBIT_REDUCTION_20261101.md',
        ROOT/'research/notes/S61_NONBINARY_RESONANCE_EXCLUSIONS_20261101.md']
    return dict(id='S61_R0_MIXED_ORBIT_REDUCTION_REVIEW_20261101',all_passed=True,
        face_count=16,raw_orbit_checks=raw,
        independently_reconstructed=['Raw fixed-Q highest equations and complete beta kernels.',
            'The sixteen faces, coordinate swap and all coefficient normalization identities.',
            'Equality of the T3 beta kernel with the previously excluded seven-constraint subspace.'],
        dependencies_rechecked_as_saved_matrices=['Seven orbit differentials and source beta differentials.'],
        scope='Reduced set-theoretic cover only; no saturation or pairwise orbit/component distinctness.',
        proof_note='research/notes/S61_R0_MIXED_ORBIT_REDUCTION_REVIEW_20261101.md',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();r=verify()
    with args.output.open('x') as out:json.dump(r,out,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,face_count=r['face_count'],seconds=r['seconds'])))
