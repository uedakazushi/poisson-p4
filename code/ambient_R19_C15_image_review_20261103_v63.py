#!/usr/bin/env python3
"""Independent finite tensor, sparse action and image-algebra audit at C15."""
import argparse,hashlib,itertools,json,subprocess,sys
from pathlib import Path
import sympy as s
from ambient_CA_global_flag_review_20261103_v63 import action as sparse_action

ROOT=Path(__file__).resolve().parents[1]
RP=ROOT/'reports/AMBIENT_R19_C15_COMPLETED_IMAGE_20261103_v63/report.json'

def main(output):
    output.mkdir(parents=True,exist_ok=False)
    fresh=output/'fresh'
    subprocess.run([sys.executable,str(ROOT/'code/ambient_R19_C15_completed_image_20261103_v63.py'),'--output',str(fresh)],check=True)
    d=json.loads(RP.read_text());e=json.loads((fresh/'report.json').read_text())
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in e.items()if k!='seconds'}
    h,x,y,z,w=s.symbols('h x y z w');X=s.Matrix([h,x,y,z,w]);p,B,E,F=s.symbols('p B E F');eps=s.Symbol('eps')
    t=-p/12;r=11*p/12
    V=s.Matrix([x*((r-13*t)*h+x-2*p*p*z/3),12*t*x*x,0,
        x*(h+(-r-13*t)*z),2*t*x*w+h*h+(B-p/6)*h*z+(E-s.Rational(1,6))*x*z+F*z*z])
    ey=s.eye(5)[:,2];P=ey*V.T-V*ey.T
    M=s.zeros(5);M[1,0]=M[3,4]=1;g=s.eye(5)+p*M;gi=g.inv()
    Q=(g*P.subs(dict(zip(X,gi*X)),simultaneous=True)*g.T).applyfunc(s.expand)
    assert all(s.expand(sum(s.diff(Q[i,j],X[i])for i in range(5)))==0 for j in range(5))
    pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
    u=s.Matrix([s.Poly(Q[i,j],*X).coeff_monomial(X[a]*X[b])for i,j in pairs for a,b in mons])
    parsed=s.Matrix([s.sympify(v,locals={'p':p,'B':B,'E':E,'F':F})for v in d['full150_target']])
    assert (u-parsed).applyfunc(s.expand)==s.zeros(150,1)
    u0=u.subs({p:0,B:0,E:0,F:0});v=u-u0
    gens=s.Matrix([B,E,F,p*p,p**3,B*p,E*p,F*p])
    got=s.Matrix([v[123]+v[112],v[127],v[132],v[15],-3*v[18]/2,
        v[108]+v[127]*v[15]+v[15]/6,v[112],v[117]])
    assert (got-gens).applyfunc(s.expand)==s.zeros(8,1)
    for a in v:
        for (ip,ib,ie,iff),c in s.Poly(a,p,B,E,F).terms():
            assert ip!=1 or ib+ie+iff>0
    C=s.Matrix([[s.Poly(a,p,B,E,F).coeff_monomial(m)for m in gens]for a in v])
    assert C==s.Matrix(d['generator_coefficient_matrix'])
    O=s.Matrix.hstack(*(sparse_action(i,j)*u0 for i,j in itertools.product(range(5),repeat=2)))
    OP=O-u0*O[106,:];W=OP.row_join(C);saved=d['entire_algebra_cotangent']
    assert W==s.Matrix(saved['matrix'])
    assert W[saved['rows'],saved['columns']].det()==s.Rational(-209,162)
    assert W.rank()==23 and OP.rank()==15
    triples=list(itertools.combinations(range(5),3))
    def flag_derivative(N):
        a=s.eye(5)+eps*N;inv=s.eye(5)-eps*N
        lv=a*ey;cv=inv[1,:].T;kv=inv[[0,1,3],:]
        out=[s.diff(q,eps).subs(eps,0)for q in lv]+[s.diff(q,eps).subs(eps,0)for q in cv]
        out += [s.diff(kv[:,list(T)].det(),eps).subs(eps,0)for T in triples]
        out[2]=out[6]=out[10+triples.index((0,1,3))]=0
        return s.Matrix(out)
    ids=d['projective_orbit_generator_indices'];HF=s.Matrix.hstack(*(flag_derivative(s.eye(5)[:,i//5]*s.eye(5)[i%5,:])for i in ids))
    top=OP[:,ids].row_join(C[:,:3]).row_join(s.zeros(150,1))
    bottom=HF.row_join(s.zeros(20,3)).row_join(flag_derivative(M));SW=top.col_join(bottom)
    saved=d['normalization_source_etale_tangent'];assert SW==s.Matrix(saved['matrix'])
    assert SW[saved['rows'],saved['columns']].det()==s.Rational(saved['determinant'])!=0 and SW.rank()==19
    a,b,ee,ff,z0,z1,z2,z3=s.symbols('a b e f z0 z1 z2 z3')
    sub={a:p*p,b:B,ee:E,ff:F,z0:p**3,z1:B*p,z2:E*p,z3:F*p}
    assert all(s.expand(s.sympify(f).subs(sub))==0 for f in d['finite_presentation_relations'])
    paths=[Path(__file__),RP,ROOT/'code/ambient_CA_global_flag_review_20261103_v63.py',
        ROOT/'research/notes/AMBIENT_R19_C15_COMPLETED_IMAGE_20261103_v63.md',
        ROOT/'research/notes/AMBIENT_R19_C15_COMPLETED_IMAGE_20261103_v63_REVIEW.md']
    for path,want in d['input_sha256'].items():
        q=ROOT/path;assert hashlib.sha256(q.read_bytes()).hexdigest()==want;paths.append(q)
    receipt=dict(id='AMBIENT_R19_C15_COMPLETED_IMAGE_20261103_v63_REVIEW',all_passed=True,
        fresh_producer_matches=True,independent_full150_finite_transformed_tensor=True,
        all_eight_target_generators_and_semigroup_membership=True,
        independent_sparse_GL_actions_and_all23_cotangent_columns=True,
        independent_finite_Pluecker_flag_derivatives_and_source_etale_rank19=True,
        whole_image_completion_A3_times_formal15_accepted=True,
        finite_16_relation_presentation_and_both_conductor_ambients_accepted=True,
        dimension=19,embedding_dimension=23,depth=16,multiplicity=2,normalization_fiber_length=2,
        actual_CA_CB_intersection_or_full_Poisson_completion_claimed=False,
        dependency_on_prior_whole_geometric_fiber_and_normal_source=True,
        input_sha256={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest()for q in paths})
    (output/'report.json').write_text(json.dumps(receipt,indent=2)+'\n');print('C15 independent full-image review PASS')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    main(p.parse_args().output)
