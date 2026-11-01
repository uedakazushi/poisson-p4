#!/usr/bin/env python3
"""Two dominant nonzero-modular minimal-nilpotent families in the unsaturated infinity locus.

All 135 leading identities are exact. This does not exhaust minimal-nilpotent infinity or
identify the saturation by the homogenizing coordinate.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61,X,J as J5
from s61_known_fiber_audit import equation_jacobian
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,J,infinitesimal_action,minor_witness,sp_basis
ROOT=Path(__file__).resolve().parents[1]
def sm(M):return [[str(v) for v in row] for row in M.tolist()]
def run(mode="M0"):
    assert mode in ("M0","M1")
    m=S61();weights=m.reduced['parameter_weights'];one=[i for i,t in enumerate(weights) if t==1];two=[i for i,t in enumerate(weights) if t==2]
    allu=s.symbols('u1:62');u=[allu[i] for i in one];x,y,z,w=X[1:];g=s.symbols('g0:9')
    mons=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3];G=sum(a*b for a,b in zip(g,mons,strict=True))
    Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-s.diff(G,w);Q[1,3]=s.diff(G,z)
    if mode=="M1":Q[0,1]=x*x/2;Q[1,2]-=x*z;Q[2,3]=x*x/2
    Q=Q-Q.T;B=X[0]**2*J5;B[1:,1:]+=Q;normalized=m.normalize(B)['parameters'];assert all(normalized[i]==0 for i in one)
    values=[allu[i] if i in one else normalized[i] for i in range(61)]
    extra_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';extra=json.loads(extra_path.read_text())['universal_new_equation']
    top=[[term for term in f if sum(i in two for i in term[0])==2] for f in m.reduced['equations']+[extra]]
    equations=[s.expand(sum(s.Rational(n,d)*s.prod(values[i] for i in mon) for mon,n,d in f)) for f in top]
    assert equations[:70]==[0]*70
    k,K,L=s.symbols('k K L');beta_sub={allu[3]:0,allu[11]:0,allu[4]:0,allu[12]:0,allu[13]:k,allu[48]:2*k+1,allu[49]:K-5*allu[34],allu[50]:5*allu[25]-L}
    a,b,c=s.symbols('a b c');shape=dict(zip(g,[a,b,0,c,0,0,0,0,0],strict=True))
    restricted=[s.factor(f.subs(beta_sub).subs(shape)) for f in equations]
    assert [i for i,f in enumerate(restricted) if f!=0]==[70,84,87]
    solution={a:-(2*K+18*(2*k+1)*allu[34])/(3*(2*k+1)*(9*k+7)),b:-(2*L-18*(2*k+1)*allu[25])/(3*(2*k+1)*(9*k+2)),c:-2/(3*(2*k+1))}
    if mode=="M1":solution={a:(-(3*k+1)*K+12*(2*k+1)*allu[34])/(3*(2*k+1)*(9*k+7)),b:((3*k+2)*L-12*(2*k+1)*allu[25])/(3*(2*k+1)*(9*k+2)),c:-(3*k+1)/(3*(2*k+1))}
    assert all(s.cancel(f.subs(solution))==0 for f in restricted)
    family_g={key:s.sympify(v).subs(solution) for key,v in shape.items()}
    assert all(s.cancel(f.subs(beta_sub).subs(family_g))==0 for f in equations)
    full_vertex=[s.expand(f.subs(beta_sub)) for f in equations[:134]]
    field=s.QQ.frac_field(k,K,L,allu[25],allu[34])
    vertex_gb=s.groebner(full_vertex,*g,order='grevlex',domain=field)
    expected=[s.cancel(gi-family_g[gi]) for gi in g]
    assert len(vertex_gb.polys)==9
    assert all(s.cancel(f-h)==0 for f,h in zip(vertex_gb,expected,strict=True))
    assert s.expand(equations[-1].subs(beta_sub).subs(shape))==0
    mod=s.Matrix([sum(s.diff(Q[i,j],X[j+1]) for j in range(4)) for i in range(4)])
    A=mod.jacobian(list(X[1:]));assert A==s.Matrix([[0,0,0,0],[-2,0,0,0],[0,0,0,0],[0,0,0,0]])
    assert A.rank()==1 and A*A==s.zeros(4) and A.T*J+J*A==s.zeros(4)
    sample_beta=dict.fromkeys(u,s.S.Zero);sample_beta.update({allu[13]:s.Integer(1),allu[48]:s.Integer(3),allu[25]:s.Integer(2),allu[34]:s.Integer(3),allu[49]:s.Integer(4),allu[50]:s.Integer(5)})
    sample_parameters={k:1,K:19,L:5,allu[25]:2,allu[34]:3}
    sample_g={key:s.cancel(v.subs(sample_parameters)) for key,v in family_g.items()}
    expected_sample=[-s.Rational(25,18),s.Rational(98,99),-s.Rational(2,9)] if mode=="M0" else [s.Rational(2,9),-s.Rational(47,99),-s.Rational(4,9)]
    assert [sample_g[g[i]] for i in [0,1,3]]==expected_sample
    point=[s.sympify(v).subs(sample_beta).subs(sample_g) for v in values];assert all(m.evaluate(f,point)==0 for f in top)
    fixed=s.Matrix(equations).jacobian(g).subs(sample_beta).subs(sample_g);fixed_minor=minor_witness(fixed);assert fixed_minor['rank']==9
    full=equation_jacobian(top,point)[:,two];full_minor=minor_witness(full)
    zz,uu,beta=primitive_beta(m);assert uu==u
    vector=s.Matrix([s.diff(beta[i,j],zz[h]) for (i,j),h in COORDS]);basis=vector.jacobian(u);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();beta_point=vector.subs(sample_beta)
    action=s.Matrix.hstack(*[left*infinitesimal_action(beta_point,aa)[rows,:] for aa in sp_basis()])
    excluded={allu[i] for i in [3,11,4,12,48]};columns=[];variables=[]
    for ui in u:
        if ui in excluded:continue
        column=s.zeros(16,1);column[u.index(ui)]=1
        if ui==allu[13]:column[u.index(allu[48])]=2
        columns.append(column);variables.append(str(ui))
    dominance=s.Matrix.hstack(*columns,action);dm=minor_witness(dominance);assert dm['rank']==16 and s.Rational(dm['determinant'])==48
    vv=s.Matrix([0,1,0,0]);line_action=s.Matrix.hstack(*[(-aa*vv)[[0,2,3],:] for aa in sp_basis()]);incidence=dominance.col_join(s.zeros(3,11).row_join(line_action));im=minor_witness(incidence)
    assert im['rank']==17 and s.Rational(im['determinant'])==20
    inputs=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_known_fiber_audit.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',m.source/'symplectic_slice_reduced_equations.json',extra_path]
    return {'id':'S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101','field':'QQ','all_passed':True,'minimal_component':mode,'G':str(G),'Q':sm(Q),'g_variables':list(map(str,g)),'beta_variables':list(map(str,u)),
      'weight_one_u_indices_zero_based':one,'weight_two_u_indices_zero_based':two,'q_coordinates':list(map(str,[normalized[i] for i in two])),
      'beta_substitution':{str(key):str(v) for key,v in beta_sub.items()},'shape':{str(key):str(v) for key,v in shape.items()},
      'nonzero_restricted_equations':{str(i):str(f) for i,f in enumerate(restricted) if f!=0},'solution':{str(key):str(v) for key,v in solution.items()},
      'generic_fixed_vertex_groebner_basis':list(map(str,vertex_gb)),'generic_fixed_vertex_field':str(field),'generic_fixed_vertex_equations':list(map(str,full_vertex[70:])),
      'denominator_open':'(2*k+1)*(9*k+7)*(9*k+2) != 0','all_135_leading_identities':True,'G6_zero_on_shape_for_vertex_beta':True,
      'modular_matrix':sm(A),'modular_partition':[2,1,1],'sample_beta':{str(key):str(v) for key,v in sample_beta.items()},'sample_g':{str(key):str(v) for key,v in sample_g.items()},
      'sample_S61_coordinates':list(map(str,point)),'fixed_v_jacobian_minor':fixed_minor,'full_infinity_jacobian_minor':full_minor,'full_projective_tangent_dimension':44-full_minor['rank'],
      'family_variables':variables,'dominance_matrix':sm(dominance),'dominance_minor':dm,'beta_line_incidence_matrix':sm(incidence),'beta_line_incidence_minor':im,
      'scope':'A nonzero-modular curve of unsaturated infinity points over sufficiently general beta. No complete M0 exhaustion, saturation, affine extension, or full S61 upper bound.',
      'proof_note':'research/notes/S61_MINIMAL_NILPOTENT_INFINITY_FAMILY_20261101.md','input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args();r=run();r['M1_family']=run('M1');args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'full_infinity_rank':r['full_infinity_jacobian_minor']['rank'],'fixed_v_rank':9,'dominance_rank':16,'joint_rank':17,'output':str(args.output)}))
if __name__=='__main__':main()
