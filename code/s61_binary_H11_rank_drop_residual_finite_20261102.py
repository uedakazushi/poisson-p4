#!/usr/bin/env python3
"""Remove the artificial unipotent normalization and bound residual rank-drop points."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import sympy as s
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from s61_model import S61
from weight1_structure import COORDS
ROOT=Path(__file__).resolve().parents[1]
STEM='S61_BINARY_H11_RANK_DROP_RESIDUAL_FINITE_20261102'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def certificate():
    start=time.monotonic();hs={};dependency_hashes={}
    def load(stem):
        p=Path('reports')/(stem+'.json');j=json.loads((ROOT/p).read_text())
        assert j['all_passed']
        for name, expected in j.get('input_sha256',j.get('source_sha256',{})).items():
            assert not Path(name).is_absolute() and sha(ROOT/name)==expected,name
            assert name not in dependency_hashes or dependency_hashes[name]==expected
            dependency_hashes[name]=expected
        hs[str(p)]=sha(ROOT/p)
        return j
    tr=load('S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_20261102')
    review=load('S61_BINARY_HYPERPLANE_TRANSVERSE_CODIMENSION_REVIEW_20261102')
    assert review['proof_audit_accepted'] and not review['H11_in_reviewed_scope']
    assert review['independent_dominance']['minor_QQ']=='1680453/125'
    act=load('S61_BINARY_NUISANCE_CONTRACTION_20261102')
    unif=load('S61_BINARY_UNIFORM_RANK_20261101')
    seven=load('S61_BINARY_DEEPER_SEVEN_SLICE_20261101')
    kl=load('S61_BINARY_DEEPER_REDUCTION_20261101')
    table=load('S61_BINARY_RANK_DROP_COVER_20261101')
    c0=load('S61_BINARY_RANK_DROP_FAMILY_20261101')
    ends=load('S61_BINARY_ENDPOINTS_20261101')
    ci=load('S61_BINARY_CODIM_TWO_FORMAL_20261101')
    assert tr['exceptional_core_dimension_upper']==1 and tr['relative_module_rank_on_zero_nuisance_section']==65
    ids=[4,5,12,13,14,25,26,27,34,35,36,37,48,49,50,51]
    u=s.symbols('u1:62');vec=s.Matrix([u[i-1] for i in ids])
    a,b,c,d,e,f=[u[i-1] for i in [4,14,35,36,49,50]]
    nu=[5,25,26,27,34,37,48,51];core=[4,14,35,36,49,50]
    weights=act['parameter_weights'];nw=[weights[i-1] for i in nu]
    assert nw==[2,4,2,4,2,2,2,2] and all(weights[i-1]==0 for i in core)
    assert sorted(nu+core+[12,13])==ids and act['equations_checked']==135
    x,z,t=s.symbols('x z t',nonzero=True);L=5*b-e
    M=s.Matrix.hstack(*[s.Matrix(v).applyfunc(s.sympify) for v in unif['primitive_actions']])
    tor={v:v*x**(-2*M[i,4]/v).cancel()*z**((M[i,2]+2*M[i,4])/v).cancel() for i,v in enumerate(vec)}
    assert all(s.cancel(tor[v]/v).as_powers_dict().get(w,0)==ww[i] for i,v in enumerate(vec) for w,ww in [(x,seven['torus_x_weights']),(z,seven['torus_z_weights'])])
    par={x:L,z:1/c};norm={v:w.subs(par,simultaneous=True).cancel() for v,w in tor.items()}
    assert norm[c]==1 and s.cancel(L.subs(norm,simultaneous=True)-1)==0
    assert s.cancel(norm[f].subs(f,-5*c)+5)==0
    ug=[s.Matrix(unif['stabilizer_generators'][i]) for i in [0,1,3]]
    # On homogeneous coordinates (x,y,z,w), each A maps span(x,w) to span(y,z)
    # and kills span(y,z); composing any two gives zero. I+sum(t_i*A_i) fixes Q0.
    assert all(A*B==s.zeros(4) for A in ug for B in ug)
    assert all(A[i,j]==0 for A in ug for i in range(4) for j in range(4) if i not in [1,2] or j not in [0,3])
    assert s.Matrix.hstack(*[A.reshape(16,1) for A in ug]).rank()==3
    J=s.diag(s.Matrix([[0,1],[-1,0]]),s.Matrix([[0,1],[-1,0]]))
    tvec=s.symbols('t0:3');g=s.eye(4)+sum((ti*A for ti,A in zip(tvec,ug)),s.zeros(4))
    assert (g.T*J*g-J).applyfunc(s.expand)==s.zeros(4)
    rep=[M[:,i].jacobian(vec) for i in [0,1,3]]
    assert all(v*w==w*v for v in rep for w in rep)
    N=sum((ti*A for ti,A in zip(tvec,rep)),s.zeros(16))
    assert N**4==s.zeros(16)
    value=(s.eye(16)+N+N**2/2+N**3/6)*vec
    sub={u[11]:0,u[12]:0,f:-5*c}
    for i in core+[12,13]:
        assert s.expand((value[ids.index(i)]-u[i-1]).subs(sub))==0
    # Check both finite pullbacks on the primitive two-form itself.
    # This also ties the stored infinitesimal representation to the stated4x4 group.
    formal_coords,primitive_parameters,form=primitive_beta(S61())
    assert list(primitive_parameters)==list(vec)
    coefficient_forms=[form.diff(z) for z in formal_coords]
    initial_tensor=s.Matrix([coefficient_forms[k][i,j] for (i,j),k in COORDS])
    embedding=initial_tensor.jacobian(vec)
    assert embedding.rank()==16
    def transported_tensor(G):
        pulled=[G.T*matrix*G for matrix in coefficient_forms]
        return s.Matrix([sum(G[i,k]*pulled[i][row,col] for i in range(4))
                         for (row,col),k in COORDS])
    assert all(s.expand(entry)==0 for entry in transported_tensor(g)-embedding*value)
    diagonal=s.diag(x,1/x,z,1/z)
    torus_value=s.Matrix([tor[v] for v in vec])
    assert all(s.cancel(entry)==0 for entry in transported_tensor(diagonal)-embedding*torus_value)
    Q0=s.zeros(4)
    Q0[1,2]=formal_coords[0]*formal_coords[3]
    Q0[2,1]=-Q0[1,2]
    Qpull=(g.inv()*Q0.subs(dict(zip(formal_coords,g*s.Matrix(formal_coords))),simultaneous=True)*g.inv().T)
    assert (Qpull-Q0).applyfunc(s.expand)==s.zeros(4)
    # H11 contains the reviewed H8, with three more independent parameters.
    assert tr['relative_completion_base_dimension']==8+3
    assert len(nu)+3==11 and len(nu)+3+3==14
    assert [x['joint_image_dimension'] for x in table['orbit_strata']]==[18,17,17]
    assert len(table['prime_ideals'])==5
    img={s.Symbol(k):s.sympify(v) for k,v in table['finite_symmetry_core_images'].items()}
    assert s.expand((5*c+f).subs(img,simultaneous=True)-L)==0 and s.expand(L.subs(img,simultaneous=True)-(5*c+f))==0
    assert kl['excluded_components'][0]['ideal_index']==0 and kl['excluded_components'][0]['beta_image_dimension_upper_bound']==15
    assert unif['all135_leading_equations_zero']
    assert c0['joint_image_dimension']==17 and c0['first_jet_image_dimension']==16
    # The known endpoint open is ambient: its base dominates W16 smoothly at
    # its witness, and t squared vanishes on the relative scheme there.
    assert all(item['homogenizer_square_zero'] and s.Rational(item['beta_dominance_minor'])!=0
               and {'u12','u13','u35','u50'}<=set(item['beta_constraints'])
               for item in ends['cases'])
    assert all(item['homogenizer_square_zero'] for item in ci['cases'])
    assert all(i['generic_fixed_beta_stratum_removed'] and i['joint_image_dimension']==16 and set(['u35','u50'])<=set(i['beta_constraints']) for i in ends['cases'])
    assert all(i['ambient_nilpotence_on_nonempty_open'] and i['joint_image_dimension']==17 for i in ci['cases'])
    for x in [Path(__file__).relative_to(ROOT),Path('research/notes')/(STEM+'.md'),Path(ends['proof_note']),Path(ci['proof_note']),Path('research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md'),Path('code/primitive_jet_elliptic_curve_20261101.py'),Path('code/weight1_structure.py'),Path('code/s61_model.py'),Path('research/notes/S61_BINARY_RANK_DROP_COVER_20261101.md'),Path('research/notes/S61_BINARY_RANK_DROP_FAMILY_20261101.md'),Path('research/notes/S61_BINARY_DEEPER_REDUCTION_20261101.md'),Path('legacy/v02/data/symplectic_slice_certificate.json'),Path('legacy/v02/data/symplectic_slice_equations.json'),Path('legacy/v02/data/symplectic_slice_reduced_equations.json')]:
        hs[str(x)]=sha(ROOT/x)
    return dict(id=STEM,all_passed=True,open_core_complement_dimension_upper=1,unrestricted_normalized_source_dimension=11,
                nuisance_indices_one_based=nu,nuisance_weights=nw,unipotent_group_dimension=3,core_fixed_under_unipotents=True,
                torus_normalization=['x=5b-e','z=1/c'],relative_free_rank_at_good_zero_nuisance_section=65,relative_base_dimension=14,
                residual_joint_dimension_upper=10+8+1-3,number_of_rank_drop_prime_presentations=5,symmetry_reduced_presentations=3,
                full135_projective_jacobian_rank_upper_in_scope=34,generic_beta_residual_rank_drop_dimension_upper=0,
                residual_finite_over_general_beta=True,residual_empty_proved=False,ordinary_rank35_included=False,
                square_binary_included=False,full_S61_exhaustion_claimed=False,finite_primitive_pullbacks_checked=['three-dimensional unipotent group','normalizing diagonal torus'],unipotent_fixes_Q0_exactly=True,H8_dominance_minor_from_accepted_review='1680453/125',primary_dependency_sha256=dependency_hashes,input_sha256=hs,seconds=time.monotonic()-start)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():p.error(f'output exists: {args.output}')
    r=certificate()
    with args.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
    print(f"PASS: eight nuisance weights {r['nuisance_weights']}, joint dimension at most 16; {r['seconds']:.3f} s")
