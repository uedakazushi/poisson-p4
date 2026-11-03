#!/usr/bin/env python3
"""Independent all-Qa source/whole-image/actual-intersection transport audit.

The producers under review are never imported. Physical raw tensors, canonical
projection, sparse GL actions, Smith products and odd modules are rebuilt here.
Only the previously audited immutable S61 coordinate adapter is reused.
"""
import argparse
from collections import defaultdict
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from sympy.polys.matrices import DomainMatrix
from s61_model import S61, X as Y, pushforward, coefficients

ROOT=Path(__file__).resolve().parents[1]
X=s.symbols('h x y z w');h,x,y,z,w=X;a=s.Symbol('a')
PAIRS=list(itertools.combinations(range(5),2))
MONS=list(itertools.combinations_with_replacement(range(5),2))
BASIS=list(itertools.product(PAIRS,MONS));INDEX={v:i for i,v in enumerate(BASIS)}
WEIGHTS=[0,1,4,2,3]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text())
def mat(values):return s.Matrix([[s.sympify(v)for v in row]for row in values])
def eq(A,B):return (A-B).applyfunc(s.expand)==s.zeros(A.rows,A.cols)
def vec(P):
    polys={ij:s.Poly(P[ij[0],ij[1]],*X)for ij in PAIRS}
    return s.Matrix([polys[ij].coeff_monomial(s.prod(X[i]for i in mm))for ij,mm in BASIS])
def canon(P,variables=X):
    d=s.Matrix([sum(s.diff(P[i,j],variables[i])for i in range(5))for j in range(5)])
    E=s.Matrix(variables)
    return (P-(E*d.T-d*E.T)/5).applyfunc(s.expand)
def raw(f,g,H,cb=False):
    P=s.zeros(5);P[1,2]=f;P[3,4]=g
    P[2,3]=-s.diff(H,w);P[2,4]=s.diff(H,z)
    if cb:P[2,3]-=s.diff(f,x)*z/2;P[2,4]-=s.diff(f,x)*w/2
    return P-P.T
def lie(P,u,v):
    out=s.zeros(150,1)
    for source,((i,j),mm)in enumerate(BASIS):
        value=P[source]
        if not value:continue
        def put(r,t,nn,c):
            if r==t:return
            if r>t:r,t,c=t,r,-c
            out[INDEX[(r,t),tuple(sorted(nn))]]+=c*value
        if i==v:put(u,j,mm,1)
        if j==v:put(i,u,mm,1)
        if u in mm:
            nn=list(mm);nn.remove(u);nn.append(v);put(i,j,nn,-mm.count(u))
    return out.applyfunc(s.expand)
def check_inputs(data,paths):
    for name,expected in data.get('input_sha256',{}).items():
        p=ROOT/name;assert sha(p)==expected,name;paths.append(p)
def rank(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ).rank()

def build(args):
    start=time.monotonic();paths=[Path(__file__),ROOT/'research/notes/Q_CA_CB_N_ALL_A_20261103_v63_REVIEW.md']
    stems=['Q_CA_CB_N_SOURCE_SMITH','Q_CA_CB_N_ALL_NONZERO_COMPLETED_MODELS','Q_CA_CB_ACTUAL_INTERSECTION_TRANSLATION']
    reports={}
    for stem,replay in zip(stems,[args.source_replay,args.completion_replay,args.translation_replay]):
        p=ROOT/'reports'/f'{stem}_20261103_v63/report.json';d=load(p);assert d['all_passed']
        actual=load(replay/'report.json')
        assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in actual.items()if k!='seconds'}
        check_inputs(d,paths);paths.append(p);reports[stem]=d
    cp=ROOT/'reports/Q_CA_CB_N_SOURCE_SMITH_20261103_v63/certificate.json.gz'
    certificate=json.loads(gzip.decompress(cp.read_bytes()));paths.append(cp)
    assert sha(cp)==reports[stems[0]]['certificate_sha256']
    assert cp.read_bytes()==(args.source_replay/'certificate.json.gz').read_bytes()
    H=s.Rational(4,3)*h*x*w-h*z*z/3-x*x*z
    Q=raw(a*h*h,-2*h*h/9,H);assert canon(Q)==Q==raw(a*h*h,-2*h*h/9,H,cb=True)
    q=vec(Q);assert q==s.Matrix(list(map(s.sympify,certificate['raw_Qa150'])))
    cubics=[s.prod((h,x,z,w)[i]for i in mm)for mm in itertools.combinations_with_replacement(range(4),3)if max(mm)>=2]
    assert list(map(str,cubics))==certificate['cubic_monomials']
    A=[raw(f,0,0)for f in [h*h,h*x,x*x]]+[raw(0,h*h,0)]+[raw(0,0,v)for v in cubics]
    B=[raw(h*h,0,0,True),raw(0,h*h,0,True),raw(a*h*x,-2*h*x/9,0,True),raw(a*x*x,-2*x*x/9,0,True)]+[raw(0,0,v,True)for v in cubics]
    N=[raw(f,0,0)for f in [h*h,h*x]]+[raw(0,h*h,0)]+[raw(0,0,v)for v in cubics]
    O=[lie(q,i,j)for i,j in itertools.product(range(5),repeat=2)]
    rw=[WEIGHTS[i]+WEIGHTS[j]-WEIGHTS[m]-WEIGHTS[n]for(i,j),(m,n)in BASIS]
    assert rw==certificate['row_weights'];smith={};R=s.QQ.poly_ring(a)
    for name,source in [('CA',A),('CB',B),('N',N)]:
        M=s.Matrix.hstack(*[vec(canon(P))for P in source],*O)
        saved=certificate['families'][name];assert M==mat(saved['full150_source_orbit_matrix'])
        assert all(not M[i,j]or rw[i]==saved['column_weights'][j]for i in range(150)for j in range(M.cols))
        all_diagonal=[];blocks=[]
        for block in saved['blocks']:
            rows,cols=block['rows'],block['columns'];small=M[rows,cols]
            assert small==mat(block['matrix'])
            assert rows==[i for i,v in enumerate(rw)if v==block['weight']]
            assert cols==[i for i,v in enumerate(saved['column_weights'])if v==block['weight']]
            U,V=mat(block['left']),mat(block['right'])
            D=(U*small*V).applyfunc(s.expand)
            diagonal=list(map(s.sympify,block['diagonal']))
            assert all(D[i,j]==(diagonal[i]if i==j else 0)for i in range(D.rows)for j in range(D.cols))
            dets=[]
            for T,key in [(U,'left_determinant'),(V,'right_determinant')]:
                det=R.to_sympy(DomainMatrix.from_Matrix(T).convert_to(R).det())
                assert det==s.sympify(block[key])and det!=0 and not det.has(a);dets.append(str(det))
            assert sum(v!=0 for v in diagonal)==block['generic_rank']
            all_diagonal+=diagonal;blocks.append(dict(weight=block['weight'],left_right_unit_determinants=dets))
        nonunits=[s.factor(v)for v in all_diagonal if v.has(a)]
        assert nonunits==([a/9]if name=='CB'else[])
        ranks={str(v):sum(s.expand(d.subs(a,v))!=0 for d in all_diagonal)for v in [0,s.Rational(2,9),1]}
        assert ranks==reports[stems[0]]['families'][name]['special_image_ranks']
        smith[name]=dict(full150_matrix_rebuilt=True,blocks=blocks,nonunits=list(map(str,nonunits)),ranks=ranks)
    print('Independent physical150 and every QQ[a] Smith product/determinant PASS',flush=True)

    # Independent physical-order permutation and actual nonzero-a frame.
    permutation=[0,1,4,2,3];old=s.zeros(5)
    substitution=dict(zip(X,[Y[i]for i in permutation]))
    for i in range(5):
        for j in range(5):old[permutation[i],permutation[j]]=Q[i,j].subs(substitution,simultaneous=True)
    assert [(i,v)for i,v in enumerate(coefficients(old))if v]==[(90,a),(105,-s.Rational(2,9)),(121,s.Rational(4,3)),(137,s.Rational(2,3)),(140,s.Integer(1))]
    frame=s.zeros(5);frame[0,0]=1;frame[1,1]=1/a;frame[2,4]=1/a**2;frame[3,2]=1/a;frame[4,3]=-9/(2*a)
    physical=(a*a*pushforward(old,frame)).applyfunc(s.cancel)
    model=S61();uQ=s.zeros(61,1);uQ[24]=-s.Rational(4,3);uQ[35]=3;uQ[37]=9*a/2
    P=model.reconstruct(uQ)
    assert canon(physical,Y)==canon(P,Y)
    assert s.Matrix(model.normalize(physical)['parameters'])==uQ
    model.verify_equations(list(uQ))
    assert reports[stems[1]]['S61_line']==list(map(str,uQ))
    # Frame entries have no parameter denominators other than a.
    for v in frame:
        den=s.Poly(s.denom(s.cancel(v)),a);assert len(den.terms())==1

    bundle_path=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json'
    bundle=load(bundle_path);paths.append(bundle_path)
    base=list(map(s.Symbol,bundle['base_parameters']));fiber=list(map(s.Symbol,bundle['fiber_parameters']))
    point=dict(zip(base,[0,0,s.Rational(4,3),s.Rational(3,2),0,0]))|dict.fromkeys(fiber,0);point[fiber[1]]=9*a/2
    for fa in bundle['families'].values():assert s.Matrix(list(map(s.sympify,fa['normalized_61_coordinates']))).subs(point)==uQ
    bi,qi=bundle['beta_indices_zero_based'],bundle['q_indices_zero_based']
    beta=s.Matrix(bundle['beta16_to_intrinsic24_adapter'])*uQ[bi,0]
    pairs4=list(itertools.combinations(range(4),2));K=s.Matrix([[beta[4*i+j]for j in range(4)]for i,ij in enumerate(pairs4)if ij!=(2,3)])
    assert K.nullspace()==[s.Matrix([0,1,0,0]),s.Matrix([0,0,0,1])]
    rr,zz,ww=s.symbols('r z w');g=s.Matrix([[1,rr,-ww,zz],[0,1,0,0],[0,zz,1,0],[0,ww,0,1]])
    full=s.eye(5);full[1:5,1:5]=g.inv();ug=s.Matrix(model.normalize(pushforward(P,full))['parameters'])
    constraints={};flag={}
    for name,fa in bundle['families'].items():
        pars=s.Matrix(fa['beta_left_inverse'])*ug[bi,0]
        shift=s.Matrix(list(map(s.sympify,fa['q_shift']))).subs(dict(zip(base,pars)),simultaneous=True)
        bc=s.Matrix(fa['beta_closed_constraints'])*ug[bi,0]
        qc=s.Matrix(fa['fiber_closed_constraints'])*(ug[qi,0]-shift)
        constraints[name]=[s.expand(v)for v in [*bc,*qc]if s.expand(v)!=0]
    constraints['N']=constraints['CA']+constraints['CB']
    for name,rows in constraints.items():
        assert list(map(str,rows))==reports[stems[1]]['all_a_flag_fibers'][name]['all_a_constraints']
        gb=list(s.groebner(rows,rr,zz,ww,a,order='lex',domain=s.QQ));assert gb==[rr,zz,ww**2]
        flag[name]=dict(ideal=list(map(str,gb)),over_polynomial_QQ_a=True,entire_projective_support_points=1,scheme_fiber_length=2)
    off=s.eye(5);off[1:5,1:5]=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
    uoff=s.Matrix(model.normalize(pushforward(P,off))['parameters'])
    witness=s.Matrix(bundle['families']['CA']['beta_closed_constraints'])*uoff[bi,0]
    assert any(witness)and not any(v.has(a)for v in witness)
    assert s.Matrix(bundle['families']['CA']['beta_embedding'])==s.Matrix(bundle['families']['CB']['beta_embedding'])
    assert list(map(str,witness))==reports[stems[1]]['off_chart_a_independent_witness']
    print('Actual physical frame/all-a flag ideal and complete off-chart exclusion PASS',flush=True)

    # Whole-image input: recompute every odd coefficient from all61 source maps.
    pinch_path=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    pinched=load(pinch_path);paths.append(pinch_path);p,sv,rv=s.symbols('p s r');odd_results={};all_u={}
    for name,fa in pinched['families'].items():
        pars=list(map(s.Symbol,fa['source_fiber_parameters']));u=s.Matrix(list(map(s.sympify,fa['full_61_coordinates'])));all_u[name]=u
        qv=u[qi,0];M=qv.jacobian(pars)
        assert not any(v.free_symbols&set(pars)for v in M)
        assert eq(qv,M*s.Matrix(pars)+qv.subs(dict.fromkeys(pars,0)))
        q0=qv.subs({p:0,sv:0});q1=qv.diff(p).subs({p:0,sv:0});M0=q0.jacobian(pars);selected=fa['actual_q_selected_rows'];L=M0[selected,:]
        assert L==mat(fa['fiber_linear_block_at_origin'])
        assert L.det()==s.sympify(fa['fiber_linear_block_determinant'])!=0
        assert s.expand(M[selected,:].det()-L.det())==0
        odd=(q1-M0*L.inv()*q1[selected,0]).applyfunc(s.expand)
        assert odd==s.Matrix(list(map(s.sympify,fa['first_odd_coefficient_at_fixed_target_K'])))
        OM=odd.jacobian(pars);assert OM==s.Matrix(fa['first_odd_coefficient_matrix'])
        desired=[v for v in pars if str(v)in ['r','k1','k3','k4','k6','k7','k8']]
        selector=s.Matrix([[int(v==c)for c in pars]for v in desired])
        assert rank(OM)==len(desired)==rank(OM.col_join(selector))
        assert odd.subs(dict.fromkeys(pars,0))==s.zeros(45,1)
        assert OM[:,[i for i,c in enumerate(pars)if c not in desired]]==s.zeros(45,len(pars)-len(desired))
        axis=dict.fromkeys(pars,0);axis[s.Symbol('k0')]=a
        assert odd.subs(axis)==s.zeros(45,1)
        assert reports[stems[1]]['whole_image_odd_ideals_on_line'][name]['pinching_codimension']==len(desired)+1
        odd_results[name]=dict(all45_independently_reconstructed=True,missing_ideal=['s',*map(str,desired)],free_fiber_coordinates=['k0','k2','k5'],all_a_axis_annihilates_odd_module=True,selected_block_determinant=str(L.det()))
    # Fixed torus has no a-dependence and preserves the pure k0 direction.
    lam,mu,BB,CC,kk=s.symbols('lambda mu B C k0',nonzero=True)
    torusH=h*(BB*x*w+CC*z*z)+kk*x*x*z
    moved=s.expand(torusH.subs({x:x/lam,z:z/mu,w:mu*w},simultaneous=True)/lam)
    expected=h*(BB*mu/lam**2*x*w+CC/(lam*mu**2)*z*z)+kk/(lam**3*mu)*x*x*z
    assert s.expand(moved-expected)==0
    C0=s.Rational(3,2);B0=s.Rational(4,3);lambda0=C0/mu**2
    assert s.cancel((B0*mu/lambda0**2-1)*C0**2)==B0*mu**5-C0**2
    assert s.cancel(C0/(lambda0*mu**2))==1
    print('All45 odd coefficients/free-k0 and fixed torus transport PASS',flush=True)

    # ACTUAL intersection: a common finite map of the same target, not an
    # inference from isomorphisms of individual completed component rings.
    builder_path=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json'
    old_path=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/report.json'
    builder,old=load(builder_path),load(old_path);paths.extend([builder_path,old_path])
    v4=s.Symbol('v4');translation=s.Symbol('translation');target=s.zeros(61,1);target[37]=1;transport={}
    for name in ['CA','CB']:
        fa=pinched['families'][name];pars=list(map(s.Symbol,fa['source_fiber_parameters']));u=all_u[name]
        velocity=s.Matrix([0,1,0,-p,0,0,p*p/3,0,0,0]);selected=fa['actual_q_selected_rows']
        assert qi[selected[4]]==37;L=u[qi,0].jacobian(pars)[selected,:];assert s.expand(L.det())==-72
        assert L*velocity==s.eye(10)[:,4]
        transformed=u.subs({c:c+translation*v for c,v in zip(pars,velocity)},simultaneous=True)
        assert eq(transformed,u+translation*target)
        # Verify the inverse finite substitution, not merely differential flow.
        reverse=transformed.subs({c:c-translation*v for c,v in zip(pars,velocity)},simultaneous=True)
        assert eq(reverse,u)
        equations=list(map(s.sympify,builder['component_presentations'][name]))
        assert len(equations)==109 and not any(f.has(v4)for f in equations)
        transport[name]=dict(full61_forward_and_inverse_translation=True,source_velocity=list(map(str,velocity)),actual109_generators_v4_independent=True)
    F=s.sympify(old['nilradical_generator']);G=s.sympify(old['target_p_multiple'])
    colon=list(map(s.sympify,old['annihilator_prime_generators']))
    assert not F.has(v4)and not G.has(v4)
    assert len(colon)==45 and not any(f.has(v4)for f in colon)
    assert 'v4'in old['annihilator_free_variables']
    print('Same-target finite61 transport + actual218/colon45 literal independence PASS',flush=True)
    dependencies=[
        'AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json',
        'AMBIENT_PINCHED_MODELS_20261103_v63_REVIEW.json',
        'AMBIENT_N_GLOBAL_R1_20261103_v63_REVIEW.json',
        'AMBIENT_R19_GEOMETRY_20261103_v63_REVIEW.json',
        'AMBIENT_Q_NORMALIZATION_20261103_v63_REVIEW/report.json',
        'AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63_REVIEW/report.json']
    for name in dependencies:
        path=ROOT/'reports'/name;assert load(path)['all_passed'];paths.append(path)
    paths+=[ROOT/'code/s61_model.py',ROOT/'legacy/v02/data/symplectic_slice_certificate.json',ROOT/'legacy/v02/data/symplectic_slice_equations.json',ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json']
    for stem in stems:paths.append(ROOT/'research/notes'/f'{stem}_20261103_v63.md')
    result=dict(id='Q_CA_CB_N_ALL_A_20261103_v63_REVIEW',all_passed=True,field='QQ[a] and exact constant torus over C',
        all_three_fresh_producer_replays_match=True,independent_physical150_and_permutation=True,source_Smith=smith,
        physical_to_S61_frame=[[str(v)for v in row]for row in frame.tolist()],physical_scalar='a^2',physical_frame_scope='a!=0 only; inverse denominators only a',
        full134_and_all150_physical_S61_identity=True,all_a_flags=flag,off_chart_exclusion_and_beta_kernel_rebuilt=True,
        whole_image_odd_modules=odd_results,whole_image_R_plus_Ip_plus_p2B_argument_accepted=True,
        all_nonzero_reduced_completed_models={'CA':'Ac8 with 19 free smooth coordinates','CB':'Ac8 with 19 free smooth coordinates','N':'Ac7 with 19 free smooth coordinates'},
        full_projective_invariants={'CA':dict(dim=28,edim=37,depth=20,mult=2),'CB':dict(dim=28,edim=37,depth=20,mult=2),'N':dict(dim=27,edim=35,depth=20,mult=2)},
        same_target_actual_intersection_transport=transport,all218_and_F_G_and45_colon_generators_v4_independent=True,
        actual_D_S61=dict(dim=18,depth=11,nilradical_square_zero=True,nilmodule_smooth_dimension=11,associated_primes=2),
        actual_D_full_projective=dict(dim=27,depth=20,nilmodule_smooth_dimension=20,associated_primes=2),
        actual_D_ring_splitting_claimed=False,physical_a0_completed_type_claimed=False,
        physical_a0_nonnormal_and_nonCM_only_by_closedness_plus_global_R1=True,
        whole_Poisson_branches_or_global_D_boundary_claimed=False,
        proof_scope='All nonzero physical Qa, including a=2/9: individual reduced CA/CB/N completions and actual CA+CB scheme-intersection transport. Source Smith alone includes a0. Physical a0 full type excluded.',
        input_sha256={str(path.relative_to(ROOT)):sha(path)for path in dict.fromkeys(paths)},seconds=time.monotonic()-start)
    with args.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print('ALL THREE PACKETS INDEPENDENT REVIEW PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-replay',type=Path,required=True)
    parser.add_argument('--completion-replay',type=Path,required=True)
    parser.add_argument('--translation-replay',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args())
