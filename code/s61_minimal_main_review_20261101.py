#!/usr/bin/env python3
"""Independent raw135 kernel, cokernel and normal-obstruction audit of finite main curves."""
import argparse, hashlib, itertools, json, time
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
from s61_model import S61, X, J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, sp_basis, infinitesimal_action

ROOT = Path(__file__).resolve().parents[1]
M = lambda x: s.Matrix(x).applyfunc(s.sympify)

def quotient_length(gb,n):
    assert gb.is_zero_dimensional
    leading=[p.LM(order=gb.order).exponents for p in gb.polys]
    counts=[]
    for degree in range(30):
        count=0
        for mon in itertools.combinations_with_replacement(range(n),degree):
            exp=tuple(mon.count(i) for i in range(n))
            count+=not any(all(a<=b for a,b in zip(lead,exp)) for lead in leading)
        counts.append(count)
        if not count:return sum(counts),counts
    raise AssertionError('length enumeration did not terminate')

def local_initial(encoded,beta,q,cols,T,rows,pivots):
    """Rebuild all quadratic and cubic residuals without the primary local helper."""
    fixed=next(i for i in range(45) if i not in cols)
    A=s.zeros(135,45); terms=[]; products={}
    for row,eq in enumerate(encoded):
        ts={}
        for mon,bmon,co in eq:
            if bmon not in products:products[bmon]=s.prod(beta[i] for i in bmon)
            co*=products[bmon]
            if len(mon)==2:
                i,j=mon
                if i!=fixed:A[row,cols.index(i)]+=co*q[j]
                if j!=fixed:A[row,cols.index(j)]+=co*q[i]
            elif len(mon)==1:A[row,44]+=co*q[mon[0]]
            if fixed not in mon:
                ij=tuple([cols.index(i) for i in mon]+[44]*(2-len(mon)))
                ts[ij]=ts.get(ij,0)+co
        terms.append(ts)
    assert A*T==s.zeros(135,T.cols)
    U=s.zeros(45,len(pivots)); inverse=A[rows,pivots].inv()
    for i,j in enumerate(pivots):U[j,:]=inverse[i,:]
    n=T.cols; mon2=list(itertools.combinations_with_replacement(range(n),2)); mon3=list(itertools.combinations_with_replacement(range(n),3)); index={m:i for i,m in enumerate(mon3)}
    H=s.zeros(135,len(mon2))
    for row,ts in enumerate(terms):
        for (i,j),co in ts.items():
            for k,(a,b) in enumerate(mon2):H[row,k]+=co*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
    W2=-U*H[rows,:]; Q2=H+A*W2; H3=s.zeros(135,len(mon3))
    for row,ts in enumerate(terms):
        for (i,j),co in ts.items():
            for a in range(n):
                for k,mon in enumerate(mon2):H3[row,index[tuple(sorted((a,)+mon))]]+=co*(T[i,a]*W2[j,k]+T[j,a]*W2[i,k])
    Q3=H3-A*U*H3[rows,:]; ar=list(Q2.T.rref()[1]); ac=list(Q2.rref()[1])
    pure=Q3-Q2[:,ac]*Q2[ar,ac].inv()*Q3[ar,:]; br=list(pure.T.rref()[1]); v=s.symbols('v0:'+str(n))
    forms=[s.expand(sum(co*s.prod(v[i] for i in mon) for co,mon in zip(row,mons))) for matrix,mons in [(Q2[ar,:],mon2),(pure[br,:],mon3)] for row in matrix.tolist()]
    return A,forms,len(ar),len(br)

def run():
    started = time.time()
    path = ROOT/'reports/S61_MINIMAL_MAIN_REMOVAL_20261101.json'
    report = json.loads(path.read_text())
    for name, digest in report['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    model = S61()
    one = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 1]
    two = [i for i,w in enumerate(model.reduced['parameter_weights']) if w == 2]
    qi = {i:j for j,i in enumerate(two)}
    ui = {i:j for j,i in enumerate(one)}
    equations = model.reduced['equations'] + [json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    encoded = [[(tuple(qi[i] for i in mon if i in qi), tuple(ui[i] for i in mon if i in ui), s.Rational(n,d)) for mon,n,d in eq] for eq in equations]
    c = s.Symbol('c'); b = s.symbols('b0:9'); variables = (c,*b)
    R,*_ = ring(','.join(map(str,variables)), s.QQ)
    poly = R.from_expr
    x,y,z,w = X[1:]; spatial = s.Matrix([x,y,z,w])

    def tensor(kind,a,bb,c):
        Q=s.zeros(4); Q[0,1]=x*x if kind=='M0' else x*x/2
        if kind=='M1': Q[2,3]=x*x/2
        Q[1,2]=-bb*x*x-c*x*z-(x*z if kind=='M1' else 0)
        Q[1,3]=a*x*x+c*x*w
        return Q-Q.T

    def coordinates(Q):
        P=X[0]**2*J; P[1:,1:]+=Q
        return s.Matrix(model.normalize(P)['parameters'])[two,:]

    def clear(vector):
        denominator=s.factor(s.lcm([s.denom(v) for v in vector]))
        return [poly(s.cancel(v*denominator)) for v in vector],poly(denominator)

    zz,bu,beta_form=primitive_beta(model)
    vec=s.Matrix([beta_form[i,j].diff(zz[k]) for (i,j),k in COORDS])
    basis=vec.jacobian(bu); rr=list(basis.T.rref()[1]); inverse=basis[rr,:].inv()
    action=s.Matrix.hstack(*[inverse*infinitesimal_action(vec,A)[rr,:] for A in sp_basis()])
    source=json.loads((ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json').read_text())
    source_beta=M(source['source_beta'])
    checks={}
    for kind,u in report['uniform'].items():
        st=time.time(); beta=M(u['beta']); q=M(u['q_point']); cols=u['chart_columns']
        fixed=next(i for i in range(45) if i not in cols)
        normalized_Q=tensor(kind,0,0,c if kind=='M0' else (c-1)/2)*(1 if kind=='M0' else 2)
        assert coordinates(normalized_Q)==q
        bp=[poly(v) for v in beta]; qp=[poly(v) for v in q]
        A=[[R.zero for j in range(45)] for i in range(135)]
        raw_terms=[]; products={}
        for row,eq in enumerate(encoded):
            value=R.zero; terms={}
            for mon,bmon,co in eq:
                if bmon not in products:
                    p=R.one
                    for i in bmon:p*=bp[i]
                    products[bmon]=p
                coefficient=R.ground_new(co)*products[bmon]
                if len(mon)==2:
                    i,j=mon; value+=coefficient*qp[i]*qp[j]
                    if i!=fixed:A[row][cols.index(i)]+=coefficient*qp[j]
                    if j!=fixed:A[row][cols.index(j)]+=coefficient*qp[i]
                elif len(mon)==1:A[row][44]+=coefficient*qp[mon[0]]
                if fixed not in mon:
                    ij=tuple([cols.index(i) for i in mon]+[44]*(2-len(mon)))
                    terms[ij]=terms.get(ij,R.zero)+coefficient
            assert not value, (kind,row,'highest')
            raw_terms.append(terms)
        T=M(u['full_kernel_frame']); L=M(u['normal_cokernel'])
        assert T[u['frame_identity_rows'],:]==s.eye(T.cols)
        cleared_T=[clear(T[:,j]) for j in range(T.cols)]
        cleared_L=[clear(L[i,:]) for i in range(L.rows)]
        for vector,den in cleared_T:
            for row in A:assert not sum((v*a for v,a in zip(vector,row)),R.zero)
        for vector,den in cleared_L:
            for j in range(45):assert not sum((v*A[i][j] for i,v in enumerate(vector) if v),R.zero)
        normal=u['normal_columns']; monomials=list(itertools.combinations_with_replacement(normal,2))
        expected=M(u['normal_minor'])
        for i,(left,ld) in enumerate(cleared_L):
            for k,(aa,bb) in enumerate(monomials):
                va,da=cleared_T[aa]; vb,db=cleared_T[bb]; actual=R.zero
                for row,lv in enumerate(left):
                    if not lv:continue
                    for (j,h),co in raw_terms[row].items():
                        actual+=lv*co*(va[j]*vb[h]+(vb[j]*va[h] if aa!=bb else R.zero))
                numerator,denominator=s.fraction(s.cancel(expected[i,k]))
                assert actual*poly(denominator)==poly(numerator)*ld*da*db,(kind,i,k)
        assert s.cancel(expected.det()-s.sympify(u['normal_determinant']))==0
        expression_A=s.Matrix([[v.as_expr() for v in row] for row in A])
        sample={c:5,**dict(zip(b,range(1,10)))}
        assert expression_A.subs(sample).rank()==45-T.cols
        beta_matrix=s.zeros(134,16)
        for row,eq in enumerate(encoded[:134]):
            for mon,bmon,co in eq:
                if len(mon)==2 and len(bmon)==1:beta_matrix[row,bmon[0]]+=co*q[mon[0]]*q[mon[1]]
        rank_certificate=report['rank_cover'][kind]
        # Scaling the M1 representative by two scales this matrix by four.
        rank_beta=beta_matrix/(1 if kind=='M0' else 4)
        determinant=s.factor(rank_beta[rank_certificate['rows'],rank_certificate['columns']].det())
        assert s.expand(determinant-s.sympify(rank_certificate['minor']))==0
        rank_roots=s.solve(determinant,c)
        assert set(rank_roots)=={s.sympify(x['c']) for x in rank_certificate['exceptions']}
        for exception in rank_certificate['exceptions']:
            assert rank_beta.subs(c,s.sympify(exception['c'])).rank()==exception['rank']
        assert (beta_matrix*beta.jacobian(b)).applyfunc(s.expand)==s.zeros(134,9)
        # Reconstruct the fixed-Q infinitesimal stabilizer, and its action on all beta coordinates.
        O=[]
        for B in sp_basis():
            dq=coordinates(sum(((B*spatial)[i]*normalized_Q.diff(spatial[i]) for i in range(4)),s.zeros(4))-B*normalized_Q-normalized_Q*B.T)
            dq-=q*dq[fixed]/q[fixed]; O.append(dq)
        O=s.Matrix.hstack(*O); geo=report['geometry'][kind]
        S=M(geo['stabilizer']); actions=M(geo['stabilizer_actions'])
        assert (O*S).applyfunc(s.cancel)==s.zeros(45,3) and S.subs(c,5).rank()==3
        assert O.subs(c,5).rank()==7
        K=beta.jacobian(b); ar=action.subs(dict(zip(bu,beta)))
        assert (K*actions-ar*S).applyfunc(s.cancel)==s.zeros(16,3)
        assert actions[6:,:]==s.Matrix([[0,-b[6],0],[0,0,b[7]],[0,0,-b[8]]])
        sp={s.Symbol(k):s.sympify(v) for k,v in geo['universal_main_source_point'].items()}
        assert (source_beta.subs(sp)-beta).applyfunc(s.cancel)==s.zeros(16,1)
        # The source convention for M1 is half the projective chart representative.
        sq=M(source['cases'][kind]['q_coordinates'])
        assert (sq.subs(sp)*(1 if kind=='M0' else 2)-q).applyfunc(s.cancel)==s.zeros(45,1)
        proper=[]
        for entry in geo['proper_images']:
            certificate=entry['certificate']; bbeta=M(certificate['beta'])
            params=list(map(s.sympify,certificate['parameters'])); differential=bbeta.jacobian(params).row_join(action.subs(dict(zip(bu,bbeta))))
            left=M(certificate['left_kernel'])
            assert (differential.T.to_DM()*left.to_DM()).is_zero_matrix
            assert left.to_DM().rank()==left.cols
            assert 16-left.cols==certificate['image_dimension_upper']
            proper.append(dict(condition=entry['condition'],image_dimension_upper=16-left.cols))
        critical=geo['noninvariant_normal_factor']
        if critical:
            P=s.sympify(critical['polynomial']); coefficient=s.diff(P,b[4])
            derivative=s.expand(sum(s.diff(P,x)*actions[i,0] for i,x in enumerate(b)))
            assert s.degree(P,b[4])==1 and s.gcd(coefficient,P.subs(b[4],0))==1
            assert s.expand(derivative-s.sympify(critical['stabilizer_derivative']))==0
            assert s.gcd(P,derivative)==1
        # Independently enumerate every irreducible numerator and denominator factor.
        seen=set()
        for key in ['first_minor','jacobian_residual_determinant','all_frame_denominator_factors','normal_determinant']:
            for part in s.fraction(s.sympify(u[key])):
                for factor,power in s.factor_list(part)[1]:
                    if factor.free_symbols:seen.add(s.Poly(factor,*variables).monic().as_expr())
        saved={s.Poly(s.sympify(e['factor']),*variables).monic().as_expr() for e in report['complete_factor_cover'][kind]}
        assert seen==saved and len(seen)==18
        # Complete finite symplectic normalization and Weyl exchange, including all a,b.
        tr=report['finite_transformations'][kind]; a,bb=s.symbols('a b')
        shear=M(tr['normalization']); exchange=M(tr['exchange']); cc=c if kind=='M0' else (c-1)/2
        assert (shear*J[1:,1:]*shear.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4)
        assert (shear.inv()*tensor(kind,a,bb,cc).subs(dict(zip(spatial,shear*spatial)),simultaneous=True)*shear.inv().T-tensor(kind,0,0,cc)).applyfunc(s.cancel)==s.zeros(4)
        assert exchange*J[1:,1:]*exchange.T==J[1:,1:]
        assert (exchange.inv()*tensor(kind,0,0,cc).subs(dict(zip(spatial,exchange*spatial)),simultaneous=True)*exchange.inv().T-tensor(kind,0,0,-c if kind=='M0' else (-c-1)/2)).applyfunc(s.expand)==s.zeros(4)
        transformed_form=exchange.T*beta_form.subs(dict(zip(zz,exchange*s.Matrix(zz))),simultaneous=True)*exchange
        transformed_vec=s.Matrix([transformed_form[i,j].diff(zz[h]) for (i,j),h in COORDS])
        transformed_beta=(inverse*transformed_vec[rr,:]).subs(dict(zip(bu,beta)))
        new_parameters=list(map(s.sympify,tr['exchanged_beta_parameters']))
        assert (transformed_beta-beta.subs(c,-c).subs(dict(zip(b,new_parameters)),simultaneous=True)).applyfunc(s.expand)==s.zeros(16,1)
        assert s.expand(new_parameters[8].subs(b[7],0))==0
        # Rational exception points: raw local equations, actual same-beta tangent, full cubic initial ideal.
        exceptional=[]; pars=list(map(s.sympify,source['source_variables']))
        for case in geo['exceptional_strata']:
            local=case['local_data']; bound=case['main_branch']['initial_bound']; point=M(local['point61']); ppbeta=point[one,:]; ppq=point[two,:]; tt=M(local['tangent_basis']); chart=local['chart_columns']
            case_sample={s.Symbol(k):s.sympify(value) for k,value in case['sample'].items()}
            if case['condition']=='b7=0':
                stratum_beta=beta.subs(b[7],0); stratum_q=q; stratum_orbit=O; stratum_parameters=[c]+[x for x in b if x!=b[7]]
            else:
                case_c=case_sample[c]; stratum_beta=beta.subs(c,case_c); stratum_q=q.subs(c,case_c); stratum_orbit=O.subs(c,case_c); stratum_parameters=list(b)
            stratum_D=stratum_beta.jacobian(stratum_parameters).row_join(action.subs(dict(zip(bu,stratum_beta))))
            stratum_joint=stratum_D.col_join(stratum_q.jacobian(stratum_parameters).row_join(stratum_orbit))
            assert stratum_D.subs(case_sample).rank()==16 and stratum_joint.subs(case_sample).rank()==16
            aa,forms,qrank,crank=local_initial(encoded,ppbeta,ppq,chart,tt,bound['jacobian_pivot_rows'],bound['jacobian_pivot_columns'])
            assert aa==M(local['jacobian'])
            saved=list(map(s.sympify,bound['initial_forms'])); assert len(forms)==len(saved) and all(s.expand(a-b)==0 for a,b in zip(forms,saved))
            assert qrank==bound['quadratic_rank'] and crank==bound['pure_cubic_rank']
            v=s.symbols('v0:'+str(tt.cols)); section=s.groebner([f.subs(v[0],0) for f in forms],*v[1:],order='grevlex'); localized=s.groebner([f.subs(v[0],1) for f in forms],*v[1:],order='grevlex')
            sl,_=quotient_length(section,tt.cols-1); gl,hf=quotient_length(localized,tt.cols-1)
            assert gl==bound['multiplicity_upper_bound'] and hf==bound['standard_degrees']
            branch=case['main_branch']; point_sub={s.Symbol(k):s.sympify(v) for k,v in branch['source_point'].items()}
            assert source_beta.subs(point_sub)==ppbeta and sq.subs(point_sub)==ppq
            BD=source_beta.jacobian(pars).subs(point_sub).row_join(action.subs(dict(zip(bu,ppbeta))))
            dq=sq.jacobian(pars).subs(point_sub); fi=next(i for i in range(45) if i not in chart); dq-=ppq*dq[fi,:]/ppq[fi]
            sample_c=s.sympify(case['sample']['c']); QQ=tensor(kind,0,0,sample_c if kind=='M0' else (sample_c-1)/2); oo=[]
            for G in sp_basis():
                value=coordinates(sum(((G*spatial)[i]*QQ.diff(spatial[i]) for i in range(4)),s.zeros(4))-G*QQ-QQ*G.T);value-=ppq*value[fi]/ppq[fi];oo.append(value)
            QD=dq[chart,:].row_join(s.Matrix.hstack(*oo)[chart,:]); direction=M(branch['source_direction']); tangent=M(branch['main_tangent'])
            assert BD*direction==s.zeros(16,1) and (QD*direction).col_join(s.zeros(1,1))==tt*tangent and tangent[0]!=0
            assert BD.rank()==16 and BD.col_join(QD).rank()==17
            exceptional.append(dict(condition=case['condition'],stratum_beta_rank=16,stratum_joint_rank=16,full_rank=aa.rank(),quadratic_rank=qrank,pure_cubic_rank=crank,initial_section_length=sl,initial_generic_length=gl,same_beta_main_tangent_verified=True))
        checks[kind]=dict(raw135_highest_verified=True,raw_J_times_T_zero=True,raw_L_times_J_zero=True,kernel_identity_rows_verified=True,normal_coefficient_count=expected.rows*expected.cols,normal_determinant=u['normal_determinant'],sample_full_rank=45-T.cols,all_H_main_graph_recovered=True,all_beta_kernel_jump_values_verified=True,stabilizer_dimension=3,proper_images=proper,primitive_linear_critical_factor_and_coprime_Lie_derivative_verified=bool(critical),complete_factor_count=len(seen),full_finite_transformations_verified=True,rational_exceptions=exceptional,seconds=time.time()-st)
        print(kind,'raw135 universal checks PASS',time.time()-st,flush=True)
    # The two conjugate exceptional values are checked over QQ(sqrt(7)), never numerically.
    field=report['algebraic_exception']; a=s.Symbol('a'); alpha=s.sqrt(7)/7; domain=s.QQ.algebraic_field(alpha)
    field_forms=[]; tangent=M(field['tangent']).subs(a,alpha); n=tangent.cols; v=s.symbols('v0:'+str(n))
    for key,degree in [('quadratic',2),('pure_cubic',3)]:
        mons=list(itertools.combinations_with_replacement(range(n),degree))
        for row in M(field[key]).subs(a,alpha).tolist():field_forms.append(s.expand(sum(co*s.prod(v[i] for i in mon) for co,mon in zip(row,mons))))
    field_jac=M(field['jacobian']).subs(a,alpha)
    assert (field_jac*tangent).applyfunc(s.simplify)==s.zeros(135,n)
    assert field_jac.to_DM().convert_to(domain).rank()==field['rank']==42
    assert M(field['fixed_c_beta_differential']).subs(a,alpha).to_DM().convert_to(domain).rank()==16
    assert M(field['main_beta_differential']).subs(a,alpha).to_DM().convert_to(domain).rank()==16
    assert s.simplify(M(field['main_tangent'])[0].subs(a,alpha))!=0
    sec=s.groebner([f.subs(v[0],0) for f in field_forms],*v[1:],order='grevlex',domain=domain); loc=s.groebner([f.subs(v[0],1) for f in field_forms],*v[1:],order='grevlex',domain=domain)
    sl,_=quotient_length(sec,n-1); gl,_=quotient_length(loc,n-1); assert sl==field['section_length'][0]==3 and gl==field['localized_length'][0]==2
    paths=[Path(__file__),path,ROOT/'research/notes/S61_MINIMAL_MAIN_REMOVAL_20261101.md',ROOT/'code/s61_minimal_main_removal_20261101.py',ROOT/'code/s61_minimal_main_uniform_data_20261101.py',ROOT/'code/s61_minimal_main_exception_field_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py']
    paths += [ROOT/name for name in report['input_sha256'] if ROOT/name not in paths]
    return dict(id='S61_MINIMAL_MAIN_REVIEW_20261101',all_passed=True,universal=checks,quadratic_field_exception=dict(field='QQ[a]/(7*a^2-1)',full_rank=42,fixed_c_beta_rank=16,main_beta_rank=16,saved_kernel_verified=True,saved_initial_section_length=sl,saved_initial_generic_length=gl,both_embeddings_covered=True),all_primary_input_hashes_verified=True,proof_review=dict(complete_finite_parameter_cover_accepted=True,all_open_same_smooth_main_curve_argument_accepted=True,noninvariant_P_intrinsic_support_argument_accepted=True,generic16D_stratum_descent_accepted=True),scope=['Universal raw135 kernel/cokernel and all ten normal quadratic coefficients independently reconstructed.','The heavy primary symbolic Schur elimination and its determinant computation are audited in source, not independently rerun.','All five rational exceptional raw local Jacobians, quadratic/cubic initial forms and same-beta branch tangents are independently reconstructed.','The quadratic-field saved kernel/ranks and initial quotient lengths are rechecked exactly; its high-order contraction is source-audited, not independently regenerated.','Prior baseline/tail/Jordan theorems remain explicit dependencies.','No projective main-pole or full affine-fiber exhaustion assertion.'],seconds=time.time()-started,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
