#!/usr/bin/env python3
"""Whole-source certificates for H + two distinct irreducible quadrics.

Initial mode: all strata of the generic pair-zero semisimple weight source.
Reuses independent sparse Fraction arithmetic, not a CAS or source producer.
"""
import hashlib,itertools,json,resource,time,sys,signal
from pathlib import Path
from classification_quadratic_cofactor_review import Ring,PAIRS



def mixed_normal_support():
    """Exact necessary support covers for the five mixed-normal supersources."""
    import sympy as s
    start=time.monotonic();signal.alarm(30)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    here=Path(__file__);root=here.resolve().parents[1];path=root/'reports/CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SYMMETRIES_20261107.json'
    source={v['name']:v for v in json.loads(path.read_text())['cases']};records=[]
    x,y,z,w=s.symbols('x y z w');nn=s.symbols('n0:13');loc={str(v):v for v in [x,y,z,w,*nn]}
    for name in source:
        c=source[name];parse=lambda v:s.sympify(v,locals=loc);P=list(map(parse,c['full_brackets']));J=list(map(parse,c['all_four_Jacobi']));F=parse(c['full_Pfaffian']);branches=[]
        if name=='double_pair_one_nil':
            assert s.rem(F,y,x,y,z,w)==0
            result='Excluded: the full Pfaffian has the linear factor y, incompatible with two irreducible quadrics.'
        elif name=='double_pair_two_nil':
            a=nn[7];restriction={nn[5]:-a/2,nn[3]:0,nn[0]:0};f=s.factor(F.subs(restriction));assert s.rem(f,y*w,x,y,z,w)==0
            assert -(2*nn[5]+a)**2 in list(map(parse,c['distinct_nonzero_Jacobi_coefficients']))
            branches.append({'name':'ss_coefficient_nonzero','necessary_substitutions':{str(k):str(v) for k,v in restriction.items()},'inverted':['n7'],
                'full_Pfaffian_after_substitution':str(f),'exclusion':'linear factors y,w; rank4 HQQ impossible'})
            branches.append({'name':'ss_coefficient_zero','necessary_substitutions':{'n7':'0','n5':'0'},'inverted':['n8'],
                'reduction':'divP=n8*N with canonical nil22-image, hence all reduced HQQ points are plane by TWO_QUADRICS_NIL22_BASE_SUPPORT.md'})
            result='Every target point reduces to the already-covered whole nil22-image plane source; its full image dimension is at most24.'
        else:
            if name=='single_pair_linear_nil':
                defs=[('all',{nn[1]:0,nn[2]:0,nn[7]:0,nn[10]:4*nn[5],nn[6]:-nn[0]*nn[8]/(2*nn[5]),nn[11]:-2*nn[9]},[0,3,4,5,8,9],['n0','n5','n8'])];normalizer=6
            elif name=='single_pair_translation':
                base={nn[1]:0,nn[2]:0,nn[10]:0,nn[11]:2*nn[6]}
                defs=[('n0_nonzero',base|{nn[12]:-nn[8],nn[7]:(nn[5]*nn[8]-nn[0]*nn[9])/nn[6]},[0,3,4,5,6,8,9],['n0','n6','n9']),
                    ('n0_zero',base|{nn[0]:0,nn[7]:nn[5]*nn[8]/nn[6]},[3,4,5,6,8,9,12],['n6','n9'])];normalizer=7
            else:
                base={nn[2]:0,nn[5]:0}
                defs=[('n1_nonzero',base|{nn[8]:-4*nn[6],nn[0]:nn[1]*nn[4]/nn[6]},[1,3,4,6,7],['n1','n6','n7']),
                    ('n1_zero',base|{nn[1]:0,nn[8]:-2*nn[6]},[0,3,4,6,7],['n0','n6','n7'])];normalizer=5
            for label,sub,free,inverted in defs:
                assert all(s.cancel(v.subs(sub,simultaneous=True))==0 for v in J)
                pp=[s.factor(v.subs(sub,simultaneous=True)) for v in P];pf=s.factor(F.subs(sub,simultaneous=True))
                assert set().union(*(v.free_symbols for v in pp))<=set([x,y,z,w,*[nn[i] for i in free]])
                assert len(free)==normalizer
                branches.append({'name':label,'necessary_substitutions':{str(k):str(v) for k,v in sub.items()},'inverted':inverted,
                    'free_parameters':[str(nn[i]) for i in free],'source_dimension_upper':len(free),'normalizer_dimension_lower':normalizer,
                    'full_brackets_after_substitution':list(map(str,pp)),'full_Pfaffian_after_substitution':str(pf),
                    'all_original_Jacobi_zero':True,'full_PGL_projective_image_dimension_upper':24+len(free)-1-normalizer})
            result='Whole target subset is covered by these necessary complete-source strata; each PGL projective image has dimension<=23 and cannot be an unknown W-dominant component of dimension>=25. No component recognition is asserted.'
        records.append({'name':name,'branches':branches,'conclusion':result})
        print(name,'SUPPORT CERTIFIED',len(branches),flush=True)
    out={'task_id':'E_reduced_nonnormal','cases':records,'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
        'producer_bytes':here.read_text(),'input_sha256':{str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()},'seconds':time.monotonic()-start,
        'scope':'Exact polynomial/rational substitutions for the necessary HQQ/cofactor support strata of all five mixed-normal supersources. Completeness of necessary restrictions and normalizer dimension argument are separate written proofs; this is a G non-dominance result, not an E component-membership assertion.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SUPPORT_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')

def mixed_normal_symmetries():
    """Finite complete supersources for joint-semisimple-rank-one normals."""
    import sympy as s
    from sympy.polys.matrices import DomainMatrix
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    x,y,z,w=s.symbols('x y z w');xx=(x,y,z,w);cc=s.symbols('c0:90');ss,tt=s.symbols('ss tt')
    mons=[s.prod(v) for degree in range(3) for v in itertools.combinations_with_replacement(xx,degree)];mon1=[s.Integer(1),*xx]
    cases=[('double_pair_one_nil',[x,y,-z,-w],[y,0,0,0]),
        ('double_pair_two_nil',[x,y,-z,-w],[y,0,w,0]),
        ('single_pair_linear_nil',[x,-y,0,0],[0,0,w,0]),
        ('single_pair_translation',[x,-y,0,0],[0,0,1,0]),
        ('single_pair_long_translation',[x,-y,0,0],[0,0,w,1])]
    records=[]
    for name,S,N in cases:
        S=s.Matrix(S);N=s.Matrix(N);entries=s.Matrix([sum(cc[15*i+j]*mon for j,mon in enumerate(mons)) for i in range(6)])
        P=s.zeros(4)
        for (i,j),v in zip(PAIRS,entries):P[i,j]=v;P[j,i]=-v
        Z=ss*S+tt*N;equations=[]
        equations += [s.Poly(sum(s.diff(P[i,j],xx[i]) for i in range(4))-Z[j],*xx).coeff_monomial(mon) for j in range(4) for mon in mon1]
        for field in [S,N]:
            for i,j in PAIRS:
                lie=sum(field[k]*s.diff(P[i,j],xx[k])-P[k,j]*s.diff(field[i],xx[k])-P[i,k]*s.diff(field[j],xx[k]) for k in range(4))
                equations +=[s.Poly(lie,*xx).coeff_monomial(mon) for mon in mons]
        M,b=s.linear_eq_to_matrix(equations,[*cc,ss,tt]);assert M.shape==(200,92) and not any(b)
        rr,piv=DomainMatrix.from_Matrix(M).convert_to(s.QQ).rref();rr=rr.to_Matrix()
        free=[i for i in range(92) if i not in piv];nn=s.symbols('n0:'+str(len(free)));coef=s.zeros(92,1)
        for j,n in zip(free,nn):coef[j]=n
        for i,j in enumerate(piv):coef[j]=-sum(rr[i,k]*coef[k] for k in free)
        assert not any(s.expand(v) for v in M*coef)
        _,pr=DomainMatrix.from_Matrix(M.T).convert_to(s.QQ).rref();det=M.extract(pr,piv).det(method='domain-ge')
        source=[s.expand(sum(coef[15*i+j]*mon for j,mon in enumerate(mons))) for i in range(6)];matrix=s.zeros(4)
        for (i,j),v in zip(PAIRS,source):matrix[i,j]=v;matrix[j,i]=-v
        jac=[s.expand(sum(matrix[a,t]*s.diff(matrix[b,d],xx[t]) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
        coefficients=[]
        for poly in jac:
            for v in s.Poly(poly,*xx).coeffs():
                v=s.factor(v)
                if v and v not in coefficients and -v not in coefficients:coefficients.append(v)
        pf=s.factor(source[0]*source[5]-source[1]*source[4]+source[2]*source[3])
        records.append({'name':name,'semisimple_field':list(map(str,S)),'nilpotent_field':list(map(str,N)),
            'full_constraint_shape':list(M.shape),'full_constraint_sparse_rows':[{str(j):str(v) for j,v in enumerate(M.row(i)) if v} for i in range(M.rows)],
            'rank':len(piv),'pivot_rows':list(pr),'pivot_columns':list(piv),'pivot_determinant':str(det),
            'free_coefficient_indices':free,'free_parameters':list(map(str,nn)),'all92_coefficients':list(map(str,coef)),
            'full_brackets':list(map(str,source)),'divergence_ss_coefficient':str(coef[90]),'divergence_nil_coefficient':str(coef[91]),
            'all_four_Jacobi':list(map(str,jac)),'distinct_nonzero_Jacobi_coefficients':list(map(str,coefficients)),
            'full_Pfaffian':str(pf)})
        print(name,'rank',len(piv),'free',len(free),'Jac generators',len(coefficients),flush=True)
    out={'task_id':'E_reduced_nonnormal','cases':records,'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
        'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'Full200-by92 supersources satisfying divP in span(S,N), LieS P=LieN P=0 for all five possible commuting semisimple/nonzero-nilpotent symmetry types. No Jacobi ideal has been solved and no factor-type/component conclusion is drawn from this entrance alone. The independent normal-Hamiltonian and factor conditions are not discarded.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SYMMETRIES_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')

def nonzero_pair_cofactor():
    """Complete90-coefficient sources and SNC limits for all lambda!=0 types."""
    import sympy as s
    from sympy.polys.matrices import DomainMatrix
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    h,x,y,z,w=s.symbols('h x y z w');xx=(x,y,z,w);cc=s.symbols('c0:90')
    mon2=[s.prod(v) for degree in range(3) for v in itertools.combinations_with_replacement(xx,degree)]
    mon3=[s.prod(v) for degree in range(4) for v in itertools.combinations_with_replacement(xx,degree)]
    mon1=[s.Integer(1),*xx]
    types=[('rank44',h*x+y*w,h*y+x*z,[(0,-1),(1,0),(1,1),(-1,-1)],[3,0,0,0,0]),
        ('rank43',h*x+y*z,y*y+x*w,[(0,-1),(s.Rational(1,2),0),(-s.Rational(1,2),-1),(1,1)],[3,0,1,0,0]),
        ('rank33',x*x+y*z,y*y+x*w,[(0,-s.Rational(1,2)),(s.Rational(1,2),0),(-s.Rational(1,2),-1),(1,s.Rational(1,2))],[3,0,0,-1,-1])]
    records=[]
    for name,qh1,qh2,characters,weights in types:
        q1,q2=qh1.subs(h,1),qh2.subs(h,1)
        entries=s.Matrix([sum(cc[15*i+j]*mon for j,mon in enumerate(mon2)) for i in range(6)])
        P=s.zeros(4)
        for (i,j),v in zip(PAIRS,entries):P[i,j]=v;P[j,i]=-v
        fields=[s.Matrix([characters[j][a]*xx[j] for j in range(4)]) for a in range(2)]
        equations=[]
        for quad,field in zip([q1,q2],fields):
            residual=-P*s.Matrix([s.diff(quad,v) for v in xx])-quad*field
            equations +=[s.Poly(v,*xx).coeff_monomial(mon) for v in residual for mon in mon3]
        Z=fields[0]+fields[1]
        equations +=[s.Poly(sum(s.diff(P[i,j],xx[i]) for i in range(4))-Z[j],*xx).coeff_monomial(mon) for j in range(4) for mon in mon1]
        for axis in range(2):
            for i,j in PAIRS:
                lie=sum(fields[axis][k]*s.diff(P[i,j],xx[k]) for k in range(4))-(characters[i][axis]+characters[j][axis])*P[i,j]
                equations +=[s.Poly(lie,*xx).coeff_monomial(mon) for mon in mon2]
        M,b=s.linear_eq_to_matrix(equations,cc);assert M.shape==(480,90)
        aug=M.row_join(b);rr,piv=DomainMatrix.from_Matrix(aug).convert_to(s.QQ).rref();rr=rr.to_Matrix()
        assert 90 not in piv
        free=[i for i in range(90) if i not in piv];nn=s.symbols('n0:'+str(len(free)));coef=s.zeros(90,1)
        for j,n in zip(free,nn):coef[j]=n
        for i,j in enumerate(piv):coef[j]=rr[i,-1]-sum(rr[i,k]*coef[k] for k in free)
        assert all(s.expand(v)==0 for v in M*coef-b)
        _,pr=DomainMatrix.from_Matrix(M.T).convert_to(s.QQ).rref()
        pivot=M.extract(pr,piv).det(method='domain-ge')
        source=[s.expand(sum(coef[15*i+j]*mon for j,mon in enumerate(mon2))) for i in range(6)]
        matrix=s.zeros(4)
        for (i,j),v in zip(PAIRS,source):matrix[i,j]=v;matrix[j,i]=-v
        jac=[s.expand(sum(matrix[a,t]*s.diff(matrix[b,d],xx[t]) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
        assert not any(jac)
        pf=s.factor(source[0]*source[5]-source[1]*source[4]+source[2]*source[3])
        kappa=s.factor(s.cancel(pf/(q1*q2)));assert kappa and not any(s.diff(kappa,v) for v in xx)
        limit=[];term_weights=[]
        for (i,j),v in zip(PAIRS,source):
            kept=s.Integer(0);ws=[]
            for exp,coefficient in s.Poly(v,*xx).terms():
                weight=(2-sum(exp))*weights[0]+sum(a*b for a,b in zip(exp,weights[1:]))-weights[i+1]-weights[j+1]
                assert weight>=0;ws.append(weight)
                if weight==0:kept+=coefficient*s.prod(a**b for a,b in zip(xx,exp))
            assert not any(s.diff(s.cancel(kept/(xx[i]*xx[j])),v) for v in xx)
            limit.append(s.expand(kept));term_weights.append(ws)
        lpf=s.factor(limit[0]*limit[5]-limit[1]*limit[4]+limit[2]*limit[3])
        assert s.expand(lpf-kappa*x*y*z*w)==0
        records.append({'name':name,'homogeneous_q1':str(qh1),'homogeneous_q2':str(qh2),
            'joint_characters':[list(map(str,v)) for v in characters],
            'cofactor_fields':[list(map(str,v)) for v in fields],'degree_at_most2_monomials':list(map(str,mon2)),
            'full_constraint_shape':list(M.shape),'full_constraint_sparse_rows':[{str(j):str(v) for j,v in enumerate(M.row(i)) if v} for i in range(M.rows)],
            'full_constraint_rhs':list(map(str,b)),'rank':len(piv),'pivot_rows':list(pr),'pivot_columns':list(piv),'pivot_determinant':str(pivot),
            'free_coefficient_indices':free,'free_parameters':list(map(str,nn)),
            'all90_coefficients':list(map(str,coef)),'full_brackets':list(map(str,source)),
            'all_four_Jacobi_zero':True,'full_Pfaffian':str(pf),'Pfaffian_multiple_of_q1q2':str(kappa),
            'cocharacter_hxyzw':weights,'all_monomial_weights':term_weights,
            'whole_source_logcanonical_limit':list(map(str,limit)),'limit_Pfaffian':str(lpf),
            'all_rank4_source_points_have_actual_SNC_orbit_limit':True})
        print(name,'rank',len(piv),'free',len(free),'kappa',kappa,'ALL SNC',flush=True)
    out={'task_id':'E_reduced_nonnormal','cases':records,'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),
        'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'All90-coefficient sources with actual normal Hamiltonians at nonzero log-bracket lambda normalized1, for the three geometrically complete joint-character quadric types. Geometry of that complete three-type reduction and the MPS theorem are separate proof inputs; all free parameters and Pfaffian-zero complements are retained.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_NONZERO_PAIR_COFACTOR_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')

def plane_even():
    """Global six-magnetic source to beta=0, with no parameter inverses."""
    import sympy as s
    start=time.monotonic();signal.alarm(30)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json'
    row=next(v for v in json.loads(path.read_text())['cases'] if v['name']=='contact4_secant_at_contact')
    h,x,y,z,w=xx=s.symbols('h x y z w');mm=s.symbols('m0:6');P=s.zeros(5)
    for (i,j),p in zip(PAIRS,row['brackets']):P[i+1,j+1]=s.sympify(p);P[j+1,i+1]=-P[i+1,j+1]
    alpha=-3*mm[2]/5;beta=2*mm[1]-4*mm[2]/3
    translation=s.eye(5);translation[3,0]=alpha;translation[4,0]=beta
    frame=s.zeros(5);frame[0,0]=frame[1,1]=frame[3,2]=frame[4,3]=1
    frame[2,4]=s.Rational(9,10);frame[4,1]=-s.Rational(9,10)*mm[0]
    change=frame*translation;assert change.det()==s.Rational(9,10)
    inverse=change.inv()*s.Matrix(xx);sub=dict(zip(xx,inverse))
    transformed=(change*P.applyfunc(lambda p:s.expand(p.subs(sub,simultaneous=True)))*change.T).applyfunc(s.expand)
    J=s.zeros(5);J[1,2]=J[3,4]=1;J[2,1]=J[4,3]=-1
    value=transformed.subs(dict(zip(xx,[1,0,0,0,0])));assert value==J
    entries=[transformed[i+1,j+1] for i,j in PAIRS]
    assert all(s.Poly(p,*xx[1:]).coeff_monomial(v)==0 for p in [q.subs(h,1) for q in entries] for v in xx[1:])
    Q=transformed-h*h*J
    def jac(A):return [s.expand(sum(A[a,t]*s.diff(A[b,d],xx[t]) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(5))) for i,j,k in itertools.combinations(range(5),3)]
    assert not any(jac(transformed)) and not any(jac(Q))
    pf=s.factor(entries[0]*entries[5]-entries[1]*entries[4]+entries[2]*entries[3])
    out={'task_id':'E_reduced_nonnormal','name':row['name'],'all_passed':True,
        'all_six_magnetic_parameters':list(map(str,mm)),'fiber_translation_alpha_beta':list(map(str,[alpha,beta])),
        'full_change_hxyzw':change.tolist(),'change_determinant':'9/10',
        'full_transformed_six_brackets':list(map(str,entries)),
        'quadratic_Q45_six_brackets':[str(Q[i+1,j+1]) for i,j in PAIRS],
        'constant_bivector_is_standard_J':True,'all_linear_coefficients_identically_zero':True,
        'all_full_Jacobi_and_compatible_Q45_conditions':True,'full_Pfaffian':str(pf),
        'conclusion':'The entire six-magnetic source is globally carried into the audited beta=0/Q45 support by a parameter-polynomial invertible projective coordinate change, hence lies in the known finite closed union from BETA_ZERO_COVERAGE. No local uniqueness or generic-only extension is used.',
        'input_sha256':{str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_EVEN_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps({'all_passed':True,'whole_source_beta_zero':True,'seconds':out['seconds']}))

def plane_graphs():
    import sympy as s
    from s61_model import S61,X,pushforward
    from classification_degenerate_quadric_q4_nil31_I_graph import darboux
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json'
    source=json.loads(path.read_text())['cases'];old=s.symbols('h x y z w');model=S61();us=s.symbols('u1:62');records=[]
    for row in source:
        if row['ambient_tangent_upper']!=25:continue
        P=s.zeros(5)
        for key,f in row['point_brackets'].items():
            i,j=map(int,key.strip('()').split(','));P[i,j]=s.sympify(f).subs(dict(zip(old,X)),simultaneous=True);P[j,i]=-P[i,j]
        a,b=next((a,b) for a in range(1,8) for b in range(8) if P.subs(dict(zip(X,(1,a,b,0,0))))[1:5,1:5].det()!=0)
        translation=s.eye(5);translation[1,0]=-a;translation[2,0]=-b
        translated=pushforward(P,translation);value=translated.subs(dict(zip(X,(1,0,0,0,0))))[1:5,1:5]
        frame=s.eye(5);frame[1:5,1:5]=darboux(value)
        normalized=model.normalize(pushforward(translated,frame));model.verify_point(normalized,check_equations=True)
        point=normalized['parameters'];membership={}
        for name,gpath in [('T39','reports/Q51_UNIVERSAL_SECTION_20261101.json'),('T40','reports/Q52_UNIVERSAL_SECTION_20261101.json')]:
            graph=json.loads((root/gpath).read_text());indices=graph['weight2_u_indices_zero_based']
            residual=[s.cancel(point[k]-s.sympify(q).subs(dict(zip(us,point)),simultaneous=True)) for k,q in zip(indices,graph['quadratic_section_q'])]
            membership[name]={'in_entire_reduced_graph':all(v==0 for v in residual),'all45_residuals':list(map(str,residual))}
        records.append({'name':row['name'],'evaluation_point':[1,a,b,0,0],
            'translation':[[str(v) for v in row] for row in translation.tolist()],
            'Darboux_frame':[[str(v) for v in row] for row in frame.tolist()],
            'perspective':list(map(str,normalized['perspective'])),'all61_parameters':list(map(str,point)),
            'all150_reconstruction_and134_equations_checked':True,'graph_membership':membership,
            'ambient_tangent_upper':25,'whole_irreducible_source_covered_if_any_graph_positive':True})
        print(row['name'],{name:v['in_entire_reduced_graph'] for name,v in membership.items()},flush=True)
    paths=[path,root/'code/s61_model.py',root/'code/classification_degenerate_quadric_q4_nil31_I_graph.py',
        root/'reports/Q51_UNIVERSAL_SECTION_20261101.json',root/'reports/Q52_UNIVERSAL_SECTION_20261101.json']
    out={'task_id':'E_reduced_nonnormal','cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'All45 defining whole-graph equations at full134-verified S61 normalizations of the two plane sources with tangent upper25. A positive graph plus known component dimension25 forces local smoothness/uniqueness and entire irreducible source containment; negative tests leave the source open.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_GRAPHS_20261107.json').write_text(json.dumps(out,indent=2)+'\n')

def nil22_rank1_support():
    import sympy as s
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    bp=root/'reports/CLASSIFICATION_TWO_QUADRICS_NIL22_BASE_20261107.json'
    op=root/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    data=json.loads(bp.read_text())['cases'][1];original=json.loads(op.read_text())['full_nilpotent_sources'][5]
    xx=s.symbols('x y z w');x,y,z,w=xx;nn=s.symbols('n0:27');names={str(v):v for v in (*xx,*nn)}
    parse=lambda p:s.sympify(p,locals=names)
    discriminant=3*nn[18]**2+4*nn[19]*nn[24]
    oldeq=list(map(parse,original['distinct_nonzero_jacobi_coefficients']))
    assert any(s.expand(p+4*discriminant/3)==0 for p in oldeq)
    assert any(s.expand(p-4*discriminant/3)==0 for p in oldeq)
    eq=list(map(parse,data['remaining_full_Jacobi_coefficients']));source=list(map(parse,data['full_brackets']))
    assert -nn[17]**2 in eq
    after17=[s.factor(p.subs(nn[17],0)) for p in eq]
    assert -4*nn[14]**2/5 in after17
    radical={nn[17]:s.Integer(0),nn[14]:s.Integer(0)}
    source=[s.expand(p.subs(radical)) for p in source]
    eq=list(dict.fromkeys(s.factor(p.subs(radical)) for p in eq if p.subs(radical)!=0))
    substitutions=dict(radical);eliminations=[]
    for parameter in nn[3:9]:
        polynomial=next(p for p in eq if s.diff(p,parameter)!=0 and s.diff(p,parameter).is_Rational)
        coefficient=s.diff(polynomial,parameter)
        solution=s.expand(-polynomial.subs(parameter,0)/coefficient)
        substitutions={k:s.expand(v.subs(parameter,solution)) for k,v in substitutions.items()};substitutions[parameter]=solution
        source=[s.expand(p.subs(parameter,solution)) for p in source]
        eq=list(dict.fromkeys(s.factor(p.subs(parameter,solution)) for p in eq if p.subs(parameter,solution)!=0))
        eliminations.append({'parameter':str(parameter),'actual_Jacobi_equation':str(polynomial),'nonzero_constant_coefficient':str(coefficient),'exact_polynomial_solution':str(solution)})
    first=nn[12]*(3*nn[16]-nn[26]-2);second=nn[12]*(3*nn[13]-nn[22])
    assert eq and all(any(s.cancel(p/q).is_Rational and s.cancel(p/q)!=0 for q in [first,second]) for p in eq)
    assert all(any(s.cancel(p/q).is_Rational and s.cancel(p/q)!=0 for p in eq) for q in [first,second])
    branches=[]
    for label,fixed in [('closed_n12_zero',{nn[12]:0}),('closed_eigen_relations',{nn[26]:3*nn[16]-2,nn[22]:3*nn[13]})]:
        entries=[s.expand(p.subs(fixed,simultaneous=True)) for p in source]
        matrix=s.zeros(4)
        for (i,j),p in zip(PAIRS,entries):matrix[i,j]=p;matrix[j,i]=-p
        jac=[s.expand(sum(matrix[a,t]*s.diff(matrix[b,d],xx[t]) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
        assert not any(jac)
        pf=s.factor(entries[0]*entries[5]-entries[1]*entries[4]+entries[2]*entries[3])
        free=sorted(set().union(*(p.free_symbols for p in entries))-set(xx),key=str)
        record={'name':label,'fixed_parameters':{str(k):str(v) for k,v in fixed.items()},'free_parameters':list(map(str,free)),
            'full_brackets':list(map(str,entries)),'all_four_Jacobi_zero':True,'actual_Pfaffian':str(pf)}
        if label=='closed_n12_zero':
            quotient=s.cancel(pf/y**2);assert quotient.is_polynomial(*xx,*nn)
            record.update(actual_Pfaffian_has_double_linear_factor='y^2',Pfaffian_divided_by_y_squared=str(s.factor(quotient)))
        else:
            assert s.diff(pf,x)==-y*y*w and s.diff(pf,z)==y**3
            assert s.expand(pf.subs(y,0))==-nn[12]**2
            record.update(Pfaffian_is_linear_in_fibers_with_coefficient='-y^2*(w*x-y*z)',
                Pfaffian_at_y_zero='-n12**2',irreducible_at_n12_nonzero_by_primitive_linear_Gauss=True,
                n12_zero_complement_is_the_first_branch=True)
        branches.append(record)
    # Exact affine GL4 witness into the independently audited old unique
    #23-dimensional H_QUARTIC_ESSENTIAL5 component (legacy v06 Q_*).
    witness={nn[12]:-1,nn[13]:0,nn[15]:1,nn[16]:1,nn[20]:0,nn[21]:0,nn[23]:0,nn[25]:4}
    transformed_old=list(map(parse,['0','y**2-y*z+2*y*w+2*y*x+w**2/2+2*w*x+x**2',
        '1+y**2+y*x','-1+y**2+y*w+y*x','y**2','-3*y*z+4*y*w+w**2+3*w*x']))
    assert all(s.expand(parse(p).subs(witness)-q)==0 for p,q in zip(branches[1]['full_brackets'],transformed_old))
    recognition={'old_source':'legacy/v06/nilpotent22.tex, subsection quarticcomponent, P_*=h^2 J+Q_*',
        'new_coordinates_in_old_order':['x4','x1','-x2','-x3'],
        'canonical_N_after_change':'y partial_x+w partial_z','canonical_base_bracket':'y^2',
        'branch_parameter_values':{str(k):str(v) for k,v in witness.items()},
        'all_six_transformed_old_brackets':list(map(str,transformed_old)),
        'all_six_branch_equal_transformed_old_source':True,
        'old_unique_component_input':'The audited old point has full Jacobi rank102 and actual orbit+vertical+scale rank48, hence a unique smooth23-dimensional ambient component H_QUARTIC_ESSENTIAL5.',
        'whole_branch_implication':'The eight-free-parameter irreducible source and its PGL saturation contain this point; every irreducible image closure through that unique-component point lies in H_QUARTIC_ESSENTIAL5. No whole-source single-orbit claim is needed.'}
    out={'task_id':'E_reduced_nonnormal','global_base_discriminant_Jacobi_coefficient':str(-4*discriminant/3),
        'reduced_square_implications':['-n17^2=0 => n17=0','-(25*n12*n17+8*n14^2)/10=0 => n14=0'],
        'global_polynomial_eliminations':eliminations,'substitutions_after_radical':{str(k):str(v) for k,v in substitutions.items()},
        'remaining_exact_Jacobi_ideal_generators':list(map(str,[first,second])),
        'full_reduced_source_before_branching':list(map(str,source)),'branches':branches,'known_quartic_component_recognition':recognition,
        'conclusion':'Any nonzero binary-base-bracket nil22-image point has a double linear Pfaffian factor or an irreducible quartic. Hence reduced H+two distinct irreducible quadrics must have zero base bracket and is essential-three-variable/plane.',
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [bp,op,root/'legacy/v06/nilpotent22.tex']},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'Complete reduced support of both nonzero binary ranks in the original whole nil22-image source, with all square/closed complements retained. The irreducible-quartic source contains the exact old unique23-component witness; geometric whole-source containment uses that audited theorem.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_NIL22_RANK1_SUPPORT_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'branch_free_dimensions':[len(v['free_parameters']) for v in branches],'seconds':out['seconds']}))

def nil22_base():
    """Global linear closure of the two nonzero binary-base-bracket strata."""
    import sympy as s
    from sympy.polys.matrices import DomainMatrix
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    path=root/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json'
    row=json.loads(path.read_text())['full_nilpotent_sources'][5]
    assert row['field_linear_matrix']==[[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]]
    xx=s.symbols('x y z w');nn=s.symbols('n0:27');names={str(v):v for v in (*xx,*nn)}
    original=list(map(lambda v:s.sympify(v,locals=names),row['brackets']));records=[]
    for label,fixed in [('rank2_yw',{nn[18]:1,nn[19]:0,nn[24]:0}),('rank1_y2',{nn[18]:0,nn[19]:0,nn[24]:-3})]:
        fixed={k:s.sympify(v) for k,v in fixed.items()}
        source=[s.expand(v.subs(fixed)) for v in original];matrix=s.zeros(4)
        for (i,j),p in zip(PAIRS,source):matrix[i,j]=p;matrix[j,i]=-p
        jac=[s.expand(sum(matrix[a,t]*s.diff(matrix[b,d],xx[t]) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
        equations=list(dict.fromkeys(s.factor(c) for p in jac for c in s.Poly(p,*xx).coeffs() if c))
        substitutions=dict(fixed);variables=[v for v in nn if v not in fixed];stages=[];empty=False
        while True:
            if not variables:
                empty=bool(equations);break
            linear=[v for v in equations if s.Poly(v,*variables).total_degree()<=1]
            if not linear:break
            M,b=s.linear_eq_to_matrix(linear,variables);aug=M.row_join(b)
            rr,pivots=DomainMatrix.from_Matrix(aug).convert_to(s.QQ).rref();rr=rr.to_Matrix()
            if M.cols in pivots:
                empty=True;stages.append({'inconsistent_constant_linear_system':True,'matrix':aug.tolist()});break
            free=[j for j in range(M.cols) if j not in pivots]
            solve={variables[p]:s.expand(rr[i,-1]-sum(rr[i,j]*variables[j] for j in free)) for i,p in enumerate(pivots)}
            stages.append({'matrix':M.tolist(),'rhs':list(b),'pivot_columns':list(pivots),'exact_affine_solution':{str(k):str(v) for k,v in solve.items()}})
            substitutions={k:s.expand(v.subs(solve,simultaneous=True)) for k,v in substitutions.items()};substitutions.update(solve)
            source=[s.expand(v.subs(solve,simultaneous=True)) for v in source]
            equations=list(dict.fromkeys(s.factor(v.subs(solve,simultaneous=True)) for v in equations if s.factor(v.subs(solve,simultaneous=True))!=0))
            variables=[variables[j] for j in free]
        pf=s.factor(source[0]*source[5]-source[1]*source[4]+source[2]*source[3])
        records.append({'stratum':label,'fixed_base_bracket_parameters':{str(k):str(v) for k,v in fixed.items()},
            'global_linear_elimination_stages':stages,'empty':empty,
            'all_parameter_substitutions':{str(k):str(v) for k,v in substitutions.items()},
            'remaining_free_parameters':list(map(str,variables)),'remaining_full_Jacobi_coefficients':list(map(str,equations)),
            'full_brackets':list(map(str,source)),'actual_Pfaffian':str(pf),
            'remaining_factorization':[{'equation':str(eq),'factors':[(str(p),int(power)) for p,power in s.factor_list(eq)[1]]} for eq in equations],
            'original_source_coefficient_count':27})
        print(label,'empty',empty,'free',len(variables),'remaining',len(equations),'Pf',pf,flush=True)
    out={'task_id':'E_reduced_nonnormal','cases':records,
        'input_sha256':{str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'Complete original nil22-image source at the two canonical nonzero binary-base-bracket types, with only globally invertible constant-coefficient linear elimination. Every nonlinear factor and all source equations remain explicit; neither component recognition nor full support decomposition is inferred.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_NIL22_BASE_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')

def plane_cubic_support():
    import sympy as s
    from sympy.polys.matrices import DomainMatrix
    from classification_rank3_mixed_component_review import COLUMNS,jacobi,vector
    from classification_three_hyperplanes_quadric_independent import ambient_matrices
    start=time.monotonic();signal.alarm(90)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    sp=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json'
    np=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_NORMALS_20261107.json'
    sources={v['name']:v for v in json.loads(sp.read_text())['cases']}
    selected=[v for v in json.loads(np.read_text())['cases'] if not v['zero_dimensional_quadratic_support']]
    xx=s.symbols('h x1:5');old=s.symbols('h x y z w');dm=lambda A:DomainMatrix.from_Matrix(A).convert_to(s.QQ)
    def matrix(data,rename=False):
        A=s.zeros(5)
        for key,value in data.items():
            i,j=map(int,key.strip('()').split(','));value=s.sympify(value)
            if rename:value=value.subs(dict(zip(old,xx)),simultaneous=True)
            A[i,j]=value;A[j,i]=-value
        return A
    def sparse(A):return {ij:dict(s.Poly(s.expand(A[ij]),*xx).terms()) for ij in itertools.combinations(range(5),2) if A[ij]!=0}
    J=lambda A:jacobi(sparse(A))
    def encode(A):return {f'{i},{j}':str(s.expand(A[i,j])) for i,j in itertools.combinations(range(5),2) if A[i,j]}
    def tensor(v):
        A=s.zeros(5)
        for ((i,j),e),a in zip(COLUMNS,v):
            if a:A[i,j]+=a*s.prod(x**p for x,p in zip(xx,e))
        return A-A.T
    records=[]
    for saved in selected:
        begin=time.monotonic();source=sources[saved['case']];P=matrix(source['point_brackets'],True)
        universal=matrix(source['full_six_magnetic_source'],True);values={s.Symbol(k):s.Rational(v) for k,v in source['parameter_values'].items()}
        directions=[universal.diff(p).subs(values) for p in values]
        D,V,G,F=ambient_matrices(P,xx,directions);pivots=list(dm(D).rref()[1]);ind=D[:,pivots]
        rows=list(dm(ind.T).rref()[1]);inverse=dm(ind[rows,:]).inv()
        RR=list(map(matrix,saved['normals']));N=s.Matrix.hstack(*[vector(sparse(R)) for R in RR]);count=len(RR)
        assert len(dm(F).rref()[1])==48 and len(dm(F.row_join(N)).rref()[1])==150-len(pivots)
        pairs=list(itertools.combinations_with_replacement(range(count),2));squares=[J(R) for R in RR]
        Q=s.Matrix.hstack(*[squares[i] if i==j else J(RR[i]+RR[j])-squares[i]-squares[j] for i,j in pairs])
        correction=-inverse.matmul(dm(Q[rows,:])).to_Matrix();TT=[]
        for k in range(len(pairs)):
            v=s.zeros(150,1)
            for i,p in enumerate(pivots):v[p]=correction[i,k]
            TT.append(tensor(v))
        qres=Q+ind*correction
        triples=list(itertools.combinations_with_replacement(range(count),3));C=s.zeros(350,len(triples))
        for i,R in enumerate(RR):
            for pair,T in zip(pairs,TT):C[:,triples.index(tuple(sorted((i,)+pair)))]+=J(R+T)-J(R)-J(T)
        A=D.row_join(Q);apivots=list(dm(A).rref()[1]);aind=A[:,apivots];arows=list(dm(aind.T).rref()[1]);ainverse=dm(aind[arows,:]).inv()
        cres=C-aind*ainverse.matmul(dm(C[arows,:])).to_Matrix()
        t=s.symbols('t0:'+str(count));qmons=[t[i]*t[j] for i,j in pairs];cmons=[s.prod(t[i] for i in inds) for inds in triples]
        qr,qp=dm(qres).rref();qr=qr.to_Matrix();cr,cp=dm(cres).rref();cr=cr.to_Matrix()
        quadrics=[s.expand(sum(qr[i,j]*qmons[j] for j in range(len(pairs)))) for i in range(len(qp))]
        cubics=[s.expand(sum(cr[i,j]*cmons[j] for j in range(len(triples)))) for i in range(len(cp))]
        gb=s.groebner(quadrics+cubics,*t,order='grevlex') if quadrics or cubics else None
        functionals=[]
        for row in dm(cres.T).rref()[1]:
            lam=s.zeros(1,350);lam[row]=1;weights=dm(aind[row,:]).matmul(ainverse).to_Matrix()
            for j,k in enumerate(arows):lam[k]-=weights[j]
            assert lam*A==s.zeros(1,A.cols)
            functionals.append({'nonzero_coefficients':{str(i):str(a) for i,a in enumerate(lam) if a},'cubic_pairing':list(map(str,list(lam*C)))})
        closed=bool(gb and gb.is_zero_dimensional)
        records.append({'case':saved['case'],'full_D_rank':len(pivots),'actual_source_rank':48,'normal_count':count,
            'true_literal_complement_verified':True,'implicit_D_columns':pivots,'implicit_D_rows':rows,
            'implicit_minor_determinant':str(dm(ind[rows,:]).det()),'quadratic_corrections':list(map(encode,TT)),
            'normal_pairs':[list(p) for p in pairs],'quadratic_initial_generators':list(map(str,quadrics)),
            'cubic_monomials':list(map(str,cmons)),
            'all_cubic_coefficient_vectors':[{str(i):str(a) for i,a in enumerate(C[:,j]) if a} for j in range(C.cols)],
            'cubic_functionals':functionals,'cubic_initial_generators':list(map(str,cubics)),
            'initial_groebner_basis':[] if gb is None else [str(g.as_expr()) for g in gb.polys],
            'zero_dimensional_initial_support':closed,'projective_component_dimension':23 if closed else None,
            'fixed_slice_only':True,'arbitrary_kernel_correction_invariance_not_required':True,'seconds':time.monotonic()-begin})
        print(saved['case'],'quadrics',quadrics,'cubics',cubics,'support0',closed,flush=True)
    inputs=[sp,np,root/'code/classification_rank3_mixed_component_review.py',root/'code/classification_three_hyperplanes_quadric_independent.py']
    out={'task_id':'E_reduced_nonnormal','cases':records,'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'One necessary cubic initial-support test in four fixed genuine source-complement slices, with all source directions retained. Stop higher-jet iteration on any nonisolated result; component use requires independent arithmetic and geometric review.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CUBIC_SUPPORT_20261107.json').write_text(json.dumps(out,indent=2)+'\n')

def plane_normals():
    import sympy as s
    from classification_two_hyperplanes_cubic_nilpotent import normal_quadrics
    start=time.monotonic();signal.alarm(90)
    here=Path(__file__);root=here.resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json'
    data=json.loads(path.read_text());old=s.symbols('h x y z w');new=s.symbols('h x1:5');records=[]
    def encode(values):
        out={}
        for key,value in values.items():
            i,j=map(int,key.strip('()').split(','));out[f'{i},{j}']=str(s.sympify(value).subs(dict(zip(old,new)),simultaneous=True))
        return out
    for c in data['cases']:
        record={'name':c['name'],'witness_brackets':encode(c['point_brackets']),
            'universal_brackets':encode(c['full_six_magnetic_source']),
            'full_source_specification':{'free_parameters':list(c['parameter_values']),'witness':c['parameter_values']}}
        result=normal_quadrics(path,0,None,record_override=record);result.pop('producer_source')
        assert result['normal_count']==c['true_normal_dimension']
        records.append(result)
        print(result['case'],'normals',result['normal_count'],'quadrics',result['quadratic_initial_generators'],'support0',result['zero_dimensional_quadratic_support'],flush=True)
    inputs=[path,root/'code/classification_two_hyperplanes_cubic_nilpotent.py',root/'code/classification_rank3_mixed_component_review.py']
    out={'task_id':'E_reduced_nonnormal','cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'First quadratic obstruction in a genuine complement to all150-variable Jacobi and all actual source directions. Zero projective support suffices for component recognition after independent review; surviving lines remain explicit obligations.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_NORMALS_20261107.json').write_text(json.dumps(out,indent=2)+'\n')

def plane_ambient():
    import sympy as s
    from sympy.polys.matrices import DomainMatrix
    from classification_three_hyperplanes_quadric_independent import ambient_matrices
    from classification_rank3_mixed_component_review import minor_certificate
    start=time.monotonic();signal.alarm(90)
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json'
    data=json.loads(path.read_text());xx=s.symbols('h x y z w');mm=s.symbols('m0:6')
    records=[]
    for c in data['cases']:
        if c['name'].startswith('shared_tangent'):continue
        P=s.zeros(5)
        for (i,j),f in zip(PAIRS,c['brackets']):P[i+1,j+1]=s.sympify(f);P[j+1,i+1]=-P[i+1,j+1]
        scale=s.ilcm(1,*[s.denom(v) for f in P for v in s.Poly(f,*xx,*mm).coeffs()]);P=P*scale
        witness=dict(zip(mm,[1,2,3,5,7,11]));point=P.subs(witness);directions=[P.diff(v) for v in mm]
        D,V,G,F=ambient_matrices(point,xx,directions)
        matrices={'D':D,'vertical':V,'orbit_scale':G,'full_source':F}
        ranks={key:len(DomainMatrix.from_Matrix(A).convert_to(s.QQ).rref()[1]) for key,A in matrices.items()}
        certificates={key:minor_certificate(A,ranks[key]) for key,A in matrices.items()}
        encode=lambda A:{str(ij):str(A[ij]) for ij in itertools.combinations(range(5),2) if A[ij]}
        records.append({'name':c['name'],'integer_scale':str(scale),'parameter_values':{str(k):str(v) for k,v in witness.items()},
            'full_six_magnetic_source':encode(P),'point_brackets':encode(point),
            'all6_parameter_directions':[encode(R) for R in directions],
            'ranks_QQ':ranks,'certificates':certificates,
            'actual_dimension_lower':ranks['full_source']-25,'ambient_tangent_upper':125-ranks['D'],
            'true_normal_dimension':150-ranks['D']-ranks['full_source'],
            'smooth_component_squeeze':ranks['D']+ranks['full_source']==150})
        print(c['name'],ranks,'normals',records[-1]['true_normal_dimension'],flush=True)
    assert len(records)==5
    inputs=[path,root/'code/classification_three_hyperplanes_quadric_independent.py',root/'code/classification_rank3_mixed_component_review.py']
    out={'task_id':'E_reduced_nonnormal','cases':records,'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start,
        'scope':'All150-variable projective Jacobi tangents, full vertical/orbit/actual6-parameter source directions at one declared point in each of five unresolved plane sources. Rank lower/upper bounds and true normals only; no component claim absent a squeeze or transverse obstruction.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json').write_text(json.dumps(out,indent=2)+'\n')

def plane_sources():
    """All six magnetic coefficients over every free plane-conic configuration."""
    import sympy as s
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,1024*1024**2))
    here=Path(__file__);root=here.resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS_20261107.json'
    data=json.loads(path.read_text());old=s.symbols('X Y Z');u=s.symbols('u')
    h,x,y,z,w=s.symbols('h x y z w');base=(h,x,y);xx=(x,y,z,w)
    names=dict(zip(('x','y','z'),old))|{'u':u};mm=s.symbols('m0:6')
    magnetic=sum(a*b for a,b in zip(mm,[h*h,h*x,h*y,x*x,x*y,y*y]))
    records=[]
    for c in data['cases']:
        if not c['has_free_degree22_basis']:continue
        H=s.sympify(c['hyperplane'],locals=names)
        rows=[[s.diff(H,v) for v in old]]
        for i in range(3):
            row=[int(i==j) for j in range(3)]
            if s.Matrix(rows+[row]).rank()>len(rows):rows.append(row)
        change=s.Matrix(rows);assert change.rows==3 and change.det()!=0
        inverse=change.inv()*s.Matrix(base);sub=dict(zip(old,inverse))
        anchor=[]
        for values in c['all_quadratic_gradient_syzygies']:
            delta=s.Matrix([s.sympify(v,locals=names) for v in values])
            multiple=s.cancel(sum(a*b for a,b in zip(delta,rows[0]))/H)
            assert s.Poly(multiple,*old).total_degree()<=1
            horizontal=delta-multiple*s.Matrix(old)
            v=(change*horizontal).applyfunc(lambda f:s.factor(f.subs(sub,simultaneous=True)))
            assert v[0]==0;anchor.append(v)
        assert len(anchor)==2;U,V=anchor
        bracket=s.Matrix([s.expand(sum(U[i]*s.diff(V[j],base[i])-V[i]*s.diff(U[j],base[i]) for i in range(3))) for j in range(3)])
        determinant=s.factor(U[1]*V[2]-U[2]*V[1]);assert determinant
        aa=s.factor(s.cancel((bracket[1]*V[2]-bracket[2]*V[1])/determinant))
        bb=s.factor(s.cancel((U[1]*bracket[2]-U[2]*bracket[1])/determinant))
        assert s.Poly(aa,*base).total_degree()<=1 and s.Poly(bb,*base).total_degree()<=1
        assert all(s.cancel(v)==0 for v in bracket-aa*U-bb*V)
        entries=[s.Integer(0),U[1],V[1],U[2],V[2],-aa*z-bb*w+magnetic]
        matrix=s.zeros(4)
        for (i,j),f in zip(PAIRS,entries):matrix[i,j]=f;matrix[j,i]=-f
        jac=[s.factor(sum(matrix[a,t]*s.diff(matrix[b,d],xx[t]) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
        assert not any(jac)
        q1=s.factor(s.sympify(c['q1'],locals=names).subs(sub,simultaneous=True))
        q2=s.factor(s.sympify(c['q2'],locals=names).subs(sub,simultaneous=True))
        pf=s.factor(entries[0]*entries[5]-entries[1]*entries[4]+entries[2]*entries[3])
        kappa=s.factor(s.cancel(pf/(q1*q2)));assert kappa and not any(s.diff(kappa,t) for t in base)
        records.append({'name':c['name'],'residue':c['residue'],'old_to_hxy_matrix':change.tolist(),
            'q1':str(q1),'q2':str(q2),'all_horizontal_anchor_U':list(map(str,U)),
            'all_horizontal_anchor_V':list(map(str,V)),
            'bracket_coefficients':list(map(str,[aa,bb])),
            'all_six_magnetic_coefficients':['m'+str(i) for i in range(6)],
            'brackets':list(map(str,entries)),'all_four_Jacobi':list(map(str,jac)),
            'full_Pfaffian':str(pf),'Pfaffian_multiple_of_q1q2':str(kappa)})
    out={'task_id':'E_reduced_nonnormal','cases':records,
        'input_sha256':{str(path.relative_to(root)):hashlib.sha256(path.read_bytes()).hexdigest()},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),
        'scope':'Every free plane H+two-conics configuration including its valid pivot exception, with all six magnetic coefficients, full Jacobi and actual factorization. Geometric completeness and component membership are separate proofs.',
        'seconds':time.monotonic()-start}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps({'all_passed':True,'cases':len(records),'seconds':out['seconds']}))

def plane_conics():
    """All quadratic gradient syzygies and Saito minors for contact candidates."""
    import sympy as s
    from sympy.polys.matrices import DomainMatrix
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    x,y,z,u=s.symbols('x y z u');xx=(x,y,z);base=x*z-y*y
    types=[('shared_tangent_generic',x,u*x*z-y*y+x*x),
           ('shared_tangent_contact3',x,base+x*y),
           ('shared_tangent_contact4',x,base+x*x),
           ('contact31_common_bitangent',x+8*y+16*z,base+x*y),
           ('contact4_other_tangent',z,base+x*x),
           ('contact31_tangent_at_simple_intersection',z,base+x*y),
           ('contact4_secant_at_contact',y,base+x*x),
           ('contact22_chord',y,x*z-u*y*y),
           ('contact31_chord',y,base+x*y)]
    mon2=[xx[i]*xx[j] for i in range(3) for j in range(i,3)]
    cc=s.symbols('a0:18');records=[]
    queue=[(name,H,q2,None) for name,H,q2 in types]
    for name,H,q2,value in queue:
        if value is not None:q2=q2.subs(u,value)
        field=s.QQ.frac_field(u) if value is None else s.QQ
        D=s.expand(H*base*q2);grad=[s.diff(D,v) for v in xx]
        delta=s.Matrix([sum(cc[6*i+j]*mon2[j] for j in range(6)) for i in range(3)])
        M,_=s.linear_eq_to_matrix(s.Poly(sum(a*b for a,b in zip(delta,grad)),*xx).coeffs(),cc)
        dm=DomainMatrix.from_Matrix(M).convert_to(field);rr,pc=dm.rref();rr=rr.to_Matrix()
        _,pr=DomainMatrix.from_Matrix(M.T).convert_to(field).rref()
        free=[i for i in range(18) if i not in pc];K=s.zeros(18,len(free))
        for j,k in enumerate(free):
            K[k,j]=1
            for i,p in enumerate(pc):K[p,j]=-rr[i,k]
        assert all(s.cancel(v)==0 for v in M*K)
        fields=[s.Matrix([sum(K[6*i+k,j]*mon2[k] for k in range(6)) for i in range(3)]) for j in range(K.cols)]
        determinants=[]
        for i,j in itertools.combinations(range(K.cols),2):
            ratio=s.factor(s.cancel(s.Matrix.hstack(s.Matrix(xx),fields[i],fields[j]).det()/D))
            assert not any(s.diff(ratio,v) for v in xx)
            determinants.append({'kernel_columns':[i,j],'Saito_determinant_divided_by_D':str(ratio)})
        pivot=s.factor(M.extract(pr,pc).det(method='domain-ge'))
        bad=s.factor(pivot*s.prod(set(s.denom(s.cancel(v)) for v in K)))
        for item in determinants:
            ratio=s.sympify(item['Saito_determinant_divided_by_D'])
            if ratio:
                numerator,denominator=s.fraction(ratio);bad=s.factor(bad*numerator*denominator)
        exceptional=s.factor_list(s.fraction(bad)[0],u)[1]
        if value is None:
            for factor,_ in exceptional:
                assert s.degree(factor,u)==1
                root=s.solve(factor,u)[0]
                if root not in (0,1):queue.append((name,H,q2,root))
        records.append({'name':name,'residue':str(value) if value is not None else 'generic u','hyperplane':str(H),'q1':str(base),'q2':str(q2),'plane_quintic':str(s.factor(D)),
            'quadratic_monomials':list(map(str,mon2)),'gradient_syzygy_matrix':M.tolist(),
            'rank':len(pc),'pivot_rows':list(pr),'pivot_columns':list(pc),'pivot_determinant':str(pivot),
            'complete_kernel':K.tolist(),'kernel_dimension':K.cols,'all_quadratic_gradient_syzygies':[list(map(str,v)) for v in fields],
            'Saito_determinants':determinants,'has_free_degree22_basis':any(v['Saito_determinant_divided_by_D']!='0' for v in determinants),
            'all_pivots_and_kernel_denominators':str(bad),'exceptional_parameter_factors':[(str(a),int(b)) for a,b in exceptional]})
    excluded=[]
    for name,H,q2,linear in [
        ('shared_tangent_contact22',x,x*z-u*y*y,[4*x,-y,-6*z]),
        ('shared_tangent_contact4',x,base+x*x,[0,x,2*y]),
        ('contact22_chord',y,x*z-u*y*y,[x,0,-z])]:
        D=H*base*q2
        assert s.expand(sum(a*s.diff(D,b) for a,b in zip(linear,xx)))==0
        excluded.append({'name':name,'hyperplane':str(H),'q1':str(base),'q2':str(q2),
            'nonzero_linear_gradient_syzygy':list(map(str,linear)),
            'obstruction_to_free_degrees22':'A free basis of degrees1,2,2 allows only Euler in degree1; a nonzero degree1 gradient syzygy cannot be an Euler multiple.'})
    here=Path(__file__);out={'task_id':'E_reduced_nonnormal','cases':records,'analytic_exclusions':excluded,
        'scope':'Every quadratic gradient syzygy and Saito determinant for the nine contact-budget candidate configurations. Completeness of the geometric candidate list, parameter exceptions and full magnetic-source/component recognition require separate proof.',
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),'seconds':time.monotonic()-start}
    (here.resolve().parents[1]/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS_20261107.json').write_text(json.dumps(out,indent=2,default=str)+'\n')
    print(json.dumps({'cases':[(c['name'],c['kernel_dimension'],c['has_free_degree22_basis'],c['exceptional_parameter_factors']) for c in records],'seconds':out['seconds']}))


def main():
    if '--mixed-normal-support' in sys.argv:return mixed_normal_support()
    if '--mixed-normal-symmetries' in sys.argv:return mixed_normal_symmetries()
    if '--nonzero-pair-cofactor' in sys.argv:return nonzero_pair_cofactor()
    if '--plane-even' in sys.argv:return plane_even()
    if '--plane-graphs' in sys.argv:return plane_graphs()
    if '--nil22-rank1-support' in sys.argv:return nil22_rank1_support()
    if '--nil22-base' in sys.argv:return nil22_base()
    if '--plane-cubic-support' in sys.argv:return plane_cubic_support()
    if '--plane-normals' in sys.argv:return plane_normals()
    if '--plane-ambient' in sys.argv:return plane_ambient()
    if '--plane-sources' in sys.argv:return plane_sources()
    if '--plane-conics' in sys.argv:return plane_conics()
    resource.setrlimit(resource.RLIMIT_CPU,(25,26))
    resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    started=time.monotonic();R=Ring(10)
    C,D,A,B,E,F,p,q,rr,inv=[R.variable('n'+str(i)) for i in range(10)]
    x,y,z,w=[R.variable(t) for t in ('x','y','z','w')]
    zero={};u,v=R.mul(x,y),R.mul(z,w)
    def sub(f,index,g):
        out={}
        for exp,coef in f.items():
            ee=list(exp);power=ee[index];ee[index]=0;term={tuple(ee):coef}
            for _ in range(power):term=R.mul(term,g)
            out=R.add(out,term)
        return out
    def subs(f,replacements):
        for i,g in replacements:f=sub(f,i,g)
        return f
    def unit_reduce(f,index):
        out={}
        for exp,coef in f.items():
            ee=list(exp);power=min(ee[index],ee[9]);ee[index]-=power;ee[9]-=power
            out[tuple(ee)]=out.get(tuple(ee),0)+coef
        return {e:c for e,c in out.items() if c}
    def eq(a,b):assert not R.add(a,R.scale(b,-1))
    def pf(entries):return R.add(R.mul(entries[0],entries[5]),R.scale(R.mul(entries[1],entries[4]),-1),R.mul(entries[2],entries[3]))
    ss=R.scale(R.add(p,q,rr),-1)
    source=[R.add(C,R.mul(A,u),R.mul(B,v)),R.mul(p,x,z),R.mul(q,x,w),R.mul(rr,y,z),R.mul(ss,y,w),R.add(D,R.mul(E,u),R.mul(F,v))]
    P=[[{} for _ in range(4)] for _ in range(4)]
    for (i,j),f in zip(PAIRS,source):P[i][j]=f;P[j][i]=R.scale(f,-1)
    jac=[R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in itertools.combinations(range(4),3)]
    v1,v2=R.add(p,q),R.add(p,rr)
    constraints=[R.mul(B,E),R.add(R.mul(v2,C),R.mul(B,D)),R.add(R.mul(E,C),R.scale(R.mul(v1,D),-1)),R.mul(B,R.add(v2,F)),R.mul(E,R.add(A,R.scale(v1,-1)))]
    def normalize(f):
        if not f:return ()
        a=f[min(f)];return tuple(sorted((e,c/a) for e,c in f.items()))
    coeffs=[]
    for f in jac:
        grouped={}
        for exp,c in f.items():grouped.setdefault(exp[-4:],{})[exp[:-4]+(0,)*4]=c
        coeffs+=list(grouped.values())
    assert {normalize(f) for f in coeffs}=={normalize(f) for f in constraints}
    totalpf=pf(source)
    eq(totalpf,R.add(R.mul(source[0],source[5]),R.mul(v1,v2,u,v)))
    branches=[]
    # B=E=0: h-simple means the whole highest Pfaffian coefficient is nonzero.
    b0=[(3,zero),(4,zero)];s0=[subs(f,b0) for f in source]
    kappa0=R.add(R.mul(A,F),R.mul(v1,v2))
    def highest(f):return {e:c for e,c in f.items() if sum(e[-4:])==2}
    eq(pf([highest(f) for f in s0]),R.mul(kappa0,u,v))
    branches.append({'stratum':'B=E=0','source':list(map(R.encode,s0)),
        'h_simple_implies_kappa_nonzero':True,'cocharacter_hxyzw':[1,0,0,0,0],
        'logcanonical_pfaffian_coefficient':R.encode(kappa0)})
    # B!=0, C!=0: irreducibility/rank4 imply v2!=0, v1=0 and A!=0.
    brep=[(4,zero),(7,R.scale(p,-1)),(1,R.scale(R.mul(v2,C,inv),-1)),(5,R.scale(v2,-1))]
    sb=[unit_reduce(subs(f,brep),3) for f in source]
    assert all(not unit_reduce(subs(f,brep),3) for f in jac)
    left=unit_reduce(R.mul(B,pf(sb)),3)
    right=R.scale(R.mul(v2,R.add(C,R.mul(A,u),R.mul(B,v)),R.add(C,R.mul(B,v))),-1)
    eq(left,right)
    # E!=0, D!=0: symmetric full source, reconstructed directly.
    erep=[(3,zero),(8,R.scale(p,-1)),(0,R.mul(v1,D,inv)),(2,v1)]
    se=[unit_reduce(subs(f,erep),4) for f in source]
    assert all(not unit_reduce(subs(f,erep),4) for f in jac)
    eq(unit_reduce(R.mul(E,pf(se)),4),R.mul(v1,R.add(D,R.mul(E,u)),R.add(D,R.mul(E,u),R.mul(F,v))))
    for label,entries,weights,kappa,nonzero in [
        ('B!=0',sb,[3,0,0,1,1],R.scale(R.mul(A,v2),-1),['B','C','A','p+r']),
        ('E!=0',se,[3,1,1,0,0],R.mul(F,v1),['E','D','F','p+q'])]:
        limit=[];allweights=[]
        for (i,j),f in zip(PAIRS,entries):
            kept={}
            for exp,c in f.items():
                degree=sum(exp[-4:]);assert degree in (0,2)
                # h^2 times the constant term, all derivative directions affine.
                weight=(2-degree)*weights[0]+sum(a*b for a,b in zip(exp[-4:],weights[1:]))-weights[i+1]-weights[j+1]
                assert weight>=0;allweights.append(weight)
                if weight==0:kept[exp]=c
            limit.append(kept)
        eq(pf(limit),R.mul(kappa,u,v))
        branches.append({'stratum':label,'source':list(map(R.encode,entries)),
            'inverted_parameter':'B' if label=='B!=0' else 'E',
            'irreducible_distinct_quadric_condition_forces_nonzero':nonzero,
            'cocharacter_hxyzw':weights,'all_monomial_weights':allweights,
            'logcanonical_limit':list(map(R.encode,limit)),
            'logcanonical_pfaffian_coefficient':R.encode(kappa),
            'full_Jacobi_and_Pfaffian_factorization_identity':True})
    here=Path(__file__);helper=Path(__import__('classification_quadratic_cofactor_review').__file__)
    result={'task_id':'E_reduced_nonnormal','all_passed':True,
        'parameter_order':['C','D','A','B','E','F','p','q','r','inverse_B_or_E'],
        'source':list(map(R.encode,source)),'all_four_Jacobi':list(map(R.encode,jac)),
        'exact_five_constraint_ideal':list(map(R.encode,constraints)),
        'full_Pfaffian':R.encode(totalpf),'branches':branches,
        'scope':'All strata of the complete generic pair-zero semisimple weight source with actual divisor h times two distinct irreducible quadrics. Extra resonance spectra, nilpotent modular fields and pure-homogeneous sources are separate obligations.',
        'method':'Exact standard-library Fraction sparse polynomial identities; full source and all Jacobi coefficients, inverse equations B*inv=1 or E*inv=1, and every cocharacter monomial retained.',
        'input_sha256':{str(helper.relative_to(here.resolve().parents[1])):hashlib.sha256(helper.read_bytes()).hexdigest()},
        'producer_sha256':hashlib.sha256(here.read_bytes()).hexdigest(),'producer_bytes':here.read_text(),
        'seconds':time.monotonic()-started}
    (here.resolve().parents[1]/'reports/CLASSIFICATION_TWO_QUADRICS_PAIR_ZERO_20261107.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'branches':len(branches),'seconds':result['seconds']}))

if __name__=='__main__':main()
