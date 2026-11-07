#!/usr/bin/env python3
"""Independent Fraction(u) reconstruction of all plane-conic Saito kernels."""
import hashlib,json,resource,time,sys
from pathlib import Path
import classification_unimodular_logarithmic_review as r


def two_lines_rank2_translation():
    """All original21 coefficients and whole nonplane support/normalizers."""
    from fractions import Fraction as F
    from itertools import combinations
    import re
    from classification_quadratic_cofactor_review import Ring,PAIRS
    import classification_quadratic_cofactor_review as helper
    import classification_two_lines_plane_cubic_component_review as q
    resource.setrlimit(resource.RLIMIT_CPU,(35,36));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    inputs=[root/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json',root/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_RANK2_TRANSLATION_SUPPORT_20261107.json']
    source=json.loads(inputs[0].read_text())['cases'][-1];assert source['case']=='j111_rank2'
    def rename(text):return re.sub(r'\b(h|u|v)\b',lambda m:{'h':'1','u':'y','v':'z'}[m[0]],text)
    def matrix(R,entries):
        M=[[{} for j in range(4)] for i in range(4)]
        for (i,j),v in zip(PAIRS,entries):M[i][j]=v;M[j][i]=R.scale(v,-1)
        return M
    def jac(R,E):
        P=matrix(R,E)
        return [R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in combinations(range(4),3)]
    def pf(R,E):return R.add(R.mul(E[0],E[5]),R.scale(R.mul(E[1],E[4]),-1),R.mul(E[2],E[3]))
    R=Ring(21);parse=lambda text:R.parse(rename(text));E=[parse(source['brackets'][str(i+1)+','+str(j+1)]) for i,j in PAIRS]
    T=R.parse('y-x*z');a,b,c=map(R.variable,['n4','n11','n18']);L=R.add(b,R.scale(R.mul(c,R.variable('x')),-1))
    A=R.add(E[3],R.scale(R.variable('y'),-1),R.mul(a,T));B=R.add(E[4],R.mul(b,T));C=R.add(E[5],R.mul(c,T));G=R.add(R.mul(R.variable('x'),C),R.scale(B,-1))
    assert all(not e[22] and not e[23] for poly in [A,B,C] for e in poly)
    cubic=R.add(R.mul(L,T),G);assert pf(R,E)==R.mul(R.variable('x'),cubic)
    I=R.add(R.mul(R.variable('x'),R.derivative(L,0)),R.mul(c,R.derivative(G,3)),R.mul(R.derivative(C,3),L))
    II=R.add(R.mul(R.variable('x'),R.derivative(G,0)),R.mul(A,L),R.scale(R.mul(C,R.derivative(G,3)),-1),
        R.scale(R.mul(R.add(R.constant(1),R.scale(a,-1),R.scale(R.derivative(C,3),-1)),G),-1))
    JJ=jac(R,E);assert not any(JJ[:3]);assert JJ[3]==R.scale(R.add(R.mul(I,T),II),-1)
    records=[]
    for packet in json.loads(inputs[1].read_text())['cases']:
        R=Ring(8);parse=lambda text:R.parse(rename(text));E=list(map(parse,packet['full_brackets']))
        assert not any(jac(R,E));assert pf(R,E)==parse(packet['full_Pfaffian'])==R.mul(R.variable('x'),parse(packet['full_irreducible_cubic']))
        A,B,C=map(parse,[packet['A'],packet['B'],packet['C']]);a,b,c=map(parse,[packet['a'],packet['b'],packet['c']]);T=R.parse('y-x*z');L=R.add(b,R.scale(R.mul(c,R.variable('x')),-1));G=R.add(R.mul(R.variable('x'),C),R.scale(B,-1))
        assert not R.add(R.mul(R.variable('x'),R.derivative(L,0)),R.mul(c,R.derivative(G,3)),R.mul(R.derivative(C,3),L))
        assert not R.add(R.mul(R.variable('x'),R.derivative(G,0)),R.mul(A,L),R.scale(R.mul(C,R.derivative(G,3)),-1),R.scale(R.mul(R.add(R.constant(1),R.scale(a,-1),R.scale(R.derivative(C,3),-1)),G),-1))
        q.N=5
        def field(*values):return list(map(q.parse,values))
        nil=[field('0','0',j,'0','0') for j in ['h','x','w']]+[field('0','0','0',j,'0') for j in ['h','x','w']]+[field('0','0','0','0',j) for j in ['h','x']]
        D1=field('0','x','y','0','0');D2=field('0','0','y','z','w');N=field('0','0','x','h','0')
        def mat(Y):return [[q.der(Y[i],j).get((0,)*5,F(0)) for j in range(5)] for i in range(5)]
        def product(A,B):return [[sum(A[i][k]*B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
        def comm(A,B):return [[a-b for a,b in zip(r,s)] for r,s in zip(product(A,B),product(B,A))]
        flat=lambda M:[v for row in M for v in row]
        NN=mat(N);nil=list(map(mat,nil));K=nil+[mat(D2)]+([mat(D1)] if packet['name']=='L_h' else []);vectors=list(map(flat,K));dim=len(K)
        assert q.rank(q.transpose(vectors))[0]==dim==packet['preserving_projective_normalizer_dimension']
        assert q.rank(q.transpose(vectors+[flat(comm(A,B)) for A,B in combinations(K,2)]))[0]==dim
        assert all(not any(flat(comm(A,NN))) for A in nil)
        assert comm(mat(D2),NN)==NN and not any(flat(comm(mat(D1),NN)))
        order=[2,3,4,0,1];assert all(not A[order[i]][order[j]] for A in nil for i in range(5) for j in range(i+1))
        upper=24+len(packet['free_parameters'])-dim;assert upper==22
        records.append({'name':packet['name'],'all_passed':True,'all_full_source_Jacobi_Pfaffian_reconstructed':True,
            'parameter_count':len(packet['free_parameters']),'nilpotent_group_dimension':8,'projective_normalizer_dimension':dim,
            'all_Lie_closure_and_cofactor_line_conditions_reconstructed':True,'whole_projective_PGL_image_dimension_upper':upper})
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'original21_parameter_source_I_II_identity_reconstructed':True,'cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs+[Path(helper.__file__),Path(q.__file__)]},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'Independent original full21-coefficient source/normal-form Jacobi identities, both complete nonplane support formulas and exact preserving normalizer matrices. Necessary polynomial-division reductions, h/x exchange and free-action dimension implications are separately written proofs. No E component membership claimed.'}
    (root/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_RANK2_TRANSLATION_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))

def reduced_divisor_normalizers():
    """Whole-source dimension groups; no tangent-point inference."""
    from fractions import Fraction as F
    from itertools import combinations
    import classification_two_lines_plane_cubic_component_review as q
    resource.setrlimit(resource.RLIMIT_CPU,(25,26));resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1];q.N=5;records=[]
    names=['h','x','y','z','w'];h,x,u,v,w=names
    def field(*data):return list(map(q.parse,data))
    common=[field('0','0','0','0',h),field('0','0','0','0',x),field('0','0',h,'0','0'),field('0','0',x,'0','0')]
    cases=[]
    for name,L,M,dim in [('p0','h+x','h',5),('h_x','h','x',6)]:
        V=field('0','0',v,L,M);K=[field('0','0',v,L,'0'),field('0','0',w,M,'0'),*common]
        cases.append((name,V,K,[2,3,4,0,1],dim))
    V=field('0','0',v,'0','h+x');K=[field('0','0',h,'0','0'),field('0','0',x,'0','0'),field('0','0',v,'0','0'),
        field('0','0',w,'h+x','0'),field('0','0','0','0',h),field('0','0','0','0',x),field('0','0','0','0',v)]
    cases.append(('translation_n6_one',V,K,[2,4,3,0,1],6))
    def matrix(Y):return [[q.der(Y[i],j).get((0,)*5,F(0)) for j in range(5)] for i in range(5)]
    def product(A,B):return [[sum(A[i][k]*B[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
    def bracket(A,B):return [[a-b for a,b in zip(r,s)] for r,s in zip(product(A,B),product(B,A))]
    def flat(A):return [v for row in A for v in row]
    for name,V,K,order,dim in cases:
        N=matrix(V);AA=list(map(matrix,K));vectors=list(map(flat,AA));rank=q.rank(q.transpose(vectors))[0]
        assert rank==len(K)
        assert all(not any(flat(bracket(A,N))) for A in AA)
        assert all(not any(A[0]) and not any(A[1]) for A in AA)
        assert all(not A[order[i]][order[j]] for A in AA for i in range(5) for j in range(i+1))
        BB=[bracket(A,B) for A,B in combinations(AA,2)]
        assert q.rank(q.transpose(vectors+list(map(flat,BB))))[0]==rank
        # Simultaneously strict upper triangular matrices integrate to a connected
        # algebraic unipotent group by the characteristic-zero exponential map.
        upper=24+dim-rank
        assert upper<=24
        records.append({'name':name,'coordinate_order':['h','x','u','v','w'],
            'cofactor_matrix':[[str(v) for v in row] for row in N],
            'all_generator_matrices':[[[str(v) for v in row] for row in A] for A in AA],
            'linearly_independent_generator_count':rank,'all_commute_with_V':True,'all_fix_h_and_x_pointwise':True,
            'all_strict_upper_triangular_in_order':order,'Lie_bracket_closure':True,
            'connected_algebraic_unipotent_group_dimension':rank,'fixed_complete_source_dimension':dim,
            'preserves_u_coefficient_and_determinant_one':True,
            'whole_projective_PGL_image_dimension_upper':upper,
            'whole_source_preservation_reason':'The complete fixed-V Jacobi support has three disjoint irreducible five-parameter components; a connected group preserves each.' if name=='p0' else
                'The full irreducible-cubic support is contained in the unique complete six-parameter source, and the group preserves V/h/x and irreducibility.' if name=='h_x' else
                'The group preserves V/h/x, determinant1 and the u coefficient of F; therefore it fixes n6=1 and preserves the unique complete six-parameter source.'})
        print(name,'WHOLE NORMALIZER PASS',rank,'dimension upper',upper,flush=True)
    paths=[root/'research/classification/TWO_HYPERPLANES_CUBIC_SUPPORT.md',
        root/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_20261107.json',
        root/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_COMPONENTS_20261107.json',
        root/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_TRANSLATION_COMPONENT_20261107.json',Path(q.__file__)]
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'Independent Fraction verification of exact algebraic unipotent groups preserving the entire complete fixed-cofactor sources. Completeness of source classification and the free-action dimension argument are separate written proof inputs. No point-orbit-rank, higher-jet or E component-membership inference.'}
    (root/'reports/CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NORMALIZER_DIMENSION_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')

def nilpotent_factor_reductions():
    """New HQQ bridges; earlier complete nonlinear reductions stay explicit inputs."""
    from fractions import Fraction as F
    from itertools import combinations,combinations_with_replacement
    from classification_quadratic_cofactor_review import Ring,PAIRS
    import classification_quadratic_cofactor_review as helper
    import classification_two_lines_plane_cubic_component_review as q
    resource.setrlimit(resource.RLIMIT_CPU,(35,36));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1];records=[];paths=[];q.N=5
    configs=[('nil31_long',['0','y','z','h','0'],['h*h','h*w','w*w','2*h*y-z*z']),
        ('nil211_short',['0','y','0','h','0'],['h*h','h*y','h*w','y*y','y*w','w*w','h*x-y*z']),
        ('nil211_long',['0','y','h','0','0'],['h*h','h*z','h*w','z*z','z*w','w*w','h*x-y*y/2'])]
    mon=q.monomials(2)
    for name,field,basis in configs:
        field=list(map(q.parse,field));basis=list(map(q.parse,basis));cols=[]
        for e in mon:
            p={e:F(1)};cols.append(q.add(*(q.mul(field[j],q.der(p,j)) for j in range(5))))
        M=[[p.get(e,F(0)) for p in cols] for e in mon];rank=q.rank(M)[0]
        assert rank+len(basis)==15 and q.rank([[p.get(e,F(0)) for p in basis] for e in mon])[0]==len(basis)
        assert all(not q.add(*(q.mul(field[j],q.der(p,j)) for j in range(5))) for p in basis)
        records.append({'name':name,'complete_quadratic_invariant_kernel':True,'rank':rank,'dimension':len(basis)})
    path=root/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json';paths.append(path)
    full=json.loads(path.read_text())['full_nilpotent_sources']
    def matrix(R,entries):
        M=[[{} for j in range(4)] for i in range(4)]
        for (i,j),v in zip(PAIRS,entries):M[i][j]=v;M[j][i]=R.scale(v,-1)
        return M
    def pf(R,v):return R.add(R.mul(v[0],v[5]),R.scale(R.mul(v[1],v[4]),-1),R.mul(v[2],v[3]))
    def substitute(R,poly,values):
        out={}
        for e,c in poly.items():
            p=R.constant(c)
            for i,d in enumerate(e):
                for k in range(d):p=R.mul(p,values.get(R.names[i],R.variable(R.names[i])))
            out=R.add(out,p)
        return out
    # Type1 canonical HQQ stratum: a full source identity, rather than a sample.
    R=Ring(13);values={f'n{i}':{} for i in [1,2,5,6,7,9,10,11]};values['n8']=R.constant(-F(1,4));values['n12']=R.constant(-F(1,2))
    E=[substitute(R,R.parse(v),values) for v in full[0]['brackets']];P=matrix(R,E)
    assert all(not R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in combinations(range(4),3))
    assert pf(R,E)==R.parse('(2*w*y-z*z)**2/48-n0*w*w/2')
    records.append({'name':'nil4_c0_HQQ_stratum','all_source_Jacobi_and_Pfaffian_identity':True,'earlier_complete_nonlinear_reduction_remains_separate_input':True})
    # Type3 exclusion holds even before Jacobi.
    R=Ring(19);E=list(map(R.parse,full[2]['brackets']));assert all(e[21]>=1 for e in pf(R,E))
    records.append({'name':'nil31_c0','unconditional_linear_factor_z':True})
    # Type9: independent full28-parameter nonlinear coordinate push and potential.
    R=Ring(28);E=list(map(R.parse,full[8]['brackets']));P=matrix(R,E);names=('x','y','w','z')
    coords=list(map(R.parse,['x-y*z','y','w','z']));inverse={'x':R.parse('x+y*z')}
    new=[]
    for a,b in PAIRS:
        v=R.add(*(R.mul(R.derivative(coords[a],i),R.derivative(coords[b],j),P[i][j]) for i in range(4) for j in range(4)))
        new.append(substitute(R,v,inverse))
    geometry=root/'reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_SHORT1_GEOMETRY_20261107.json';paths.append(geometry)
    data=json.loads(geometry.read_text())['nil_short_translation_geometry'];parse=lambda v:R.parse(v.replace('u','x'))
    assert new==list(map(parse,data['transformed_brackets']))
    Phi=parse(data['potential']);V=list(map(parse,data['vector_field']));axes=[0,1,3]
    assert [new[0],new[1],new[3]]==[R.derivative(Phi,3),R.scale(R.derivative(Phi,1),-1),R.derivative(Phi,0)]
    assert [new[2],new[4],new[5]]==[R.scale(v,-1) for v in V]
    assert R.add(*(R.derivative(v,a) for v,a in zip(V,axes)))==R.constant(-1)
    VPhi=R.add(*(R.mul(v,R.derivative(Phi,a)) for v,a in zip(V,axes)));Q=R.add(VPhi,Phi)
    NP=matrix(R,new);logical=[0,1,3,2]
    JJ=[R.add(*(R.mul(NP[a][t],R.derivative(NP[b][d],logical[t])) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in combinations(range(4),3)]
    assert JJ==[{},R.derivative(Q,3),R.scale(R.derivative(Q,1),-1),R.derivative(Q,0)]
    assert pf(R,new)==R.scale(VPhi,-1)
    u2={e[:28]+(0,0,0,0):c for e,c in Q.items() if e[28:]==(2,0,0,0)}
    assert u2==R.parse('-n20*n25')
    records.append({'name':'nil211_short','full28_parameter_coordinate_push_and_potential_reconstructed':True,
        'all_Jacobi_exactly_three_derivatives_of_VPhi_plus_Phi':True,'actual_u_squared_coefficient_minus_n20n25':True,
        'quadratic_coefficient_c_rank_at_most2_on_Jacobi_by_two_explicit_factor_branches':True})
    # Type10 whole three-parameter family and its three-dimensional group action.
    R=Ring(6);E=list(map(R.parse,['(x-y*y/2)/2','-y*z/2','-2*n1*z+n2*y','-z/2','n2','n0*(x-y*y/2)']));P=matrix(R,E)
    assert all(not R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in combinations(range(4),3))
    assert pf(R,E)==R.parse('n0*(x-y*y/2)**2/2+n1*z*z')
    # Use diagonal g,d and translation f; clear inverses by comparing pullbacks.
    coords=list(map(R.parse,['x','y','n3*z','n4*w+n5']));target_params={'n0':R.parse('n3*n4*n0'),'n1':R.parse('n4*n1'),'n2':R.parse('n4*n2')}
    # For a'=d*a/g, multiply only its affected bracket by g before comparison.
    for idx,(a,b) in enumerate(PAIRS):
        pushed=R.add(*(R.mul(R.derivative(coords[a],i),R.derivative(coords[b],j),P[i][j]) for i in range(4) for j in range(4)))
        if idx==2:expected=R.parse('-2*n4*n1*z+n4*n2*y')
        else:expected=substitute(R,substitute(R,E[idx],target_params),dict(zip(('x','y','z','w'),coords)))
        assert pushed==expected
    records.append({'name':'nil211_long','whole3_parameter_Jacobi_Pfaffian_and_3dim_normalizer_action_reconstructed':True,
        'delta_prime':'g*d*delta','a_prime':'d*a/g','kappa_prime':'d*kappa','projective_PGL_image_dimension_upper':24})
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[Path(helper.__file__),Path(q.__file__)]},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'New HQQ invariant-kernel/potential/full canonical-family arithmetic bridges. Earlier complete nonlinear nil4/nil31long reductions are explicit separate proof inputs; no claim that their ideal-membership calculations were newly independently reconstructed. Type7 remains separate.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PURE_NILPOTENT_FACTOR_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))

def mixed_normal_support():
    """Independent all-row sources and complete-support stratum arithmetic."""
    import ast
    from fractions import Fraction as F
    from itertools import combinations, combinations_with_replacement
    from classification_quadratic_cofactor_review import Ring,PAIRS
    import classification_quadratic_cofactor_review as helper
    import classification_two_lines_plane_cubic_component_review as q
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    paths=[root/'reports'/('CLASSIFICATION_TWO_QUADRICS_'+p+'_20261107.json') for p in ['MIXED_NORMAL_SYMMETRIES','MIXED_NORMAL_SUPPORT']]
    sources=json.loads(paths[0].read_text())['cases'];supports={v['name']:v for v in json.loads(paths[1].read_text())['cases']}
    R=Ring(13);xx=[R.variable(v) for v in ('x','y','z','w')];zero=F(0);records=[]
    # Only nonzero source parameters from the written support proof are inverted.
    def inverse_monomial(poly):
        assert len(poly)==1
        e,c=next(iter(poly.items()));assert not any(e[13:])
        return {tuple(-v for v in e):1/c}
    def parse(expression):
        def walk(node):
            if isinstance(node,ast.Constant) and type(node.value) is int:return R.constant(node.value)
            if isinstance(node,ast.Name):return R.variable(node.id)
            if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return R.scale(walk(node.operand),-1)
            if isinstance(node,ast.BinOp):
                a,b=walk(node.left),walk(node.right)
                if isinstance(node.op,ast.Add):return R.add(a,b)
                if isinstance(node.op,ast.Sub):return R.add(a,R.scale(b,-1))
                if isinstance(node.op,ast.Mult):return R.mul(a,b)
                if isinstance(node.op,ast.Div):return R.mul(a,inverse_monomial(b))
                if isinstance(node.op,ast.Pow):
                    assert set(b)<={R.zero};power=b.get(R.zero,zero);assert power.denominator==1
                    if power<0:a=inverse_monomial(a);power=-power
                    out=R.constant(1)
                    for k in range(int(power)):out=R.mul(out,a)
                    return out
            raise AssertionError(ast.dump(node))
        return walk(ast.parse(expression,mode='eval').body)
    def substitute(poly,mapping):
        out={}
        for e,c in poly.items():
            term=R.constant(c)
            for i,degree in enumerate(e):
                assert degree>=0
                value=mapping.get(R.names[i],R.variable(R.names[i]))
                for k in range(degree):term=R.mul(term,value)
            out=R.add(out,term)
        return out
    def mons(degree):
        out=[]
        for d in range(degree+1):
            for indexes in combinations_with_replacement(range(4),d):
                e=[0]*4
                for j in indexes:e[j]+=1
                out.append(tuple(e))
        return out
    mon2,mon1=mons(2),mons(1)
    def coefficients(poly,mon):
        assert all(not any(e[:13]) for e in poly)
        return [poly.get((0,)*13+e,zero) for e in mon]
    def matrix(entries):
        M=[[{} for j in range(4)] for i in range(4)]
        for (i,j),v in zip(PAIRS,entries):M[i][j]=v;M[j][i]=R.scale(v,-1)
        return M
    def pf(v):return R.add(R.mul(v[0],v[5]),R.scale(R.mul(v[1],v[4]),-1),R.mul(v[2],v[3]))
    for source in sources:
        begin=time.monotonic();name=source['name'];S=list(map(parse,source['semisimple_field']));N=list(map(parse,source['nilpotent_field']))
        def lhs(entries,ss,tt):
            P=matrix(entries);rows=[]
            for j in range(4):rows+=coefficients(R.add(*(R.derivative(P[i][j],i) for i in range(4)),R.scale(R.mul(ss,S[j]),-1),R.scale(R.mul(tt,N[j]),-1)),mon1)
            for field in [S,N]:
                for i,j in PAIRS:
                    lie=R.add(*(R.add(R.mul(field[k],R.derivative(P[i][j],k)),R.scale(R.mul(P[k][j],R.derivative(field[i],k)),-1),R.scale(R.mul(P[i][k],R.derivative(field[j],k)),-1)) for k in range(4)))
                    rows+=coefficients(lie,mon2)
            assert len(rows)==200
            return rows
        columns=[]
        for col in range(92):
            entries=[{} for j in range(6)];ss=tt={}
            if col<90:entries[col//15]={(0,)*13+mon2[col%15]:F(1)}
            elif col==90:ss=R.constant(1)
            else:tt=R.constant(1)
            columns.append(lhs(entries,ss,tt))
        rows=[{j:columns[j][i] for j in range(92) if columns[j][i]} for i in range(200)]
        assert rows==[{int(j):F(v) for j,v in row.items()} for row in source['full_constraint_sparse_rows']]
        coef=list(map(parse,source['all92_coefficients']));free=source['free_coefficient_indices'];count=len(free)
        assert len(coef)==92 and count==len(source['free_parameters'])==92-source['rank']
        assert [coef[j] for j in free]==[R.variable('n'+str(i)) for i in range(count)]
        assert all(not R.add(*(R.scale(coef[j],a) for j,a in row.items())) for row in rows)
        det=q.determinant([[rows[i].get(j,zero) for j in source['pivot_columns']] for i in source['pivot_rows']])
        assert det==F(source['pivot_determinant']) and det and len(source['pivot_columns'])==92-count
        entries=[R.add(*(R.mul(coef[15*i+j],{(0,)*13+e:F(1)}) for j,e in enumerate(mon2))) for i in range(6)]
        assert entries==list(map(parse,source['full_brackets']));P=matrix(entries)
        Jac=[R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in combinations(range(4),3)]
        assert Jac==list(map(parse,source['all_four_Jacobi'])) and pf(entries)==parse(source['full_Pfaffian'])
        branches=[]
        for branch in supports[name]['branches']:
            mapping={j:parse(v) for j,v in branch['necessary_substitutions'].items()}
            actual=[substitute(v,mapping) for v in entries];actual_pf=pf(actual)
            if 'full_Pfaffian_after_substitution' in branch:assert actual_pf==parse(branch['full_Pfaffian_after_substitution'])
            if 'full_brackets_after_substitution' in branch:
                assert actual==list(map(parse,branch['full_brackets_after_substitution']))
                assert all(not substitute(v,mapping) for v in Jac)
                free_names=set(branch['free_parameters']);assert len(free_names)==branch['source_dimension_upper']==branch['normalizer_dimension_lower']
                assert branch['full_PGL_projective_image_dimension_upper']==23
                assert all(all(not d or R.names[i] in free_names or i>=13 for i,d in enumerate(e)) for v in actual for e in v)
            elif name=='double_pair_two_nil' and branch['name']=='ss_coefficient_nonzero':
                assert all(e[14]>=1 and e[16]>=1 for e in actual_pf)
            branches.append({'name':branch['name'],'all_original_stratum_Jacobi_and_Pfaffian_arithmetic_passed':True,
                'source_dimension_upper':branch.get('source_dimension_upper'),'normalizer_dimension_lower':branch.get('normalizer_dimension_lower')})
        if name=='double_pair_one_nil':assert all(e[14]>=1 for e in pf(entries))
        records.append({'name':name,'all_passed':True,'full200_by92_matrix_reconstructed':True,'complete_kernel_and_rank_minor_reconstructed':True,
            'full_four_Jacobi_and_Pfaffian_reconstructed':True,'rank':source['rank'],'free_parameters':count,'branches':branches,'seconds':time.monotonic()-begin})
        print(name,'FULL200/92 AND LAURENT SUPPORT PASS',flush=True)
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[Path(helper.__file__),Path(q.__file__)]},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'Independent Fraction sparse matrix/full kernel/rank/Jacobi/Pfaffian and Laurent stratum arithmetic for all five supersources. The whole-support normal-Hamiltonian restrictions, finite normal-form completeness, and free normalizer-action dimension implication are separate proofs in TWO_QUADRICS_JOINT_SS_RANK1_SUPPORT.md. No E component recognition asserted.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SUPPORT_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))

def nonzero_pair_cofactor():
    """Independent all-row arithmetic and exhaustive joint-character matching."""
    from fractions import Fraction as F
    from itertools import combinations, combinations_with_replacement, permutations
    from classification_quadratic_cofactor_review import Ring,PAIRS,complete_solution
    import classification_two_lines_plane_cubic_component_review as q
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_NONZERO_PAIR_COFACTOR_20261107.json'
    saved=json.loads(path.read_text());records=[];R=Ring(1);xx=[R.variable(v) for v in ('x','y','z','w')]
    zero=F(0);e0=(0,)*4
    def mons(degree):
        out=[]
        for d in range(degree+1):
            for indexes in combinations_with_replacement(range(4),d):
                e=[0]*4
                for j in indexes:e[j]+=1
                out.append(tuple(e))
        return out
    mon2,mon3,mon1=mons(2),mons(3),mons(1)
    def coeffs(poly,mon):
        assert all(not e[0] for e in poly)
        return [poly.get((0,)+e,zero) for e in mon]
    def polynomial(e):return {(0,)+e:F(1)}
    def matrix(entries):
        M=[[{} for j in range(4)] for i in range(4)]
        for (i,j),v in zip(PAIRS,entries):M[i][j]=v;M[j][i]=R.scale(v,-1)
        return M
    expected={
        'rank44':('x+y*w','y+x*z',[(0,-1),(1,0),(1,1),(-1,-1)]),
        'rank43':('x+y*z','y*y+x*w',[(0,-1),(F(1,2),0),(-F(1,2),-1),(1,1)]),
        'rank33':('x*x+y*z','y*y+x*w',[(0,-F(1,2)),(F(1,2),0),(-F(1,2),-1),(1,F(1,2))])}
    for c in saved['cases']:
        begin=time.monotonic();name=c['name'];qt1,qt2,chars=expected[name];quads=list(map(R.parse,[qt1,qt2]));chars=[tuple(map(F,t)) for t in chars]
        assert chars==[tuple(map(F,t)) for t in c['joint_characters']]
        fields=[[R.scale(xx[j],chars[j][a]) for j in range(4)] for a in range(2)]
        assert fields==[list(map(R.parse,v)) for v in c['cofactor_fields']]
        def lhs(entries):
            P=matrix(entries);rows=[]
            for quad in quads:
                for i in range(4):rows+=coeffs(R.scale(R.add(*(R.mul(P[i][j],R.derivative(quad,j)) for j in range(4))),-1),mon3)
            for j in range(4):rows+=coeffs(R.add(*(R.derivative(P[i][j],i) for i in range(4))),mon1)
            for a in range(2):
                for i,j in PAIRS:
                    rows+=coeffs(R.add(*(R.mul(fields[a][k],R.derivative(P[i][j],k)) for k in range(4)),R.scale(P[i][j],-chars[i][a]-chars[j][a])),mon2)
            assert len(rows)==480
            return rows
        rhs=[]
        for quad,field in zip(quads,fields):
            for v in field:rhs+=coeffs(R.mul(quad,v),mon3)
        for j in range(4):rhs+=coeffs(R.add(fields[0][j],fields[1][j]),mon1)
        rhs += [zero]*180
        columns=[]
        for i in range(90):
            entries=[{} for j in range(6)];entries[i//15]=polynomial(mon2[i%15]);columns.append(lhs(entries))
        rows=[{j:columns[j][i] for j in range(90) if columns[j][i]} for i in range(480)]
        assert rows==[{int(j):F(v) for j,v in row.items()} for row in c['full_constraint_sparse_rows']]
        assert rhs==list(map(F,c['full_constraint_rhs'])) and c['full_constraint_shape']==[480,90]
        free,piv,solution=complete_solution(rows,rhs)
        assert free==c['free_coefficient_indices'] and piv==c['pivot_columns'] and len(free)==1 and len(piv)==89
        coefficients=list(map(R.coefficient,solution));assert coefficients==list(map(R.parse,c['all90_coefficients']))
        det=q.determinant([[rows[i].get(j,zero) for j in c['pivot_columns']] for i in c['pivot_rows']])
        assert det==F(c['pivot_determinant']) and det
        entries=[R.add(*(R.mul(coefficients[15*i+j],polynomial(e)) for j,e in enumerate(mon2))) for i in range(6)]
        assert entries==list(map(R.parse,c['full_brackets']));P=matrix(entries)
        jac=[R.add(*(R.mul(P[a][t],R.derivative(P[b][d],t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4))) for i,j,k in combinations(range(4),3)]
        assert not any(jac)
        pf=lambda v:R.add(R.mul(v[0],v[5]),R.scale(R.mul(v[1],v[4]),-1),R.mul(v[2],v[3]))
        kappa=R.parse(c['Pfaffian_multiple_of_q1q2']);assert kappa and all(not any(e[1:]) for e in kappa)
        assert pf(entries)==R.parse(c['full_Pfaffian'])==R.mul(kappa,*quads)
        weights=c['cocharacter_hxyzw'];limit=[]
        for (i,j),v in zip(PAIRS,entries):
            kept={}
            for e,a in v.items():
                degree=sum(e[1:]);weight=(2-degree)*weights[0]+sum(d*b for d,b in zip(e[1:],weights[1:]))-weights[i+1]-weights[j+1]
                assert weight>=0
                if not weight:kept[e]=a
            target=tuple(int(t in (i,j)) for t in range(4))
            assert all(e[1:]==target for e in kept);limit.append(kept)
        assert limit==list(map(R.parse,c['whole_source_logcanonical_limit']))
        assert pf(limit)==R.parse(c['limit_Pfaffian'])==R.mul(kappa,*xx)
        records.append({'name':name,'all_passed':True,'full480_by90_reconstructed':True,'complete_affine_source_reconstructed':True,
            'rank':89,'free_coefficients':1,'nonzero_rank_minor':str(det),'all_Jacobi_Pfaffian_and_actual_SNC_limit_reconstructed':True,'seconds':time.monotonic()-begin})
        print(name,'FULL480/90 FRACTION PASS',flush=True)
    # Every nondegenerate symmetric rank-r weight block has an involution matching
    # on r coordinate weights: a self-dual block gives fixed points, opposite
    # weight blocks pair. Enumerate all such matchings, with h coordinate zero.
    def involutions(items):
        if not items:yield [];return
        i,*tail=items
        for rest in involutions(tail):yield [(i,i)]+rest
        for j in tail:
            for rest in involutions([k for k in tail if k!=j]):yield [(i,j)]+rest
    def patterns(rank):return [v for subset in combinations(range(5),rank) for v in involutions(list(subset))]
    def solve(pattern1,pattern2):
        eq=[]
        for a,pattern in enumerate([pattern1,pattern2]):
            target=(F(0),F(-1)) if a==0 else (F(1),F(0))
            for i,j in pattern:
                row=[F(int(k==i)+int(k==j)) for k in range(1,5)]+list(target);eq.append(row)
        eq.append([F(1)]*4+[F(1),F(-1)])
        row=0;piv=[]
        for col in range(4):
            k=next((k for k in range(row,len(eq)) if eq[k][col]),None)
            if k is None:continue
            eq[row],eq[k]=eq[k],eq[row];d=eq[row][col];eq[row]=[v/d for v in eq[row]]
            for k in range(len(eq)):
                if k!=row and eq[k][col]:
                    d=eq[k][col];eq[k]=[a-d*b for a,b in zip(eq[k],eq[row])]
            piv.append(col);row+=1
        if any(not any(v[:4]) and any(v[4:]) for v in eq):return None
        assert len(piv)==4, ('unclassified positive-dimensional weights',pattern1,pattern2,eq)
        return [(F(0),F(0))]+[tuple(eq[piv.index(i)][4:]) for i in range(4)]
    enumeration=[]
    for ranks in [(4,4),(4,3),(3,3)]:
        count=0;orbits=set()
        for p1 in patterns(ranks[0]):
            for p2 in patterns(ranks[1]):
                weights=solve(p1,p2)
                if weights is None:continue
                count+=1;orbits.add(tuple(sorted(weights[1:])))
        target=tuple(sorted(expected['rank'+''.join(map(str,ranks))][2]))
        assert orbits=={target}
        enumeration.append({'ranks':list(ranks),'matching_pairs_checked':len(patterns(ranks[0]))*len(patterns(ranks[1])),
            'consistent_pairs':count,'joint_weight_multisets':[[list(map(str,v)) for v in orbit] for orbit in sorted(orbits)],
            'all_solutions_unique_with_four_distinct_nonzero_weights':True})
        print('JOINT CHARACTER',ranks,'EXHAUSTIVE MATCHING PASS',count,flush=True)
    inputs=[path,Path(Ring.__module__.replace('.','/')+'.py')]
    import classification_quadratic_cofactor_review as helper
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,'full_joint_character_matching_classification':enumeration,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,Path(helper.__file__),Path(q.__file__)]},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'Independent all480 rows/all90 coefficients/rank/full Jacobi/Pfaffian/all-parameter SNC for the three sources, plus exhaustive symmetric-weight involution matchings including repeated weights and h=0; geometric reductions and the MPS theorem remain separately stated proof inputs.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_NONZERO_PAIR_COFACTOR_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))

def plane_cubic_support():
    import classification_two_lines_plane_cubic_component_review as q
    resource.setrlimit(resource.RLIMIT_CPU,(85,86));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1];q.N=6
    paths=[root/'reports'/name for name in ['CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json',
        'CLASSIFICATION_TWO_QUADRICS_PLANE_NORMALS_20261107.json','CLASSIFICATION_TWO_QUADRICS_PLANE_COMPONENT_REVIEW_20261107.json',
        'CLASSIFICATION_TWO_QUADRICS_PLANE_CUBIC_SUPPORT_20261107.json']]
    sources={v['name']:v for v in json.loads(paths[0].read_text())['cases']};normals={v['case']:v for v in json.loads(paths[1].read_text())['cases']}
    checked=json.loads(paths[2].read_text());assert checked['all_passed']
    assert checked['input_sha256'][str(paths[0].relative_to(root))]==hashlib.sha256(paths[0].read_bytes()).hexdigest()
    assert checked['input_sha256'][str(paths[1].relative_to(root))]==hashlib.sha256(paths[1].read_bytes()).hexdigest()
    checks={v['name']:v for v in checked['cases']};saved=json.loads(paths[3].read_text())['cases'];records=[]
    def plus(a,b):return {ij:q.add(a.get(ij,{}),b.get(ij,{})) for ij in set(a)|set(b)}
    J=lambda a:q.jvec(q.jacobi(a))
    def mixed(a,b):return [v-x-y for v,x,y in zip(J(plus(a,b)),J(a),J(b))]
    def vector_dict(v):return {str(i):str(a) for i,a in enumerate(v) if a}
    for packet in saved:
        begin=time.monotonic();name=packet['case'];point=q.decode(sources[name]['point_brackets']);R=[q.decode(v) for v in normals[name]['normals']];count=len(R)
        assert checks[name]['all_passed'] and checks[name]['complete_kernel_equals_source_plus_true_normals']
        assert checks[name]['true_normal_count']==count and checks[name]['full_actual_source_rank']==48
        Dcols=[]
        for ij,e in q.COLS:
            p={key:dict(value) for key,value in point.items()};p[ij]=q.add(p.get(ij,{}),{e+(0,):q.F(1)});p[ij[::-1]]=q.scale(p[ij],-1);Dcols.append(J(p))
        D=q.transpose(Dcols);rows,cols=packet['implicit_D_rows'],packet['implicit_D_columns']
        determinant=q.determinant([[D[i][j] for j in cols] for i in rows])
        assert determinant and determinant==q.F(packet['implicit_minor_determinant']) and len(cols)==102-count
        pairs=[tuple(v) for v in packet['normal_pairs']]
        assert pairs==([(0,0)] if count==1 else [(0,0),(0,1),(1,1)])
        Q=[J(R[i]) if i==j else mixed(R[i],R[j]) for i,j in pairs]
        T=[q.decode(v) for v in packet['quadratic_corrections']];vectors=[q.vec(v) for v in T]
        assert len(T)==len(Q)
        residual=[]
        for correction,quad in zip(vectors,Q):
            assert all(not a for i,a in enumerate(correction) if i not in cols)
            value=[a+sum(row[j]*correction[j] for j in cols) for a,row in zip(quad,D)]
            assert all(not value[i] for i in rows);residual.append(value)
        qrank=q.rank(q.transpose(residual))[0]
        assert qrank==checks[name]['quadratic_initial_space_rank']
        mixed_terms={(i,j):mixed(R[i],T[j]) for i in range(count) for j in range(len(T))}
        if count==1:C=[mixed_terms[0,0]]
        else:C=[mixed_terms[0,0],[a+b for a,b in zip(mixed_terms[1,0],mixed_terms[0,1])],
            [a+b for a,b in zip(mixed_terms[1,1],mixed_terms[0,2])],mixed_terms[1,2]]
        assert [vector_dict(v) for v in C]==packet['all_cubic_coefficient_vectors']
        pairings=[]
        for function in packet['cubic_functionals']:
            left={int(i):q.F(v) for i,v in function['nonzero_coefficients'].items()}
            assert all(not sum(a*D[i][j] for i,a in left.items()) for j in range(150))
            assert all(not sum(a*col[i] for i,a in left.items()) for col in Q)
            actual=[sum(a*col[i] for i,a in left.items()) for col in C]
            assert actual==list(map(q.F,function['cubic_pairing']));pairings.append(actual)
        crank=q.rank(pairings)[0] if pairings else 0
        isolated=False;reduced_certificate=None
        if name=='contact4_other_tangent':
            assert count==1 and qrank==0 and crank==1 and pairings[0][0]
            isolated=True;reduced_certificate={'nonzero_cubic_coefficient':str(pairings[0][0]),'reduced_normal_support':'origin'}
        elif name=='contact4_secant_at_contact':
            assert count==1 and qrank==0 and crank==0;reduced_certificate={'remaining_normal_support':'whole line at this order','stop_higher_jets':True}
        elif name=='contact31_tangent_at_simple_intersection':
            assert count==2 and qrank==2
            assert all(not sum(a*b for a,b in zip(row,[25**3,25**2,25,1])) for row in pairings)
            reduced_certificate={'quadratic_support_line':'t0=25*t1','all_cubic_pairings_vanish_on_that_line':True,'stop_higher_jets':True}
        elif name=='contact31_chord':
            assert count==2 and qrank==0 and crank==2
            gcd=()
            for row in pairings:gcd=r.ug(gcd,tuple(reversed(row)))
            assert gcd==(q.F(1),q.F(1)) and any(row[0] for row in pairings)
            reduced_certificate={'cubic_support_line':'t0=-t1','exact_monic_gcd_of_all_cubic_pairings':['1','1'],'no_infinity_point':True,'stop_higher_jets':True}
        else:raise AssertionError(name)
        assert isolated==packet['zero_dimensional_initial_support']
        records.append({'name':name,'all_passed':True,'uses_previously_independent_full_source_kernel_review':True,
            'all150_derivative_columns_and_selected_minor_reconstructed':True,'fixed_true_slice_quadratic_corrections_reconstructed':True,
            'all_cubic_coefficient_vectors_reconstructed':True,'all_cubic_functionals_annihilate_D_and_all_raw_normal_quadrics':True,
            'quadratic_initial_space_rank':qrank,'cubic_initial_subspace_rank':crank,'cubic_pairing_matrix':[[str(v) for v in row] for row in pairings],
            'reduced_support_certificate':reduced_certificate,'full23_component_certified':isolated,
            'component_implication':'The verified nonzero cubic in the genuine one-normal slice forces zero-dimensional transverse reduced support; Krull gives ambient lift dimension<=48, matching the full actual source. Its projective closure is a23-dimensional component.' if isolated else 'No component or membership conclusion; a structural full-source deformation is required after this single cubic test.',
            'seconds':time.monotonic()-begin})
        print(name,'INDEPENDENT CUBIC PASS','new23',isolated,flush=True)
    assert len(records)==4
    inputs=paths+[Path(q.__file__),Path(r.__file__)]
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'method':'Independent Fraction sparse arithmetic, using the previously independently verified full source/kernel only; no producer/CAS import. Every original derivative column, fixed-slice correction, cubic coefficient and actual cokernel pairing is reconstructed.',
        'scope':'One new23-dimensional component from cubic support; all three nonisolated outcomes are explicitly kept open with the exact remaining line/plane and no higher-jet iteration.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CUBIC_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))

def plane_components():
    import classification_two_lines_plane_cubic_component_review as q
    from itertools import product
    resource.setrlimit(resource.RLIMIT_CPU,(110,111));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    source_path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json'
    normal_path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_NORMALS_20261107.json'
    sources={v['name']:v for v in json.loads(source_path.read_text())['cases']};normals=json.loads(normal_path.read_text())['cases']
    records=[]
    for saved in normals:
        begin=time.monotonic();source=sources[saved['case']]
        rename=lambda text:text.replace('m0','r0').replace('m1','r1').replace('m2','r2').replace('m3','r3').replace('m4','r4').replace('m5','r5')
        q.N=11;universal=q.decode({ij:rename(v) for ij,v in source['full_six_magnetic_source'].items()})
        assert all(not v for v in q.rawjac(universal).values())
        values=[q.atom(1,j) for j in range(5)]+[q.atom(q.F(source['parameter_values']['m'+str(i)])) for i in range(6)]
        P={ij:q.substitute(poly,values) for ij,poly in universal.items()}
        directions=[{ij:q.substitute(q.der(poly,5+j),values) for ij,poly in universal.items()} for j in range(6)]
        reduce=lambda A:{ij:{e[:5]+(0,):v for e,v in poly.items()} for ij,poly in A.items() if poly}
        P=reduce(P);directions=list(map(reduce,directions));q.N=6
        assert P==q.decode(source['point_brackets']) and directions==[q.decode(v) for v in source['all6_parameter_directions']]
        R=[q.decode(v) for v in saved['normals']];count=len(R)
        Dcols=[]
        for ij,e in q.COLS:
            v={key:dict(value) for key,value in P.items()};v[ij]=q.add(v.get(ij,{}),{e+(0,):q.F(1)});v[ij[::-1]]=q.scale(v[ij],-1)
            Dcols.append(q.jvec(q.jacobi(v)))
        D=q.transpose(Dcols);VV=[];OO=[]
        for a,b in product(range(5),repeat=2):
            vertical={};orbit={}
            for i,j in q.PAIRS:
                vertical[i,j]=q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
                orbit[i,j]=q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
            VV.append(q.vec(vertical));OO.append(q.vec(orbit))
        Family=q.transpose(VV+OO+[q.vec(P)]+[q.vec(v) for v in directions])
        assert q.rank(q.transpose(VV))[0]==24 and q.rank(q.transpose(VV+OO+[q.vec(P)]))[0]==47
        assert q.rank(Family)[0]==48
        rv=[q.vec(v) for v in R];K=[row+[v[i] for v in rv] for i,row in enumerate(Family)]
        assert q.rank(K)[0]==48+count
        columns=q.transpose(K)
        for row in D:
            sparse=[(j,a) for j,a in enumerate(row) if a]
            assert all(not sum(a*col[j] for j,a in sparse) for col in columns)
        cert=source['certificates']['D'];rr,cc=cert['rows'],cert['columns']
        determinant=q.determinant([[D[i][j] for j in cc] for i in rr])
        assert determinant and determinant==q.F(cert['determinant']) and len(rr)+48+count==150
        squares=[q.jvec(q.jacobi(v)) for v in R]
        Qcols=[]
        for i in range(count):
            for j in range(i,count):
                if i==j:Qcols.append(squares[i])
                else:
                    total=q.jvec(q.jacobi({ij:q.add(R[i].get(ij,{}),R[j].get(ij,{})) for ij in set(R[i])|set(R[j])}))
                    Qcols.append([v-a-b for v,a,b in zip(total,squares[i],squares[j])])
        pairings=[]
        for function in saved['quadratic_cokernel_functionals']:
            left={int(i):q.F(v) for i,v in function['nonzero_coefficients'].items()}
            assert all(not sum(a*D[i][j] for i,a in left.items()) for j in range(150))
            actual=[sum(a*col[i] for i,a in left.items()) for col in Qcols]
            assert actual==list(map(q.F,function['quadratic_pairing']));pairings.append(actual)
        pr=q.rank(pairings)[0] if pairings else 0
        support0=(count==1 and pr==1) or (count==2 and pr==3)
        if count==2 and pr==2:
            _,basis=q.rank(q.transpose(pairings));a,b=[pairings[i] for i in basis]
            resultant=q.determinant([[a[0],a[1],a[2],0],[0,a[0],a[1],a[2]],[b[0],b[1],b[2],0],[0,b[0],b[1],b[2]]])
            support0=bool(resultant)
        assert support0==saved['zero_dimensional_quadratic_support']
        literal_axes=[i for i,v in enumerate(squares) if not any(v)]
        assert literal_axes==saved['literal_Poisson_normal_axes']
        extra_line=None
        if saved['case']=='contact31_tangent_at_simple_intersection':
            line={ij:q.add(q.scale(R[0].get(ij,{}),25),R[1].get(ij,{})) for ij in set(R[0])|set(R[1])}
            residual=q.jvec(q.jacobi(line));extra_line={'coefficients':[25,1],'literal_Poisson_pencil':not any(residual),'nonzero_residual_count':sum(bool(v) for v in residual)}
        records.append({'name':saved['case'],'all_passed':True,'all150_Jacobi_columns_reconstructed':True,
            'full_actual_source_rank':48,'full_Jacobi_rank':len(rr),'complete_kernel_equals_source_plus_true_normals':True,
            'true_normal_count':count,'all_cokernel_pairings_reconstructed':True,'quadratic_initial_space_rank':pr,
            'zero_projective_quadratic_support':support0,'full23_component_certified':support0,
            'literal_Poisson_normal_axes':literal_axes,'surviving_line_test':extra_line,
            'component_implication':'A genuine48-cut slice has an invertible Jacobi block and zero reduced normal support, so Krull bounds ambient dimension by48, matching actual source. Projectivization removes25 vertical/scale dimensions.' if support0 else 'Transverse support remains open; no component or local uniqueness is claimed.',
            'seconds':time.monotonic()-begin})
        print(saved['case'],'FULL FRACTION PASS','new23',support0,'axes',literal_axes,'line',extra_line,flush=True)
    inputs=[source_path,normal_path,Path(q.__file__)]
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'method':'Independent standard-library Fraction reconstruction from all original150 coefficients and all actual6 magnetic directions, complete true normals, exact derivative minors and every cokernel pairing; no producer/CAS import.',
        'scope':'Full ambient/kernel/obstruction arithmetic for all5 fixed plane sources. A23-dimensional component is certified only where the actual projective quadratic normal support is empty; distinctness and full-source membership use separate geometric proofs.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_COMPONENT_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))

def plane_sources():
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json'
    kernel_path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS_20261107.json'
    data=json.loads(path.read_text());kernel=json.loads(kernel_path.read_text())
    r.NAMES=('x','y','z','w','h')+tuple('m'+str(i) for i in range(6))
    x,y,z,w=(r.atom(1,i) for i in range(4));h=r.atom(1,4);base=[h,x,y]
    def combine(matrix,vector):return [r.add(*(r.scale(v,c) for c,v in zip(row,vector))) for row in matrix]
    def invert(a):
        determinant=r.det(a);assert determinant
        return [[((-1)**(i+j))*r.det([[a[k][l] for l in range(3) if l!=i] for k in range(3) if k!=j])/determinant for j in range(3)] for i in range(3)]
    records=[]
    for c in data['cases']:
        r.U=r.RF((r.F(0),r.F(1))) if c['residue']=='generic u' else r.RF(r.F(c['residue']))
        original=next(v for v in kernel['cases'] if (v['name'],v['residue'])==(c['name'],c['residue']))
        change=r.matrix(c['old_to_hxy_matrix']);inverse=invert(change);old_sub=combine(inverse,base)
        assert r.add(*(r.scale(r.atom(1,i),a) for i,a in enumerate(change[0])))==r.parse(original['hyperplane'])
        def substitute(poly):
            out={}
            for exp,coef in poly.items():
                assert not any(exp[3:]);term=r.atom(coef)
                for p,power in zip(old_sub,exp[:3]):
                    for _ in range(power):term=r.mul(term,p)
                out=r.add(out,term)
            return out
        anchors=[]
        for values in original['all_quadratic_gradient_syzygies']:
            transformed=combine(change,[substitute(r.parse(v)) for v in values]);quotient={}
            for exp,coef in transformed[0].items():
                assert exp[4]>0;ee=list(exp);ee[4]-=1;quotient[tuple(ee)]=coef
            anchor=[r.add(a,r.scale(r.mul(quotient,b),-1)) for a,b in zip(transformed,base)]
            assert not anchor[0];anchors.append(anchor)
        U,V=anchors
        assert U==list(map(r.parse,c['all_horizontal_anchor_U'])) and V==list(map(r.parse,c['all_horizontal_anchor_V']))
        q1,q2=[substitute(r.parse(original[key])) for key in ['q1','q2']]
        assert q1==r.parse(c['q1']) and q2==r.parse(c['q2'])
        aa,bb=map(r.parse,c['bracket_coefficients'])
        bracket=[r.add(*(r.add(r.mul(U[i],r.der(V[j],k)),r.scale(r.mul(V[i],r.der(U[j],k)),-1)) for i,k in enumerate([4,0,1]))) for j in range(3)]
        assert bracket==[r.add(r.mul(aa,a),r.mul(bb,b)) for a,b in zip(U,V)]
        magnetic=r.add(*(r.mul(r.atom(1,5+i),mon) for i,mon in enumerate([r.mul(h,h),r.mul(h,x),r.mul(h,y),r.mul(x,x),r.mul(x,y),r.mul(y,y)])))
        P=[{},U[1],V[1],U[2],V[2],r.add(r.scale(r.mul(aa,z),-1),r.scale(r.mul(bb,w),-1),magnetic)]
        assert P==list(map(r.parse,c['brackets']))
        matrix={}
        for (i,j),p in zip(r.PAIRS,P):matrix[i,j]=p;matrix[j,i]=r.scale(p,-1)
        for i,j,k in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]:
            jac=r.add(*(r.mul(matrix.get((a,t),{}),r.der(matrix.get((b,d),{}),t)) for a,b,d in [(i,j,k),(j,k,i),(k,i,j)] for t in range(4)))
            assert not jac
        pf=r.pf(P);ratio=r.scalar(c['Pfaffian_multiple_of_q1q2']);assert ratio
        assert pf==r.parse(c['full_Pfaffian'])==r.scale(r.mul(q1,q2),ratio)
        records.append({'name':c['name'],'residue':c['residue'],'full_kernel_anchors_independently_transformed':True,
            'all_six_magnetic_coefficients_and_full_source_reconstructed':True,'all_four_Jacobi_and_actual_factorization':True,'all_passed':True})
        print(c['name'],c['residue'],'ALL SIX MAGNETIC PASS',flush=True)
    assert len(records)==8
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'method':'Standard-library Fraction(u) sparse reconstruction from the complete gradient kernels, independent coordinate-matrix inversion, horizontalization, all6 magnetic coefficients, all4 Jacobi and actual Pfaffian.',
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,kernel_path,Path(r.__file__)]},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start,
        'scope':'All eight full magnetic sources (seven configurations plus one valid parameter exception). Geometric completeness and component recognition remain separate proofs.'}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))


def main():
    if '--two-lines-rank2-translation' in sys.argv:return two_lines_rank2_translation()
    if '--reduced-divisor-normalizers' in sys.argv:return reduced_divisor_normalizers()
    if '--nilpotent-factor-reductions' in sys.argv:return nilpotent_factor_reductions()
    if '--mixed-normal-support' in sys.argv:return mixed_normal_support()
    if '--nonzero-pair-cofactor' in sys.argv:return nonzero_pair_cofactor()
    if '--plane-cubic-support' in sys.argv:return plane_cubic_support()
    if '--plane-components' in sys.argv:return plane_components()
    if '--plane-sources' in sys.argv:return plane_sources()
    resource.setrlimit(resource.RLIMIT_CPU,(55,56));resource.setrlimit(resource.RLIMIT_AS,(768*1024**2,768*1024**2))
    start=time.monotonic();root=Path(__file__).resolve().parents[1]
    path=root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS_20261107.json';data=json.loads(path.read_text())
    expected={
        'shared_tangent_generic':('x','u*x*z-y*y+x*x'),
        'shared_tangent_contact3':('x','x*z-y*y+x*y'),
        'shared_tangent_contact4':('x','x*z-y*y+x*x'),
        'contact31_common_bitangent':('x+8*y+16*z','x*z-y*y+x*y'),
        'contact4_other_tangent':('z','x*z-y*y+x*x'),
        'contact31_tangent_at_simple_intersection':('z','x*z-y*y+x*y'),
        'contact4_secant_at_contact':('y','x*z-y*y+x*x'),
        'contact22_chord':('y','x*z-u*y*y'),
        'contact31_chord':('y','x*z-y*y+x*y')}
    assert {c['name'] for c in data['cases']}==set(expected)
    extra={(c['name'],c['residue']) for c in data['cases'] if c['residue']!='generic u'}
    assert extra=={('shared_tangent_generic','1/2'),('contact22_chord','-1')}
    records=[]
    for c in data['cases']:
        begin=time.monotonic();r.U=r.RF((r.F(0),r.F(1))) if c['residue']=='generic u' else r.RF(r.F(c['residue']))
        x,y,z=(r.atom(1,i) for i in range(3));xx=[x,y,z]
        H,q2=map(r.parse,expected[c['name']]);q1=r.parse('x*z-y*y');D=r.mul(H,r.mul(q1,q2))
        assert H==r.parse(c['hyperplane']) and q1==r.parse(c['q1']) and q2==r.parse(c['q2']) and D==r.parse(c['plane_quintic'])
        mon2=[r.mul(xx[i],xx[j]) for i,j in r.combinations_with_replacement(range(3),2)]
        assert mon2==list(map(r.parse,c['quadratic_monomials']))
        columns=[r.mul(mon,r.der(D,i)) for i in range(3) for mon in mon2]
        mon6=sorted(set().union(*(p.keys() for p in columns)),reverse=True)
        M=r.transpose([[p.get(e,r.RF()) for e in mon6] for p in columns])
        assert M==r.matrix(c['gradient_syzygy_matrix'])
        K=r.matrix(c['complete_kernel']);dim=c['kernel_dimension'];rank=c['rank']
        assert r.annihilate(M,K) and r.rank(K)==dim and dim+rank==18
        pivot=r.det(r.submatrix(M,c['pivot_rows'],c['pivot_columns']))
        assert pivot and pivot==r.scalar(c['pivot_determinant']) and len(c['pivot_rows'])==rank
        fields=[]
        for j in range(dim):
            fields.append([r.add(*(r.scale(mon,K[6*i+k][j]) for k,mon in enumerate(mon2))) for i in range(3)])
        assert fields==[list(map(r.parse,v)) for v in c['all_quadratic_gradient_syzygies']]
        ratios=[]
        for item in c['Saito_determinants']:
            i,j=item['kernel_columns'];a,b=fields[i],fields[j]
            determinant=r.add(r.mul(x,r.add(r.mul(a[1],b[2]),r.scale(r.mul(a[2],b[1]),-1))),
                r.scale(r.mul(y,r.add(r.mul(a[0],b[2]),r.scale(r.mul(a[2],b[0]),-1))),-1),
                r.mul(z,r.add(r.mul(a[0],b[1]),r.scale(r.mul(a[1],b[0]),-1))))
            ratio=r.scalar(item['Saito_determinant_divided_by_D']);ratios.append(ratio)
            assert determinant==r.scale(D,ratio)
        assert len(ratios)==dim*(dim-1)//2 and bool(any(ratios))==c['has_free_degree22_basis']
        if c['residue']=='generic u':
            roots={r.F(d['residue']) for d in data['cases'] if d['name']==c['name'] and d['residue']!='generic u'}
            allowed=roots|{r.F(0),r.F(1)};found=set()
            for v in [pivot]+[v for v in ratios if v]:found|=r.roots_supported(v.p,allowed)|r.roots_supported(v.q,allowed)
            for row in K:
                for v in row:found|=r.roots_supported(v.q,allowed)
            assert found-{r.F(0),r.F(1)}==roots
        records.append({'name':c['name'],'residue':c['residue'],'complete_kernel_dimension':dim,
            'all18_coefficients_and_full_matrix_reconstructed':True,'all_Saito_determinants_reconstructed':True,
            'free_degree22_basis_exists':bool(any(ratios)),'all_passed':True,'seconds':time.monotonic()-begin})
        print(c['name'],c['residue'],'PASS',dim,bool(any(ratios)),flush=True)
    assert len(records)==11
    r.U=r.RF((r.F(0),r.F(1)))
    excluded={
        'shared_tangent_contact22':('x','x*z-u*y*y',['4*x','-y','-6*z']),
        'shared_tangent_contact4':('x','x*z-y*y+x*x',['0','x','2*y']),
        'contact22_chord':('y','x*z-u*y*y',['x','0','-z'])}
    assert {c['name'] for c in data['analytic_exclusions']}==set(excluded)
    for c in data['analytic_exclusions']:
        H,q2,v=excluded[c['name']];H,q2=r.parse(H),r.parse(q2);v=list(map(r.parse,v))
        D=r.mul(H,r.mul(r.parse('x*z-y*y'),q2))
        assert v==list(map(r.parse,c['nonzero_linear_gradient_syzygy'])) and any(v)
        assert not r.add(*(r.mul(a,r.der(D,i)) for i,a in enumerate(v)))
    out={'task_id':'E_reduced_nonnormal','all_passed':True,'cases':records,
        'method':'Standard-library Fraction(u) sparse polynomial matrix/kernel/minor and full Saito determinant reconstruction, with all exceptional parameter values; no producer/CAS imported.',
        'all_three_global_linear_gradient_obstructions_reconstructed':True,
        'scope':'Nine complete contact-budget kernels, two non-excluded pivot exceptions and three global linear-gradient exclusions covering the tenth candidate. Geometric completeness, full magnetic-source reduction and component membership remain separate proofs.',
        'input_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [path,Path(r.__file__)]},
        'reviewer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    (root/'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS_REVIEW_20261107.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'seconds':out['seconds']}))

if __name__=='__main__':main()
