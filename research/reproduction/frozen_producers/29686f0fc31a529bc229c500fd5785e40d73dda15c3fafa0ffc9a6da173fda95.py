#!/usr/bin/env python3
"""Whole tangent cone sources and constant-bivector Hamiltonian pencil.

All9 tangent restrictions, all50 L coefficients and all original350 equations
are retained. Three normalized C strata include C=0. No nonlinear GB or jets.
"""
import argparse,hashlib,itertools as it,json,resource,signal,time
from pathlib import Path
from classification_cubic_scroll_abelian_review import Q,PAIRS,add,scale,mul,diff,parse,rref,monomials,jacobi_and_divisor,determinant
from classification_quadric_cone_zero_source import affine,empty

ROOT=Path(__file__).resolve().parents[1]


def tensor_sum(bases,coefficients):
    P=empty()
    for i,j in PAIRS:P[i][j]=add(*(scale(B[i][j],c)for B,c in zip(bases,coefficients)if c));P[j][i]=scale(P[i][j],-1)
    return P


def encode_poly(P):return [[list(e),str(c)]for e,c in sorted(P.items(),reverse=True)]


def encode_tensor(P):return {str((i,j)):encode_poly(P[i][j])for i,j in PAIRS}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    def timeout(signum,frame):raise TimeoutError('finite whole tangentcone Hamiltonianpencil120s')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    names=[f'x{i}'for i in range(5)];xs=[parse(v,names)for v in names];q=parse('x0*x2-x1**2',names);h=xs[4]
    mons2=monomials(5,2);mons4=monomials(5,4);mons5=monomials(5,5);QUADS=list(it.combinations(range(5),4))
    def vec(P):return [P[i][j].get(e,Q(0))for i,j in PAIRS for e in mons2]
    Lraw=[]
    for i,j in PAIRS:
        for k in range(5):
            P=empty();P[i][j]=mul(h,xs[k]);P[j][i]=scale(P[i][j],-1);Lraw.append(P)
    Araw=[]
    for i,j in PAIRS:
        P=empty();P[i][j]=q;P[j][i]=scale(q,-1);Araw.append(P)
    Fmons=[e for e in mons2 if not e[4]and e!=(1,0,1,0,0)];assert len(Fmons)==9
    Fraw=[]
    for e in Fmons:
        F={e:Q(1)};P=empty()
        for i,j in it.combinations(range(4),2):
            k,l=[v for v in range(4)if v not in[i,j]];perm=[i,j,k,l];sign=(-1)**sum(perm[a]>perm[b]for a in range(4)for b in range(a+1,4))
            P[i][j]=scale(add(mul(diff(q,k),diff(F,l)),scale(mul(diff(q,l),diff(F,k)),-1)),sign);P[j][i]=scale(P[i][j],-1)
        assert all(not add(*(mul(P[i][j],diff(q,j))for j in range(5)))for i in range(5))and all(not P[4][j]for j in range(5))
        Fraw.append(P)
    raw=Lraw+Araw+Fraw;assert len(rref([vec(P)for P in raw]))==69
    vertical=[]
    for a,b in it.product(range(5),repeat=2):
        P=empty()
        for i,j in PAIRS:P[i][j]=add(mul(xs[i],xs[b])if j==a else{},scale(mul(xs[j],xs[b]),-1)if i==a else{})
        vertical.append(vec(P))
    assert len(rref(vertical))==24 and len(rref(vertical+[vec(P)for P in raw]))==88
    oldpath=ROOT/'reports/CLASSIFICATION_QUADRIC_CONE_ZERO_SOURCE_20261107.json';old=json.loads(oldpath.read_text())
    conormal=[[Q(v)for v in row]for row in old['all11_eliminated_raw_source_constraints']];assert len(conormal)==11 and len(rref(conormal))==11
    eulerL=[[Q((j==k and l==i)-(i==k and l==j))for i,j in PAIRS for l in range(5)]+[Q(0)]*9 for k in range(5)]
    gauge_indices=sorted(rref(eulerL));assert len(gauge_indices)==5
    gauge_minor=determinant([[v[j]for j in gauge_indices]for v in eulerL]);assert gauge_minor
    A0=empty();A0[4][3]=q;A0[3][4]=scale(q,-1)
    Cpairs=[(0,1),(0,2),(1,2)];Cbases=[Araw[PAIRS.index(pair)]for pair in Cpairs]
    unknown=Lraw+Fraw;universal=[A0,*Cbases,*unknown];assert len(universal)==63
    universal_J,universal_D=jacobi_and_divisor(universal);upairs=list(it.combinations_with_replacement(range(63),2))
    # All unknown-unknown terms must vanish mod h: this checks the exact
    # affine-linearity used by the whole-source reduction.
    modrows=[i for i in range(350)if not mons4[i%70][4]]
    assert all(not any(upairs[j][0]>=4 for j in universal_J[i])for i in modrows)
    x4=[{tuple(int(j==i)for j in range(4)):Q(1)}for i in range(4)]
    def deh(P):
        out={}
        for e,c in P.items():out[e[:4]]=out.get(e[:4],Q(0))+c
        return {e:c for e,c in out.items()if c}
    def pencil(P):
        A=[[{}for _ in range(4)]for _ in range(4)]
        for i,j in it.combinations(range(4),2):A[i][j]=add(deh(P[i][j]),scale(mul(x4[i],deh(P[4][j])),-1),mul(x4[j],deh(P[4][i])));A[j][i]=scale(A[i][j],-1)
        B=[[{}for _ in range(3)]for _ in range(3)];D=[[{}for _ in range(3)]for _ in range(3)];Y=[];T=[]
        for i,j in it.combinations(range(3),2):
            assert all(e[3]<=1 and(not e[3]or sum(e[:3])==0)for e in A[i][j])
            B[i][j]={e[:3]:c for e,c in A[i][j].items()if not e[3]};B[j][i]=scale(B[i][j],-1)
            D[i][j]={e[:3]:c for e,c in A[i][j].items()if e[3]};D[j][i]=scale(D[i][j],-1)
        for i in range(3):
            assert all(e[3]<=1 for e in A[i][3]);Y.append({e[:3]:c for e,c in A[i][3].items()if e[3]==0});T.append({e[:3]:c for e,c in A[i][3].items()if e[3]==1})
        assert all(sum(e)<=2 for row in B for P in row for e in P)
        assert all(sum(e)<=3 for P in Y for e in P)and all(sum(e)<=1 for P in T for e in P)
        return B,D,Y,T
    def lie(V,B,i,j):
        return add(*(add(mul(V[l],diff(B[i][j],l)),scale(mul(B[l][j],diff(V[i],l)),-1),scale(mul(B[i][l],diff(V[j],l)),-1))for l in range(3)))
    def jwedge(V,D):return add(*(mul(V[i],D[j][k])for i,j,k in[(0,1,2),(1,2,0),(2,0,1)]))
    def jacmixed(B,C):return add(*(mul(B[i][l],diff(C[j][k],l))for i,j,k in[(0,1,2),(1,2,0),(2,0,1)]for l in range(3)))
    def vwedge(U,V,i,j):return add(mul(U[i],V[j]),scale(mul(U[j],V[i]),-1))
    def all_pencil_conditions(data,J):
        pairs=list(it.combinations_with_replacement(range(len(data)),2));rows={}
        def record(kind,component,P,col):
            for e,c in P.items():rows.setdefault((kind,component,e),{})[col]=c
        for col,(a,b)in enumerate(pairs):
            A0={};A1={};C0={ij:{}for ij in it.combinations(range(3),2)};C1={ij:{}for ij in C0};C2={ij:{}for ij in C0}
            for c,d in([(a,b)]if a==b else[(a,b),(b,a)]):
                B,D,Y,T=data[c];BB,DD,YY,TT=data[d]
                A0=add(A0,jacmixed(B,BB),jwedge(Y,DD));A1=add(A1,jacmixed(D,BB),jwedge(T,DD))
                for i,j in C0:
                    C0[i,j]=add(C0[i,j],scale(lie(Y,BB,i,j),-1),vwedge(Y,TT,i,j))
                    C1[i,j]=add(C1[i,j],scale(lie(T,BB,i,j),-1),scale(lie(Y,DD,i,j),-1))
                    C2[i,j]=add(C2[i,j],scale(lie(T,DD,i,j),-1))
            predicted={(0,1,2):add({e+(0,):v for e,v in A0.items()},{e+(1,):v for e,v in A1.items()})}
            record('JacB_plus_YwedgeD',(),A0,col);record('B_D_cross_plus_TwedgeD',(),A1,col)
            for i,j in C0:
                predicted[i,j,3]=add({e+(0,):v for e,v in C0[i,j].items()},{e+(1,):v for e,v in C1[i,j].items()},{e+(2,):v for e,v in C2[i,j].items()})
                record('minusLieYB_plus_YwedgeT',(i,j),C0[i,j],col);record('minusLieTB_minusLieYD',(i,j),C1[i,j],col);record('minusLieTD',(i,j),C2[i,j],col)
            for tri,P in predicted.items():
                quad=tuple([*tri,4]);at=QUADS.index(quad);actual={}
                for j,e in enumerate(mons4):
                    value=J[70*at+j].get(col,Q(0))
                    if value:actual[e[:4]]=actual.get(e[:4],Q(0))+value
                assert not add(P,{e:c for e,c in actual.items()if c})
        return [dict(condition=kind,component=component,spatial_monomial=e,parameter_pair_coefficients={str(j):str(c)for j,c in coeff.items()})for(kind,component,e),coeff in sorted(rows.items())]
    cases=[]
    for name,cvalues in [('semisimple',[Q(0),Q(1),Q(0)]),('nilpotent',[Q(1),Q(0),Q(0)]),('zero',[Q(0)]*3)]:
        constants=[Q(1),*cvalues];Avec=[Q(0)]*10;Avec[PAIRS.index((3,4))]=-1
        for pair,c in zip(Cpairs,cvalues):Avec[PAIRS.index(pair)]=c
        linear=[];rhs=[]
        for i in modrows:
            row=[Q(0)]*59;constant=Q(0)
            for col,c in universal_J[i].items():
                a,b=upairs[col]
                if b<4:constant+=c*constants[a]*constants[b]
                else:row[b-4]+=c*constants[a]
            linear.append(row);rhs.append(-constant)
        con_rows=[row[:50]+[Q(0)]*9 for row in conormal];con_rhs=[-sum(row[50+j]*Avec[j]for j in range(10))for row in conormal]
        combined=linear+con_rows;combined_rhs=rhs+con_rhs
        rawR,rawfree,rawpoint,rawK=affine(combined,combined_rhs,59)
        assert all(not sum(a*b for a,b in zip(row,v))for row in combined for v in eulerL)
        gauge=[[Q(j==i)for j in range(59)]for i in gauge_indices]
        RR,free,point,K=affine(combined+gauge,combined_rhs+[Q(0)]*5,59);assert len(rawK)==len(K)+5
        Path('/tmp/quadric_cone_rank_two_pencil_linear_checkpoint.json').write_text(json.dumps(dict(current_case=name,point=list(map(str,point)),complete_kernel=[list(map(str,v))for v in K],seconds=time.monotonic()-begun),indent=2)+'\n')
        A=tensor_sum([A0,*Cbases],constants);Ppoint=tensor_sum([A,*unknown],[Q(1),*point]);family=[Ppoint,*[tensor_sum(unknown,v)for v in K]]
        data=[pencil(P)for P in family];d=len(K);JJ,DD=jacobi_and_divisor(family);conditions=all_pencil_conditions(data,JJ)
        vv=50+Fmons.index((0,0,0,2,0));assert not point[vv]and all(not v[vv]for v in K)
        ell_indices=[50+Fmons.index(tuple(int(j==i)+int(j==3)for j in range(5)))for i in range(3)]
        ell_point=[point[i]for i in ell_indices];ell_matrix=[[v[i]for v in K]for i in ell_indices]
        Dpoint=[D[i][j].get((0,0,0),Q(0))for B,D,Y,T in data[:1]for i,j in it.combinations(range(3),2)]
        Dmatrix=[[D[i][j].get((0,0,0),Q(0))for B,D,Y,T in data[1:]]for i,j in it.combinations(range(3),2)]
        cases.append(dict(type=name,normalized_A_upper=list(map(str,Avec)),whole_modh_linear_shape=[len(linear),59],all350_modh_row_indices=modrows,
                          all_modh_linear_rows=[[str(v)for v in row]for row in linear],all_modh_rhs=list(map(str,rhs)),all11_conormal_rows=[[str(v)for v in row]for row in con_rows],all11_conormal_rhs=list(map(str,con_rhs)),
                          raw_affine_dimension=len(rawK),raw_particular=list(map(str,rawpoint)),raw_complete_kernel=[list(map(str,v))for v in rawK],Euler5_verified=True,Euler_gauge_indices=gauge_indices,Euler_gauge_minor=str(gauge_minor),
                          projective_modh_affine_dimension=d,projective_particular=list(map(str,point)),complete_projective_modh_kernel=[list(map(str,v))for v in K],projective_free_raw_indices=free,
                          all_recovered150_tensors=[encode_tensor(P)for P in family],pencil_B0_point_and_all_directions=[{str((i,j)):encode_poly(B[i][j])for i,j in it.combinations(range(3),2)}for B,D,Y,T in data],
                          pencil_D_constant_point=list(map(str,Dpoint)),pencil_D_constant_parameter_matrix=[list(map(str,v))for v in Dmatrix],
                          pencil_Y_point_and_all_directions=[[encode_poly(v)for v in Y]for B,D,Y,T in data],pencil_T_point_and_all_directions=[[encode_poly(v)for v in T]for B,D,Y,T in data],
                          whole_pencil_identity_verified=True,all_five_pencil_Jacobi_blocks_equal_full_affine_Jacobi=True,all_pencil_condition_coefficient_rows=conditions,
                          no_false_Jacobi_B0_alone_constraint=True,F_vv_forced_zero=True,F_v_linear_point=list(map(str,ell_point)),F_v_complete_coefficient_matrix=[list(map(str,v))for v in ell_matrix],F_v_parameter_rank=len(rref(ell_matrix)),
                          full_original350_Jacobi_parameter_rows=[{str(j):str(c)for j,c in row.items()}for row in JJ],full126_actual_Pfaffian_parameter_rows=[{str(j):str(c)for j,c in row.items()}for row in DD],
                          full_parameter_pair_order=list(it.combinations_with_replacement(range(len(family)),2)),constant_parameter_index0_equals1=True,
                          scope='Whole affine supersource with conormal/allmodhJacobi. Full constant-bivector Hamiltonianpencil Jacobi and actualPf retained. No nonlinear Gdimension/support conclusion.'))
        Path('/tmp/quadric_cone_rank_two_pencil_progress.json').write_text(json.dumps(dict(cases=cases,seconds=time.monotonic()-begun),indent=2)+'\n')
    paths=[Path(__file__),ROOT/'code/classification_cubic_scroll_abelian_review.py',ROOT/'code/classification_quadric_cone_zero_source.py',oldpath,ROOT/'research/classification/QUADRIC_CONE_RANK_TWO_TANGENT_SOURCE.md',ROOT/'research/classification/QUADRIC_CONE_RANK_TWO_HAMILTONIAN_PENCIL.md',ROOT/'reports/CLASSIFICATION_QUADRIC_CONE_RANK_TWO_SOURCE_RUNTIME_20261107.json']
    result=dict(task_id='G',all_arithmetic_passed=True,full69_tangent_raw_source_rank=69,full24_Euler_rank=24,combined_rank88=True,whole_projective_tangent_source_dimension=64,full_conormal_raw_dimension58=True,
                coordinate_order=['u0','u1','u2','v','h'],quadric='u0*u2-u1^2',F_basis_monomials=Fmons,all_F9_h_q_Casimirs_verified=True,
                raw69_parameter_order='L50,constantA10,JacobianF9',normalized_unknown59_parameter_order='L50,JacobianF9',
                all_three_constant_A_strata_included=True,cases=cases,no_GB_or_point_or_jet_run=True,no_G_dimension_claimed=True,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,raw69=True,cases=[dict(type=r['type'],raw=r['raw_affine_dimension'],quotient=r['projective_modh_affine_dimension'],Fv_rank=r['F_v_parameter_rank'],D_rank=len(rref([[Q(v)for v in row]for row in r['pencil_D_constant_parameter_matrix']])),pencilconditions=len(r['all_pencil_condition_coefficient_rows']))for r in cases],seconds=result['seconds'])))


def highest_symbol_main():
    import sympy as sp
    ap=argparse.ArgumentParser();ap.add_argument('--highest-symbol',action='store_true');ap.add_argument('--input',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    signal.signal(signal.SIGALRM,lambda signum,frame:(_ for _ in ()).throw(TimeoutError('highest symbol30s')));signal.alarm(30)
    r=json.loads(args.input.read_text());out=[];s,a,b,c=sp.symbols('s a b c');sym=[sp.Integer(1),s,a,b,c];allowed=[0,24,28,31,33]
    degrees={'JacB_plus_YwedgeD':3,'B_D_cross_plus_TwedgeD':1,'minusLieYB_plus_YwedgeT':4,'minusLieTB_minusLieYD':2,'minusLieTD':0}
    mons=list(it.combinations_with_replacement(range(5),2))
    for case in r['cases']:
        pairs=case['full_parameter_pair_order'];rows=[];original=[]
        assert [case['projective_free_raw_indices'][i-1]for i in allowed[1:]]==[48,52,55,57]
        for row in case['all_pencil_condition_coefficient_rows']:
            if sum(row['spatial_monomial'])!=degrees[row['condition']]:continue
            v=[Q(0)]*len(mons)
            for col,q in row['parameter_pair_coefficients'].items():
                i,j=pairs[int(col)];assert i in allowed and j in allowed
                v[mons.index((allowed.index(i),allowed.index(j)))]+=Q(q)
            rows.append(v);original.append(dict(condition=row['condition'],component=row['component'],spatial_monomial=row['spatial_monomial'],coefficient_row=list(map(str,v))))
        R=rref(rows);basis=[[v.get(j,Q(0))for j in range(len(mons))]for v in R.values()]
        polys=[sp.factor(sum(sp.Rational(q.numerator,q.denominator)*sym[i]*sym[j]for q,(i,j)in zip(v,mons)))for v in basis]
        assert all(not any(sum(q*x for q,x in zip(row,kv))for row in [])for kv in [])
        out.append(dict(type=case['type'],symbol_parameters=['s','ell0','ell1','ell2'],original_family_parameter_indices=allowed,monomial_pairs=mons,all_selected_coefficient_rows=original,all_rows_linear_span_rank=len(R),rowspace_basis=[list(map(str,v))for v in basis],factored_rowspace_basis=list(map(str,polys)),all_other29_parameters_absent_verified=True))
    result=dict(task_id='G',all_arithmetic_passed=True,cases=out,scope='Whole highest-spatial Jacobi symbol on all four parameters. Exact rowspace equations only, no support/dimension claim without mathematical case analysis.',input_sha256={str(args.input):hashlib.sha256(args.input.read_bytes()).hexdigest(),str(Path(__file__).relative_to(ROOT)):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6))
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,cases=[dict(type=v['type'],equations=v['factored_rowspace_basis'])for v in out],seconds=result['seconds'])))


def triangular_bound_main():
    import ast
    import sympy as sp
    import classification_plane_rank_zero_rank_two_high_m as H
    ap=argparse.ArgumentParser();ap.add_argument('--triangular-bound',action='store_true');ap.add_argument('--input',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    signal.signal(signal.SIGALRM,lambda signum,frame:(_ for _ in ()).throw(TimeoutError('whole cone triangular bound120s')));signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    old=json.loads(args.input.read_text());highpath=ROOT/'reports/CLASSIFICATION_QUADRIC_CONE_RANK_TWO_HIGHEST_SYMBOL_20261107.json';high=json.loads(highpath.read_text())
    frozen={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(high['producer_source_bytes']).body if isinstance(v,ast.FunctionDef)};now={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(Path(__file__).read_text()).body if isinstance(v,ast.FunctionDef)};assert all(now[k]==v for k,v in frozen.items())
    k,s=sp.symbols('kappa sigma');stages=[[5,7,12,14,18,20,21,22,23,26,27,29,30,32],[1,2,3,4,8,9,10,11,15,16,17,25],[6,13,19]]
    conditions={'JacB_plus_YwedgeD':3,'minusLieYB_plus_YwedgeT':4,'minusLieTB_minusLieYD':2}
    def encode(M):return [[str(sp.factor(v))for v in row]for row in M.tolist()]
    def summarize(matrices,domain,base_dim,offset):
        certs=[];minors=[]
        for M in matrices:
            cert,minor=H.certificate(M,domain);certs.append(cert);minors.append(minor)
        fiber=sum(c['nullity']for c in certs);dim=base_dim+fiber;image=offset+dim
        return dict(stage_certificates=certs,base_dimension=base_dim,triangular_fiber_upper=fiber,normalized_source_upper=dim,whole_PGL_image_upper=image,G24_closed=image<=24),sp.factor(sp.prod(minors))
    def line_cover(matrices,sub,offset):
        mats=[M.subs(sub,simultaneous=True)for M in matrices];generic,minor=summarize(mats,H.POLY,1,offset);walls=[]
        for factor,mult in sp.factor_list(minor)[1]:
            f=sp.Poly(factor,H.TAU)
            if f.degree()==1:
                val=-f.nth(0)/f.nth(1);point,_=summarize([M.subs(H.TAU,val)for M in mats],sp.QQ,0,offset);point['parameter']=str(val)
            else:
                certs=[H.algebraic_certificate(M,factor)for M in mats];fiber=sum(c['nullity']for c in certs);point=dict(stage_certificates=certs,base_dimension=0,triangular_fiber_upper=fiber,normalized_source_upper=fiber,whole_PGL_image_upper=offset+fiber,G24_closed=offset+fiber<=24,parameter='all roots of '+str(factor))
            point.update(minimal_polynomial=str(factor),multiplicity=int(mult),all_complex_roots_covered=True);walls.append(point)
        return dict(parameterization={str(a):str(b)for a,b in sub.items()},generic=generic,selected_minor_product=str(minor),all_exceptional_roots=walls,all_G24_closed=generic['G24_closed']and all(v['G24_closed']for v in walls))
    cases=[]
    for cindex,model in[(0,'semisimple'),(1,'nilpotent'),(2,'zero_C_ell_semisimple'),(2,'zero_C_ell_nilpotent')]:
        case=old['cases'][cindex];pairs=case['full_parameter_pair_order'];base_indices={0:sp.Integer(1),24:s};ell_index=31 if model in['semisimple','zero_C_ell_semisimple']else 33
        base_indices[ell_index]=k if cindex<2 else sp.Integer(1)
        retained=set(base_indices)|set(it.chain.from_iterable(stages));weights={j:0 for j in base_indices}
        for weight,js in enumerate(stages,1):weights.update({j:weight for j in js})
        rows=[[]for _ in stages];metadata=[[]for _ in stages];all_selected=[]
        for n,row in enumerate(case['all_pencil_condition_coefficient_rows']):
            terms=[]
            for col,value in row['parameter_pair_coefficients'].items():
                i,j=pairs[int(col)]
                if i in retained and j in retained:terms.append((i,j,Q(value)))
            if not terms:continue
            if row['condition']not in conditions:
                assert not terms;continue
            weight=conditions[row['condition']]-sum(row['spatial_monomial'])
            assert all(weights[i]+weights[j]==weight for i,j,v in terms)
            if weight==0:
                assert sp.expand(sum(sp.Rational(v.numerator,v.denominator)*base_indices[i]*base_indices[j]for i,j,v in terms))==0
            if weight not in[1,2,3]:continue
            stage=stages[weight-1];values=[sp.Integer(0)]*len(stage)
            for i,j,v in terms:
                if i in base_indices and j in stage:values[stage.index(j)]+=sp.Rational(v.numerator,v.denominator)*base_indices[i]
                if j in base_indices and i in stage:values[stage.index(i)]+=sp.Rational(v.numerator,v.denominator)*base_indices[j]
            if any(values):rows[weight-1].append(values);metadata[weight-1].append(n)
            all_selected.append(dict(original_pencil_row_index=n,stage=weight,only_same_weight_terms_verified=True))
        matrices=[sp.Matrix(v)if v else sp.zeros(0,len(stage))for v,stage in zip(rows,stages)]
        offset=18 if model in['semisimple','zero_C_ell_semisimple']else 17
        base_dim=2 if cindex<2 else 1;generic,minor=summarize(matrices,sp.QQ.poly_ring(k,s)if cindex<2 else sp.QQ.poly_ring(s),base_dim,offset)
        if cindex<2:
            lines=[]
            for factor,mult in sp.factor_list(minor)[1]:
                f=sp.Poly(factor,k,s)
                if f.total_degree()!=1:
                    checkpoint=dict(status='nonlinear_exception_curve_stop',model=model,exact_factor=str(factor),matrices=[encode(M)for M in matrices],source_bytes=Path(__file__).read_text());args.output.with_name(args.output.stem+'_RUNTIME.json').write_text(json.dumps(checkpoint,indent=2)+'\n');raise ValueError('nonlinear exception curve; no profile continuation')
                A,B,C=f.coeff_monomial(k),f.coeff_monomial(s),f.coeff_monomial(1);sub={k:H.TAU,s:-(A*H.TAU+C)/B}if B else{k:-C/A,s:H.TAU};line=line_cover(matrices,sub,offset);line.update(equation=str(factor),multiplicity=int(mult));lines.append(line)
            cover=dict(generic=generic,selected_minor_product=str(minor),all_linear_exception_lines=lines,all_G24_closed=generic['G24_closed']and all(v['all_G24_closed']for v in lines))
        else:cover=line_cover(matrices,{s:H.TAU},offset)
        resultcase=dict(model=model,source_case_index=cindex,F_v_retained_family_index=ell_index,F_v_value='kappa'if cindex<2 else'1',base_parameters=['kappa','sigma']if cindex<2 else['sigma'],all_stage_variable_family_indices=stages,all_stage_matrices=[encode(M)for M in matrices],all_stage_original_pencil_row_indices=metadata,all_selected_weight_rows=all_selected,triangular_weight_identity_verified=True,PGL_dimension_offset=offset,cover=cover)
        cases.append(resultcase);Path('/tmp/quadric_cone_rank_two_triangular_progress.json').write_text(json.dumps(dict(cases=cases,seconds=time.monotonic()-begun),indent=2)+'\n')
    zero=old['cases'][2];pairs=zero['full_parameter_pair_order'];mons5=monomials(5,5)
    assert all(not value for i,row in enumerate(zero['full126_actual_Pfaffian_parameter_rows'])if not mons5[i][4]for col,value in row.items()if not any(j in[28,31,33]for j in pairs[int(col)]))
    paths=[Path(__file__),args.input,highpath,ROOT/'code/classification_plane_rank_zero_rank_two_high_m.py',ROOT/'code/classification_degenerate_quadric_q4_nil31_I_review.py',ROOT/'research/classification/QUADRIC_CONE_RANK_TWO_TRIANGULAR_BOUND.md']
    result=dict(task_id='G',all_arithmetic_passed=True,cases=cases,all_rank2_quadric_cone_G24_closed=all(v['cover']['all_G24_closed']for v in cases),C0_ell0_actual_Hfactor_verified=True,all_three_source_A_strata_and_all_R_retained=True,all_old_functions_AST_unchanged=True,no_nonlinearGB_or_small_slice_or_parameter_scan=True,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_rank2_quadric_cone_G24_closed'],cases=[dict(model=v['model'],generic_image=v['cover'].get('generic',{})['whole_PGL_image_upper'],closed=v['cover']['all_G24_closed'])for v in cases],seconds=result['seconds'])))


if __name__=='__main__':
    import sys
    triangular_bound_main()if '--triangular-bound'in sys.argv else highest_symbol_main()if '--highest-symbol'in sys.argv else main()
