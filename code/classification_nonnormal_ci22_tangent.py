#!/usr/bin/env python3
"""Full150 tangent/conormal sources on every non-normal integral CI22 orbit."""
import argparse,hashlib,itertools as it,json,resource,signal,time
from pathlib import Path
from classification_cubic_scroll_abelian_review import Q,PAIRS,add,scale,mul,diff,parse,rref,monomials,jacobi_and_divisor,determinant
from classification_quadric_cone_zero_source import empty,kernel
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();begun=time.monotonic()
    signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in ()).throw(TimeoutError('five whole nonnormalCI22 sources120s')));signal.alarm(120)
    resource.setrlimit(resource.RLIMIT_CPU,(115,118));resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,)*2)
    names=[f'x{i}'for i in range(5)];xs=[parse(v,names)for v in names];mons2=monomials(5,2);mons3=monomials(5,3);mons4=monomials(5,4);mons5=monomials(5,5);ix4={e:i for i,e in enumerate(mons4)};ix5={e:i for i,e in enumerate(mons5)};cols=[(ij,e)for ij in PAIRS for e in mons2]
    models=[('regular22_1','2*x0*x1+2*x2*x3+x4**2','x1**2+x3**2+x4**2',3),('nil32','2*x0*x2+x1**2+2*x3*x4','2*x1*x2+x4**2',4),('singularM1','x0*x1+x3**2','x0*x2+x4**2',4),('cone_node','2*x0*x1+x2**2+x3**2','x1**2+x2**2+2*x3**2',5),('cone_cusp','2*x0*x2+x1**2+x3**2','2*x1*x2+x3**2',5)]
    def rowpoly(P,index):return {index[e]:v for e,v in P.items()}
    def remainder(P,R,index):
        row=rowpoly(P,index)
        for pivot,v in R.items():
            if pivot in row:row=add(row,scale(v,-row[pivot]))
        return row
    def decode(v):
        P=empty()
        for value,((i,j),e)in zip(v,cols):
            if value:P[i][j][e]=value;P[j][i][e]=-value
        return P
    def encode_tensor(P):return {str((i,j)):[[list(e),str(v)]for e,v in sorted(P[i][j].items(),reverse=True)]for i,j in PAIRS}
    def sparse(rows):return [{str(j):str(v)for j,v in row.items()}for row in rows]
    E=[]
    for a,b in it.product(range(5),repeat=2):
        P=empty()
        for i,j in PAIRS:P[i][j]=add(mul(xs[i],xs[b])if j==a else{},scale(mul(xs[j],xs[b]),-1)if i==a else{})
        E.append([P[i][j].get(e,Q(0))for(i,j),e in cols])
    RE=rref(E);gauge_indices=sorted(RE);assert len(RE)==24;independent_E=[list(v.get(j,Q(0))for j in range(150))for v in RE.values()];gauge_minor=determinant([[v[i]for i in gauge_indices]for v in independent_E]);assert gauge_minor
    cases=[]
    for model,q1text,q2text,aut in models:
        q1,q2=parse(q1text,names),parse(q2text,names);I4=[mul(q,{e:Q(1)})for q in[q1,q2]for e in mons2];RI4=rref([rowpoly(v,ix4)for v in I4]);assert len(RI4)==29
        I25=[mul(mul(a,b),x)for a,b in[(q1,q1),(q1,q2),(q2,q2)]for x in xs];RI25=rref([rowpoly(v,ix5)for v in I25]);assert len(RI25)==15
        tangency={};conormal={}
        for column,((i,j),e)in enumerate(cols):
            mon={e:Q(1)};P=empty();P[i][j]=mon;P[j][i]=scale(mon,-1)
            for qi,q in enumerate([q1,q2]):
                V=[add(*(mul(P[a][b],diff(q,b))for b in range(5)))for a in range(5)]
                for a,b in PAIRS:
                    rem=remainder(add(mul(xs[a],V[b]),scale(mul(xs[b],V[a]),-1)),RI4,ix4)
                    for r,v in rem.items():tangency.setdefault((qi,a,b,r),{})[column]=v
            h12=add(*(mul(mul(P[a][b],diff(q1,a)),diff(q2,b))for a in range(5)for b in range(5)))
            h10=add(*(mul(diff(q1,a),P[a][0])for a in range(5)));h02=add(*(mul(P[0][b],diff(q2,b))for b in range(5)))
            numerator=add(mul(xs[0],h12),scale(mul(q2,h10),-2),scale(mul(q1,h02),-2));rem=remainder(numerator,RI25,ix5)
            for r,v in rem.items():conormal.setdefault(r,{})[column]=v
        tl=sorted(tangency);cl=sorted(conormal);T=[tangency[v]for v in tl];C=[conormal[v]for v in cl];M=T+C
        assert all(not sum(c*v[j]for j,c in row.items())for row in M for v in independent_E)
        RT,tf,KT=kernel(T,150);R,free,K=kernel(M,150)
        gauge=[{i:Q(1)}for i in gauge_indices];RR,frees,Kproj=kernel(M+gauge,150);assert len(K)==len(Kproj)+24
        assert all(not sum(c*v[j]for j,c in row.items())for row in M+gauge for v in Kproj)
        d=len(Kproj);assert d<=24
        family=[decode(v)for v in Kproj];JJ,DD=jacobi_and_divisor(family);pairs=list(it.combinations_with_replacement(range(d),2));RJ=rref(JJ)
        leading=[pairs[j]for j in RJ];loops=sorted({a for a,b in leading if a==b});used=set(loops);edges=[]
        for a,b in leading:
            if a!=b and a not in used and b not in used:edges.append((a,b));used.update([a,b])
        height=len(loops)+len(edges);upper=d-height;image=23+upper-aut
        allD=True;residual=[]
        for i,row in enumerate(DD):
            rem=dict(row)
            for pivot,v in RJ.items():
                if pivot in rem:rem=add(rem,scale(v,-rem[pivot]))
            if rem:allD=False;residual.append(dict(row=i,remainder={str(j):str(v)for j,v in rem.items()}))
        linear_factors=[i for i in range(5)if all(e[i]or not row for e,row in zip(mons5,DD))]
        qfactors=[]
        for text,q in[(q1text,q1),(q2text,q2)]:
            RQ=rref([rowpoly(mul(q,{e:Q(1)}),ix5)for e in mons3]);assert len(RQ)==35
            if all(not remainder({e:DD[i].get(col,Q(0))for i,e in enumerate(mons5)if DD[i].get(col,Q(0))},RQ,ix5)for col in range(len(pairs))):qfactors.append(text)
        emptybranch=allD or bool(linear_factors)or bool(qfactors);closed=emptybranch or image<=24
        cases.append(dict(model=model,q1=q1text,q2=q2text,full150_column_order=[dict(pair=ij,monomial=e)for ij,e in cols],I_degree4_rank29=True,I_square_degree5_rank15=True,
             all_tangency_row_labels=tl,full_tangency_matrix=sparse(T),tangency_raw_kernel_dimension=len(KT),tangency_projective_dimension=len(KT)-24,
             all_conormal_remainder_monomial_indices=cl,full_conormal_matrix=sparse(C),whole_tangent_conormal_raw_dimension=len(K),whole_projective_source_dimension=d,
             Euler24_verified=True,Euler_gauge_indices=gauge_indices,Euler_gauge_minor=str(gauge_minor),complete_projective_source_kernel=[list(map(str,v))for v in Kproj],complete_projective_free_original_indices=frees,all150_source_tensors=[encode_tensor(P)for P in family],
             full350_Jacobi_rows=sparse(JJ),full126_actualPf_rows=sparse(DD),parameter_pair_order=pairs,full_Jacobi_coefficient_rowspace_basis=sparse(list(RJ.values())),Jacobi_rowspace_rank=len(RJ),
             exact_lex_leading_monomials=leading,mandatory_square_variable_indices=loops,selected_pairwise_disjoint_edges=edges,initial_ideal_height_lower=height,whole_Poisson_source_dimension_upper=upper,
             actual_source_preserving_group_dimension_lower=aut,whole_projective_PGL_image_upper=image,all_actual_D_coefficients_in_original_Jacobi_rowspace=allD,actual_D_rowspace_nonmembership_residuals=residual,
             all_source_actual_D_fixed_coordinate_factors=linear_factors,all_source_actual_D_fixed_quadric_factors=qfactors,irreducible_quintic_branch_empty=emptybranch,G24_closed=closed))
        Path('/tmp/nonnormal_ci22_tangent_progress.json').write_text(json.dumps(dict(cases=cases,seconds=time.monotonic()-begun),indent=2)+'\n')
    paths=[Path(__file__),ROOT/'code/classification_cubic_scroll_abelian_review.py',ROOT/'code/classification_quadric_cone_zero_source.py',ROOT/'research/classification/NONNORMAL_CI22_RANK_TWO_SOURCE_ENTRANCE.md',ROOT/'research/classification/REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_NONDOMINANCE.md',ROOT/'research/classification/ALL_DEGENERATE_CI22_ABELIAN_CONORMAL_EXCLUSION.md']
    result=dict(task_id='G',all_arithmetic_passed=True,all_five_nonnormal_CI22_orbits_G24_closed=all(v['G24_closed']for v in cases),cases=cases,scope='Full150 tangent/conormal source on all five nonnormalCI22 orbits; no anticanonicalnormalization restriction ansatz. Exact necessary source/initial monomial bound only; failed Gbounds remain explicit.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=result['all_five_nonnormal_CI22_orbits_G24_closed'],cases=[dict(model=v['model'],tangent=v['tangency_projective_dimension'],source=v['whole_projective_source_dimension'],height=v['initial_ideal_height_lower'],sourceupper=v['whole_Poisson_source_dimension_upper'],image=v['whole_projective_PGL_image_upper'],empty=v['irreducible_quintic_branch_empty'])for v in cases],seconds=result['seconds'])))
def support_completion_main():
    import ast
    import sympy as sp
    ap=argparse.ArgumentParser();ap.add_argument('--support-completion',action='store_true');ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();begun=time.monotonic()
    signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in ()).throw(TimeoutError('whole five support completion60s')));signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(1024*1024**2,)*2)
    src=json.loads(args.input.read_text());frozen={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(src['producer_source_bytes']).body if isinstance(v,ast.FunctionDef)};now={v.name:ast.dump(v,include_attributes=False)for v in ast.parse(Path(__file__).read_text()).body if isinstance(v,ast.FunctionDef)};assert all(now[k]==v for k,v in frozen.items())
    x=sp.symbols('x0:5');mons5=monomials(5,5)
    def setup(case):
        p=sp.symbols('p0:'+str(case['whole_projective_source_dimension']));pairs=case['parameter_pair_order']
        def ev(row,sub):return sp.factor(sum(sp.Rational(v)*sub[p[i]]*sub[p[j]]for col,v in row.items()for i,j in[pairs[int(col)]]))
        E=[sp.factor(sum(sp.Rational(v)*p[i]*p[j]for col,v in row.items()for i,j in[pairs[int(col)]]))for row in case['full_Jacobi_coefficient_rowspace_basis']]
        return p,ev,E
    def test_map(case,p,ev,E,sub):
        assert set(sub)==set(p)
        return dict(parameter_map={str(k):str(v)for k,v in sub.items()},all_equivalent_original_Jacobi_residuals=[str(sp.factor(v.subs(sub,simultaneous=True)))for v in E])
    regular=src['cases'][0];p,ev,E=setup(regular);reg=[]
    for name,keep in[('A',[0,3,4,5,6]),('B',[3,5,6,7,9])]:
        sub={v:(v if i in keep else sp.Integer(0))for i,v in enumerate(p)};assert all(ev(row,sub)==0 for row in regular['full350_Jacobi_rows'])
        coeff=[ev(row,sub)for row in regular['full126_actualPf_rows']];D=sp.expand(sum(value*sp.prod(v**n for v,n in zip(x,e))for value,e in zip(coeff,mons5)));factor=sp.factor_list(D);records=[]
        for f,m in factor[1]:
            degrees={sum(e[:5])for e,c in sp.Poly(f,*x,*p).terms()};assert len(degrees)==1;records.append(dict(factor=str(f),multiplicity=int(m),spatial_degree=next(iter(degrees))))
        assert sp.expand(D-factor[0]*sp.prod(f**m for f,m in factor[1]))==0
        positive=[v for v in records if v['spatial_degree']];reducible=(D==0)or sum(v['multiplicity']for v in positive)>=2
        reg.append(dict(branch=name,complete_linear_source_indices=keep,all_original350_zero_verified=True,actual_D=str(D),rational_factorization_scalar=str(factor[0]),all_factors=records,actual_D_zero_or_positive_degree_product=reducible))
    nil=src['cases'][1];p,ev,E=setup(nil);A,B,G,Z,H,F=sp.symbols('A B G Z H F');nil_maps=[]
    def base():return {v:sp.Integer(0)for v in p}
    sub=base();sub.update({p[10]:A,p[8]:B,p[11]:G,p[2]:Z,p[7]:H,p[4]:3*A/2,p[3]:3*B/2,p[5]:B**2/(8*A)})
    sub[p[1]]=-(G+sub[p[5]])*B/(6*A);sub[p[9]]=-2*Z-G*B/A-sub[p[5]]*B/(2*A);sub[p[0]]=(-12*sub[p[1]]*G-6*sub[p[1]]*sub[p[5]]+8*(G+sub[p[5]])*Z+10*(G+sub[p[5]])*sub[p[9]])/(12*A)
    nil_maps.append(dict(case='p10_nonzero',inverted='A',source_upper=5,**test_map(nil,p,ev,E,sub)))
    sub=base();sub.update({p[4]:A,p[8]:B,p[5]:F,p[11]:-F,p[7]:H,p[9]:Z,p[3]:2*B,p[2]:-2*Z});nil_maps.append(dict(case='p10_zero_p4_nonzero',inverted='A',source_upper=5,**test_map(nil,p,ev,E,sub)))
    sub=base();sub.update({p[8]:B,p[3]:2*B,p[5]:F,p[11]:-F,p[7]:H,p[9]:Z,p[2]:-2*Z});nil_maps.append(dict(case='core_line0_nonzero',inverted='B',source_upper=4,**test_map(nil,p,ev,E,sub)))
    for line,t,r,p2,p9,p0 in[(1,4*B,-3*B,2*A,-A/2,A**2/(2*B)),(2,5*B/2,-3*B/2,7*A/3,-4*A/3,10*A**2/(9*B))]:
        sub=base();sub.update({p[8]:B,p[3]:t,p[6]:r,p[1]:A,p[2]:p2,p[9]:p9,p[0]:p0,p[7]:H});nil_maps.append(dict(case='core_line'+str(line)+'_nonzero',inverted='B',source_upper=3,**test_map(nil,p,ev,E,sub)))
    sub=base();sub.update({p[5]:F,p[11]:-F,p[0]:A,p[2]:B,p[7]:H,p[9]:Z});nil_maps.append(dict(case='core_zero_f_nonzero_g_minus_f',inverted='F',source_upper=5,**test_map(nil,p,ev,E,sub)))
    sub=base();sub.update({p[11]:G,p[1]:B,p[0]:A,p[9]:Z,p[2]:(6*B-5*Z)/4,p[7]:-(6*B**2-3*B*Z)/(4*G)});nil_maps.append(dict(case='core_zero_f0_g_nonzero',inverted='G',source_upper=4,**test_map(nil,p,ev,E,sub)))
    for line in[0,1]:
        sub=base();sub.update({p[0]:A,p[2]:B,p[7]:H,p[9]:Z}if line==0 else{p[0]:A,p[2]:B,p[7]:H,p[1]:Z,p[9]:2*Z});nil_maps.append(dict(case='all_core_zero_final_factor'+str(line),inverted=None,source_upper=4,**test_map(nil,p,ev,E,sub)))
    node=src['cases'][3];JJ=[{int(j):Q(v)for j,v in row.items()}for row in node['full350_Jacobi_rows']];nonzero=[i for i,v in enumerate(JJ)if v];ncols=len(node['parameter_pair_order']);R=rref([[row.get(j,Q(0))for j in range(ncols)]+[Q(i==k)for k in nonzero]for i in nonzero for row in[JJ[i]]]);nodecert=[]
    for ri,row in enumerate(node['full126_actualPf_rows']):
        target={int(j):Q(v)for j,v in row.items()};rem=dict(target);coeff=[Q(0)]*len(nonzero)
        for pivot,v in R.items():
            if pivot>=ncols:continue
            if pivot in rem:
                c=rem[pivot];rem=add(rem,scale({j:q for j,q in v.items()if j<ncols},-c));coeff=[a+c*v.get(ncols+j,Q(0))for j,a in enumerate(coeff)]
        assert not rem and all(sum(c*JJ[i].get(j,Q(0))for c,i in zip(coeff,nonzero))==target.get(j,Q(0))for j in range(ncols))
        if target:nodecert.append(dict(actual_D_row=ri,original350_functional={str(i):str(c)for i,c in zip(nonzero,coeff)if c}))
    cusp=src['cases'][4];p,ev,E=setup(cusp);sub0={v:(sp.Integer(0)if i in[13,15]else v)for i,v in enumerate(p)};assert all(ev(row,sub0)==0 for row in cusp['full126_actualPf_rows']);cusp_maps=[]
    sub={v:sp.Integer(0)for v in p};sub.update({p[6]:F,p[9]:B,p[11]:G,p[12]:A,p[13]:H,p[4]:2*A,p[8]:-F,p[14]:-B});assert all(ev(row,sub)==0 for row in cusp['full350_Jacobi_rows']);cusp_maps.append(dict(open='p13_nonzero',source_upper=5,**test_map(cusp,p,ev,E,sub)))
    sub={v:sp.Integer(0)for v in p};sub.update({p[3]:F,p[10]:B,p[12]:A,p[14]:G,p[15]:H,p[1]:2*A/3,p[4]:4*A/3,p[7]:(-12*B+9*F)/2});assert all(ev(row,sub)==0 for row in cusp['full350_Jacobi_rows']);cusp_maps.append(dict(open='p15_nonzero',source_upper=5,**test_map(cusp,p,ev,E,sub)))
    paths=[Path(__file__),args.input.resolve(),ROOT/'research/classification/NONNORMAL_CI22_REDUCED_SUPPORT_REDUCTION.md'];closed=all(v['actual_D_zero_or_positive_degree_product']for v in reg)
    result=dict(task_id='G',all_arithmetic_passed=True,regular22_largest_two_whole_source_factors=reg,regular22_remaining_support_dimension_upper=3 if closed else 5,nil32_complete_branch_cover_maps=nil_maps,nil32_source_dimension_upper_from_whole_caseproof=5,singularM1_source_dimension_upper_from_whole_caseproof=5,cone_node_all_D_original350_functionals=nodecert,cone_cusp_both_core_coefficients_zero_actualD0=True,cone_cusp_complete_two_open_maps=cusp_maps,cone_cusp_rank4_source_dimension_upper=5,all_five_nonnormal_CI22_G24_closed=closed,
        scope='Exact fullbranch factor identities/whole rational supersource maps and original functionals. Mathematical completeness/dimension of all inverse branch cases is the referenced handwritten proof, not inferred from map substitution alone.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},producer_source_bytes=Path(__file__).read_text(),seconds=round(time.monotonic()-begun,6),maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with args.output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_arithmetic_passed=True,all_G24_closed=closed,regular=[dict(branch=v['branch'],factors=v['all_factors'])for v in reg],seconds=result['seconds'])))

if __name__=='__main__':
    import sys
    support_completion_main()if '--support-completion'in sys.argv else main()
