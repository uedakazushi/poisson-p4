#!/usr/bin/env python3
"""ONE cached-source2 E4 quadratic/linear-Schur calculation, no new jets."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as sy
from ambient_B18_top_graph_telescoping_review_20261104_v63 import solve,closure
from ambient_B18_T3_zero_structural_20261104_v63 import MON,POS,qmul,eval2,insert_basis
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63'
ZERO=(F(0),)*55
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def enc(a):
    if isinstance(a,F):return str(a)
    if isinstance(a,(list,tuple)):return [enc(x) for x in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return a
def compose(a,sub):
    out=[a[0]]+[F(0)]*54
    for i in range(9):
        if a[1+i]:
            for k,c in enumerate(sub[i]):out[k]+=a[1+i]*c
    for j,m in enumerate(MON[10:],10):
        if a[j]:
            prod=qmul(sub[m[0]],sub[m[1]])
            for k,c in enumerate(prod):out[k]+=a[j]*c
    return tuple(out)
def row_basis(M):
    b=[];p=[]
    for i,v in enumerate(M):
        if insert_basis(b,v) is not None:p.append(i)
    return p,b
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    modelp=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json'
    gates=[modelp.parent/'report.json',ROOT/'reports/AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    S,J,E=map(read,[sp,jp,modelp]);nodes=S['DAG_scalar_nodes'];rows=S['physical_zero_base_rows9']
    x=[next(i for i,n in enumerate(nodes) if n==['v',f'x{r}']) for r in rows]
    targets=[d['physical_row'] for d in E['all12_zero_physical_difference_columns']]
    assert targets==[18,19,28,29,60,67,72,63,64,68,73,74]
    data={};common={}
    for family in ['CA','CB']:
        cache=J['both57_source2'][family];rec=S['families'][family]
        bind={i:tuple(map(F,p)) for i,p in cache['all57_source2jets']}
        assert len(bind)==57 and all(len(p)==55 for p in bind.values())
        pair=rec['complete_normalized_physical_w_coefficients'][0]
        roots=cache['equations57']+[rec['actual_tracefree_nu']]+[i for r in targets for i in pair[r]]
        needed=sorted(closure(nodes,roots));assert {i for i in needed if nodes[i][0]=='v'}<=set(bind)
        val=eval2(nodes,needed,bind)
        desired=[ZERO]*48
        for k in range(9):a=[F(0)]*55;a[1+k]=1;desired.append(tuple(a))
        assert all(val[r]==desired[i] for i,r in enumerate(cache['equations57']))
        X=[bind[i] for i in x];L=[list(p[1:10]) for p in X]
        Li,det9=solve(L,[[F(i==j) for j in range(9)] for i in range(9)]);assert det9
        lin=[(F(0),*r,*([F(0)]*45)) for r in Li]
        q=[compose((F(0),*([F(0)]*9),*p[10:]),lin) for p in X]
        sub=[]
        for i in range(9):
            a=list(lin[i])
            for k in range(10,55):a[k]=-sum((Li[i][j]*q[j][k] for j in range(9)),F(0))
            sub.append(tuple(a))
        check=[compose(p,sub) for p in X]
        assert all(p==tuple([F(0)]+[F(i==j) for j in range(9)]+[F(0)]*45) for i,p in enumerate(check))
        result={r:[compose(val[k],sub) for k in pair[r]] for r in targets}
        nu=compose(val[rec['actual_tracefree_nu']],sub)
        common[family]={'outputs':result,'nu':nu}
        data[family]={'cached_source57_unknowns':cache['unknown57'],
          'cached_source2_all57':cache['all57_source2jets'],'all57_equations_all55_coefficients_zero_or_coordinate_checked':True,
          'own_F9_coordinate_order':cache['F9_nodes'],'physical_R0_rows9':rows,
          'physical_R0_readback_jets_in_own_F9':X,'R0_by_F9_unit_Jacobian':L,'unit9_determinant':det9,
          'fixed_inverse9':Li,'inverse_F9_in_common_R0_through2':sub,
          'all9_inverse_composition_identities_through2_checked':True,
          'nu_in_common_R0_through2':nu,'all12_trace_pairs_in_common_R0_through2':[[r,result[r]] for r in targets]}
    vars11=[f'x{r}' for r in rows]+['tau_CA','tau_CB']
    mon11=[()]+[(i,) for i in range(11)]+list(itertools.combinations_with_replacement(range(11),2));p11={m:i for i,m in enumerate(mon11)}
    def lift(a):
        r=[F(0)]*len(mon11)
        for m,c in zip(MON,a):r[p11[m]]=c
        return r
    equations=[];definitions=[]
    for family,t in [('CA',9),('CB',10)]:
        a=[-c for c in lift(common[family]['nu'])];a[p11[(t,t)]]+=1
        equations.append(a);definitions.append(f'tau_{family}^2-nu_{family}')
    for r in targets:
        a,b=common['CA']['outputs'][r];c,d=common['CB']['outputs'][r]
        q=lift(tuple(x-y for x,y in zip(a,c)))
        assert b[0]==d[0]==0
        for i in range(9):q[p11[(i,9)]]+=b[1+i];q[p11[(i,10)]]-=d[1+i]
        equations.append(q);definitions.append(f'physical{r}_CA_minus_CB')
    assert all(q[0]==0 for q in equations)
    linear=[q[1:12] for q in equations];selected,lb=row_basis(linear);rank=len(lb)
    assert rank<=11
    pivotcols,_=row_basis([list(c) for c in zip(*[linear[i] for i in selected])])
    assert len(pivotcols)==rank
    free=[i for i in range(11) if i not in pivotcols]
    nq=len(free);assert nq==9,('unexpected_embedding_dimension',rank,nq)
    A=[[linear[i][j] for j in pivotcols] for i in selected]
    inverse,det=solve(A,[[F(i==j) for j in range(rank)] for i in range(rank)])
    sub11=[]
    for i in range(11):
        a=[F(0)]*55
        if i in free:a[1+free.index(i)]=1
        else:
            k=pivotcols.index(i)
            for j,f in enumerate(free):a[1+j]=-sum((inverse[k][h]*linear[e][f] for h,e in enumerate(selected)),F(0))
        sub11.append(tuple(a))
    def compose11(poly,sub):
        out=[F(0)]*55
        for m,c in zip(mon11,poly):
            if not c:continue
            if not m:out[0]+=c
            else:
                p=sub[m[0]] if len(m)==1 else qmul(sub[m[0]],sub[m[1]])
                for k,v in enumerate(p):out[k]+=c*v
        return tuple(out)
    raw=[compose11(p,sub11) for p in equations]
    for k,variable in enumerate(pivotcols):
        a=list(sub11[variable])
        for j in range(10,55):a[j]=-sum((inverse[k][h]*raw[e][j] for h,e in enumerate(selected)),F(0))
        sub11[variable]=tuple(a)
    after=[compose11(p,sub11) for p in equations]
    assert all(not any(q[:10]) for q in after)
    assert all(not any(after[e]) for e in selected)
    quad=[list(p[10:]) for p in after];pivots,qbasis=row_basis(quad);qr=len(qbasis)
    # Degree2 rational factor arithmetic only, not a Groebner calculation.
    zz=sy.symbols('z0:9');polys=[sum((sy.Rational(c.numerator,c.denominator)*zz[i]*zz[j] for (i,j),c in zip(MON[10:],q)),sy.Integer(0)) for q in quad]
    factor_data=[];common_factors=[]
    if pivots:
        unit,factors=sy.factor_list(polys[pivots[0]],*zz)
        for fac,power in factors:
            if sy.Poly(fac,*zz).total_degree()!=1:continue
            quotients=[]
            for p in polys:
                try:v=sy.Poly(p,*zz).exquo(sy.Poly(fac,*zz)).as_expr()
                except sy.polys.polyerrors.ExactQuotientFailed:break
                quotients.append(v)
            if len(quotients)==len(polys):
                lc=[F(sy.Poly(fac,*zz).coeff_monomial(z)) for z in zz]
                ll=[[F(sy.Poly(v,*zz).coeff_monomial(z)) for z in zz] for v in quotients]
                qi,qb=row_basis(ll)
                common_factors.append({'factor_linear_coefficients':lc,'linear_quotient_rows':ll,
                  'linear_quotient_rank':len(qb),'linear_quotient_independent_rows':qi})
        for i in pivots:factor_data.append({'equation_index':i,'factorization':str(sy.factor(polys[i]))})
    proof={'source2_monomials':list(map(list,MON)),'families':data,'variables11':vars11,
      'monomials11_through2':list(map(list,mon11)),'all14_equation_definitions':definitions,
      'all14_raw_equations_through2':equations,'linear_Jacobian14_by11':linear,
      'linear_rank':rank,'linear_pivot_equations':selected,'linear_pivot_variables':pivotcols,
      'linear_Schur_minor':det,'remaining_variable_indices':free,
      'remaining_variable_names':[vars11[i] for i in free],
      'all11_variables_as_remaining9_through2':sub11,
      'all14_equations_after_exact_linear_Schur_through2':after,
      'full_quadratic_initial_span_rank':qr,'quadratic_independent_equation_indices':pivots,
      'independent_quadric_factorizations':factor_data,'all_common_rational_linear_factors':common_factors,
      'full_ring_socle_or_depth_not_inferred_from_quadrics':True}
    ep=out/'cached-source2-and-exact-E4-quadrics.json';text=json.dumps(enc(proof),separators=(',',':'))+'\n'
    assert len(text.encode())<32*1024*1024;ep.write_text(text)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',sp,jp,modelp,*gates,
      ROOT/'code/ambient_B18_T3_zero_structural_20261104_v63.py',ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    rep={'id':STEM,'all_passed':True,'status':'actual_E4_full_quadratic_initial_after_exact_linear_Schur_only',
      'cached_source2_coefficients_each_family':55,'full_physical_R0_inverse_coordinate_correction_used':True,
      'all57_source_equations_checked_in_both_caches':True,'defining_equations':14,'ambient_variables_before_Schur':11,
      'linear_rank':rank,'embedding_dimension':nq,'full_quadratic_initial_rank':qr,
      'common_rational_linear_factor_count':len(common_factors),
      'full_ring_socle_or_depth_claimed':False,'new_source_jets_48_IFT_or_old_r_g_composition_used':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps({k:rep[k] for k in ['all_passed','linear_rank','embedding_dimension','full_quadratic_initial_rank','common_rational_linear_factor_count','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
