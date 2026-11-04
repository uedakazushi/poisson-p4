#!/usr/bin/env python3
"""Independent dictionary microjets of cached FULL9; no source solve or jet raise."""
import argparse, collections, hashlib, itertools, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
AUTHOR='AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63';ID=AUTHOR+'_REVIEW'
MON=[()]+[(i,) for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2))
MON11=[()]+[(i,) for i in range(11)]+list(itertools.combinations_with_replacement(range(11),2))
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b):
    c=dict(a)
    for m,v in b.items():
        c[m]=c.get(m,F(0))+v
        if not c[m]:del c[m]
    return c
def scale(a,c):return {m:c*v for m,v in a.items() if c*v}
def mul(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            if len(m)+len(n)>2:continue
            k=tuple(sorted(m+n));c[k]=c.get(k,F(0))+x*y
    return {m:v for m,v in c.items() if v}
def inv(a):
    c=a.get((),F(0));assert c
    h=scale(a,1/c);h.pop((),None)
    return scale(plus(plus({():F(1)},scale(h,-1)),mul(h,h)),1/c)
def decode(v,mons=MON):return {m:F(q) for m,q in zip(mons,v) if F(q)}
def vector(p,mons=MON):return [str(p.get(m,F(0))) for m in mons]
def coordinate(i):return {(i,):F(1)}
def composition(p,sub):
    out={}
    for mon,c in p.items():
        term={():c}
        for i in mon:term=mul(term,sub[i])
        out=plus(out,term)
    return out
def ancestors(nodes,roots):
    todo=list(roots);seen=set()
    while todo:
        i=todo.pop()
        if i in seen:continue
        seen.add(i);op,*a=nodes[i]
        if op in ('a','m'):todo.extend(a)
        elif op in ('n','i'):todo.append(a[0])
        else:assert op in ('c','v')
    return sorted(seen)
def evaluate(nodes,needed,bind):
    out={};inversions=0
    for i in needed:
        op,*a=nodes[i]
        if op=='v':p=bind[i]
        elif op=='c':p={():F(*a)} if a[0] else {}
        elif op=='a':p=plus(out[a[0]],out[a[1]])
        elif op=='n':p=scale(out[a[0]],-1)
        elif op=='m':p=mul(out[a[0]],out[a[1]])
        else:assert op=='i';p=inv(out[a[0]]);inversions+=1
        out[i]=p
    return out,inversions
def rank(rows):
    basis=[]
    for row in rows:
        v={i:F(x) for i,x in enumerate(row) if x}
        for k,w in basis:
            a=v.get(k,F(0))
            for i,x in w.items():
                v[i]=v.get(i,F(0))-a*x
                if not v[i]:v.pop(i,None)
        if v:
            k=min(v);a=v[k];basis.append((k,{i:x/a for i,x in v.items()}))
    return len(basis)
def determinant(A):
    a=[list(r) for r in A];d=F(1)
    for i in range(len(a)):
        j=next(j for j in range(i,len(a)) if a[j][i])
        if j!=i:a[j],a[i]=a[i],a[j];d=-d
        p=a[i][i];d*=p
        for j in range(i+1,len(a)):
            q=a[j][i]/p
            for k in range(i+1,len(a)):a[j][k]-=q*a[i][k]
    return d
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False);base=ROOT/'reports'/AUTHOR;hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);assert d['all_passed'];bind(p)
        for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;bind(p.parent/n)
        return d
    pr=audit(base/'report.json');assert pr['cached_source2_coefficients_each_family']==55
    for name in ('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW',
                 'AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW',
                 'AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW'):
        audit(ROOT/'reports'/name/'report.json')
    sourcep=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    cachep=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    modelp=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json'
    source,cache,model=map(read,[sourcep,cachep,modelp]);nodes=source['DAG_scalar_nodes']
    cert=read(base/'cached-source2-and-exact-E4-quadrics.json')
    assert cert['source2_monomials']==list(map(list,MON))
    xrows=source['physical_zero_base_rows9'];assert xrows==[27,108,109,112,113,117,118,119,61]
    variables={n[1]:i for i,n in enumerate(nodes) if n[0]=='v'};xids=[variables[f'x{r}'] for r in xrows]
    targets=[p['physical_row'] for p in model['all12_zero_physical_difference_columns']]
    assert targets==[18,19,28,29,60,67,72,63,64,68,73,74]
    common={};details={}
    for family in ('CA','CB'):
        rec=source['families'][family];saved=cache['both57_source2'][family];cc=cert['families'][family]
        expected_vars=[i for p in rec['coefficient_unknown_pairs'].values() for i in p]
        expected_vars += [rec['p_quadratic_relation'][k] for k in ('s','t')]+xids
        gen=rec['actual_zero_odd_I0_generators']+rec['actual_positive_odd_J0_generators']
        f9=[gen[i][1] for i in [0,1,2,3,4,5,6,8,9]]
        eq=rec['all48_zero_IFT_quadratic_coefficient_equations']+f9
        assert saved['unknown57']==expected_vars and saved['equations57']==eq and saved['F9_nodes']==f9
        assert cc['cached_source57_unknowns']==expected_vars and cc['own_F9_coordinate_order']==f9
        assert cc['cached_source2_all57']==saved['all57_source2jets']
        bindjet={i:decode(v) for i,v in saved['all57_source2jets']}
        assert set(bindjet)==set(expected_vars) and len(bindjet)==57
        assert all(() not in p for p in bindjet.values())
        pairs=rec['complete_normalized_physical_w_coefficients'][0]
        roots=eq+[rec['actual_tracefree_nu']]+[k for r in targets for k in pairs[r]]
        need=ancestors(nodes,roots);assert {i for i in need if nodes[i][0]=='v'}==set(bindjet)
        values,ninv=evaluate(nodes,need,bindjet)
        assert all(not values[r] if i<48 else values[r]==coordinate(i-48) for i,r in enumerate(eq))
        X=[bindjet[i] for i in xids]
        assert [vector(p) for p in X]==cc['physical_R0_readback_jets_in_own_F9']
        A=[[p.get((j,),F(0)) for j in range(9)] for p in X]
        assert [[str(x) for x in row] for row in A]==cc['R0_by_F9_unit_Jacobian']
        det=determinant(A);assert det==F(cc['unit9_determinant'])==({'CA':F(3,128),'CB':F(96)}[family])
        Ai=[list(map(F,row)) for row in cc['fixed_inverse9']]
        for P,Q in ((A,Ai),(Ai,A)):
            assert all(sum(P[i][k]*Q[k][j] for k in range(9))==int(i==j) for i in range(9) for j in range(9))
        inverse=[decode(v) for v in cc['inverse_F9_in_common_R0_through2']]
        # Both compositions, in a different sparse dictionary monomial algebra.
        assert [composition(p,inverse) for p in X]==[coordinate(i) for i in range(9)]
        assert [composition(p,X) for p in inverse]==[coordinate(i) for i in range(9)]
        nu=composition(values[rec['actual_tracefree_nu']],inverse)
        actual={r:[composition(values[k],inverse) for k in pairs[r]] for r in targets}
        assert vector(nu)==cc['nu_in_common_R0_through2']
        assert [[r,[vector(p) for p in actual[r]]] for r in targets]==cc['all12_trace_pairs_in_common_R0_through2']
        assert all(not p[1].get((),F(0)) for p in actual.values())
        common[family]=(nu,actual)
        details[family]=dict(original_ancestor_nodes=len(need),unit_inverse_nodes=ninv,
            all57_equations_all55_coefficients_checked=True,R0_by_F9_determinant=str(det),
            both9_inverse_compositions_checked=True,all12_common_physical_trace_pairs_checked=True)
    equations=[]
    for family,t in (('CA',9),('CB',10)):
        equations.append(plus(scale(common[family][0],-1),{(t,t):F(1)}))
    for r in targets:
        a,b=common['CA'][1][r];c,d=common['CB'][1][r]
        p=plus(a,scale(c,-1))
        # The two trace bases remain distinct. Only b1 contributes through2.
        p=plus(p,{(i,9):b.get((i,),F(0)) for i in range(9) if b.get((i,),F(0))})
        p=plus(p,{(i,10):-d.get((i,),F(0)) for i in range(9) if d.get((i,),F(0))})
        equations.append(p)
    assert [vector(p,MON11) for p in equations]==cert['all14_raw_equations_through2']
    linear=[[p.get((i,),F(0)) for i in range(11)] for p in equations]
    assert rank(linear)==2 and [[str(x) for x in r] for r in linear]==cert['linear_Jacobian14_by11']
    elim=[0,8];free=[i for i in range(11) if i not in elim]
    assert cert['linear_pivot_equations']==[0,1] and cert['linear_pivot_variables']==elim
    A=[[linear[i][j] for j in elim] for i in range(2)];det=A[0][0]*A[1][1]-A[0][1]*A[1][0]
    assert det==F(-20,9)==F(cert['linear_Schur_minor'])
    Ai=[[A[1][1]/det,-A[0][1]/det],[-A[1][0]/det,A[0][0]/det]]
    sub=[{} for _ in range(11)]
    for j,i in enumerate(free):sub[i]=coordinate(j)
    for k,i in enumerate(elim):
        sub[i]={(j,):-sum(Ai[k][h]*linear[h][f] for h in range(2)) for j,f in enumerate(free)}
        sub[i]={m:v for m,v in sub[i].items() if v}
    raw=[composition(p,sub) for p in equations]
    for k,i in enumerate(elim):
        sub[i]=plus(sub[i],total(scale({m:v for m,v in raw[h].items() if len(m)==2},-Ai[k][h]) for h in range(2)))
    # The producer leaves literal integer1 in free-coordinate jets while
    # its other coefficients are rational strings. Compare in QQ, not JSON types.
    assert sub==[decode(p) for p in cert['all11_variables_as_remaining9_through2']]
    after=[composition(p,sub) for p in equations]
    assert not after[0] and not after[1] and all(all(len(m)==2 for m in p) for p in after)
    assert [vector(p) for p in after]==cert['all14_equations_after_exact_linear_Schur_through2']
    quadmatrix=[[p.get(m,F(0)) for m in MON[10:]] for p in after]
    assert rank(quadmatrix)==5==cert['full_quadratic_initial_span_rank']
    L={(2,):F(5),(4,):F(-20),(7,):F(24),(8,):F(6)}
    expect=[mul(L,coordinate(i)) for i in (0,1,3,5,6)]
    expectedmatrix=[[p.get(m,F(0)) for m in MON[10:]] for p in expect]
    assert rank(expectedmatrix)==rank(expectedmatrix+quadmatrix)==5
    factors=cert['all_common_rational_linear_factors'];assert len(factors)==1
    fac=factors[0];assert list(map(F,fac['factor_linear_coefficients']))==[L.get((i,),F(0)) for i in range(9)]
    assert len(fac['linear_quotient_rows'])==14
    for p,row in zip(after,fac['linear_quotient_rows']):assert p==mul(L,{(i,):F(v) for i,v in enumerate(row) if F(v)})
    assert cert['remaining_variable_indices']==free
    assert cert['remaining_variable_names']==['x108','x109','x112','x113','x117','x118','x119','tau_CA','tau_CB']
    proof=dict(independent_sparse_dictionary_degree2_only=True,cached_source_equations=114,
        cached_monomials_per_unknown=55,source_details=details,all14_true_equations_independently_reassembled=True,
        distinct_family_trace_bases_preserved=True,linear_rank=2,Schur_minor=str(det),
        eliminated_physical_rows=[27,61],all_quadratic_coordinate_corrections_retained=True,
        full_quadratic_initial_rank=5,remaining_embedding_dimension=9,
        exact_quadratic_span='(5*x112-20*x117+24*tau_CA+6*tau_CB)*(x108,x109,x113,x118,x119)',
        quadratic_initial_is_complete_after_unit_linear_elimination=True,
        no_allorders_factorization_socle_depth_nilindex_or_D0_injectivity_claim=True,
        source_solver_author_worker_newjet_point_prime_GB_used=False)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/f'research/notes/{ID}.md',sourcep,cachep,modelp):bind(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,
        evidence_sha256={p.name:sha(p) for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Actual E4 embedding dimension9 and complete quadratic initial span5 only; no all-order ring, D0 injection, socle or depth assertion.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
def total(ps):
    q={}
    for p in ps:q=plus(q,p)
    return q
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
