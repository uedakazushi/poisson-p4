#!/usr/bin/env python3
"""Independent cached-source2 composition, sparseAD and dictionary cubic jets."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_zero_normalizer_quadratic_review_20261104_v63 import decode,composition,coordinate,ancestors,determinant,rank
from s61_sparse_generic5_full134_normal_circuit_review_20261104_v63 import plus as freeadd,times as freemul,neg as freeneg,const,symbol
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_CORE4_QUADRATIC_20261104_v63_REVIEW'
M9=[()]+[(i,)for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2))
M4=[()]+[(i,)for i in range(4)]+list(itertools.combinations_with_replacement(range(4),2))
M3=M4+list(itertools.combinations_with_replacement(range(4),3))
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
    c=dict(a)
    for m,v in b.items():
        c[m]=c.get(m,F(0))+v
        if not c[m]:del c[m]
    return c
def scale(a,c):return{m:c*v for m,v in a.items()if c*v}
def mul(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            if len(m)+len(n)<=3:
                k=tuple(sorted(m+n));c[k]=c.get(k,F(0))+x*y
    return{m:v for m,v in c.items()if v}
def invert(a):
    c=a.get((),F(0));assert c
    h=scale(a,1/c);h.pop((),None);h2=mul(h,h)
    return scale(add(add(add({():F(1)},scale(h,-1)),h2),scale(mul(h2,h),-1)),1/c)
def vector(a,mons=M3):return list(map(str,(a.get(m,F(0))for m in mons)))
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    base=ROOT/'reports';hashes={}
    def bind(p):hashes[str(p.relative_to(ROOT))]=sha(p)
    def audit(p):
        d=read(p);assert d['all_passed'];bind(p)
        for n,h in d.get('input_sha256',{}).items():assert sha(ROOT/n)==h,n;hashes[n]=h
        for n,h in d.get('evidence_sha256',{}).items():assert sha(p.parent/n)==h,n;bind(p.parent/n)
        return d
    for n in ['AMBIENT_B18_CORE4_QUADRATIC_PREFLIGHT_20261104_v63','AMBIENT_B18_CORE4_QUADRATIC_INITIAL_20261104_v63',
              'AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63_REVIEW','AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63_REVIEW',
              'AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW']:
        audit(base/n/'report.json')
    P=read(base/'AMBIENT_B18_CORE4_QUADRATIC_PREFLIGHT_20261104_v63/literal-adjoint98-core4-preflight.json')
    Q=read(base/'AMBIENT_B18_CORE4_QUADRATIC_INITIAL_20261104_v63/seven-actual-core-quadrics.json')
    C=read(ROOT/P['original_scalar_DAG_file']);nodes=C['DAG_scalar_nodes']
    E=read(base/'AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json')
    J=read(base/'AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json')
    S=read(base/'AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json')
    prefix=S['DAG_scalar_nodes'];assert nodes[:len(prefix)]==prefix
    unknown=C['unit_schur_unknown_nodes98'];equations=C['unit_schur_equations98'];tau=C['tau_nodes']
    assert len(unknown)==len(equations)==98 and unknown==P['actual98_unknown_nodes'] and equations==P['actual98_equation_nodes']
    expected=[];eq96=[]
    for family in ('CA','CB'):
        rec=S['families'][family];expected += [i for p in rec['coefficient_unknown_pairs'].values()for i in p]+[rec['p_quadratic_relation'][k]for k in ('s','t')]
        eq96+=rec['all48_zero_IFT_quadratic_coefficient_equations']
    variables={n[1]:i for i,n in enumerate(nodes)if n[0]=='v'}
    assert unknown==expected+[variables['x27'],variables['x61']] and equations[:96]==eq96
    rows=[18,19,28,29,60,67,72,63,64,68,73,74];assert C['actual_all12_difference_rows']==rows
    # Literal source-root binding of both monic equations and all12 outputs.
    memo={}
    def actual(i):
        if i in memo:return memo[i]
        op,*a=nodes[i]
        if i<len(prefix):v=const(F(*a))if op=='c'else symbol(('old',i))
        elif op=='v':v=symbol(('tau',i))
        elif op=='c':v=const(F(*a))
        elif op=='a':v=freeadd(actual(a[0]),actual(a[1]))
        elif op=='n':v=freeneg(actual(a[0]))
        else:assert op=='m';v=freemul(actual(a[0]),actual(a[1]))
        memo[i]=v;return v
    for i,family in enumerate(('CA','CB')):
        t=actual(tau[family]);nu=actual(S['families'][family]['actual_tracefree_nu'])
        assert actual(equations[96+i])==freeadd(freemul(t,t),freeneg(nu))
    for row,root in zip(rows,C['all12_actual_difference_nodes']):
        pair={f:S['families'][f]['complete_normalized_physical_w_coefficients'][0][row]for f in ('CA','CB')}
        lhs=freeadd(actual(pair['CA'][0]),freemul(actual(pair['CA'][1]),actual(tau['CA'])))
        rhs=freeadd(actual(pair['CB'][0]),freemul(actual(pair['CB'][1]),actual(tau['CB'])))
        assert actual(root)==freeadd(lhs,freeneg(rhs))
    selected=[5,6,7,8,9,10,11]
    expr=[[(C['all12_actual_difference_nodes'][i],F(1))]+([(C['all12_actual_difference_nodes'][4],F(-1))]if i==6 else[])for i in selected]
    assert [[(int(r),F(c))for r,c in row]for row in P['seven_residual_linear_combinations']]==expr
    assert P['residuals_physical_rows']==[67,72,63,64,68,73,74]
    # The zero2jet premise is on all nine Schur variables, before core restriction.
    for i in selected:
        a=list(map(F,E['all14_equations_after_exact_linear_Schur_through2'][i+2]))
        if i==6:a=[x-F(y)for x,y in zip(a,E['all14_equations_after_exact_linear_Schur_through2'][6])]
        assert not any(a)
    keep=[2,4,7,8];linear=[coordinate(keep.index(j))if j in keep else{}for j in range(9)]
    sub=[composition(decode(a,M9),linear)for a in E['all11_variables_as_remaining9_through2']]
    bound={}
    for family in ('CA','CB'):
        cache=J['both57_source2'][family];f=E['families'][family]
        own=[composition(decode(a,M9),sub[:9])for a in f['inverse_F9_in_common_R0_through2']]
        proof=P['source2_cache_reuse'][family]
        assert own==[decode(a,M4)for a in proof['family_F9_in_core4_through2']]
        data={i:composition(decode(a,M9),own)for i,a in cache['all57_source2jets']}
        assert data=={i:decode(a,M4)for i,a in proof['all57_source2_in_core4']}
        for i,a in data.items():
            if i in bound:assert bound[i]==a
            bound[i]=a
    params=C['actual9_parameter_nodes'];assert params==[variables[f'x{r}']for r in [108,109,112,113,117,118,119]]+[tau['CA'],tau['CB']]
    for i,node in enumerate(params):
        if node in bound:assert bound[node]==linear[i]
        bound[node]=linear[i]
    assert bound=={i:decode(a,M4)for i,a in P['all_needed_source_and_parameter2_bindings']} and len(bound)==107
    roots=equations+[r for row in expr for r,c in row];needed=ancestors(nodes,roots)
    assert needed==P['source_DAG_ancestor_nodes'] and len(needed)==17085
    # Independent forward sparse first derivatives at the literal origin.
    values={};gradient={};U=set(unknown)
    for i in needed:
        op,*a=nodes[i]
        if op=='c':v=F(*a);g={}
        elif op=='v':v=F(0);g={i:F(1)}if i in U else{}
        elif op=='a':v=values[a[0]]+values[a[1]];g=add(gradient[a[0]],gradient[a[1]])
        elif op=='n':v=-values[a[0]];g=scale(gradient[a[0]],-1)
        elif op=='m':v=values[a[0]]*values[a[1]];g=add(scale(gradient[a[0]],values[a[1]]),scale(gradient[a[1]],values[a[0]]))
        else:assert op=='i'and values[a[0]];v=1/values[a[0]];g=scale(gradient[a[0]],-v*v)
        values[i]=v;gradient[i]=g
    A=[[gradient[r].get(v,F(0))for v in unknown]for r in equations]
    B=[[sum(c*gradient[r].get(v,F(0))for r,c in row)for v in unknown]for row in expr]
    assert A==[list(map(F,r))for r in P['actual98_origin_Jacobian']]
    assert B==[list(map(F,r))for r in P['seven_raw_unknown_gradients']]
    det=determinant(A);assert det==F(P['unit98_determinant'])==F(-12524124635136,78125)
    lam=[list(map(F,r))for r in P['adjoint98_rows_lambda']];assert len(lam)==7 and all(len(r)==98 for r in lam)
    assert all(sum(lam[i][k]*A[k][j]for k in range(98)if lam[i][k])==B[i][j]for i in range(7)for j in range(98))
    assert all(not any(r)for r in lam) and all(not any(r)for r in B)
    # Four-variable dictionary jets; omitted cubic source unknowns are not
    # interpreted as exact solutions. Their influence is canceled by lambda J=B.
    jets={}
    for i in needed:
        op,*a=nodes[i]
        if op=='c':v={():F(*a)}if a[0]else{}
        elif op=='v':v=bound[i]
        elif op=='a':v=add(jets[a[0]],jets[a[1]])
        elif op=='n':v=scale(jets[a[0]],-1)
        elif op=='m':v=mul(jets[a[0]],jets[a[1]])
        else:assert op=='i';v=invert(jets[a[0]])
        jets[i]=v
    eqjets=[jets[r]for r in equations]
    assert all(all(len(m)==3 for m in p)for p in eqjets)
    assert [vector(p)for p in eqjets]==Q['all98_unit_equations_evaluated_through3']
    L={(0,):F(5),(1,):F(-20),(2,):F(24),(3,):F(6)};quadrics=[]
    for i,row in enumerate(expr):
        raw={}
        for r,c in row:raw=add(raw,scale(jets[r],c))
        assert all(len(m)==3 for m in raw)and vector(raw)==Q['all7_raw_residuals_through3'][i]
        corrected=dict(raw)
        for j,c in enumerate(lam[i]):
            if c:corrected=add(corrected,scale(eqjets[j],-c))
        assert vector(corrected)==Q['all7_after_exact_adjoint98_correction'][i]
        quad=decode(Q['all7_core_quadrics_padded_through3'][i],M3);assert all(len(m)==2 for m in quad)
        assert mul(L,quad)==corrected;quadrics.append(quad)
    qmat=[[q.get(m,F(0))for m in M4[5:]]for q in quadrics]
    assert qmat==[list(map(F,r))for r in Q['quadratic_coefficient_matrix7_by10']]and rank(qmat)==2
    assert all(not q for q in quadrics[2:])and quadrics[0]and quadrics[1]
    rr=Q['independent_quadratic_rows'];cc=Q['minor_columns'];minor=determinant([[qmat[i][j]for j in cc]for i in rr])
    assert rr==[0,1]and cc==[0,1]and minor==F(Q['nonzero_minor'])==F(703,82944)
    proof=dict(actual_source_prefix_and_all12_output_difference_nodes_bound=True,actual96_plus2_monic_equations_and98_unknowns_bound=True,
        all107_fourvariable_source2_bindings_rebuilt_from_both57_cache=True,source_ancestor_nodes=17085,
        independent_sparse_AD98_Jacobian_and686_lambda_identities=True,unit98_determinant=str(det),
        all7_lambda_and_raw_unknown_gradients_are_actually_zero=True,unknown_source3_independence_proved_not_assumed=True,
        all98_equation_jets_zero_through2=True,all98_saved_cubic_residuals_reproduced=True,
        all7_raw_and_corrected_cubic_jets_reproduced=True,all245_L1_product_coefficients_checked=True,
        true7_core_quadratic_initial_span_rank=2,only_first_two_quadric_nonzero=True,nonzero_minor=str(minor),
        no_L2_L3_or_five_G_second_acceleration_needed=True,no_higher_source_solve_fullideal_socle_or_depth_claim=True)
    (out/'independent-proof.json').write_text(json.dumps(proof,indent=2)+'\n')
    for p in (Path(__file__),ROOT/'code/ambient_B18_zero_normalizer_quadratic_review_20261104_v63.py',ROOT/'code/s61_sparse_generic5_full134_normal_circuit_review_20261104_v63.py',ROOT/f'research/notes/{ID}.md'):bind(p)
    report=dict(id=ID,all_passed=True,input_sha256=hashes,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start,
        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,scope='Actual core4 complete quadratic initial space of rank2 only; no higher-order ideal, nilpotence, socle or depth conclusion.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(all_passed=True,seconds=report['seconds'],rss=report['max_rss_kib'],report_sha256=sha(out/'report.json'))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.output)
