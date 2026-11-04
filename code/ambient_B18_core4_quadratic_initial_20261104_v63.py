#!/usr/bin/env python3
"""One exact cubic evaluation, seven core quadrics, no new source solve."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as sy
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CORE4_QUADRATIC_INITIAL_20261104_v63'
MON=[()]+[(i,) for i in range(4)]+list(itertools.combinations_with_replacement(range(4),2))+list(itertools.combinations_with_replacement(range(4),3))
POS={m:i for i,m in enumerate(MON)}
CONV=[[(i,j) for i,m in enumerate(MON) for j,n in enumerate(MON) if tuple(sorted(m+n))==r] for r in MON]
ZERO=(F(0),)*35
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x
def mul(a,b):return tuple(sum((a[i]*b[j] for i,j in terms if a[i] and b[j]),F(0)) for terms in CONV)
def eval3(nodes,need,bind):
    values={}
    for i in need:
        n=nodes[i];o=n[0]
        if o=='c':a=(F(n[1],n[2]),)+ZERO[1:]
        elif o=='v':a=bind[i]
        elif o=='n':a=tuple(-v for v in values[n[1]])
        elif o=='a':a=tuple(x+y for x,y in zip(values[n[1]],values[n[2]]))
        elif o=='m':a=mul(values[n[1]],values[n[2]])
        elif o=='i':
            b=values[n[1]];b0=b[0];assert b0
            h=(F(0),)+b[1:];h2=mul(h,h);h3=mul(h2,h)
            a=(1/b0,)+tuple(-h[k]/b0**2+h2[k]/b0**3-h3[k]/b0**4 for k in range(1,35))
        else:raise ValueError(o)
        values[i]=a
    return values
def divide_linear(poly):
    L=[F(5),F(-20),F(24),F(6)];remainder={m:poly[i] for i,m in enumerate(MON) if len(m)==3 and poly[i]};quot={}
    while any(3 in m for m in remainder):
        m=max((m for m in remainder if 3 in m),key=lambda m:(m.count(3),m));c=remainder[m]/6
        q=list(m);q.remove(3);q=tuple(q);quot[q]=quot.get(q,F(0))+c
        for j,lj in enumerate(L):
            target=tuple(sorted(q+(j,)));remainder[target]=remainder.get(target,F(0))-c*lj
            if not remainder[target]:del remainder[target]
    assert not remainder,('nonzero exact L remainder',remainder)
    result=tuple(quot.get(m,F(0)) if len(m)==2 else F(0) for m in MON)
    ljet=(F(0),*L,*([F(0)]*30))
    assert mul(ljet,result)==poly
    return result
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    pp=ROOT/'reports/AMBIENT_B18_CORE4_QUADRATIC_PREFLIGHT_20261104_v63/literal-adjoint98-core4-preflight.json'
    pr=pp.parent/'report.json';P=read(pp);assert read(pr)['all_passed']
    cp=ROOT/P['original_scalar_DAG_file'];C=read(cp);nodes=C['DAG_scalar_nodes']
    bind={int(i):tuple(map(F,a))+ZERO[15:] for i,a in P['all_needed_source_and_parameter2_bindings']}
    V=eval3(nodes,P['source_DAG_ancestor_nodes'],bind)
    equations=[V[i] for i in P['actual98_equation_nodes']]
    assert all(not any(a[:15]) for a in equations)
    lam=[[F(c) for c in row] for row in P['adjoint98_rows_lambda']]
    raw=[];corrected=[];quadrics=[]
    for i,expr in enumerate(P['seven_residual_linear_combinations']):
        f=tuple(sum(F(c)*V[int(r)][j] for r,c in expr) for j in range(35))
        assert not any(f[:15])
        g=tuple(f[j]-sum(lam[i][k]*equations[k][j] for k in range(98)) for j in range(35))
        assert not any(g[:15]);q=divide_linear(g)
        raw.append(f);corrected.append(g);quadrics.append(q)
    qmat=[[a[i] for i,m in enumerate(MON) if len(m)==2] for a in quadrics]
    M=sy.Matrix([[sy.Rational(c.numerator,c.denominator) for c in row] for row in qmat]);rank=M.rank()
    cols=list(M.rref()[1]);rows=list(M.T.rref()[1]);minor=F(M[rows,cols].det()) if rank else F(1)
    zz=sy.symbols('x112 x117 tauA tauB');expressions=[]
    for q in quadrics:
        expression=sum((sy.Rational(c.numerator,c.denominator)*sy.prod(zz[j] for j in m) for c,m in zip(q,MON) if c),sy.Integer(0))
        expressions.append(str(sy.factor(expression)))
    proof={'monomials4_through3':list(map(list,MON)),'four_coordinate_order':list(map(str,zz)),
      'all98_unit_equations_evaluated_through3':equations,
      'all98_equations_zero_through2_verified':True,
      'all7_raw_residuals_through3':raw,'all7_after_exact_adjoint98_correction':corrected,
      'all7_core_quadrics_padded_through3':quadrics,
      'all7_exact_L1_product_identities_all35_coefficients_verified':True,
      'quadratic_monomials':list(map(list,MON[5:15])),'quadratic_coefficient_matrix7_by10':qmat,
      'quadratic_span_rank':rank,'independent_quadratic_rows':rows,'minor_columns':cols,'nonzero_minor':minor,
      'exact_factored_quadric_expressions':expressions,
      'no_unknown_source3_set_to_zero_in_result':'all unknown cubic corrections canceled by verified adjoint98',
      'no_source3_solve_Npair3_corecubic_source4_GB_or_depth_result':True}
    p=out/'seven-actual-core-quadrics.json';p.write_text(json.dumps(enc(proof),separators=(',',':'))+'\n');assert p.stat().st_size<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pp,pr,cp]
    rep={'id':STEM,'all_passed':True,'status':'seven_actual_core_quadratic_initial_forms_computed_from_cached_source2_and_adjoint98',
      'four_variables':4,'seven_core_functions':7,'quadratic_span_rank':rank,
      'quadratic_expressions':expressions,'exact_L_division_identities':7,'coefficient_checks_each_identity':35,
      'new_source3_solve_or_higher_sourcejet_used':False,'socle_depth_or_full_core_ideal_determined':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
