#!/usr/bin/env python3
"""ONE145-slot core-cubic computation from cached source2; no source3/4 solve."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as sy
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CORE_CUBIC_INITIAL_20261104_v63'
MON=[()]+[(i,) for i in range(4)]+list(itertools.combinations_with_replacement(range(4),2))+list(itertools.combinations_with_replacement(range(4),3))+list(itertools.combinations_with_replacement(range(4),4))
POS={m:i for i,m in enumerate(MON)}
CONV70=[[(i,j) for i,m in enumerate(MON) for j,n in enumerate(MON) if tuple(sorted(m+n))==r] for r in MON]
CONV15=[[(i,j) for i,m in enumerate(MON[:15]) for j,n in enumerate(MON[:15]) if tuple(sorted(m+n))==r] for r in MON[:15]]
ZERO=(F(0),)*70
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x
def conv(a,b,table):return tuple(sum((a[i]*b[j] for i,j in row if a[i] and b[j]),F(0)) for row in table)
def mul145(a,b):
    base=conv(a[:70],b[:70],CONV70);dual=[]
    for j in range(5):
        k=70+15*j;x=conv(a[:15],b[k:k+15],CONV15);y=conv(a[k:k+15],b[:15],CONV15)
        dual.extend(u+v for u,v in zip(x,y))
    return base+tuple(dual)
def inverse145(a):
    assert a[0];r=[1/a[0]]
    for k in range(1,70):r.append(-sum((a[i]*r[j] for i,j in CONV70[k] if i and a[i] and r[j]),F(0))/a[0])
    invsq=conv(r[:15],r[:15],CONV15);dual=[]
    for j in range(5):dual.extend(-v for v in conv(a[70+15*j:85+15*j],invsq,CONV15))
    return tuple(r)+tuple(dual)
def eval145(nodes,need,bind):
    V={}
    for i in need:
        n=nodes[i];o=n[0]
        if o=='c':v=(F(n[1],n[2]),)+(F(0),)*144
        elif o=='v':v=bind[i]
        elif o=='a':v=tuple(x+y for x,y in zip(V[n[1]],V[n[2]]))
        elif o=='n':v=tuple(-x for x in V[n[1]])
        elif o=='m':v=mul145(V[n[1]],V[n[2]])
        elif o=='i':v=inverse145(V[n[1]])
        else:raise ValueError(o)
        V[i]=v
    return V
def homogeneous(a,d):return tuple(c if len(m)==d else F(0) for c,m in zip(a,MON))
def divL(a,d):
    assert all(not c or len(m)==d for c,m in zip(a,MON))
    residual={m:c for c,m in zip(a,MON) if c};q={};lc=[F(5),F(-20),F(24),F(6)]
    while any(3 in m for m in residual):
        m=max((m for m in residual if 3 in m),key=lambda m:(m.count(3),m));co=residual[m]/6
        mm=list(m);mm.remove(3);mm=tuple(mm);q[mm]=q.get(mm,F(0))+co
        for j,c in enumerate(lc):
            target=tuple(sorted(mm+(j,)));residual[target]=residual.get(target,F(0))-co*c
            if not residual[target]:del residual[target]
    assert not residual,('exact L1 division failed',d,residual)
    quotient=tuple(q.get(m,F(0)) for m in MON);L=(F(0),*lc,*((F(0),)*65))
    assert conv(L,quotient,CONV70)==a
    return quotient
def binarymul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:out[i+j]+=x*y
    return out
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    pp=ROOT/'reports/AMBIENT_B18_CORE_CUBIC_ADJOINT_PREFLIGHT_20261104_v63/literal-lambda1-dual5-preflight.json';pr=pp.parent/'report.json'
    P=read(pp);assert read(pr)['all_passed'];cp=ROOT/P['scalar_DAG_source'];C=read(cp);nodes=C['DAG_scalar_nodes']
    A=[[F(c) for c in r] for r in P['unit98_Jacobian']]
    l0=[[F(c) for c in r] for r in P['selected5_lambda0_rows']]
    l1=[[[F(c) for c in r] for r in a] for a in P['remaining5_lambda1_rows_by_core_direction']]
    g0=P['own_reverse_linear_selected_gradients'];g1=P['own_reverse_linear_remaining_gradients']
    assert all(not F(c[0]) for row in g1 for c in row)
    assert all(sum(l0[i][k]*A[k][j] for k in range(98))==F(g0[i][j][0]) for i in range(5) for j in range(98))
    assert all(sum(l1[i][d][k]*A[k][j] for k in range(98))==F(g1[i][j][d+1]) for i in range(5) for d in range(4) for j in range(98))
    b0={int(i):tuple(map(F,v)) for i,v in P['all_four_variable_source2_bindings']}
    bd={int(i):[[F(c) for c in a] for a in rows] for i,rows in P['all_five_dual_source2_bindings']}
    assert b0.keys()==bd.keys()
    bind={i:v+(F(0),)*55+tuple(c for row in bd[i] for c in row) for i,v in b0.items()}
    assert all(len(a)==145 for a in bind.values())
    V=eval145(nodes,P['scalar_ancestor_nodes'],bind)
    equations=[V[i] for i in P['unit98_equations']]
    assert all(not any(a[:15]) for a in equations)
    assert all(not any(a[70+15*d:75+15*d]) for a in equations for d in range(5))
    E3=[homogeneous(a[:70],3) for a in equations]
    selectedraw=[V[i] for i in P['selected5_F_roots']]
    selectedcorrect=[];U2=[]
    diag=list(map(F,P['selected5_Diagonal_G1_coefficients']))
    for i,a in enumerate(selectedraw):
        assert not any(a[:15])
        r=homogeneous(a[:70],3);correct=tuple(r[j]-sum(l0[i][k]*E3[k][j] for k in range(98)) for j in range(70))
        q=divL(correct,3);selectedcorrect.append(correct);U2.append(tuple(-c/diag[i] for c in q))
    raw5=[V[i] for i in P['remaining5_F_roots']];mixed=[];correct4=[];core3=[]
    for i,a in enumerate(raw5):
        assert not any(a[:35])
        raw=homogeneous(a[:70],4);adj=[F(0)]*70
        for d in range(4):
            for k,c in enumerate(l1[i][d]):
                if not c:continue
                for j,m in enumerate(MON):
                    if len(m)==3 and E3[k][j]:adj[POS[tuple(sorted(m+(d,)))]]+=c*E3[k][j]
        mixed_i=[];acc=[F(0)]*70
        for d in range(5):
            dual=a[70+15*d:85+15*d];assert not any(dual[:5])
            q=tuple(c if len(m)==2 else F(0) for c,m in zip(dual,MON[:15]))+(F(0),)*55
            mixed_i.append(q);prod=conv(q,U2[d],CONV70)
            acc=[x+y for x,y in zip(acc,prod)]
        numerator=tuple(raw[j]-adj[j]+acc[j] for j in range(70));q=divL(numerator,4)
        mixed.append(mixed_i);correct4.append(numerator);core3.append(q)
    geo=ROOT/'reports/AMBIENT_B18_CORE_NODAL_QUARTIC_20261104_v63/fixed-quadrics-nodal-quartic.json';G=read(geo)
    phi=[[F(c) for c in row] for row in G['binary_form_coefficient_matrix4_by5']]
    binary=[]
    for a in core3:
        r=[F(0)]*13
        for c,m in zip(a,MON):
            if not c:continue
            assert len(m)==3;v=[F(1)]
            for j in m:v=binarymul(v,phi[j])
            for j,x in enumerate(v):r[j]+=c*x
        binary.append(r)
    upper=any(any(r) for r in binary)
    cubicmons=[m for m in MON if len(m)==3];cmat=[[q[POS[m]] for m in cubicmons] for q in core3]
    MM=sy.Matrix([[sy.Rational(c.numerator,c.denominator) for c in r] for r in cmat]);rank=MM.rank()
    zz=sy.symbols('x112 x117 tauA tauB');uu,ww=sy.symbols('u w')
    cexpr=[str(sy.factor(sum(sy.Rational(c.numerator,c.denominator)*sy.prod(zz[j] for j in m) for c,m in zip(a,MON) if c))) for a in core3]
    bexpr=[str(sy.factor(sum(sy.Rational(c.numerator,c.denominator)*uu**(12-j)*ww**j for j,c in enumerate(r)))) for r in binary]
    rgates=[ROOT/'reports/AMBIENT_B18_R_ZERO_SLICE_IMAGE_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in rgates)
    termcount=sum(bool(c) for a in equations+selectedcorrect+U2+correct4+core3 for c in a)+sum(bool(c) for a in binary for c in a)
    assert termcount<=150000
    proof={'four_variable_monomials_through4':MON,'dual5_core_monomials_through2':MON[:15],
      'all98_residual_jets145':equations,'all98_base_residuals_zero_through2':True,'all98_dual_residuals_zero_through_core1':True,
      'all2450_adjoint_identities_rechecked':True,'selected5_raw145':selectedraw,'selected5_corrected_cubics70':selectedcorrect,
      'selected5_eliminated_U2_quadrics70':U2,'remaining5_raw145':raw5,'remaining5_mixed_dual_quadrics':mixed,
      'remaining5_corrected_quartic_numerators70':correct4,'all5_actual_core_cubics70':core3,
      'all10_L_division_identities_all70coefficients_verified':True,
      'remaining5_physical_rows':P['remaining5_physical_rows'],'cubic_monomials':cubicmons,'actual_core_cubic_matrix5_by20':cmat,
      'actual_core_cubic_span_rank':rank,'actual_core_cubic_factored_expressions':cexpr,
      'binarydegree12_normalization_pullback_coefficients':binary,'binarydegree12_factored_expressions':bexpr,
      'some_binary_section_nonzero':upper,'actual_core_dimension_at_most1_if_nonzero_section':upper,
      'actual_R_lower_bound_source_proof':'E4 finite integral over D0 has a point above the actual R DVR prime; that prime misses Npair projection by norm lead -81600T4, so L is invertible and the prime is on S/J. Hence dim(core)>=1.',
      'actual_core_dimension_exactly1':upper,'source3_4_not_solved_or_assumed_zero':True,
      'no_H0_CM_socle_E4_D0_depth_or_full_ideal_claim':True,'expanded_nonzero_coefficient_term_count':termcount}
    p=out/'five-actual-core-cubics-and-binary12.json';p.write_text(json.dumps(enc(proof),separators=(',',':'))+'\n');assert p.stat().st_size<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',pp,pr,cp,geo,*rgates]
    report={'id':STEM,'all_passed':True,'status':'first5_actual_core_cubics_with_IFT_and_adjoint_corrections_computed',
      'cubic_span_rank':rank,'actual_core_cubic_expressions':cexpr,'nonzero_binary12_sections':[bool(any(r)) for r in binary],
      'actual_core_dimension_upper1_proved':upper,'actual_R_source_proves_core_dimension_lower1':True,
      'actual_core_dimension_exactly1':upper,'source3_or4_solve_used':False,'socle_CM_or_E4_D0_depth_proved':False,
      'expanded_nonzero_coefficient_terms':termcount,'output_bytes':p.stat().st_size,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
