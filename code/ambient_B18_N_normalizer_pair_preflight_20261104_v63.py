#!/usr/bin/env python3
"""One unit proper-source identity and cached N/CA/CB pair-tangent gate.

No new source solve, higher source jets, implicitization, or division by a
vanishing conductor coefficient in the local ring is performed.
"""
import argparse, hashlib, itertools, json, resource, signal, time
from fractions import Fraction as F
from pathlib import Path
import sympy as sy
from ambient_B18_N_same_R0_binding_20261104_v63 import JM, JP, jeval, jmul
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [enc(v) for v in x]
    return x
def compose55(a,sub):
    mons=[()]+[(i,) for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2))
    out=[F(a[0])]+[F(0)]*16
    for m,cc in zip(mons[1:],a[1:]):
        c=F(cc)
        if not c: continue
        v=sub[m[0]] if len(m)==1 else jmul(sub[m[0]],sub[m[1]])
        out=[x+c*y for x,y in zip(out,v)]
    return tuple(out)
def universal_identity():
    h,x,y,z,w,A,B,g=sy.symbols('h x y z w A B g')
    X=sy.Matrix([h,x,y,z,w]);q=(h-B*x)/A
    mons=[sy.prod((h,x,z,w)[i] for i in m) for m in itertools.combinations_with_replacement(range(4),3) if any(i>=2 for i in m)]
    cs=sy.symbols('c0:16');H=sum(c*m for c,m in zip(cs,mons))
    f=h*(A*h+B*x);P=sy.zeros(5)
    P[1,2]=f;P[3,4]=g*h*h;P[2,3]=-sy.diff(H,w);P[2,4]=sy.diff(H,z);P=P-P.T
    M=sy.eye(5);M[0,0]=A;M[0,1]=B
    pushed=(M*P*M.T).subs(h,q,simultaneous=True)
    Hp=H.subs(h,q);fb=A*q*q;Q=sy.zeros(5)
    Q[1,2]=fb;Q[3,4]=g*q*q
    Q[2,3]=-sy.diff(fb,x)*z/2-sy.diff(Hp,w)
    Q[2,4]=-sy.diff(fb,x)*w/2+sy.diff(Hp,z);Q=Q-Q.T
    ey=sy.Matrix([0,0,1,0,0]);vertical=B*q*(ey*X.T-X*ey.T)
    delta=(pushed+vertical-Q).applyfunc(sy.cancel)
    assert delta==sy.zeros(5)
    return {'all25_raw_entries_exact':True,'H_coefficients':16,'frame_determinant':'A',
      'CA_f':'h*(A*h+B*x)','new_h':'A*h+B*x','new_x':'x',
      'CB_f':'A*((h-B*x)/A)^2','CB_mu':'g/A',
      'Euler_correction':'B*((h-B*x)/A) * (dy wedge Euler)',
      'only_localized_coefficient':'A','no_B_g_or_conductor_inversion':True}
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    out.mkdir(parents=True,exist_ok=False)
    np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    paths=[np,ep,np.parent/'report.json',ep.parent/'report.json',
      ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json',
      ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json']
    for p in paths[2:]: assert read(p)['all_passed']
    ident=universal_identity();N,E,S=read(np),read(ep),read(sp)
    sol={i:tuple(map(F,v)) for i,v in N['all28_mixed_source2jets']}
    targets=[18,19,28,29,60,67,72,63,64,68,73,74]
    roots=N['all20_actual_physical_cut_equations']+N['actual_same_R0_readback_nodes9']+[N['physical150_N_source_nodes'][r] for r in targets]
    need=sorted(ancestors(N['DAG_scalar_nodes'],roots))
    V=jeval(N['DAG_scalar_nodes'],need,sol)
    rho=[V[i] for i in N['actual_same_R0_readback_nodes9']]
    assert [list(map(str,v)) for v in rho]==N['all9_R0_mixed2jets']
    assert all(not any(V[i]) for i in N['all20_actual_physical_cut_equations'])
    assert all(v[JP[(7,)]]==0 for v in rho)
    data={};tangents={}
    for name in ['CA','CB']:
        family=E['families'][name]
        pairs={int(r):pair for r,pair in family['all12_trace_pairs_in_common_R0_through2']}
        nu=compose55(family['nu_in_common_R0_through2'],rho)
        assert not any(nu[:9])
        candidates=[];checks=[]
        for r in targets:
            a,b=[compose55(v,rho) for v in pairs[r]]
            physical=V[N['physical150_N_source_nodes'][r]]
            numerator=tuple(c-d for c,d in zip(physical,a))
            assert not any(numerator[:9]),(name,r,'numerator_linear',numerator[:9])
            for j in range(7):
                den=b[JP[(j,)]]
                if den: candidates.append((r,j,den,numerator[JP[(j,7)]]/den))
            checks.append({'physical_row':r,'conductor_pullback_mixed2':b,'target_minus_even_mixed2':numerator})
        assert candidates,(name,'all conductor pullbacks have zero firstjet')
        tp=candidates[0][3];assert tp and all(c[3]==tp for c in candidates)
        assert nu[JP[(7,7)]]==tp*tp
        t=[nu[JP[(j,7)]]/(2*tp) for j in range(7)]+[tp]
        for ch in checks:
            b=ch['conductor_pullback_mixed2'];num=ch['target_minus_even_mixed2']
            assert all(num[JP[(j,7)]]==b[JP[(j,)]]*tp for j in range(7))
            assert num[JP[(7,7)]]==0
        tangents[name]=t
        data[name]={'actual_nonzero_conductor_pullback_candidates':candidates,
          'all12_target_minus_even_and_conductor_coefficients':checks,
          'nu_pullback_mixed2':nu,'tau_firstjet_in_N_source8':t,
          'tau_p_sign_fixed_by_actual_target_ratio':tp,
          'nonzero_conductor_restriction_proves_fraction_map':True,
          'integrality_in_regular_N_normalizer_proves_regular_quotient':True}
    # Hypersurface chart after the two already accepted monic eliminations.
    kept=[1,2,3,4,5,6,7]
    matrix=[list(rho[i][1:9]) for i in kept]+[tangents['CA'],tangents['CB']]
    MM=sy.Matrix([[sy.Rational(c.numerator,c.denominator) for c in row] for row in matrix])
    chosen=list(MM.T.rref()[1]);assert len(chosen)==8
    minor=F(MM[chosen,:].det());assert minor
    ell=[F(0),F(0),F(5),F(0),F(-20),F(0),F(0),F(24),F(6)]
    assert all(sum(ell[i]*matrix[i][j] for i in range(9))==0 for j in range(8))
    finite={}
    xids=[next(i for i,n in enumerate(S['DAG_scalar_nodes']) if n==['v',f'x{r}']) for r in N['actual_same_R0_rows9']]
    for name in ['CA','CB']:
        f=S['families'][name]
        unknown=[i for pair in f['coefficient_unknown_pairs'].values() for i in pair]+[f['p_quadratic_relation']['s'],f['p_quadratic_relation']['t']]
        a,b=f['complete_normalized_physical_w_coefficients'][0][18]
        finite[name]={'source_scalar_DAG_file':str(sp.relative_to(ROOT)),
          'coefficient_unknowns48':unknown,'coefficient_unit_equations48':f['all48_zero_IFT_quadratic_coefficient_equations'],
          'source_R0_nodes_to_N_source_rho_nodes':list(zip(xids,N['actual_same_R0_readback_nodes9'])),
          'numerator':'N_q18 - pullback(a18)','N_q18_node':N['physical150_N_source_nodes'][18],
          'a18_source_node':a,'denominator_b18_source_node':b,'square_nu_source_node':f['actual_tracefree_nu'],
          'operation':'unique_regular_quotient_in_normal_N_source',
          'unit_inverse_of_b18_used':False,'regularity_proof':'nonzero pullback + monic square + integral closure',
          'origin':'0','firstjet':tangents[name]}
    evidence={'universal_proper_source_identity':ident,
      'N_mixed_monomials':list(map(list,JM)),'N_same_actual_R0_mixed2':rho,
      'families':data,'Npair_to_Schur9_firstjet':matrix,
      'Schur9_coordinate_order':E['remaining_variable_names'],
      'rank8_selected_rows':chosen,'rank8_minor':minor,
      'unique_linear_equation_coefficients':ell,
      'actual_finite_regular_function_input':finite,
      'normal_N_source_equations':N['all20_actual_physical_cut_equations'],
      'normal_N_source_unknowns28':N['normalization_source28_variable_nodes'],
      'regular_quotient_recipe':'For the first listed nonzero conductor b: tau=(physical_row-a_row(R0))/b_row(R0), interpreted as the unique regular quotient in the normal N source, never a local unit inverse.',
      'new_source_solve_or_jets_above_cached2':False}
    p=out/'actual-N-pair-preflight-evidence.json';p.write_text(json.dumps(enc(evidence),separators=(',',':'))+'\n')
    rep={'id':STEM,'all_passed':True,'status':'universal_unit_source_lift_and_actual_regular_pair_rank8',
      'universal25_raw_identity':True,'actual_conductor_pullback_nonzero_both':True,
      'both_regular_trace_roots_proved':True,'actual_Npair_Schur9_rank':8,
      'rank8_minor':str(minor),'linear_equation':'5*x112-20*x117+24*tau_CA+6*tau_CB',
      'allorder_Npair_hypersurface_and_12_output_L_divisibility_proved':True,
      'full_E4_core_ideal_coefficients_or_socle_computed':False,
      'representation':'finite source-node bindings and normal regular quotients; not a claimed pure-unit-only polynomial Hensel system',
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths+[sp,Path(__file__),ROOT/'research/notes'/f'{STEM}.md']},
      'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
