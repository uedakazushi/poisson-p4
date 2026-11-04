#!/usr/bin/env python3
"""Actual source-bound first derivative interface; no new source jet solve.

Only first implicit differentiation of the frozen two48 unit systems and
constant source/group matrices.  Full nonlinear coefficient DAG extraction
is deliberately a separate representation gate.
"""
import argparse,gzip,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63'
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def em(M):return [[str(c) for c in r] for r in M.tolist()]
def mat(x):return s.Matrix([[s.Rational(c) for c in r] for r in x])
def dependencies(nodes,roots):
    todo=list(roots);done=set()
    while todo:
        i=todo.pop()
        if i in done:continue
        done.add(i)
        if nodes[i][0] not in ['c','v']:todo.extend(nodes[i][1:])
    return sorted(done)
def ad(nodes,roots,inputs):
    pos={v:i for i,v in enumerate(inputs)};v={};der={}
    for i in dependencies(nodes,roots):
        q=nodes[i];op=q[0]
        if op=='c':v[i]=F(q[1],q[2]);der[i]={}
        elif op=='v':v[i]=F(0);der[i]={pos[i]:F(1)} if i in pos else {}
        elif op=='n':v[i]=-v[q[1]];der[i]={j:-c for j,c in der[q[1]].items()}
        elif op in ['a','m']:
            a,b=q[1:];da,db=der[a],der[b]
            if op=='a':v[i]=v[a]+v[b];d={j:da.get(j,F(0))+db.get(j,F(0)) for j in da.keys()|db.keys()}
            else:v[i]=v[a]*v[b];d={j:da.get(j,F(0))*v[b]+v[a]*db.get(j,F(0)) for j in da.keys()|db.keys()}
            der[i]={j:c for j,c in d.items() if c}
        elif op=='i':
            a=q[1];assert v[a];v[i]=1/v[a];der[i]={j:-c/v[a]**2 for j,c in der[a].items()}
        else:raise ValueError((i,q))
    return v,der
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json'
    pp=ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json'
    mp=ROOT/'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63/actual-proper-source14-equations17.json'
    hp=ROOT/'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json'
    qp=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    rp=ROOT/'reports/AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63/report.json'
    op=ROOT/'reports/AMBIENT_B18_ACTUAL_INTERSECTION_ORBIT_PRODUCT_20261103_v63/certificate.json.gz'
    E,P,M,H,Q,S,R,O=map(read,[ep,pp,mp,hp,qp,sp,rp,op])
    nodes=E['DAG_scalar_nodes'];rnodes=E['R0_variable_nodes'];oldmap={int(k):v for k,v in E['original_node_to_pruned_node']}
    cols=Q['remaining_variable_names'];assert cols==['x108','x109','x112','x113','x117','x118','x119','tau_CA','tau_CB']
    base=s.zeros(9,9)
    for j in range(7):base[j+1,j]=1
    tau=[s.zeros(1,9),s.zeros(1,9)];tau[0][0,7]=1;tau[1][0,8]=1
    orbit=mat(O['physical150_by19_scalar_then_GL'])
    chosen=[tuple(x) for x in P['actual_G15_generator_positions']];orbits=mat(P['all150_GL4_orbit16'])
    allgl=[tuple(x) for x in P['GL4_generator_positions16']]
    orb15=s.Matrix.hstack(*[orbits[:,allgl.index(q)] for q in chosen])
    A=mat(P['all150_CA_fixed_flag19']);rows=P['zero_coordinate_rows36'];ai=mat(P['CA_left_inverse19']);ar=P['CA_constant_left_inverse_rows19']
    al=ai*mat([[1 if j==rows[k] else 0 for j in range(150)] for k in ar])
    fixed=P['CA_fiber10_fixed_coefficient_indices'];free=P['CA_free_fiber9_coefficient_indices']
    Jg=mat(P['source_gauge15_flag5_then_fiber10_Jacobian']);assert Jg.det()==s.Rational(-2,3)
    bad=[(1,0),(1,3),(1,4),(0,3),(0,4)]
    physical=[];flagframes=[];family_records=[]
    for family,tau_i in zip(['CA','CB'],tau):
        sys=next(x for x in E['normalizer_sources48_each'] if x['family']==family)
        unknown=sys['unknown48'];equations=sys['equations48'];inputs=unknown+rnodes
        values,der=ad(nodes,equations,inputs);assert all(values[q]==0 for q in equations)
        J=s.Matrix([[s.Rational(der[q].get(j,F(0))) for j in range(48)] for q in equations])
        B=s.Matrix([[s.Rational(der[q].get(48+j,F(0))) for j in range(9)] for q in equations])
        # Stored block determinant uses (a23,t,b23,s) variables and
        # even24/odd24 equations; current literal order is paired a,b then s,t.
        colperm=list(range(0,46,2))+[47]+list(range(1,46,2))+[46]
        rowperm=list(range(0,48,2))+list(range(1,48,2))
        Jdet=J.det()
        assert Jdet==-s.Rational(sys['unit48_determinant'])
        assert J[rowperm,colperm].det()==s.Rational(sys['unit48_determinant'])
        sol=-J.inv()*B;assert J*sol+B==s.zeros(48,9)
        Dsol=sol*base;unkrows={node:i for i,node in enumerate(unknown)}
        f=S['families'][family];names=f['zero_hensel_unknowns23_order']
        shifts={name:Dsol[unkrows[oldmap[f['coefficient_unknown_pairs'][name][0]]],:] for name in names}
        mean=Dsol[unkrows[oldmap[f['p_quadratic_relation']['s']]],:]
        pfirst=tau_i+mean/2
        outer=s.zeros(5,5);gam=s.zeros(19,9);k=s.zeros(9,9)
        for name,v in shifts.items():
            if name.startswith('gamma'):
                n=int(name[5:]);gam[n,:]=v
            else:k[int(name[1:]),:]=v
        # H has 5x5 first-derivative matrices. gamma0 is scalar, not a flag matrix.
        positions=O['GL18_generator_indices']
        Hfirst=[s.zeros(5,5) for _ in range(9)]
        for n in range(1,19):
            i,j=divmod(positions[n-1],5)
            for z in range(9):Hfirst[z][i,j]+=gam[n,z]
        coreM=mat(R['families'][family]['source_kernel_Lie_matrix'])
        for z in range(9):Hfirst[z]+=pfirst[0,z]*coreM
        physical_first=orbit*gam+mat(R['families'][family]['linear_J'])[:,19:]*k
        assert physical_first[E['common_cut_zero_rows15'],:]==s.zeros(15,9)
        assert physical_first[E['actual_R0_rows9'],:]==base
        flags=s.Matrix([[-Hfirst[z][i,j] for z in range(9)] for i,j in bad])
        physical.append(physical_first);flagframes.append((Hfirst,flags))
        family_records.append({'family':family,'implicit_unknown48':unknown,'implicit_equations48':equations,
          'actual_J48':em(J),'actual_literal_order_J48_determinant':str(Jdet),
          'block_row_permutation':rowperm,'block_column_permutation':colperm,
          'actual_input_rhs48_by9':em(B),'solution48_by9':em(sol),
          'physical150_by9_first':em(physical_first),'normalizer_flag5_by9_first':em(flags),
          'p_first':em(pfirst),'tau_first':em(tau_i),'source_frame5_by5_by9_first':[em(z) for z in Hfirst],
          'all48_linear_residuals_zero':True})
    HA,flagsA=flagframes[0]
    correction=s.Matrix.hstack(*[orbits[:,allgl.index(q)] for q in bad])*flagsA
    # Correcting the own flag frame adds its five flag matrix entries.
    fiber_before=al*(physical[0]+correction)
    constraints=flagsA.col_join(fiber_before[fixed,:])
    gauge=-Jg.inv()*constraints;assert Jg*gauge+constraints==s.zeros(15,9)
    PAfirst=physical[0]+orb15*gauge
    proper_fiber=al*PAfirst
    assert proper_fiber[fixed,:]==s.zeros(10,9)
    gmat=[s.zeros(5,5) for _ in range(9)]
    for j,(i,k) in enumerate(chosen):
        for z in range(9):gmat[z][i,k]=gauge[j,z]
    HB=flagframes[1][0]
    relativeflags=s.Matrix([[-HB[z][i,j]-gmat[z][i,j] for z in range(9)] for i,j in bad])
    actual14=proper_fiber[free,:].col_join(relativeflags)
    # M uses fiber9 followed by a,b,c,d,e.
    retain=[9,12,2,3,4,5,6,7,8]
    J9=actual14[retain,:];det9=s.factor(J9.det())
    # Compute the exact target-difference row adapter at the origin.
    ann=mat(P['CB_annihilator17_by36']);targets=[q['physical_row'] for q in E['all12_zero_physical_difference_columns']]
    R17=ann[:,[rows.index(q) for q in targets]]
    J17=mat(M['Jacobian17_by14']);piv=[0,1,10,11,13];sel=list(range(5));rest=list(range(5,17))
    JS=J17[sel,piv];assert JS.det()==s.Rational(8,3)
    L=J17[rest,piv]*JS.inv()
    norms=[s.Rational(q['constant_denominator'],q['primitive_content_unit']) for q in H['all12_actual_remaining_and_divided_records']]
    R12=s.diag(*norms)*(R17[rest,:]-L*R17[sel,:]);det12=s.factor(R12.det())
    # Corrected14 retained source coordinates pull the same quotient tangent.
    residual_first=J17*actual14
    assert residual_first==s.zeros(17,9)
    assert det9 and det12
    # Typed finite operator.  Matrix action is the same literal coefficient
    # action as in the pinned150-coordinate source builders, not an unnamed
    # pullback across an abstract quotient isomorphism.
    typed_operator={
      'coefficient_ring':{'base':'same physical R0 nine','basis':E['free4_basis'],
        'monics':['tauA^2=nuA','tauB^2=nuB'],'nu_nodes':[E['nu_CA_node'],E['nu_CB_node']]},
      'input_reference_files':{'two48_coefficients':str(ep.relative_to(ROOT)),
        'source150_core_and_groups':str(rp.relative_to(ROOT)),
        'own_source_variable_names':str(sp.relative_to(ROOT)),
        'fixed19_left_inverse_and_gauge':str(pp.relative_to(ROOT)),
        'literal_source17':str(mp.relative_to(ROOT)),
        'first5_exact_roots_and_row_units':str(hp.relative_to(ROOT))},
      'operations':[
        {'id':'parameter_pairs','op':'read_two48_pointed_roots','unknown_and_equation_records':E['normalizer_sources48_each'],
          'original_to_pruned_variable_map_reference':'two48_coefficients.original_node_to_pruned_node',
          'parameter_name_pair_map_reference':'own_source_variable_names.families[family].coefficient_unknown_pairs'},
        {'id':'p_i','op':'affine_trace_shift','formula':'tau_i + s_i/2'},
        {'id':'H_i','op':'matrix_product','formula':'G_gamma_i * core_group_i(p_i)',
          'G_gamma_rule':'I+sum_gamma_nonzero gamma_j*E_(GL18[j-1]); gamma0 is omitted here and is a separate scalar',
          'core_group_reference':'source150_core_and_groups.families[family].core_group'},
        {'id':'P_i','op':'literal150_coefficient_action','formula':'(1+gamma0_i) * Action(G_gamma_i, core150_i(p_i,k_i))',
          'core150_reference':'source150_core_and_groups.families[family].core_full150_coordinates','k2':'zero'},
        {'id':'G','op':'matrix_linear_section','formula':'I+sum_(j=0)^14 t_j*E_positions[j]',
          'positions':chosen,'unknowns15':[f't{j}' for j in range(15)]},
        {'id':'flag_frame','op':'normalized_covector_frame',
          'definition':['ell=row_x(H^-1)/row_x(H^-1)[x]',
            'k0=row_h(H^-1)-row_h(H^-1)[x]*ell','k=k0/k0[h]',
            'F rows(h,x,y,z,w)=(k,ell,e_y,e_z,e_w)'],
          'all_normalizing_units_origin':1},
        {'id':'CA_own_flag','op':'flag_frame','argument':'G*H_A'},
        {'id':'CA_own_coeff19','op':'fixed_source_left_inverse',
          'formula':'CA_left_inverse19 * selected19(Action(CA_own_flag*G,P_A))',
          'rows19':[rows[q] for q in ar],'matrix19':P['CA_left_inverse19']},
        {'id':'gauge15','op':'pointed_unit_system','unknowns':'t0..t14',
          'equations':[{'op':'flag_coefficient','frame':'CA_own_flag','coefficient':q} for q in ['a','b','c','d','e']]
             +[{'op':'fiber_coefficient_minus_center','coefficient':j,'center':P['actual_CA_source19_center'][j]} for j in fixed],
          'scalar_expansion':'t_j=t_j_even+t_j_odd*tauA; reduce tauA^2=nuA',
          'scalar_equations_order':'15 even then15 odd','scalar_unknown_order':'15 even then15 odd',
          'origin_Jacobian15':em(Jg),'origin_Jacobian15_determinant':'-2/3',
          'origin_scalar30_determinant':'4/9'},
        {'id':'actual14','op':'source_readback','fiber9':'CA_own_coeff19[free9]-center[free9]',
          'free9':free,'flag5':'normalized_covector_frame(G*H_B)[a,b,c,d,e]'},
        {'id':'kappa9','op':'select','input':'actual14','indices':retain,'origin_Jacobian':em(J9),'origin_determinant':str(det9)},
        {'id':'R17','op':'coefficient_action_matrix',
          'formula':'Ann_CB * restrict36 * Action(flag_frame(G*H_B)*G) * reconstruct150 * inject_F12',
          'Ann_CB':P['CB_annihilator17_by36'],'coordinate_rows36':rows,'F12_physical_rows':targets,
          'reconstruct150_reference':'two48_coefficients.four_divergence_relations40/eliminated_four_coordinate_reconstruction',
          'origin':em(R17)},
        {'id':'root5','op':'evaluate_explicit_rational_roots','input':'kappa9','roots':H['five_explicit_rational_roots'],
          'root_variable_order':H['remaining9_variable_order'],'pivot_order':['v0','v1','b','c','e']},
        {'id':'A5_B12','op':'finite_monomial_telescoping',
          'input_polynomials_reference':'literal_source17.equations14variable_sparse',
          'source_variables':[f'v{i}' for i in range(9)]+['a','b','c','d','e'],
          'actual_variables':'actual14','pivot_variables':['v0','v1','b','c','e'],
          'pivot_basepoint':'root5(kappa9)',
          'finite_quotient_rule':'(X^e-Y^e)/(X-Y)=sum_(j=0)^(e-1) X^(e-1-j)Y^j',
          'A5_rows':list(range(5)),'B12_rows':list(range(5,17)),
          'A5_inverse':'adj(A5)/det(A5); det origin8/3','origin_A5':em(JS)},
        {'id':'R12','op':'matrix_schur_and_unit_scaling',
          'formula':'diag(N_i)*(R17_rest-B12*A5^-1*R17_selected)',
          'N_i':[{'QQ_factor':str(q),'D_power':rec['D_denominator_power']} for q,rec in zip(norms,H['all12_actual_remaining_and_divided_records'])],
          'D':'1-kappa_d^2','identity':'dH12(kappa9)=R12*F12',
          'origin_matrix':em(R12),'origin_determinant':str(det12),
          'inverse':'adj(R12)/det(R12); known origin unit'}],
      'operator_has_no_unspecified_coefficients_or_nonunit_quotients':True,
      'operator_is_matrix_typed_not_expanded_scalar_DAG':True,
      'scalar_four_R0_component_extraction_performed':False}
    evidence={'physical_regular9_order':cols,'physical_monic_pivots':['x27','x61'],
      'physical_monic2_minor':'-20/9','physical_R0_first9_by9':em(base),
      'families':family_records,'actual_group15_positions':chosen,'group15_unit_Jacobian':em(Jg),
      'group15_first9':em(gauge),'source14_order':[f'v{i}' for i in range(9)]+['a','b','c','d','e'],
      'source14_first_by_physical9':em(actual14),'retained_source9_indices':retain,
      'retained_source9_order':['a','d','v2','v3','v4','v5','v6','v7','v8'],
      'actual_source9_by_physical9_Jacobian':em(J9),'actual_source9_unit_determinant':str(det9),
      'source17_by_physicalF12_origin':em(R17),'source_first5_pivot_Jacobian':em(JS),
      'source_remaining12_subtract_first5_origin':em(L),
      'source_dH12_by_physicalF12_origin':em(R12),'source_dH12_row_unit_determinant':str(det12),
      'all17_source_equations_first_residual_zero':True,
      'no_source_jet_or_nonunit_division_used':True,
      'source_bound_finite_matrix_typed_operator':typed_operator,
      'typed_nonlinear_recipe':{
        'pure_coefficients':'the frozen96 equations/unknown roots in families, no extra scalar solve',
        'trace_source_p':'p_i=tau_i+s_i(R0)/2; all gamma/k source functions are their own a+b*p_i',
        'source_matrix':'H_i=G_gamma_i*core_group_i(p_i), outer scalar1+gamma0_i',
        'flag_extraction':'ell=row_x(H_i^-1)/row_x(H_i^-1)[x]; k=(row_h-row_h[x]*ell)/coefficient_h',
        'group15':'G(t)=I+sum_j t_j E_chosen_j; flags_A(G(t)H_A)=0 and own fixed10 fiber readbacks=center10',
        'group30_scalar_unit_determinant':'(-2/3)^2; rank2 CA trace expansion of the15 equations',
        'source14_readbacks':'CAfree9 coefficients after own group15 gauge; five normalized flags of G(t)H_B',
        'source9_map':'retain indices9,12,2,3,4,5,6,7,8 of actual14 readbacks',
        'source17_identity':'Ann_CB * Action(F_CB*G)(P_A-P_B); physical differences outsideF12 are zero',
        'row12_identity':'subtract the exact finite first5-pivot Hadamard multiples; scale by known rational/D units to obtain dH12',
        'all_unit_inverse_origins_pinned':{'source96':['84934656/625','331776/625'],'flag_x_and_h_normalizers':1,'group15':'-2/3','first5_source':'8/3','two_monics':'-20/9','source9':str(det9),'row12':str(det12)},
        'scalar_DAG_of_nonlinear_group30_and_row12_serialized':False,
        'coefficient_extraction_still_requires_literal_nonlinear_root_and_row_nodes':True}}
    p=out/'actual-source-physical-firstjet-interface.json';p.write_text(json.dumps(evidence,indent=2)+'\n');assert p.stat().st_size<32*1024*1024
    gates=[ROOT/'reports'/f'{x}_20261104_v63_REVIEW/report.json' for x in ['AMBIENT_B18_PROPER_POLYNOMIAL_MODEL','AMBIENT_B18_PROPER_RATIONAL_ELIMINATION','AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP','AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL']]
    assert all(read(q)['all_passed'] for q in gates)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',ep,pp,mp,hp,qp,sp,rp,op,*gates]
    report={'id':STEM,'all_passed':True,'status':'actual96source_first_differentiation_group15_source9_and_relation12_unit_interfaces',
      'actual_source9_unit_minor':str(det9),'actual_relation12_unit_minor':str(det12),
      'physical_source96_unit_jacobians_checked':True,'first_source14_and_physical12_identity_bound':True,
      'nonlinear_unit_operator_recipe_proved':True,'literal_nonlinear_coefficient_DAG_or128_R0_syzygies_claimed':False,
      'D0_K_depth_or_torsion_claimed':False,'new_source_jets_or_GB':False,
      'input_sha256':{str(q.relative_to(ROOT)):sha(q) for q in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','actual_source9_unit_minor','actual_relation12_unit_minor','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
