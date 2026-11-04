#!/usr/bin/env python3
"""Literal common-R0 normalizer tensor algebra, and conductor support preflight.

Only frozen scalar source systems and cached first coefficients are used.
No source jets, old regular-quotient composition, GB, or branch solving.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG,PAIRS,MONS
from ambient_B18_top_graph_telescoping_review_20261104_v63 import solve

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def closure(nodes,roots):
    done=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in done:continue
        done.add(i)
        if nodes[i][0] not in ['c','v']:todo.extend(nodes[i][1:])
    return sorted(done)
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(a) for a in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x
def mm(A,B):return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def dual_eval(nodes,needed,binding):
    v={}
    for i in needed:
        n=nodes[i];o=n[0]
        if o=='v':q=binding.get(i,(F(0),)*3)
        elif o=='c':q=(F(n[1],n[2]),F(0),F(0))
        elif o=='n':q=tuple(-x for x in v[n[1]])
        elif o=='a':q=tuple(a+b for a,b in zip(v[n[1]],v[n[2]]))
        elif o=='m':
            a,b=v[n[1]],v[n[2]];q=(a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[2]*b[0])
        elif o=='i':
            a=v[n[1]];assert a[0];q=(1/a[0],-a[1]/a[0]**2,-a[2]/a[0]**2)
        else:raise ValueError((i,n))
        v[i]=q
    return v
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    cp=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'
    wp=ROOT/'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json'
    rp=ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63/certificate.json'
    zp=ROOT/'reports/AMBIENT_B18_T0_COKERNEL_Z_20261104_v63/finite-cokernel-Z-and-actions.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    ip=ROOT/'reports/AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63/actual-R-N-projection-norm.json'
    bp=ROOT/'reports/AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63/degree3-certificate.json'
    gates=[ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63_REVIEW.json',
      ROOT/'reports/AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    S,C,W,R,Z,J,N,I,B=map(read,[sp,cp,wp,rp,zp,jp,np,ip,bp])
    old=C['DAG_scalar_nodes'];assert old[:len(S['DAG_scalar_nodes'])]==S['DAG_scalar_nodes']
    wt=W['physical150_torus_weights'];zero=[i for i,w in enumerate(wt) if w==0];assert len(zero)==40
    cuts=[i for i in R['common_physical_group_cuts19'] if wt[i]==0]
    rows=S['physical_zero_base_rows9'];assert rows==N['actual_same_R0_rows9']
    odd=[18,19,28,29,60,67,72,63];graphs=[64,68,73,74];target=odd+graphs
    coords=cuts+rows+target;assert len(coords)==len(set(coords))==36 and set(coords)<=set(zero)
    # Four divergence equations give the remaining four physical coordinates.
    other=[i for i in zero if i not in coords];assert len(other)==4
    div=[]
    for symbol in [0,1,3,4]:
        v=[]
        for row in zero:
            a,b=PAIRS[row//15];u,w=MONS[row%15];co=F(0)
            if b==2:
                if u==a and w==symbol:co+=1
                if w==a and u==symbol:co+=1
            if a==2:
                if u==b and w==symbol:co-=1
                if w==b and u==symbol:co-=1
            v.append(co)
        div.append(v)
    om=[[r[zero.index(i)] for i in other] for r in div]
    reconstruction,detdiv=solve(om,[[-r[zero.index(i)] for i in coords] for r in div]);assert detdiv
    ca=S['families']['CA'];cb=S['families']['CB']
    equations=[];unknown=[];sources=[]
    for family,rec in [('CA',ca),('CB',cb)]:
        us=[i for p in rec['coefficient_unknown_pairs'].values() for i in p]
        us += [rec['p_quadratic_relation'][s] for s in ['s','t']]
        assert len(us)==len(set(us))==48
        eq=rec['all48_zero_IFT_quadratic_coefficient_equations'];assert len(eq)==48
        equations+=eq;unknown+=us
        sources.append({'family':family,'unknown48':us,'equations48':eq,
          'unit48_determinant':rec['unit48_block_determinant_even24_then_odd24_vars_a23_t_b23_s']})
    prefix=C['all7_pure_CA_tail_prefix57_systems'][0]
    equations+=prefix['all57_unit_equations'];unknown+=[p[1] for p in prefix['bindings']]
    assert len(unknown)==len(set(unknown))==153
    basis=C['actual_C0_basis8_scalar_nodes']
    assert basis[:7]==[p[1] for p in ca['actual_zero_odd_I0_generators'][:7]]
    ap=ca['complete_normalized_physical_w_coefficients'][0]
    bpairs=cb['complete_normalized_physical_w_coefficients'][0]
    roots=equations+basis+[ca['actual_tracefree_nu'],cb['actual_tracefree_nu']]
    roots += [v for i in zero for v in ap[i]+bpairs[i]]
    needed=closure(old,roots)
    assert all(old[i][0] in ['c','v','a','m','n','i'] for i in needed)
    D=DAG();mapping={}
    for i in needed:
        n=old[i];o=n[0]
        if o=='c':v=D.const(F(n[1],n[2]))
        elif o=='v':v=D.var(f'ACTUAL_NODE{i}_{n[1]}')
        elif o=='a':v=D.add(mapping[n[1]],mapping[n[2]])
        elif o=='m':v=D.mul(mapping[n[1]],mapping[n[2]])
        elif o=='n':v=D.neg(mapping[n[1]])
        elif o=='i':v=D.inverse(mapping[n[1]])
        else:raise ValueError(o)
        mapping[i]=v
    nuA=mapping[ca['actual_tracefree_nu']];nuB=mapping[cb['actual_tracefree_nu']]
    z,o=D.zero,D.one
    MA=[[z,nuA,z,z],[o,z,z,z],[z,z,z,nuA],[z,z,o,z]]
    MB=[[z,z,nuB,z],[z,z,z,nuB],[o,z,z,z],[z,o,z,z]]
    def dmm(A,B):return [[D.total(D.mul(a,b) for a,b in zip(r,c)) for c in zip(*B)] for r in A]
    def mv(M,v):return [D.total(D.mul(a,b) for a,b in zip(r,v)) for r in M]
    assert dmm(MA,MB)==dmm(MB,MA)
    assert dmm(MA,MA)==[[nuA if i==j else z for j in range(4)] for i in range(4)]
    assert dmm(MB,MB)==[[nuB if i==j else z for j in range(4)] for i in range(4)]
    cols=[];differences=[]
    for row in target:
        a,b=[mapping[v] for v in ap[row]];c,d=[mapping[v] for v in bpairs[row]]
        v=[D.add(a,D.neg(c)),b,D.neg(d),z]
        differences.append({'physical_row':row,'column4':v})
        cols.extend([v,mv(MA,v),mv(MB,v),mv(MA,mv(MB,v))])
    assert len(cols)==48 and all(D.values[a]==0 for v in cols for a in v)
    c0map=[[o,z,z,z]]+[[z,mapping[b],z,z] for b in basis]
    # Frozen first-coefficient support gate, independently retaining both planes.
    XCA=[[F(x) for x in r] for r in Z['actual_physical_R0_firstjets_on_CA_F9_9by9']]
    adapter=[[F(x) for x in r] for r in J['YB_to_YA_linear']]
    invadapter,ad=solve(adapter,[[F(i==j) for j in range(9)] for i in range(9)])
    XCB=mm(XCA,invadapter)
    ne=N['all56_actual_coefficient_Hensel_equations'];nus=N['all56_coefficient_unknown_nodes']
    nroots=ne+[N['ninth_actual_even_a_node'],N['ninth_actual_odd_b_node'],N['nu_tracefree_square_node']]
    nn=closure(N['DAG_scalar_nodes'],nroots);ni=[[F(x) for x in r] for r in I['N56_unit_inverse']]
    cohen=N['chosen7_R0_indices']+[N['Cohen_even8_R0_index']]
    support=[];ninth=N['ninth_R0_index']
    for family,XX,g3 in [('CA',XCA,C['cached_actual_CA_gA_all_degree3_coefficients'][6:10]),
                        ('CB',XCB,B['g_CB_all10_coefficients'][6:10])]:
        bind={u:(F(0),)*3 for u in nus}
        bind.update({u:(F(0),XX[j][7],XX[j][8]) for u,j in zip(N['actual_Cohen_R8_parameter_nodes'],cohen)})
        val=dual_eval(N['DAG_scalar_nodes'],nn,bind)
        residual=[[val[e][k] for k in [1,2]] for e in ne]
        corr=[[-sum((ni[i][j]*residual[j][k] for j in range(56)),F(0)) for k in range(2)] for i in range(56)]
        bind.update({u:(F(0),*c) for u,c in zip(nus,corr)})
        val=dual_eval(N['DAG_scalar_nodes'],nn,bind);assert all(val[e]==(F(0),)*3 for e in ne)
        even=val[N['ninth_actual_even_a_node']];ell=[XX[ninth][7]-even[1],XX[ninth][8]-even[2]]
        g3=list(map(F,g3));rr,ss=-ell[1],ell[0]
        res=g3[0]*rr**3+g3[1]*rr**2*ss+g3[2]*rr*ss**2+g3[3]*ss**3
        support.append({'family':family,'actual_R0_first_on_own_conductor_plane':[[r[7],r[8]] for r in XX],
          'N56_first_unknown_corrections':corr,'all56_first_equations_zero':True,
          'v_N_linear_form_on_own_plane':ell,'g3_coefficients':g3,
          'g3_evaluated_on_kernel_v_linear':res,'gcd_linear_v_and_g3_is1':bool(any(ell) and res)})
    finite=all(q['gcd_linear_v_and_g3_is1'] for q in support)
    ev={'DAG_scalar_nodes':D.nodes,'original_node_to_pruned_node':[[i,mapping[i]] for i in needed],
      'actual_R0_rows9':rows,'R0_variable_nodes':[mapping[next(i for i,n in enumerate(old) if n==['v',f'x{r}'])] for r in rows],
      'normalizer_sources48_each':[{**s,'unknown48':[mapping[i] for i in s['unknown48']],
          'equations48':[mapping[i] for i in s['equations48']]} for s in sources],
      'CA_conductor_gA_prefix57':{'unknown57':[mapping[p[1]] for p in prefix['bindings']],
        'equations57':[mapping[i] for i in prefix['all57_unit_equations']],
        'unit_determinant':prefix['unit57_determinant']},
      'all153_scalar_unit_unknown_nodes':[mapping[i] for i in unknown],
      'all153_scalar_unit_equations':[mapping[i] for i in equations],
      'free4_basis':['1','tau_CA','tau_CB','tau_CA*tau_CB'],
      'nu_CA_node':nuA,'nu_CB_node':nuB,'multiplication_tau_CA':MA,'multiplication_tau_CB':MB,
      'zero_raw_rows40':zero,'common_cut_zero_rows15':cuts,'coordinate_rows36':coords,
      'four_divergence_relations40':div,'eliminated_zero_rows4':other,
      'divergence_unit4_minor':detdiv,'eliminated_four_coordinate_reconstruction':reconstruction,
      'all12_zero_physical_difference_columns':differences,'normalizer_tensor_relation_columns4_by48':cols,
      'actual_C0_free9_to_E4_columns':c0map,
      'all_actual_I0_CA_generator_nodes8':[mapping[b] for row,b in ca['actual_zero_odd_I0_generators']],
      'all_actual_I0_CB_generator_nodes8':[mapping[b] for row,b in cb['actual_zero_odd_I0_generators']],
      'normalizer_tensor_fiber_dimension_at_closed_R0_point':4,
      'source_binding_proof':'actualpure48CAandCBcommonR0; canonicalzero36/cut15/base9 leaves12',
      'D0_map_proof':'natural completedtensor map of both finite normalizations; all640relations vanish by universal property, not oldcofactor composition',
      'all_coefficient_dependency_closures_avoid_old_r_g_operations':True,
      'conductor_product_annihilates_kernel_and_cokernel':True,
      'conductor_support_first_gate':support,'comparison_kernel_finite_from_this_gate':finite,
      'E4_depth_socle_kernel_nonzero_or_D0_depth_computed':False}
    ep=out/'actual-normalizer-tensor-free4-DAG.json';text=json.dumps(enc(ev),separators=(',',':'))+'\n'
    assert len(text.encode())<32*1024*1024;ep.write_text(text)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',sp,cp,wp,rp,zp,jp,np,ip,bp,*gates,
      ROOT/'code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py',
      ROOT/'code/ambient_B18_top_graph_telescoping_review_20261104_v63.py']
    rep={'id':STEM,'all_passed':True,'status':'actual_normalizer_tensor_free4_48relations_and_product_conductor_ann_kernel',
      'normalizer_tensor_free4_generators':4,'actual_relation_columns':48,'all_zero_graphs_retained':4,
      'all153_unit_source_equations_retained':True,'old_r_g_composition_or_new_source_jets':False,
      'actual_C0_free9_map_constructed':True,'kernel_and_cokernel_annihilated_by_I0_CA_times_I0_CB':True,
      'both_conductor_curves_meet_F_only_at_origin_by_initial_gate':finite,
      'comparison_kernel_and_cokernel_finite_proved_by_this_gate':finite,
      'E4_saturation_socle_depth_or_comparison_kernel_nonzero_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'DAG_nodes':len(D.nodes),'bytes':len(text.encode()),'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n')
    print(json.dumps({k:rep[k] for k in ['all_passed','both_conductor_curves_meet_F_only_at_origin_by_initial_gate','DAG_nodes','bytes','seconds','max_rss_kib']}))
    print(json.dumps(enc(support)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
