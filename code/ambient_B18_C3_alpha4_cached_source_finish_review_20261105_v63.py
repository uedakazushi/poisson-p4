#!/usr/bin/env python3
"""Finish the independent alpha4 SOURCE audit from its frozen bound checkpoint."""
import ast,hashlib,importlib.util,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_C3_ALPHA4_CACHED_SOURCE_20261105_v63'
ID=BASE+'_REVIEW'
SOURCE_PIN='3c4909ff3908664e2986679aabeb99c592ecb7debc107b8faed4dbbdb9e86082'
CORE_SHA='430cc27fc33b263cdcaec0fd616b87e88824bef702a3ff49fbc7136eaf428313'
CHECKPOINT_SHA='3146169fb3eab4007554271fa5e116408150fc821c2299df38cb0b41055639c0'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def main():
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;assert out.exists() and not(out/'report.json').exists()
    corepath=ROOT/'code/ambient_B18_C3_alpha4_cached_source_review_20261105_v63.py';assert sha(corepath)==CORE_SHA
    spec=importlib.util.spec_from_file_location('own_source_core',corepath);own=importlib.util.module_from_spec(spec);spec.loader.exec_module(own)
    pins={};totalbytes=0
    def pin(p,h=None):
        nonlocal totalbytes
        p=pathlib.Path(p);s=sha(p)
        if h is not None:assert s==h,(str(p),s,h)
        key=str(p.relative_to(ROOT))
        if key not in pins:totalbytes+=p.stat().st_size
        pins[key]=s;return read(p) if p.suffix=='.json' else None
    pin(corepath,CORE_SHA)
    checkpoint=pin(out/'source-majorant-and-frame-checkpoint.json',CHECKPOINT_SHA)
    assert checkpoint['whole_operational_prebound']==61 and checkpoint['no_positive_order_requested_root_or_alpha_value_computed']
    pr=pin(ROOT/'reports'/BASE/'report.json',SOURCE_PIN);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/BASE/p,h)
    c=read(ROOT/'reports'/BASE/'compiled-eight-even-roots-and-provided-joint105-plane3-cache.json')
    profile=read(ROOT/'reports'/BASE/'alpha4-source-callback-binding-and-operational-height-profile.json')
    supplement=pin(ROOT/'reports'/BASE/'operational-height-scope-companion.json',
        '996f041e53427dd1387a439892f553ca99200855ffd0cf32b861f1749587141f')
    pin(ROOT/'code/ambient_B18_C3_alpha4_operational_scope_20261105_v63.py',supplement['code_sha256'])
    assert supplement['all_temporary_scope_bound_including_original_DAG']==79
    assert checkpoint['input_conversion_operation_prebound']<=supplement['input_conversion_operation_bound']
    assert checkpoint['inverse_diagnostic_prebound']<=max(v[3] for v in supplement['inverse_diagnostic_bounds'])
    assert checkpoint['final_delta_subtraction_prebound']<=max(v[2] for v in supplement['delta_subtraction_bounds'])
    peer=pin(ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
        '2c2c157f217212c14c7cc46a94f252701d9b6909e1c3f13fdc5b470351c9ecd1');assert peer['all_passed']
    jointpeer=read(ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json')
    assert jointpeer['independent_sparse_dict_jets_all105_equations_and_unknowns_through3']
    assert jointpeer['literal_finite105_DAG_and_every_binding_reconstructed_from_actual_source']
    jp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    assert own.field(jp,'CA_F9_targets')==['0']*7+['r','s']
    assert Q(own.field(jp,'det_joint105'))==Q(1202315964973056,390625)
    op=ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json'
    assert own.field(op,'CA_target_positions')==list(range(48,57))
    assert own.field(op,'actual_target_pure_ell')=='z7=0,rA=0,sA=ell/5'
    # All40 expected lower values are old, independently accepted data.
    old=own.field(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json','all_requested164_root_values')
    old={i:dict(v) for i,v in old}
    basis=own.field(ROOT/'reports/AMBIENT_B18_ELL2_RAW33_QUOTIENT_20261105_v63/actual-ell2-product33-and-raw45-actions.json','actual33_basis')
    idx=[basis.index(m) for m in [[0,0,0],[0,1,0],[0,0,1],[0,1,1]]]
    ell=dict(own.field(ROOT/'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json','all164_actual_pure_ell_root_values3'))
    roots=[r[k] for r in c['root_pairs_physical4'] for k in ['aA_node','aB_node']]
    expected=[[i,[old[i].get(j,'0') for j in idx]+[ell[i][2]]] for i in roots]
    assert expected==c['expected8_lower_B6_coefficients_order_1_r_ell_rell_ell2']
    lower=dict(expected)
    for r in c['root_pairs_physical4']:
        assert list(map(Q,lower[r['aA_node']]))==list(map(Q,lower[r['aB_node']]))
    complete=read(ROOT/'reports/AMBIENT_B18_C3_COMPLETE_CUBIC_LINE_SOURCE_20261105_v63/complete155-source-cubic-operator-with-four-formal-alphas.json')
    callbacks=complete['remaining_four_alpha_source_functionals'];assert len(callbacks)==4
    for callback,meta in zip(callbacks,c['root_pairs_physical4']):
        assert callback['literal_roots']==meta and callback['value']=='UNAVAILABLE_FORMAL_NOT_ZERO'
    # Exact coefficient extraction of r*ell2 under s=r+ell/5.
    extraction={}
    for nr,ns in c['provided_monomials10']:
        value=Q(0)
        for k in range(ns+1):
            if (nr+ns-k,k)==(1,2):
                import math
                value+=Q(math.comb(ns,k),5**k)
        if value:extraction[(nr,ns)]=value
    assert extraction=={(1,2):Q(1,25),(0,3):Q(3,25)}
    custody=pin(ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUE_CUSTODY_20261105_v63/report.json',
        '5318a5b3e08ed15b5dd5e55a8350682ff53979b01acfff465ac1b3b2a05501cf')
    for p,h in custody['input_sha256'].items():pin(ROOT/p,h)
    assert custody['eight_root_or_alpha_values_never_evaluated']
    assert checkpoint['all45_constants_exact_numerator_denominator']==custody['literal_numerator_denominator_constant_fixtures']
    assert sum(v[2]!=1 for v in custody['literal_numerator_denominator_constant_fixtures'])==32
    path=ROOT/'code/ambient_B18_C3_alpha4_cached_values_20261105_v63.py';code=path.read_text();tree=ast.parse(code)
    funcs={n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)}
    # Literal order is essential: BEFORE each potentially growing Fraction
    # operation. Zero shortcuts and negation have no growth.
    def statements(name):return funcs[name].body
    assert isinstance(statements('times')[-2],ast.Expr)
    assert ast.unparse(statements('times')[-2].value).startswith('checked_height(')
    assert ast.unparse(statements('times')[-1])=='return note(a * b)'
    assert ast.unparse(statements('add')[-2].value).startswith('checked_height(')
    assert ast.unparse(statements('add')[-1])=='return note(a + b)'
    assert ast.unparse(statements('parse')[-2].value).startswith('checked_height(')
    assert ast.unparse(statements('parse')[-1])=='return note(Q(s))'
    assert ast.unparse(statements('neg')[-1])=='return note(-a)'
    for name in ['conv','inverse']:
        calls={n.func.id for n in ast.walk(funcs[name]) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)}
        assert calls<=({'Q','enumerate','add','times'} if name=='conv' else {'Q','range','enumerate','add','neg','times'})
        assert not any(isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Mult,ast.Div,ast.Pow))
                       and not(isinstance(n.op,ast.Mult) and isinstance(n.left,ast.List))
                       for n in ast.walk(funcs[name]))
    assert "elif op=='c':v=[parse(str(r[2])+'/'+str(r[3]))]+[Q(0)]*5" in code
    for literal in [
       'signal.alarm(60)','resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))',
       'assert owned_bound<=70000<100000','assert raw_inputs+150000<=1000000',
       'assert peak==184','ROOT_ONLY_STOP_CHECKPOINT.json','assert MAX_PRECHECK<=profile',
       "assert conv(v,cache[r[2]])==one==conv(cache[r[2]],v)",
       'assert delta[:5]==[Q(0)]*5','assert overlap_count==630',
       "out.mkdir(exist_ok=False)","parser.add_argument('--source-peer-report',required=True)"]:
        assert literal in code,literal
    assert not (ROOT/'reports/AMBIENT_B18_C3_ALPHA4_CACHED_VALUES_20261105_v63').exists()
    # Reconstruct the worker's full owned-cell formula before any values.
    own_report_reserve=2000
    def bools(x):
        if isinstance(x,bool):return 1
        if isinstance(x,list):return sum(map(bools,x))
        if isinstance(x,dict):return sum(map(bools,x.values()))
        return 0
    # Its conservative int counter includes bool control slots too.
    future_inputs=sum(own.nc(x)+bools(x) for x in [pr,c,profile])+own_report_reserve
    worker_upper=future_inputs+2940+735+7945+1116+20000
    assert worker_upper==custody['whole_future_owned_including2000peerreserve']==59866
    assert worker_upper<70000<=75000<100000
    assert totalbytes+200000<32<<20
    # An unrelated nilpotent unit fixture, not any source coefficient.
    M=own.MON;MI=own.MI
    def prod(a,b):
        d={}
        for m,c0 in a.items():
            for n,c1 in b.items():
                k=tuple(i+j for i,j in zip(m,n))
                if k in MI:d[k]=d.get(k,Q(0))+c0*c1
        return {m:c for m,c in d.items() if c}
    def plus(a,b,s=1):
        d=dict(a)
        for m,c in b.items():d[m]=d.get(m,Q(0))+s*c
        return {m:c for m,c in d.items() if c}
    u={(1,0):Q(2,3),(0,1):Q(-5,7),(1,1):Q(4,9),(0,2):Q(-1,11),(1,2):Q(3,13)}
    inv=plus(plus({(0,0):Q(1)},u,-1),prod(u,u));inv=plus(inv,prod(prod(u,u),u),-1)
    unit=plus({(0,0):Q(1)},u)
    assert prod(unit,inv)==prod(inv,unit)=={(0,0):Q(1)}
    overflow_limbs=(1<<2047).bit_length()+Q(4).numerator.bit_length()
    assert overflow_limbs==2051 and overflow_limbs>2048
    assert ast.unparse(funcs['checked_height'].body[-1])=='assert n <= 2048'
    proof={'all_eight_original_root_ids_and_complete_four_callback_bindings_verified':True,
      'all7945_original_DAG_records_and_minimal_ancestor_closure_verified':True,
      'same105_original_unknown_and_equation_ids_and_actual_r_s_targets':True,
      'actual_joint_inverse_Jacobian_is_QQ_unit':'1202315964973056/390625',
      'all630_input_lower_overlaps_rebuilt_without_requested_readbacks':True,
      'B6_exact_ring_homomorphism':'s=r+ell/5, then quotient by r2,ell3; unit inverses commute with this map',
      'alpha_input_coefficient_formula':'(delta_(r*s2)+3delta_(s3))/25; no delta is paired in this SOURCE audit',
      'all40_expected_lower_root_values_rebound_to_C2_and_pureell_peers':True,
      'all45_n_d_constants_and32_nonintegral_fixtures_verified':True,
      'complete_per_node_absolute_denominator_majorant_hash_rebuilt':True,
      'lastuse_peak':184,'coefficient_height_bound':39,'published_all_temporary_bound':79,
      'own_tighter_before_operation_bounds':{'DAG':61,'input_conversion':40,'inverse_checks':35,'delta_subtraction':22,'whole':61},
      'strict_worker_prechecks_before_add_mul_parse_and_every_convolution_verified':True,
      'both_inverse_product_and_all40_prior_overlap_checks_in_worker':True,
      'unrelated_B6_geometric_inverse_and_before_product_stop_fixtures':True,
      'future_whole_numeric_with2000_peer_reserve':59866,'static_cap':70000,'parent_static_cap':75000,'hard_cap':100000,
      'source_review_whole_static_numeric_bound':125000,'source_hard_cap':150000,
      'fresh_RAW_hard_bytes':32<<20,'worker_RAW_input_output_cap':1000000,'hard_bits':2048,
      'hard_seconds':60,'hard_address_bytes':2<<30,
      'runtime_stops_preserve_only_completed_roots_no_cap_raise':True,
      'all_four_alpha_and_all_positive_order_requested_roots_still_uncomputed':True,
      'no_source_solve_gauge_new_physical_IFT_or_jet':True,
      'conditional_ONE_existing_cache_worker_ready':True,
      'scope':'Any nonzero alpha gives exact HF3=3 and upper19; allzero gives exact HF3=4 and lower19. Neither determines full C3 length or all torsion exponents.'}
    ep=out/'independent-source-frame-lifetime-code-and-ONE-worker-readiness.json';save(ep,proof)
    result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,
      'input_sha256':pins,'core_code_sha256':CORE_SHA,'finish_code_sha256':sha(__file__),
      'evidence_sha256':{'source-majorant-and-frame-checkpoint.json':CHECKPOINT_SHA,ep.name:sha(ep)},
      'prepared_worker_sha256':'3067f4f71e9e2eaa097680e696a3c91d4f36ac9d7187ebaeebfc8a7f9c90e048',
      'conditional_ONE_cache_emission_GO_ready':True,'no_alpha_values_evaluated':True,'elapsed_finish_seconds':time.monotonic()-started}
    save(out/'report.json',result);save(out/'phase.json',{'status':'FINAL_SOURCE_PASS','value_emission_not_done_by_reviewer':True})
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}));signal.alarm(0)
if __name__=='__main__':main()
