#!/usr/bin/env python3
"""RIGHT2 post-root bound and pending prepared-worker independent audit."""
import argparse,ast,hashlib,importlib.util,json,math,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_C3_RIGHT2_LAMBDA3_CACHE_SOURCE_20261105_v63'
ID=BASE+'_REVIEW'
SOURCE='ae2b9c8162f53a396dcfb52288edbb53035569742b384de9b34286cf8cdc0d53'
CORE='f24ba57a8fa691231d150d05f9783569b3bdd221d5707ada02b4f495a3ea1717'
CHECK='32fc381eacf048dc3e97b2ac7124e4133aa09a0484a5fea5b48fbe77318dae26'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,d):pathlib.Path(p).write_text(json.dumps(d,indent=2)+'\n')
def ceilabs(q):return (abs(q.numerator)+q.denominator-1)//q.denominator
def post_bound_check():
    path=ROOT/'reports/AMBIENT_B18_C3_RIGHT2_WHOLE_VALUE_GATE_20261105_v63/report.json'
    assert sha(path)=='9fc7c8fffd965ebd5fc0cbaac7399d203e960b671fb69707c802849295f1a0bb'
    gate=read(path)
    op=read(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json')
    c2=read(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    c2={r['physical_row']:dict(r['delta']) for r in c2['physical12_coefficients']}
    L={r:(U,D) for r,U,D in gate['lambda3_positive_coefficient_majorants']}
    assert L=={18:(1680,3456000),60:(726,1728000),67:(1094,6912000)}
    def bound(terms):
        U=0;D=1;parse=0
        for (kind,row,index),co in terms:
            c=Q(co)
            if kind=='lambda':u,d=L[row]
            elif kind=='cache':q=Q(c2[row].get(index,'0'));u,d=ceilabs(q),q.denominator
            elif kind=='constant':u,d=1,1
            else:raise AssertionError(kind)
            U+=ceilabs(c)*u;D=math.lcm(D,c.denominator*d)
            parse=max(parse,abs(c.numerator).bit_length(),c.denominator.bit_length())
        return {'U':str(U),'D':str(D),'height':max((U*D).bit_length(),D.bit_length()),'coefficient_parse_bound':parse,'term_count':len(terms)}
    bounds=[[bound(terms) for terms in row] for row in op['matrix2_rows_A4_B4_columns_HB0_HB1_linear_leaf_operators']]
    assert bounds==gate['four_affine_entry_majorants']
    H=max(x['height'] for row in bounds for x in row);assert H==56
    assert 12*H+3==gate['all_postroot_affine_minor_inverse_rank1_operation_bound']==675<2048
    # Independently reconstruct the two affine factors without reading any
    # prospective lambda, using our own coherent symbolic source operator.
    hp=ROOT/'code/ambient_B18_C3_HB_right2_source_review_20261105_v63.py'
    assert sha(hp)=='d0f891265be64f5029ff2a6ec1f447f478ad1188b8d05450353964645ed505c9'
    s=importlib.util.spec_from_file_location('own_HB_symbolic',hp);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    affine=m.prepare_operator()['HB2_quartic_NF_affine_lambda_const_18_60_67']
    A0,B0,A1,B1=[list(map(Q,v)) for v in [affine[0][0],affine[0][1],affine[1][0],affine[1][1]]]
    C0=[b+a/4 for a,b in zip(A0,B0)];C1=[b+a/4 for a,b in zip(A1,B1)]
    F=[10*a+b for a,b in zip(A0,A1)]
    assert all(a+10*b==0 for a,b in zip(C1,C0))
    assert C0==[Q(-209831,2160),Q(0),Q(80,3),Q(40,3)]
    assert F==[Q(61,6),Q(100,3),-Q(40,3),Q(50)]
    factors=[]
    for v in [C0,F]:
        terms=[]
        for i,c in enumerate(v):
            if c:terms.append([['constant',0,0] if i==0 else ['lambda',[18,60,67][i-1],-1],str(c)])
        factors.append(terms)
    assert factors==gate['two_affine_factor_forms']
    assert [bound(t) for t in factors]==gate['two_affine_factor_majorants']
    # Universal determinant factorization on formal affine coefficients.
    def product(a,b):
        out={}
        for i,c in enumerate(a):
            for j,d in enumerate(b):
                key=tuple(sorted([i,j]));out[key]=out.get(key,Q(0))+c*d
        return out
    det=product(A0,B1);other=product(A1,B0);rhs=product(C0,F)
    for k,c in other.items():det[k]=det.get(k,Q(0))-c
    assert all(det.get(k,Q(0))==-rhs.get(k,Q(0)) for k in set(det)|set(rhs))
    return {'whole_gate_sha256':sha(path),'all_four_affine_U_D_height_and_parse_bounds_verified':True,
      'four_entry_heights':[[b['height'] for b in row] for row in bounds],
      'whole_postroot_operation_bound':675,'formal_factorization':'det=-C0*(10A0+A1)',
      'both_factor_forms_independently_reconstructed':True,
      'no_lambda_readback_or_actual_RIGHT2_value_computed':True,
      '2x2_inverse_bound_reason':'entry<=H, det<=4H+1, inverse<=5H+1; multiplying entry by inverse and adding two terms gives<=12H+3; rank1 ratios/factors fit<=6H+1.'}
def finalize():
    started=time.monotonic();out=ROOT/'reports'/ID;assert not(out/'report.json').exists();pins={}
    def pin(p,h=None):
        p=pathlib.Path(p);s=sha(p)
        if h is not None:assert s==h,(str(p),s,h)
        pins[str(p.relative_to(ROOT))]=s
        return read(p) if p.suffix=='.json' else None
    core=ROOT/'code/ambient_B18_C3_right2_lambda3_cache_source_review_20261105_v63.py';pin(core,CORE)
    cp=pin(out/'six-root-source-bound-checkpoint.json',CHECK)
    pp=pin(out/'postroot-affine-bound-and-formal-factor-checkpoint.json','4c0121bc4e1d16242238a9581e4e80ef438e61d5353998675616504ade5512fe')
    assert pp['whole_postroot_operation_bound']==675 and pp['no_lambda_readback_or_actual_RIGHT2_value_computed']
    pr=pin(ROOT/'reports'/BASE/'report.json',SOURCE);assert pr['all_passed']
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/BASE/p,h)
    custody=pin(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_VALUE_CUSTODY_20261105_v63/report.json',
        'fed89e1e40d200870350d1ece68f0af8c5c14087d2dc6d91fc1d10568b05aaad')
    for p,h in custody['input_sha256'].items():pin(ROOT/p,h)
    assert cp['literal_n_d_fixtures']==custody['all_literal_c_n_d_fixtures']
    assert len(cp['literal_n_d_fixtures'])==45 and sum(v[2]!=1 for v in cp['literal_n_d_fixtures'])==32
    c=read(ROOT/'reports'/BASE/'compiled-six-even-roots-and-provided-joint105-plane3-cache.json')
    for (row,ids),meta in zip(custody['exact_three_callback_all_four_source_root_bindings'],c['root_pairs_physical3']):
        assert row==meta['physical_row'] and ids==[meta[k] for k in ['aA_node','aB_node','bA_node','bB_node']]
    p=ROOT/'code/ambient_B18_C3_right2_quartic_cached_values_20261105_v63.py';code=p.read_text()
    pin(p,'a77b8597abf8b4a5af1b62bd42ec7983178badd8881ec13e5d5ea57cec01f740')
    f={n.name:n for n in ast.parse(code).body if isinstance(n,ast.FunctionDef)}
    for name,last in [('times','return note(a * b)'),('add','return note(a + b)'),('parse','return note(Q(s))')]:
        assert ast.unparse(f[name].body[-1])==last
        assert ast.unparse(f[name].body[-2].value).startswith('checked_height(')
    assert ast.unparse(f['checked_height'].body[-1])=='assert n <= 2048'
    for name in ['conv','inverse','mm']:
        calls={n.func.id for n in ast.walk(f[name]) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)}
        assert calls<={'Q','enumerate','range','zip','add','times','neg'}
        assert not any(isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Mult,ast.Div,ast.Pow))
           and not(isinstance(n.op,ast.Mult) and isinstance(n.left,ast.List)) for n in ast.walk(f[name]))
    assert ast.unparse(f['divide'].body[-1])=='return times(a, note(Q(b.denominator, b.numerator)))'
    assert ast.unparse(f['divide'].body[-2].value).startswith('checked_height(')
    required=['assert owned_bound<=75000<100000','assert raw_inputs+150000<=1000000',
      'signal.alarm(60)','resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))',
      "out.mkdir(exist_ok=False)","parser.add_argument('--source-peer-report',required=True)",
      'assert overlap_count==630','assert peak==162','assert delta[:5]==[Q(0)]*5',
      "elif op=='c':v=[parse(str(r[2])+'/'+str(r[3]))]+[Q(0)]*5",
      'assert MAX_PRECHECK<=79','assert MAX_PRECHECK<=gate[',
      'ROOT_ONLY_STOP_CHECKPOINT.json','ROOT_ONLY_POST_STOP_CHECKPOINT.json',
      'assert conv(v,cache[r[2]])==one==conv(cache[r[2]],v)',
      'assert determinant==neg(times(c0,f))','assert left==eye==right',
      'assert all(times(column[t],coefficients[k])==matrix[t][k]for t in range(2)for k in range(2))']
    for t in required:assert t in code,t
    assert code.index('assert MAX_PRECHECK<=79')<code.index('post=quartic_result(')
    assert code.index("phase1.write_text(")<code.index('post=quartic_result(')
    assert not(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_CACHED_VALUES_20261105_v63').exists()
    hp=ROOT/'code/ambient_B18_C3_alpha4_cached_source_review_20261105_v63.py'
    assert sha(hp)=='430cc27fc33b263cdcaec0fd616b87e88824bef702a3ff49fbc7136eaf428313'
    s=importlib.util.spec_from_file_location('own_counter',hp);own=importlib.util.module_from_spec(s);s.loader.exec_module(own)
    def bools(x):
        if isinstance(x,bool):return 1
        if isinstance(x,list):return sum(map(bools,x))
        if isinstance(x,dict):return sum(map(bools,x.values()))
        return 0
    def count(x):return own.nc(x)+bools(x)
    profile=read(ROOT/'reports'/BASE/'lambda3-cache-and-right2-whole-source-profile.json')
    gate=read(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_WHOLE_VALUE_GATE_20261105_v63/report.json')
    mathpeer=read(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63_REVIEW/report.json')
    op=read(ROOT/'reports/AMBIENT_B18_C3_RIGHT2_QUARTIC_SOURCE_20261105_v63/two-coherent-HB-quartic-sources-and-unpaired-affine-operator.json')
    c2=read(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    owned=sum(map(count,[pr,c,profile,gate,mathpeer,op,c2]))+2000+2940+735+7095+984+22000
    assert owned==custody['whole_future_owned_including2000_sourcepeer_reserve']==61309
    assert owned<75000<100000 and gate['source_whole_numeric_envelope_including_boundcache']<150000
    assert gate['whole_root_and_post_operation_bound']==675<2048
    # Independent tiny fixtures for both inverse directions and rank1 factors.
    def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
    A=[[Q(2),Q(3)],[Q(5),Q(7)]];iv=[[-Q(7),Q(3)],[Q(5),-Q(2)]]
    eye=[[Q(1),Q(0)],[Q(0),Q(1)]];assert mm(A,iv)==mm(iv,A)==eye
    B=[[Q(2),Q(6)],[Q(5),Q(15)]];assert all([Q(2),Q(5)][i]*[Q(1),Q(3)][j]==B[i][j] for i in range(2) for j in range(2))
    raw=sum((ROOT/p).stat().st_size for p in pins);assert raw+200000<32<<20
    evidence={'all7095_original_nodes_sixroots_and_same105frame_verified_by_fixed_checkpoint':True,
      'all630_input_and30_prior_root_overlaps_verified':True,'lastuse_peak':162,
      'all45_exact_n_d_fixtures_32_nonintegral_match':True,
      'full_coherent_RIGHT2_source_peer_and_all24_oldleaf_3lambda_weights_pinned':True,
      'root_phase_whole_QQ_bound':79,'four_affine_entry_heights':[48,49,54,56],
      'whole_postroot_bound':675,'hard_bits':2048,
      'before_add_mul_parse_divide_and_every_convolution_guards_verified':True,
      'fixed2x2det_both_inverse_products_rank1_all4_factor_checks_verified':True,
      'root_prefix_and_postroot_phase_checkpoint_no_rerun_semantics_verified':True,
      'whole_VALUE_owned_with2000_peerreserve':owned,'static_VALUE_cap':75000,'hard_VALUE_cap':100000,
      'SOURCE_whole_upper':134546,'hard_SOURCE_cap':150000,
      'hard_seconds':60,'hard_address_bytes':2<<30,'worker_RAW_cap':1000000,'hard_RAW_bytes':32<<20,
      'all3_lambda_values_and_actual_RIGHT2_matrix_rank_still_uncomputed':True,
      'no_source_solve_IFТ_gauge_physical_or_new_jet':True,
      'conditional_SAME_ONE_combined_cache_and_fixedquartic_worker_ready':True,
      'scope':'Rank2 proves exactC3length19/HF36640. Rank0or1 concerns only this two-source subset, never proves length20 or fullquarticrank.'}
    ep=out/'independent-sixroot-whole-code-lifetime-and-ONE-worker-readiness.json';save(ep,evidence)
    result={'id':ID,'all_passed':True,'frozen':True,'source_only':True,
      'input_sha256':pins,'core_code_sha256':CORE,'finish_code_sha256':sha(__file__),
      'evidence_sha256':{'six-root-source-bound-checkpoint.json':CHECK,
          'postroot-affine-bound-and-formal-factor-checkpoint.json':'4c0121bc4e1d16242238a9581e4e80ef438e61d5353998675616504ade5512fe',ep.name:sha(ep)},
      'prepared_worker_sha256':'a77b8597abf8b4a5af1b62bd42ec7983178badd8881ec13e5d5ea57cec01f740',
      'conditional_ONE_combined_worker_GO_ready':True,'no_lambda_values_evaluated':True,
      'elapsed_finish_seconds':time.monotonic()-started}
    save(out/'report.json',result);save(out/'phase.json',{'phase':'FINAL_SOURCE_CODE_BUDGET_PASS','no_value_by_reviewer':True})
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--finalize',action='store_true');a=parser.parse_args()
    signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    if a.finalize:finalize()
    else:
        out=ROOT/'reports'/ID;p=out/'postroot-affine-bound-and-formal-factor-checkpoint.json';assert not p.exists()
        save(p,post_bound_check());print(json.dumps({'postroot_SOURCE_bound_passed':True,'lambda_values_unread':True,'final_worker_custody_pending':True}))
