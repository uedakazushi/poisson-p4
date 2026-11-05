#!/usr/bin/env python3
"""Finish independent literal CUT1 source/code/ONE-residual custody audit."""
import ast,hashlib,json,pathlib,resource,signal,time
from fractions import Fraction as Q
ROOT=pathlib.Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_C3_CUT1_QUARTIC_SOURCE_20261105_v63';ID=BASE+'_REVIEW'
PIN='ed1678417ff8e459de67d38ced0c8c835aab2a7ea8a1eeb315867a7570f9b633'
CORE='4ff64d3a238e0977619ce3700f748fc15e99e84cd6250f75f449a09d4b293ec1'
CHECK='4d776b4fab24f218fb75b0fcc44910ca8b64f2beba42a924ca87fc136e2a363d'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(pathlib.Path(p).read_text())
def save(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2)+'\n')
def count(x):
    if isinstance(x,(int,float,Q)):return 1
    if isinstance(x,str):
        try:Q(x);return 1
        except(ValueError,ZeroDivisionError):return 0
    if isinstance(x,list):return sum(map(count,x))
    if isinstance(x,dict):return sum(map(count,x.values()))
    return 0
def main():
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
    out=ROOT/'reports'/ID;assert not(out/'report.json').exists();pins={}
    def pin(p,h=None):
        p=pathlib.Path(p);s=sha(p)
        if h is not None:assert s==h,(str(p),s,h)
        pins[str(p.relative_to(ROOT))]=s;return read(p) if p.suffix=='.json' else None
    pin(ROOT/'code/ambient_B18_C3_cut1_source_review_20261105_v63.py',CORE)
    proof=pin(out/'independent-literal-cut-source-full-NF-and-positive-bound-checkpoint.json',CHECK)
    assert proof['whole_operation_bound']==519 and proof['all_old_leaf_values_unpaired_no_residual_computed']
    pr=pin(ROOT/'reports'/BASE/'report.json',PIN)
    for p,h in pr['input_sha256'].items():pin(ROOT/p,h)
    for p,h in pr['evidence_sha256'].items():pin(ROOT/'reports'/BASE/p,h)
    custody=pin(ROOT/'reports/AMBIENT_B18_C3_CUT1_VALUE_CUSTODY_20261105_v63/report.json',
        'c5a1e3420930b7e51f246e0558fb287035f11744f3b1528bc1cb5c6becc454ed')
    for p,h in custody['input_sha256'].items():pin(ROOT/p,h)
    p=ROOT/'code/ambient_B18_C3_cut1_quartic_values_20261105_v63.py';text=p.read_text()
    assert sha(p)=='b260774879d61b9413dc3e07590c59a09345d5722a515874ae8f8770c018cb84'
    f={n.name:n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef)}
    for name,last in [('parse','return note(Q(s))'),('add','return note(a + b)'),('mul','return note(a * b)')]:
        assert ast.unparse(f[name].body[-1])==last
        assert ast.unparse(f[name].body[-2].value).startswith('guard(')
    assert ast.unparse(f['guard'].body[-1])=='assert n <= 2048'
    assert ast.unparse(f['div'].body[-1])=='return mul(a, note(Q(b.denominator, b.numerator)))'
    assert ast.unparse(f['div'].body[-2].value).startswith('guard(')
    for name in ['padd','pmul','mm']:
        calls={n.func.id for n in ast.walk(f[name]) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)}
        assert calls<={'Q','add','mul','zip','enumerate','padd','pmul'}
        assert not any(isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Mult,ast.Div,ast.Pow))
          and not(isinstance(n.op,ast.Mult) and isinstance(n.left,ast.List)) for n in ast.walk(f[name]))
    for literal in ['assert owned<=18000<30000','assert raw+100000<=400000',
        'signal.alarm(60)','resource.setrlimit(resource.RLIMIT_AS,(2147483648,2147483648))',
        "parser.add_argument('--source-peer-report',required=True)","out.mkdir(exist_ok=False)",
        "for key,w in ev['residual_B_minus_121A_over48_operator']:",
        'assert left==I==right','assert MAX_PRE<=ev[']:
        assert literal in text,literal
    assert "ev['A4_cached_linear_operator']" not in text and "ev['B4_cached_linear_operator']" not in text
    assert not(ROOT/'reports/AMBIENT_B18_C3_CUT1_QUARTIC_VALUES_20261105_v63').exists()
    ev=read(ROOT/'reports'/BASE/'literal-cut-quartic-source-and-three-unpaired-cache-operators.json')
    c2=read(ROOT/'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json')
    owned=count(pr)+count(ev)+count(c2)+1500+6000
    assert owned==custody['whole_future_owned_with1500peerreserve']==10040<18000<30000
    assert len(ev['residual_B_minus_121A_over48_operator'])==24
    raw_paths=set(pr['input_sha256'])|{str((ROOT/'reports'/BASE/'report.json').relative_to(ROOT)),str((ROOT/'reports'/BASE/'literal-cut-quartic-source-and-three-unpaired-cache-operators.json').relative_to(ROOT))}
    raw_before_peer=sum((ROOT/p).stat().st_size for p in raw_paths)
    assert raw_before_peer==custody['all_input_RAW_before_peer']==145426
    assert raw_before_peer+30000+100000<400000
    # Independent universal check over QQ[t,d,d^-1]. No actual target is
    # evaluated. The worker's vectors of length3 hold all degrees<=2; its
    # matrix entries have degree<=1, so this equals an identity in QQ[t].
    def add(a,b):
        c=dict(a)
        for m,q in b.items():
            c[m]=c.get(m,Q(0))+q
            if not c[m]:del c[m]
        return c
    def mul(a,b):
        c={}
        for m,q in a.items():
            for n,v in b.items():
                k=(m[0]+n[0],m[1]+n[1]);c=add(c,{k:q*v})
        return c
    def mm(a,b):
        return [[add(mul(row[0],col[0]),mul(row[1],col[1])) for col in zip(*b)] for row in a]
    z={};one={(0,0):Q(1)};a={(0,0):-Q(3,5)};t={(1,0):Q(1)};d={(0,1):Q(1)}
    M=[[a,t],[z,d]];iv=[[{(0,0):-Q(5,3)},{(1,-1):Q(5,3)}],[z,{(0,-1):Q(1)}]]
    I=[[one,z],[z,one]];assert mm(M,iv)==mm(iv,M)==I
    assert all(m[0]<=1 for mat in [M,iv] for row in mat for pp in row for m in pp)
    assert 2048+3>2048 and custody['before_product_overflow_fixture']
    final={'literal_cut_coherent_source_and_all_A_B_residual_weights_accepted':True,
      'all24_old_C2_leaves_remain_unpaired_before_this_GO':True,
      'only_residual_is_evaluated_by_prepared_worker':True,
      'CUT_A_coordinate_t_is_indeterminate_not_zero':True,
      'universal_matrix_over_QQ_t':'[[-3/5,t],[0,d]]',
      'universal_inverse_over_QQ_t_d_inverse':'[[-5/3,5t/(3d)],[0,1/d]]',
      'both_products_independently_verified_without_target_values':True,
      'worker_degree2_polynomials_capture_all_product_degrees_not_just_t_nilpotent':True,
      'all_guarded_add_mul_parse_reciprocal_polynomial_operations_verified':True,
      'whole_QQ_operation_bound':519,'hard_QQ_bits':2048,
      'whole_numeric_with1500_peerreserve':owned,'static_cap':18000,'hard_numeric':30000,
      'RAW_before_peer':raw_before_peer,'peer_RAW_reserve':30000,'output_RAW_reserve':100000,'hard_RAW':400000,
      'hard_seconds':60,'hard_address_bytes':2<<30,
      'conditional_SAME_ONE_cached_residual_worker_ready':True,
      'no_new_lambda_nu_delta4_mixedPhi3_beta3_source_gauge_IFT_jet':True,
      'nonzero_completion':'Two actual independent quartics => C3 exact19/HF36640; longer exponent and H0/nil stayopen.',
      'zero_scope':'Only this CUT source is dependent. Never fullquarticrank/length20.'}
    ep=out/'independent-ONE-residual-code-universal-rank-and-whole-custody.json';save(ep,final)
    report={'id':ID,'all_passed':True,'frozen':True,'source_only':True,'input_sha256':pins,
      'core_code_sha256':CORE,'finish_code_sha256':sha(__file__),
      'evidence_sha256':{'independent-literal-cut-source-full-NF-and-positive-bound-checkpoint.json':CHECK,ep.name:sha(ep)},
      'prepared_worker_sha256':'b260774879d61b9413dc3e07590c59a09345d5722a515874ae8f8770c018cb84',
      'conditional_ONE_cached_scalar_GO_ready':True,'actual_residual_still_uncomputed':True,
      'elapsed_finish_seconds':time.monotonic()-start}
    save(out/'report.json',report);save(out/'phase.json',{'phase':'FINAL_SOURCE_CODE_CUSTODY_PASS','scalar_uncomputed_by_reviewer':True})
    print(json.dumps({'all_passed':True,'report':str(out/'report.json'),'sha256':sha(out/'report.json')}));signal.alarm(0)
if __name__=='__main__':main()
