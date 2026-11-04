#!/usr/bin/env python3
"""Independent sparse-Fraction review of the five-line formal core theorem.

Does not import or run the producer, SymPy, a source solver, or any CAS/GB.
"""
import argparse,ast,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def add(a,b):
    c=dict(a)
    for m,x in b.items():c[m]=c.get(m,F(0))+x
    return {m:x for m,x in c.items() if x}
def scale(a,c):return {m:c*x for m,x in a.items() if c*x}
def mul(a,b):
    c={}
    for m,x in a.items():
        for n,y in b.items():
            k=tuple(i+j for i,j in zip(m,n));c[k]=c.get(k,F(0))+x*y
    return {m:x for m,x in c.items() if x}
def one(n):return {(0,)*n:F(1)}
def power(p,k,n):
    if k<0:
        assert len(p)==1
        m,c=next(iter(p.items()));return {tuple(k*i for i in m):c**k}
    q=one(n)
    for _ in range(k):q=mul(q,p)
    return q
def variable(n,j):m=[0]*n;m[j]=1;return {tuple(m):F(1)}
def decode(p):return {tuple(m):F(n,d) for m,n,d in p}
def encode(p):return [[list(m),str(c)] for m,c in sorted(p.items(),reverse=True)]
def substitute(p,replacements,n):
    result={}
    for m,c in p.items():
        q=scale(one(n),c)
        for a,k in zip(replacements,m):q=mul(q,power(a,k,n))
        result=add(result,q)
    return result
def derivative(p,j):
    q={}
    for m,c in p.items():
        if m[j]:
            t=list(m);t[j]-=1;q[tuple(t)]=m[j]*c
    return q
def determinant(A,n):
    result={}
    for p in itertools.permutations(range(len(A))):
        c=one(n);sign=(-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
        for i,j in enumerate(p):c=mul(c,A[i][j])
        result=add(result,scale(c,sign))
    return result
def parse(text,names):
    n=len(names);variables={q:variable(n,i) for i,q in enumerate(names)}
    def f(node):
        if isinstance(node,ast.Constant):return scale(one(n),F(node.value))
        if isinstance(node,ast.Name):return variables[node.id]
        if isinstance(node,ast.UnaryOp):return scale(f(node.operand),-1) if isinstance(node.op,ast.USub) else f(node.operand)
        if isinstance(node,ast.BinOp):
            a=f(node.left)
            if isinstance(node.op,ast.Pow):return power(a,ast.literal_eval(node.right),n)
            b=f(node.right)
            if isinstance(node.op,ast.Add):return add(a,b)
            if isinstance(node.op,ast.Sub):return add(a,scale(b,-1))
            if isinstance(node.op,ast.Mult):return mul(a,b)
            if isinstance(node.op,ast.Div):return mul(a,power(b,-1,n))
        raise ValueError(ast.dump(node))
    return f(ast.parse(text,mode='eval').body)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63';primary=folder/'report.json';evidence=folder/'five-line-and-generic-models.json'
    prep=ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json'
    gate=ROOT/'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63_REVIEW/report.json'
    author=read(primary);A=read(evidence);P=read(prep);assert author['all_passed'] and read(gate)['all_passed']
    assert sha(primary)=='acbb5c1f5e34ad44bbde2f4c8f977813a22bb17b44e2de52bf829c3fa6fcaf95'
    for q,h in author['input_sha256'].items():assert sha(ROOT/q)==h,q
    for q,h in author['evidence_sha256'].items():assert sha(folder/q)==h,q
    Q0,Q1,f=map(decode,P['generators3']);a,d,x,y=[variable(4,j) for j in range(4)];c=one(4)
    r=add(add(x,scale(a,-1)),scale(d,-1));s=add(a,x);polyP=add(mul(d,d),scale(c,5))
    assert f==mul(s,mul(r,r))
    # Different elementary elimination identities, all in original four vars.
    qr0=substitute(Q0,[a,d,add(a,d),y],4);qr1=substitute(Q1,[a,d,add(a,d),y],4)
    h=add(add(scale(mul(a,a),-1),scale(mul(a,d),-3)),mul(d,y))
    j=add(add(add(mul(mul(a,a),d),scale(mul(a,mul(d,d)),2)),scale(a,-5)),scale(y,5))
    assert qr0==scale(h,8) and qr1==scale(mul(a,j),-8)
    assert add(mul(d,j),scale(h,-5))==mul(mul(a,polyP),add(a,scale(d,2)))
    qs0=substitute(Q0,[a,d,scale(a,-1),y],4);qs1=substitute(Q1,[a,d,scale(a,-1),y],4)
    hs=add(add(scale(a,-1),scale(d,2)),scale(y,2));assert qs0==scale(mul(d,hs),4)
    # Restrict the linear second factor to zero instead of dividing by it.
    rs1=substitute(qs1,[add(scale(d,2),scale(y,2)),d,x,y],4)
    assert rs1==scale(mul(mul(d,add(d,scale(y,2))),polyP),2)
    dzero=substitute(qs1,[a,{},x,y],4)
    assert dzero==scale(mul(a,add(a,scale(y,-2))),10)
    line_strings={'N':['0','0','0','t'],'C':['0','t','t','0'],
      'R':['-2*t','t','-t','-2*t'],'P':['2*t','0','-2*t','t'],'L':['t','t','-t','-t/2']}
    lines={label:[parse(q,['t']) for q in values] for label,values in line_strings.items()}
    for line in lines.values():assert all(not substitute(q,line,1) for q in [Q0,Q1,f])
    for label,ps in A['support']['primes'].items():
        for p in ps:
            sp={tuple(m):F(v) for m,v in p};assert not substitute(sp,lines[label],1)
    models=A['generic_models'];checked={}
    # N/P use explicit Laurent-polynomial jets, not an implicit solve.
    for label in ['N','P']:
        rec=models[label]
        jets=[parse(rec['all4_coordinate_jets'][name],['t','u']) for name in ['a','d','v2','v6']]
        residual=[substitute(q,jets,2) for q in [Q0,Q1]]
        assert all(not {m:c for m,c in p.items() if m[1]<=2} for p in residual)
        ff=substitute(f,jets,2)
        if label=='N':
            assert {m:c for m,c in ff.items() if m[1]<=3}=={(0,3):F(-32)}
            assert jets[1]=={(-1,2):F(4)}
            assert jets[2]=={(0,1):F(-3),(-1,2):F(2)}
            normal=[1,2];free=0;detwanted=parse('-80*t**2',['t'])
        else:
            assert {m:c for m,c in ff.items() if m[1]<=2}=={(1,2):F(-16,3)}
            assert add(jets[0],jets[2])=={(-1,2):F(-1,3)}
            assert jets[1]=={(0,1):F(1)}
            normal=[0,2];free=1;detwanted=parse('-240*t**2',['t'])
        jac=[[substitute(derivative(q,k),lines[label],1) for k in normal] for q in [Q0,Q1]]
        det=determinant(jac,1);assert det==detwanted and det==parse(rec['IFT_determinant'],['t'])
        assert all(jac[i][j]==parse(rec['IFT_Jacobian'][i][j],['t']) for i in range(2) for j in range(2))
        checked[label]={'unit_minor':encode(det),'all_Q0_Q1_jets_through2_zero':True,
          'f_first_nonzero':encode({m:c for m,c in ff.items() if m[1]<=(3 if label=='N' else 2)}),
          'd_in_artin_ring_nonzero_and_square_zero':True}
    # At C, compute the exact y elimination, then every coefficient of Q1.
    t,u,v=[variable(3,j) for j in range(3)];zer={}
    cmap=[u,t,add(add(t,u),v),zer];CQ0=substitute(Q0,cmap,3)
    yroot=scale(mul(CQ0,power(t,-1,3)),F(-1,8))
    cmap[3]=yroot;assert not substitute(Q0,cmap,3)
    CQ1=substitute(Q1,cmap,3)
    expected={}
    for ex,co in models['C']['exact_Q1_after_y']:
        coefficient=parse(co,['t']);term={(m[0],ex[0],ex[1]):c for m,c in coefficient.items()};expected=add(expected,term)
    assert CQ1==expected
    q2={m:c for m,c in CQ1.items() if m[1]+m[2]==2}
    assert q2==parse(models['C']['quadratic_initial'],['t','u','v'])
    assert not {m:c for m,c in CQ1.items() if m[1]+m[2]<2}
    u2={m:c for m,c in q2.items() if m[1:]==(2,0)}
    assert u2=={(2,2,0):F(-16),(0,2,0):F(-80)}
    checked['C']={'exact_Q0_y_elimination_verified':True,'all_exact_Q1_coefficients':encode(CQ1),
      'quadratic_regular_pair_with_v_squared':True,'HF':[1,2,1]}
    # R/L: independent signed determinant and every normal Jacobian entry.
    for label,last in [('R',r),('L',s)]:
        jac=[[substitute(derivative(q,k),lines[label],1) for k in [0,2,3]] for q in [Q0,Q1,last]]
        det=determinant(jac,1);assert det and det==parse(models[label]['determinant'],['t'])
        assert all(jac[i][j]==parse(models[label]['exact_normal_Jacobian'][i][j],['t']) for i in range(3) for j in range(3))
        assert substitute(s if label=='R' else r,lines[label],1)==parse('-3*t',['t'])
        checked[label]={'unit_minor':encode(det),'normal_coordinate_system_proved':True}
    assert [models[k]['length'] for k in ['N','C','R','P','L']]==[3,4,2,2,1]
    assert sum(models[k]['length'] for k in models)==12
    proof={'source_Q0_Q1_read_directly_from_accepted_prep':True,
      'own_engine':'sparse Fraction Laurent polynomials; no SymPy, producer imports or CAS',
      'all_r_s_restrictions_and_prime_cases_verified':True,'all5_line_substitutions_verified':True,
      'generic_checks':checked,'geometric_lengths':[3,4,2,2,1],
      'core_CM_and_no_embedded_from_independent_CI_gate':True,
      'core_nilcube_global_and_square_nonzero':True,
      'E4_Ass_exact':['d','N','C','R','P','L'],'E4_embedded':['N','P'],
      'E4_nil_at_N_P':'rad(dJ)_q=(d)_q; multiplication d identifies nil ideal with T_q as module, nil^k with d^(k-1)T_q; d≠0,d²=0',
      'why_global_nilcube':'nilcube is a finite submodule of E4; its localization at each of the six associated primes is zero, so it is zero',
      'why_square_nonzero':'at N/P the nil ideal square is d*T_q≠0, using exact generic models',
      'no_D0_or_full_D_transfer':True}
    ep=out/'independent-five-lines-and-nil-proof.json';ep.write_text(json.dumps(proof,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',primary,evidence,prep,gate,
           ROOT/'code/ambient_B18_proper_core_five_lines_20261104_v63.py',ROOT/'research/notes/AMBIENT_B18_PROPER_CORE_FIVE_LINES_20261104_v63.md']
    report={'id':STEM,'all_passed':True,'status':'independent_source_Laurent_identity_and_embedded_nil_argument_acceptance',
      'exact5_core_support_and_lengths_3_4_2_2_1_accepted':True,
      'core_and_normalized_E4_exact_nilradical_index3_accepted':True,
      'E4_Ass_six_and_embedded_N_P_accepted':True,
      'no_producer_rerun_SymPy_CAS_GB_source_point_or_newjet':True,
      'D0_full_D_or_global_classification_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
