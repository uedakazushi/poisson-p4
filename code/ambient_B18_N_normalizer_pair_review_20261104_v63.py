#!/usr/bin/env python3
"""Independent sparse-Laurent and quotient-jet audit of the actual N pair."""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63';STEM=BASE+'_REVIEW'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def enc(a):
    if isinstance(a,F):return str(a)
    if isinstance(a,(list,tuple)):return [enc(x) for x in a]
    if isinstance(a,dict):return {str(k):enc(v) for k,v in a.items()}
    return a
def plus(a,b):
    v=dict(a)
    for m,c in b.items():
        q=v.get(m,F(0))+c
        if q:v[m]=q
        else:v.pop(m,None)
    return v
def neg(a):return {m:-c for m,c in a.items()}
def scale(a,c):return {m:c*v for m,v in a.items() if c*v}
def universal_identity():
    n=24;zero=(0,)*n;one={zero:F(1)}
    def var(i):a=list(zero);a[i]=1;return {tuple(a):F(1)}
    def mul(a,b):
        v={}
        for m,c in a.items():
            for q,d in b.items():
                e=tuple(x+y for x,y in zip(m,q));v[e]=v.get(e,F(0))+c*d
        return {m:c for m,c in v.items() if c}
    def power(a,k):
        v=one
        for j in range(k):v=mul(v,a)
        return v
    def derivative(a,i):
        v={}
        for m,c in a.items():
            if m[i]:
                e=list(m);e[i]-=1;v[tuple(e)]=c*m[i]
        return v
    h,x,y,z,w,A,B,g=[var(i) for i in range(8)];ainv=list(zero);ainv[5]=-1;ainv={tuple(ainv):F(1)}
    q=mul(plus(h,neg(mul(B,x))),ainv)
    mons=[m for m in itertools.combinations_with_replacement([0,1,3,4],3) if 3 in m or 4 in m]
    assert len(mons)==16
    H={}
    for i,mon in enumerate(mons):
        p=var(8+i)
        for j in mon:p=mul(p,var(j))
        H=plus(H,p)
    def substitute_h(p):
        out={}
        for m,c in p.items():
            e=list(m);k=e[0];e[0]=0
            out=plus(out,mul({tuple(e):c},power(q,k)))
        return out
    hp=substitute_h(H);f=mul(h,plus(mul(A,h),mul(B,x)))
    old=[[{} for j in range(5)] for i in range(5)]
    for i,j,v in [(1,2,f),(3,4,mul(g,power(h,2))),(2,3,neg(derivative(H,4))),(2,4,derivative(H,3))]:old[i][j]=v;old[j][i]=neg(v)
    M=[[one if i==j else {} for j in range(5)] for i in range(5)];M[0][0]=A;M[0][1]=B
    pushed=[[{} for j in range(5)] for i in range(5)]
    for i in range(5):
        for j in range(5):
            for a in range(5):
                for b in range(5):pushed[i][j]=plus(pushed[i][j],mul(mul(M[i][a],substitute_h(old[a][b])),M[j][b]))
    X=[h,x,y,z,w];bq=mul(B,q)
    for i in range(5):
        for j in range(5):
            if i==2:pushed[i][j]=plus(pushed[i][j],mul(bq,X[j]))
            if j==2:pushed[i][j]=plus(pushed[i][j],neg(mul(bq,X[i])))
    fB=mul(A,power(q,2));target=[[{} for j in range(5)] for i in range(5)]
    for i,j,v in [(1,2,fB),(3,4,mul(g,power(q,2))),
      (2,3,plus(neg(scale(mul(derivative(fB,1),z),F(1,2))),neg(derivative(hp,4)))),
      (2,4,plus(neg(scale(mul(derivative(fB,1),w),F(1,2))),derivative(hp,3)))]:target[i][j]=v;target[j][i]=neg(v)
    assert all(pushed[i][j]==target[i][j] for i in range(5) for j in range(5))
    assert all(sum(e[:5])==2 and all(e[k]>=0 for k in list(range(5))+list(range(6,24))) for row in pushed for p in row for e in p)
    return {'independent_full25_raw_Laurent_identities':True,'all150_quadratic_physical_coefficients_equal':True,
      'H_coefficients':16,'only_parameter_inverse':'A','denominator_B_g_or_conductor_inverted':False}

JM=[()]+[(i,) for i in range(8)]+[(i,7) for i in range(7)]+[(7,7)];ALLOWED=set(JM)
def jmul(a,b):
    v={}
    for m,c in a.items():
        for n,d in b.items():
            k=tuple(sorted(m+n))
            if k in ALLOWED:v[k]=v.get(k,F(0))+c*d
    return {m:c for m,c in v.items() if c}
def jet_eval(nodes,roots,binding):
    todo=list(roots);need=set()
    while todo:
        i=todo.pop()
        if i in need:continue
        need.add(i)
        if nodes[i][0] not in ['c','v']:todo.extend(nodes[i][1:])
    v={}
    for i in sorted(need):
        n=nodes[i];o=n[0]
        if o=='c':q={():F(n[1],n[2])} if n[1] else {}
        elif o=='v':q=binding[i]
        elif o=='a':q=plus(v[n[1]],v[n[2]])
        elif o=='n':q=neg(v[n[1]])
        elif o=='m':q=jmul(v[n[1]],v[n[2]])
        elif o=='i':
            a=v[n[1]];c=a.get((),F(0));assert c
            t=dict(a);t.pop((),None)
            q=plus({():1/c},plus(scale(t,-1/c**2),scale(jmul(t,t),1/c**3)))
        else:raise ValueError(o)
        v[i]=q
    return v,len(need)
def compose(a,sub):
    m=[()]+[(i,) for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2));out={}
    for mon,c in zip(m,map(F,a)):
        if not c:continue
        v={():c}
        for i in mon:v=jmul(v,sub[i])
        out=plus(out,v)
    return out
def array(p):return [p.get(m,F(0)) for m in JM]
def rref(M):
    A=[list(r) for r in M];piv=[];row=0
    for j in range(len(A[0])):
        k=next((i for i in range(row,len(A)) if A[i][j]),None)
        if k is None:continue
        A[row],A[k]=A[k],A[row];q=A[row][j];A[row]=[x/q for x in A[row]]
        for i in range(len(A)):
            if i!=row and A[i][j]:
                q=A[i][j];A[i]=[x-q*y for x,y in zip(A[i],A[row])]
        piv.append(j);row+=1
        if row==len(A):break
    return A,piv
def det(M):
    A=[list(r) for r in M];d=F(1)
    for j in range(len(A)):
        k=next(i for i in range(j,len(A)) if A[i][j])
        if k!=j:A[j],A[k]=A[k],A[j];d=-d
        p=A[j][j];d*=p
        for i in range(j+1,len(A)):
            q=A[i][j]/p
            for h in range(j,len(A)):A[i][h]-=q*A[j][h]
    return d
def main(out):
    start=time.monotonic();signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    folder=ROOT/'reports'/BASE;ap=folder/'actual-N-pair-preflight-evidence.json';ar=folder/'report.json'
    np=ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    gates=[ar,ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',
      ROOT/'reports/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json',
      ROOT/'reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW.json']
    assert all(read(p)['all_passed'] for p in gates)
    A,N,E,S=map(read,[ap,np,ep,sp]);universal=universal_identity()
    bind={i:{m:F(c) for m,c in zip(JM,p) if F(c)} for i,p in N['all28_mixed_source2jets']}
    roots=N['frame28_equation_nodes']+N['actual_same_R0_readback_nodes9']+N['physical150_N_source_nodes']
    v,count=jet_eval(N['DAG_scalar_nodes'],roots,bind)
    frame=N['frame28_equation_nodes'];want=[{}]*20+[{(i,):F(1)} for i in range(8)]
    assert [v[i] for i in frame]==want
    rho=[v[i] for i in N['actual_same_R0_readback_nodes9']]
    assert [array(p) for p in rho]==[[F(x) for x in p] for p in N['all9_R0_mixed2jets']]
    assert all(p.get((7,),F(0))==0 for p in rho)
    found={};tau={}
    for family in ['CA','CB']:
        pairs={r:pp for r,pp in E['families'][family]['all12_trace_pairs_in_common_R0_through2']}
        nu=compose(E['families'][family]['nu_in_common_R0_through2'],rho)
        assert all(not nu.get(m,F(0)) for m in JM[:9])
        checks=[];ratios=[]
        for row in [18,19,28,29,60,67,72,63,64,68,73,74]:
            even,b=map(lambda p:compose(p,rho),pairs[row])
            num=plus(v[N['physical150_N_source_nodes'][row]],neg(even))
            assert not any(num.get(m,F(0)) for m in JM[:9])
            for j in range(7):
                if b.get((j,),F(0)):ratios.append(num.get((j,7),F(0))/b[(j,)])
            checks.append({'row':row,'b_pullback17':array(b),'numerator17':array(num)})
        assert ratios and all(t==ratios[0] for t in ratios) and ratios[0]
        t=ratios[0];assert nu.get((7,7),F(0))==t*t
        tau[family]=[nu.get((j,7),F(0))/(2*t) for j in range(7)]+[t]
        for c in checks:
            beta=dict(zip(JM,c['b_pullback17']));num=dict(zip(JM,c['numerator17']))
            assert all(num[(j,7)]==beta[(j,)]*t for j in range(7)) and num[(7,7)]==0
        author=A['families'][family]
        assert tau[family]==list(map(F,author['tau_firstjet_in_N_source8']))
        assert array(nu)==list(map(F,author['nu_pullback_mixed2']))
        assert [[c['row'],c['b_pullback17'],c['numerator17']] for c in checks]==[[c['physical_row'],list(map(F,c['conductor_pullback_mixed2'])),list(map(F,c['target_minus_even_mixed2']))] for c in author['all12_target_minus_even_and_conductor_coefficients']]
        finite=A['actual_finite_regular_function_input'][family];sf=S['families'][family]
        xid=[next(i for i,n in enumerate(S['DAG_scalar_nodes']) if n==['v',f'x{r}']) for r in N['actual_same_R0_rows9']]
        assert finite['source_R0_nodes_to_N_source_rho_nodes']==[[i,j] for i,j in zip(xid,N['actual_same_R0_readback_nodes9'])]
        assert [finite['a18_source_node'],finite['denominator_b18_source_node']]==sf['complete_normalized_physical_w_coefficients'][0][18]
        assert finite['square_nu_source_node']==sf['actual_tracefree_nu'] and not finite['unit_inverse_of_b18_used']
        found[family]={'nu_pullback17':array(nu),'all12_checked_target_and_conductor_jets':checks,
          'tau_firstjet':tau[family],'b18_nonzero_X108_coefficient':checks[0]['b_pullback17'][1]}
    assert found['CA']['b18_nonzero_X108_coefficient']==2
    assert found['CB']['b18_nonzero_X108_coefficient']==F(-1,2)
    matrix=[[p.get((j,),F(0)) for j in range(8)] for p in rho[1:8]]+[tau['CA'],tau['CB']]
    assert matrix==[[F(c) for c in r] for r in A['Npair_to_Schur9_firstjet']]
    rr,piv=rref([list(c) for c in zip(*matrix)]);assert len(piv)==8
    minor=det([matrix[i] for i in piv]);assert minor==F(-1,2)
    free=next(i for i in range(9) if i not in piv);ell=[F(0)]*9;ell[free]=1
    for i,j in enumerate(piv):ell[j]=-rr[i][free]
    ell=[x*6/ell[8] for x in ell]
    assert ell==[F(0),F(0),F(5),F(0),F(-20),F(0),F(0),F(24),F(6)]
    evidence={'independent_universal_source_identity':universal,'own_sparse_N_source_ancestor_count':count,
      'all28_cached_frame_identities17_verified':True,'all150_N_physical_outputs_rebuilt':True,
      'families':found,'own_Npair_firstjet9_by8':matrix,'own_rank8_selected_rows':piv,'own_minor8':minor,
      'own_normal_linear_equation':ell,
      'normal_regular_quotient_proof':'actual zero N map by nonnegative torus degeneration; b18!=0 defines rational lift on normalizer-equals-image open; monic square makes root integral in normal B_N, hence regular',
      'all_order_factor_proof':'same physical outputs agree on dense conductor open then on whole normal source; rank8 makes image a smooth hypersurface, so every defining F_j belongs its principal ideal',
      'residual4_core_consequence':'I_E4=L_formal*(u5,h1..h7(core4)) after formal coordinate change; retain all L_formal*u5 equations',
      'core_functions_full_E4_socle_or_D0_depth_computed':False}
    op=out/'independent-N-pair-and-root-proof.json';op.write_text(json.dumps(enc(evidence),separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',ap,ar,np,ep,sp,*gates,
      ROOT/'code/ambient_B18_N_normalizer_pair_preflight_20261104_v63.py',ROOT/'research/notes'/f'{BASE}.md']
    rep={'id':STEM,'all_passed':True,'independent_Laurent150_and_sparse_N17_arithmetic':True,
      'actual_generic_conductor_open_and_integral_closure_root_proof_accepted':True,
      'actual_N_pair_smooth_hypersurface_and_all12_L_divisions_accepted':True,
      'rank8_minor':'-1/2','no_producer_rerun_new_source_jets_GB_or_denominator_unit_assumption':True,
      'full_E4_core_socle_or_actual_D_depth_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{op.name:sha(op)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','rank8_minor','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
