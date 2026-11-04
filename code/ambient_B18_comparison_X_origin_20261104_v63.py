#!/usr/bin/env python3
"""Actual constant full comparison X, using only saved source data and C0 length4.
No fresh source/implicit solve, jet, point or GB. Native inverses are unit
finite-algebra inverses of the already specified rational source operators.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_CA_critical_V1_relative_tangent_review_20261104_v63 import action,PAIRS,MONS,PI,MI
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63'
Z=(F(0),)*4;ONE=(F(1),F(0),F(0),F(0));TA=(F(0),F(1),F(0),F(0));TB=(F(0),F(0),F(1),F(0))
INVERSES=[]
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def c(q):return (F(q),F(0),F(0),F(0))
def add(*ps):return tuple(sum(x,F(0)) for x in zip(*ps)) if ps else Z
def neg(p):return tuple(-x for x in p)
def scale(p,q):return tuple(q*x for x in p)
def mul(*ps):
    out=ONE
    for p in ps:
        q=[F(0)]*4
        for i,x in enumerate(out):
            if x:
                for j,y in enumerate(p):
                    if y and not(i&j):q[i^j]+=x*y
        out=tuple(q)
    return out
def power(p,n):
    if n<0:return power(inv(p),-n)
    out=ONE
    for _ in range(n):out=mul(out,p)
    return out
def inv(p):
    assert p[0],('NONUNIT',p)
    n=scale(add(p,c(-p[0])),-1/p[0]);v=scale(add(ONE,n,mul(n,n)),1/p[0]);assert mul(p,v)==ONE
    INVERSES.append([list(map(str,p)),list(map(str,v))]);return v
def ident(n):return [[ONE if i==j else Z for j in range(n)] for i in range(n)]
def mm(A,B):
    out=[[Z for _ in B[0]] for _ in A]
    for i,row in enumerate(A):
        for k,p in enumerate(row):
            if p!=Z:
                for j,q in enumerate(B[k]):
                    if q!=Z:out[i][j]=add(out[i][j],mul(p,q))
    return out
def mi(A):
    n=len(A);a=[r[:]+s for r,s in zip(A,ident(n))]
    for j in range(n):
        k=next(i for i in range(j,n) if a[i][j][0])
        if k!=j:a[k],a[j]=a[j],a[k]
        u=inv(a[j][j]);a[j]=[mul(x,u) for x in a[j]]
        for i in range(n):
            if i!=j and a[i][j]!=Z:
                t=a[i][j];a[i]=[add(x,neg(mul(t,y))) for x,y in zip(a[i],a[j])]
    assert [r[:n] for r in a]==ident(n);out=[r[n:] for r in a];assert mm(A,out)==mm(out,A)==ident(n);return out
def constmat(A):return [[c(x) for x in row] for row in A]
def evalpoly(terms,values):
    out=Z
    for rec in terms:
        if len(rec)==3:m,n,d=rec;q=F(n,d)
        else:m,q=rec;q=F(q)
        v=c(q)
        for p,k in zip(values,m):
            if k:v=mul(v,power(p,k))
        out=add(out,v)
    return out
def derivative(terms,j):
    result=[]
    for m,n,d in terms:
        if m[j]:
            q=list(m);q[j]-=1;result.append([q,n*m[j],d])
    return result
def telescope(terms,values,indices,roots):
    # Ordered polynomial difference quotients, valid even when Xi=Yi.
    result=[Z]*len(indices);position={j:i for i,j in enumerate(indices)}
    for ex,n,d in terms:
        for k,j in enumerate(indices):
            if not ex[j]:continue
            t=Z
            for h in range(ex[j]):t=add(t,mul(power(values[j],ex[j]-1-h),power(roots[k],h)))
            t=scale(t,F(n,d))
            for q,e in enumerate(ex):
                if q==j or not e:continue
                p=roots[position[q]] if q in position and position[q]>k else values[q]
                t=mul(t,power(p,e))
            result[k]=add(result[k],t)
    return result
def frame(H):
    I=mi(H);ell=[mul(x,inv(I[1][1])) for x in I[1]]
    k=[add(x,neg(mul(I[0][1],y))) for x,y in zip(I[0],ell)];den=inv(k[0]);k=[mul(x,den) for x in k]
    Fm=ident(5);Fm[0]=k;Fm[1]=ell;return Fm
def physical_push(v,G,Gi):
    out=[Z]*150
    for t,p in enumerate(v):
        if p==Z:continue
        a,b=PAIRS[t//15];u,w=MONS[t%15];pol={}
        for i,j in itertools.product(range(5),repeat=2):
            q=mul(p,Gi[u][i],Gi[w][j])
            if q!=Z:
                m=tuple(sorted((i,j)));pol[m]=add(pol.get(m,Z),q)
        for i,j in PAIRS:
            f=add(mul(G[i][a],G[j][b]),neg(mul(G[i][b],G[j][a])))
            if f!=Z:
                for m,q in pol.items():out[15*PI[i,j]+MI[m]]=add(out[15*PI[i,j]+MI[m]],mul(f,q))
    return out
def serial(A):return [[list(map(str,p)) for p in row] for row in A]

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'interface':'reports/AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63/actual-source-physical-firstjet-interface.json',
      'proper':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json',
      'equations':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63/actual-proper-source14-equations17.json',
      'first5':'reports/AMBIENT_B18_PROPER_RATIONAL_ELIMINATION_20261104_v63/actual12-flag-divided-polynomials.json',
      'second5':'reports/AMBIENT_B18_PROPER_CORE_ELIMINATION_20261104_v63/actual4variable7core-polynomials.json',
      'membership':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/membership-and-source-cofactors.json',
      'preparation':'reports/AMBIENT_B18_PROPER_CORE_CI_MEMBERSHIP_20261104_v63/preparation.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'firstjets':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'source105':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
      'normalizer':'reports/AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63/actual-normalizer-tensor-free4-DAG.json',
      'fullV':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
      'Qoperator':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
      'quartic_firsts':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json'}
    D={k:read(ROOT/p) for k,p in paths.items()};I,P,Q,H,K,C,prep,B,C0,F9,S,E,V,Op,quartic=[D[k] for k in paths]
    hashes={p:sha(ROOT/p) for p in paths.values()}
    (out/'phase.json').write_text(json.dumps({'phase':'all_saved_source_data_loaded','input_sha256':hashes,'no_new_source_or_implicit_solve':True},indent=2)+'\n')
    # Stage1 native source values in C0=QQ[tA,tB]/(tA²,tB²).
    fam={f['family']:f for f in I['families']};HA=ident(5);HB=ident(5)
    for r,s in itertools.product(range(5),repeat=2):
        HA[r][s]=add(HA[r][s],scale(TA,F(fam['CA']['source_frame5_by5_by9_first'][7][r][s])))
        HB[r][s]=add(HB[r][s],scale(TB,F(fam['CB']['source_frame5_by5_by9_first'][8][r][s])))
    gauge=ident(5)
    for (i,j),row in zip(I['actual_group15_positions'],I['group15_first9']):gauge[i][j]=add(gauge[i][j],scale(TA,F(row[7])))
    fA=frame(mm(gauge,HA));assert fA==ident(5)
    fB=frame(mm(gauge,HB));flags=[fB[1][0],fB[1][3],fB[1][4],fB[0][3],fB[0][4]]
    point=list(map(F,P['actual_P0_physical150']));ga=[[F(0)]*5 for _ in range(5)]
    for (i,j),row in zip(I['actual_group15_positions'],I['group15_first9']):ga[i][j]+=F(row[7])
    velocity=[F(0)]*150
    for i,j in itertools.product(range(5),repeat=2):
        if ga[i][j]:velocity=[x+ga[i][j]*y for x,y in zip(velocity,action(point,i,j))]
    moved=[add(c(x),scale(TA,y)) for x,y in zip(point,velocity)]
    rows=P['zero_coordinate_rows36'];ar=P['CA_constant_left_inverse_rows19'];ai=constmat(P['CA_left_inverse19'])
    own=mm(ai,[[moved[rows[j]]] for j in ar]);own=[r[0] for r in own]
    center=list(map(F,P['actual_CA_source19_center']));free=P['CA_free_fiber9_coefficient_indices'];fixed=P['CA_fiber10_fixed_coefficient_indices']
    assert all(own[j]==c(center[j]) for j in fixed)
    v=[add(own[j],c(-center[j])) for j in free]
    actual=v+flags;assert len(actual)==14
    for j,p in enumerate(actual):
        assert p[:3]==(F(0),F(I['source14_first_by_physical9'][j][7]),F(I['source14_first_by_physical9'][j][8]))
    assert all(evalpoly(p,actual)==Z for p in Q['equations14variable_sparse'])
    kappa=[actual[j] for j in I['retained_source9_indices']]
    (out/'phase1-native-values.json').write_text(json.dumps({'actual14_C0':serial([actual])[0],'kappa9_C0':serial([kappa])[0],'gauge':serial(gauge),'CB_flag_frame':serial(fB),'source_equations17_all_zero':True},indent=2)+'\n')
    # R17 is the literal finite coefficient action of flag_frame(G*HB)*G.
    transform=mm(fB,gauge);itransform=mi(transform);targets=B['physical_target_rows12']
    R17=[[Z]*12 for _ in range(17)];ann=constmat(P['CB_annihilator17_by36'])
    for j,target in enumerate(targets):
        inp=[Z]*150;inp[target]=ONE
        for k,row in enumerate(E['eliminated_zero_rows4']):
            inp[row]=c(E['eliminated_four_coordinate_reconstruction'][k][rows.index(target)])
        movedcolumn=physical_push(inp,transform,itransform)
        val=mm(ann,[[movedcolumn[i]] for i in rows])
        for k in range(17):R17[k][j]=val[k][0]
    Dunit=add(ONE,neg(mul(kappa[1],kappa[1])));Punit=add(c(5),mul(kappa[1],kappa[1]))
    firstnames=['v0','v1','b','c','e'];firstindices=[0,1,10,11,13]
    firstroots=[]
    for name in firstnames:
        r=H['five_explicit_rational_roots'][name]
        firstroots.append(mul(evalpoly(r['numerator'],kappa),power(Dunit,-r['D_denominator_power'])))
    assert [actual[j] for j in firstindices]==firstroots
    firstder=[[evalpoly(derivative(p,j),actual) for j in firstindices] for p in Q['equations14variable_sparse']]
    ai5=mi(firstder[:5]);ba=mm(firstder[5:],ai5)
    R12=[[add(R17[i+5][j],neg(add(*(mul(ba[i][k],R17[k][j]) for k in range(5))))) for j in range(12)] for i in range(12)]
    for i,r in enumerate(H['all12_actual_remaining_and_divided_records']):
        unit=scale(power(Dunit,r['D_denominator_power']),F(r['constant_denominator'],r['primitive_content_unit']))
        R12[i]=[mul(unit,p) for p in R12[i]]
    # Source H12 -> regular8+zero4 adapter, without any d inversion.
    hpolys=[r['divided_primitive_polynomial'] for r in H['all12_actual_remaining_and_divided_records']]
    core=[kappa[j] for j in [0,1,2,6]];secondnames=K['selected5_pivot_order']
    secondindices=[K['source9_variable_order'].index(name) for name in secondnames]
    secondroots=[]
    for name in secondnames:
        r=K['five_actual_rational_roots'][name]
        secondroots.append(mul(evalpoly(r['numerator'],core),power(Dunit,-r['D_power']),power(Punit,-r['P_power'])))
    coefficients=[telescope(p,kappa,secondindices,secondroots) for p in hpolys]
    order=K['selected5_source_H_rows']+[j for j in range(12) if j not in K['selected5_source_H_rows']]
    assert order==Op['source_adapter_T']['permutation']
    A5=[coefficients[j] for j in order[:5]];B7=[coefficients[j] for j in order[5:]];Ainv=mi(A5);B7A=mm(B7,Ainv)
    T1=[[Z]*12 for _ in range(12)]
    for i in range(5):
        for j in range(5):T1[i][order[j]]=Ainv[i][j]
    for i in range(7):
        T1[i+5][order[i+5]]=ONE
        for j in range(5):T1[i+5][order[j]]=neg(B7A[i][j])
    T=T1
    for i,r in enumerate(K['remaining7_actual_core_polynomials']):
        lam=scale(mul(power(Dunit,-r['D_denominator_power']),power(Punit,-r['P_denominator_power'])),F(r['QQ_unit_content_numerator'],r['QQ_unit_content_denominator']))
        T[i+5]=[mul(inv(lam),p) for p in T[i+5]]
    E3=ident(12);E3[5][5]=inv(Dunit);E3[6][6]=inv(Dunit)
    rr=add(core[2],neg(core[0]),neg(core[1]));ss=add(core[0],core[2])
    E3[7][5]=neg(mul(rr,inv(mul(Punit,Dunit))));E3[7][7]=scale(inv(mul(Punit,Dunit)),-2)
    T=mm(E3,T);E4=ident(12)
    for j in range(4):
        factor=power(Dunit,[2,2,2,3][j])
        for i in range(3):E4[8+j][5+i]=neg(mul(factor,evalpoly(C['cofactor3_by4'][i][j],core)))
    T=mm(E4,T);N=mm(T,R12)
    gs=[add(kappa[j],neg(root)) for j,root in zip(secondindices,secondroots)]
    gs += [evalpoly(p,core) for p in prep['generators3']]
    hvalues=[[evalpoly(p,kappa)] for p in hpolys]
    assert mm(T,hvalues)==[[p] for p in gs]+[[Z] for _ in range(4)]
    assert all(mul(kappa[1],p)==Z for p in gs)
    # All native source adapters are units; certify their scalar augmentations.
    Ninv=mi(N);assert mm(N,Ninv)==ident(12)
    (out/'phase2-native-adapter.json').write_text(json.dumps({'R17':serial(R17),'R12':serial(R12),'T':serial(T),'N':serial(N),'N_inverse':serial(Ninv),'G8_C0':list(map(lambda p:list(map(str,p)),gs)),'all_source_H12_to_G8_identities':True,'d_not_inverted':True},separators=(',',':'))+'\n')
    Snative=[[Z]*32 for _ in range(12)]
    for i,j,g,sgn in Op['native_S12_by32']['koszul28_in_first8_rows']:Snative[i][j]=scale(gs[g],F(sgn))
    for i,j,value in Op['native_S12_by32']['free4_zero_source_columns']:Snative[i][j]=c(value)
    native=mm([list(row) for row in zip(*N)],Snative)
    syz=[[F(0)]*128 for _ in range(48)]
    for i,j in itertools.product(range(12),range(32)):
        for b in range(4):
            eb=tuple(F(k==b) for k in range(4));v=mul(native[i][j],eb)
            for k in range(4):syz[4*i+k][4*j+b]=v[k]
    # Constants of every chain scalar: reuse saved regular-function caches,
    # never evaluate a vanishing numerator divided by a vanishing denominator.
    constants={int(i):F(p[0]) for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA'] for i,p in F9[key].items()}
    constants.update({i:F(dict(p).get(0,'0')) for i,p in S['all_requested_original_source_root_values']})
    constants.update({i:F(0) for i in C0['actual_C0_basis8_scalar_nodes']})
    for r in C0['all35_regular_coordinate_divisions']:
        first=F9['both57_source2']['CA']['pure_function2jets'][str(r['original_pure_CA_source_function'])]
        for j,q in enumerate(r['all7_regular_coordinate_divisions']):
            assert F(q['origin'])==F(first[1+j]);constants[q['node']]=F(first[1+j])
    oldg=quartic['original_regular_function_nodes_unmodified']['old_CB'];newg=quartic['original_regular_function_nodes_unmodified']['new_CA']
    constants.update({i:F(0) for i in oldg+newg})
    # Both g plane functions begin at degree3.
    gb=next(r['gB'] for r in B['accepted_actual_regular_function_identities'] if r['kind']=='accepted_actual_CB_f63_tail')
    constants[gb]=F(0)
    full=C0['DAG_scalar_nodes'];memo=dict(constants)
    def original(i):
        if i in memo:return memo[i]
        n=full[i];tag=n[0]
        if tag=='c':v=F(n[1],n[2])
        elif tag=='a':v=original(n[1])+original(n[2])
        elif tag=='n':v=-original(n[1])
        elif tag=='m':
            a=original(n[1]);v=F(0) if not a else a*original(n[2])
        else:raise AssertionError(('unbound_original_regular_constant',i,n))
        memo[i]=v;return v
    ops=B['operator_scalar_nodes'];scalar=[]
    for n in ops:
        tag=n[0]
        if tag=='c':v=F(n[1],n[2])
        elif tag=='actual_C0_scalar_function':v=original(n[1])
        elif tag=='sum':v=sum((scalar[j]*F(a,b) for j,a,b in n[1]),F(0))
        elif tag=='product':
            v=F(n[1],n[2])
            for j in n[3]:v*=scalar[j]
        else:raise AssertionError(n)
        scalar.append(v)
    X=[[F(0)]*768 for _ in range(56)]
    for j,rec in enumerate(B['B48_by640_sparse_columns']):
        for i,k in rec:X[i][j]=scalar[k]
    assert C0['D0_relation_origin_rank']==0
    for i in range(48):X[i][640:]=syz[i]
    m0=[[F(x) for x in row] for row in V['all24_by56_m0']]
    for row in m0:
        for j in range(768):assert sum((row[i]*X[i][j] for i in range(56) if row[i] and X[i][j]),F(0))==0
    selected=V['fixed_selected24_columns'];unselected=[i for i in range(56) if i not in selected]
    residual=[X[i] for i in unselected];pivots={};basiscols=[];receipts=[];combinations=[]
    for j in range(768):
        v={i:residual[i][j] for i in range(32) if residual[i][j]};comb={j:F(1)};steps=[]
        while v and min(v) in pivots:
            k=min(v);q=v[k];p,pc=pivots[k];steps.append([k,str(q)])
            for i,a in p.items():v[i]=v.get(i,F(0))-q*a
            v={i:a for i,a in v.items() if a}
            for i,a in pc.items():comb[i]=comb.get(i,F(0))-q*a
            comb={i:a for i,a in comb.items() if a}
        rec={'source_column':j,'subtract':steps}
        if v:
            k=min(v);q=v[k];v={i:a/q for i,a in v.items()};comb={i:a/q for i,a in comb.items()};pivots[k]=(v,comb);basiscols.append(j)
            rec.update(pivot=k,divide=str(q),normalized=[[i,str(a)] for i,a in sorted(v.items())])
        else:rec['zero']=True
        rec['constant_original_source_combination']=[[i,str(a)] for i,a in sorted(comb.items())];receipts.append(rec)
        # Exact full56 readback of every transformed zero/pivot column.
        actual=[sum((X[i][k]*a for k,a in comb.items()),F(0)) for i in range(56)]
        if rec.get('zero'):assert not any(actual)
        else:assert [actual[i] for i in unselected]==[v.get(k,F(0)) for k in range(32)]
    rank=len(pivots);nonzero=sum(bool(x) for row in X for x in row)
    assert nonzero+sum(bool(x) for row in syz for x in row)<150000
    evidence={'C0_native_algebra':'QQ[tauA,tauB]/(tauA^2,tauB^2)',
      'full_X0_56_by768':[[str(x) for x in row] for row in X],
      'all48_by128_SyzQ_origins':[[str(x) for x in row] for row in syz],
      'all48_by640_chain_B_origins':[[str(x) for x in row[:640]] for row in X[:48]],
      'honest_P8_origins_zero':True,'all24_by768_m0_X0_zero':True,
      'free_kernel_rows32':unselected,'rank_X0':rank,'selected_independent_columns':basiscols,
      'all768_fixed_elimination_receipts_and_constant_combinations':receipts,
      'mu_T_lower_bound':32-rank,'T_length_lower_bound':32-rank,
      'mu_D_K_lower_bound':max(0,28-rank),'depth19_proved_from_X0_only':rank<=23,
      'actual_X1_not_evaluated':True,'native_unit_inverse_certificates':INVERSES,
      'source_function_constant_refs':[[i,str(original(i))] for i in sorted(n[1] for n in ops if n[0]=='actual_C0_scalar_function')]}
    ep=out/'actual-X0-and-unit-cancellation.json';ep.write_text(json.dumps(evidence,separators=(',',':'))+'\n')
    report={'id':STEM,'all_passed':True,'actual_X0_rank':rank,'mu_Rbar_T_lower_bound':32-rank,
      'mu_D_K_lower_bound':max(0,28-rank),'depth19_proved_from_X0_only':rank<=23,
      'all640_B_P_origins_and128_Res4_NtS_origins_bound':True,
      'no_new_source_implicit_solve_point_jet_GB_or_X1':True,
      'coefficient_nonzeros':nonzero,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['id','all_passed','actual_X0_rank','mu_Rbar_T_lower_bound','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
