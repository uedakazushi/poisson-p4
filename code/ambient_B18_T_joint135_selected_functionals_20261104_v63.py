#!/usr/bin/env python3
"""NEW selected third implicit functionals; no full105/gauge3 jet vector."""
import argparse,hashlib,itertools,json,re,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_comparison_X_origin_20261104_v63 as H
from ambient_B18_comparison_T_constant_rank_20261104_v63 import rec,encode,addscaled
from ambient_B18_native128_first_20261104_v63 import PATHS as OLD_PATHS,Algebra as Dual
from ambient_B18_T_joint135_functional_tools_20261104_v63 import MON,INDEX,Finite,Proper,source_evaluate,source_first_AD
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63'
THIRD=[(7,7,7),(0,7,7),(1,7,7),(4,7,7),(5,7,7)];BM=[INDEX[m] for m in THIRD]
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x
def nnz(x):
    # Sparse dictionary coefficients and dense tuple zeros are BOTH counted.
    if isinstance(x,F):return 1
    if isinstance(x,str) and re.fullmatch(r'-?\d+(?:/\d+)?',x):return 1
    if isinstance(x,dict):return sum(nnz(v) for v in x.values())
    if isinstance(x,(list,tuple)):return sum(nnz(v) for v in x)
    return 0
def rmul(v,rows):
    z={}
    for i,a in v.items():addscaled(z,rows[i],a)
    return z
def put(v,k,a):
    z=v.get(k,F(0))+a
    if z:v[k]=z
    else:v.pop(k,None)
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False);saved=0;maxcombined=0;filecounts={}
    paths={k:OLD_PATHS[k] for k in ['source','interface','proper','second5','preparation','Q_operator']}
    paths.update({
      'method':'reports/AMBIENT_B18_JOINT135_SELECTED_FUNCTIONAL_METHOD_PREFLIGHT_20261104_v63/report.json',
      'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
      'common3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'source20':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
      'gauge20':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase2-canonical270-zz-gauge.json',
      'proper20':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/phase3-actual-proper-G8-base20.json',
      'J30':'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/phase2-fixed30-Jacobian.json',
      'first_preflight':'reports/AMBIENT_B18_NATIVE128_FIRST_SOURCE_PREFLIGHT_20261104_v63/actual-native128-first-input-and-mixed-index-profile.json',
      'native0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/phase2-native-adapter.json',
      'native2':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
      'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
      'nativeS0':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
      'oldhigh':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/phase4-all619-high-old-source-sensitivities.json'})
    paths['MS_no_extra_requests']='reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_HIGHER_SOURCE_PREFLIGHT_20261104_v63/report.json'
    def save(n,x,live=()):
        nonlocal saved,maxcombined
        count=nnz(x);old=filecounts.get(n,0);maxcombined=max(maxcombined,saved-old+count+nnz(live));assert maxcombined<150000,('QQ_GUARD',n,maxcombined)
        p=out/n;p.write_text(json.dumps(enc(x),separators=(',',':'))+'\n');saved+=count-old;filecounts[n]=count;assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
        (out/'phase.json').write_text(json.dumps({'phase':n,'seconds':time.monotonic()-start,'saved_QQ':saved,'max_saved_live_QQ':maxcombined})+'\n');print(json.dumps({'phase':n,'seconds':time.monotonic()-start}),flush=True)
    hashes={p:H.sha(ROOT/p) for p in paths.values()};save('source-manifest.json',{'input_sha256':hashes,'source_helpers_sha256':H.sha(Path(__file__).with_name('ambient_B18_T_joint135_functional_tools_20261104_v63.py'))})
    D={k:H.read(ROOT/p) for k,p in paths.items()};assert D['method']['all_passed'] and D['MS_no_extra_requests']['all_passed']
    S,I,P,K,prep,OP=[D[k] for k in ['source','interface','proper','second5','preparation','Q_operator']]
    J=D['joint'];unknown=J['all105_unknown_original_node_order'];uindex={i:j for j,i in enumerate(unknown)};eq=J['all105_equation_original_node_order'];nodes=S['DAG_scalar_nodes']
    c20={i:rec(v) for i,v in D['source20']['all105_actual_source2_values']};phi0={i:p.get(0,F(0)) for i,p in c20.items()};assert set(phi0)==set(unknown)
    roots=[i for f in S['families'].values() for p in f['complete_normalized_physical_w_coefficients'][0] for i in p]+[f['actual_tracefree_nu'] for f in S['families'].values()]
    sourceAD,adpeak=source_first_AD(nodes,roots,phi0)
    pairnodes={f:q['complete_normalized_physical_w_coefficients'][0] for f,q in S['families'].items()};nunodes=[f['actual_tracefree_nu'] for f in S['families'].values()]
    ev0=[F(p[0]) for p in D['gauge20']['actual_even15_base20']];od0=[F(p[0]) for p in D['gauge20']['actual_odd15_base20']]
    proper=Proper(S,I,P,K,prep,D['first_preflight']['GL18_source_parameter_positions']);Gbase=[tuple(map(F,p[:4])) for p in D['proper20']['actualG8_C80']]
    jc=[rec([[j,c] for j,c in enumerate(row)]) for row in J['fixed105_Jacobian']];jci=[rec([[j,c] for j,c in enumerate(row)]) for row in J['fixed105_Jacobian_inverse']]
    jg_expected=[rec([[j,c] for j,c in enumerate(row)]) for row in D['J30']['J30_inverse']]
    # Constant source/gauge Jacobians: nil monics vary with phi, never freeze them.
    Hg=[{} for _ in range(30)];Jg=[{} for _ in range(30)];Dt=[{} for _ in range(34)]
    save('preflight.json',{'kind':'NEW SELECTED THIRD IMPLICIT FUNCTIONAL RESEARCH','unknown_order':['original105','gauge_even15','gauge_odd15'],
      'requested_monomials':THIRD,'full_general_source3_or_gauge3_vector_permitted':False,
      'constant_Jacobian_and_target_map_dense_ceiling':135*34+30*135+2*105*105+2*30*30,
      'RHS3_only_known_base2_and270gauge':True,'source_AD_peak_QQ':adpeak,'cap':150000,
      'higherN3_drops_by_separate_nativeS0_c0_zero':True,'m1_m2_removed_by_degree3_start_sourcebridge':True})
    for direction in range(135):
        leaf=unknown[direction] if direction<105 else None
        nud=[sourceAD[n][1].get(leaf,F(0)) if leaf is not None else F(0) for n in nunodes]
        a=Dual(1,2,[[F(0),d] for d in nud]);phi={i:(phi0[i],F(i==leaf)) for i in unknown}
        pairs={f:[[(sourceAD[x][0],sourceAD[x][1].get(leaf,F(0)) if leaf is not None else F(0)),(sourceAD[y][0],sourceAD[y][1].get(leaf,F(0)) if leaf is not None else F(0))] for x,y in ps] for f,ps in pairnodes.items()}
        ev=[(x,F(direction==105+i)) for i,x in enumerate(ev0)];od=[(x,F(direction==120+i)) for i,x in enumerate(od0)]
        con,gs,kap=proper.evaluate(a,phi,pairs,ev,od)
        assert all(q[:4]==(F(0),)*4 and not q[6] and not q[7] for q in con)
        assert [tuple(q[:4]) for q in gs]==Gbase
        col=[q[4] for q in con]+[q[5] for q in con]
        for i,c in enumerate(col):
            if c:(Hg if direction<105 else Jg)[i][direction if direction<105 else direction-105]=c
        for i,q in enumerate(gs):
            for b in range(4):
                if q[4+b]:Dt[4*i+b][direction]=q[4+b]
        for i,c in enumerate(nud):
            if c:Dt[32+i][direction]=c
        if direction%15==14:save('phase1-Jacobian-checkpoint.json',{'completed_input_columns':direction+1,'Jg_phi_rows':list(map(encode,Hg)),'Jg_gamma_rows':list(map(encode,Jg)),'target_G_nu_rows':list(map(encode,Dt))})
    # The checkpoint is retained; accounting replaces its prior contents.
    assert all(rmul(row,Jg)=={i:F(1)} for i,row in enumerate(jg_expected))
    assert all(rmul(row,jc)=={i:F(1)} for i,row in enumerate(jci))
    save('phase1-actual-J135-blocks-and-target-Jacobian.json',{'physical105_order':unknown,'canonical_gauge30_order':['even0..14','odd0..14'],
      'Jg_phi30_by105':list(map(encode,Hg)),'Jg_gamma30':list(map(encode,Jg)),'Jg_gamma_inverse':list(map(encode,jg_expected)),
      'properG32_and_nu2_target_Jacobian34_by135':list(map(encode,Dt)),
      'all_base_G8_and_gauge_constraints_checked':True,'source_dependent_monic_carry_included':True,'J135_formula':'[[J105,0],[Jg_phi,J30]]'})
    # Known cubic residuals with ALL unknown cubic source/gauge coefficients zero.
    a0=Finite();bind={i:tuple([c20[i].get(j,F(0)) for j in range(20)]+[F(0)]*5) for i in unknown}
    vals,closure,peak=source_evaluate(nodes,list(dict.fromkeys(roots+eq)),bind,a0)
    pairs={f:[[vals[x],vals[y]] for x,y in ps] for f,ps in pairnodes.items()};nu=[vals[i] for i in nunodes]
    ev=[tuple(map(F,p))+(F(0),)*5 for p in D['gauge20']['actual_even15_base20']];od=[tuple(map(F,p))+(F(0),)*5 for p in D['gauge20']['actual_odd15_base20']]
    a=Finite(2,nu)
    dense_transient_bound=10*25*100+3*150*100+5000
    assert saved+nnz([a.table,Hg,Jg,Dt,jc,jci,jg_expected,bind,pairs,ev,od])+dense_transient_bound<150000
    save('phase2-whole-live-gate.json',{'base25_monomials':MON,'third_order_by_tuple':THIRD,
      'method_ring_permutation_bound_by_monomial_tuples':True,'transient_matrix_tensor_QQ_bound':dense_transient_bound,
      'saved_plus_fixed_live_QQ':saved+nnz([a.table,Hg,Jg,Dt,jc,jci,jg_expected,bind,pairs,ev,od]),
      'whole_less_than150000':True,'source_specific_off_source_extension_same_for_J_and_RHS':True,
      'no_Tnil_highest_callbacks_requested_MS_receipt_pinned':True})
    con,gs,kap=proper.evaluate(a,bind,pairs,ev,od)
    assert all(not q[4*b+t] for q in con for b in range(20) for t in range(4))
    assert all(not q[4*b+t] for q in con for b in BM for t in [2,3])
    rhsE=[[vals[i][b] for i in eq] for b in BM];rhsG=[[q[4*b] for q in con]+[q[4*b+1] for q in con] for b in BM]
    for j,i in enumerate(eq):
        wanted=[F(0)]*25
        if 48<=j<57:
            v=j-48
            if v in [0,1,4,5]:wanted[INDEX[(v,)]]=F(1)
            elif v>=7:wanted[INDEX[(7,)]]=F(1)
        assert list(vals[i][:20])==wanted[:20]
    known3={i:list(map(F,p))[3] for i,p in D['common3']['actual_common105_rline_through3']}
    assert all(rhsE[0][i]+sum((c*known3[unknown[j]] for j,c in row.items()),F(0))==0 for i,row in enumerate(jc))
    save('phase2-known-third-RHS-and-proper-response.json',{'monomials5':THIRD,'source105_RHS3_by_monomial':rhsE,'gauge30_RHS3_by_monomial':rhsG,
      'known_properG8_C100':gs,'known_nu2_base25':nu,'all_base20_source_and_gauge_constraints_verified':True,
      'all105_pure_r3_cached_source_crosschecks':True,'gauge30_r3_not_claimed_previously_cached':True,'source_DAG_closure':closure},live=[Hg,Jg,Dt,bind,pairs])
    del vals,kap,pairs,bind,ev,od,con
    # Build the requested native derivative covectors through32G+2nu only.
    high=D['high'];mons3=list(itertools.combinations_with_replacement(range(8),3));mi3={m:i for i,m in enumerate(mons3)}
    gamma=rec(high['Gamma4_56_times120']);grows=[{src:a for index,a in gamma.items() for src,mm in [divmod(index,120)] if mm==mi3[m] and src<48} for m in THIRD]
    native0=D['native0'];N0=[[[F(c) for c in p] for p in row] for row in native0['N']]
    source0=[{row:F(dict(p).get(0,'0')) for row,p in col if F(dict(p).get(0,'0'))} for col in D['nativeS0']['actual_S_source_columns128_rows48_sparse_B11']]
    caseN=[]
    for row in high['all619']:
        c={j-640:v for j,v in rec(row['actual_original768_constant_combination']).items() if j>=640};z={}
        for j,v in c.items():addscaled(z,source0[j],v)
        assert not z;caseN.append(c)
    used=sorted({j for c in caseN for j in c});prim=[[{} for j in range(128)] for m in THIRD]
    for var in range(34):
        nud=[F(var==32),F(var==33)];ad=Dual(1,2,[[F(0),v] for v in nud]);ad.install()
        gg=[]
        for i,p in enumerate(Gbase):gg.append(tuple(p)+(tuple(F(var==4*i+b) for b in range(4))))
        sn=[[ad.Z]*32 for _ in range(12)]
        for i,j,k,sign in OP['native_S12_by32']['koszul28_in_first8_rows']:sn[i][j]=ad.scale(gg[k],F(sign))
        for i,j,c in OP['native_S12_by32']['free4_zero_source_columns']:sn[i][j]=ad.c(c)
        for j in used:
            jn,b=divmod(j,4);unit=[F(0)]*8;unit[b]=F(1);values={}
            for outrow in range(12):
                q=ad.Z
                for k in range(12):
                    if sn[k][jn]!=ad.Z:q=ad.add(q,ad.mul(tuple(N0[k][outrow])+(F(0),)*4,sn[k][jn]))
                q=ad.mul(q,tuple(unit))
                for t in range(4):
                    if q[4+t]:values[4*outrow+t]=q[4+t]
            for m,row in enumerate(grows):
                v=sum((c*values.get(i,F(0)) for i,c in row.items()),F(0))
                if v:prim[m][j][var]=v
    save('phase3-native-requested-functional-factors.json',{'target_order':'G0mask0..3,...,G7mask0..3,nuA,nuB',
      'native_source_derivative_by_monomial_and_column':[[encode(p) for p in row] for row in prim],
      'all619_native_S0_c0_exact_zero':True,'N3_component_omitted_only_after_this_zero':True,
      'monic_nu3_scalar_reduction_derivative_included':True,'m1_m2_absent_by_degree3_start_619_bridge':True},live=[Dt,Hg,Jg,gs,nu])
    # Stream sparse block adjoints; do not construct any source/gauge cubic vector.
    oldhigh=D['oldhigh']['all619'];requests=[];allvalues=[];support=0;nonzero=0
    for number,(case,orow) in enumerate(zip(caseN,oldhigh)):
        rowvalues=[]
        for mon,m in enumerate(THIRD):
            alpha={}
            for j,c in case.items():addscaled(alpha,prim[mon][j],-c)
            v=rmul(alpha,Dt)
            oldphi={}
            for (kind,leaf,mm),c in orow['old_high_source3_sensitivity']:
                if kind=='physical' and mm==mi3[m]:put(oldphi,uindex[leaf],F(c))
            addscaled(v,oldphi,F(1));vp={i:c for i,c in v.items() if i<105};vg={i-105:c for i,c in v.items() if i>=105}
            rg=rmul(vg,jg_expected);res=dict(vp);addscaled(res,rmul(rg,Hg),F(-1));rp=rmul(res,jci)
            chk=rmul(rp,jc);addscaled(chk,rmul(rg,Hg),F(1));assert chk==vp and rmul(rg,Jg)==vg
            value=-sum((c*rhsE[mon][i] for i,c in rp.items()),F(0))-sum((c*rhsG[mon][i] for i,c in rg.items()),F(0))
            rowvalues.append(value)
            if v:
                nonzero+=1;support+=len(alpha)+len(oldphi)+len(rp)+len(rg);assert saved+support+nnz([Dt,Hg,Jg,prim,gs,nu])+10000<150000
                requests.append({'kernel619_index':number,'monomial_original8':m,'native_gradient34':encode(alpha),
                  'old_physical_gradient105':encode(oldphi),'rho_phi105':encode(rp),'rho_gauge30':encode(rg),
                  'all135_rhoJ_equals_requested_gradient':True,'actual_implicit_third_functional_value':str(value)})
        allvalues.append(rowvalues)
    save('phase4-selected-joint135-third-adjoint-values.json',{'all_nonzero_requested_adjoints':requests,
      'all619_actual_implicit_third_functional_values_by5':allvalues,'all135_components_checked_for_each_adjoint':True,
      'rho_formula':'rho_gamma=v_gamma*J30inv;rho_phi=(v_phi-rho_gamma*H)*J105inv',
      'no_full_source105_or_gauge30_cubic_vector_solved_or_saved':True,
      'these_are_requested_implicit_contributions_not_complete619_quartic_values':True},live=[Dt,Hg,Jg,prim])
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'kind':'NEW SELECTED JOINT135 THIRD IMPLICIT FUNCTIONALS','actual_nonzero_requested_adjoints':nonzero,
      'full_general_source105_gauge30_cubic_jet_constructed':False,'all_adjoint_unit_source_constraints_verified':True,
      'all105_pure_r3_cached_crosschecks':True,'Tnil_specific_highest_callbacks_added':False,
      'actual619_complete_quartic_values_HF4_or_length_computed':False,'source6_or_unrequested_source4_computed':False,
      'max_saved_live_QQ':maxcombined,'input_sha256':hashes,'code_sha256':H.sha(Path(__file__)),
      'helper_code_sha256':H.sha(Path(__file__).with_name('ambient_B18_T_joint135_functional_tools_20261104_v63.py')),
      'note_sha256':H.sha(note),'evidence_sha256':{p.name:H.sha(p) for p in out.iterdir() if p.name!='phase.json'},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'independent_review':'PENDING'}
    save('report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual_nonzero_requested_adjoints','max_saved_live_QQ','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
