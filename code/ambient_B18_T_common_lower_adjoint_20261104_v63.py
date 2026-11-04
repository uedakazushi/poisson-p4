#!/usr/bin/env python3
"""One common D56x36 for619 lower quartic terms from existing caches."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_COMMON_LOWER_ADJOINT_20261104_v63'
def put(v,k,a):
    c=v.get(k,F(0))+a
    if c:v[k]=c
    else:v.pop(k,None)
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def main(out):
    st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    def dump(n,x):
        (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths={
      'preflight':'reports/AMBIENT_B18_T_COMMON_LOWER_ADJOINT_SOURCE_PREFLIGHT_20261104_v63/report.json',
      'lambda':'reports/AMBIENT_B18_T_COMMON_LOWER_ADJOINT_SOURCE_PREFLIGHT_20261104_v63/all121-common-lower-lambda-source-functionals.json',
      'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'composite':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'P3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
      'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
      'gates':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
      'X0':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
      'pi':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
      'star':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json',
      'L37':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'low':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/actual-linear-quadratic-T-certificate.json',
      'syzygies':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_PURE_R3_20261104_v63/actual-linear-syzygies-and-pure-r3-source-values.json',
      'mu':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'm0':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
      'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};dump('phase.json',{'phase':'frozen_common_lower_inputs_before_cache_values','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};PF,S,AD,CV,N,P3,MS,G2,G,X,P,ST,A,L,LOW,Y,M,V,HI,B,C=[D[k] for k in paths]
    assert PF['all_passed'] and PF['all_required_values_already_fixed']
    mons=list(itertools.combinations_with_replacement(range(8),2));mi={m:i for i,m in enumerate(mons)}
    mons3=list(itertools.combinations_with_replacement(range(8),3));gam=rec(S['Gamma4_56_times120'])
    allowed=[(0,7),(1,7),(4,7),(5,7),(7,7)];pB=dict(zip(allowed,[7,8,9,10,2]));pAD={p:j for j,p in enumerate(allowed)}
    td={i:(rec(g),rec(q)) for i,v,g,f,q in AD['node_data']};cv={i:(rec(v),rec(g)) for i,v,g in CV['all108_composite_values_and_opaque_second_gradients']}
    nv=[{row:rec(v) for row,v in col} for col in N['source_columns128_rows48_sparse_B11_effective_second']]
    jets={i:list(map(F,p)) for i,p in P3['raw_source_roots_with_missing_quotient2_3_set_zero_ONLY_for_known_part']}
    jac={i:(rec(a),rec(b)) for i,a,b in MS['source_root_Jacobians0_and_r']};g2={i:F(a) for i,a in G2['actual_eight_regular_g_second_values']}
    def value(label):
        kind=label[0];pair=mons[label[-1]];b=pB[pair]
        if kind==2:return nv[label[1]].get(label[2],{}).get(b,F(0))
        if kind==1:assert not cv[label[1]][1];return cv[label[1]][0].get(b,F(0))
        g,q=td[label[1]];assert all(not cv[i][1] for i in g)
        return q.get(pAD[pair],F(0))+sum((a*cv[i][0].get(b,F(0)) for i,a in g.items()),F(0))
    labels={tuple(k) for r in S['all121'] for f in r['lambda_by8variables_X2_functionals'] for k,c in f}
    bound=PF['source_QQ_support']+len(labels)+sum(len(v) for col in nv for v in col.values())+sum(len(c) for c in LOW['quadratic_certificate']['basis_normalized_columns'])+30000
    assert bound<150000;dump('preflight.json',{'whole_forecast_QQ_support':bound,'cap':150000,'all2025_shared_labels_fixed':True,'no_new_source_coefficient_or_solve':True})
    vals={k:value(k) for k in labels};lb={};ll={}
    for row in S['all121']:
        z={}
        for var,f in enumerate(row['lambda_by8variables_X2_functionals']):
            v=sum((F(a)*vals[tuple(k)] for k,a in f),F(0))
            if var==7:
                for i,a in row['additional_rawP3_functional_at_variable_r']:
                    assert not jac[i][0] and not jac[i][1];v+=F(a)*jets[i][3]
            if v:z[var]=v
        target=lb if row['shared_source_kind']=='star' else ll
        for var,v in z.items():target[8*row['index']+var]=v
    xp1=list(map(rec,P['actual_pivot_X1_by27_sparse56_times8']));xp0=list(map(rec,P['actual_pivot_X0_columns27_sparse56']));la={}
    for index,a in gam.items():
        src,m=divmod(index,120);mon=mons3[m]
        for var in set(mon):
            pair=list(mon);pair.remove(var);m2=mi[tuple(pair)]
            for p,first in enumerate(xp1):put(la,36*p+m2,a*first.get(8*src+var,F(0)))
    xprr=[]
    for p,record in enumerate(G['normalized_Xp27_source55_rr']):
        v=F(0)
        for root,a in record['rawP8_root_coefficients']:
            v+=F(a)*(jets[root][2]+sum((b*g2[i] for i,b in jac[root][0].items()),F(0)))
        xprr.append(v);put(la,36*p+35,-v)
    free=X['free_kernel_rows32'];receipts=[r for r in X['all768_fixed_elimination_receipts_and_constant_combinations'] if 'pivot' in r]
    pb={r['pivot']:rec(r['normalized']) for r in receipts};pindex={r['pivot']:i for i,r in enumerate(receipts)};zeta=[{} for p in range(27)]
    for j,src in enumerate(free):
        v={j:F(1)}
        for p in sorted(pb):
            a=v.get(p,F(0))
            if a:put(zeta[pindex[p]],src,a);addscaled(v,pb[p],-a)
    assert all(sum((a*xp0[j].get(src,F(0)) for src,a in zeta[i].items()),F(0))==F(i==j) for i in range(27) for j in range(27))
    delta={}
    for index,a in la.items():
        p,m=divmod(index,36)
        for src,b in zeta[p].items():put(delta,36*src+m,a*b)
    star1=list(map(rec,ST['actual_unitstar84_first_sparse56_times8']));lin1=list(map(rec,A['actual_selected37_X1_sparse56_by8']))
    def subtract_first(lam,firsts,dual):
        r=dict(lam)
        for p,first in enumerate(firsts):
            for var in range(8):
                a=sum((b*dual.get(36*(i//8)+mi[tuple(sorted((var,i%8)))],F(0)) for i,b in first.items()),F(0))
                put(r,8*p+var,-a)
        return r
    lb1=subtract_first(lb,star1,delta);ll1=subtract_first(ll,lin1,delta)
    syz=[rec(r['coefficients296_sparse']) for r in Y['all119_linear_coefficient_syzygies']]
    checks=[sum((a*ll1.get(i,F(0)) for i,a in c.items()),F(0)) for c in syz]
    dump('phase1-actual-common-lambda-values.json',{'lambda_a2_27_times36':encode(la),'lambda_b1_84_times8':encode(lb),'lambda_l1_37_times8':encode(ll),
      'zeta27_by56_rows':list(map(encode,zeta)),'delta_from_a2_source56_times36':encode(delta),
      'lambda_b1_after_a2':encode(lb1),'lambda_l1_after_a2':encode(ll1),'all119_kernel_checks':list(map(str,checks)),
      'all_shared_X2_used_after_separate_native_zero_origin_combinations':True})
    assert not any(checks),'common lambda_l is not well-defined on N177: sign/source gate failed'
    norm=[]
    for j,r in enumerate(LOW['quadratic_certificate']['fixed_elimination_receipts'][:296]):
        c={j:F(1)}
        for p,a in r['subtract']:addscaled(c,norm[p],-F(a))
        if 'pivot_row' in r:norm.append({i:a/F(r['divide']) for i,a in c.items()})
    assert len(norm)==177
    basis=list(map(rec,LOW['quadratic_certificate']['basis_normalized_columns']));lead={min(v):j for j,v in enumerate(basis)}
    mus=[sum((a*ll1.get(i,F(0)) for i,a in c.items()),F(0)) for c in norm];mu={}
    for i in range(180):
        v={i:F(1)};a=F(0)
        for p in sorted(lead):
            if p in v:
                c=v[p];a+=c*mus[lead[p]];addscaled(v,basis[lead[p]],-c)
        if a:mu[i]=a
    assert all(sum((a*mu.get(i,F(0)) for i,a in b.items()),F(0))==v for b,v in zip(basis,mus))
    qnorm=list(map(rec,L['normalized84_unitcycle_Q1_columns_sparse']));lb2=subtract_first(lb1,qnorm,mu)
    pirows=[{P['free_kernel_rows32'][j]:a for j,a in rec(p).items()} for p in P['pi_free32_rows']]
    cov=dict(delta)
    for index,a in mu.items():
        t,m=divmod(index,36)
        for src,b in pirows[t].items():put(cov,36*src+m,a*b)
    wb=list(map(rec,M['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):j for j,v in enumerate(wb)};wc=[]
    for w in range(84):
        v={w:F(1)};c={}
        while v:
            i=min(v);p=wl[i];a=v[i];put(c,p,a);addscaled(v,wb[p],-a)
        wc.append(c)
    kap={}
    for w,(sector,pair) in enumerate(M['literal_wedge_basis84']):
        i,j=pair
        for var in range(8):
            a=sum((b*lb2.get(8*p+var,F(0)) for p,b in wc[w].items()),F(0))
            if not a:continue
            terms=[(8*sector+7,i,1)] if j==7 else [(8*sector+j,i,1),(8*sector+i,j,-1)]
            for row,k,sgn in terms:
                mon=tuple(sorted((var,k)));zd=sum(u<7 for u in mon);den=zd if j==7 else zd+1
                put(kap,36*row+mi[mon],a*F(sgn*mon.count(k),den))
    m0=[list(map(F,r)) for r in V['all24_by56_m0']]
    for index,a in kap.items():
        row,m=divmod(index,36)
        for src,b in enumerate(m0[row]):
            if b:put(cov,36*src+m,a*b)
    needed=sorted({mons[i%36] for i in cov})
    dump('phase2-common-D56-by36.json',{'D_source56_times36':encode(cov),'mu177_extension180':encode(mu),
      'lambda_b1_after_a2_and_l1':encode(lb2),'m0_kappa2_adjoint24_times36':encode(kap),'required_second_pairs':needed,
      'all119_lambda_l_kernel_zero':True,'fixed177_lift_not_new_rank':True,
      'formula':'lower619 = D(C2), C2=Xc2-Xp1*a1; all619 share D','no_new_source_coefficient':True})
    # Fold D against actual619 C2. Evaluate only if EVERY dependency fixed.
    ops=B['operator_scalar_nodes'];bc=list(map(dict,B['B48_by640_sparse_columns']));nodes=C['DAG_scalar_nodes']
    def key(j,src,m):
        if j>=640:return (2,j-640,src,m) if src<48 else None
        if src<48:
            root=bc[j].get(src)
            if root is None or ops[root][0]=='c':return None
            op=ops[root];return (1,op[1],m) if op[0]=='actual_C0_scalar_function' else (0,root,m)
        root=C['honest_D0_relation_columns'][j][src-48+1]
        return None if nodes[root][0]=='c' else (1,root,m)
    fs=[];allkeys=set();support=PF['source_QQ_support']+len(cov)+len(delta)+len(mu)+len(lb)+len(ll)+len(la)+len(labels)
    for r in HI['all619']:
        c=rec(r['actual_original768_constant_combination']);a1=rec(r['actual_Xp27_linear_a1']);f={};known=F(0)
        for index,a in cov.items():
            src,m=divmod(index,36);u,v=mons[m]
            for j,b in c.items():
                label=key(j,src,m)
                if label is not None:put(f,label,a*b)
            for p,first in enumerate(xp1):
                known-=a*a1.get(8*p+u,F(0))*first.get(8*src+v,F(0))
                if u!=v:known-=a*a1.get(8*p+v,F(0))*first.get(8*src+u,F(0))
        fs.append((f,known));allkeys.update(f);support+=len(f)+bool(known);assert support<150000
    missingpairs=sorted(set(needed)-set(allowed));typed={k[1] for k in allkeys if k[0]==0};direct={k[1] for k in allkeys if k[0]==1}
    missingtyped=typed-td.keys();atoms={i for r in typed&td.keys() for i in td[r][0]};missingdirect=(direct|atoms)-cv.keys();opaque={i:cv[i][1] for i in (direct|atoms)&cv.keys() if cv[i][1]}
    ready=not any([missingpairs,missingtyped,missingdirect,opaque]);actual=[]
    if ready:
        lookup={k:value(k) for k in allkeys}
        actual=[known+sum((a*lookup[k] for k,a in f.items()),F(0)) for f,known in fs]
    dump('phase3-all619-common-D-contractions.json',{'all619_source_functionals':[{'kernel619_index':i,'X2_source_functional':enc(f),'known_Xp1_a1':str(k)} for i,(f,k) in enumerate(fs)],
      'all_required_values_fixed':ready,'actual619_lower_contribution_values':list(map(str,actual)) if ready else None,
      'missing_pairs':missingpairs,'missing_typed_roots':sorted(missingtyped),'missing_direct_roots':sorted(missingdirect),
      'opaque_sensitivities_retained':[[i,encode(p)] for i,p in sorted(opaque.items())],
      'not_complete619_quartic_values':True,'no_new_source_jet_or_solve':True})
    dump('phase.json',{'phase':'FINAL','D_nonzeros':len(cov),'all619_lower_fixed':ready})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'common_D_source_nonzeros':len(cov),'required_second_pairs':needed,
      'all119_internal_kernel_checks_zero':True,'all619_lower_contributions_fixed':ready,
      'actual619_lower_nonzero_count':sum(bool(a) for a in actual) if ready else None,'complete_quartic_values_or_HF4_computed':False,
      'new_source_jet_coefficients_or_solve':False,'saved_source_QQ_support':support,'input_sha256':hs,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump('report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','common_D_source_nonzeros','required_second_pairs','all619_lower_contributions_fixed','actual619_lower_nonzero_count','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
