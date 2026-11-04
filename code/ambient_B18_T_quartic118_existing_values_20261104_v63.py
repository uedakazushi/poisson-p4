#!/usr/bin/env python3
"""Restricted118 pure-quartic channels, using only frozen P3 and B11 X2.

All scalar roots are source-bound before value evaluation. No new source
coefficient, jet, normalizer inverse, or coefficient solve is performed.
"""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63'
def dump(out,name,x):
    (out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
def put(v,i,a):
    c=v.get(i,F(0))+a
    if c:v[i]=c
    else:v.pop(i,None)
def enc(v):return [[list(i) if isinstance(i,tuple) else i,str(a)] for i,a in sorted(v.items()) if a]
def main(out,evaluate,theorem):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'cases':'reports/AMBIENT_B18_T_QUARTIC118_NATIVE_FIRST_GATE_20261104_v63/actual118-old-native-first-separation.json',
      'theta2':'reports/AMBIENT_B18_T_CONSTANT_QUADRATIC_KERNEL_R3_ADJOINT_20261104_v63/phase1-theta-and-pure-r3-second-adjoint.json',
      'offset':'reports/AMBIENT_B18_T_QUADRATIC_L_SYZYGY_R4_20261104_v63/actual-quadratic-L-syzygy-scalar-offset-certificate.json',
      'rewrite':'reports/AMBIENT_B18_T_KNOWN_QUADRATIC_CUBIC_CLOSURE_20261104_v63/linear-module-rewrite-certificate.json',
      'gates':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
      'linear':'reports/AMBIENT_B18_COMPARISON_T_LINEAR_QUADRATIC_20261104_v63/phase1-actual-five-target-linear-relations.json',
      'mu':'reports/AMBIENT_B18_COMPARISON_T_CONSTANT_RANK_20261104_v63/actual-linear-cycle-muT-certificate.json',
      'pi':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
      'm0':'reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json',
      'source620':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
      'rank1':'reports/AMBIENT_B18_T_ACTUAL_PROJECTED_QUADRATIC_RANK_20261104_v63/actual-three-by620-quadratic-image-certificate.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'composite':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
      'P3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
      'source_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json'}
    if theorem:paths['restricted_theorem']=str(theorem)
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'frozen_sources_before_contraction','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()}
    C,T,O,NF,G,L,M,P,V,S,R,B,H,AD,CV,N,MS,P3,G2,SP=[D[k] for k in list(paths)[:20]]
    assert SP['all_passed'] and sha(ROOT/paths['source_peer'])=='808637882f3d8d28f393a7b0f4abfeec59e076bd1269034bd67d8ac6dd8c706a'
    mons2=list(itertools.combinations_with_replacement(range(8),2));mi2={m:i for i,m in enumerate(mons2)}
    mons3=list(itertools.combinations_with_replacement(range(8),3));mi3={m:i for i,m in enumerate(mons3)}
    assert NF['monomials3_120']==[list(m) for m in mons3]
    k={i:F(c) for i,lift,c in O['all600_chosen_source_lifts_and_scalar_offsets'] if F(c)}
    assert k=={239:-F(55495,20736),359:F(27785,1536)}
    nf0={i:rec(v)[0] for i,step,v in NF['all600_rewrite_steps'] if rec(v).get(0)}
    assert nf0=={359:F(1,4),479:F(1)}
    theta2=rec(T['theta_extension180']);adj2=rec(T['theta_plus_Ystar2_beta_quadratic_source_adjoint56_times36'])
    star=[F(r['known_rr']) for r in G['normalized_Ystar84_source55_rr']]
    qnorm=list(map(rec,L['normalized84_unitcycle_Q1_columns_sparse']))
    wb=list(map(rec,M['complete_rank_certificate']['basis_normalized_columns']));wl={min(v):j for j,v in enumerate(wb)};wc=[]
    for i in range(84):
        v={i:F(1)};c={}
        while v:
            q=min(v);j=wl[q];a=v[q];put(c,j,a);addscaled(v,wb[j],-a)
        wc.append(c)
    m0=[list(map(F,row)) for row in V['all24_by56_m0']]
    pirows=[{P['free_kernel_rows32'][j]:a for j,a in rec(p).items()} for p in P['pi_free32_rows']]
    def gamma(theta,include_star):
        beta={};z={}
        for p,q in enumerate(qnorm):
            for m,pair in enumerate(mons2):
                a=-sum((b*theta.get(120*(i//8)+mi3[tuple(sorted(pair+(i%8,)))],F(0)) for i,b in q.items()),F(0))
                if include_star and pair==(7,7):a+=star[p]
                if a:beta[36*p+m]=a
        for i,a in theta.items():
            t,m=divmod(i,120)
            for src,b in pirows[t].items():put(z,120*src+m,a*b)
        for w,(sector,pair) in enumerate(M['literal_wedge_basis84']):
            i,j=pair
            for m,pair2 in enumerate(mons2):
                a=sum((b*beta.get(36*p+m,F(0)) for p,b in wc[w].items()),F(0))
                if not a:continue
                terms=[(8*sector+7,i,1)] if j==7 else [(8*sector+j,i,1),(8*sector+i,j,-1)]
                for row,var,sgn in terms:
                    mon=tuple(sorted(pair2+(var,)));zd=sum(v<7 for v in mon);den=zd if j==7 else zd+1
                    c=a*F(sgn*mon.count(var),den)
                    for src,b in enumerate(m0[row]):
                        if b:put(z,120*src+mi3[mon],c*b)
        return z
    gam0=gamma(k,True);gamb=gamma(nf0,False)
    expected={120*(i//36)+mi3[mons2[i%36]+(7,)]:a for i,a in adj2.items()}
    assert gam0==expected and gamb=={120*29+119:F(1,4),120*30+119:F(1)}
    xp1=list(map(rec,P['actual_pivot_X1_by27_sparse56_times8']))
    ops=B['operator_scalar_nodes'];bcols=list(map(dict,B['B48_by640_sparse_columns']));hnode=H['DAG_scalar_nodes']
    def key(j,src,m):
        if j>=640:return (2,j-640,src,m) if src<48 else None
        if src<48:
            root=bcols[j].get(src)
            if root is None or ops[root][0]=='c':return None
            op=ops[root];return (1,op[1],m) if op[0]=='actual_C0_scalar_function' else (0,root,m)
        root=H['honest_D0_relation_columns'][j][src-48+1]
        return None if hnode[root][0]=='c' else (1,root,m)
    def c1_functional(c1,a2,gam):
        f={};known=F(0)
        for index,c in gam.items():
            src,m=divmod(index,120);mon=mons3[m]
            for v in set(mon):
                p=list(mon);p.remove(v);m2=mi2[tuple(p)]
                for j,a in c1[v].items():
                    label=key(j,src,m2)
                    if label is not None:put(f,label,-c*a)
                for p,first in enumerate(xp1):known+=c*first.get(8*src+v,F(0))*a2.get(36*p+m2,F(0))
        return f,known
    def rawP(c):
        z={}
        for j,a in c.items():
            if j<640:
                root=H['honest_D0_relation_columns'][j][8]
                if hnode[root][0]!='c':put(z,root,a)
        return z
    def xp_correction(a,rrstride):
        z={}
        for p in range(27):
            c=a.get(rrstride*p+(7 if rrstride==8 else 35),F(0))
            if c:
                for root,b in G['normalized_Xp27_source55_rr'][p]['rawP8_root_coefficients']:put(z,root,-c*F(b))
        return z
    # Canonical qB lift = selected620[563]/(-250/3) minus the N177 lift.
    pivot=563;scale=-F(250,3);s=S['all620'][pivot];c=rec(s['combination_in_actual768']);a1=rec(s['linear_Xp27_correction_row8j_plus_variable'])
    p3b=rawP(c);xpb=xp_correction(a1,8);fb={};kb=F(0)
    for index,a in adj2.items():
        src,m=divmod(index,36);u,v=mons2[m]
        for p,first in enumerate(xp1):
            kb+=a*a1.get(8*p+u,F(0))*first.get(8*src+v,F(0))
            if u!=v:kb+=a*a1.get(8*p+v,F(0))*first.get(8*src+u,F(0))
        for j,b in c.items():
            label=key(j,src,m)
            if label is not None:put(fb,label,-a*b)
    candidates=[];all_labels=set(fb);raw_roots=set(p3b);xp_roots=set(xpb);support=len(fb)+len(p3b)+len(xpb)+len(adj2)+len(gam0)+len(gamb)
    for row in C['all118']:
        c1=list(map(rec,row['linear_original768_combinations_by8variables']));a2=rec(row['known_normalized_Xp27_quadratic_a2'])
        p3=rawP(c1[7]);xp=xp_correction(a2,36);f0,k0=c1_functional(c1,a2,gam0);f1,k1=c1_functional(c1,a2,gamb)
        support+=sum(map(len,c1))+len(a2)+len(p3)+len(xp)+len(f0)+len(f1)+bool(k0)+bool(k1)
        assert support<150000
        all_labels.update(f0);all_labels.update(f1);raw_roots.update(p3);xp_roots.update(xp)
        candidates.append((row,p3,xp,f0,f1,k0,k1))
    allowed=[(0,7),(1,7),(4,7),(5,7),(7,7)];pairB=dict(zip(allowed,[7,8,9,10,2]));pairAD={p:j for j,p in enumerate(allowed)}
    requested={mons2[label[-1]] for label in all_labels};assert requested<=set(allowed),('UNCACHED_PAIR',requested-set(allowed))
    td={i:(rec(g),rec(q)) for i,v,g,f,q in AD['node_data']}
    cv={i:(rec(v),rec(g)) for i,v,g in CV['all108_composite_values_and_opaque_second_gradients']}
    native=[{row:rec(v) for row,v in col} for col in N['source_columns128_rows48_sparse_B11_effective_second']]
    typed={q[1] for q in all_labels if q[0]==0};direct={q[1] for q in all_labels if q[0]==1}
    assert typed<=td.keys() and direct<=cv.keys()
    atoms={i for root in typed for i in td[root][0]};assert atoms<=cv.keys()
    opaque={i:cv[i][1] for i in atoms|direct if cv[i][1]}
    assert not opaque,('UNCANCELLED_OPAQUE_X2',opaque)
    jets={i:list(map(F,p)) for i,p in P3['raw_source_roots_with_missing_quotient2_3_set_zero_ONLY_for_known_part']}
    jac={i:(rec(a),rec(b)) for i,a,b in MS['source_root_Jacobians0_and_r']}
    assert raw_roots<=jets.keys() and raw_roots<=jac.keys()
    assert all(not jac[i][0] and not jac[i][1] for i in raw_roots),('OPAQUE_P3_SENSITIVITY',[(i,jac[i]) for i in raw_roots if any(jac[i])])
    # Named quotient nodes with a literal zero value remain defined sources.
    # rec() is a sparse-vector reader and deliberately drops their keys.
    g2={i:F(a) for i,a in G2['actual_eight_regular_g_second_values']}
    xpvalues={i:sum((a*g2[j] for j,a in jac[i][0].items()),F(0)) for i in xp_roots}
    assert all(not jac[i][1] for i in xp_roots)
    for i in xp_roots:assert jets[i][2]==0,'known part of each gdiff root must be0'
    assert support+sum(len(v) for col in native for v in col.values())+sum(len(td[i][0]) for i in typed)+10000<150000
    pre={'actual118_count':len(candidates),'Gamma0_nonzeros':len(gam0),'Gammab_nonzeros':len(gamb),
      'Gamma0_equals_frozen31_adjoint_times_r':True,'Gammab_literal_source29_over4_plus_source30_pure_r3':True,
      'canonical_qB_pivot620':pivot,'canonical_qB_pivot_value':str(scale),'root_X2_labels':len(all_labels),
      'typed_B_roots':len(typed),'direct_C0_roots':len(direct),'rawP3_roots':len(raw_roots),'XP2_roots':len(xp_roots),
      'requested_pairs':sorted(requested),'all_X2_opaque_second_sensitivities_zero':True,
      'all_rawP3_opaque_second_and_third_sensitivities_zero':True,'actual_gdiff2_source5_already_frozen':True,
      'visible_QQ_support_bound':support+sum(len(v) for col in native for v in col.values())+sum(len(td[i][0]) for i in typed)+10000,
      'cap':150000,'no_new_source_coefficient_jet_normalizer_or_solve':True}
    dump(out,'preflight.json',pre)
    dump(out,'phase1-Gamma-and-b3-source.json',{'monomials2_36':mons2,'monomials3_120':mons3,
      'k600':encode(k),'NF0_600':encode(nf0),'Gamma0_56_times120':encode(gam0),'Gammab_56_times120':encode(gamb),
      'b3_source':{'selected620':pivot,'divide':str(scale),'rawP3':encode(p3b),'XP2':encode(xpb),'X2':enc(fb),'known':str(kb)},
      'Theta3_formula':'k600 + b3*NF0; b3 is corrected F620[563]/(-250/3)'})
    dump(out,'phase2-all118-source-functionals.json',{'all118':[{
      'kernel118_index':r['kernel118_index'],'combination_in119':r['combination_in119'],
      'rawP3':encode(p3),'XP2':encode(xp),'X2_constant':enc(f0),'X2_times_b3':enc(f1),
      'known_constant':str(k0),'known_times_b3':str(k1)} for r,p3,xp,f0,f1,k0,k1 in candidates],
      'formula':'rawP3 + XP2 + X2_constant + known_constant + b3*(X2_times_b3 + known_times_b3)',
      'native_effective_second_used_only_after_separate_native_zero_origin_combos':True})
    dump(out,'phase.json',{'phase':'SOURCE_PREFLIGHT_COMPLETE','all_required_values_cached':True})
    if not evaluate:
        dump(out,'report.json',{'id':STEM+'_PREFLIGHT','all_passed':True,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
          'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start});print(json.dumps(pre));return
    assert theorem and D['restricted_theorem']['all_passed'],'restricted-source theorem must be frozen before actual values'
    def x2val(label):
        kind=label[0];pair=mons2[label[-1]];b=pairB[pair]
        if kind==2:return native[label[1]].get(label[2],{}).get(b,F(0))
        if kind==1:return cv[label[1]][0].get(b,F(0))
        g,q=td[label[1]];return q.get(pairAD[pair],F(0))+sum((a*cv[i][0].get(b,F(0)) for i,a in g.items()),F(0))
    values={label:x2val(label) for label in all_labels}
    def val(p3,xp,f,k):return k+sum((a*jets[i][3] for i,a in p3.items()),F(0))+sum((a*xpvalues[i] for i,a in xp.items()),F(0))+sum((a*values[i] for i,a in f.items()),F(0))
    b3raw=val(p3b,xpb,fb,kb);b3=b3raw/scale
    actual=[]
    for row,p3,xp,f0,f1,k0,k1 in candidates:
        a=val(p3,xp,f0,k0);b=k1+sum((c*values[i] for i,c in f1.items()),F(0));actual.append(a+b3*b)
    nonzero=[i for i,a in enumerate(actual) if a]
    dump(out,'actual-b3-and-all118-pure-quartic-values.json',{'canonical_qB_pivot620':pivot,'F620_corrected_pure_cubic':str(b3raw),
      'canonical_qB_pure_cubic_b3':str(b3),'all_X2_scalar_root_pair_values':[[list(k),str(v)] for k,v in sorted(values.items())],
      'all118_actual_corrected_pure_quartic_values':list(map(str,actual)),'nonzero_indices':nonzero,
      'actual_T_total_length_if_any_nonzero':11 if nonzero else None,'actual_HF4_if_any_nonzero':0 if nonzero else None,
      'if_allzero_scope':'Only this118 linear-syzygy cubic-kernel family; not actualHF4 or totalT length',
      'no_new_source_coefficient_jet_or_solve':True})
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL','nonzero118_values':len(nonzero),'b3':str(b3)})
    report={'id':STEM,'all_passed':True,'actual118_nonzero_corrected_quartic_values':len(nonzero),'canonical_qB_pure_cubic':str(b3),
      'actual_T_total_length':11 if nonzero else None,'actual_HF4':0 if nonzero else None,
      'allzero_means_only_this_family_negative':not nonzero,'new_source_coefficients_computed':0,
      'independent_source_and_theorem_review_status':'PENDING','input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual118_nonzero_corrected_quartic_values','canonical_qB_pure_cubic','actual_T_total_length','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--evaluate',action='store_true');p.add_argument('--theorem',type=Path);a=p.parse_args();main(a.output,a.evaluate,a.theorem)
