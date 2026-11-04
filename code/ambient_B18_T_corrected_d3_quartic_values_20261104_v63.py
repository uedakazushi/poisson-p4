#!/usr/bin/env python3
"""Canonical c4 and168 corrected d3 pure quartics from frozen caches only."""
import argparse,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_CORRECTED_D3_QUARTIC_VALUES_20261104_v63'
def dump(out,name,x):
    (out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
def put(v,i,a):
    c=v.get(i,F(0))+a
    if c:v[i]=c
    else:v.pop(i,None)
def enc(v):return [[list(i) if isinstance(i,tuple) else i,str(a)] for i,a in sorted(v.items()) if a]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'theorem':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/report.json',
      'source168':'reports/AMBIENT_B18_T_CORRECTED_D3_QUARTIC_PREFLIGHT_20261104_v63/actual168-corrected-d3-quartic-source-preflight.json',
      'gamma':'reports/AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63/phase1-Gamma-and-b3-source.json',
      'gamma_peer':'reports/AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63_REVIEW/report.json',
      'b3':'reports/AMBIENT_B18_T_QUARTIC118_EXISTING_VALUES_20261104_v63/actual-b3-and-all118-pure-quartic-values.json',
      'PiC':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase1-cubic-Pi2-and-source-adjoint.json',
      'c1_119':'reports/AMBIENT_B18_T_LINEAR_SYZYGY_CUBIC_ADJOINT_20261104_v63/phase2-linear-syzygy-source-root-contractions.json',
      'a2_119':'reports/AMBIENT_B18_T_CUBIC_X2_ONLY_PREFLIGHT_20261104_v63/actual119-cubic-X2-only-source-operator-preflight.json',
      'gates':'reports/AMBIENT_B18_T_CUBIC_BOUNDARY_SOURCE3_PREFLIGHT_20261104_v63/phase2-known-Schur-second-r-gates.json',
      'pi':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'composite':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json',
      'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
      'P3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
      'source_peer':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63_REVIEW/report.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'sources_pinned_before168_contractions','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};TH,S,GA,GP,B3,PC,C119,A119,G,P,B,H,AD,CV,N,MS,P3,G2,SP=[D[k] for k in paths]
    assert TH['all_passed'] and GP['all_passed'] and SP['all_passed']
    assert sha(ROOT/paths['theorem'])=='95a6207dad460626caf6c80492650d8d7e62ff43166c4f7c0ed476efd4bb9462'
    assert sha(ROOT/paths['gamma_peer'])=='a73e2516ad0aabcb47b2d644a95e79bde76f849bc603d2f11387a2d02531fe71'
    assert sha(ROOT/paths['source_peer'])=='808637882f3d8d28f393a7b0f4abfeec59e076bd1269034bd67d8ac6dd8c706a'
    mons2=list(itertools.combinations_with_replacement(range(8),2));mi2={m:i for i,m in enumerate(mons2)}
    mons3=list(itertools.combinations_with_replacement(range(8),3))
    assert GA['monomials2_36']==PC['monomials2_36']==[list(m) for m in mons2]
    assert GA['monomials3_120']==PC['monomials3_120']==[list(m) for m in mons3]
    b3=F(B3['canonical_qB_pure_cubic_b3']);assert b3==-F(42215,3456)
    gam3=rec(GA['Gamma0_56_times120']);addscaled(gam3,rec(GA['Gammab_56_times120']),b3)
    gamC=rec(PC['two_Xsource3_adjoint_rows56_times120'][0]);assert PC['two_quotient_dual_rows600'][0]==S['Theta4_cubic600']['times_c4']
    xp1=list(map(rec,P['actual_pivot_X1_by27_sparse56_times8']))
    ops=B['operator_scalar_nodes'];bcols=list(map(dict,B['B48_by640_sparse_columns']));nodes=H['DAG_scalar_nodes']
    def key(j,src,m):
        if j>=640:return (2,j-640,src,m) if src<48 else None
        if src<48:
            root=bcols[j].get(src)
            if root is None or ops[root][0]=='c':return None
            op=ops[root];return (1,op[1],m) if op[0]=='actual_C0_scalar_function' else (0,root,m)
        root=H['honest_D0_relation_columns'][j][src-48+1]
        return None if nodes[root][0]=='c' else (1,root,m)
    def c1_functional(c1,a2,gam):
        f={};known=F(0)
        for index,c in gam.items():
            src,m=divmod(index,120);mon=mons3[m]
            for v in set(mon):
                pair=list(mon);pair.remove(v);m2=mi2[tuple(pair)]
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
                if nodes[root][0]!='c':put(z,root,a)
        return z
    def xp_correction(a2):
        z={}
        for p in range(27):
            c=a2.get(36*p+35,F(0))
            if c:
                for root,a in G['normalized_Xp27_source55_rr'][p]['rawP8_root_coefficients']:put(z,root,-c*F(a))
        return z
    def functional(c1,a2):
        p=rawP(c1[7]);xp=xp_correction(a2);f0,k0=c1_functional(c1,a2,gam3);fC,kC=c1_functional(c1,a2,gamC)
        return p,xp,f0,fC,k0,kC
    pivot=71;scale=-F(251099375,1051776)
    assert S['c4_source_normalization']['actual119_pivot']==pivot and F(S['c4_source_normalization']['divide'])==scale
    c1pivot=list(map(rec,C119['all119_linear_coefficient_original768_combinations'][pivot]));a2pivot=rec(A119['all119_actual_a2_sparse27_times36'][pivot])
    pivot_functional=functional(c1pivot,a2pivot)
    cases=[];labels=set();proots=set();xproots=set();support=S['QQ_support_bound']
    def absorb(f):
        nonlocal support
        p,xp,f0,fC,k0,kC=f;labels.update(f0);labels.update(fC);proots.update(p);xproots.update(xp)
        support+=len(p)+len(xp)+len(f0)+len(fC)+bool(k0)+bool(kC)
    absorb(pivot_functional)
    for row in S['all168']:
        assert row['all56_x2_equals_Xp0_a2'] and row['native_and_old_source_origin_parts_separately_zero']
        assert all(i%8!=7 for i,a in row['b1_U0_inverse_D3linear84_times8'])
        c1=list(map(rec,row['linear_original768_combinations_by8variables']));a2=rec(row['actual_a2_sparse27_times36'])
        f=functional(c1,a2);cases.append(f);absorb(f)
    assert len(cases)==168 and support<150000
    def serial(f):
        p,xp,f0,fC,k0,kC=f
        return {'rawP3':encode(p),'XP2':encode(xp),'X2_Theta3':enc(f0),'X2_times_c4':enc(fC),'known_Theta3':str(k0),'known_times_c4':str(kC)}
    dump(out,'phase1-c4-and168-source-functionals.json',{'monomials2_36':mons2,'monomials3_120':mons3,
      'Gamma3_56_times120':encode(gam3),'GammaC_56_times120':encode(gamC),'canonical_qB_pure_cubic':str(b3),
      'canonical_c4_source':{'pivot119':pivot,'divide':str(scale),'lower598_corrected_functional':serial(pivot_functional)},
      'all168':[{**serial(f),'literal_d3_source':i,'corrected_d3_is_NEGATIVE_source_functional':True} for i,f in enumerate(cases)],
      'ordinary_source_formula':'rawP3+XP2+X2_Theta3+known_Theta3+c4*(X2_times_c4+known_times_c4)',
      'no_new_source_coefficient_jet_or_inverse':True})
    allowed=[(0,7),(1,7),(4,7),(5,7),(7,7)];pairB=dict(zip(allowed,[7,8,9,10,2]));pairAD={p:j for j,p in enumerate(allowed)}
    pairs={mons2[k[-1]] for k in labels};extra=sorted(pairs-set(allowed))
    pre={'actual168_count':168,'all_source_signs_negative_ordinary_functional':True,'Theta4':'Theta3+c4*PiC; PiC annihilates qB*r',
      'canonical_c4_pivot119':pivot,'canonical_c4_divide':str(scale),'source_labels':len(labels),'required_pairs':sorted(pairs),
      'pairs_outside_B11':extra,'rawP3_roots':len(proots),'XP2_roots':len(xproots),'source_support':support,'cap':150000}
    dump(out,'preflight.json',pre)
    if extra:
        dump(out,'report.json',{'id':STEM+'_SOURCE_PREFLIGHT','all_passed':True,'all_required_values_cached':False,
          'missing_X2_pairs':extra,'actual168_values_evaluated':False,'no_new_source_coefficients':True,'input_sha256':hashes,
          'code_sha256':sha(Path(__file__)),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start})
        print(json.dumps(pre));return
    td={i:(rec(g),rec(q)) for i,v,g,f,q in AD['node_data']}
    cv={i:(rec(v),rec(g)) for i,v,g in CV['all108_composite_values_and_opaque_second_gradients']}
    native=[{row:rec(v) for row,v in col} for col in N['source_columns128_rows48_sparse_B11_effective_second']]
    typed={k[1] for k in labels if k[0]==0};direct={k[1] for k in labels if k[0]==1};atoms={a for i in typed for a in td[i][0]}
    assert typed<=td.keys() and direct<=cv.keys() and atoms<=cv.keys()
    opaque={i:cv[i][1] for i in atoms|direct if cv[i][1]};assert not opaque,('OPAQUE_SECOND_SENSITIVITY',opaque)
    jets={i:list(map(F,p)) for i,p in P3['raw_source_roots_with_missing_quotient2_3_set_zero_ONLY_for_known_part']}
    jac={i:(rec(a),rec(b)) for i,a,b in MS['source_root_Jacobians0_and_r']}
    assert proots<=jets.keys() and proots<=jac.keys()
    assert all(not jac[i][0] and not jac[i][1] for i in proots),('OPAQUE_P3_SENSITIVITY',[(i,jac[i]) for i in proots if any(jac[i])])
    g2={i:F(a) for i,a in G2['actual_eight_regular_g_second_values']}
    xpvalues={i:sum((a*g2[j] for j,a in jac[i][0].items()),F(0)) for i in xproots}
    assert all(not jac[i][1] and not jets[i][2] for i in xproots)
    whole=support+sum(len(v) for col in native for v in col.values())+sum(len(td[i][0]) for i in typed)+len(labels)+10000
    assert whole<150000
    pre.update({'whole_visible_QQ_support_bound':whole,'all_required_values_already_cached':True,'all_unavailable_second_third_sensitivities_zero':True})
    dump(out,'preflight.json',pre);dump(out,'phase.json',{'phase':'WHOLE_SOURCE_CACHE_PREFLIGHT_PASS_BEFORE_VALUES','support':whole})
    def value(label):
        kind=label[0];p=mons2[label[-1]];b=pairB[p]
        if kind==2:return native[label[1]].get(label[2],{}).get(b,F(0))
        if kind==1:return cv[label[1]][0].get(b,F(0))
        g,q=td[label[1]];return q.get(pairAD[p],F(0))+sum((a*cv[i][0].get(b,F(0)) for i,a in g.items()),F(0))
    values={label:value(label) for label in labels}
    def components(f):
        p,xp,f0,fC,k0,kC=f
        a=k0+sum((c*jets[i][3] for i,c in p.items()),F(0))+sum((c*xpvalues[i] for i,c in xp.items()),F(0))+sum((c*values[k] for k,c in f0.items()),F(0))
        b=kC+sum((c*values[k] for k,c in fC.items()),F(0));return a,b
    pivotraw,pivot_extra=components(pivot_functional);c4=pivotraw/scale
    actual=[-(a+c4*b) for a,b in map(components,cases)];nz=[i for i,a in enumerate(actual) if a]
    dump(out,'actual-c4-and168-corrected-d3-quartic-values.json',{'canonical119_pivot':pivot,'pivot_lower598_corrected_pure_quartic':str(pivotraw),
      'canonical_qC_pure_quartic_c4':str(c4),'all_X2_scalar_root_pair_values':[[list(k),str(v)] for k,v in sorted(values.items())],
      'all168_actual_corrected_d3_pure_quartic_values':list(map(str,actual)),'nonzero_indices':nz,
      'all168_source_signs_checked_as_negative_ordinary_functional':True,
      'actual_HF4_if_nonzero':0 if nz else None,'actual_T_total_length_if_nonzero':11 if nz else None,
      'allzero_scope':'Only the corrected168 d3 family; actualHF4 and total length remain open',
      'new_source_coefficients_computed':0})
    dump(out,'phase.json',{'phase':'FINAL','nonzero168_values':len(nz),'canonical_c4':str(c4)})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_nonzero_corrected168_quartic_values':len(nz),'canonical_qC_pure_quartic':str(c4),
      'actual_HF4':0 if nz else None,'actual_T_total_length':11 if nz else None,'allzero_means_only_this_family_negative':not nz,
      'new_source_coefficient_jet_inverse_or_solve':False,'independent_source_value_review':'PENDING','input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual_nonzero_corrected168_quartic_values','canonical_qC_pure_quartic','actual_T_total_length','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
