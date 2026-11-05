#!/usr/bin/env python3
"""Common105 source2 by saved inverse composition; oldB AD stays symbolic."""
import argparse,collections,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63'
ACTIVE=[0,1,4,5,7]
MON=[()]+[(v,) for v in ACTIVE]+[p for p in itertools.combinations_with_replacement(ACTIVE,2) if p!=(1,1)]
IDX={m:i for i,m in enumerate(MON)}
def put(v,k,a):
    c=v.get(k,F(0))+a
    if c:v[k]=c
    else:v.pop(k,None)
def add(a,b,c=F(1)):
    a=dict(a);addscaled(a,b,c);return a
def mul(a,b):
    z={}
    for i,c in a.items():
        for j,d in b.items():
            k=IDX.get(tuple(sorted(MON[i]+MON[j])))
            if k is not None:put(z,k,c*d)
    return z
def inv(a):
    c=a.get(0,F(0));assert c;h={i:v for i,v in a.items() if i};z={0:1/c};addscaled(z,h,-1/c**2);addscaled(z,mul(h,h),1/c**3)
    assert mul(a,z)=={0:F(1)};return z
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def main(out,resume=None):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    def dump(n,x):
        (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths={
      'source2':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'inverse2':'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json',
      'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
      'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
      'copy2':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
      'composite':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
      'lower':'reports/AMBIENT_B18_T_COMMON_LOWER_ADJOINT_20261104_v63/phase3-all619-common-D-contractions.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json'}
    if resume:paths['completed_source_cache']=str(resume)
    hs={p:sha(ROOT/p) for p in paths.values()};dump('phase.json',{'phase':'source_inverse_composition_inputs_bound','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};J,Q,S,JC,B11,C,CP,CP2,CV,AF,LO,CH=[D[k] for k in list(paths)[:12]]
    ownmons=[tuple(m) for m in J['source2_monomials']];assert ownmons==[tuple(m) for m in Q['source2_monomials']]
    vars8=[{IDX[(i,)]:F(1)} if i in ACTIVE else {} for i in range(8)];YA=vars8[:7]+[vars8[7],vars8[7]]
    def compose(p,args):
        z={}
        for m,c in zip(ownmons,p):
            c=F(c)
            if not c:continue
            v={0:c}
            for i in m:v=mul(v,args[i])
            addscaled(z,v,F(1))
        return z
    if not resume:
        qa=Q['families']['CA'];qb=Q['families']['CB'];assert qa['physical_R0_rows9']==qb['physical_R0_rows9']
        xa=[compose(p,YA) for p in qa['physical_R0_readback_jets_in_own_F9']]
        YB=[compose(p,xa) for p in qb['inverse_F9_in_common_R0_through2']]
        assert [compose(p,YB) for p in qb['physical_R0_readback_jets_in_own_F9']]==xa
        assert [compose(p,xa) for p in qa['inverse_F9_in_common_R0_through2']]==YA
        ca={i:compose(p,YA) for i,p in J['both57_source2']['CA']['all57_source2jets']}
        cb={i:compose(p,YB) for i,p in J['both57_source2']['CB']['all57_source2jets']}
        assert len(set(ca)&set(cb))==9 and all(ca[i]==cb[i] for i in set(ca)&set(cb))
        common=ca|cb;assert len(common)==105
    else:
        cached=D['completed_source_cache'];assert cached['base20_monomials_original8']==[list(m) for m in MON]
        assert cached['all105_equations_verified'] and cached['all57_copied_equations_verified'] and cached['all105_and57_B11_projections_exact']
        common={i:rec(v) for i,v in cached['all105_actual_source2_values']};copyvals={i:rec(v) for i,v in cached['all57_actual_copy_source2_values']}
        YB=list(map(rec,cached['actual_YB9_nonlinear_through2']))
    Bmons=[tuple(i for i,n in enumerate(m) for _ in range(n)) for m in B11['B11_basis_monomials8']]
    bi={m:i for i,m in enumerate(Bmons)}
    def project(p):return {bi[MON[i]]:a for i,a in p.items() if MON[i] in bi}
    if not resume:
        oldB={i:rec(p) for i,p in B11['all105_actual_B11_root_values']}
        assert set(common)==set(oldB) and all(project(p)==oldB[i] for i,p in common.items())
        bindcopy=dict(CP['explicit57_source_copy_bindings']);copyargs=[{}]*7+YB[7:9]
        copyvals={bindcopy[i]:compose(p,copyargs) for i,p in J['both57_source2']['CB']['all57_source2jets']}
        oldcopy={i:rec(p) for i,p in CP2['all57_copied_unknown_values_through2']}
        assert all(project(p)==oldcopy[i] for i,p in copyvals.items())
    nodes=C['DAG_scalar_nodes'];assert nodes[:len(S['DAG_scalar_nodes'])]==S['DAG_scalar_nodes']
    eq105=JC['all105_equation_original_node_order'];eqcopy=CP['all57_copied_equation_indices']
    def closure(roots,stops):
        ns=set();todo=list(roots);missing={}
        while todo:
            i=todo.pop()
            if i in ns:continue
            ns.add(i)
            if i in stops:continue
            n=nodes[i]
            if n[0] in ['a','m']:todo.extend(n[1:3])
            elif n[0] in ['n','i']:todo.append(n[1])
            elif n[0]!='c':missing[i]=n
        return sorted(ns),missing
    def evaluate(roots,bind,opaque=None):
        opaque=opaque or {};order,missing=closure(roots,set(bind)|set(opaque));assert not missing,('MISSING_SOURCE_BINDING',missing)
        uses=collections.Counter()
        for i in order:
            if i in bind or i in opaque:continue
            n=nodes[i]
            if n[0] in ['a','m','n','i']:uses.update(n[1:])
        uses.update(roots);val={};grad={};live=0;peak=0
        for i in order:
            n=nodes[i];o=n[0]
            if i in bind:v=bind[i];g={}
            elif i in opaque:v=opaque[i];g={i:F(1)}
            elif o=='c':v={0:F(n[1],n[2])} if n[1] else {};g={}
            elif o=='a':v=add(val[n[1]],val[n[2]]);g=add(grad[n[1]],grad[n[2]])
            elif o=='n':v={k:-a for k,a in val[n[1]].items()};g={k:-a for k,a in grad[n[1]].items()}
            elif o=='m':
                v=mul(val[n[1]],val[n[2]]);g={};addscaled(g,grad[n[1]],val[n[2]].get(0,F(0)));addscaled(g,grad[n[2]],val[n[1]].get(0,F(0)))
            elif o=='i':v=inv(val[n[1]]);g={k:-a/val[n[1]][0]**2 for k,a in grad[n[1]].items()}
            else:raise AssertionError((i,n))
            val[i]=v;grad[i]=g;live+=len(v)+len(g);peak=max(peak,live);assert peak+20000<150000
            if i not in bind and i not in opaque:
                for j in n[1:] if o in ['a','m','n','i'] else []:
                    uses[j]-=1
                    if not uses[j]:live-=len(val[j])+len(grad[j]);del val[j];del grad[j]
        return {i:val[i] for i in roots},{i:grad[i] for i in roots},len(order),peak
    if not resume:
        eqvals,_,ec,ep=evaluate(eq105+eqcopy,common|copyvals)
        for j,i in enumerate(eq105):assert eqvals[i]==(YA[j-48] if 48<=j<57 else {})
        assert all(not eqvals[i] for i in eqcopy)
        dump('phase1-actual105-source2-and-copy.json',{'base20_monomials_original8':MON,'actual_YA9':list(map(encode,YA)),
          'actual_YB9_nonlinear_through2':list(map(encode,YB)),'all105_actual_source2_values':[[i,encode(v)] for i,v in sorted(common.items())],
          'all57_actual_copy_source2_values':[[i,encode(v)] for i,v in sorted(copyvals.items())],
          'all105_equations_verified':True,'all57_copied_equations_verified':True,'all105_and57_B11_projections_exact':True,
          'no_implicit_or_gauge_solve':True,'equation_closure_nodes':ec,'peak_live_QQ':ep})
    else:dump('phase1-resumed-source-cache.json',{'source_path':str(resume),'sha256':sha(ROOT/resume),'recomposed_source105':False})
    dump('phase.json',{'phase':'ALL105_SOURCE2_COMPOSITION_VERIFIED_BEFORE_OLDB_AD'})
    # Original honest quotient second coefficients stay independent symbols.
    opaque={}
    def fromfirst(p):
        p=list(map(F,p));z={0:p[0]} if p[0] else {}
        for v in ACTIVE:
            a=(p[8]+p[9]) if v==7 else p[1+v]
            if a:z[IDX[(v,)]]=a
        return z
    for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
        for i,p in J[key].items():opaque[int(i)]=fromfirst(p)
    for row in AF['all35_ordered_CA_prefix_firsts']:
        p=rec(row['first9_bar']);opaque[row['node']]={0:p[0]} if p.get(0) else {}
        for v in ACTIVE:
            if p.get(1+v):opaque[row['node']][IDX[(v,)]]=p[1+v]
    for i,p in AF['all8_g_quotient_overrides']:
        p=rec(p);opaque[i]={0:p[0]} if p.get(0) else {}
        for v in ACTIVE:
            if p.get(1+v):opaque[i][IDX[(v,)]]=p[1+v]
    bind=common|copyvals
    for v,i in enumerate(C['actual_C0_basis8_scalar_nodes']):bind[i]=vars8[v] if v<7 else {}
    tnodes=CH['operator_scalar_nodes'];oldlabels={tuple(k) for row in LO['all619_source_functionals'] for k,a in row['X2_source_functional'] if k[0]<2}
    troots={k[1] for k in oldlabels if k[0]==0};croots={k[1] for k in oldlabels if k[0]==1};tc=set();todo=list(troots)
    while todo:
        i=todo.pop()
        if i in tc:continue
        tc.add(i);n=tnodes[i]
        if n[0]=='actual_C0_scalar_function':croots.add(n[1])
        elif n[0]=='sum':todo.extend(j for j,a,b in n[1])
        elif n[0]=='product':todo.extend(n[3])
        elif n[0] in ['a','m','n','i']:todo.extend(n[1:])
        else:assert n[0]=='c',('UNKNOWN_TYPED_SOURCE',i,n)
    oldcv={i:(rec(v),rec(g)) for i,v,g in CV['all108_composite_values_and_opaque_second_gradients']};croots.update(oldcv)
    corder,missing=closure(croots,set(bind)|set(opaque))
    profile={'oldB_typed_root_count':len(troots),'typed_closure_nodes':len(tc),'C0_root_count':len(croots),'C0_closure_nodes':len(corder),
      'missing_source_leaves':[[i,n] for i,n in sorted(missing.items())],'no_name_heuristic_or_unknown_gradient_invention':True,
      'all_native_extra_second_entries_remain_symbolic':True,'no_new_gauge_coefficient_or_source3_or6':True}
    dump('phase2-oldB-source-profile.json',profile)
    if missing:
        result={'oldB_AD_complete':False,'missing_source_leaves':profile['missing_source_leaves'],'source105_composition_complete':True}
    else:
        cv,cg,cc,peak=evaluate(croots,bind,opaque)
        for i,(v,g) in oldcv.items():assert project(cv[i])==v and cg[i]==g,('OLD_B11_COMPOSITE_MISMATCH',i)
        tv={};tg={};typed_support=0
        for i in sorted(tc):
            n=tnodes[i];o=n[0]
            if o=='c':v={0:F(n[1],n[2])} if n[1] else {};g={}
            elif o=='actual_C0_scalar_function':v=cv[n[1]];g=cg[n[1]]
            elif o=='sum':
                v={};g={}
                for j,a,b in n[1]:addscaled(v,tv[j],F(a,b));addscaled(g,tg[j],F(a,b))
            elif o=='product':
                v={0:F(n[1],n[2])};g={}
                for j in n[3]:v=mul(v,tv[j])
                for k,j in enumerate(n[3]):
                    c=F(n[1],n[2])
                    for l,h in enumerate(n[3]):
                        if l!=k:c*=tv[h].get(0,F(0))
                    addscaled(g,tg[j],c)
            elif o=='a':v=add(tv[n[1]],tv[n[2]]);g=add(tg[n[1]],tg[n[2]])
            elif o=='n':v={k:-a for k,a in tv[n[1]].items()};g={k:-a for k,a in tg[n[1]].items()}
            elif o=='m':v=mul(tv[n[1]],tv[n[2]]);g={};addscaled(g,tg[n[1]],tv[n[2]].get(0,F(0)));addscaled(g,tg[n[2]],tv[n[1]].get(0,F(0)))
            elif o=='i':v=inv(tv[n[1]]);g={k:-a/tv[n[1]][0]**2 for k,a in tg[n[1]].items()}
            else:raise AssertionError((i,n))
            tv[i]=v;tg[i]=g;typed_support+=len(v)+len(g)
            assert typed_support+sum(len(v)+len(cg[i]) for i,v in cv.items())+20000<150000
        dump('phase3-oldB-common-source2-values-and-symbols.json',{'base20_monomials_original8':MON,
          'C0_composite_values_and_opaque_second_gradients':[[i,encode(cv[i]),encode(cg[i])] for i in sorted(croots)],
          'typedB_root_values_and_opaque_second_gradients':[[i,encode(tv[i]),encode(tg[i])] for i in sorted(troots)],
          'all108_old_B11_values_and_gradients_reproduced':True,'unknown_second_symbols_are_true_honest_source_leaves':True})
        pairs=list(itertools.combinations_with_replacement(range(8),2));actual=[];native=[];sens=[]
        for row in LO['all619_source_functionals']:
            z=F(row['known_Xp1_a1']);g={};nn={}
            for k,a in row['X2_source_functional']:
                k=tuple(k);a=F(a);m=pairs[k[-1]];bi20=IDX[m]
                if k[0]==2:put(nn,k,a);continue
                vv,gg=(tv[k[1]],tg[k[1]]) if k[0]==0 else (cv[k[1]],cg[k[1]])
                z+=a*vv.get(bi20,F(0))
                for leaf,b in gg.items():put(g,(leaf,k[-1]),a*b)
            actual.append(z);native.append(nn);sens.append(g)
        dump('phase4-all619-lower-knownparts-and-symbols.json',{'all619_oldB_known_parts':list(map(str,actual)),
          'all619_opaque_honest_second_sensitivities':list(map(enc,sens)),
          'all619_uncomputed_native_second_source_functionals':list(map(enc,native)),
          'all_ordinary_honest_second_sensitivities_zero':not any(sens),
          'all_native_second_values_still_symbolic':True,'not_complete619_lower_or_quartic_values':True})
        result={'oldB_AD_complete':True,'oldB_opaque_nonzero_source_count':sum(bool(g) for g in sens),
          'opaque_QQ_entries':sum(map(len,sens)),'native_uncomputed_source_count':sum(bool(n) for n in native),
          'native_uncomputed_QQ_entries':sum(map(len,native)),'typed_closure_nodes':len(tc),'C0_closure_nodes':cc,'C0_peak_live_QQ':peak}
    dump('phase.json',{'phase':'FINAL','source105_complete':True,**result})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM+('_OLDB_AD' if resume else ''),'all_passed':True,'all105_actual_seconds_composed_without_solve':True,'all105_equations_and57_copy_and_B11_verified':True,
      **result,'new_physical_jet_or_gauge_solve':False,'actual619_quartic_HF4_or_length_computed':False,'input_sha256':hs,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump('report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','all105_actual_seconds_composed_without_solve','oldB_AD_complete','seconds']}));print(json.dumps(result))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--resume-source-cache',type=Path);a=p.parse_args();main(a.output,a.resume_source_cache)
