#!/usr/bin/env python3
"""Source-only0/r Jacobian inventory for old619 high quartic functionals."""
import argparse,collections,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63'
Z=(F(0),F(0));ONE=(F(1),F(0))
def da(a,b):return (a[0]+b[0],a[1]+b[1])
def dm(a,b):return (a[0]*b[0],a[0]*b[1]+a[1]*b[0])
def ds(a,c):return (a[0]*c,a[1]*c)
def di(a):assert a[0];return (1/a[0],-a[1]/a[0]**2)
def gadd(a,b,c=ONE):
    a=dict(a)
    for k,v in b.items():
        z=da(a.get(k,Z),dm(c,v))
        if z!=Z:a[k]=z
        else:a.pop(k,None)
    return a
def put(a,k,v):
    z=a.get(k,F(0))+v
    if z:a[k]=z
    else:a.pop(k,None)
def enc(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def gdump(v):return [[list(k),list(map(str,a))] for k,a in sorted(v.items())]
def main(out):
    st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    def dump(n,d):
        (out/n).write_text(json.dumps(d,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths={
      'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'source2':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
      'own':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
      'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
      'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'common3':'reports/AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63/phase1-existing-source3-and-copy-composition.json',
      'nativeN2':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/report.json',
      'lower':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/report.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};dump('phase.json',{'phase':'joint_high_old_source_inventory_inputs_persisted','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};HI,C,B,S,J,CP,AF,OLD,C3,NN,LO=[D[k] for k in paths]
    assert NN['all_passed'] and LO['all_passed'];nodes=C['DAG_scalar_nodes'];ops=B['operator_scalar_nodes']
    mon20=[tuple(m) for m in S['base20_monomials_original8']];rindex=mon20.index((7,))
    common={i:(rec(v).get(0,F(0)),rec(v).get(rindex,F(0))) for i,v in S['all105_actual_source2_values']}
    copyfirst={i:(rec(v).get(0,F(0)),rec(v).get(rindex,F(0))) for i,v in S['all57_actual_copy_source2_values']}
    cb=J['both57_source2']['CB'];ownmons=[tuple(m) for m in J['source2_monomials']];mi={m:i for i,m in enumerate(ownmons)}
    byown={i:list(map(F,p)) for i,p in cb['all57_source2jets']};callbacks={new:old for old,new in CP['explicit57_source_copy_bindings']};yb=cb['F9_nodes'][7:9]
    yfirst=[tuple(map(F,p)) for i,p in C3['actual_YB_last2_through3']];assert [p[1] for p in yfirst]==[F(17,12),F(7,12)]
    opaque={}
    for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
        for i,p in J[key].items():opaque[int(i)]=(F(p[0]),F(p[8])+F(p[9]))
    for q in AF['all35_ordered_CA_prefix_firsts']:
        p=rec(q['first9_bar']);opaque[q['node']]=(p.get(0,F(0)),p.get(8,F(0)))
    for i,q in AF['all8_g_quotient_overrides']:
        p=rec(q);opaque[i]=(p.get(0,F(0)),p.get(8,F(0)))
    coords={i:Z for i in C['actual_C0_basis8_scalar_nodes']}
    troots={};croots=set();p4roots=set();xp3roots=set()
    for row in HI['all619']:
        p4roots.update(i for i,a in row['rawP55_pure_r4_roots']);xp3roots.update(i for i,a in row['Xp55_pure_r3_times_a1r_roots'])
        for k,a in row['actual_X3_high_order_source_functional']:
            if k[0]==0:troots[k[1]]=True
            elif k[0]==1:croots.add(k[1])
    croots|=p4roots|xp3roots;tc=set();todo=list(troots)
    while todo:
        i=todo.pop()
        if i in tc:continue
        tc.add(i);n=ops[i]
        if n[0]=='sum':todo.extend(j for j,a,b in n[1])
        elif n[0]=='product':todo.extend(n[3])
        elif n[0]=='actual_C0_scalar_function':croots.add(n[1])
        else:assert n[0]=='c'
    needed=set();todo=list(croots);missing={}
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i)
        if i in common or i in opaque or i in coords:continue
        if i in callbacks:todo.extend(yb);continue
        n=nodes[i]
        if n[0] in ['a','m','n','i']:todo.extend(n[1:])
        elif n[0]!='c':missing[i]=n
    dump('phase1-source-closure-and-bindings.json',{'C0_roots':sorted(croots),'typed_roots':sorted(troots),'C0_closure_nodes':len(needed),'typed_closure_nodes':len(tc),
      'missing_source_leaves':[[i,n] for i,n in sorted(missing.items())],
      'copied57_source_rule':'ownCB zero7 inverse psi(0^7,YB7,YB8); differentiate savedpsi1/psi2 coherently',
      'gA_rule':'actual zero7 conductor extension, known ownr/s source3/4; never fullf63',
      'source105_higher_coeffs_and_honest_higher_coeffs_remain_formal':True})
    assert not missing,('EXACT_MISSING_SOURCE_BINDING',missing)
    def children(i):
        if i in common or i in opaque or i in coords:return []
        if i in callbacks:return yb
        n=nodes[i];return n[1:] if n[0] in ['a','m','n','i'] else []
    uses=collections.Counter(j for i in needed for j in children(i));uses.update(croots);v={};g={};live=0;peak=0
    for i in sorted(needed):
        n=nodes[i];o=n[0]
        if i in common:x=common[i];grad={('physical',i):ONE}
        elif i in coords:x=coords[i];grad={}
        elif i in opaque:x=opaque[i];grad={('honest',i):ONE}
        elif i in callbacks:
            p=byown[callbacks[i]];x=copyfirst[i];grad={}
            for k,root in enumerate(yb):
                d0=p[mi[(7+k,)]];d1=F(0)
                for m,a in zip(ownmons,p):
                    if len(m)==2 and 7+k in m and all(t>=7 for t in m):
                        rest=list(m);rest.remove(7+k);d1+=a*m.count(7+k)*v[yb[rest[0]-7]][1]
                grad=gadd(grad,g[root],(d0,d1))
            assert x[1]==sum((p[mi[(7+k,)]]*v[root][1] for k,root in enumerate(yb)),F(0))
        elif o=='c':x=(F(n[1],n[2]),F(0));grad={}
        elif o=='a':x=da(v[n[1]],v[n[2]]);grad=gadd(g[n[1]],g[n[2]])
        elif o=='n':x=ds(v[n[1]],F(-1));grad={k:ds(a,F(-1)) for k,a in g[n[1]].items()}
        elif o=='m':x=dm(v[n[1]],v[n[2]]);grad=gadd({},g[n[1]],v[n[2]]);grad=gadd(grad,g[n[2]],v[n[1]])
        elif o=='i':x=di(v[n[1]]);grad={k:ds(dm(a,dm(x,x)),F(-1)) for k,a in g[n[1]].items()}
        else:raise AssertionError((i,n))
        v[i]=x;g[i]=grad;live+=2+2*len(grad);peak=max(peak,live);assert peak+20000<150000
        for j in children(i):
            uses[j]-=1
            if not uses[j]:live-=2+2*len(g[j]);del g[j];del v[j]
    dump('phase2-C0-source-Jacobians0-r.json',{'all_C0_roots':[[i,list(map(str,v[i])),gdump(g[i])] for i in sorted(croots)],
      'leaf_kind_physical_is_original_common105_unknown':True,'leaf_kind_honest_is_true_opaque_regular_function':True,
      'copy57_sources_not_independent_variables':True,'peak_live_QQ':peak})
    usesT=collections.Counter()
    for i in tc:
        n=ops[i]
        if n[0]=='sum':usesT.update(j for j,a,b in n[1])
        elif n[0]=='product':usesT.update(n[3])
    usesT.update(troots);tv={};tg={};tlive=0;tpeak=0;old={i:(F(a),rec(f)) for i,a,gg,f,q in OLD['node_data']};checks=0
    for i in sorted(tc):
        n=ops[i];o=n[0]
        if o=='c':x=(F(n[1],n[2]),F(0));grad={}
        elif o=='actual_C0_scalar_function':x=v[n[1]];grad=g[n[1]]
        elif o=='sum':
            x=Z;grad={}
            for j,a,b in n[1]:x=da(x,ds(tv[j],F(a,b)));grad=gadd(grad,tg[j],(F(a,b),F(0)))
        elif o=='product':
            x=(F(n[1],n[2]),F(0));grad={}
            for j in n[3]:x=dm(x,tv[j])
            for k,j in enumerate(n[3]):
                c=(F(n[1],n[2]),F(0))
                for l,h in enumerate(n[3]):
                    if l!=k:c=dm(c,tv[h])
                grad=gadd(grad,tg[j],c)
        else:raise AssertionError(n)
        if i in old:assert x==(old[i][0],old[i][1].get(7,F(0))),('OLD_TYPED_FIRST_MISMATCH',i);checks+=1
        tv[i]=x;tg[i]=grad;tlive+=2+2*len(grad);tpeak=max(tpeak,tlive);assert peak+tpeak+20000<150000
        ch=[j for j,a,b in n[1]] if o=='sum' else n[3] if o=='product' else []
        for j in ch:
            usesT[j]-=1
            if not usesT[j]:tlive-=2+2*len(tg[j]);del tv[j];del tg[j]
    dump('phase3-typed-high-old-source-Jacobians0-r.json',{'typed_roots':[[i,list(map(str,tv[i])),gdump(tg[i])] for i in sorted(troots)],
      'old_typed_origin_first_checks':checks,'peak_live_QQ':tpeak})
    records=[];active3=set();active4=set();support=peak+tpeak
    for row in HI['all619']:
        a3={};a4={}
        for root,c in row['rawP55_pure_r4_roots']:
            for leaf,d in g[root].items():put(a4,leaf,F(c)*d[0]);put(a3,leaf+(119,),F(c)*d[1])
        for root,c in row['Xp55_pure_r3_times_a1r_roots']:
            for leaf,d in g[root].items():put(a3,leaf+(119,),F(c)*d[0])
        for key,c in row['actual_X3_high_order_source_functional']:
            if key[0]==2:continue
            gg=tg[key[1]] if key[0]==0 else g[key[1]]
            for leaf,d in gg.items():put(a3,leaf+(key[-1],),F(c)*d[0])
        active3.update(a3);active4.update(a4);support+=len(a3)+len(a4);assert support+20000<150000
        records.append({'kernel619_index':row['kernel619_index'],'old_high_source3_sensitivity':enc(a3),'old_high_source4_sensitivity':enc(a4)})
    inventory=[]
    for kind,i in sorted({k[:2] for k in active3}|active4):
        inventory.append({'kind':kind,'source_node':i,'literal_node':nodes[i],
          'source3_monomials':[mons3 for mons3 in sorted({key[2] for key in active3 if key[:2]==(kind,i)})],
          'source4_r4_active':(kind,i) in active4,
          'common105_r3_already_cached':kind=='physical','ordinary_honest_higher_still_unknown':kind=='honest'})
    dump('phase4-all619-high-old-source-sensitivities.json',{'all619':records,'active_source_inventory':inventory,
      'source3_monomials_order':list(itertools.combinations_with_replacement(range(8),3)),
      'native_G3_and_monic3_source_sensitivities_NOT_combined_yet':True,'lower_honest_second_products_NOT_evaluated_here':True,
      'no_new_source3_4_6_or_gauge3_values_computed':True})
    dump('phase.json',{'phase':'FINAL_SOURCE_ONLY','source3_entries':len(active3),'source4_entries':len(active4)})
    note=ROOT/'research/notes'/f'{STEM}.md'
    r={'id':STEM,'all_passed':True,'old_high_source3_active_pairs':len(active3),'old_high_source4_active_leaves':len(active4),
      'old_source3_nonzero_output_count':sum(bool(q['old_high_source3_sensitivity']) for q in records),
      'old_source4_nonzero_output_count':sum(bool(q['old_high_source4_sensitivity']) for q in records),
      'new_source_or_gauge_coefficients':0,'native_source3_not_combined':True,'actual619_values_or_HF4_computed':False,
      'total_live_saved_QQ_bound':support+20000,'input_sha256':hs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    dump('report.json',r);print(json.dumps({k:r[k] for k in ['id','all_passed','old_high_source3_active_pairs','old_high_source4_active_leaves','old_source3_nonzero_output_count','old_source4_nonzero_output_count','total_live_saved_QQ_bound','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
