#!/usr/bin/env python3
"""Fold fixed three-class source-root labels; old B atom AD at saved values.

P8 functions remain opaque unless literally the same saved C0 atom.  No
uncached C0 source derivative is invented. Native X2 entries are kept typed.
"""
import argparse,itertools,json,resource,signal,sys,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63'
CAP=150000
def put(d,k,a):
    b=d.get(k,F(0))+a
    if b:d[k]=b
    else:d.pop(k,None)
def sparse(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'sens':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/projected-source-X2-sensitivity.json',
      'comb':'reports/AMBIENT_B18_COMPARISON_T_NEXT3_SOURCE_PREFLIGHT_20261104_v63/zero-linear620-source-combinations.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'atoms':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'proof':'reports/AMBIENT_B18_T_QUADRATIC_M1_CANCELLATION_PREFLIGHT_20261104_v63/actual-projected-three-class-X2-operator-preflight.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()}
    dump(out,'phase.json',{'phase':'source_hashes_persisted_before_root_fold','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};S,C,B,A,H,V=[D[k] for k in paths]
    mons=list(itertools.combinations_with_replacement(range(8),2));sens=list(map(rec,S['three_rows_sparse']))
    pairs=sorted({mons[i%36] for row in sens for i in row});assert pairs==[(0,7),(1,7),(4,7),(5,7),(7,7)]
    pairid={p:i for i,p in enumerate(pairs)};ops=B['operator_scalar_nodes'];cols=[dict(col) for col in B['B48_by640_sparse_columns']]
    atom={i:rec(v) for i,v in A['all79_atom_values']};atom0={i:v.get(0,F(0)) for i,v in atom.items()}
    xp1=list(map(rec,V['actual_pivot_X1_by27_sparse56_times8']));assert len(xp1)==27
    folded=[];native=[];known=[];roots=set();saved_support=0;oldonly=[]
    for number,row in enumerate(C['all620']):
        c=rec(row['combination_in_actual768']);alpha=rec(row['linear_Xp27_correction_row8j_plus_variable'])
        fs=[];ns=[];ks=[]
        for functional in sens:
            f={};n={};k=F(0)
            for index,a in functional.items():
                src,m=divmod(index,36);u,v=mons[m];p=pairid[(u,v)]
                # Known Xp1*alpha1 coefficient, with ordinary polynomial
                # coefficients rather than doubled symmetric Hessians.
                for pivot,x in enumerate(xp1):
                    k-=a*alpha.get(8*pivot+u,F(0))*x.get(8*src+v,F(0))
                    if u!=v:k-=a*alpha.get(8*pivot+v,F(0))*x.get(8*src+u,F(0))
                for j,b in c.items():
                    if j>=640:
                        # Native128 bottom eight rows are literally zero in
                        # the complete source operator, at all orders.
                        if src<48:put(n,(j-640,src,p),a*b)
                        continue
                    if src<48:
                        root=cols[j].get(src)
                        if root is None:continue
                        node=ops[root]
                        if node[0]=='c':continue
                        if node[0]=='actual_C0_scalar_function':key=(1,node[1],p)
                        else:key=(0,root,p);roots.add(root)
                    else:
                        root=H['honest_D0_relation_columns'][j][src-48+1]
                        if H['DAG_scalar_nodes'][root][0]=='c':continue
                        key=(1,root,p)
                    put(f,key,a*b)
            fs.append(f);ns.append(n);ks.append(k)
        folded.append(fs);native.append(ns);known.append(ks)
        if not any(ns):oldonly.append(number)
        saved_support+=sum(map(len,fs))+sum(map(len,ns))+sum(bool(k) for k in ks)
    assert saved_support<CAP
    dump(out,'phase1-folded-root-functionals.json',{'needed_pairs':pairs,'620_rows':[
      {'index':i,'three_oldB_P_root_functionals':list(map(sparse,fs)),
       'three_native_X2_functionals':list(map(sparse,native[i])),
       'three_known_Xp1alpha1_corrections':list(map(str,known[i]))}
      for i,fs in enumerate(folded)],'root_kind0_typed_B_kind1_original_C0':True,
      'all_P8_entries_retained_unless_literal_constant':True,'folded_QQ_support':saved_support})
    # Closure and a tight structural derivative support upper bound. Scalars
    # are evaluated ONLY as polynomials in the79 already saved atom origins.
    needed=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i);node=ops[i]
        if node[0]=='sum':todo.extend(j for j,a,b in node[1])
        elif node[0]=='product':todo.extend(node[3])
    val={};deps={};first={};qknown={}
    def product_partial_coeff(indices,omit):
        z=F(1)
        for k,j in enumerate(indices):
            if k not in omit:z*=val[j]
        return z
    for i in sorted(needed):
        n=ops[i]
        if n[0]=='c':val[i]=F(n[1],n[2]);deps[i]=set()
        elif n[0]=='actual_C0_scalar_function':val[i]=atom0[n[1]];deps[i]={n[1]}
        elif n[0]=='sum':
            val[i]=sum((F(a,b)*val[j] for j,a,b in n[1]),F(0));deps[i]=set().union(*(deps[j] for j,a,b in n[1] if a))
        elif n[0]=='product':
            js=n[3];val[i]=F(n[1],n[2])*product_partial_coeff(js,set());deps[i]=set()
            for k,j in enumerate(js):
                if product_partial_coeff(js,{k}):deps[i].update(deps[j])
        else:raise AssertionError(n)
    root_ad_bound=0
    for fs in folded:
        for f in fs:
            symbols=set()
            for (kind,i,p),a in f.items():
                if kind:symbols.add((i,p))
                else:symbols.update((j,p) for j in deps[i])
            root_ad_bound+=len(symbols)
    forecast=saved_support+sum(map(len,deps.values()))+root_ad_bound+14*len(needed)+2000
    pre={'typed_closure_nodes':len(needed),'structural_gradient_slots':sum(map(len,deps.values())),
      'folded_QQ_support':saved_support,'atom_functional_support_bound':root_ad_bound,
      'typed_0_first8_known_pairs5_upper':14*len(needed),'TOTAL_forecast':forecast,'cap':CAP,
      'P8_nonatom_roots_remain_opaque':True,'no_uncached_C0_variable_derivative':True}
    dump(out,'preflight.json',pre);assert forecast<CAP,('AD_FORECAST_STOP',forecast)
    dump(out,'phase.json',{'phase':'root_fold_and_strict_AD_forecast_pass','TOTAL_forecast':forecast})
    grad={}
    for i in sorted(needed):
        n=ops[i]
        if n[0]=='c':g={};f={};q={}
        elif n[0]=='actual_C0_scalar_function':g={n[1]:F(1)};f={k-1:a for k,a in atom[n[1]].items() if k};q={}
        elif n[0]=='sum':
            g={};f={};q={}
            for j,a,b in n[1]:
                addscaled(g,grad[j],F(a,b));addscaled(f,first[j],F(a,b));addscaled(q,qknown[j],F(a,b))
        elif n[0]=='product':
            js=n[3];g={};f={};q={};c=F(n[1],n[2])
            for k,j in enumerate(js):
                t=c*product_partial_coeff(js,{k})
                addscaled(g,grad[j],t);addscaled(f,first[j],t);addscaled(q,qknown[j],t)
            for k in range(len(js)):
                for l in range(k+1,len(js)):
                    t=c*product_partial_coeff(js,{k,l})
                    if not t:continue
                    for p,(u,v) in enumerate(pairs):
                        a=first[js[k]].get(u,F(0))*first[js[l]].get(v,F(0))
                        if u!=v:a+=first[js[k]].get(v,F(0))*first[js[l]].get(u,F(0))
                        if a:put(q,p,t*a)
        assert set(g)<=deps[i]
        grad[i]=g;first[i]=f;qknown[i]=q
    adrows=[];knownfinal=[];active=set();constant_only=[];final_support=0
    for number,fs in enumerate(folded):
        arow=[];krow=[]
        for t,f in enumerate(fs):
            a={};k=known[number][t]
            for (kind,i,p),c in f.items():
                if kind:put(a,(i,p),c)
                else:
                    for j,b in grad[i].items():put(a,(j,p),c*b)
                    k+=c*qknown[i].get(p,F(0))
            arow.append(a);krow.append(k);active.update(a)
        adrows.append(arow);knownfinal.append(krow)
        if not any(arow) and not any(native[number]):constant_only.append({'index':number,'known_three_class_vector':list(map(str,krow))})
        final_support+=sum(map(len,arow))+sum(bool(k) for k in krow)
    live=saved_support+sum(len(grad[i])+len(first[i])+len(qknown[i])+bool(val[i]) for i in needed)+final_support
    assert live<CAP
    # Function source types/order are recorded literally. Composite C0 roots
    # cannot be assigned uncached variable gradients by treating h as a free
    # coordinate; their active second coefficient remains a genuine source gate.
    prefix={x['node']:x for x in A['all35_ordered_CA_prefix_firsts']};gover={i for i,p in A['all8_g_quotient_overrides']}
    inventory=[]
    for i in sorted({i for i,p in active}):
        node=H['DAG_scalar_nodes'][i];ps=sorted(p for j,p in active if j==i)
        kind=node[0];order=5 if i in gover or kind=='g' else 3 if i in prefix or kind=='r' else None
        inventory.append({'actual_C0_node':i,'literal_source_node':node,'active_pairs':[pairs[p] for p in ps],
          'is_original_B_atom79':i in atom,'is_saved_ordered_prefix':i in prefix,
          'is_saved_g_override':i in gover,'underlying_numerator_order_upper_if_named_quotient':order,
          'composite_source_second_gate_retained':order is None})
    dump(out,'phase2-atom-second-adjoint.json',{'pairs':pairs,'620_rows':[
      {'index':i,'three_atom_or_opaque_P_second_functionals':list(map(sparse,adrows[i])),
       'three_native_X2_functionals':list(map(sparse,native[i])),
       'known_first_product_and_Xp_correction_three_vector':list(map(str,knownfinal[i]))}
      for i in range(620)],'active_actual_C0_function_inventory':inventory,
      'old_only_combinations':oldonly,'entirely_unknown_free_combinations':constant_only,
      'known_nonzero_class_vectors_without_source2_input':[x for x in constant_only if any(F(c) for c in x['known_three_class_vector'])],
      'new_source_second_values_computed':False,'expanded_live_QQ_support':live})
    dump(out,'oldB-AD-fixed-node-values-and-gradients.json',{'all_needed_nodes':sorted(needed),
      'atom_source_sha256':sha(ROOT/paths['atoms']),
      'node_data':[[i,str(val[i]),encode(grad[i]),encode(first[i]),encode(qknown[i])] for i in sorted(needed)],
      'atom2_zero_used_only_to_split_known_first_products':True})
    nonzero=[x for x in constant_only if any(F(c) for c in x['known_three_class_vector'])]
    dump(out,'phase.json',{'phase':'FINAL','active_atom_or_P_roots':len(inventory),'all_unknown_free_combos':len(constant_only),'known_nonzero_class_vectors':len(nonzero)})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_combinations':620,'folded_root_QQ_support':saved_support,
      'typed_AD_nodes':len(needed),'active_atom_or_opaque_P_function_count':len(inventory),
      'active_atom_pair_count':len(active),'final_functional_QQ_support':final_support,'expanded_live_QQ_support':live,
      'unknown_free_combination_count':len(constant_only),'unknown_free_nonzero_class_count':len(nonzero),
      'all_P8_roots_retained':True,'no_uncached_functionvalues_gradients_sourcejets_rank_GB':True,
      'T_length_and_depth_not_computed':True,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','active_atom_or_opaque_P_function_count','active_atom_pair_count','unknown_free_combination_count','unknown_free_nonzero_class_count','expanded_live_QQ_support','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
