#!/usr/bin/env python3
"""Exact copied-CB57 leaf binding and existing-cache composite second AD.

No source solve or source jet: original105 H20 values and the already saved
CB57 degree-two inverse in its OWN YB coordinates are the only function data.
"""
import argparse,collections,itertools,json,resource,signal,sys,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63'
CAP=150000
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def put(v,k,a):
    b=v.get(k,F(0))+a
    if b:v[k]=b
    else:v.pop(k,None)
def enc_keys(v):return [[list(k) if isinstance(k,tuple) else k,str(a)] for k,a in sorted(v.items()) if a]

def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));sys.set_int_max_str_digits(0)
    out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'adj':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/phase2-atom-second-adjoint.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'full':'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json',
      'source':'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json',
      'first':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'cache':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
      'H20':'reports/AMBIENT_B18_CONDUCTOR_PRODUCT_SPECIAL_FIBER20_20261104_v63/actual-special-fiber20-coefficients.json',
      'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'source_hashes_persisted_before_label_binding','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};A,C,Full,S,J,Cache,H,AF=[D[k] for k in paths]
    nodes=C['DAG_scalar_nodes'];original=S['DAG_scalar_nodes'];assert nodes[:len(original)]==original
    inst=Full['all16_pointed_unit57_inverse_instances'][14];assert inst['retained_Y_indices']==[7,8]
    bind=dict(inst['bindings']);CB=J['both57_source2']['CB'];assert set(bind)==set(CB['unknown57'])
    assert CB['equations57']==Full['original_unit48_equations_CB']+Full['CB_coordinate_Y9_nodes']
    assert original[11]==['v','CA_trace_p'] and S['families']['CA']['p_quadratic_relation']['s']==11
    assert 11 in Cache['joint105_unknown_order']
    # Reconstruct the copied expression indices using ONLY fixed source DAG
    # syntax and the serialized binding, not a name comparison or a builder.
    table={tuple(n):i for i,n in enumerate(nodes)};zero=table[('c',0,1)];one=table[('c',1,1)]
    def cc(q):return table[('c',q.numerator,q.denominator)]
    def neg(a):
        n=nodes[a]
        if a==zero:return a
        if n[0]=='c':return cc(-F(n[1],n[2]))
        if n[0]=='n':return n[1]
        return table[('n',a)]
    def plus(a,b):
        if a==zero:return b
        if b==zero:return a
        if nodes[a]==['n',b] or nodes[b]==['n',a]:return zero
        if nodes[a][0]==nodes[b][0]=='c':return cc(F(nodes[a][1],nodes[a][2])+F(nodes[b][1],nodes[b][2]))
        return table[('a',min(a,b),max(a,b))]
    def times(a,b):
        if a==zero or b==zero:return zero
        if a==one:return b
        if b==one:return a
        if nodes[a][0]==nodes[b][0]=='c':return cc(F(nodes[a][1],nodes[a][2])*F(nodes[b][1],nodes[b][2]))
        return table[('m',min(a,b),max(a,b))]
    def inverse(a):
        if nodes[a][0]=='c':return cc(1/F(nodes[a][1],nodes[a][2]))
        return table[('i',a)]
    need=set();todo=CB['equations57'][:]
    while todo:
        i=todo.pop()
        if i in need:continue
        need.add(i);n=original[i]
        if n[0] in ['a','m']:todo.extend(n[1:3])
        elif n[0] in ['n','i']:todo.append(n[1])
    copied={}
    for i in sorted(need):
        n=original[i];op=n[0]
        if op=='v':copied[i]=bind[i]
        elif op=='c':copied[i]=i
        elif op=='a':copied[i]=plus(copied[n[1]],copied[n[2]])
        elif op=='m':copied[i]=times(copied[n[1]],copied[n[2]])
        elif op=='n':copied[i]=neg(copied[n[1]])
        elif op=='i':copied[i]=inverse(copied[n[1]])
        else:raise AssertionError(op)
    eq=copied_eq=[copied[i] for i in CB['equations57'][:48]]
    eq+= [plus(copied[i],neg(i if k in [7,8] else zero)) for k,i in enumerate(CB['F9_nodes'])]
    assert eq==inst['all57_pointed_equations']
    dump(out,'phase1-exact-leaf-label-transport.json',{'original_source_prefix_literal_identical':True,
      'instance15_retains_YB_indices':[7,8],'explicit57_source_copy_bindings':inst['bindings'],
      'all57_copied_equation_indices':eq,'all57_copied_equations_literal_match':True,
      'CA_trace_p_node11_definition':{'original_literal':original[11],'quadratic_rule':S['families']['CA']['p_quadratic_relation'],
       'cache105_unknown_index':Cache['joint105_unknown_order'].index(11)},
      'no_name_heuristic_or_new_source_inverse':True})
    # B11 is the degree<=2 directional quotient with r² and four zi*r kept.
    e=[0]*8;Bmons=[tuple(e)]
    def mon(**kw):
        x=[0]*8
        for k,v in kw.items():x[int(k)]=v
        return tuple(x)
    for k,p in [(7,1),(7,2),(0,1),(1,1),(4,1),(5,1)]:
        x=[0]*8;x[k]=p;Bmons.append(tuple(x))
    for k in [0,1,4,5]:x=[0]*8;x[k]=x[7]=1;Bmons.append(tuple(x))
    bm={m:i for i,m in enumerate(Bmons)};assert len(bm)==11
    mulindex={(i,j):bm.get(tuple(a+b for a,b in zip(m,n))) for i,m in enumerate(Bmons) for j,n in enumerate(Bmons)}
    def mul(a,b):
        z={}
        for i,c in a.items():
            for j,d in b.items():
                k=mulindex[i,j]
                if k is not None:put(z,k,c*d)
        return z
    def add(a,b,c=F(1)):
        z=dict(a);addscaled(z,b,c);return z
    def inv(a):
        c=a.get(0,F(0));assert c
        h={i:v for i,v in a.items() if i};hh=mul(h,h);z={0:1/c};addscaled(z,h,-1/c**2);addscaled(z,hh,1/c**3)
        assert mul(a,z)=={0:F(1)};return z
    vars8=[{bm[tuple(int(j==k) for j in range(8))]:F(1)} if tuple(int(j==k) for j in range(8)) in bm else {} for k in range(8)]
    images=[{bm[tuple(m)]:F(1)} if tuple(m) in bm else {} for m in H['quotient_basis_monomials']]
    def phi(p):
        v={}
        for i,a in p.items():addscaled(v,images[i],a)
        return v
    for k,action in enumerate(H['multiplication8_column_sparse']):
        for j,p in enumerate(action):assert mul(vars8[k],images[j])==phi(rec(p))
    assert Cache['A20_basis']==H['quotient_basis_monomials']
    cached={i:phi(rec(p)) for i,p in Cache['all_requested_original_source_root_values']}
    cached.update({i:phi(rec(p)) for i,p in zip(Cache['joint105_unknown_order'],Cache['all105_unknown_final_values'])})
    # Source variables are fixed before evaluation; no new57/105 solve.
    def poly_closure(roots,stops):
        needed=set();todo=list(roots)
        while todo:
            i=todo.pop()
            if i in needed:continue
            needed.add(i)
            if i in stops:continue
            n=nodes[i]
            if n[0] in ['a','m']:todo.extend(n[1:3])
            elif n[0] in ['n','i']:todo.append(n[1])
            elif n[0]!='c':raise AssertionError(('MISSING_CACHED_LEAF',i,n))
        return needed
    def eval_known(roots,cache):
        ns=poly_closure(roots,set(cache));uses=collections.Counter()
        for i in ns:
            if i in cache:continue
            n=nodes[i]
            if n[0] in ['a','m']:uses.update(n[1:3])
            elif n[0] in ['n','i']:uses.update(n[1:])
        uses.update(roots);live=set();peak=0
        for i in sorted(ns):
            live.add(i);peak=max(peak,len(live));n=nodes[i]
            if i not in cache:
                if n[0] in ['a','m']:children=n[1:3]
                elif n[0] in ['n','i']:children=n[1:]
                else:children=[]
                for j in children:
                    uses[j]-=1
                    if not uses[j]:live.discard(j)
        assert 11*peak+len(cache)*11<75000
        vals={};uses=collections.Counter()
        for i in ns:
            if i in cache:continue
            n=nodes[i]
            if n[0] in ['a','m']:uses.update(n[1:3])
            elif n[0] in ['n','i']:uses.update(n[1:])
        uses.update(roots)
        for i in sorted(ns):
            n=nodes[i];o=n[0]
            if i in cache:v=cache[i]
            elif o=='c':v={0:F(n[1],n[2])} if n[1] else {}
            elif o=='a':v=add(vals[n[1]],vals[n[2]])
            elif o=='m':v=mul(vals[n[1]],vals[n[2]])
            elif o=='n':v={j:-c for j,c in vals[n[1]].items()}
            elif o=='i':v=inv(vals[n[1]])
            else:raise AssertionError((i,n))
            vals[i]=v
            if i not in cache:
                children=n[1:] if o in ['n','i'] else n[1:3] if o in ['a','m'] else []
                for j in children:
                    uses[j]-=1
                    if not uses[j]:vals.pop(j,None)
        return {i:vals[i] for i in roots},len(ns),peak
    yroots=CB['F9_nodes'][7:9];yb,ynodes,ypeak=eval_known(yroots,cached)
    args=[{} for _ in range(7)]+[yb[i] for i in yroots]
    ownmons=[tuple(m) for m in J['source2_monomials']]
    copyvalues={}
    for original_var,p in CB['all57_source2jets']:
        v={}
        for monomial,c in zip(ownmons,p):
            c=F(c)
            if not c:continue
            term={0:c}
            for k in monomial:term=mul(term,args[k])
            addscaled(v,term,F(1))
        copyvalues[bind[original_var]]=v
    copiedcache=cached|copyvalues
    equations_values,enodes,epeak=eval_known(eq,copiedcache);assert all(not v for v in equations_values.values())
    dump(out,'phase2-existing-cache-B11-copy-values.json',{'B11_monomial_order':Bmons,'H20_action_intertwining_identities160':True,
      'actual_YB_last2_from_saved105':[[i,encode(yb[i])] for i in yroots],
      'all57_copied_unknown_values_through2':[[i,encode(v)] for i,v in sorted(copyvalues.items())],
      'all57_copied_equations_values':[[i,encode(equations_values[i])] for i in eq],
      'YB_source_closure_nodes':ynodes,'YB_peak_live_nodes':ypeak,'copied_equation_closure_nodes':enodes,
      'copied_equation_peak_live_nodes':epeak,'all57_equations_exactly_zero_in_B11':True,
      'copy_leaf_values_are_saved_CB57_2jets_composed_with_actual_YB7_8':True})
    dump(out,'phase.json',{'phase':'all_copy57_bindings_and_B11_equations_verified'})
    # Exact first caches for honest regular divisions. Their second values
    # are NOT read from H20 or assigned zero as an actual-source assertion.
    unknownfirst={}
    def full9bar(p):return [F(x) for x in p[:8]]+[F(p[8])+F(p[9])]
    for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
        for i,p in J[key].items():unknownfirst[int(i)]=full9bar(p)
    for p in AF['all35_ordered_CA_prefix_firsts']:
        q=rec(p['first9_bar']);unknownfirst[p['node']]=[q.get(i,F(0)) for i in range(9)]
    for i,p in AF['all8_g_quotient_overrides']:
        q=rec(p);unknownfirst[i]=[q.get(k,F(0)) for k in range(9)]
    # The actual C0 coordinates have literal value z_i; gA starts in degree3
    # by the pinned source3 certificate and is therefore zero in B11.
    for k,i in enumerate(C['actual_C0_basis8_scalar_nodes']):copiedcache[i]=vars8[k] if k<7 else {}
    active=[x['actual_C0_node'] for x in A['active_actual_C0_function_inventory']]
    stops=set(copiedcache)|set(unknownfirst);closure=poly_closure(active,stops)
    opaque=[i for i in closure if i in unknownfirst and i not in copiedcache]
    # Forecast only possible sets of opaque second inputs, no invented
    # variable gradients. All coefficients arise from true polynomial rules.
    val={};grad={};inverses=[]
    structural_support_bound=len(closure)*(11+len(opaque))+15000
    # Actual expected closure is small; use a structural dependency bound.
    deps={}
    for i in sorted(closure):
        if i in copiedcache:deps[i]=set()
        elif i in unknownfirst:deps[i]={i}
        else:
            n=nodes[i]
            deps[i]=set().union(*(deps[j] for j in (n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])))
    forecast=11*len(closure)+sum(map(len,deps.values()))+15000
    dump(out,'preflight.json',{'composite_closure_nodes':len(closure),'opaque_cached_regular_division_leaves':len(opaque),
      'known_B11_slots_upper':11*len(closure),'opaque_gradient_slots_upper':sum(map(len,deps.values())),
      'TOTAL_forecast':forecast,'cap':CAP,'no_uncached_function_leaf':True})
    assert forecast<CAP
    for i in sorted(closure):
        n=nodes[i];o=n[0]
        if i in copiedcache:v=copiedcache[i];g={}
        elif i in unknownfirst:
            p=unknownfirst[i];v={0:p[0]} if p[0] else {}
            for k in range(8):addscaled(v,vars8[k],p[k+1])
            g={i:F(1)}
        elif o=='c':v={0:F(n[1],n[2])} if n[1] else {};g={}
        elif o=='a':v=add(val[n[1]],val[n[2]]);g=add(grad[n[1]],grad[n[2]])
        elif o=='n':v={j:-c for j,c in val[n[1]].items()};g={j:-c for j,c in grad[n[1]].items()}
        elif o=='m':
            v=mul(val[n[1]],val[n[2]]);g={};addscaled(g,grad[n[1]],val[n[2]].get(0,F(0)));addscaled(g,grad[n[2]],val[n[1]].get(0,F(0)))
        elif o=='i':
            v=inv(val[n[1]]);g={j:-c/val[n[1]][0]**2 for j,c in grad[n[1]].items()};inverses.append([i,n[1],encode(val[n[1]]),encode(v)])
        else:raise AssertionError(('MISSING_CACHED_COMPOSITE',i,n))
        assert set(g)<=deps[i];val[i]=v;grad[i]=g
    # Compose with all620*3 original atom/root adjoints. B11 pair indices are
    # [7,8,9,10,2] for [z0r,z1r,z4r,z5r,r²].
    pair_to_B=[7,8,9,10,2];rows=[];active_quotients=set();new_nz=0;known_nonzero=[]
    for row in A['620_rows']:
        af=[];kv=[]
        for t,f in enumerate(row['three_atom_or_opaque_P_second_functionals']):
            g={};known=F(row['known_first_product_and_Xp_correction_three_vector'][t])
            for (i,p),c in f:
                c=F(c);known+=c*val[i].get(pair_to_B[p],F(0))
                for j,a in grad[i].items():put(g,(j,p),c*a)
            af.append(g);kv.append(known);active_quotients.update(g)
        item={'index':row['index'],'three_unknown_division_second_functionals':list(map(enc_keys,af)),
          'three_known_oldB_P_and_pivot_correction_values':list(map(str,kv)),
          'three_native_X2_functionals':row['three_native_X2_functionals']}
        rows.append(item);new_nz+=sum(map(len,af))+sum(bool(x) for x in kv)
        if not any(af) and not any(row['three_native_X2_functionals']) and any(kv):known_nonzero.append(item['index'])
    total=sum(len(val[i])+len(grad[i]) for i in closure)+new_nz
    assert total<CAP
    dump(out,'phase3-composite-second-source-functional.json',{'B11_monomial_order':Bmons,
      'all108_composite_values_and_opaque_second_gradients':[[i,encode(val[i]),encode(grad[i])] for i in active],
      'all_unit_inverse_checks':inverses,'620_projected_source_functionals':rows,
      'active_honest_division_second_pairs':[[i,p,nodes[i]] for i,p in sorted(active_quotients)],
      'known_nonzero_classes_without_unknown_source':known_nonzero,'new_source_jet_or_solve':False,
      'all_unknown_second_values_zero_only_in_known_part_split':True,'expanded_live_QQ_support':total})
    dump(out,'phase.json',{'phase':'FINAL','active_division_second_pairs':len(active_quotients),
      'known_nonzero_without_unknown_source':known_nonzero})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'all48_FULLINV15_leaves_bound_by_exact_original_ID':True,
      'all57_copy_equations_structural_and_B11_pass':True,'CA_trace_p_literal_node11_bound':True,
      'H20_to_B11_all160_action_identities':True,'active_composite_count':len(active),
      'active_division_second_pair_count':len(active_quotients),'active_division_second_function_count':len({i for i,p in active_quotients}),
      'known_nonzero_class_without_unknown_source_count':len(known_nonzero),
      'actual_composite_expanded_QQ_support':total,'no_new_functionvalues_sourcejet_IFTsolve_GB_point':True,
      'T_length_and_depth_not_computed':True,'input_sha256':hashes,'code_sha256':sha(Path(__file__)),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','active_division_second_pair_count','active_division_second_function_count','known_nonzero_class_without_unknown_source_count','actual_composite_expanded_QQ_support','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
