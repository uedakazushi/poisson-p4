#!/usr/bin/env python3
"""Preflight / existing-cache ordinary source3 evaluation before any source5.

All missing g2 stay symbolic via the frozen sensitivity matrix. We only
evaluate its annihilating combinations as actual relations.
"""
import argparse,collections,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled,sparse_rank
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_PURE619_CACHED_THIRD_ECONOMY_20261104_v63'
ZERO=(F(0),)*4
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def mul(a,b):return tuple(sum((a[i]*b[k-i] for i in range(k+1)),F(0)) for k in range(4))
def inv(a):
    assert a[0];q=[1/a[0]]
    for k in range(1,4):q.append(-sum((a[i]*q[k-i] for i in range(1,k+1)),F(0))/a[0])
    assert mul(a,tuple(q))==(F(1),F(0),F(0),F(0));return tuple(q)
def bit(q):return max(q.numerator.bit_length(),q.denominator.bit_length())
def main(out,evaluate):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'functional':'reports/AMBIENT_B18_T_CONSTANT_QUADRATIC_KERNEL_R3_ADJOINT_20261104_v63/phase2-all619-pure-cubic-source-functionals.json',
      'missing':'reports/AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63/actual619-missing-quotient-jet-sensitivities.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'joint':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json',
      'joint_peer':'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63_REVIEW/report.json',
      'own4':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
      'own4_peer':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63_REVIEW/report.json',
      'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
      'copy2':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
      'hfirst':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
      'oldAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'composites':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json',
      'native':'reports/AMBIENT_B18_T_NATIVE_S2_B11_20261104_v63/actual-native128-effective-S2-B11.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'all_cache_source_hashes_persisted','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};S,MS,H,J,JP,K,KP,CP,CP2,FJ,AF,AD,CV,N=[D[k] for k in paths]
    assert JP['all_passed'] and KP['all_passed'] and not MS['all_third_unknown_nodes_remaining']
    nodes=H['DAG_scalar_nodes'];common={}
    for i,p in J['all105_unknown_jets']:
        q=[F(0)]*4
        for m,c in zip(J['monomials'],p):q[sum(m)]+=F(c)
        common[i]=tuple(q)
    assert len(common)==105
    cb=K['families']['CB'];copybind=dict(CP['explicit57_source_copy_bindings']);assert set(copybind)==set(cb['actual_unknown57_order'])
    yroots=[i for i,p in CP2['actual_YB_last2_from_saved105']]
    eq105=J['all105_equation_original_node_order'];copyeq=CP['all57_copied_equation_indices']
    opaque={}
    for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
        for i,p in FJ[key].items():opaque[int(i)]=(F(p[0]),F(p[8])+F(p[9]),F(0),F(0))
    for p in AF['all35_ordered_CA_prefix_firsts']:
        q=rec(p['first9_bar']);opaque[p['node']]=(q.get(0,F(0)),q.get(8,F(0)),F(0),F(0))
    for i,p in AF['all8_g_quotient_overrides']:
        q=rec(p);opaque[i]=(q.get(0,F(0)),q.get(8,F(0)),F(0),F(0))
    proots={i for r in S['all619'] for key in ['rawP8_last_cubic_roots','Xp_second_r2_roots'] for i,c in r[key]}
    def closure(roots,stops):
        seen=set();todo=list(roots)
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i)
            if i in stops:continue
            n=nodes[i]
            if n[0] in ['a','m']:todo.extend(n[1:3])
            elif n[0] in ['n','i']:todo.append(n[1])
            else:assert n[0]=='c',('MISSING_CACHE_SOURCE_LEAF',i,n)
        return sorted(seen)
    def live_profile(order,roots,stops):
        uses=collections.Counter()
        for i in order:
            if i in stops:continue
            n=nodes[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
        uses.update(roots);live=set();peak=0
        for i in order:
            live.add(i);peak=max(peak,len(live))
            if i not in stops:
                n=nodes[i]
                for j in n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else []:
                    uses[j]-=1
                    if not uses[j]:live.discard(j)
        return peak
    roots1=eq105+yroots;order1=closure(roots1,set(common));p1=live_profile(order1,roots1,set(common))
    finalstops=set(common)|set(copybind.values())|set(opaque)|set(H['actual_C0_basis8_scalar_nodes'])
    roots2=copyeq+sorted(proots);order2=closure(roots2,finalstops);p2=live_profile(order2,roots2,finalstops)
    bound=4*(max(p1,p2)+len(common)+len(copybind)+len(opaque)+len(roots2))+20000
    assert bound<150000
    pre={'stage1_original105_equations_and_YB_tail':{'source_DAG_nodes':len(order1),'peak_live_nodes':p1},
      'stage2_copy57_equations_and_rawP3':{'source_DAG_nodes':len(order2),'peak_live_nodes':p2},
      'series_coefficients_per_node':4,'all105_common_unknowns_source3_cached':True,
      'own_CB57_bivariate_source3_in_source4_cache':True,'no_new_source_jet_or_IFTsolve':True,
      'all_missing_g2_remain_symbolic_via_sensitivity_matrix':True,'TOTAL_live_QQ_slots_bound':bound,
      'cap':150000,'artifact_bytes_cap':32*1024*1024,'address_space_bytes_cap':2000000000,
      'cached_common_source3_max_bits':max(bit(c) for p in common.values() for c in p)}
    dump(out,'preflight.json',pre)
    if not evaluate:
        dump(out,'report.json',{'id':STEM+'_PREFLIGHT','all_passed':True,'actual_source3_composition_or_619_values_executed':False,
          'input_sha256':hashes,'preflight_sha256':sha(out/'preflight.json'),'code_sha256':sha(Path(__file__)),'seconds':time.monotonic()-start})
        print(json.dumps(pre));return
    def ev(order,roots,bind):
        uses=collections.Counter()
        for i in order:
            if i in bind:continue
            n=nodes[i];uses.update(n[1:3] if n[0] in ['a','m'] else n[1:] if n[0] in ['n','i'] else [])
        uses.update(roots);vals={}
        for i in order:
            n=nodes[i];op=n[0]
            if i in bind:v=bind[i]
            elif op=='c':v=(F(n[1],n[2]),F(0),F(0),F(0))
            elif op=='a':v=tuple(a+b for a,b in zip(vals[n[1]],vals[n[2]]))
            elif op=='n':v=tuple(-a for a in vals[n[1]])
            elif op=='m':v=mul(vals[n[1]],vals[n[2]])
            elif op=='i':v=inv(vals[n[1]])
            else:raise AssertionError((i,n))
            vals[i]=v
            if i not in bind:
                for j in n[1:3] if op in ['a','m'] else n[1:] if op in ['n','i'] else []:
                    uses[j]-=1
                    if not uses[j]:vals.pop(j,None)
        return {i:vals[i] for i in roots}
    v1=ev(order1,roots1,common)
    for j,i in enumerate(eq105):assert v1[i]==((F(0),F(1),F(0),F(0)) if j in [55,56] else ZERO)
    ys=[v1[i] for i in yroots];copyvalues={}
    for original,p in cb['all57_source4_jets']:
        q=ZERO
        for m,c in zip(K['monomials_bivariate_through4'],p):
            if sum(m)>3 or not F(c):continue
            term=(F(c),F(0),F(0),F(0))
            for j,k in enumerate(m):
                for unused in range(k):term=mul(term,ys[j])
            q=tuple(a+b for a,b in zip(q,term))
        copyvalues[copybind[original]]=q
    finalbind=common|copyvalues|opaque
    gAcoeff=[F(0)]*4
    assert J['g_CA_all10_coefficients']==H['cached_actual_CA_gA_all_degree3_coefficients']
    for m,c in zip(J['monomials'],J['g_CA_all10_coefficients']):gAcoeff[sum(m)]+=F(c)
    for j,i in enumerate(H['actual_C0_basis8_scalar_nodes']):finalbind[i]=tuple(gAcoeff) if j==7 else ZERO
    v2=ev(order2,roots2,finalbind);assert all(v2[i]==ZERO for i in copyeq)
    dump(out,'phase1-existing-source3-and-copy-composition.json',{'actual_common105_rline_through3':[[i,list(map(str,p))] for i,p in sorted(common.items())],
      'actual_YB_last2_through3':[[i,list(map(str,v1[i]))] for i in yroots],
      'all105_source_equations_through3_verified':True,'copy57_unknowns_through3':[[i,list(map(str,p))] for i,p in sorted(copyvalues.items())],
      'all57_copy_equations_through3_verified':True,'actual_gA_rline_through3':list(map(str,gAcoeff)),
      'raw_source_roots_with_missing_quotient2_3_set_zero_ONLY_for_known_part':[[i,list(map(str,v2[i]))] for i in sorted(proots)]})
    td={i:(rec(g),rec(q)) for i,v,g,f,q in AD['node_data']};cv={i:(rec(v),rec(g)) for i,v,g in CV['all108_composite_values_and_opaque_second_gradients']}
    native=[{row:rec(v) for row,v in col} for col in N['source_columns128_rows48_sparse_B11_effective_second']]
    allowed=[(0,7),(1,7),(4,7),(5,7),(7,7)];pairB=dict(zip(allowed,[7,8,9,10,2]));pairAD={p:j for j,p in enumerate(allowed)}
    pairs=[tuple(p) for p in S['pairs36_order']]
    def x2value(key):
        p=pairs[key[-1]];assert p in pairB;b=pairB[p]
        if key[0]==2:return native[key[1]].get(key[2],{}).get(b,F(0))
        if key[0]==1:assert not cv[key[1]][1];return cv[key[1]][0].get(b,F(0))
        g,q=td[key[1]];assert all(not cv[i][1] for i in g)
        return q.get(pairAD[p],F(0))+sum((a*cv[i][0].get(b,F(0)) for i,a in g.items()),F(0))
    knownpart=[]
    for row in S['all619']:
        q=F(row['known_Xp1_a1_contraction'])
        q+=sum((F(c)*v2[i][3] for i,c in row['rawP8_last_cubic_roots']),F(0))
        q+=sum((F(c)*v2[i][2] for i,c in row['Xp_second_r2_roots']),F(0))
        q+=sum((F(c)*x2value(key) for key,c in row['quadratic_source_root_correction']),F(0))
        knownpart.append(q)
    missing_nodes=[i for i,n in MS['all_second_unknown_nodes_remaining']];sens=[]
    for row in MS['all619']:
        c=rec(row['opaque_second_r2_sensitivity']);assert not row['opaque_third_r3_sensitivity']
        sens.append({missing_nodes.index(i):a for i,a in c.items()})
    rank=sparse_rank(sens,8);norm=[];free=[]
    for j,op in enumerate(rank['fixed_elimination_receipts']):
        c={j:F(1)}
        for p,a in op['subtract']:addscaled(c,norm[p],-F(a))
        if 'pivot_row' in op:norm.append({i:a/F(op['divide']) for i,a in c.items()})
        else:
            check={}
            for i,a in c.items():addscaled(check,sens[i],a)
            assert not check
            v=sum((a*knownpart[i] for i,a in c.items()),F(0));free.append({'source619':j,'combination619':encode(c),'actual_missing_independent_pure_r3_value':str(v)})
    nonzero=[r for r in free if F(r['actual_missing_independent_pure_r3_value'])]
    dump(out,'phase2-known-parts-and-missing-independent-projection.json',{'all619_known_part_missing_coefficients_split':list(map(str,knownpart)),
      'missing_second_primitive_nodes':missing_nodes,'missing_sensitivity_rank':rank,'all_missing_annihilator_combinations':free,
      'nonzero_actual_pure_r3_relation_witnesses':nonzero,'missing_g2_values_computed':False,
      'all_source3_values_are_existing_cache_compositions':True})
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL','missing_sensitivity_rank':rank['rank'],'nonzero_missing_independent_relations':len(nonzero)})
    report={'id':STEM,'all_passed':True,'all105_and_copy57_source_equations_verified_through3':True,
      'missing_second_sensitivity_rank':rank['rank'],'missing_independent_relation_count':len(free),
      'nonzero_actual_missing_independent_cubic_relation_count':len(nonzero),
      'm3_T_zero_proved_primary':bool(nonzero),'T_total_length_if_proved':10 if nonzero else None,
      'no_new_source3_or5_jet_or_g2_solve':True,'independent_source_review_status':'PENDING',
      'input_sha256':hashes,'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','missing_second_sensitivity_rank','missing_independent_relation_count','nonzero_actual_missing_independent_cubic_relation_count','T_total_length_if_proved','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--evaluate',action='store_true');a=p.parse_args();main(a.output,a.evaluate)
