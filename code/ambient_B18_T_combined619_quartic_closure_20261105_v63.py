#!/usr/bin/env python3
"""Combined619 pure-quartic cache closure with honest-second LIN/QUAD symbols."""
import argparse,collections,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
from ambient_B18_T_quartic_symbolic_tools_20261104_v63 import MON,INDEX,pa,ja,js,jm,ji,scalar,known,count,encode_poly,encode_jet
from ambient_B18_T_joint135_functional_tools_20261104_v63 import MON as MON25,Finite
import ambient_B18_comparison_X_origin_20261104_v63 as H
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_T_COMBINED619_QUARTIC_CLOSURE_20261105_v63'
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    def dump(n,x):
        (out/n).write_text(json.dumps(x,separators=(',',':'))+'\n');assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    paths={
      'high':'reports/AMBIENT_B18_T_CONSTANT619_QUARTIC_HIGH_ORDER_ADJOINT_20261104_v63/all619-high-order-root-adjoints.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'chain':'reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json',
      'source2':'reports/AMBIENT_B18_T_SOURCE2_CACHE_COMPOSITION_20261104_v63/phase1-actual105-source2-and-copy.json',
      'own':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'own4':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
      'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase1-exact-leaf-label-transport.json',
      'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
      'g2':'reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_20261104_v63/actual-gdiff-second-and-complete619-pure-cubic-values.json',
      'g3':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63/actual-eight-g-third-values.json',
      'CA6':'reports/AMBIENT_B18_T_GDIFF_THIRD_UNIVARIATE6_20261104_v63/CA-source6-and-regular-quotient-third.json',
      'phi4':'reports/AMBIENT_B18_T_SOURCE4_ADJOINT9_20261104_v63/actual-nine-fourth-source-functionals.json',
      'highJac':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63/phase4-all619-high-old-source-sensitivities.json',
      'highJac_peer':'reports/AMBIENT_B18_T_HIGH_OLD_JACOBIAN_PREFLIGHT_20261104_v63_REVIEW/report.json',
      'implicit':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase4-selected-joint135-third-adjoint-values.json',
      'implicit_peer':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63_REVIEW/report.json',
      'response':'reports/AMBIENT_B18_T_JOINT135_SELECTED_FUNCTIONALS_20261104_v63/phase2-known-third-RHS-and-proper-response.json',
      'N2':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/phase2-actual-N-B11.json',
      'Soperator':'reports/AMBIENT_B18_COMPLETE_NORMALIZER_SYZYGY_OPERATOR_20261104_v63/actual-complete-Res4-syzygy-operator.json',
      'B11':'reports/AMBIENT_B18_T_NATIVE_B11_SOURCE_PREFLIGHT_20261104_v63/actual-B11-source-and-gauge-preflight.json',
      'lower':'reports/AMBIENT_B18_T_NATIVE_ZZ_GAUGE_CACHED_20261104_v63/all619-lower-contribution-values.json',
      'native_partial':'reports/AMBIENT_B18_T_NATIVE_N2_CACHED_20261104_v63/actual-L2-requested-values-and619-native-contributions.json'}
    hs={p:sha(ROOT/p) for p in paths.values()};dump('phase.json',{'phase':'combined_cache_sources_bound_before_symbolic_values','input_sha256':hs})
    D={k:read(ROOT/p) for k,p in paths.items()};HI,C,B,S,J,O4,CP,AF,G2,G3,CA6,P4,HJ,HJP,IMP,IPP,RESP,N2,OP,B11,LOW,NAT=[D[k] for k in paths]
    assert HJP['all_passed'] and IPP['all_passed']
    assert sha(ROOT/paths['implicit_peer'])=='554d2b46684dc1f044f8cf98f658db34e896f7c7f43fc50c785575b54cf4fc36'
    nodes=C['DAG_scalar_nodes'];ops=B['operator_scalar_nodes'];mons2=list(__import__('itertools').combinations_with_replacement(range(8),2));mons3=list(__import__('itertools').combinations_with_replacement(range(8),3))
    phi4={i:F(a) for i,a in P4['actual_nine_fourth_values']};g2={i:F(a) for i,a in G2['actual_eight_regular_g_second_values']};g3={i:F(a) for i,a in G3['actual_eight_common_regular_g_third_values']}
    assert len(phi4)==9 and len(g3)==8
    bind={i:known(v) for i,v in S['all105_actual_source2_values']}
    for i,a in phi4.items():
        if a:bind[i][INDEX[(7,7,7,7)]]={():a}
    cb=J['both57_source2']['CB'];callbacks={new:old for old,new in CP['explicit57_source_copy_bindings']};yb=cb['F9_nodes'][7:9]
    cb4=O4['families']['CB'];own4={i:list(map(F,p)) for i,p in cb4['all57_source4_jets']};mon4=O4['monomials_bivariate_through4'];assert set(own4)==set(callbacks.values())
    opaque={}
    def first(p):
        p=list(map(F,p));v=scalar(p[0])
        for k in [0,1,4,5,7]:
            a=p[8]+p[9] if k==7 else p[1+k]
            if a:v[INDEX[(k,)]]={():a}
        return v
    for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
        for i,p in J[key].items():opaque[int(i)]=first(p)
    for row in AF['all35_ordered_CA_prefix_firsts']:
        p=rec(row['first9_bar']);v=scalar(p.get(0,F(0)))
        for k in [0,1,4,5,7]:
            if p.get(k+1):v[INDEX[(k,)]]={():p[k+1]}
        opaque[row['node']]=v
    for i,p in AF['all8_g_quotient_overrides']:
        p=rec(p);v=scalar(p.get(0,F(0)))
        for k in [0,1,4,5,7]:
            if p.get(k+1):v[INDEX[(k,)]]={():p[k+1]}
        opaque[i]=v
    parameter={};registry=[]
    def parameter_for(leaf,m):
        k=(leaf,m)
        if k not in parameter:parameter[k]=len(registry);registry.append({'parameter':parameter[k],'actual_honest_source_node':leaf,'literal_node':nodes[leaf],'ordinary_second_monomial':m})
        return parameter[k]
    def opaque_value(i):
        v={k:dict(p) for k,p in opaque[i].items()}
        for m in mons2:
            if m not in INDEX:continue
            b=INDEX[m]
            if m==(7,7) and i in g2:
                if g2[i]:v[b]={():g2[i]}
            else:v[b]={(parameter_for(i,m),):F(1)}
        if i in g3 and g3[i]:v[INDEX[(7,7,7)]]={():g3[i]}
        return v
    coords={}
    for k,i in enumerate(C['actual_C0_basis8_scalar_nodes']):
        if k<7:coords[i]={INDEX[(k,)]:{():F(1)}} if (k,) in INDEX else {}
        else:
            gp=list(map(F,CA6['g_through6']));assert gp[3]==F(125,24) and not gp[4];coords[i]={INDEX[(7,7,7)]:{():gp[3]}}
    rootC=set();rootT=set()
    for row in HI['all619']:
        rootC.update(i for key in ['rawP55_pure_r4_roots','Xp55_pure_r3_times_a1r_roots'] for i,c in row[key])
        for key in ['actual_X3_high_order_source_functional','actual_Xp2_times_a1_source_functional']:
            for k,a in row[key]:
                if k[0]==0:rootT.add(k[1])
                elif k[0]==1:rootC.add(k[1])
    needT=set();todo=list(rootT)
    while todo:
        i=todo.pop()
        if i in needT:continue
        needT.add(i);n=ops[i]
        if n[0]=='sum':todo.extend(j for j,a,b in n[1])
        elif n[0]=='product':todo.extend(n[3])
        elif n[0]=='actual_C0_scalar_function':rootC.add(n[1])
        else:assert n[0]=='c'
    def child(i):
        if i in bind or i in coords or i in opaque:return []
        if i in callbacks:return yb
        n=nodes[i]
        if n[0] in ['a','m','n','i']:return n[1:]
        assert n[0]=='c',('MISSING_LITERAL_SOURCE',i,n);return []
    needC=set();todo=list(rootC)
    while todo:
        i=todo.pop()
        if i in needC:continue
        needC.add(i);todo.extend(child(i))
    dump('preflight.json',{'kind':'EXISTING_CACHE_COMBINED_QUARTIC_EVALUATION','C0_root_count':len(rootC),'C0_closure':len(needC),'typed_roots':len(rootT),'typed_closure':len(needT),
      'ordinary_honest_seconds':'ALL retained base20 secondmonomials; LIN andQUADRATIC products explicit',
      'old_source3_known_part':'allcommon105thirds zero; selected135 adds theiractualprojection exactlyonce',
      'source4_known_part':'nineactualtargets; allotherphysicalfourths annihilated bypinnedhighJac',
      'honest3_known_part':'eightactualgthirds; otherhonest3/4 annihilated bypinnedhighJac',
      'N3_omission':'all619 S0c0zero fromselected135; no arbitraryrowgauge',
      'cap_QQ':150000,'numeric_native_source_and_arithmetic_reserve_QQ':50000,'cap_bytes':32*1024*1024,'no_new_source_or_implicit_solve':True})
    uses=collections.Counter(j for i in needC for j in child(i));uses.update(rootC);cv={};live=0;peak=0
    for i in sorted(needC):
        n=nodes[i];o=n[0]
        if i in bind:v=bind[i]
        elif i in coords:v=coords[i]
        elif i in opaque:v=opaque_value(i)
        elif i in callbacks:
            v={};powers=[[scalar(1)] for k in range(2)]
            for k in range(2):
                for p in range(1,5):powers[k].append(jm(powers[k][-1],cv[yb[k]]))
            for m,a in zip(mon4,own4[callbacks[i]]):
                if a:v=ja(v,jm(powers[0][m[0]],powers[1][m[1]]),a)
        elif o=='c':v=scalar(F(n[1],n[2]))
        elif o=='a':v=ja(cv[n[1]],cv[n[2]])
        elif o=='n':v=js(cv[n[1]],F(-1))
        elif o=='m':v=jm(cv[n[1]],cv[n[2]])
        elif o=='i':v=ji(cv[n[1]])
        else:raise AssertionError((i,n))
        cv[i]=v;live+=count(v);peak=max(peak,live);assert peak+50000<150000,('OLD_SYMBOLIC_LIVE_CAP',i,peak)
        for j in child(i):
            uses[j]-=1
            if not uses[j]:live-=count(cv[j]);del cv[j]
    dump('phase1-C0-symbolic-quartic-roots.json',{'R0_monomials':MON,'honest_second_parameter_registry':registry,
      'all_C0_roots':[[i,encode_jet(cv[i])] for i in sorted(rootC)],'peak_live_QQ':peak,
      'copy57_own_bivariate_order4_composed_coherently':True})
    usesT=collections.Counter()
    for i in needT:
        n=ops[i]
        if n[0]=='sum':usesT.update(j for j,a,b in n[1])
        elif n[0]=='product':usesT.update(n[3])
    usesT.update(rootT);tv={};tlive=0;tpeak=0
    for i in sorted(needT):
        n=ops[i];o=n[0]
        if o=='c':v=scalar(F(n[1],n[2]))
        elif o=='actual_C0_scalar_function':v=cv[n[1]]
        elif o=='sum':
            v={}
            for j,a,b in n[1]:v=ja(v,tv[j],F(a,b))
        elif o=='product':
            v=scalar(F(n[1],n[2]))
            for j in n[3]:v=jm(v,tv[j])
        else:raise AssertionError(n)
        tv[i]=v;tlive+=count(v);tpeak=max(tpeak,tlive);assert peak+tpeak+50000<150000,('TYPED_SYMBOLIC_CAP',i,peak,tpeak)
        ch=[j for j,a,b in n[1]] if o=='sum' else n[3] if o=='product' else []
        for j in ch:
            usesT[j]-=1
            if not usesT[j]:tlive-=count(tv[j]);del tv[j]
    dump('phase2-typed-symbolic-quartic-roots.json',{'all_typed_roots':[[i,encode_jet(tv[i])] for i in sorted(rootT)],
      'peak_live_QQ':tpeak,'all_unknown_parameter_products_retained':True})
    # Numeric known native response only; unknown third inputs are phase4implicit.
    nu=[tuple(map(F,p)) for p in RESP['known_nu2_base25']];gs=[tuple(map(F,p)) for p in RESP['known_properG8_C100']];a=Finite(2,nu);a.install()
    b11=[tuple(i for i,n in enumerate(m) for _ in range(n)) for m in B11['B11_basis_monomials8']]
    nmat=[]
    for row in N2['N_C_B44']:
        nr=[]
        for p in row:
            v=[F(0)]*100
            for j,m in enumerate(b11):
                for t in range(4):v[4*MON25.index(m)+t]=F(p[4*j+t])
            nr.append(tuple(v))
        nmat.append(nr)
    sn={}
    for i,j,k,sign in OP['native_S12_by32']['koszul28_in_first8_rows']:sn[i,j]=a.scale(gs[k],F(sign))
    for i,j,c in OP['native_S12_by32']['free4_zero_source_columns']:sn[i,j]=a.c(c)
    native_cache={};native_j_done=set()
    all_native={tuple(k) for row in HI['all619'] for key in ['actual_X3_high_order_source_functional','actual_Xp2_times_a1_source_functional'] for k,c in row[key] if k[0]==2}
    req_by_j=collections.defaultdict(set)
    for k in all_native:req_by_j[k[1]].add(k)
    def native_value(k):
        j=k[1]
        if j not in native_j_done:
            jn,b=divmod(j,4);unit=[F(0)]*100;unit[b]=F(1)
            sourcecol={}
            for row in range(12):
                q=a.Z
                for source in range(12):
                    if (source,jn) in sn:q=a.add(q,a.mul(nmat[source][row],sn[source,jn]))
                sourcecol[row]=a.mul(q,tuple(unit))
            for label in req_by_j[j]:
                _,jj,src,order,m=label;mon=mons3[m] if order==3 else mons2[m];row,t=divmod(src,4)
                native_cache[label]=sourcecol[row][4*MON25.index(mon)+t]
            native_j_done.add(j)
        return native_cache[k]
    lower=list(map(F,LOW['all619_actual_lower_contribution_values']));implicit=[sum(map(F,row),F(0)) for row in IMP['all619_actual_implicit_third_functional_values_by5']]
    part=list(map(F,NAT['all619_native_L2_contribution_values']))
    priority=sorted(range(619),key=lambda i:(-(4*bool(lower[i])+2*bool(implicit[i])+bool(part[i])),i))
    records=[];unbound=[];witness=None;support=peak+tpeak
    for number in priority:
        row=HI['all619'][number];poly={}
        for root,c in row['rawP55_pure_r4_roots']:poly=pa(poly,cv[root].get(INDEX[(7,7,7,7)],{}),F(c))
        for root,c in row['Xp55_pure_r3_times_a1r_roots']:poly=pa(poly,cv[root].get(INDEX[(7,7,7)],{}),F(c))
        native=F(0)
        for key in ['actual_X3_high_order_source_functional','actual_Xp2_times_a1_source_functional']:
            for k,c in row[key]:
                k=tuple(k);c=F(c)
                if k[0]==2:native+=c*native_value(k);continue
                order=k[-2] if k[0]<2 else k[-2];m=k[-1];mon=mons3[m] if order==3 else mons2[m]
                value=(tv if k[0]==0 else cv)[k[1]].get(INDEX[mon],{});poly=pa(poly,value,c)
        symbols={m:c for m,c in poly.items() if m};constant=poly.get((),F(0));total=constant+native+implicit[number]+lower[number]
        support+=len(poly)+bool(native)+bool(total);assert support+len(native_cache)+50000<150000
        record={'kernel619_index':number,'old_known_part_with_honest_second_polynomial':encode_poly(poly),
          'honest_second_linear_and_quadratic_terms_cancel':not symbols,'native_known_response':str(native),
          'selected135_implicit_contribution':str(implicit[number]),'actual_common_D_lower':str(lower[number]),
          'actual_corrected_pure_quartic_if_bound':str(total) if not symbols else None}
        records.append(record)
        if symbols:unbound.append(number)
        elif total:
            witness={**record,'actual_original768_constant_combination':row['actual_original768_constant_combination'],
              'all_source_inputs_fixed':True,'canonical_pure_class':'lasttarget(source55)*rA^4',
              'actual_HF0_through4':[5,3,2,1,0],'actual_T_total_length':11,
              'proof':'nonzeroactualrelationkillssoleknown gr4candidate; m4T=m5T, Nakayama=>m4T=0'}
            dump('actual-nonzero-quartic-witness.json',witness);break
    dump('phase3-combined619-evaluated-cases.json',{'priority_order':priority,'evaluated_cases':records,
      'unbound_evaluated_indices':unbound,'honest_second_parameter_registry':registry,
      'actual_native_known_source_values':[[list(k),str(v)] for k,v in sorted(native_cache.items())],
      'witness_found':bool(witness),'all619_evaluated':len(records)==619,
      'allzero_without_completeness_means_only_family_negative':True,
      'known_N2_zz_omitted_only_because_requested3monomials_have_atmostonez':True,
      'N3_actual_components_drop_by_all619_S0c0zero_pinned_selected135':True,
      'source3_restored_exactly_once_by_selected135':True,'new_implicit_or_source_jet_or_coefficient':False})
    dump('phase.json',{'phase':'FINAL','witness_found':bool(witness),'evaluated_cases':len(records),'unbound_cases':len(unbound)})
    note=ROOT/'research/notes'/f'{STEM}.md'
    report={'id':STEM,'all_passed':True,'actual_nonzero_corrected_quartic_witness':bool(witness),
      'actual_T_total_length':11 if witness else None,'actual_HF0_through4':[5,3,2,1,0] if witness else None,
      'actual_evaluated_cases':len(records),'unbound_cases':unbound,'all619_evaluated':len(records)==619,
      'allzero_is_family_negative_without_separate_completeness':not witness,'new_source_implicit_coefficients':0,
      'max_symbolic_C0_live_QQ':peak,'max_symbolic_typed_live_QQ':tpeak,'all_saved_QQ_bound':support+len(native_cache)+50000,
      'input_sha256':hs,'code_sha256':sha(Path(__file__)),'helper_code_sha256':sha(Path(__file__).with_name('ambient_B18_T_quartic_symbolic_tools_20261104_v63.py')),
      'note_sha256':sha(note),'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,
      'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'independent_source_value_review':'PENDING'}
    dump('report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','actual_nonzero_corrected_quartic_witness','actual_T_total_length','actual_evaluated_cases','all619_evaluated','all_saved_QQ_bound','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
