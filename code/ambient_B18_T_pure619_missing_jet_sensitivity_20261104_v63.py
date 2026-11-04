#!/usr/bin/env python3
"""Keep all honest quotient second/third coefficients symbolic before source3."""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_comparison_T_constant_rank_20261104_v63 import read,sha,rec,encode,addscaled
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T_PURE619_MISSING_JET_SENSITIVITY_20261104_v63'
def dump(out,name,x):(out/name).write_text(json.dumps(x,separators=(',',':'))+'\n')
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out=out.resolve();out.mkdir(parents=True,exist_ok=False)
    paths={
      'source':'reports/AMBIENT_B18_T_CONSTANT_QUADRATIC_KERNEL_R3_ADJOINT_20261104_v63/phase2-all619-pure-cubic-source-functionals.json',
      'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
      'hfirst':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
      'cache':'reports/AMBIENT_B18_RAW_COKERNEL_ARTIN20_EVALUATION_20261104_v63/actual-source105-Artin20-values.json',
      'copy':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase2-existing-cache-B11-copy-values.json',
      'atomfirst':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
      'typedAD':'reports/AMBIENT_B18_T_NEXT3_ROOT_ATOM_ADJOINT_20261104_v63/oldB-AD-fixed-node-values-and-gradients.json',
      'composites':'reports/AMBIENT_B18_T_NEXT3_EXISTING_CACHE_TRANSPORT_20261104_v63/phase3-composite-second-source-functional.json'}
    hashes={p:sha(ROOT/p) for p in paths.values()};dump(out,'phase.json',{'phase':'source_refs_before_missing_jet_AD','input_sha256':hashes})
    D={k:read(ROOT/p) for k,p in paths.items()};S,H,J,C,CP,A,TD,CV=[D[k] for k in paths];nodes=H['DAG_scalar_nodes']
    basis=C['A20_basis'];const=basis.index([0]*8);rmon=[0]*8;rmon[7]=1;rp=basis.index(rmon)
    def first(p):
        v=rec(p);return v.get(const,F(0)),v.get(rp,F(0))
    known={i:first(p) for i,p in C['all_requested_original_source_root_values']}
    known.update({i:first(p) for i,p in zip(C['joint105_unknown_order'],C['all105_unknown_final_values'])})
    known.update({i:(rec(p).get(0,F(0)),rec(p).get(1,F(0))) for i,p in CP['all57_copied_unknown_values_through2']})
    for i in H['actual_C0_basis8_scalar_nodes']:known[i]=(F(0),F(0))
    opaque={}
    for key in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
        for i,p in J[key].items():opaque[int(i)]=(F(p[0]),F(p[8])+F(p[9]))
    for p in A['all35_ordered_CA_prefix_firsts']:
        q=rec(p['first9_bar']);opaque[p['node']]=(q.get(0,F(0)),q.get(8,F(0)))
    for i,p in A['all8_g_quotient_overrides']:
        q=rec(p);opaque[i]=(q.get(0,F(0)),q.get(8,F(0)))
    roots={i for row in S['all619'] for key in ['rawP8_last_cubic_roots','Xp_second_r2_roots'] for i,c in row[key]}
    needed=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in needed:continue
        needed.add(i)
        if i in known or i in opaque:continue
        n=nodes[i]
        if n[0] in ['a','m']:todo.extend(n[1:3])
        elif n[0] in ['n','i']:todo.append(n[1])
        else:assert n[0]=='c',('MISSING_CACHED_FIRST_LEAF',i,n)
    # Existing X2 functions are separately all fixed by old762/108 caches;
    # they cannot add an unknown second/third quotient sensitivity here.
    tdata={i:(rec(g),rec(q)) for i,v,g,f,q in TD['node_data']};cdata={i:(rec(v),rec(g)) for i,v,g in CV['all108_composite_values_and_opaque_second_gradients']}
    typed={q[1] for row in S['all619'] for q,a in row['quadratic_source_root_correction'] if q[0]==0}
    direct={q[1] for row in S['all619'] for q,a in row['quadratic_source_root_correction'] if q[0]==1}
    assert typed<=tdata.keys() and direct<=cdata.keys()
    atoms={i for j in typed for i in tdata[j][0]};assert atoms<=cdata.keys() and all(not cdata[i][1] for i in atoms|direct)
    v={};deps0={};deps1={}
    for i in sorted(needed):
        n=nodes[i];op=n[0]
        if i in known:v[i]=known[i];deps0[i]=set();deps1[i]=set()
        elif i in opaque:v[i]=opaque[i];deps0[i]={i};deps1[i]=set()
        elif op=='c':v[i]=(F(n[1],n[2]),F(0));deps0[i]=set();deps1[i]=set()
        elif op in ['a','n']:
            js=n[1:3] if op=='a' else n[1:]
            v[i]=tuple(sum((v[j][k] for j in js),F(0))*(1 if op=='a' else -1) for k in range(2))
            deps0[i]=set().union(*(deps0[j] for j in js));deps1[i]=set().union(*(deps1[j] for j in js))
        elif op=='m':
            a,b=n[1:3];a0,a1=v[a];b0,b1=v[b];v[i]=(a0*b0,a0*b1+a1*b0)
            deps0[i]=(deps0[a] if b0 else set())|(deps0[b] if a0 else set())
            deps1[i]=(deps1[a] if b0 else set())|(deps0[a] if b1 else set())|(deps1[b] if a0 else set())|(deps0[b] if a1 else set())
        elif op=='i':
            a=n[1];a0,a1=v[a];assert a0;v[i]=(1/a0,-a1/a0**2);deps0[i]=deps0[a];deps1[i]=deps1[a]|(deps0[a] if a1 else set())
        else:raise AssertionError((i,n))
    bound=2*len(needed)+sum(len(deps0[i])+len(deps1[i]) for i in needed)+619*100+10000
    assert bound<150000
    dump(out,'preflight.json',{'cached_first_DAG_closure':len(needed),'structural_Jacobian0_plus_r_slots':sum(len(deps0[i])+len(deps1[i]) for i in needed),
      'TOTAL_support_bound':bound,'cap':150000,'all_X2_sources_fixed_from_preexisting_caches':True})
    g0={};g1={}
    for i in sorted(needed):
        n=nodes[i];op=n[0]
        if i in known or op=='c':a={};b={}
        elif i in opaque:a={i:F(1)};b={}
        elif op=='a':a=dict(g0[n[1]]);b=dict(g1[n[1]]);addscaled(a,g0[n[2]],F(1));addscaled(b,g1[n[2]],F(1))
        elif op=='n':a={j:-c for j,c in g0[n[1]].items()};b={j:-c for j,c in g1[n[1]].items()}
        elif op=='m':
            x,y=n[1:3];x0,x1=v[x];y0,y1=v[y];a={};b={}
            addscaled(a,g0[x],y0);addscaled(a,g0[y],x0)
            addscaled(b,g1[x],y0);addscaled(b,g0[x],y1);addscaled(b,g1[y],x0);addscaled(b,g0[y],x1)
        elif op=='i':
            x=n[1];x0,x1=v[x];a={j:-c/x0**2 for j,c in g0[x].items()};b={}
            addscaled(b,g1[x],-1/x0**2);addscaled(b,g0[x],2*x1/x0**3)
        else:raise AssertionError((i,n))
        assert set(a)<=deps0[i] and set(b)<=deps1[i];g0[i]=a;g1[i]=b
    checks=[];bad2=set();bad3=set();nonzero2=[];nonzero3=[]
    for row in S['all619']:
        s2={};s3={}
        for root,c in row['rawP8_last_cubic_roots']:
            addscaled(s2,g1[root],F(c));addscaled(s3,g0[root],F(c))
        for root,c in row['Xp_second_r2_roots']:addscaled(s2,g0[root],F(c))
        if s2:nonzero2.append(row['constant_quadratic_kernel_index619']);bad2.update(s2)
        if s3:nonzero3.append(row['constant_quadratic_kernel_index619']);bad3.update(s3)
        checks.append({'kernel619_index':row['constant_quadratic_kernel_index619'],
          'opaque_second_r2_sensitivity':encode(s2),'opaque_third_r3_sensitivity':encode(s3)})
    support=sum(len(g0[i])+len(g1[i])+2 for i in needed)+sum(len(r['opaque_second_r2_sensitivity'])+len(r['opaque_third_r3_sensitivity']) for r in checks)
    assert support<150000
    dump(out,'actual619-missing-quotient-jet-sensitivities.json',{'all619':checks,
      'source_root_Jacobians0_and_r':[[i,encode(g0[i]),encode(g1[i])] for i in sorted(roots)],
      'all_second_unknown_nodes_remaining':[[i,nodes[i]] for i in sorted(bad2)],
      'all_third_unknown_nodes_remaining':[[i,nodes[i]] for i in sorted(bad3)],
      'four_gdiff2_and_ALL_ordinary_h_g_second_and_third_coefficients_retained_before_test':True,
      'all_X2_correction_unknown_second_gradients_individually_zero':True,
      'all_missing_second_and_third_sensitivities_cancel':not bad2 and not bad3,
      'no_new_source3_or5_jet_or_actual619_values':True,'expanded_QQ_support':support})
    note=ROOT/'research/notes'/f'{STEM}.md';dump(out,'phase.json',{'phase':'FINAL','remaining_unknown_second_nodes':len(bad2),'remaining_unknown_third_nodes':len(bad3)})
    report={'id':STEM,'all_passed':True,'all619_missing_second_and_third_coefficients_cancel':not bad2 and not bad3,
      'nonzero_second_sensitivity_relation_count':len(nonzero2),'nonzero_third_sensitivity_relation_count':len(nonzero3),
      'remaining_opaque_second_node_count':len(bad2),'remaining_opaque_third_node_count':len(bad3),
      'no_new_source_jet_values_or_actual619_values_computed':True,'input_sha256':hashes,
      'code_sha256':sha(Path(__file__)),'note_sha256':sha(note),
      'evidence_sha256':{p.name:sha(p) for p in out.iterdir()},'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
    dump(out,'report.json',report);print(json.dumps({k:report[k] for k in ['id','all_passed','all619_missing_second_and_third_coefficients_cancel','remaining_opaque_second_node_count','remaining_opaque_third_node_count','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);main(p.parse_args().output)
