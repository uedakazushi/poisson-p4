#!/usr/bin/env python3
"""Read-only D0 firstjet nuisance/support preflight; no source4 job."""
import argparse
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import resource,signal,time
from ambient_B18_T3_zero_structural_20261104_v63 import insert_basis
from ambient_B18_top_graph_telescoping_review_20261104_v63 import closure
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T0_FIRSTJET_SENSITIVITY_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main(out):
    st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    cp=ROOT/'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json'
    sp=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fp=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    jointp=ROOT/'reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/joint-degree3-certificate.json'
    c,s,f,j=map(read,[cp,sp,fp,jp]);nodes=c['DAG_scalar_nodes'];old=s['DAG_scalar_nodes'];N=len(old)
    adapter=[list(map(F,a)) for a in j['YB_to_YA_linear']];binding={}
    for name in ['CA','CB']:
        for i,a in j['both57_source2'][name]['all57_source2jets']:
            d=list(map(F,a[1:10]))
            if name=='CB':d=[sum((d[k]*adapter[k][q] for k in range(9)),F(0)) for q in range(9)]
            v=[F(0)]+d+[F(0)]*8
            if i in binding:assert binding[i]==v
            binding[i]=v
    def product(a,b):return [a[0]*b[0]]+[a[0]*b[q]+b[0]*a[q] for q in range(1,18)]
    def inverse(a):assert a[0];return [1/a[0]]+[-x/a[0]**2 for x in a[1:]]
    original=[]
    for i,n in enumerate(old):
        o=n[0]
        if o=='c':v=[F(n[1],n[2])]+[F(0)]*17
        elif o=='v':v=binding[i]
        elif o=='a':v=[x+y for x,y in zip(original[n[1]],original[n[2]])]
        elif o=='n':v=[-x for x in original[n[1]]]
        elif o=='m':v=product(original[n[1]],original[n[2]])
        elif o=='i':v=inverse(original[n[1]])
        else:raise ValueError(o)
        original.append(v)
    known={int(i):list(map(F,a))+[F(0)]*8 for i,a in j['all784_CB_h_firstjets_in_YA'].items()}
    known.update({int(i):list(map(F,a))+[F(0)]*8 for i,a in j['all792_CA_h_firstjets'].items()})
    sourceg=[q['node'] for q in f['all4_known_regular_g_divisions']]+[q['node'] for q in c['all4_known_regular_gA_divisions']]
    for k,i in enumerate(sourceg):v=[F(0)]*18;v[10+k]=1;known[i]=v
    known[f['CB_g_restriction_node']]=[F(0)]*18
    known[c['actual_C0_basis8_scalar_nodes'][-1]]=[F(0)]*18
    mon=list(map(tuple,j['source2_monomials']));pos={m:i for i,m in enumerate(mon)}
    CA=j['both57_source2']['CA']['pure_function2jets']
    for rec in c['all35_regular_coordinate_divisions']:
        p=list(map(F,CA[str(rec['original_pure_CA_source_function'])]))
        for k,q in enumerate(rec['all7_regular_coordinate_divisions']):
            a=[p[1+k]]+[p[pos[tuple(sorted((k,t)))]] if t<=k or t>=7 else F(0) for t in range(9)]
            assert a[0]==F(q['origin']);known[q['node']]=a+[F(0)]*8
    A=c['honest_D0_relation_columns'];needed=set();stack=[i for col in A for i in col]
    while stack:
        i=stack.pop()
        if i in needed:continue
        needed.add(i)
        if i<N or i in known:continue
        n=nodes[i]
        if n[0] in ['a','m']:stack+=n[1:]
        elif n[0] in ['n','i']:stack.append(n[1])
        elif n[0]!='c':raise ValueError((i,n))
    V={}
    for i in sorted(needed):
        if i<N:v=original[i]
        elif i in known:v=known[i]
        else:
            n=nodes[i];o=n[0]
            if o=='c':v=[F(n[1],n[2])]+[F(0)]*17
            elif o=='a':v=[x+y for x,y in zip(V[n[1]],V[n[2]])]
            elif o=='n':v=[-x for x in V[n[1]]]
            elif o=='m':v=product(V[n[1]],V[n[2]])
            elif o=='i':v=inverse(V[n[1]])
            else:raise ValueError(o)
        V[i]=v
    fixed=[];knowncols=[];records=[];effects=[]
    for index,col in enumerate(A):
        assert all(V[i][0]==0 for i in col)
        a=[V[i][1+k] for i in col for k in range(9)];knowncols.append(a)
        sens=[[V[i][10+q] for i in col] for q in range(8)]
        if not any(any(a) for a in sens):insert_basis(fixed,a)
        for q,coeff in enumerate(sens):
            for t in range(2):
                direction=adapter[7+t] if q<4 else [F(int(k==7+t)) for k in range(9)]
                v=[coef*direction[k] for coef in coeff for k in range(9)]
                if any(v):effects.append([index,q,t,v])
        if any(any(a) for a in sens):records.append({'relation':index,'sensitivities8_by9':[[str(x) for x in a] for a in sens]})
    enlarged=[(p,dict(b)) for p,b in fixed]
    for a in knowncols:insert_basis(enlarged,a)
    for _,_,_,a in effects:insert_basis(enlarged,a)
    outside=[]
    for index,q,t,a in effects:
        b=[(p,dict(v)) for p,v in fixed]
        if insert_basis(b,a) is not None:outside.append([index,q,t])
    # Source4 cost/size only: reuse the CA57 inverse and saved joint105 source3.
    ca=s['families']['CA'];Y=j['both57_source2']['CA']['F9_nodes'];f63=j['both57_source2']['CA']['f63']
    zero_rows=[64,68,73,74];broots=[ca['complete_normalized_physical_w_coefficients'][0][r][1] for r in zero_rows]
    source_roots=ca['all48_zero_IFT_quadratic_coefficient_equations']+Y+[f63]+broots
    source_needed=closure(old,source_roots)
    allowed=set(j['both57_source2']['CA']['unknown57'])
    assert {i for i in source_needed if old[i][0]=='v'}<=allowed
    joint=read(jointp);assert allowed<=set(i for i,_ in joint['all105_unknown_jets'])
    g3=list(map(F,c['cached_actual_CA_gA_all_degree3_coefficients']));assert g3[6]*g3[9]==F(1,4)
    ev={'all640_known_firstjet81_columns':[[str(x) for x in a] for a in knowncols],
      'all_nonzero_opaque_function_sensitivities':records,
      'unknown_function_order':['oldCBg0','oldCBg1','oldCBg2','oldCBg3','newCAH0','newCAH1','newCAH2','newCAH3'],
      'each_function_has_two_unknown_first_derivatives':True,
      'fixed_relation_firstspan_rank':len(fixed),'max_enlarged_span_rank':len(enlarged),
      'unknown_effects_not_in_fixed_span':outside,
      'CA57_source4_ancestor_count':len(source_needed),'CA57_source4_required_unknowns':57,
      'CA57_source4_new_coefficients':57*5,'CA57_source4_bivariate_monomials_through4':15,
      'new_quotient_firsts_from_Bbar4_equals_gA3_H1':True,
      'quartic_gate_2_by2_endpoint_minor':'1/4',
      'no_source4_coefficients_computed':True,'no_socle_or_depth_conclusion':True}
    ep=out/'formal-sensitivities-and-source4-size.json';ep.write_text(json.dumps(ev,separators=(',',':'))+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',cp,sp,fp,jp,jointp]
    r={'id':STEM,'all_passed':True,'status':'firstjet_nuisance_and_exact_source4_size_preflight_only',
      'relations_with_nonzero_unknown_function_first_sensitivity':len(records),
      'new_CA_H_first_sensitivity_cancels_from_all_entries':not any(any(F(x) for x in a['sensitivities8_by9'][q]) for a in records for q in range(4,8)),
      'fixed_relation_firstspan_rank':len(fixed),'max_enlarged_span_rank':len(enlarged),
      'new_CA_H_effects_outside_fixed_span':sum(q>=4 for _,q,_ in outside),
      'old_CB_g_effects_outside_fixed_span':sum(q<4 for _,q,_ in outside),
      'source4_job_launched':False,'T0_socle_or_depth_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps({k:r[k] for k in ['all_passed','relations_with_nonzero_unknown_function_first_sensitivity','fixed_relation_firstspan_rank','max_enlarged_span_rank','new_CA_H_effects_outside_fixed_span','old_CB_g_effects_outside_fixed_span','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
