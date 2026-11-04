#!/usr/bin/env python3
"""Read-only cost/source preparation; no218 Jacobian test or dual77 solve."""
import argparse,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];ID='AMBIENT_B18_NIL_DERIVATION_EXTRACTOR_PREFLIGHT_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    bodyp=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    bp=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json';ep=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/report.json';builderp=ep.parent/'builder.json'
    body,bundle,old,builder=map(read,[bodyp,bp,ep,builderp]);assert body['all_passed'] and old['all_passed']
    gates=[ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63_REVIEW/report.json']
    assert all(read(p)['all_passed'] for p in gates)
    names=builder['ambient_variables'];vars=s.symbols(' '.join(names));vm=dict(zip(names,vars));p,ss,r=s.symbols('p s r');ks=s.symbols('k0:9')
    restricted={p:0,ss:0,**{ks[i]:0 for i in [3,4,6,7,8]}}
    source=[s.expand(s.sympify(x).subs(restricted)) for x in body['families']['CA']['full_61_coordinates']]
    q=[source[i] for i in bundle['q_indices_zero_based']];images={n:s.Integer(0) for n in names}
    for j,i in enumerate(builder['common_q_selected_rows']):images[f'v{j}']=q[i]
    for j,i in enumerate(builder['remaining_q_rows']):images[f'z{j}']=q[i]
    free=['v1','v4','v5','v7'];KtoV={ks[0]:vm['v4'],ks[1]:-vm['v1'],ks[2]:vm['v5']/2,ks[5]:vm['v7']/3}
    for n in names:
        at0=s.expand(images[n].subs(r,0).subs(KtoV))
        assert at0==(vm[n] if n in free else 0)
    candidate=[s.expand(s.diff(images[n],r).subs(r,0).subs(KtoV)) for n in names]
    F0=s.sympify(old['nilradical_generator'],locals=vm)
    normalizer=s.expand(sum(s.diff(F0,vm[n])*v for n,v in zip(names,candidate)))
    equations=builder['component_presentations']['CA']+builder['component_presentations']['CB'];assert len(equations)==218
    support=[i for i,v in enumerate(candidate) if v];counts=[];upper=0
    for text in equations:
        z=s.Poly(s.sympify(text,locals=vm),*vars)
        counts.append({'terms':len(z.terms()),'total_degree':z.total_degree()})
        upper+=sum(sum(mon[i] for i in support) for mon,c in z.terms())
    # This is input preparation, deliberately NOT the tangent-identity calculation.
    record={'actual_CA_source_specialization':{'p':0,'s':0,'zero_k':[3,4,6,7,8]},'kept_source_variables':['r','k0','k1','k2','k5'],
        'actual_old49_map':{n:str(images[n]) for n in names},
        'nil_support_parameter_identification':{'v1':'-k1','v4':'k0','v5':'2*k2','v7':'3*k5'},
        'candidate_derivation49':[str(v) for v in candidate],'candidate_nonzero_coordinates':{names[i]:str(candidate[i]) for i in support},
        'candidate_is_constant':all(not v.free_symbols for v in candidate),
        'dF_on_candidate':str(normalizer),'F':str(F0),'candidate_nonzero_scalar_normalization':bool(normalizer and not normalizer.free_symbols),
        'all218_relation_term_degree_counts':counts,'sum_relation_terms':sum(z['terms'] for z in counts),'maximum_relation_degree':max(z['total_degree'] for z in counts),
        'upper_bound_selected_partial_term_operations':upper,'actual218_Jacobian_dot_NOT_computed':True,
        'proposed_next_exact_test':'QQ[v1,v4,v5,v7] all218 Jac dot candidate/dF after oldann45 specialization; no point, GB or source jet',
        'proposed_positive_followthrough':'Compose old20cut nilcurve and actualN56coefficients via physicalR8, replace curveparameter by z1=b61, save77pointed unit equations and77linear dual equations with source seed candidate/dF; no expandedformalseries or77rationalGaussian.',
        'dual77_center_unit_from_accepted_blocks':'2*(-80)*298598400=-47775744000 up to declared row/column ordering',
        'budgets':{'seconds':60,'address_space_bytes':2000000000,'artifact_bytes':32*1024*1024},
        'estimated_input_scalar_nodes_before_dual':70000,'safe_dual_DAG_strategy':'two-copy AD and implicit linear77unit input; do not expand adjugate/inverse matrices',
        'no_delta_ringmap_extractor_or_extension_coefficients_proved':True}
    out.mkdir(parents=True,exist_ok=False);cert=out/'preflight.json';cert.write_text(json.dumps(record,indent=2)+'\n')
    paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),bodyp,bp,ep,builderp,*gates]
    report={'id':ID,'all_passed':True,'source_binding_candidate_and_finite_cost_prepared':True,'all218_tangent_identity_or_delta_extractor_proved':False,
        'candidate_nonzero_coordinates':record['candidate_nonzero_coordinates'],'dF_on_candidate':record['dF_on_candidate'],
        'relation_terms':record['sum_relation_terms'],'partial_term_operations_upper_bound':upper,
        'next_coefficient_job_requires_separate_parent_authorization':True,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{cert.name:sha(cert)},
        'seconds':time.monotonic()-st,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
    print(json.dumps({k:report[k] for k in ['all_passed','candidate_nonzero_coordinates','dF_on_candidate','relation_terms','partial_term_operations_upper_bound','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
