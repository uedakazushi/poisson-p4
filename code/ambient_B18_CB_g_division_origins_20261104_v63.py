#!/usr/bin/env python3
"""ONE degree-three, two-variable jet of the actual CB unit57 system.

This computes four origin values of quotients already proved regular.  It does
not infer full membership from jets, construct new points, or solve a colon.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
from fractions import Fraction as Q
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors, forward, linear_solve

ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CB_G_DIVISION_ORIGINS_20261104_v63'
MON=[(d-i,i) for d in range(4) for i in range(d+1)]
POS={a:i for i,a in enumerate(MON)}
CONV=[[(i,j) for i,a in enumerate(MON) for j,b in enumerate(MON)
       if (a[0]+b[0],a[1]+b[1])==c] for c in MON]
Z=Q(0); O=Q(1); ZERO=(Z,)*10
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(M):return [[str(x) for x in row] for row in M]

def add(a,b):
    if a is ZERO:return b
    if b is ZERO:return a
    return tuple(x+y for x,y in zip(a,b))
def neg(a):return ZERO if a is ZERO else tuple(-x for x in a)
def mul(a,b):
    if a is ZERO or b is ZERO:return ZERO
    return tuple(sum((a[i]*b[j] for i,j in terms if a[i] and b[j]),Z) for terms in CONV)
def inv(a):
    assert a[0]
    c=1/a[0]; ans=[c]+[Z]*9
    for k in range(1,10):
        ans[k]=-c*sum((a[i]*ans[j] for i,j in CONV[k] if i and a[i] and ans[j]),Z)
    return tuple(ans)

def evaluate(nodes,needed,bindings):
    val={}
    for i in needed:
        n=nodes[i]; op=n[0]
        if op=='c':v=(Q(n[1],n[2]),)+(Z,)*9 if n[1] else ZERO
        elif op=='v':v=bindings[i]
        elif op=='a':v=add(val[n[1]],val[n[2]])
        elif op=='n':v=neg(val[n[1]])
        elif op=='m':v=mul(val[n[1]],val[n[2]])
        elif op=='i':v=inv(val[n[1]])
        else:raise ValueError(op)
        val[i]=v
    return val

def main(out):
    started=time.monotonic(); stage='read inputs'
    out.mkdir(parents=True,exist_ok=False)
    def timeout(_signum,_frame):raise TimeoutError('authorized single60second cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    source=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    fullpath=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json'
    fullpeer=ROOT/'reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63_REVIEW/report.json'
    cbmu=ROOT/'reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/report.json'
    cbmua=ROOT/'reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63/report.json'
    data=read(source); nodes=data['DAG_scalar_nodes']; cb=data['families']['CB']
    assert read(fullpeer)['all_passed'] and read(cbmu)['all_passed']
    assert read(cbmua)['actual_odd_ideal_minimal_generator_count']==10
    assert read(cbmua)['exact_conductor_multiplicity4']
    full=read(fullpath); gd=full['all4_known_regular_g_divisions']
    assert len(gd)==4 and all(v['origin_value']=='not evaluated; no unit use' for v in gd)
    varids={n[1]:i for i,n in enumerate(nodes) if n[0]=='v'}
    xids=[varids[f'x{i}'] for i in data['physical_zero_base_rows9']]
    zids=[i for pair in cb['coefficient_unknown_pairs'].values() for i in pair]
    zids += [cb['p_quadratic_relation']['s'],cb['p_quadratic_relation']['t']]
    variables=zids+xids; assert len(variables)==57 and len(set(variables))==57
    eq48=cb['all48_zero_IFT_quadratic_coefficient_equations']
    gens=cb['actual_zero_odd_I0_generators']+cb['actual_positive_odd_J0_generators']
    selected=[0,1,2,3,4,5,6,8,9]
    F9=[gens[i][1] for i in selected]; g=gens[7][1]
    assert [gens[i][0] for i in selected]==[18,19,28,29,60,67,72,46,95]
    assert gens[7][0]==63
    brows=[64,68,73,74]
    assert [full['DAG_scalar_nodes'][record['node']][3] for record in gd]==[
        f'graph{row}:known_I0_membership' for row in brows]
    bnodes=[cb['complete_normalized_physical_w_coefficients'][0][r][1] for r in brows]
    roots=eq48+F9+[g]+bnodes
    needed=sorted(ancestors(nodes,roots))
    assert {i for i in needed if nodes[i][0]=='v'}<=set(variables)
    stage='fixed rational Jacobian57'
    vals,grads=forward(nodes,set(variables)); assert all(vals[i]==0 for i in roots)
    eq=eq48+F9
    J=[[Q(str(grads[i].get(v,F(0)))) for v in variables] for i in eq]
    Jinv,det=linear_solve(J,[[Q(int(i==j)) for j in range(57)] for i in range(57)])
    assert det==Q(-3456,625)
    assert all(sum((J[i][k]*Jinv[k][j] for k in range(57) if J[i][k] and Jinv[k][j]),Z)==int(i==j)
               for i in range(57) for j in range(57))
    targets=[ZERO for i in range(57)]
    rr=[Z]*10;rr[POS[(1,0)]]=O
    ss=[Z]*10;ss[POS[(0,1)]]=O
    targets[55]=tuple(rr);targets[56]=tuple(ss)
    sol={i:ZERO for i in variables}; stages=[]
    # Only the fixed rational inverse is used at all three homogeneous orders.
    for degree in range(1,4):
        stage=f'Newton homogeneous order{degree}'
        v=evaluate(nodes,needed,sol)
        cols=[k for k,a in enumerate(MON) if sum(a)==degree]
        residual=[[v[row][k]-targets[e][k] for k in cols] for e,row in enumerate(eq)]
        correction=[[-sum((Jinv[i][k]*residual[k][j] for k in range(57) if Jinv[i][k] and residual[k][j]),Z)
                     for j in range(len(cols))] for i in range(57)]
        for j,var in enumerate(variables):
            current=list(sol[var])
            for k,col in enumerate(cols):current[col]+=correction[j][k]
            sol[var]=tuple(current)
        stages.append({'order':degree,'monomials':[list(MON[k]) for k in cols],
                       'correction57_by_order_plus1':enc(correction)})
    stage='final degree3 equations and regular quotient constants'
    v=evaluate(nodes,needed,sol)
    assert all(v[row]==targets[e] for e,row in enumerate(eq))
    gj=v[g]; assert all(gj[k]==0 for k,a in enumerate(MON) if sum(a)<3)
    leading=next(k for k,a in enumerate(MON) if sum(a)==3 and gj[k])
    quotient=[]
    for row,b in zip(brows,bnodes):
        bj=v[b]; assert all(bj[k]==0 for k,a in enumerate(MON) if sum(a)<3)
        h=bj[leading]/gj[leading]
        assert bj==tuple(h*x for x in gj)
        quotient.append({'physical_graph_row':row,'old_source_odd_coefficient_node':b,
                         'known_regular_full_DAG_g_quotient_node':gd[len(quotient)]['node'],
                         'all10_bbar_coefficients':list(map(str,bj)),
                         'quotient_origin_value':str(h)})
    evidence={'coefficient_ring':'QQ[r,s]/(r,s)^4','ordered_monomials':list(map(list,MON)),
              'single_actual_original_source_DAG_hash':sha(source),'original_source_needed_node_ids':needed,
              'single57_unknown_source_node_ids':variables,'single57_unknown_source_variable_names':[nodes[i][1] for i in variables],
              'original48_equation_nodes':eq48,'selected_F9_original_nodes':F9,
              'F9_targets':['0']*7+['r','s'],
              'fixed_J57':enc(J),'fixed_J57_inverse':enc(Jinv),
              'fixed_J57_determinant':str(det),'homogeneous_Newton_corrections':stages,
              'all57_unknown_jets':[[var,list(map(str,sol[var]))] for var in variables],
              'all57_equations_verified_through_degree3':True,
              'g_CB_f63_original_node':g,'g_CB_all10_coefficients':list(map(str,gj)),
              'lex_first_nonzero_g3_monomial':list(MON[leading]),'lex_first_nonzero_g3_coefficient':str(gj[leading]),
              'all4_regular_quotient_origin_values':quotient,
              'ord_g3_and_regular_divisibility_are_independently_proved_input_not_inferred_from_jet':True,
              'no_g_division_node_evaluated_until_final_exact_constant_ratio':True}
    encoded=json.dumps(evidence,separators=(',',':'))+'\n'; assert len(encoded.encode())<32*1024*1024
    (out/'degree3-certificate.json').write_text(encoded)
    paths=[Path(__file__),source,fullpath,fullpeer,cbmu,cbmua]
    report={'id':STEM,'all_passed':True,'status':'four_exact_regular_quotient_origin_constants_determined',
            'single_pointed57_system_only':True,'fixed_J57_determinant':str(det),
            'homogeneous_orders_computed':[1,2,3],'g_CB_order_exactly3_verified':True,
            'all4_bbar_lower_orders_zero_and_degree3_proportional_to_g':True,
            'origin_values_by_physical_row':{str(q['physical_graph_row']):q['quotient_origin_value'] for q in quotient},
            'known_regular_quotient_full_membership_not_inferred_from_jets':True,
            'no_depth_colon_nil_full_membership_or_higher_order_claim':True,
            'no_GB_new_point_source_regeneration_ungauged_core_or_p_inverse':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{'degree3-certificate.json':sha(out/'degree3-certificate.json')},
            'needed_original_source_DAG_nodes':len(needed),'source_DAG_nodes':len(nodes),
            'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-started,
            'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','origin_values_by_physical_row','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    try:main(args.output)
    except (TimeoutError,MemoryError) as error:
        (args.output/'blocked.json').write_text(json.dumps({'all_passed':False,'reason':str(error),'no_automatic_order_raise_or_retry':True})+'\n')
        raise
