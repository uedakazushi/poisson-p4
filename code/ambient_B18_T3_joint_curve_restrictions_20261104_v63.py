#!/usr/bin/env python3
"""ONE actual CAcurve57+CBsource48 unit system and first degree<=3 restrictions."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import sympy as sy
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors,forward,linear_solve,replay
from ambient_B18_CB_g_division_origins_20261104_v63 import MON,POS,Z,ZERO,evaluate,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
    started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
    out.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/actual-source-coefficients.json'
    annp=ROOT/'reports/AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_PREFLIGHT_20261104_v63/report.json'
    cap=ROOT/'reports/AMBIENT_B18_T3_CURVE_ANNIHILATOR_20261104_v63_REVIEW/report.json'
    sourcepeer=ROOT/'reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63_REVIEW/report.json'
    data=read(source);assert read(annp)['all_passed'] and read(cap)['all_passed'] and read(sourcepeer)['all_passed']
    nodes=data['DAG_scalar_nodes'];ca,cb=[data['families'][f] for f in ['CA','CB']]
    varids={n[1]:i for i,n in enumerate(nodes) if n[0]=='v'}
    xids=[varids[f'x{i}'] for i in data['physical_zero_base_rows9']]
    def zs(f):return [i for p in f['coefficient_unknown_pairs'].values() for i in p]+[f['p_quadratic_relation']['s'],f['p_quadratic_relation']['t']]
    az,bz=zs(ca),zs(cb); variables=az+xids+bz
    assert len(variables)==105 and len(set(variables))==105
    agens=ca['actual_zero_odd_I0_generators']+ca['actual_positive_odd_J0_generators']
    sel=[0,1,2,3,4,5,6,8,9];F9=[agens[i][1] for i in sel]
    g=agens[7][1]; assert agens[7][0]==63
    brows=[r for r,i in cb['actual_zero_odd_I0_generators']]
    bnodes=[i for r,i in cb['actual_zero_odd_I0_generators']]
    eq=ca['all48_zero_IFT_quadratic_coefficient_equations']+F9+cb['all48_zero_IFT_quadratic_coefficient_equations']
    roots=eq+[g]+bnodes;needed=sorted(ancestors(nodes,roots))
    assert {i for i in needed if nodes[i][0]=='v'}<=set(variables)
    # A literal, finite105-unit algebraic input, independent of the jet computation.
    D=DAG();original=replay(D,nodes,{})
    assert all(original[i]==i for i in range(len(nodes)))
    R,S=D.var('JOINT_CURVE_r'),D.var('JOINT_CURVE_s')
    binding={i:D.var('JOINT_INV_'+nodes[i][1]) for i in variables}
    mapped=replay(D,nodes,binding,set(needed))
    algebraic_eq=[mapped[i] for i in eq]
    algebraic_eq[55]=D.add(algebraic_eq[55],D.neg(R));algebraic_eq[56]=D.add(algebraic_eq[56],D.neg(S))
    assert all(D.values[i]==0 for i in algebraic_eq)
    vals,grads=forward(nodes,set(variables));assert all(vals[i]==0 for i in roots)
    J=[[grads[i].get(v,F(0)) for v in variables] for i in eq]
    Jinv,det=linear_solve(J,[[F(int(i==j)) for j in range(105)] for i in range(105)])
    _,detCA=linear_solve([r[:57] for r in J[:57]],[[F(int(i==j)) for j in range(57)] for i in range(57)])
    _,detCB=linear_solve([r[57:] for r in J[57:]],[[F(int(i==j)) for j in range(48)] for i in range(48)])
    assert det==detCA*detCB and detCA==F(-3623878656,625) and detCB==F(-331776,625)
    assert all(J[i][j]==0 for i in range(57) for j in range(57,105))
    # Block inverse identity is also saved/rechecked, not merely determinant output.
    assert all(sum((J[i][k]*Jinv[k][j] for k in range(105) if J[i][k] and Jinv[k][j]),Z)==int(i==j)
               for i in range(105) for j in range(105))
    targets=[ZERO for i in range(105)];rv=[Z]*10;rv[POS[(1,0)]]=F(1);sv=[Z]*10;sv[POS[(0,1)]]=F(1)
    targets[55]=tuple(rv);targets[56]=tuple(sv)
    sol={i:ZERO for i in variables};stages=[]
    for degree in range(1,4):
        v=evaluate(nodes,needed,sol);cols=[i for i,a in enumerate(MON) if sum(a)==degree]
        res=[[v[row][j]-targets[e][j] for j in cols] for e,row in enumerate(eq)]
        cor=[[-sum((Jinv[i][k]*res[k][j] for k in range(105) if Jinv[i][k] and res[k][j]),Z)
              for j in range(len(cols))] for i in range(105)]
        for i,var in enumerate(variables):
            a=list(sol[var])
            for j,col in enumerate(cols):a[col]+=cor[i][j]
            sol[var]=tuple(a)
        stages.append({'degree':degree,'monomials':[list(MON[i]) for i in cols],'all105_corrections':enc(cor)})
    v=evaluate(nodes,needed,sol);assert all(v[row]==targets[e] for e,row in enumerate(eq))
    gj=v[g];assert all(gj[i]==0 for i,a in enumerate(MON) if sum(a)<3) and any(gj[i] for i,a in enumerate(MON) if sum(a)==3)
    restrictions=[{'CB_physical_zero_row':row,'old_source_node':node,'all_degree_le3_coefficients':list(map(str,v[node]))}
                  for row,node in zip(brows,bnodes)]
    # Only the leading homogeneous forms give this honest initial-ideal gate.
    r,s,t=sy.symbols('r s t')
    def poly(jet,degree):return sum(sy.Rational(a.numerator,a.denominator)*r**m[0]*s**m[1] for a,m in zip(jet,MON) if sum(m)==degree)
    g3=poly(gj,3);forms=[g3];initial=[]
    for row,node in zip(brows,bnodes):
        jet=v[node];degree=next((d for d in range(4) if any(a for a,m in zip(jet,MON) if sum(m)==d)),None)
        p=poly(jet,degree) if degree is not None else sy.Integer(0)
        if degree is not None:forms.append(p)
        initial.append({'CB_row':row,'first_detected_order':degree,'leading_form':str(p)})
    G=sy.Poly(forms[0].subs({r:1,s:t}),t,domain=sy.QQ)
    for p in forms[1:]:G=sy.gcd(G,sy.Poly(p.subs({r:1,s:t}),t,domain=sy.QQ))
    common_infinity=all(p.subs({r:0,s:1})==0 for p in forms)
    mprimary=G.degree()==0 and not common_infinity
    literal={'DAG_scalar_nodes':D.nodes,'formal_curve_parameter_nodes':[R,S],
             'all105_unknown_bindings':[[a,b] for a,b in binding.items()],
             'all105_pointed_unit_equations':algebraic_eq,'unit_Jacobian_determinant':str(det),
             'actual_CA_curve_g_node':mapped[g],
             'actual_CB_zero8_restriction_nodes':[[row,mapped[i]] for row,i in zip(brows,bnodes)],
             'actual_common_physical_base_rows':data['physical_zero_base_rows9'],
             'annihilator_of_T3':'(g_CA,all8 actual CB zero restrictions) in QQ[[r,s]]',
             'no_p_inverse_arbitrary_coefficient_or_truncated_membership_definition':True}
    jets={'coefficient_ring':'QQ[r,s]/(r,s)^4','monomials':list(map(list,MON)),
          'all105_unknown_original_node_order':variables,'all105_equation_original_node_order':eq,
          'CA_F9_targets':['0']*7+['r','s'],'fixed105_Jacobian':enc(J),'fixed105_Jacobian_inverse':enc(Jinv),
          'det_CA57':str(detCA),'det_CB48':str(detCB),'det_joint105':str(det),
          'homogeneous_corrections':stages,'all105_unknown_jets':[[i,list(map(str,sol[i]))] for i in variables],
          'all105_equations_zero_through_degree3':True,'g_CA_all10_coefficients':list(map(str,gj)),
          'g_CA_degree3':str(g3),'all8_CB_zero_restrictions':restrictions,
          'leading_forms':initial,'all_leading_forms_common_finite_slope_polynomial':str(G.as_expr()),
          'all_leading_forms_common_infinity':bool(common_infinity),'initial_annihilator_m_primary_certified':bool(mprimary)}
    for name,obj in [('joint105-algebraic-DAG.json',literal),('joint-degree3-certificate.json',jets)]:
        encoded=json.dumps(obj,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
        (out/name).write_text(encoded)
    paths=[Path(__file__),source,annp,cap,sourcepeer,ROOT/'code/ambient_B18_CB_g_division_origins_20261104_v63.py']
    report={'id':STEM,'all_passed':True,'status':'actual_joint105_and_degree3_ideal_restrictions',
            'literal_joint105_unit_algebraic_system_constructed':True,'joint105_determinant':str(det),
            'all105_equations_verified_through_degree3':True,'g_CA_degree3':str(g3),
            'CB_zero_restriction_initial_forms':initial,'initial_common_finite_slope_factor':str(G.as_expr()),
            'initial_common_infinity':bool(common_infinity),'m_primary_annihilator_proved_by_initial_forms':bool(mprimary),
            'no_nonfiniteness_inferred_if_initial_gate_fails':True,
            'no_T3_nonzero_colon_depth_nil_or_module_fiber_inference':True,
            'no_GB_new_point_ungauged_core_p_inverse_or_order_above3':True,
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'evidence_sha256':{name:sha(out/name) for name in ['joint105-algebraic-DAG.json','joint-degree3-certificate.json']},
            'original_DAG_ancestor_count':len(needed),'literal_DAG_node_count':len(D.nodes),
            'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['all_passed','g_CA_degree3','CB_zero_restriction_initial_forms','initial_common_finite_slope_factor','initial_common_infinity','m_primary_annihilator_proved_by_initial_forms','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
