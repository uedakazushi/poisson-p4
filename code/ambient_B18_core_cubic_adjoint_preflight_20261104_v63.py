#!/usr/bin/env python3
"""Preflight only: source2 reuse, reverse-linear AD, adjoint98 cancellation.

No cubic/quartic source evaluation or source3/4 solve is executed.
"""
import argparse,hashlib,itertools,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors,linear_solve
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_CORE_CUBIC_ADJOINT_PREFLIGHT_20261104_v63'
M9=[()]+[(i,) for i in range(9)]+list(itertools.combinations_with_replacement(range(9),2));P9={m:i for i,m in enumerate(M9)}
M4=[()]+[(i,) for i in range(4)]+list(itertools.combinations_with_replacement(range(4),2));P4={m:i for i,m in enumerate(M4)}
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    return x
def mul9(a,b):
    r=[a[0]*b[0]]+[a[0]*b[1+i]+b[0]*a[1+i] for i in range(9)]
    for i,j in M9[10:]:
        c=a[0]*b[P9[i,j]]+b[0]*a[P9[i,j]]+a[1+i]*b[1+j]
        if i!=j:c+=a[1+j]*b[1+i]
        r.append(c)
    return tuple(r)
def compose9(a,sub):
    r=[F(a[0])]+[F(0)]*54
    for m,c in zip(M9[1:],a[1:]):
        c=F(c)
        if not c:continue
        v=sub[m[0]] if len(m)==1 else mul9(sub[m[0]],sub[m[1]])
        for j,x in enumerate(v):r[j]+=c*x
    return tuple(r)
Z1=(F(0),)*5
def add1(a,b):return tuple(x+y for x,y in zip(a,b))
def neg1(a):return tuple(-x for x in a)
def mul1(a,b):return (a[0]*b[0],)+tuple(a[0]*b[j]+b[0]*a[j] for j in range(1,5))
def eval_linear(nodes,need,bind):
    V={}
    for i in need:
        n=nodes[i];o=n[0]
        if o=='c':v=(F(n[1],n[2]),)+Z1[1:]
        elif o=='v':v=bind[i][:5]
        elif o=='a':v=add1(V[n[1]],V[n[2]])
        elif o=='n':v=neg1(V[n[1]])
        elif o=='m':v=mul1(V[n[1]],V[n[2]])
        elif o=='i':
            x=V[n[1]];assert x[0];v=(1/x[0],)+tuple(-a/x[0]**2 for a in x[1:])
        else:raise ValueError(o)
        V[i]=v
    return V
def reverse_linear(nodes,need,values,root,variables):
    adj={root:(F(1),)+Z1[1:]}
    def add(i,v):
        if any(v):adj[i]=add1(adj.get(i,Z1),v)
    for i in reversed(need):
        a=adj.get(i)
        if a is None or not any(a):continue
        n=nodes[i];o=n[0]
        if o=='a':add(n[1],a);add(n[2],a)
        elif o=='n':add(n[1],neg1(a))
        elif o=='m':add(n[1],mul1(a,values[n[2]]));add(n[2],mul1(a,values[n[1]]))
        elif o=='i':add(n[1],neg1(mul1(a,mul1(values[i],values[i]))))
        elif o not in ('v','c'):raise ValueError(o)
    return [adj.get(v,Z1) for v in variables]
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
    cp=ROOT/'reports/AMBIENT_B18_NORMALIZER_CORE4_REDUCTION_20261104_v63/typed-core4-input.json'
    ep=ROOT/'reports/AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63/cached-source2-and-exact-E4-quadrics.json'
    jp=ROOT/'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json'
    pp=ROOT/'reports/AMBIENT_B18_CORE4_QUADRATIC_PREFLIGHT_20261104_v63/literal-adjoint98-core4-preflight.json'
    qp=ROOT/'reports/AMBIENT_B18_CORE4_QUADRATIC_INITIAL_20261104_v63/seven-actual-core-quadrics.json'
    C,E,J,P,Q=map(read,[cp,ep,jp,pp,qp]);nodes=C['DAG_scalar_nodes']
    unknown=C['unit_schur_unknown_nodes98'];eq=C['unit_schur_equations98'];A=[[F(c) for c in r] for r in P['actual98_origin_Jacobian']]
    assert len(unknown)==len(eq)==98
    core=[2,4,7,8];udir=[0,1,3,5,6];bind9={}
    rho=[tuple(map(F,a)) for a in E['all11_variables_as_remaining9_through2'][:9]]
    for name in ['CA','CB']:
        f9=[compose9(a,rho) for a in E['families'][name]['inverse_F9_in_common_R0_through2']]
        for node,jt in J['both57_source2'][name]['all57_source2jets']:
            q=compose9(jt,f9)
            if node in bind9:assert bind9[node]==q
            bind9[node]=q
    for j,node in enumerate(C['actual9_parameter_nodes']):
        v=(F(0),)+tuple(F(i==j) for i in range(9))+(F(0),)*45
        if node in bind9:assert bind9[node]==v
        bind9[node]=v
    four={};dual={}
    for node,jt in bind9.items():
        base=[F(0)]*15;db=[[F(0)]*15 for j in range(5)]
        for m,c in zip(M9,jt):
            if not c:continue
            us=[i for i in m if i in udir]
            if not us:base[P4[tuple(core.index(i) for i in m)]]+=c
            elif len(us)==1 and len(m)<=2:
                rest=list(m);rest.remove(us[0]);db[udir.index(us[0])][P4[tuple(core.index(i) for i in rest)]]+=c
        four[node]=tuple(base);dual[node]=db
    assert all(n in four for n in unknown)
    selected=[0,1,4,2,3];remaining=[7,8,9,10,11]
    Fselected=[C['all12_actual_difference_nodes'][i] for i in selected]
    Ffive=[C['all12_actual_difference_nodes'][i] for i in remaining]
    roots=eq+Fselected+Ffive;need=sorted(ancestors(nodes,roots))
    assert all(i in four for i in need if nodes[i][0]=='v')
    vals=eval_linear(nodes,need,four)
    gsel=[reverse_linear(nodes,need,vals,r,unknown) for r in Fselected]
    gf=[reverse_linear(nodes,need,vals,r,unknown) for r in Ffive]
    assert all(not any(g[0] for g in row) for row in gf),('remaining raw center gradient nonzero')
    assert all(not any(F(c) for c in a) for a in Q['all7_core_quadrics_padded_through3'][2:])
    assert all(not any(F(c) for c in a) for a in Q['all7_after_exact_adjoint98_correction'][2:])
    # Selected five use lambda0; five remaining use lambda1(C), with twenty
    # rational right-hand sides. All are exact constant-matrix solves.
    rhs=[[g[0] for g in row] for row in gsel]
    rhs += [[g[j+1] for g in row] for row in gf for j in range(4)]
    solT,det=linear_solve([list(c) for c in zip(*A)],[list(c) for c in zip(*rhs)])
    sol=[list(c) for c in zip(*solT)];lambda0=sol[:5];lambda1=[sol[5+4*i:9+4*i] for i in range(5)]
    assert all(sum(sol[i][k]*A[k][j] for k in range(98))==rhs[i][j] for i in range(25) for j in range(98))
    M70=[()]+[(i,) for i in range(4)]+list(itertools.combinations_with_replacement(range(4),2))+list(itertools.combinations_with_replacement(range(4),3))+list(itertools.combinations_with_replacement(range(4),4))
    conv70=sum(len(a)+len(b)<=4 for a in M70 for b in M70);conv15=sum(len(a)+len(b)<=2 for a in M4 for b in M4)
    assert len(M70)==70 and conv70==495 and conv15==45
    nmult=sum(nodes[i][0]=='m' for i in need);ninv=sum(nodes[i][0]=='i' for i in need)
    maxpairs=nmult*(conv70+10*conv15)+ninv*(4*conv70+15*conv15)
    proof={'scalar_DAG_source':str(cp.relative_to(ROOT)),'scalar_ancestor_nodes':need,
      'unit98_unknowns':unknown,'unit98_equations':eq,'unit98_Jacobian':A,'unit98_determinant':det,
      'full9_source2_monomials':M9,'all107_full9_source2_bindings':list(bind9.items()),
      'core4_order':['x112','x117','tau_CA','tau_CB'],'U5_order':['x108','x109','x113','x118','x119'],
      'all_four_variable_source2_bindings':list(four.items()),'all_five_dual_source2_bindings':list(dual.items()),
      'selected5_F_roots':Fselected,'remaining5_F_roots':Ffive,'remaining5_physical_rows':[63,64,68,73,74],
      'own_reverse_linear_selected_gradients':gsel,'own_reverse_linear_remaining_gradients':gf,
      'all5_remaining_center_unknown_gradients_zero':True,
      'selected5_lambda0_rows':lambda0,'remaining5_lambda1_rows_by_core_direction':lambda1,
      'all2450_adjoint_coefficient_identities_verified':True,
      'selected5_Diagonal_G1_coefficients':['1/12','1/12','1/24','1/12','1/12'],
      'all5_remaining_core_G2_zero_from_actual_cached_quadrics':True,
      'source3_not_solved_or_assumed_zero':'u3=-J0^-1 E3; its quartic contribution is canceled by lambda1*E3',
      'source4_drops':'dFi/du(0)=0 for all5; no u4 term',
      'mixed_dual_unknown_source3_drops':'the unprovided dual u3 coefficients multiply the same zero center gradient',
      'selectedU2_formula':'U2_j=-D5_jj^-1 * (selectedFj3(Y2)-lambda0_j*E98_3(Y2))/L1',
      'remainingK3_formula':'(Fi4(Y2)-lambda1_i(C)*E98_3(Y2)+sum_j(dUj Fi3)|core*U2_j)/L1',
      'L2_Npair3_not_required':'Gi1=0 and Gi2|core=0',
      'proposed_base70_plus_dual75_coefficient_count':145,'base70_convolution_pairs':conv70,'dual15_convolution_pairs':conv15,
      'multiplication_nodes':nmult,'unit_inverse_nodes':ninv,'conservative_scalar_pair_budget':maxpairs,
      'coefficient_slots_if_dense':len(need)*145,
      'finite_completion_targets':['five actual homogeneous core cubics','all70 coefficient L1 division checks',
        'binarydegree12 pullbacks through frozen nodalquartic normalization: a nonzero pullback proves that actual core residual is nonzero in the integral two-function surface, hence core dimension<=1'],
      'actual_cubic_or_quartic_coefficient_evaluation_run':False,'no_auto_raise':True}
    p=out/'literal-lambda1-dual5-preflight.json';p.write_text(json.dumps(enc(proof),separators=(',',':'))+'\n');assert p.stat().st_size<32*1024*1024
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md',cp,ep,jp,pp,qp,
      ROOT/'reports/AMBIENT_B18_CORE4_QUADRATIC_20261104_v63_REVIEW/report.json']
    report={'id':STEM,'all_passed':True,'status':'literal_source2_reverse_linear_lambda1_and_dual5_cancellation_preflight_only',
      'all5_actual_raw_center_unknown_gradients_zero':True,'source3_and_source4_solves_needed':False,
      'all2450_adjoint_identities_checked':True,'coefficient_slots_per_node_future':145,
      'source_ancestor_nodes':len(need),'conservative_scalar_pair_budget':maxpairs,
      'proposed_five_actual_core_cubics_computed':False,'no_higher_coefficient_job_performed':True,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{p.name:sha(p)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
