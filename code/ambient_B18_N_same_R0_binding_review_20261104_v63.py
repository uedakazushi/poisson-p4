#!/usr/bin/env python3
"""Independent actual N/R0 source, restricted QQ2jets and literal Koszul cone."""
import argparse
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import resource
import signal
import time
import sympy as s
from ambient_B18_actual_top_coefficient_preflight_review_20261104_v63 import DAGReview,QuadraticReplay,det
from ambient_W1_S61_intersection_recognition_review_20261103_v63 import rref,inverse

ROOT=Path(__file__).resolve().parents[1]
BASE='AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63';CONE='AMBIENT_B18_N_LITERAL_MODULE_CONE_20261104_v63'
def read(p):
    import gzip
    return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):
    import hashlib
    return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(A):return [list(map(F,row)) for row in A]
def mm(A,B):return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def ancestors(nodes,roots):
    found=set();todo=list(roots)
    while todo:
        i=todo.pop()
        if i in found:continue
        found.add(i);n=nodes[i]
        if n[0] not in ['c','v']:todo.extend(n[1:])
    return sorted(found)

# Own sparse quotient k[X0..6,kappa]/(XiXj, all cubics), only17 monomials.
E0=(0,)*8
LINEAR=[tuple(int(i==j) for i in range(8)) for j in range(8)]
MIXED=[tuple(int(i==j)+int(i==7) for i in range(8)) for j in range(7)]+[tuple(2*int(i==7) for i in range(8))]
ALLOWED={E0,*LINEAR,*MIXED}
def jadd(a,b):
    r=a.copy()
    for m,c in b.items():r[m]=r.get(m,F(0))+c
    return {m:c for m,c in r.items() if c}
def jscale(a,c):return {m:v*c for m,v in a.items() if v*c}
def jmul(a,b):
    r={}
    for m,c in a.items():
        for n,v in b.items():
            e=tuple(x+y for x,y in zip(m,n))
            if e in ALLOWED:r[e]=r.get(e,F(0))+c*v
    return {m:c for m,c in r.items() if c}
def jinverse(a):
    c=a.get(E0,F(0));assert c
    h=jadd(a,{E0:-c})
    return jadd({E0:1/c},jadd(jscale(h,-1/c**2),jscale(jmul(h,h),1/c**3)))

# Independent universal polynomial ring on z1..z8,nu,v; monomial is sorted indices.
def padd(a,b):
    r=a.copy()
    for m,c in b.items():r[m]=r.get(m,0)+c
    return {m:c for m,c in r.items() if c}
def pscale(a,c):return {m:v*c for m,v in a.items() if v*c}
def pmul(a,b):
    r={}
    for m,c in a.items():
        for n,v in b.items():e=tuple(sorted(m+n));r[e]=r.get(e,0)+c*v
    return {m:c for m,c in r.items() if c}
def pvar(i):return {(i,):1}
def coladd(col,i,p):col[i]=padd(col.get(i,{}),p)
def compose(cols,col):
    ans={}
    for j,p in col.items():
        for i,v in cols[j].items():coladd(ans,i,pmul(p,v))
    return {i:p for i,p in ans.items() if p}
def decode_col(col):return {int(i):{tuple(m):int(c) for c,m in ts if c} for i,ts in col}

def main(output):
    started=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
    def timeout(*_):raise TimeoutError('independent N binding/cone60s cap')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
    base=ROOT/'reports'/BASE;cb=ROOT/'reports'/CONE
    c=read(base/'actual-N-R0-binding-DAG.json');nc=read(cb/'actual-N-literal-cone-DAG.json')
    gates=[ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW/report.json',ROOT/'reports/AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63_REVIEW.json',ROOT/'reports/AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63.json']
    checked={}
    for rp in [base/'report.json',cb/'report.json',*gates]:
        r=read(rp);assert r['all_passed']
        for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h,n;checked[n]=h
        for n,h in r.get('evidence_sha256',{}).items():p=rp.parent/n;assert sha(p)==h,n;checked[str(p.relative_to(ROOT))]=h
    assert nc['DAG_scalar_nodes'][:len(c['DAG_scalar_nodes'])]==c['DAG_scalar_nodes']
    D=DAGReview(nc['DAG_scalar_nodes']);unknown=c['normalization_source28_variable_nodes'];assert len(unknown)==28
    tail=read(ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63/literal-tail-binding-DAG.json')
    body=read(ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json')
    bundle=read(ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json')
    builder=read(ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json')
    actual=c['oldN_source11_parameter_values'];names=['p','s']+[f'k{i}' for i in range(9)]
    assert [D.variables['N_SOURCE_'+n] for n in names]==unknown[:11]
    for n in names:assert actual[n]==D.add(D.variables['N_SOURCE_'+n],D.const(1 if n=='k0' else F(1,3) if n=='k5' else 0))
    def poly(text):
        z=s.sympify(text);vs=sorted(z.free_symbols,key=str)
        if not vs:return D.const(z)
        terms=[]
        for mon,co in s.Poly(z,*vs).terms():
            v=D.const(co)
            for n,k in zip(vs,mon):
                for _ in range(k):v=D.mul(v,actual[str(n)])
            terms.append(v)
        return D.total(terms)
    u=[poly(z) for z in body['families']['N']['full_61_coordinates']];q45=[u[i] for i in bundle['q_indices_zero_based']]
    pnode=actual['p'];snode=actual['s'];old49={'ss':snode,'aa':D.mul(pnode,pnode),'bb':D.mul(D.mul(pnode,pnode),pnode),'cc':D.mul(pnode,snode)}
    for j,i in enumerate(builder['common_q_selected_rows']):old49[f'v{j}']=q45[i]
    for j,i in enumerate(builder['remaining_q_rows']):old49[f'z{j}']=q45[i]
    binding={tail['old49_centered_variable_nodes'][n]:D.add(v,D.const(-1 if n in ['v4','v7'] else 0)) for n,v in old49.items()}
    binding.update(zip(tail['actual_group17_parameter_nodes'],unknown[11:]))
    cache={}
    def ts(i):
        if i in cache:return cache[i]
        n=tail['DAG_scalar_nodes'][i]
        if n[0]=='v':v=binding[i]
        elif n[0]=='c':v=D.const(F(n[1],n[2]))
        elif n[0]=='n':v=D.neg(ts(n[1]))
        elif n[0]=='i':v=D.inverse(ts(n[1]))
        else:v={'a':D.add,'m':D.mul}[n[0]](ts(n[1]),ts(n[2]))
        cache[i]=v;return v
    cuts=[ts(i) for i in tail['all20_common_cut_equation_nodes']];rho=[ts(i) for i in tail['same_physical_R0_readback_nodes9']];physical=[ts(i) for i in tail['actual_full_physical150_nodes']]
    assert cuts==c['all20_actual_physical_cut_equations'] and rho==c['actual_same_R0_readback_nodes9'] and physical==c['physical150_N_source_nodes']
    assert all(D.values[i]==0 for i in cuts+rho)
    C20=[[D.grad[i].get(v,F(0)) for v in unknown] for i in cuts];assert C20==mat(c['cut20_Jacobian'])
    piv=c['cut20_pivot_columns'];free=c['cut20_free_columns'];Ci,cdet=inverse([[row[i] for i in piv] for row in C20]);assert cdet==-96
    Z=[[-x for x in row] for row in mm(Ci,[[r[i] for i in free] for r in C20])]
    def on8(root):return [D.grad[root].get(unknown[i],F(0))+sum(D.grad[root].get(unknown[p],F(0))*Z[k][j] for k,p in enumerate(piv)) for j,i in enumerate(free)]
    rfirst=[on8(i) for i in rho];assert rfirst==mat(c['all9_R0_firstjets_on_source8']) and len(rref(rfirst)[1])==7
    chosen=c['chosen7_R0_indices'];even=c['Cohen_even8_R0_index'];ninth=c['ninth_R0_index']
    assert chosen==[1,2,3,4,5,6,7] and even==0 and ninth==8
    frame=cuts+[rho[i] for i in chosen]+[pnode];assert frame==c['frame28_equation_nodes']
    J28=[[D.grad[i].get(v,F(0)) for v in unknown] for i in frame];Ji,jdet=inverse(J28)
    assert J28==mat(c['frame28_Jacobian']) and jdet==23040
    crit=unknown[11:]+[D.variables['N_SOURCE_'+n] for n in ['k0','k2','k5']]
    assert crit==c['nonnormal_quotient_free20_source_parameters'] and det([[D.grad[i].get(v,F(0)) for v in crit] for i in cuts])==-480
    # Verify the stored restricted jets, not a fresh Newton or higher jet computation.
    jm=[tuple(sum(int(k==j) for k in mon) for j in range(8)) for mon in c['mixed_jet_monomials']]
    assert set(jm)==ALLOWED and len(jm)==17
    jetbinding={i:{m:F(v) for m,v in zip(jm,vals) if F(v)} for i,vals in c['all28_mixed_source2jets']}
    assert set(jetbinding)==set(unknown)
    for i,var in enumerate(unknown):assert [jetbinding[var].get(m,F(0)) for m in LINEAR]==Ji[i][20:28]
    values={}
    for i in ancestors(D.nodes,frame+rho):
        n=D.nodes[i]
        if n[0]=='c':v={E0:F(n[1],n[2])} if n[1] else {}
        elif n[0]=='v':v=jetbinding[i]
        elif n[0]=='a':v=jadd(values[n[1]],values[n[2]])
        elif n[0]=='m':v=jmul(values[n[1]],values[n[2]])
        elif n[0]=='n':v=jscale(values[n[1]],-1)
        else:assert n[0]=='i';v=jinverse(values[n[1]])
        values[i]=v
    assert all(values[node]==({LINEAR[i-20]:F(1)} if i>=20 else {}) for i,node in enumerate(frame))
    rjets=[[values[i].get(m,F(0)) for m in jm] for i in rho]
    assert rjets==mat(c['all9_R0_mixed2jets'])
    kp=MIXED[-1];fold=values[rho[even]].get(kp,F(0));nquad=values[rho[ninth]].get(kp,F(0))
    assert fold==F(-3,4) and nquad==F(-3,2)
    mixed=[values[rho[ninth]].get(m,F(0))-nquad/fold*values[rho[even]].get(m,F(0)) for m in MIXED[:7]]
    assert mixed==[0,0,F(-1,2),0,2,0,0]==list(map(F,c['ninth_odd_linear_on7_parameters']))
    # Exact quadratic coefficient system, actual source E and all150 functions.
    cohen=c['actual_Cohen_R8_parameter_nodes'];S,T=c['quadratic_trace_S_T_nodes'];Qring=QuadraticReplay(D,S,T)
    one=(D.one,D.zero);pp=(D.zero,D.one);zero=(D.zero,D.zero)
    us=c['all56_coefficient_unknown_nodes'];others=[v for v in unknown if v!=pnode];pairs={pnode:pp}
    for j,var in enumerate(others):pairs[var]=(us[j],D.add(us[27+j],D.const(Ji[unknown.index(var)][27])))
    E=cuts+[D.add(rho[i],D.neg(cohen[j])) for j,i in enumerate(chosen)]+[D.add(rho[even],D.neg(cohen[7]))]
    qm={}
    for i in ancestors(D.nodes,E+rho+physical):
        n=D.nodes[i]
        if n[0]=='c':v=(D.const(F(n[1],n[2])),D.zero)
        elif n[0]=='v':v=pairs[i] if i in pairs else (i,D.zero)
        elif n[0]=='a':v=Qring.add(qm[n[1]],qm[n[2]])
        elif n[0]=='m':v=Qring.mul(qm[n[1]],qm[n[2]])
        elif n[0]=='n':v=Qring.neg(qm[n[1]])
        else:assert n[0]=='i';v=Qring.inverse(qm[n[1]])
        qm[i]=v
    equations=[v for root in E for v in qm[root]]
    assert equations==c['all56_actual_coefficient_Hensel_equations'] and len(equations)==56
    assert all(D.values[i]==0 for i in equations)
    A56=[[D.grad[i].get(v,F(0)) for v in us] for i in equations];Ai,adet=inverse(A56)
    assert A56==mat(c['unit56_Jacobian']) and adet==298598400
    B56=[[-D.grad[i].get(v,F(0)) for v in cohen] for i in equations];implicit=mm(Ai,B56)
    def over8(root):return [D.grad[root].get(v,F(0))+sum(D.grad[root].get(w,F(0))*implicit[k][j] for k,w in enumerate(us)) for j,v in enumerate(cohen)]
    p150=[qm[i] for i in physical]
    assert [list(qm[i]) for i in rho]==c['all9_R0_quadratic_p_readback_pairs']
    assert [list(v) for v in p150]==nc['all150_actual_quadratic_p_pairs']
    trace=[(D.add(a,D.mul(D.const(F(1,2)),D.mul(b,S))),b) for a,b in p150]
    assert [list(v) for v in trace]==nc['all150_actual_tracefree_pairs_even_odd']
    nu=D.add(T,D.mul(D.const(F(1,4)),D.mul(S,S)));assert nu==c['nu_tracefree_square_node']==nc['nu_actual_node']
    assert over8(nu)[7]==F(-4,3) and D.values[nu]==0
    assert trace[61][0]==c['ninth_actual_even_a_node'] and trace[61][1]==c['ninth_actual_odd_b_node']
    assert all(D.values[b]==0 for a,b in trace)
    brow=[over8(b) for a,b in trace];assert brow==mat(c['all150_actual_target_odd_linear_on_Cohen8'])
    sel=c['chosen8_odd_target_physical_rows'];assert sel==[61,18,19,28,29,45,46,60] and det([brow[i] for i in sel])==F(-4,9)
    assert brow[61]==[0,0,F(-1,2),0,2,0,0,0]
    z=[trace[i][1] for i in sel];assert z==c['chosen8_odd_coefficient_nodes']==nc['z1_through_z8_actual_odd_coefficient_nodes']
    assert [D.add(a,D.neg(b)) for a,b in zip(z,c['new_m8_coordinate_nodes_z'])]==c['unit8_R8_coordinate_change_equations']
    assert nc['actual_R8_physical_rows8']==[108,109,112,113,117,118,119,27]
    rho61=nc['actual_independent_R0_rho61_node'];vnode=D.add(rho61,D.neg(trace[61][0]));assert vnode==nc['actual_v_node']
    # Own36 source relations and84 chains, derived directly from Koszul/cone.
    pairs8=list(itertools.combinations(range(8),2));triples=list(itertools.combinations(range(8),3));pairidx={p:i for i,p in enumerate(pairs8)}
    d1=[]
    for i,j in pairs8:
        col={}
        if i==0:coladd(col,0,pscale(pmul(pvar(j),pvar(9)),-1));coladd(col,j,pvar(0))
        else:coladd(col,j,pvar(i));coladd(col,i,pscale(pvar(j),-1))
        d1.append(col)
    for i in range(8):
        col={0:pscale(pmul(pvar(8),pmul(pvar(0),pvar(i))),-1)}
        if i==0:coladd(col,0,pmul(pvar(9),pvar(9)))
        else:coladd(col,i,pvar(9))
        d1.append(col)
    d2=[]
    for a,b,d0 in triples:
        d2.append({pairidx[b,d0]:pvar(a),pairidx[a,d0]:pscale(pvar(b),-1),pairidx[a,b]:pvar(d0)})
    for idx,(i,j) in enumerate(pairs8):d2.append({idx:pvar(9),28+i:pvar(j),28+j:pscale(pvar(i),-1)})
    assert d1==[decode_col(v) for v in nc['d1_universal_sparse_polynomial_columns']]
    assert d2==[decode_col(v) for v in nc['d2_universal_sparse_polynomial_columns']]
    assert all(not compose(d1,col) for col in d2)
    def badd(a,b):return padd(a[0],b[0]),padd(a[1],b[1])
    def bmul(a,b):return padd(pmul(a[0],b[0]),pmul(pvar(8),pmul(a[1],b[1]))),padd(pmul(a[0],b[1]),pmul(a[1],b[0]))
    def peval(p):
        result=({},{});variables=[(pvar(i),{}) for i in range(9)]+[({},pvar(0))]
        for mon,c0 in p.items():
            value=({():c0},{})
            for j in mon:value=bmul(value,variables[j])
            result=badd(result,value)
        return result
    gens=[({():1},{})]+[({},pvar(i)) for i in range(1,8)]
    for col in d1:
        result=({}, {})
        for i,p in col.items():result=badd(result,bmul(peval(p),gens[i]))
        assert result==({}, {})
    actualvars=z+[nu,vnode]
    def at(p):
        terms=[]
        for mon,c0 in sorted(p.items()):
            value=D.const(c0)
            for j in mon:value=D.mul(value,actualvars[j])
            terms.append(value)
        return D.total(terms)
    assert [{str(i):at(p) for i,p in col.items()} for col in d1]==nc['R0_sparse_d1_8_by36_columns']
    assert [{str(i):at(p) for i,p in col.items()} for col in d2]==nc['R0_sparse_d2_36_by84_columns']
    bettiR=[9]+[math.comb(8,i+1) for i in range(1,8)]
    betti=[bettiR[i] if i<len(bettiR) else 0 for i in range(9)]
    for i in range(1,9):betti[i]+=bettiR[i-1]
    betti[0]-=1;betti[1]-=1
    assert betti==nc['N_full_minimal_free_resolution_ranks']==[8,36,84,126,126,84,36,9,1]
    assert all(all(mon for mon in p) for col in d1+d2 for p in col.values())
    output.mkdir(parents=True,exist_ok=False)
    evidence={'cut20_minor':'-96','actual_Cohen_rows8':[108,109,112,113,117,118,119,27],
        'rho9_firstjet_rank':7,'frame28_minor':'23040','quotient20_minor':'-480',
        'all17_monomial_source2jet_equations_verified':True,'fold_kappa2':'-3/4','rho61_kappa2':'-3/2',
        'unit56_minor':'298598400','nu_V_derivative':'-4/3','odd8_minor':'-4/9',
        'b61_firstjet':['0','0','-1/2','0','2','0','0','0'],
        'all150_source_p_and_trace_pairs_exact_and_odd_centers_zero':True,
        'actual_same_R0_N_equals_R8_plus_m8_tau':True,
        'z1_nilcurve_parameter_derivative_in_rho117':2,
        'all36_source_relations_and84_chains_rederived':True,
        'minimal_Betti':betti,'exactness_by_Koszul_monic_cone_not_only_d1d2_zero':True,
        'D4_extension_or_T3_or_origin_depth_solved':False}
    ep=output/'independent-proof.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
    for p in [Path(__file__),ROOT/'research/notes'/(BASE+'_REVIEW.md'),base/'report.json',base/'actual-N-R0-binding-DAG.json',cb/'report.json',cb/'actual-N-literal-cone-DAG.json',*gates,
        ROOT/'code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py',ROOT/'code/ambient_W1_S61_intersection_recognition_review_20261103_v63.py']:
        checked[str(p.relative_to(ROOT))]=sha(p)
    common={'all_passed':True,'actual28_source_and20cuts_same_R0_binding_independent':True,
        'restricted_QQ2jet_all17_monomials_and_true56_Hensel_pair_system_verified':True,
        'all150_actual_pairs_and_m8_odd_image_Nakayama_accepted':True,
        'actual_same_R0_N_R8_plus_m8_tau_and_v_z1tau_accepted':True,
        'literal36_84_and_full_minimal_Betti_cone_accepted':True,
        'D4_extension_T3_presentation_or_origin_depth_claim':False,
        'no_author_solver_new_source_jet_raise_GB_or_newpoint':True,
        'input_sha256':checked,'evidence_sha256':{ep.name:sha(ep)},
        'seconds':time.monotonic()-started,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'frozen':True}
    (output/'report.json').write_text(json.dumps({'id':BASE+'_REVIEW',**common},indent=2)+'\n')
    signal.alarm(0);print(json.dumps({k:common[k] for k in ['all_passed','seconds','max_rss_kib']}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
