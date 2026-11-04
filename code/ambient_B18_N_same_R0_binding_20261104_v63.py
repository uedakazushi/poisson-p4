#!/usr/bin/env python3
"""ONE actual N source/cuts/Cohen fold/odd-unit gate and finite quadratic input."""
import argparse
from fractions import Fraction as F
import gzip,hashlib,json
from pathlib import Path
import resource,signal,time
import sympy as sy
from ambient_B18_actual_top_coefficient_preflight_20261104_v63 import DAG,Quadratic
from ambient_B18_top_graph_telescoping_20261104_v63 import ancestors,forward,linear_solve,independent_rows
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63'
JM=[()]+[(i,) for i in range(8)]+[(i,7) for i in range(7)]+[(7,7)];JP={m:i for i,m in enumerate(JM)}
CONV=[[(i,j) for i,a in enumerate(JM) for j,b in enumerate(JM) if tuple(sorted(a+b))==m] for m in JM]
ZERO=(F(0),)*17
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(M):return [[str(x) for x in r] for r in M]
def jmul(a,b):return tuple(sum((a[i]*b[j] for i,j in terms if a[i] and b[j]),F(0)) for terms in CONV)
def jeval(nodes,needed,binding):
    v={}
    for i in needed:
        n=nodes[i];o=n[0]
        if o=='c':p=(F(n[1],n[2]),)+(F(0),)*16
        elif o=='v':p=binding[i]
        elif o=='a':p=tuple(a+b for a,b in zip(v[n[1]],v[n[2]]))
        elif o=='n':p=tuple(-a for a in v[n[1]])
        elif o=='m':p=jmul(v[n[1]],v[n[2]])
        elif o=='i':
            a=v[n[1]];assert a[0];h=(-F(1)/a[0]**2);hh=jmul((F(0),)+a[1:],(F(0),)+a[1:])
            p=(1/a[0],)+tuple(h*a[k]+hh[k]/a[0]**3 for k in range(1,17))
        else:raise ValueError(o)
        v[i]=p
    return v
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    tailp=ROOT/'reports/AMBIENT_B18_W1_TAIL_BINDING_20261104_v63/literal-tail-binding-DAG.json'
    tailr=tailp.parent/'report.json';bodyp=ROOT/'reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json'
    bundlep=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json';builderp=ROOT/'reports/AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63/builder.json'
    tail,body,bundle,builder=map(read,[tailp,bodyp,bundlep,builderp]);assert read(tailr)['all_passed']
    D=DAG();names=['p','s']+[f'k{i}' for i in range(9)];nv={n:D.var('N_SOURCE_'+n) for n in names}
    actual={n:D.add(nv[n],D.const(1 if n=='k0' else F(1,3) if n=='k5' else 0)) for n in names}
    group=[D.var(f'N_GROUP17_{i}') for i in range(17)];unknown=list(nv.values())+group;assert len(unknown)==28
    def poly(text):
        expr=sy.sympify(text);syms=sorted(expr.free_symbols,key=str)
        if not syms:return D.const(expr)
        terms=[]
        for mon,co in sy.Poly(expr,*syms).terms():
            v=D.const(co)
            for n,k in zip(syms,mon):
                for j in range(k):v=D.mul(v,actual[str(n)])
            terms.append(v)
        return D.total(terms)
    u=list(map(poly,body['families']['N']['full_61_coordinates']));q=[u[i] for i in bundle['q_indices_zero_based']]
    pp,ss=actual['p'],actual['s'];source49={'ss':ss,'aa':D.mul(pp,pp),'bb':D.mul(D.mul(pp,pp),pp),'cc':D.mul(pp,ss)}
    for j,i in enumerate(builder['common_q_selected_rows']):source49[f'v{j}']=q[i]
    for j,i in enumerate(builder['remaining_q_rows']):source49[f'z{j}']=q[i]
    binding={tail['old49_centered_variable_nodes'][n]:D.add(v,D.const(-1 if n in ['v4','v7'] else 0)) for n,v in source49.items()}
    binding.update(dict(zip(tail['actual_group17_parameter_nodes'],group)))
    oldnodes=tail['DAG_scalar_nodes'];roots=tail['all20_common_cut_equation_nodes']+tail['same_physical_R0_readback_nodes9']+tail['actual_full_physical150_nodes']
    needed=sorted(ancestors(oldnodes,roots));mp={}
    for i in needed:
        n=oldnodes[i];o=n[0]
        if o=='v':assert i in binding;v=binding[i]
        elif o=='c':v=D.const(F(n[1],n[2]))
        elif o=='n':v=D.neg(mp[n[1]])
        elif o=='a':v=D.add(mp[n[1]],mp[n[2]])
        elif o=='m':v=D.mul(mp[n[1]],mp[n[2]])
        elif o=='i':v=D.inverse(mp[n[1]])
        else:raise ValueError(o)
        mp[i]=v
    cuts=[mp[i] for i in tail['all20_common_cut_equation_nodes']];rho=[mp[i] for i in tail['same_physical_R0_readback_nodes9']]
    physical=[mp[i] for i in tail['actual_full_physical150_nodes']];assert all(D.values[i]==0 for i in cuts+rho)
    vals,grad=forward(D.nodes,set(unknown));C=[[grad[i].get(v,F(0)) for v in unknown] for i in cuts]
    piv=independent_rows([list(c) for c in zip(*C)]);assert len(piv)==20
    free=[i for i in range(28) if i not in piv];assert len(free)==8
    Z,det20=linear_solve([[r[i] for i in piv] for r in C],[[-r[i] for i in free] for r in C])
    def on8(node):return [grad[node].get(unknown[i],F(0))+sum(grad[node].get(unknown[p],F(0))*Z[k][j] for k,p in enumerate(piv)) for j,i in enumerate(free)]
    rfirst=[on8(i) for i in rho];chosen7=independent_rows(rfirst);assert len(chosen7)==7
    frame=cuts+[rho[i] for i in chosen7]+[nv['p']]
    J=[[grad[i].get(v,F(0)) for v in unknown] for i in frame]
    Ji,det28=linear_solve(J,[[F(int(i==j)) for j in range(28)] for i in range(28)]);assert det28
    target=[ZERO]*28
    for j in range(8):v=[F(0)]*17;v[JP[(j,)]]=1;target[20+j]=tuple(v)
    sol={i:ZERO for i in unknown};jneeded=sorted(ancestors(D.nodes,frame+rho))
    for degree in [1,2]:
        vv=jeval(D.nodes,jneeded,sol);cols=[i for i,m in enumerate(JM) if len(m)==degree]
        res=[[vv[i][j]-target[k][j] for j in cols] for k,i in enumerate(frame)]
        corr=[[-sum(Ji[i][k]*res[k][j] for k in range(28)) for j in range(len(cols))] for i in range(28)]
        for i,var in enumerate(unknown):
            a=list(sol[var])
            for j,col in enumerate(cols):a[col]+=corr[i][j]
            sol[var]=tuple(a)
    vv=jeval(D.nodes,jneeded,sol);assert all(vv[i]==target[j] for j,i in enumerate(frame))
    remain=[i for i in range(9) if i not in chosen7];k2=JP[(7,7)];even=next(i for i in remain if vv[rho[i]][k2]);ninth=next(i for i in remain if i!=even)
    c=vv[rho[even]][k2];d=vv[rho[ninth]][k2]
    blinear=[vv[rho[ninth]][JP[(j,7)]]-d/c*vv[rho[even]][JP[(j,7)]] for j in range(7)]
    assert any(blinear)
    # Nonnormal quotient's old free3 source block is killed transversely by20cuts.
    critvars=group+[nv[n] for n in ['k0','k2','k5']]
    critJ=[[grad[i].get(v,F(0)) for v in critvars] for i in cuts]
    _,critdet=linear_solve(critJ,[[F(int(i==j)) for j in range(20)] for i in range(20)]);assert critdet
    # Actual rank2 coefficient Hensel system over eight selected R0 coordinates.
    base=[D.var(f'N_COHEN_X{i}') for i in range(7)]+[D.var('N_COHEN_V')]
    E=cuts+[D.add(rho[i],D.neg(base[j])) for j,i in enumerate(chosen7)]+[D.add(rho[even],D.neg(base[7]))]
    Q=Quadratic(D,'N_EXACT_FOLD');others=[i for i in unknown if i!=nv['p']]
    pairs={nv['p']:Q.p};Avec=[];Bvec=[]
    for var in others:
        a=D.var(f'N_COEF_A_{var}');b=D.var(f'N_COEF_B_{var}')
        pairs[var]=(a,D.add(b,D.const(sol[var][JP[(7,)]])));Avec.append(a);Bvec.append(b)
    us=Avec+Bvec+[Q.S,Q.T];assert len(us)==56
    qroots=E+rho+physical; qneeded=sorted(ancestors(D.nodes[:],qroots)); qmap={}
    for i in qneeded:
        n=D.nodes[i];o=n[0]
        if o=='c':z=Q.scalar(D.const(F(n[1],n[2])))
        elif o=='v':z=pairs[i] if i in pairs else Q.scalar(i)
        elif o=='a':z=Q.add(qmap[n[1]],qmap[n[2]])
        elif o=='n':z=Q.neg(qmap[n[1]])
        elif o=='m':z=Q.mul(qmap[n[1]],qmap[n[2]])
        elif o=='i':z=Q.inv(qmap[n[1]])
        else:raise ValueError(o)
        qmap[i]=z
    eq56=[a for i in E for a in qmap[i]];assert len(eq56)==56 and all(D.values[i]==0 for i in eq56)
    nu=D.add(Q.T,D.mul(D.const(F(1,4)),D.mul(Q.S,Q.S)))
    q9=qmap[rho[ninth]];traceeven=D.add(q9[0],D.mul(D.const(F(1,2)),D.mul(q9[1],Q.S)))
    _,g56=forward(D.nodes,set(us+base))
    A=[[g56[i].get(v,F(0)) for v in us] for i in eq56];B=[[-g56[i].get(v,F(0)) for v in base] for i in eq56]
    ux,det56=linear_solve(A,B);assert det56
    def bd(node):return [g56[node].get(v,F(0))+sum(g56[node].get(w,F(0))*ux[k][j] for k,w in enumerate(us)) for j,v in enumerate(base)]
    brow=[bd(qmap[i][1]) for i in physical];assert all(D.values[qmap[i][1]]==0 for i in physical)
    assert brow[tail['same_physical_R0_rows9'][ninth]][:7]==blinear
    first=tail['same_physical_R0_rows9'][ninth];basisrows=[first]
    for i,row in enumerate(brow):
        if i not in basisrows and len(independent_rows([brow[j] for j in basisrows]+[row]))>len(basisrows):basisrows.append(i)
        if len(basisrows)==8:break
    assert len(basisrows)==8
    _,bdet=linear_solve([brow[i] for i in basisrows],[[F(int(i==j)) for j in range(8)] for i in range(8)]);assert bdet
    bnodes=[qmap[physical[i]][1] for i in basisrows];z=[D.var(f'N_FINAL_z{i+1}') for i in range(8)]
    zchange=[D.add(b,D.neg(v)) for b,v in zip(bnodes,z)]
    data={'DAG_scalar_nodes':D.nodes,'normalization_source28_variable_nodes':unknown,
          'oldN_source11_parameter_values':actual,'all20_actual_physical_cut_equations':cuts,'physical150_N_source_nodes':physical,
          'actual_same_R0_rows9':tail['same_physical_R0_rows9'],'actual_same_R0_readback_nodes9':rho,
          'cut20_Jacobian':enc(C),'cut20_pivot_columns':piv,'cut20_free_columns':free,'cut20_minor':str(det20),
          'all9_R0_firstjets_on_source8':enc(rfirst),'chosen7_R0_indices':chosen7,'Cohen_even8_R0_index':even,'ninth_R0_index':ninth,
          'frame28_equation_nodes':frame,'frame28_Jacobian':enc(J),'frame28_minor':str(det28),
          'mixed_jet_monomials':list(map(list,JM)),'all28_mixed_source2jets':[[i,list(map(str,sol[i]))] for i in unknown],
          'all9_R0_mixed2jets':[list(map(str,vv[i])) for i in rho],
          'fold_quadratic_coefficient':str(c),'ninth_quadratic_coefficient':str(d),'ninth_odd_linear_on7_parameters':list(map(str,blinear)),
          'nonnormal_quotient_free20_source_parameters':critvars,'nonnormal_quotient_cut20_minor':str(critdet),
          'actual_Cohen_R8_parameter_nodes':base,'quadratic_trace_S_T_nodes':[Q.S,Q.T],
          'all56_coefficient_unknown_nodes':us,'all56_actual_coefficient_Hensel_equations':eq56,
          'unit56_Jacobian':enc(A),'unit56_minor':str(det56),
          'all9_R0_quadratic_p_readback_pairs':[qmap[i] for i in rho],
          'nu_tracefree_square_node':nu,'ninth_actual_even_a_node':traceeven,'ninth_actual_odd_b_node':q9[1],
          'all150_actual_target_odd_linear_on_Cohen8':enc(brow),'chosen8_odd_target_physical_rows':basisrows,
          'chosen8_odd_coefficient_nodes':bnodes,'chosen8_odd_linear_determinant':str(bdet),
          'new_m8_coordinate_nodes_z':z,'unit8_R8_coordinate_change_equations':zchange,
          'actual_ninth_R0_v_after_even_shift':'rho_ninth-a(R8)=z1*tau, tau=p-S/2',
          'exact_actual_N_target_image':'R8 + m8*tau, tau^2=nu; followsalloddcenters0+unit8oddJacobian Nakayama',
          'same_R0_identification':'R0=R8[[v]], v=z1*tau, R8actual8ofR0followedbyunit8change',
          'no_T3_presentation_extension_lift_or_B18_depth_solution':True}
    encoded=json.dumps(data,separators=(',',':'))+'\n';assert len(encoded.encode())<32*1024*1024
    (out/'actual-N-R0-binding-DAG.json').write_text(encoded)
    paths=[Path(__file__),tailp,tailr,bodyp,bundlep,builderp]
    rep={'id':STEM,'all_passed':True,'status':'actual_N_same_R0_Cohen8_fold_and_m8_odd_image_bound',
         'normalization28_cuts20_rank':20,'R0_nine_firstjet_rank':7,'fold_quadratic_unit':str(c),
         'ninth_actual_odd_linear_nonzero':True,'unit56_minor':str(det56),'unit8_odd_minor':str(bdet),
         'chosen_Cohen_R0_rows8':[tail['same_physical_R0_rows9'][i] for i in chosen7+[even]],
         'ninth_R0_physical_row':tail['same_physical_R0_rows9'][ninth],
         'actual_R0_model_R8_v_and_N_R8_m8_tau_proved':True,
         'nilcurve_extension_lifts_Jinf_syzygy_T3_or_depth_solved':False,
         'no_GB_new_points_source_jets_above2_or_opaque_branch_pullback':True,
         'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},
         'evidence_sha256':{'actual-N-R0-binding-DAG.json':sha(out/'actual-N-R0-binding-DAG.json')},
         'DAG_nodes':len(D.nodes),'certificate_bytes':len(encoded.encode()),'seconds':time.monotonic()-start,
         'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps({k:rep[k] for k in ['all_passed','R0_nine_firstjet_rank','fold_quadratic_unit','unit56_minor','unit8_odd_minor','chosen_Cohen_R0_rows8','ninth_R0_physical_row','DAG_nodes','certificate_bytes','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
